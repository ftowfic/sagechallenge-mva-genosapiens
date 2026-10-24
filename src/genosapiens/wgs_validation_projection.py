"""Operation-owned acceptance queries over complete normalized scientific evidence.

Only source/member joins are spilled. Hypothesis groups are traversed through
complete orders and scalar outcomes; no complete rich context or H-sized JSON
collection is constructed. Evaluator truth is supplied only to target_recovery,
when a completed execution has already been authenticated.
"""
from collections import Counter
from contextlib import ExitStack
from dataclasses import dataclass
import json

from genosapiens.evidence_model import _sha
from genosapiens.evidence_serialization import _canonical
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_scientific_artifact import ScientificArtifactReader
from genosapiens.wgs_scientific_products import product_descriptors, validate_scientific_products
from genosapiens.wgs_validation_audit import iter_validation_audit

MAX_GROUP_BYTES = 64 * 1024**2
MAX_RECOVERY_METADATA_BYTES = 8 * 1024**2


@dataclass(frozen=True, slots=True, repr=False)
class ValidationRecoveryBinding:
    scientific_manifest_sha256: str
    audit_manifest_sha256: str
    order_manifest_sha256: str
    selection_manifest_sha256: str
    population_selector: tuple[str, str, str | None] | None
    schema: str = "genosapiens.wgs-normalized-validation-recovery"
    version: int = 1

    def __post_init__(self):
        for value in (self.scientific_manifest_sha256, self.audit_manifest_sha256,
                self.order_manifest_sha256, self.selection_manifest_sha256):
            _sha(value)
        selector = self.population_selector
        if (self.schema != "genosapiens.wgs-normalized-validation-recovery"
                or type(self.version) is not int or self.version != 1
                or selector is not None and (type(selector) is not tuple or len(selector) != 3
                    or any(type(v) is not str or not v for v in selector[:2])
                    or selector[2] is not None and (type(selector[2]) is not str or not selector[2]))):
            raise EvidenceStoreError("invalid normalized validation recovery binding")


@dataclass(frozen=True, slots=True)
class TargetRecovery:
    admitted_member_count: int
    generated_member_count: int
    exact_generated: bool
    generated_rank: int | None
    clinical_generated: bool
    clinical_rank: int | None
    population_excluded: bool
    route_recall: tuple[tuple[str, int], ...]


class NormalizedValidationProjection:
    def __init__(self, reader, products, audit, *, audit_directory, epcr_assignment, population_selector=None):
        from genosapiens.track1_validation_runner import (
            MODEL_COVERAGE_KEYS, POPULATION_LOOKUP_KEYS, _member_keys,
            _population_frequency_observations, _population_lookup_counts,
        )
        if type(reader) is not ScientificArtifactReader:
            raise EvidenceStoreError("verified normalized scientific owner required")
        reader._require_open()
        validate_scientific_products(reader, product_descriptors(products),
            expected_prefilter_manifest_sha256=products.primary.prefilter_manifest_sha256)
        if reader.artifact != products.scientific:
            raise EvidenceStoreError("foreign validation scientific owner")
        self.reader, self.products, self.connection = reader, products, reader.reader.connection
        self.tables, self.closed, self.stack = [], False, ExitStack()
        self.binding = ValidationRecoveryBinding(products.scientific.sha256, audit.sha256,
            products.order.sha256, products.selection.sha256, population_selector)
        self.coverage = dict.fromkeys(MODEL_COVERAGE_KEYS, 0)
        self.lookups = dict.fromkeys(POPULATION_LOOKUP_KEYS, 0)
        self.associations, self.outcomes = Counter(), Counter()
        self.prefix = "normalized_validation"
        try:
            self.stack.enter_context(reader.reader.read_scope())
            if (epcr_assignment.selection_sha256 != products.selection.sha256
                    or tuple((r.hypothesis_id, r.source_rank, r.source_score) for r in epcr_assignment.values)
                    != tuple((r.primary_hypothesis_id, r.source_rank, r.source_score) for r in products.selection.selected)):
                raise EvidenceStoreError("validation selected EPCR binding changed")
            for suffix, definition in (
                ("inputs", "id TEXT PRIMARY KEY"),
                ("members", "allele TEXT,member TEXT,PRIMARY KEY(allele,member)"),
                ("admitted", "allele TEXT PRIMARY KEY"),
                ("required", "context TEXT PRIMARY KEY"),
                ("contexts", "context TEXT PRIMARY KEY,ordinal INTEGER UNIQUE"),
                ("annotations", "context TEXT PRIMARY KEY,clinical INTEGER"),
                ("excluded", "allele TEXT PRIMARY KEY"),
                ("routes", "route TEXT,member TEXT,PRIMARY KEY(route,member)"),
                ("frequencies", "identity TEXT,frequency TEXT,an INTEGER,PRIMARY KEY(identity,frequency,an)"),
                ("unresolved", "kind TEXT,id TEXT,PRIMARY KEY(kind,id)"),
            ):
                table = self.prefix + "_" + suffix
                self.connection.execute(f"CREATE TABLE {table} ({definition})")
                self.tables.append(table)
            self.connection.execute(f"CREATE INDEX {self.prefix}_member_key ON {self.prefix}_members(member,allele)")
            parents = dict(products.scientific.parent_bindings)
            audit_parents = {key: parents[key] for key in ("prepared_case", "policy_bindings", "execution_config")}
            audit_parents["scientific"] = products.scientific.sha256
            audit_parents.update({k:v for k,v in parents.items() if k in {
                'analysis_eligibility_policy','analysis_eligibility_manifest','analysis_eligibility_rescue_manifest'}})
            for row in iter_validation_audit(audit, directory=audit_directory,
                    execution_config=products.execution_config, expected_parents=audit_parents, scientific_reader=reader):
                if row.kind == "admission":
                    for value in row.admission.inputs:
                        self.connection.execute(f"INSERT INTO {self.prefix}_inputs VALUES (?)", (value.subject_id,))
                        if value.allele is not None:
                            self.connection.executemany(f"INSERT INTO {self.prefix}_members VALUES (?,?)",
                                ((value.allele.allele_id, member) for member in _member_keys(value.allele)))
                    for decision in row.admission.decisions:
                        if not self.connection.execute(f"SELECT 1 FROM {self.prefix}_inputs WHERE id=?", (decision.subject_id,)).fetchone():
                            raise EvidenceStoreError("validation decision has no admission input")
                        self.associations[decision.association_state] += 1
                        self.outcomes[decision.outcome] += 1
                        if decision.outcome == "admitted":
                            if not decision.context_ids:
                                raise EvidenceStoreError("validation admission omitted its contexts")
                            self.connection.execute(f"INSERT OR IGNORE INTO {self.prefix}_admitted VALUES (?)", (decision.subject_id,))
                            self.connection.executemany(f"INSERT OR IGNORE INTO {self.prefix}_required VALUES (?)",
                                ((key,) for key in decision.context_ids))
                else:
                    annotation = row.annotation
                    self.connection.execute(f"INSERT INTO {self.prefix}_annotations VALUES (?,?)",
                        (annotation.context_id, int(annotation.clinical_completed)))
                    for batch in reader.annotation_reports(annotation):
                        from types import SimpleNamespace
                        report = SimpleNamespace(source_fields=batch["population"]["source_fields"])
                        for key, count in _population_lookup_counts(SimpleNamespace(report=report)):
                            self.lookups[key] += count
            for row in reader.case.contexts():
                if row["member_count"]:
                    self.connection.execute(f"INSERT INTO {self.prefix}_contexts VALUES (?,?)", (row["context_id"], row["ordinal"]))
            for left, right in (("required", "contexts"), ("contexts", "required"), ("contexts", "annotations"), ("annotations", "contexts")):
                if self.connection.execute(f"SELECT 1 FROM {self.prefix}_{left} l LEFT JOIN {self.prefix}_{right} r "
                        "ON l.context=r.context WHERE r.context IS NULL LIMIT 1").fetchone():
                    raise EvidenceStoreError("validation admission, annotation and scientific context closure differs")
            for row in reader.iter_member_records():
                allele = row["member"].allele
                keys = _member_keys(allele)
                saved = tuple(value for value, in self.connection.execute(
                    f"SELECT member FROM {self.prefix}_members WHERE allele=? ORDER BY member", (allele.allele_id,)))
                if keys != saved:
                    raise EvidenceStoreError("validation scientific source members changed")
                self.connection.executemany(f"INSERT OR IGNORE INTO {self.prefix}_routes VALUES (?,?)",
                    ((route, member) for route in row["seed"].routes for member in keys))
                self.connection.executemany(f"INSERT OR IGNORE INTO {self.prefix}_frequencies VALUES (?,?,?)",
                    _population_frequency_observations((allele,), population_selector))
            # Exclusions are global and include only actually emitted hypotheses.
            # They must close before search-level unresolved/contradiction counts.
            for value in reader.iter_scalar_records():
                subjects = reader.case.identities.subject_ids(value["ordering"]["locator"])
                if len(subjects) != len(value["population_excluded_members"]):
                    raise EvidenceStoreError("validation scalar subject coverage changed")
                self.connection.executemany(f"INSERT OR IGNORE INTO {self.prefix}_excluded VALUES (?)",
                    ((key,) for key, excluded in zip(subjects, value["population_excluded_members"], strict=True) if excluded))
                self.coverage["completed-population-lookup-unknown-rescue"] += value["unknown_rescue_members"]
                if (len(subjects) == 2 and value["model_state"] != "contradicted"
                        and "unresolved" in value["member_states"] and not any(value["clinical_support"])):
                    self.coverage["clinical-unasserted-exploratory-pair"] += 1
                if (value["completeness"] == "unknown" and value["model_state"] != "contradicted"
                        and all(state is True for state in value["alternate_present"])):
                    self.coverage["unknown-count-hypothesis"] += 1
            self._search_coverage()
            self.eligible, self.clinical_eligible = [], []
            for group in products.selection.selected:
                ordered = reader.case.ordered_at(group.source_rank)
                scalar = reader.case.scalar_facts.get(ordered["locator"])
                value = self.group(ordered["locator"])
                self.eligible.append(value)
                if self._clinical(scalar):
                    self.clinical_eligible.append(value)
            self.eligible, self.clinical_eligible = tuple(self.eligible), tuple(self.clinical_eligible)
            self.epcr = tuple(value.epcr for value in epcr_assignment.values)
        except BaseException:
            self.close()
            raise

    def _root(self, value, kind):
        root = self.reader.reader.reference(value["name"], kind=kind)
        if root.sha256 != value["sha256"]:
            raise EvidenceStoreError("validation scientific dependency changed")
        return root

    def _excluded(self, allele):
        return self.connection.execute(f"SELECT 1 FROM {self.prefix}_excluded WHERE allele=?", (allele,)).fetchone() is not None

    def _search_coverage(self):
        evidence = self.reader.reader
        self.completed = 0
        self.coverage["population-annotation-completed"] = 1
        self.coverage["clinical-annotation-completed"] = 1
        for row in self.reader.case.contexts():
            if not row["member_count"]:
                continue
            context = evidence.reconstruct(self._root(row["context"], "scientific-case-context-v1"))
            universe = evidence.reconstruct(self._root(context["universe"], "complete-hypothesis-universe-v1"))
            search = evidence.reconstruct(self._root(universe["search"], "complete-partner-search-v1"))
            if search["enumeration_status"] != "complete":
                raise EvidenceStoreError("normalized validation search is incomplete")
            self.completed += 1
            clinical, = self.connection.execute(f"SELECT clinical FROM {self.prefix}_annotations WHERE context=?",
                (row["context_id"],)).fetchone()
            self.coverage["clinical-annotation-completed"] &= clinical
            for member in evidence.iter_array(self._root(search["member_qualifications"], "search-member-qualifications-v1")):
                if member["state"] == "unresolved" and self._excluded(member["allele_id"]):
                    self.coverage["population-excluded-unresolved-partner"] += 1
            for kind, values in (
                    ("source", evidence.iter_array(self._root(search["unresolved"], "search-unresolved-sources-v1"))),
                    ("relevance", iter(search["unresolved_relevance_source_allele_ids"]))):
                self.connection.executemany(f"INSERT OR IGNORE INTO {self.prefix}_unresolved VALUES (?,?)",
                    ((kind, key) for key in values))
            for locator, value in self.reader.case.scalar_facts.pair_records(row["ordinal"]):
                if value["model_state"] == "contradicted" and "nonqualifying" in value["member_states"]:
                    subjects = self.reader.case.identities.subject_ids((row["ordinal"], *locator))
                    if any(self._excluded(key) for key in subjects):
                        self.coverage["independent-non-population-contradiction"] += 1
            availability, availability_bytes = {}, 0
            # Only conditional singleton subjects are needed in memory. The
            # original singleton records have the existing bounded row limit.
            for value in evidence.iter_array(self._root(universe["singletons"], "constructed-singletons-v1")):
                if value["conditional_singleton_search_satisfied"]:
                    hypothesis = value["hypothesis"]
                    for subject in (*hypothesis["allele_refs"], *hypothesis["locus_call_refs"]):
                        key = subject["allele_id"]
                        if key not in availability:
                            availability_bytes += len(key) * 4 + 128
                            if availability_bytes > MAX_GROUP_BYTES:
                                raise EvidenceStoreError("validation singleton window exceeds capacity", capacity=True)
                            availability[key] = False
            if availability and universe["pair_evidence"] is not None:
                pairs = evidence.reconstruct(self._root(universe["pair_evidence"], "complete-pair-evidence-v1"))
                for value in evidence.iter_array(self._root(pairs["availability"], "pair-partner-availability-v1")):
                    if value["allele_id"] in availability:
                        availability[value["allele_id"]] = value["has_not_nonqualifying_partner"]
            self.coverage["unresolved-partner-singleton-violation"] += sum(availability.values())
        if self.completed != self.products.scientific.positive_context_count:
            raise EvidenceStoreError("validation completed context count changed")

    @staticmethod
    def _clinical(value):
        return bool(value["member_states"]) and all(state == "qualifying" and supported
            for state, supported in zip(value["member_states"], value["clinical_support"], strict=True))

    def group(self, locator):
        self._require_open()
        values = set()
        size = 0
        for allele in self.reader.case.identities.subject_ids(locator):
            for member, in self.connection.execute(f"SELECT member FROM {self.prefix}_members WHERE allele=? ORDER BY member", (allele,)):
                if member not in values:
                    size += len(member) * 4 + 128
                    if size > MAX_GROUP_BYTES:
                        raise EvidenceStoreError("validation member group exceeds capacity", capacity=True)
                    values.add(member)
        if not values:
            raise EvidenceStoreError("scientific hypothesis has no evaluator source identity")
        return tuple(sorted(values))

    def target_recovery(self, target):
        self._require_open()
        if type(target) is not tuple or not target or target != tuple(sorted(set(target))):
            raise EvidenceStoreError("invalid evaluator target members")
        if len(_canonical(target)) > MAX_GROUP_BYTES:
            raise EvidenceStoreError("evaluator target exceeds capacity", capacity=True)
        target_set, generated = set(target), set()
        first = clinical_first = None
        excluded = False
        clinical_rank = 0
        for value in self.reader.iter_scalar_records():
            group = self.group(value["ordering"]["locator"])
            generated.update(target_set.intersection(group))
            clinical = self._clinical(value)
            clinical_rank += int(clinical)
            if group == target:
                if first is None:
                    first = value["ordering"]["rank"]
                if clinical and clinical_first is None:
                    clinical_first = clinical_rank
                excluded |= value["population_outcome"] == "population_excluded"
        admitted = sum(self.connection.execute(f"SELECT 1 FROM {self.prefix}_members m JOIN {self.prefix}_admitted a "
            "ON m.allele=a.allele WHERE m.member=? LIMIT 1", (key,)).fetchone() is not None for key in target)
        routes = []
        for route, in self.connection.execute(f"SELECT DISTINCT route FROM {self.prefix}_routes ORDER BY route"):
            recall = sum(self.connection.execute(f"SELECT 1 FROM {self.prefix}_routes WHERE route=? AND member=?",
                (route, key)).fetchone() is not None for key in target)
            routes.append((route, recall))
        return TargetRecovery(admitted, len(generated), first is not None, first,
            clinical_first is not None, clinical_first, excluded, tuple(routes))

    def frequency_observations(self, identities=None):
        self._require_open()
        if identities is None:
            yield from self.connection.execute(f"SELECT * FROM {self.prefix}_frequencies ORDER BY identity,frequency,an")
        else:
            for key in sorted(set(identities)):
                _sha(key)
                yield from self.connection.execute(f"SELECT * FROM {self.prefix}_frequencies WHERE identity=? ORDER BY frequency,an", (key,))

    def metrics(self):
        self._require_open()
        count = self.products.scientific.positive_context_count
        def unresolved(kind):
            return None if not count else self.connection.execute(
                f"SELECT count(*) FROM {self.prefix}_unresolved WHERE kind=?", (kind,)).fetchone()[0]
        value = dict(status="succeeded", recovery_binding=self.binding,
            eligible_member_groups=self.eligible, eligible_epcr=self.epcr,
            clinical_only_eligible_member_groups=self.clinical_eligible,
            candidate_burden=count, hypothesis_count=self.products.scientific.hypothesis_count,
            required_enumerations=count, completed_enumerations=self.completed,
            association_state_counts=tuple(sorted(self.associations.items())),
            admission_outcome_counts=tuple(sorted(self.outcomes.items())),
            model_coverage_counts=tuple(sorted(self.coverage.items())), source_inventory_status="complete",
            broader_unresolved_source_count=unresolved("source"),
            broader_unresolved_relevance_source_count=unresolved("relevance"),
            population_lookup_counts=tuple(sorted(self.lookups.items())))
        if len(_canonical(value)) > MAX_RECOVERY_METADATA_BYTES:
            raise EvidenceStoreError("validation metadata exceeds capacity", capacity=True)
        return value

    def _require_open(self):
        if self.closed:
            raise EvidenceStoreError("normalized validation projection is closed")
        self.reader._require_open()

    def close(self):
        for table in reversed(self.tables):
            self.connection.execute("DROP TABLE " + table)
        self.tables.clear()
        self.stack.close()
        self.closed = True

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, *_):
        self.close()
