"""Immutable diagnosis-guided construction output and source-specific admission.

Raw VCF/index/reference identities are construction attestations. Synthesis
rehashes only the declared small companions; it never opens original sources.
The content bindings detect drift, not malicious fabrication by an authorized
builder/reviewer. Scientific support is evaluated separately from integrity.
"""
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path

from genosapiens.json_schema import validate_json_schema
from genosapiens.track2_input import AnalysisOrigin, Track2InputError
from genosapiens.track2_mechanism import assess_mechanism_readiness
from genosapiens.track2_published_genotype import canonical, sha256, _source_file

ROOT = Path(__file__).resolve().parents[2]
SCHEMA_PATH = ROOT / "schemas/track2/diagnosis-guided-package.schema.json"
SCHEMA = "genosapiens.track2-diagnosis-guided-package/v5"
ROLES = ("scope", "selection", "structural", "variants", "scenarios", "phenotype", "evidence")


def input_acceptance(manifest, parts, timestamp):
    from genosapiens.track2_adjudication import automatic_decision, AUTOMATIC_CHECKS
    validate_components(parts,classification=manifest['classification'],evidence_origin=manifest['evidence_origin'])
    evidence={'artifacts':manifest['artifacts'],'scope_policy':parts['scope']['policy_sha256'],
        'selection_policy':parts['selection']['coverage']['filter_policy']}
    return automatic_decision(domain='input',subject=manifest,evidence=evidence,
        checks={k:{'state':'accepted','reasons':['validated_'+k]} for k in AUTOMATIC_CHECKS['input']},
        execution={'code_sha256':manifest['builder']['code_sha256'],'image':manifest['builder']['runtime']['image']},
        timestamp=timestamp)


def _require(condition, reason):
    if not condition:
        raise Track2InputError(reason)


def _digest(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def _fields(value, fields):
    _require(isinstance(value, dict) and set(value) == set(fields.split()),
             "diagnosis package record fields differ from the current contract")


def _index(rows, key):
    _require(isinstance(rows, list) and all(isinstance(r, dict) and r.get(key) for r in rows),
             "diagnosis package record identities are missing")
    indexed = {row[key]: row for row in rows}
    _require(len(indexed) == len(rows), "diagnosis package record identities are duplicated")
    return indexed


def scenario_identity(scenario):
    """Bind the actual genotype and its evidence; HPO is contextual only."""
    return "diagnosis-genotype:" + sha256(canonical({
        key: value for key, value in scenario.items() if key not in {"scenario_id", "mechanism_inputs"}
    }))[:24]


def validate_components(parts, *, classification, evidence_origin):
    _require(set(parts) == set(ROLES), "diagnosis package companions are incomplete")
    scope, selection = parts["scope"], parts["selection"]
    _fields(scope, "schema original_diagnosis resolved_diagnosis reference_sha256 resource_bindings "
                   "policy_sha256 genes intervals assertions mapping_evidence limitations")
    _require(scope["schema"] == "genosapiens.track2-diagnosis-scope/v1"
             and scope["original_diagnosis"] and scope["resolved_diagnosis"]
             and _digest(scope["reference_sha256"]) and _digest(scope["policy_sha256"])
             and scope["resource_bindings"] and all(_digest(h) for h in scope["resource_bindings"].values())
             and scope["genes"] and scope["intervals"] and scope["assertions"] and scope["mapping_evidence"],
             "diagnosis scope or essential source bindings are missing")
    _fields(selection, "schema classification evidence_origin source_sha256 sample reference_sha256 "
                       "index_sha256 scope_sha256 subset_sha256 structural_sha256 strategy coverage selected_alleles vcf_policy_sha256 vcf_evidence")
    from genosapiens.vcf_contracts import require_current_policy
    require_current_policy(selection['vcf_policy_sha256'])
    from genosapiens.track2_vcf_evidence import validate as validate_vcf_evidence
    validate_vcf_evidence(selection['vcf_evidence'],selection,parts['structural'])
    from genosapiens.track2_vcf_evidence import validate_calls
    validate_calls(parts)
    _require(selection["schema"] == "genosapiens.track2-source-selection/v4"
             and selection["classification"] == classification and selection["evidence_origin"] == evidence_origin
             and _digest(selection["source_sha256"]) and selection["sample"]
             and selection["reference_sha256"] == scope["reference_sha256"]
             and selection["scope_sha256"] == sha256(canonical(scope))
             and _digest(selection["subset_sha256"])
             and selection["strategy"] in {"indexed", "streaming"}
             and (selection["index_sha256"] is None or _digest(selection["index_sha256"]))
             and (selection["strategy"] != "indexed" or selection["index_sha256"] is not None),
             "source selection does not bind diagnosis, sample, reference or subset")
    coverage = selection["coverage"]
    _fields(coverage, "state selected_allele_count source_record_count source_alternate_class_counts source_accounting unsupported_classes limitations "
                      "filter_policy filter_exclusions filter_excluded_record_count filter_excluded_alternate_count retention deferred_inventory")
    from genosapiens.track2_retention import validate_summary
    validate_summary(coverage['retention'], selection, scope)
    from genosapiens.track2_vcf_selection import ACCOUNTING
    accounting=coverage['source_accounting']
    _require(set(accounting)=={'unit',*ACCOUNTING} and accounting['unit']=='source_ALT_occurrences'
        and all(type(accounting[k]) is int and accounting[k]>=0 for k in ACCOUNTING)
        and sum(accounting[k] for k in ACCOUNTING)==sum(coverage['source_alternate_class_counts'].values())
        and accounting['quality_excluded']==coverage['filter_excluded_alternate_count']
        and accounting['sequence_scoped']==sum(coverage['retention']['occurrences'][k] for k in ('retained','deferred')),
        'Source ALT accounting does not reconcile with retained and deferred selection.')
    from genosapiens.track2_deferred import validate_descriptor as validate_deferred
    validate_deferred(coverage['deferred_inventory'], selection)
    from genosapiens.track2_vcf_selection import POLICY
    from genosapiens.track2_exclusions import validate_descriptor
    _require(coverage['filter_policy']==POLICY, 'FILTER selection accounting differs.')
    validate_descriptor(coverage['filter_exclusions'],selection=selection)
    _require(coverage["state"] == "complete_within_declared_scope"
             and type(coverage["selected_allele_count"]) is int and coverage["selected_allele_count"] >= 0
             and type(coverage["source_record_count"]) is int and coverage["source_record_count"] >= 0
             and isinstance(coverage["source_alternate_class_counts"], dict)
             and all(type(v) is int and v >= 0 for v in coverage["source_alternate_class_counts"].values())
             and set(coverage["unsupported_classes"]) >= {"SV", "CNV", "repeat", "symbolic"}
             and coverage["limitations"], "essential selection coverage is incomplete")
    selected = _index(selection["selected_alleles"], "allele_id")
    _require(len(selected) == coverage["selected_allele_count"], "selection coverage count disagrees")
    _require(len({r['source_record_id'] for r in selected.values()}) <= coverage['source_record_count']
             and len(selected) <= sum(coverage['source_alternate_class_counts'].values()),
             "selection coverage cannot exceed observed source records")
    for row in selected.values():
        _fields(row, "allele_id source_record_id source_allele_id")
        _require(row["source_record_id"] and row["source_allele_id"], "source allele provenance is missing")
    _require(selection["subset_sha256"] == sha256(canonical(selection["selected_alleles"])),
             "selected allele subset changed")
    _require(selection['structural_sha256'] == sha256(canonical(parts['structural'])),
             'structural evidence checksum changed')
    from genosapiens.track2_structural import validate_component
    try:
        validate_component(parts['structural'], selection, scope)
    except (ValueError, KeyError, TypeError, IndexError, StopIteration) as error:
        raise Track2InputError('structural evidence is invalid') from error
    variants = _index(parts["variants"], "allele_id")
    _require(set(variants) <= set(selected), "interpreted allele is outside source selection")
    evidence = _index(parts["evidence"], "evidence_id")
    for row in evidence.values():
        _fields(row, "evidence_id source_identity lineage_key assertion evidence_role")
        _require(all(isinstance(v, str) and v for v in row.values()), "evidence provenance is missing")
        _require(row["evidence_role"] in {"source_observation", "prediction", "curated", "synthetic"},
                 "evidence role is undefined")
        _require(row["evidence_role"] != "synthetic" or evidence_origin == "independently_synthetic",
                 "synthetic evidence cannot support actual case science")
    for row in variants.values():
        _fields(row, "allele_id source_record_id source_allele_id gene_ids assessment evidence_ids")
        _require(all(row[k] == selected[row["allele_id"]][k] for k in selected[row["allele_id"]])
                 and set(row["gene_ids"]) <= set(scope["genes"])
                 and row["assessment"] in {"admitted", "uncertain", "excluded"}
                 and row["evidence_ids"] and set(row["evidence_ids"]) <= set(evidence),
                 "interpreted allele lacks source-bound evidence")
    scenarios = _index(parts["scenarios"], "scenario_id")
    for row in scenarios.values():
        _fields(row, "scenario_id disease_id gene_id genotype_id allele_ids evidence_ids mechanism_inputs limitations")
        _require(row["scenario_id"] == scenario_identity(row) and row["genotype_id"]
                 and row["gene_id"] in scope["genes"] and row["disease_id"]
                 and row["allele_ids"] and len(set(row["allele_ids"])) == len(row["allele_ids"])
                 and set(row["allele_ids"]) <= set(variants)
                 and all(row["gene_id"] in variants[a]["gene_ids"]
                         and variants[a]["assessment"] != "excluded" for a in row["allele_ids"])
                 and row["evidence_ids"] and set(row["evidence_ids"]) <= set(evidence),
                 "genotype scenario is unbound or outside selected scope")
        facts = row["mechanism_inputs"]
        readiness = assess_mechanism_readiness(facts)
        _require(readiness.scenario_id == row["scenario_id"]
                 and {a.allele_id for a in readiness.allele_evidence} == set(row["allele_ids"])
                 and all(set(a.evidence_ids) <= set(evidence) for a in readiness.allele_evidence),
                 "mechanism facts disagree with the source genotype")
        # The existing readiness kernel owns all scientific thresholds. Unknown
        # direction/phase/partner is admissible and remains non-passing there.
    phenotype = parts["phenotype"]
    _fields(phenotype, "state source_sha256 observations family_context limitations")
    _require(phenotype["state"] in {"assessed", "not_assessed"}
             and isinstance(phenotype["observations"], list) and isinstance(phenotype["family_context"], list),
             "optional phenotype state is invalid")
    if phenotype["state"] == "not_assessed":
        _require(phenotype["source_sha256"] is None and not phenotype["observations"]
                 and not phenotype["family_context"], "absent HPO cannot invent patient observations")
    else:
        _require(_digest(phenotype["source_sha256"]), "supplied HPO source binding is missing")


def write_diagnosis_package(output, *, components, classification, evidence_origin, runtime, exclusion_source,
                            builder_code_sha256, policy_sha256, controlled_store=None, deferred_rows=()):
    """Materialize the builder's small source-bound records, without ready flags.

    The VCF builder owns source authentication and scientific construction. This
    serialization boundary neither reads a VCF nor creates a review on its behalf.
    """
    validate_components(components, classification=classification, evidence_origin=evidence_origin)
    artifacts = {role: {"path": role + ".json", "sha256": sha256(canonical(components[role]))}
                 for role in ROLES}
    manifest = dict(schema=SCHEMA, input_origin="diagnosis_guided_vcf", classification=classification,
                    evidence_origin=evidence_origin, artifacts=artifacts,
                    builder=dict(code_sha256=builder_code_sha256, policy_sha256=policy_sha256,
                                 schema_sha256=sha256(SCHEMA_PATH.read_bytes()), runtime=runtime))
    manifest["package_id"] = "diagnosis-package:" + sha256(canonical(manifest))[:24]
    validate_json_schema(manifest, json.loads(SCHEMA_PATH.read_bytes()))
    from genosapiens.track2_storage import output_storage
    with output_storage(output, classification=classification, controlled_store=controlled_store,
                        kind='track2-diagnosis-package') as target:
        target.mkdir(parents=True, exist_ok=False)
        from genosapiens.track2_exclusions import copy_exclusions, DIRECTORY
        copy_exclusions(exclusion_source,target/DIRECTORY,components['selection']['coverage']['filter_exclusions'],
            selected_record_ids=(r['source_record_id'] for r in components['selection']['selected_alleles']))
        from genosapiens.track2_deferred import write_inventory, DIRECTORY as DEFERRED_DIRECTORY
        write_inventory(target/DEFERRED_DIRECTORY, deferred_rows, components['selection'])
        for role in ROLES:
            (target / artifacts[role]["path"]).write_bytes(canonical(components[role]))
        (target / "package.json").write_bytes(canonical(manifest))
        (target / 'acceptance.json').write_bytes(canonical(input_acceptance(manifest,components,
            datetime.now(timezone.utc).isoformat())))
    return manifest


@dataclass(frozen=True)
class DiagnosisPackage:
    data: dict
    components: dict
    input_sha256: str
    source_hashes: dict
    review: dict | None
    review_sha256: str | None
    input_report_ready: bool
    automatic_acceptance: dict | None = None
    automatic_acceptance_sha256: str | None = None
    root: Path | None = None

    @property
    def human_review_ready(self):
        return bool(self.review and self.review['decision']=='accepted'
            and self.review['reviewer_type'] in {'human','synthetic_fixture'})

    def admitted_origins(self):
        _require(self.input_report_ready, "diagnosis package requires source-bound acceptance")
        return self.report_origins()

    def report_origins(self):
        """Identify constructed scenarios, retaining absent review as null."""
        return tuple(AnalysisOrigin("diagnosis_guided_vcf", row["scenario_id"], self.input_sha256,
                                    self.review_sha256) for row in self.components["scenarios"])


def render_diagnosis_dossier(bundle, *, exclusion_path='package/filter-exclusions/records.jsonl',
                             deferred_path='package/deferred-alleles/records.jsonl'):
    """Minimal package view only; no raw input or public acquisition access."""
    from genosapiens.track2_phenotype import render_phenotype_context
    p=bundle.components
    lines=['# Diagnosis-guided input dossier','',
        'Input origin: diagnosis_guided_vcf. Classification: '+bundle.data['classification']+'.',
        'Evidence origin: '+bundle.data['evidence_origin']+'.',
        'Package validity and review do not establish a causal genotype or treatment.','',
        'Diagnosis: '+p['scope']['resolved_diagnosis']+'.',
        f"Genes in the frozen scope: {len(p['scope']['genes'])}. Gene and regulatory intervals: {len(p['scope']['intervals'])}.",
        f"Selected sequence alleles: {len(p['selection']['selected_alleles'])}. Mechanism alternatives: {len(p['scenarios'])}.",
        'Structural interpretation and outside-scope partners remain unassessed.','',
        '| Scenario | Gene | Genotype completeness | Phase | Direction |','|---|---|---|---|---|']
    for row in p['scenarios']:
        facts=row['mechanism_inputs']
        lines.append('| '+ ' | '.join((row['scenario_id'],row['gene_id'],facts['genotype_completeness'],
            facts['phase'],facts['direction']))+' |')
    if not p['scenarios']:lines+=['','The builder constructed no genotype/mechanism hypothesis within the declared scope.']
    coverage=p['selection']['coverage']
    from genosapiens.track2_selection_report import render_selection
    lines+=['',render_selection(p['selection'],deferred_path=deferred_path)]
    lines+=['','## Automatic input acceptance','',
        'The machine checked input integrity and source bindings. This decision does not establish a causal genotype.',
        'Policy: '+bundle.automatic_acceptance['policy']+'.',
        'Outcome: '+bundle.automatic_acceptance['outcome']+'.',
        'Input review authority (exact source value): `'+(bundle.review['reviewer_type'] if bundle.review is not None else 'none')+'`.','',
        f"Exact FILTER exclusions: {coverage['filter_excluded_record_count']} source records and "
        f"{coverage['filter_excluded_alternate_count']} alternate entries.",
        'Counts cover every excluded source occurrence, including duplicate lines.',
        'The [complete exclusion artifact]('+exclusion_path+') retains the source identities and reasons.',
        'The record retains ExcessHet without a new threshold.','',*construction_details(p)]
    from genosapiens.track2_vcf_evidence import render as render_vcf_evidence
    lines+=['',*render_vcf_evidence(p)]
    from genosapiens.track2_public_matching import render_matching
    from genosapiens.track2_structural import render_structural
    return ('\n'.join(lines)+'\n\n'+render_structural(p['structural'])+'\n'
        +render_phenotype_context(p['phenotype'])+'\n'+render_matching(p))


def construction_details(parts):
    """Render only source-bound copy, quality, and transcript records."""
    def cell(value):
        return str(value if value is not None else 'unresolved').replace('|','\\|').replace('\n',' ')
    evidence={e['evidence_id']:e for e in parts['evidence']}
    lines=['## Copy requirements and sample quality','',
        'Observed copies and required affected copies have different meanings.',
        'A diploid expectation does not supply a missing allele.','',
        '| Scenario | Allele | Observed ploidy | ALT copies | Expected ploidy | Required copies | Derivation | Quality |',
        '|---|---|---|---|---|---|---|---|']
    transcripts=[];measurements=[]
    for scenario in parts['scenarios']:
        proofs=[json.loads(evidence[k]['assertion']) for k in scenario['evidence_ids'] if k.startswith('scoped-model:')]
        if not proofs: continue
        proof=proofs[0];genotype=proof['genotype'];requirements=genotype['hypothesis']['context']['requirements']
        mapping=proof['normalized_to_package_allele_ids']
        alleles={r['allele_id']:r for k in scenario['evidence_ids'] if k.startswith('scoped-allele:')
                 for r in [json.loads(evidence[k]['assertion'])]}
        for member in genotype['model']['members']:
            key=member['allele']['allele_id'];prop=member['genotype']
            qualities=[c['source']['source_fields'].get('automatic_quality',{}) for c in alleles[key]['calls'] or ()]
            from genosapiens.genotype_model import _LITERAL_RULES, _INHERITANCE_RULES
            rules=[_LITERAL_RULES.get(r['allelic_requirement']) or _INHERITANCE_RULES.get(r['inheritance'])
                   for r in requirements]
            expected=sorted({n for rule in rules if rule and rule[1]==member['locus']['kind'] for n in rule[2]}) or 'unresolved'
            values=(scenario['scenario_id'],mapping[key],prop['reported_ploidies'],
                [prop['minimum_alt_copies'],prop['maximum_alt_copies']],expected,
                [r['required_allele_count'] for r in requirements],
                [d for r in requirements for d in r['diagnostics'] if d.startswith('affected_copies:')],
                [q.get('state','unresolved') for q in qualities])
            lines.append('| '+' | '.join(cell(v) for v in values)+' |')
            for q in qualities:
                measurements.extend(['',f"Allele `{cell(mapping[key])}`: GT `{cell(q.get('genotype'))}`, "
                    f"GQ `{cell(q.get('genotype_quality'))}`, DP `{cell(q.get('depth'))}`, "
                    f"balance `{cell(q.get('allele_balance'))}`.",
                    'Quality reasons (exact identifiers): '+', '.join(q.get('reasons',['unresolved']))+'.'])
        for key in scenario['evidence_ids']:
            if not key.startswith('scoped-prediction:'): continue
            p=json.loads(evidence[key]['assertion'])
            for tx in p['transcripts']:
                transcripts.append((scenario['scenario_id'],mapping[p['allele_id']],tx['transcript_id'],p['selection_reason'],
                    tx['source']['source_release'],tx['consequences'],tx['sift'],tx['polyphen'],tx['direction'],
                    tx['predictor_disagreement']))
    lines+=[*measurements,'','## Transcript predictions','',
        'These predictions do not establish protein absence, residual activity, pathogenicity, or treatment response.',
        'The structured source records retain lower-priority transcripts and missing mappings.','',
        '| Scenario | Allele | Transcript | Selection | Release | Consequence | SIFT | PolyPhen | Predicted direction | Predictor disagreement |',
        '|---|---|---|---|---|---|---|---|---|---|',
        *('| '+' | '.join(cell(v) for v in row)+' |' for row in transcripts)]
    return lines


def load_diagnosis_package(root, *, require_review=True, allow_synthetic=False,
                           expected_input_sha256=None, expected_review_sha256=None):
    root = Path(root)
    raw = _source_file(root, "package.json").read_bytes()
    data = json.loads(raw)
    validate_json_schema(data, json.loads(SCHEMA_PATH.read_bytes()))
    _require(raw == canonical(data), "diagnosis package serialization is not canonical")
    _require(data["package_id"] == "diagnosis-package:" + sha256(canonical({
        k: v for k, v in data.items() if k != "package_id"}))[:24], "diagnosis package identity changed")
    _require(data["builder"]["schema_sha256"] == sha256(SCHEMA_PATH.read_bytes()),
             "diagnosis package uses an obsolete contract")
    _require(allow_synthetic or data["evidence_origin"] != "independently_synthetic",
             "synthetic package requires explicit synthetic input selection")
    _require((data["classification"] == "controlled") == (data["evidence_origin"] == "controlled_case"),
             "case source classification cannot be substituted")
    parts, hashes = {}, {}
    for role, record in data["artifacts"].items():
        _require(record["path"] == role + ".json", "unexpected diagnosis companion path")
        content = _source_file(root, record["path"]).read_bytes()
        _require(sha256(content) == record["sha256"], "diagnosis source companion checksum mismatch")
        parts[role] = json.loads(content)
        _require(canonical(parts[role]) == content, "diagnosis companion serialization is not canonical")
        hashes[role] = record["sha256"]
    validate_components(parts, classification=data["classification"], evidence_origin=data["evidence_origin"])
    from genosapiens.track2_exclusions import check_exclusions, DIRECTORY
    check_exclusions(root/DIRECTORY,parts['selection']['coverage']['filter_exclusions'],
        selected_record_ids=(r['source_record_id'] for r in parts['selection']['selected_alleles']))
    from genosapiens.track2_deferred import check_inventory, DIRECTORY as DEFERRED_DIRECTORY
    check_inventory(root/DEFERRED_DIRECTORY, parts['selection'])
    input_hash = sha256(raw)
    acceptance_raw=_source_file(root,'acceptance.json').read_bytes()
    acceptance=json.loads(acceptance_raw)
    _require(acceptance_raw==canonical(acceptance) and acceptance==input_acceptance(data,parts,acceptance['timestamp']),
             'Automatic input acceptance differs from its validated sources.')
    _require(expected_input_sha256 is None or expected_input_sha256 == input_hash,
             "frozen diagnosis package changed")
    review, review_hash, ready = None, None, False
    if (root / "review.json").exists():
        review_raw = _source_file(root, "review.json").read_bytes()
        review, review_hash = json.loads(review_raw), sha256(review_raw)
        _fields(review, "schema reviewer reviewer_type independence reviewed_at input_sha256 "
                        "source_hashes scenario_ids decision rationale limitations")
        _require(review["schema"] == "genosapiens.track2-diagnosis-review/v1"
                 and review["input_sha256"] == input_hash and review["source_hashes"] == hashes
                 and review["scenario_ids"] == [s["scenario_id"] for s in parts["scenarios"]]
                 and review["reviewer"] and review["rationale"] and review["limitations"]
                 and review["reviewer_type"] in {"human", "ai_assisted_extraction", "synthetic_fixture"}
                 and review["independence"] in {"independent", "self_review_disclosed"}
                 and review["decision"] in {"accepted", "pending_human_review", "rejected"}
                 and datetime.fromisoformat(review["reviewed_at"].replace("Z", "+00:00")).tzinfo is not None,
                 "diagnosis source review is missing, stale or invalid")
        synthetic_review = review["reviewer_type"] == "synthetic_fixture"
        _require(not synthetic_review or allow_synthetic and data["evidence_origin"] == "independently_synthetic",
                 "synthetic review cannot admit actual case evidence")
        ready = review["decision"] == "accepted" and (review["reviewer_type"] == "human" or synthetic_review)
        _require(review["decision"] != "accepted" or ready, "AI review cannot impersonate human acceptance")
    _require(expected_review_sha256 is None or expected_review_sha256 == review_hash,
             "frozen diagnosis review changed")
    _require(not require_review or ready, "diagnosis package requires actual source-bound review")
    automatic_ready=acceptance['outcome']=='accepted' and (review is None or review['decision']!='rejected')
    return DiagnosisPackage(data, parts, input_hash, hashes, review, review_hash,ready or automatic_ready,
                            acceptance,sha256(acceptance_raw),root)
