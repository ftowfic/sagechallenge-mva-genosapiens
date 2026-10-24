"""Bounded projection of authenticated LongTR source records into native views."""

from dataclasses import dataclass, replace
import hashlib
from itertools import groupby
import json

from genosapiens import longtr_admission as longtr
from genosapiens.evidence_model import SourcedRecord
from genosapiens.variant_normalization import QuarantineReason
from genosapiens.vcf_qc import _structured_attributes, _parse_info
from genosapiens.wgs_admission import (
    LedgerAdmissionInput, SourceTechnicalObligation, WgsAdmissionError,
    RULE_IDS, _occurrence, _source_record, _technical,
)
from genosapiens.wgs_ledger import canonical, iter_ledger
from genosapiens.wgs_source_routing import (
    WgsSourceRoutingError,
    classify_source_occurrence,
)


@dataclass(frozen=True, slots=True)
class _Record:
    source: SourcedRecord
    qc: object


@dataclass(frozen=True, slots=True)
class _Mapping:
    normalized_allele_id: str


def _header(header):
    lines = tuple(line.rstrip("\r\n") for line in header.raw_lines)
    definitions = {}
    try:
        for line in lines:
            if not line.startswith(("##INFO=<", "##FORMAT=<")):
                continue
            namespace = line[2:].split("=", 1)[0]
            attributes = _structured_attributes(line.split("=<", 1)[1][:-1])
            key = (namespace, attributes["ID"])
            if not line.endswith(">") or key in definitions:
                raise ValueError
            definitions[key] = (attributes["Number"], attributes["Type"])
        longtr.validate_header(lines, definitions)
    except (ValueError, TypeError, KeyError, IndexError):
        raise WgsAdmissionError("bounded LongTR source header is invalid") from None


def _retained_projection_rows(connection, source_manifest, working_manifest, execution_config):
    cursor = connection.execute("SELECT o.record_number,o.allele_number,o.source_allele_id,"
        "p.payload_json FROM occurrences o JOIN source_payloads p USING(record_number) "
        "ORDER BY o.record_number,o.allele_number")
    observed = 0
    for record_number, group in groupby(cursor, key=lambda row: row[0]):
        rows = tuple(group)
        observed += len(rows)
        if (not rows or type(record_number) is not int or record_number < 1
                or record_number > source_manifest.summary.total_records
                or any(row[1] != ordinal or row[3] != rows[0][3]
                       for ordinal, row in enumerate(rows, 1))):
            raise WgsAdmissionError("bounded LongTR retained locus is incomplete")
        payload_json = rows[0][3]
        retained = tuple(row[2] for row in rows)
        try:
            expected = tuple(row["source_allele_id"] for row in
                json.loads(payload_json)["normalization"]["alleles"])
        except (ValueError, TypeError, KeyError):
            raise WgsAdmissionError("bounded LongTR retained locus is invalid") from None
        if (retained != expected or len(payload_json.encode("ascii")) * len(retained)
                > execution_config.annotation.max_response_row_bytes):
            raise WgsAdmissionError("bounded LongTR retained locus is incomplete")
        yield retained, payload_json
    if observed != working_manifest.summary.retained_occurrences:
        raise WgsAdmissionError("bounded LongTR retained coverage is incomplete")


def iter_longtr_admission_inputs(connection, *, source_root, source_manifest,
                                working_manifest, routing, header, execution_config,
                                from_retained_projection=False):
    """Project one complete source record at a time, preserving every retained ALT.

    The existing LongTR kernel supplies interpretation and parent/member identity.
    The ledger and working-set replay supply original-input and retention authority.
    A locus with unsupported member correspondence keeps its parent and explicit
    unresolved source records; a reference locus remains a native call.
    """
    _header(header)
    contract = routing.contract
    profile = {
        "schema": "genosapiens.supplied-vcf-format",
        "input_sha256": source_manifest.source.input_sha256,
        "reference_role": contract.reference_role.value,
        "reference_identity_id": contract.reference_identity_id,
        "reference_manifest_sha256": contract.reference_manifest_sha256,
        "mitochondrial_contigs": list(contract.mitochondrial_contigs),
        "repeat_identity": longtr.FORMAT,
    }
    profile_sha256 = hashlib.sha256(canonical(profile)).hexdigest()
    required = tuple(rule for rule in RULE_IDS if rule in (RULE_IDS[2], RULE_IDS[3], RULE_IDS[4]))
    ordinal = 0
    def source_ledger_rows():
        for row in iter_ledger(source_root, source_manifest):
            retained = tuple(source_id for source_id in row.source_allele_ids if connection.execute(
                "SELECT 1 FROM occurrences WHERE source_allele_id=?", (source_id,)).fetchone())
            if retained:
                if retained != row.source_allele_ids:
                    raise WgsAdmissionError("bounded LongTR retained locus is incomplete")
                yield retained, row.payload_json

    source_rows = (_retained_projection_rows(connection, source_manifest,
                                           working_manifest, execution_config)
                   if from_retained_projection else source_ledger_rows())
    for retained, payload_json in source_rows:
        if not retained:
            continue
        if len(payload_json.encode("ascii")) * len(retained) > execution_config.annotation.max_response_row_bytes:
            raise WgsAdmissionError("bounded LongTR record exceeds capacity")
        values = tuple(_occurrence(payload_json, source_allele_id=source_id,
            source_manifest=source_manifest, contract=contract) for source_id in retained)
        try:
            routes = tuple(classify_source_occurrence(
                value, contract=contract, renderer=routing.renderer,
            ) for value in values)
        except WgsSourceRoutingError:
            raise WgsAdmissionError("bounded LongTR source routing is invalid") from None
        if any(
            route.preparation_disposition != "native_class"
            or route.normalized_allele_id is not None
            or "repeat" not in route.native_class_tags
            for route in routes
        ):
            raise WgsAdmissionError("bounded LongTR retained route is inconsistent")
        first = values[0]
        qc = replace(first.record, sample_calls=(first.sample,))
        blocked = tuple(dict.fromkeys(value.normalization.quarantine_reason.value for value in values
            if value.normalization.quarantine_reason in (
                QuarantineReason.REFERENCE_MISMATCH, QuarantineReason.REFERENCE_UNRESOLVABLE)))
        views = ()
        reasons = {source_id: blocked for source_id in retained}
        try:
            if blocked or qc.raw_line is None:
                raise ValueError
            columns = tuple(qc.raw_line.rstrip("\r\n").split("\t"))
            source = SourcedRecord("supplied-vcf", "VCFv4.1", qc.source_record_id,
                source_manifest.source.input_sha256, None, {
                    **_source_record(first, _technical(first)).source_fields,
                    "source_payload_sha256": first.payload_sha256,
                    "columns": columns,
                    "info": _parse_info(columns[7]),
                    "format_keys": tuple(columns[8].split(":")),
                    "sample_values": dict(zip(header.sample_names, columns[9:], strict=True)),
                })
            normalized = {value.normalization.normalized_allele.normalized_allele_id:
                value.normalization.normalized_allele for value in values
                if value.normalization.normalized_allele is not None}
            mappings = {value.source.source_allele_id:
                _Mapping(value.normalization.normalized_allele.normalized_allele_id)
                for value in values if value.normalization.normalized_allele is not None}
            views, links, reasons = longtr.build_views(_Record(source, qc),
                tuple(value.source for value in values), profile=profile,
                profile_sha256=profile_sha256, role=contract.reference_role,
                reference_id=contract.reference_identity_id,
                reference_sha256=contract.reference_manifest_sha256,
                normalized=normalized, mappings=mappings)
            retained_set = set(retained)
            view_ids = {allele.allele_id for allele in views}
            if (
                set(links) != retained_set
                or set(reasons) != retained_set
                or len(view_ids) != len(views)
                or any(
                    set(identity_ids) - view_ids
                    for identity_ids in links.values()
                )
            ):
                raise ValueError
            for allele in views:
                expected_sources = tuple(
                    source_id for source_id in retained
                    if allele.allele_id in links[source_id]
                )
                actual_sources = tuple(
                    call.source.source_fields.get("source_allele_id")
                    for call in allele.calls or ()
                )
                if not expected_sources or actual_sources != expected_sources:
                    raise ValueError
        except (ValueError, TypeError, KeyError, IndexError, OverflowError):
            views = ()
            reasons = {source_id: blocked or ("longtr-source-interpretation-unresolved",)
                       for source_id in retained}
        for allele in views:
            obligations = tuple(SourceTechnicalObligation(
                call.source.source_fields["source_allele_id"], call.sample.sample_name,
                tuple(call.source.source_fields["technical_flags"]),
                call.source.source_fields["read_validation_required"], bool(call.filters),
            ) for call in allele.calls)
            ordinal += 1
            result = LedgerAdmissionInput(ordinal, allele.allele_id,
                tuple(value.source_allele_id for value in obligations), "repeat_catalog", allele,
                None, required, obligations, source_manifest.sha256,
                working_manifest.sha256, contract.sha256, "pending", ())
            result.as_bytes()
            yield result
        for value in values:
            limits = tuple(dict.fromkeys(reasons[value.source.source_allele_id]))
            if not limits:
                continue
            ordinal += 1
            result = LedgerAdmissionInput(ordinal, value.source.source_allele_id,
                (value.source.source_allele_id,), "unresolved", None, value.source,
                required, (_technical(value),), source_manifest.sha256,
                working_manifest.sha256, contract.sha256, "unresolved", limits)
            result.as_bytes()
            yield result
