"""Source-specific G0 admission into one scientific identity contract.

An AnalysisOrigin identifies records; only the corresponding source loader can
admit them. In particular, constructing an identity never supplies a G0 pass.
"""

from dataclasses import dataclass
import re

SCENARIO_PREFIXES = {
    "published_patient_genotype": "published-genotype:",
    "track1_case_handoff": "track1-hypothesis:",
    "diagnosis_guided_vcf": "diagnosis-genotype:",
}


class Track2InputError(ValueError):
    """Source admission failed before scientific synthesis."""


@dataclass(frozen=True, slots=True)
class AnalysisOrigin:
    input_origin: str
    scenario_id: str
    source_sha256: str
    review_sha256: str | None

    def __post_init__(self) -> None:
        prefix = SCENARIO_PREFIXES.get(self.input_origin)
        if prefix is None or not isinstance(self.scenario_id, str) or not self.scenario_id.startswith(prefix):
            raise Track2InputError("scenario identity disagrees with input origin")
        if not self.scenario_id[len(prefix):]:
            raise Track2InputError("scenario identity is required")
        hashes = (self.source_sha256, self.review_sha256)
        if self.input_origin == "diagnosis_guided_vcf" and self.review_sha256 is None:
            hashes = (self.source_sha256,)
        if any(not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None
               for value in hashes):
            raise Track2InputError("source and review SHA-256 bindings are required")


def load_analysis_input(root, *, input_origin, require_review=True, allow_synthetic=False,
                        expected_input_sha256=None, expected_review_sha256=None):
    """G0: validate the selected authority, never infer origin from filenames."""
    options = dict(require_review=require_review, allow_synthetic=allow_synthetic,
                   expected_input_sha256=expected_input_sha256,
                   expected_review_sha256=expected_review_sha256)
    if input_origin == "published_patient_genotype":
        from genosapiens.track2_published_genotype import load_published_genotype
        return load_published_genotype(root, **options)
    if input_origin == "diagnosis_guided_vcf":
        from genosapiens.track2_diagnosis_package import load_diagnosis_package
        return load_diagnosis_package(root, **options)
    if input_origin == "track1_case_handoff":
        raise Track2InputError("Track 1 minimal-view and immutable-selection dependency is unmet")
    raise Track2InputError("unknown scientific input origin")


def admitted_origins(bundle):
    return (bundle.admitted_origin(),) if bundle.data['input_origin']=='published_patient_genotype' else bundle.admitted_origins()


def assessed_report_origins(bundle, *, draft=False, assessment=None):
    """Draft diagnosis reports identify sources without claiming reviewed admission."""
    if assessment is not None:
        from genosapiens.track2_mechanism_assessment import validate_assessment
        validate_assessment(bundle, assessment)
        if not draft and not bundle.input_report_ready:
            raise Track2InputError('diagnosis package requires actual source-bound review')
        return tuple(AnalysisOrigin('diagnosis_guided_vcf', row['scenario_id'],
            bundle.input_sha256, bundle.review_sha256) for row in assessment['scenarios'])
    return report_origins(bundle, draft=draft)


def input_evidence(bundle):
    return bundle.data['evidence'] if bundle.data['input_origin']=='published_patient_genotype' else bundle.components['evidence']


def assessed_mechanism_inputs(bundle, supplied, *, draft=False, assessment=None):
    if assessment is None:
        return source_mechanism_inputs(bundle, supplied, draft=draft)
    origins = assessed_report_origins(bundle, draft=draft, assessment=assessment)
    expected = {s['scenario_id']:s['mechanism_inputs'] for s in assessment['scenarios']}
    if supplied != expected:
        raise Track2InputError('mechanism assessments differ from the reviewed diagnosis package')
    return origins


def report_origins(bundle, *, draft=False):
    """Draft diagnosis reports identify sources without claiming reviewed admission."""
    if draft and bundle.data['input_origin'] == 'diagnosis_guided_vcf':
        return bundle.report_origins()
    return admitted_origins(bundle)


def source_mechanism_inputs(bundle, supplied, *, draft=False):
    """Apply each origin's actual genotype authority before shared synthesis."""
    origins=report_origins(bundle,draft=draft)
    if set(supplied)!={o.scenario_id for o in origins}:
        raise Track2InputError('mechanism scenarios differ from the admitted input')
    if bundle.data['input_origin']=='published_patient_genotype':
        facts=supplied[origins[0].scenario_id]
        if 'prediction_acceptance' in facts or facts.get('joint_evidence')=='prediction_accepted':
            raise Track2InputError('Automatic prediction acceptance requires diagnosis construction evidence.')
        phase='not_applicable' if bundle.data['genotype']['zygosity']=='homozygous' else 'confirmed_trans'
        if (facts['scenario_id']!=origins[0].scenario_id or facts['genotype_completeness']!='complete'
            or facts['phase']!=phase or {a['allele_id'] for a in facts['allele_evidence']}
            !={a['allele_id'] for a in bundle.data['genotype']['alleles']}):
            raise Track2InputError('mechanism genotype disagrees with the reviewed source')
    else:
        expected={s['scenario_id']:s['mechanism_inputs'] for s in bundle.components['scenarios']}
        if supplied!=expected:
            raise Track2InputError('mechanism assessments differ from the reviewed diagnosis package')
    return origins
