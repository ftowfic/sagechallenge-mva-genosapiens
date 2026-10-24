"""Explicit case parents and once-per-case scientific resource ownership."""
from contextlib import ExitStack
import hashlib

from genosapiens.evidence_fusion import FusionPolicy, policy_as_json
from genosapiens.evidence_model import _ids, _sha, _typed
from genosapiens.evidence_serialization import _canonical
from genosapiens.genotype_evaluation import prepare_genotype_policies
from genosapiens.genotype_hypotheses import HypothesisPolicy, hypothesis_policy_bytes
from genosapiens.track1_sensitivity import DisclosedDiagnosis, PriorPolicy
from genosapiens.track1_submission import SubmissionReferenceBinding
from genosapiens.wgs_evidence_replay import EvidenceReplay
from genosapiens.wgs_evidence_store import EvidenceStoreError
from genosapiens.wgs_gene_evidence import VerifiedCaseGeneEvidence
from genosapiens.wgs_phenotype_evidence import VerifiedCasePhenotypeEvidence, _hash
from genosapiens.wgs_population_evidence import VerifiedPopulationPolicy

MAX_CASE_BINDING_BYTES = 8 * 1024**2


def _require_policy_bindings(hypothesis_policy, fusion_policy, *,
        admission_policy_sha256, hypothesis_policy_sha256, fusion_policy_sha256,
        outside_context_partner_uncertainty):
    if (hypothesis_policy.maximum_pairs_per_context is not None
            or type(outside_context_partner_uncertainty) is not bool
            or hashlib.sha256(hypothesis_policy_bytes(hypothesis_policy,
                admission_policy_sha256=admission_policy_sha256)).hexdigest() != hypothesis_policy_sha256
            or hashlib.sha256(policy_as_json(fusion_policy)).hexdigest() != fusion_policy_sha256):
        raise EvidenceStoreError("scientific case policy binding changed")


def require_saved_case_policy_bindings(value):
    """Check small saved policies before restoring large case resources."""
    _require_policy_bindings(
        _typed(value['hypothesis_policy'], HypothesisPolicy, decode=True),
        _typed(value['fusion_policy'], FusionPolicy, decode=True),
        **{key: value[key] for key in ('admission_policy_sha256', 'hypothesis_policy_sha256',
            'fusion_policy_sha256', 'outside_context_partner_uncertainty')})


def bound_reference(reader, value, kind, *, extra=()):
    if (type(value) is not dict or set(value) != {"name","kind","sha256",*extra}
            or value["kind"] != kind):
        raise EvidenceStoreError("invalid case dependency reference")
    reference = reader.reference(value["name"],kind=kind)
    if reference.sha256 != value["sha256"]:
        raise EvidenceStoreError("case dependency digest changed")
    return reference


def component_reference(owner, kind):
    return {"name":owner.name,"kind":kind,"sha256":owner.sha256}


class ScientificCaseResources:
    def __init__(self, writer, *, genes, phenotypes, population, hypothesis_policy,
                 hypothesis_policy_sha256, admission_policy_sha256, fusion_policy,
                 fusion_policy_sha256, prior_policy, diagnosis, sample_name,
                 source_manifest_sha256, reference, parent_bindings,
                 outside_context_partner_uncertainty=True, analysis_eligibility=None, analysis_profile=None, name):
        self.writer,self.connection,self.name = writer,writer.connection,name
        self.genes,self.phenotypes,self.population = genes,phenotypes,population
        try:
            genes.require(writer,genes.prepared.evidence)
            population._require_open()
            if any(o.connection is not self.connection for o in (phenotypes,population)):
                raise EvidenceStoreError("foreign scientific case resources")
            _ids((sample_name,))
            _sha(source_manifest_sha256)
            _sha(admission_policy_sha256)
            if (type(parent_bindings) is not dict or not parent_bindings
                    or any(type(k) is not str or not k for k in parent_bindings)):
                raise EvidenceStoreError("invalid scientific case parents")
            for digest in parent_bindings.values():
                _sha(digest)
            self.genotype_policies = prepare_genotype_policies(
                _typed(hypothesis_policy,HypothesisPolicy), ('loci','requirements','clinical'))
            self.hypothesis_policy = self.genotype_policies.parent
            self.fusion_policy = _typed(fusion_policy,FusionPolicy)
            self.prior_policy = _typed(prior_policy,PriorPolicy)
            self.diagnosis = _typed(diagnosis,DisclosedDiagnosis | None)
            self.reference_binding = _typed(reference,SubmissionReferenceBinding)
            _require_policy_bindings(self.hypothesis_policy, self.fusion_policy,
                admission_policy_sha256=admission_policy_sha256,
                hypothesis_policy_sha256=hypothesis_policy_sha256, fusion_policy_sha256=fusion_policy_sha256,
                outside_context_partner_uncertainty=outside_context_partner_uncertainty)
            self.hypothesis_policy_sha256 = hypothesis_policy_sha256
            self.fusion_policy_sha256 = fusion_policy_sha256
            self.admission_policy_sha256 = admission_policy_sha256
            self.sample_name,self.source_manifest_sha256 = sample_name,source_manifest_sha256
            self.observation_ids = tuple(o.observation_id for o in phenotypes.view.observations)
            self.outside_context_partner_uncertainty = outside_context_partner_uncertainty
            self.parent_bindings = dict(parent_bindings)
            self.analysis_eligibility=analysis_eligibility
            component=None
            if analysis_eligibility is not None:
                from genosapiens.wgs_analysis_eligibility import eligibility_component_parents
                expected=eligibility_component_parents(analysis_eligibility)
                if (any(self.parent_bindings.get(k)!=v for k,v in expected.items())
                        or any(dict(analysis_eligibility['input_dispositions'])[key] for key in ('defer_refcall','exclude_common'))
                        and not outside_context_partner_uncertainty):
                    raise EvidenceStoreError('scientific eligibility parent or partner uncertainty changed')
                digest=writer.add('case-analysis-eligibility','analysis-eligibility-component-v1',analysis_eligibility)
                component=dict(name='case-analysis-eligibility',kind='analysis-eligibility-component-v1',sha256=digest)
            elif any(k.startswith('analysis_eligibility_') for k in parent_bindings):
                raise EvidenceStoreError('missing scientific eligibility component')
            self.analysis_profile=analysis_profile
            profile_component=None
            if analysis_profile is not None:
                from genosapiens.wgs_profile_projection import profile_component_parents
                expected=profile_component_parents(analysis_profile)
                if (analysis_eligibility is None or any(self.parent_bindings.get(k)!=v for k,v in expected.items())
                        or not outside_context_partner_uncertainty):
                    raise EvidenceStoreError('scientific profile parent or partner uncertainty changed')
                digest=writer.add('case-analysis-profile','analysis-profile-component-v1',analysis_profile)
                profile_component=dict(name='case-analysis-profile',kind='analysis-profile-component-v1',sha256=digest)
            elif any(k.startswith('analysis_profile_') for k in parent_bindings):
                raise EvidenceStoreError('missing scientific profile component')
            from genosapiens.wgs_analysis_eligibility import scientific_assessment_domain
            self.assessment_domain=scientific_assessment_domain(analysis_eligibility,self.parent_bindings,analysis_profile)
            self.assessment_domain_sha256=_hash(self.assessment_domain)
            view_sha256, profile_sha256 = phenotypes.logical_hashes(writer)
            self.header = {
                "schema":"genosapiens.wgs-scientific-case-resources","version":4,
                "assessment_domain":self.assessment_domain,
                "analysis_eligibility":component,
                "analysis_profile":profile_component,
                "analysis_domain_status":'not_bound' if analysis_eligibility is None else analysis_eligibility['analysis_domain_status'],
                "rescue_assessment_status":'not_performed' if analysis_eligibility is None else analysis_eligibility['rescue_assessment_status'],
                "genes":genes.reference,"phenotypes":phenotypes.reference,"population":population.reference,
                "run_manifest_sha256":phenotypes.run_manifest_sha256,
                "view_sha256":view_sha256,"profile_sha256":profile_sha256,
                "ic_sha256":phenotypes.ic_sha256,"hypothesis_policy":self.hypothesis_policy,
                "hypothesis_policy_sha256":hypothesis_policy_sha256,"admission_policy_sha256":admission_policy_sha256,
                "fusion_policy":self.fusion_policy,"fusion_policy_sha256":fusion_policy_sha256,
                "prior_policy":self.prior_policy,"diagnosis":self.diagnosis,"sample_name":sample_name,
                "source_manifest_sha256":source_manifest_sha256,"reference":self.reference_binding,
                "parent_bindings":self.parent_bindings,
                "outside_context_partner_uncertainty":outside_context_partner_uncertainty,
            }
            if len(_canonical(self.header)) > MAX_CASE_BINDING_BYTES:
                raise EvidenceStoreError("scientific case bindings exceed capacity",capacity=True)
            self.sha256 = writer.add(name,"scientific-case-resources-v1",self.header)
        except BaseException:
            writer.failed = True
            raise

    @property
    def reference(self):
        return component_reference(self,"scientific-case-resources-v1")

    @property
    def gene_logical_sha256(self):
        return self.header['genes']['logical_sha256']

    @property
    def profile_evidence_id(self):
        return None if self.phenotypes.profile is None else self.phenotypes.profile.evidence_id

    @property
    def population_reference(self):
        return self.header['population']

    def require(self, writer):
        if self.connection is None or writer.connection is not self.connection:
            raise EvidenceStoreError("scientific case resource owner is closed or foreign")
        if type(writer) is EvidenceReplay:
            writer.reader._check_files()
        if getattr(writer,"failed",False):
            raise EvidenceStoreError("scientific case resource operation has failed")
        self.genes.require(writer,self.genes.prepared.evidence)
        self.population._require_open()


class CompletedScientificCaseResources(ScientificCaseResources):
    """Bound resource metadata; full data is loaded only for selected fallback."""
    def __init__(self, reader, reference, *, expected_parent_bindings, eligibility_root=None,
                 profile_root=None, progress=None, authenticated_resources=None, _live_resources=None):
        from genosapiens.wgs_completed_science import require, fields
        from genosapiens.wgs_analysis_eligibility import scientific_assessment_domain, verify_analysis_eligibility_component
        from genosapiens.wgs_profile_projection import verify_analysis_profile_component
        require(reference.kind == 'scientific-case-resources-v1')
        value = reader.reconstruct(reference, max_bytes=MAX_CASE_BINDING_BYTES)
        fields(value, ('schema', 'version', 'genes', 'phenotypes', 'population', 'run_manifest_sha256',
            'view_sha256', 'profile_sha256', 'ic_sha256', 'hypothesis_policy', 'hypothesis_policy_sha256',
            'admission_policy_sha256', 'fusion_policy', 'fusion_policy_sha256', 'prior_policy', 'diagnosis',
            'sample_name', 'source_manifest_sha256', 'reference', 'parent_bindings',
            'outside_context_partner_uncertainty', 'analysis_eligibility', 'analysis_domain_status',
            'rescue_assessment_status', 'assessment_domain', 'analysis_profile'))
        require(value['schema'] == 'genosapiens.wgs-scientific-case-resources'
            and type(value['version']) is int and value['version'] == 4 and value['parent_bindings'] == expected_parent_bindings)
        require_saved_case_policy_bindings(value)
        self.reader, self.connection, self.name, self.sha256 = reader, reader.connection, reference.name, reference.sha256
        self.header, self.parent_bindings = value, dict(expected_parent_bindings)
        self._loaded, self.closed = None, False
        self._restore_arguments = dict(expected_parent_bindings=expected_parent_bindings,
            eligibility_root=eligibility_root, profile_root=profile_root, progress=progress,
            authenticated_resources=authenticated_resources)
        self.analysis_eligibility = self.analysis_profile = None
        if value['analysis_eligibility'] is not None:
            self.analysis_eligibility = reader.reconstruct(bound_reference(reader, value['analysis_eligibility'], 'analysis-eligibility-component-v1'),
                max_bytes=MAX_CASE_BINDING_BYTES)
            verify_analysis_eligibility_component(reader.directory.parent, self.analysis_eligibility, full=False,
                scratch_directory=reader.scratch_directory, eligibility_root=eligibility_root)
            require(all(value[k] == self.analysis_eligibility[k] for k in ('analysis_domain_status', 'rescue_assessment_status')))
        else:
            require(not any(k.startswith('analysis_eligibility_') for k in expected_parent_bindings)
                and value['analysis_domain_status'] == 'not_bound' and value['rescue_assessment_status'] == 'not_performed')
        if value['analysis_profile'] is not None:
            self.analysis_profile = reader.reconstruct(bound_reference(reader, value['analysis_profile'], 'analysis-profile-component-v1'),
                max_bytes=MAX_CASE_BINDING_BYTES)
            verify_analysis_profile_component(reader.directory.parent, self.analysis_profile, full=False,
                scratch_directory=reader.scratch_directory, profile_root=profile_root)
        else:
            require(not any(k.startswith('analysis_profile_') for k in expected_parent_bindings))
        self.assessment_domain = scientific_assessment_domain(self.analysis_eligibility, expected_parent_bindings, self.analysis_profile)
        require(value['assessment_domain'] == self.assessment_domain)
        self.assessment_domain_sha256 = _hash(self.assessment_domain)
        for key, cls in (('hypothesis_policy', HypothesisPolicy), ('fusion_policy', FusionPolicy),
                         ('prior_policy', PriorPolicy), ('diagnosis', DisclosedDiagnosis | None)):
            setattr(self, key, _typed(value[key], cls, decode=True))
        self.reference_binding = _typed(value['reference'], SubmissionReferenceBinding, decode=True)
        for key in ('hypothesis_policy_sha256', 'admission_policy_sha256', 'fusion_policy_sha256',
                    'source_manifest_sha256', 'view_sha256', 'profile_sha256', 'ic_sha256', 'run_manifest_sha256'):
            _sha(value[key])
            setattr(self, key, value[key])
        self.sample_name = value['sample_name']
        _ids((self.sample_name,))
        self.outside_context_partner_uncertainty = value['outside_context_partner_uncertainty']
        bound_reference(reader, value['genes'], 'case-gene-evidence-v1', extra=('logical_sha256',))
        _sha(self.gene_logical_sha256)
        bound_reference(reader, value['population'], 'population-policy-evidence-v1', extra=('policy_id', 'policy_sha256'))
        phenotypes = reader.reconstruct(bound_reference(reader, value['phenotypes'], 'case-phenotype-evidence-v1'))
        require(phenotypes['schema'] == 'genosapiens.wgs-case-phenotype-evidence'
            and type(phenotypes['version']) is int and phenotypes['version'] == 3
            and all(phenotypes[k] == value[k] for k in ('run_manifest_sha256', 'view_sha256', 'profile_sha256', 'ic_sha256')))
        profile = phenotypes['profile']
        if _live_resources is not None:
            require(type(_live_resources) is ScientificCaseResources
                and _live_resources.reference == self.reference
                and reader._live_prefix_source is not None
                and (_live_resources.writer is reader._live_prefix_source
                    or _live_resources.writer._copied_prefix_source is reader._live_prefix_source))
            _live_resources.require(_live_resources.writer)
            self._profile_evidence_id = _live_resources.profile_evidence_id
            require((profile is None) == (self._profile_evidence_id is None))
        elif profile is not None:
            profile = reader.reconstruct(bound_reference(reader, profile, 'case-phenotype-profile-v1'))
            require(profile['schema'] == 'genosapiens.wgs-case-phenotype-profile'
                and type(profile['version']) is int
                and profile['version'] == 2)
            self._profile_evidence_id = profile['evidence_id']
        else:
            self._profile_evidence_id = None
        require(self._profile_evidence_id is None or type(self._profile_evidence_id) is str and self._profile_evidence_id)

    @property
    def profile_evidence_id(self):
        return self._profile_evidence_id

    def require(self, writer):
        if self.closed or self.connection is None or writer.connection is not self.connection:
            raise EvidenceStoreError('completed scientific resource owner is closed or foreign')
        self.reader._check_files()

    def restore_for_selected(self):
        self.require(self.reader)
        if self._loaded is None:
            from genosapiens.wgs_runtime import analysis_progress
            analysis_progress(self._restore_arguments['progress'], 'publication', 'selected_detail_resources',
                counts=dict(scientific_resources_restored=0), force=True)
            self._loaded = VerifiedScientificCaseResources(self.reader,
                self.reader.reference(self.name, kind='scientific-case-resources-v1'), **self._restore_arguments)
            analysis_progress(self._restore_arguments['progress'], 'publication', 'selected_detail',
                counts=dict(scientific_resources_restored=1), force=True)
        return self._loaded

    def close(self):
        if self._loaded is not None:
            self._loaded.close()
        self.connection, self.closed = None, True

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


class VerifiedScientificCaseResources(ScientificCaseResources):
    """Replay resources from retained bytes, with independently supplied parents."""
    def __init__(self, reader, reference, *, expected_parent_bindings, eligibility_root=None, profile_root=None,
                 progress=None, authenticated_resources=None):
        try:
            with reader.read_scope():
                self._restore(reader,reference,expected_parent_bindings=expected_parent_bindings,
                    eligibility_root=eligibility_root,profile_root=profile_root,progress=progress,
                    authenticated_resources=authenticated_resources)
        except BaseException:
            if hasattr(self,"stack"):
                self.close()
            raise

    def _restore(self, reader, reference, *, expected_parent_bindings, eligibility_root=None, profile_root=None,
                 progress=None, authenticated_resources=None):
        from genosapiens.wgs_runtime import analysis_progress
        self.stack = ExitStack()
        try:
            if reference.kind != "scientific-case-resources-v1":
                raise EvidenceStoreError("foreign scientific case resource root")
            value = reader.reconstruct(reference,max_bytes=MAX_CASE_BINDING_BYTES)
            keys = {"schema","version","genes","phenotypes","population","run_manifest_sha256",
                "view_sha256","profile_sha256","ic_sha256","hypothesis_policy","hypothesis_policy_sha256",
                "admission_policy_sha256","fusion_policy","fusion_policy_sha256","prior_policy",
                "diagnosis","sample_name","source_manifest_sha256","reference","parent_bindings",
                "outside_context_partner_uncertainty","analysis_eligibility","analysis_domain_status",
                "rescue_assessment_status","assessment_domain","analysis_profile"}
            eligibility=None
            if value.get('analysis_eligibility') is not None:
                eligibility=reader.reconstruct(bound_reference(reader,value['analysis_eligibility'],'analysis-eligibility-component-v1'),
                    max_bytes=MAX_CASE_BINDING_BYTES)
                from genosapiens.wgs_analysis_eligibility import verify_analysis_eligibility_component
                verify_analysis_eligibility_component(reader.directory.parent,eligibility,full=False,
                    scratch_directory=reader.scratch_directory,eligibility_root=eligibility_root)
                if any(value[k]!=eligibility[k] for k in ('analysis_domain_status','rescue_assessment_status')):
                    raise EvidenceStoreError('scientific assessment domain changed')
            profile=None
            if value.get('analysis_profile') is not None:
                profile=reader.reconstruct(bound_reference(reader,value['analysis_profile'],'analysis-profile-component-v1'),
                    max_bytes=MAX_CASE_BINDING_BYTES)
                from genosapiens.wgs_profile_projection import verify_analysis_profile_component
                verify_analysis_profile_component(reader.directory.parent,profile,full=False,
                    scratch_directory=reader.scratch_directory,profile_root=profile_root)
            if (set(value) != keys or value["schema"] != "genosapiens.wgs-scientific-case-resources"
                    or type(value["version"]) is not int or value["version"]!=4
                    or value["parent_bindings"] != expected_parent_bindings):
                raise EvidenceStoreError("scientific case parent binding changed")
            from genosapiens.wgs_analysis_eligibility import scientific_assessment_domain
            if value['assessment_domain']!=scientific_assessment_domain(eligibility,expected_parent_bindings,profile):
                raise EvidenceStoreError('scientific assessment domain changed')
            require_saved_case_policy_bindings(value)
            analysis_progress(progress, 'publication', 'gene_resources', force=True)
            genes = VerifiedCaseGeneEvidence(reader,bound_reference(reader,value["genes"],"case-gene-evidence-v1",
                extra=("logical_sha256",)),expected_logical_sha256=value["genes"]["logical_sha256"])
            analysis_progress(progress, 'publication', 'phenotype_resources', force=True)
            phenotypes = VerifiedCasePhenotypeEvidence(reader,
                bound_reference(reader,value["phenotypes"],"case-phenotype-evidence-v1"),
                run_manifest_sha256=value["run_manifest_sha256"],expected_view_sha256=value["view_sha256"],
                expected_profile_sha256=value["profile_sha256"],expected_ic_sha256=value["ic_sha256"],
                ic_store=None if authenticated_resources is None else authenticated_resources.ic_store,
                authenticated_resources=authenticated_resources)
            analysis_progress(progress, 'publication', 'population_resources', force=True)
            population = self.stack.enter_context(VerifiedPopulationPolicy(reader,
                bound_reference(reader,value["population"],"population-policy-evidence-v1",
                    extra=("policy_id","policy_sha256")),
                expected_policy_sha256=value["population"]["policy_sha256"]))
            super().__init__(EvidenceReplay(reader),genes=genes,phenotypes=phenotypes,population=population,
                hypothesis_policy=_typed(value["hypothesis_policy"],HypothesisPolicy,decode=True),
                fusion_policy=_typed(value["fusion_policy"],FusionPolicy,decode=True),
                prior_policy=_typed(value["prior_policy"],PriorPolicy,decode=True),
                diagnosis=_typed(value["diagnosis"],DisclosedDiagnosis | None,decode=True),
                reference=_typed(value["reference"],SubmissionReferenceBinding,decode=True),
                name=reference.name,analysis_eligibility=eligibility,analysis_profile=profile,**{k:value[k] for k in ("hypothesis_policy_sha256","admission_policy_sha256",
                    "fusion_policy_sha256","sample_name","source_manifest_sha256","parent_bindings",
                    "outside_context_partner_uncertainty")})
        except BaseException:
            self.close()
            raise

    def close(self):
        self.stack.close()
        self.connection = None

    def __enter__(self):
        return self

    def __exit__(self,*_):
        self.close()
