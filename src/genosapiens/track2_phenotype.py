"""Optional, sourced HPO review context; never a genomic selection input."""
from dataclasses import dataclass
import json

from genosapiens.bounded_json import read_bounded_bytes
from genosapiens.evidence_model import Observation
from genosapiens.evidence_serialization import _plain
from genosapiens.ontology_bundle import OBSOLETE_RESOLUTION_POLICY, ParsedOntology
from genosapiens.phenotype_ingestion import ingest_phenotypes, MAX_HPO_LIST_BYTES
from genosapiens.phenotype_scoring import ObservationProvenance, ScoringObservation, score_phenotypes
from genosapiens.track2_published_genotype import canonical, sha256


class OptionalPhenotypeError(ValueError):
    """Static failure; never log case observations or input paths."""


def require(ok, message='optional phenotype input is invalid or changed'):
    if not ok:
        raise OptionalPhenotypeError(message)


@dataclass(frozen=True)
class OptionalPhenotype:
    component: dict
    observations: tuple[Observation, ...]


def ingest_optional_phenotype(path=None, *, expected_input_sha256=None,
        hpo_ontology=None, ontology_manifest_sha256=None, input_format='hpo-list'):
    if path is None:
        require(expected_input_sha256 is None)
        return OptionalPhenotype({'state':'not_assessed','source_sha256':None,'observations':[],
            'family_context':[],'limitations':['No HPO supplied; patient manifestation relevance is not assessed.',
                'Unlisted features remain unknown.']},())
    require(type(hpo_ontology) is ParsedOntology and isinstance(expected_input_sha256,str)
        and len(expected_input_sha256)==64 and isinstance(ontology_manifest_sha256,str)
        and len(ontology_manifest_sha256)==64)
    if input_format=='hpo-list':
        observations=ingest_phenotypes(path,expected_input_sha256=expected_input_sha256,
            hpo_ontology=hpo_ontology,ontology_manifest_sha256=ontology_manifest_sha256).observations
    else:
        require(input_format=='structured-observations')
        raw=read_bounded_bytes(path,maximum_bytes=MAX_HPO_LIST_BYTES)
        require(sha256(raw)==expected_input_sha256)
        try:
            payload=json.loads(raw)
            require(set(payload)=={'schema','records'} and payload['schema']=='genosapiens.track2-hpo-observations/v1'
                and isinstance(payload['records'],list) and payload['records'])
            observations=[]
            optional={'subject_id','onset','age_years','sex','domain','dependency_cluster','quality'}
            for index,row in enumerate(payload['records'],1):
                require(isinstance(row,dict) and {'source_hpo_id','scope','status'}<=set(row)
                    and set(row)<={'source_hpo_id','scope','status',*optional})
                provenance=ObservationProvenance(expected_input_sha256,'record:'+str(index),'structured')
                resolution=hpo_ontology.resolve_term(row['source_hpo_id'])
                identity=sha256(canonical({'input':expected_input_sha256,'record':index,
                    'resolution':_plain(resolution),'ontology':ontology_manifest_sha256}))
                scoring=ScoringObservation('observation:'+identity[:24],row['source_hpo_id'],row['scope'],row['status'],
                    provenance,**{k:row[k] for k in optional if k in row})
                observations.append(Observation(scoring,resolution,ontology_manifest_sha256,
                    OBSOLETE_RESOLUTION_POLICY,'typed-structured',context_provenance=(provenance,),source_fields=row))
            observations=tuple(observations)
        except (ValueError,TypeError,KeyError,AttributeError):
            raise OptionalPhenotypeError('structured phenotype input cannot be resolved') from None
    family=tuple(o for o in observations if o.scoring.scope=='family_history')
    patient=tuple(o for o in observations if o.scoring.scope!='family_history')
    limitations=['HPO is review context only; it cannot change disease/gene/interval scope, variant membership or genotype order.',
        'Unlisted features remain unknown; family context is separate from patient observations.',
        'Phenotype compatibility does not establish causal direction, reversibility or drug benefit.']
    if any(o.resolution.canonical_id is None for o in observations):
        limitations.append('Supplied unresolved HPO observations remain sourced and are not negative findings.')
    return OptionalPhenotype({'state':'assessed','source_sha256':expected_input_sha256,
        'observations':_plain(patient),'family_context':_plain(family),'limitations':limitations},observations)


def contextual_phenotype_scores(phenotype, profile_ids, *, ic_store=None, profile_store=None):
    """One shared comparison per disease profile, never one score per variant."""
    require(type(phenotype) is OptionalPhenotype)
    profiles=tuple(sorted(set(profile_ids)))
    if phenotype.component['state']=='not_assessed':
        return {'state':'not_assessed','profiles':{},'ordering_effect':'none',
                'endpoint_relevance':'patient_manifestations_not_assessed'}
    if ic_store is None or profile_store is None:
        return {'state':'not_assessed','profiles':{},'ordering_effect':'none',
                'endpoint_relevance':'observations_retained; comparison_sources_unavailable'}
    patient=tuple(o.scoring for o in phenotype.observations if o.scoring.scope!='family_history')
    family=tuple(o.scoring for o in phenotype.observations if o.scoring.scope=='family_history')
    rows={key:{'patient':None if not patient else _plain(score_phenotypes(patient,key,
        ic_store=ic_store,profile_store=profile_store)),
        'family_context':None if not family else _plain(score_phenotypes(family,key,
        ic_store=ic_store,profile_store=profile_store))} for key in profiles}
    return {'state':'assessed','profiles':rows,'ordering_effect':'none',
        'endpoint_relevance':'context_for_review_only; benefit_and_reversibility_unassessed'}


def render_phenotype_context(component):
    if component['state']=='not_assessed':
        return '## Phenotype context\n\nThe input supplies no Human Phenotype Ontology (HPO) terms. Relevance to patient features is **not assessed**. Disease-level objectives remain available for review.\n'
    observations=component['observations']
    present=sum(o['scoring']['status']=='present' for o in observations)
    absent=sum(o['scoring']['status']=='explicitly_absent' for o in observations)
    unresolved=sum(o['resolution']['canonical_id'] is None for o in observations)
    return ('## Phenotype context\n\n'
        f'{present} present and {absent} explicitly absent patient observations. The report retains {unresolved} unresolved identifiers. '
        f"{len(component['family_context'])} family observations appear separately.\n\n"
        'HPO provides review context only. It does not alter genomic scope, variant membership or genotype order. '
        'Comparisons against the same disease profile cannot distinguish its variants. '
        'Unlisted features, patient-specific endpoint relevance and benefit remain unassessed unless separately supported.\n')
