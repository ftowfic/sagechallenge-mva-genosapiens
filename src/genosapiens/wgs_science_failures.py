"""Classify failures only after an independent context owns its output."""
from genosapiens.population_filter_candidate import OptionalPopulationFilterError, PopulationCandidateBindingError
from genosapiens.wgs_science_coordinator import ContextEvaluationFailure, ContextFailure
from genosapiens.wgs_science_runtime import ScienceBudget, ScienceClock, ScienceDeadlineExpired
from genosapiens.wgs_scientific_artifact import _ScienceFailureObservation

LOCAL_STAGES = frozenset({'member_enrichment', 'member_read', 'population_query', 'population_bind',
    'clinical_query', 'clinical_bind', 'member_encode', 'genotype_facts', 'hypothesis_evidence',
    'population', 'phenotype', 'scores', 'representation', 'eligibility', 'ordering'})
LOCAL_QUERIES = frozenset({'population_query', 'clinical_index_query', 'context_population_query', 'context_clinical_query'})


def observe_setup_failure(error):
    """Emit bounded installed-code locations, never exception text or values."""
    import json
    try:
        frames = []
        seen = set()
        truncated = False
        while error is not None and id(error) not in seen and len(seen) < 4:
            seen.add(id(error))
            current, limited = _ScienceFailureObservation._frames(error)
            frames.extend(current)
            truncated |= limited
            error = error.__cause__ if error.__cause__ is not None else error.__context__
        record = dict(schema='genosapiens.science-setup-failure/v1', diagnostic_only=True,
            source_frames=frames[-32:], source_frames_truncated=truncated or len(frames)>32 or error is not None)
        print(json.dumps(record, sort_keys=True), flush=True)
    except BaseException:
        # Diagnostic failure cannot replace the original failure or authorize reuse.
        pass


class IndependentScienceFailureObservation(_ScienceFailureObservation):
    def __init__(self, root, run_manifest_sha256, progress, *, budget=None):
        if budget is not None and type(budget) not in (ScienceBudget, ScienceClock):
            raise ValueError('scientific failure timing owner is invalid')
        super().__init__(root, run_manifest_sha256, progress, terminal=False)
        self.classified = None
        self.budget = budget

    def capture(self, error):
        if self.classified is not None:
            return self.classified
        primary = super().capture(error)
        overall_deadline = (type(error) is ScienceDeadlineExpired
            and (self.budget is None or self.budget.remaining() <= 0))
        if (overall_deadline or type(error) in (OptionalPopulationFilterError, ContextEvaluationFailure)
                or isinstance(error, (KeyboardInterrupt, SystemExit, GeneratorExit))):
            self.classified = error
            return error
        scope = 'unknown'
        if type(error) is PopulationCandidateBindingError or primary.category in {'authentication', 'storage'}:
            scope = 'shared'
        elif self.subphase in LOCAL_STAGES:
            if primary.category == 'response_capacity':
                scope = 'context'
            elif primary.category in {'missing_response', 'timeout'} and primary.operation in LOCAL_QUERIES:
                scope = 'context'
        self.classified = ContextEvaluationFailure(ContextFailure(self.subphase, primary.operation, primary.failure_code, scope))
        return self.classified
