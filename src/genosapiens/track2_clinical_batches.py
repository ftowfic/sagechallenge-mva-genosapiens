"""Keep Track 2 clinical batches within the existing public query limits."""

from genosapiens.admission_links import _inputs
from genosapiens.wgs_clinical_session import _canonical
from genosapiens.wgs_failures import WgsFailure


_CAPACITY_OPERATIONS = {
    'clinical_index': 'clinical_index_query',
    'clinical_replay': 'clinical_response_bind',
    'clinical_adaptation': 'clinical_response_adapt',
}


def _enrich(session, batch):
    """Publish results only after both operations succeed for the same batch."""
    associations = None
    try:
        associations = session.association_batch(batch)
        annotation = session.annotate(batch)
    except WgsFailure as error:
        if (len(batch) <= 1 or error.category != 'response_capacity'
                or error.operation != _CAPACITY_OPERATIONS.get(error.origin)
                or error.exception_kind != 'project_failure'
                or error.secondary_failures or error.secondary_failures_truncated):
            raise
    else:
        yield associations, annotation
        return
    # Source proofs and clinical evidence can exceed the input byte estimate.
    # Release the failed batch's results before retrying. Single-row failures
    # remain failures, including rows that exceed the unchanged row limit.
    associations = None
    middle = len(batch) // 2
    yield from _enrich(session, batch[:middle])
    yield from _enrich(session, batch[middle:])


def clinical_batches(session, alleles):
    """Bound input bytes without changing the session's existing row order."""
    limits = session.execution_config.public_queries
    # The session sorts each count-bounded batch by allele ID. Sort before
    # splitting so new byte boundaries cannot change the published order.
    for offset in range(0, len(alleles), limits.batch_rows):
        rows = _inputs(alleles[offset:offset + limits.batch_rows])
        batch, size = [], 0
        for allele in rows:
            row_size = len(_canonical(allele))
            if batch and size + row_size > limits.max_response_bytes:
                yield from _enrich(session, tuple(batch))
                batch, size = [], 0
            batch.append(allele)
            size += row_size
        if batch:
            yield from _enrich(session, tuple(batch))
