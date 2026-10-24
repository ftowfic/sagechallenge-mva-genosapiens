"""One bounded offline stage session for population and clinical resources."""

from __future__ import annotations

from contextlib import AbstractContextManager, ExitStack
import os
from pathlib import Path

from genosapiens.common_allele_public_indexes import (
    MAX_ROW_BYTES,
    CommonAlleleIndexError,
    CommonAllelePublicIndexSession,
    PublicIndexQuery,
    PublicIndexQueryResult,
    VerifiedPublicIndexes,
)
from genosapiens.bounded_wgs_config import PublicQueryLimits
from genosapiens.wgs_query_limits import QueryBudget
from genosapiens.population_annotation import (
    GnomadAllSitesReference,
    GnomadJointReference,
    PopulationAnnotationResult,
    PopulationQueryLimits,
    PopulationReferenceCompatibility,
    PopulationRequest,
    PopulationResourceSession,
)
from genosapiens.wgs_failures import static_failure

def _identity(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)


class PopulationClinicalResourceSession(
    AbstractContextManager["PopulationClinicalResourceSession"]
):
    """Reuse authenticated resources across bounded source-ledger batches.

    The caller supplies a qualified public index whose canonical projection and
    SQLite accelerator reconciled at publication. This session pins their file
    metadata while population contracts are checked once and selected published
    data files are opened lazily under the same metadata-trust boundary.
    """

    def __init__(
        self, *, joint: GnomadJointReference, all_sites: GnomadAllSitesReference,
        reference: PopulationReferenceCompatibility, runtime_image: str,
        limits: PopulationQueryLimits, public_query_limits: PublicQueryLimits,
        clinical_indexes: VerifiedPublicIndexes,
    ) -> None:
        if (
            type(clinical_indexes) is not VerifiedPublicIndexes
            or type(public_query_limits) is not PublicQueryLimits
            or (
                limits.max_alleles,
                limits.max_output_bytes,
                limits.timeout_seconds,
            ) != (
                public_query_limits.batch_rows,
                public_query_limits.max_response_bytes,
                public_query_limits.timeout_seconds,
            )
        ):
            raise CommonAlleleIndexError("common-allele public index validation failed")
        self._owner_pid = os.getpid()
        self._population_arguments = dict(
            joint=joint, all_sites=all_sites, reference=reference,
            runtime_image=runtime_image, limits=limits,
        )
        self._limits = limits
        self._clinical_limits = public_query_limits
        self._verified = clinical_indexes
        self._stack: ExitStack | None = None
        self._population: PopulationResourceSession | None = None
        self._clinical: CommonAllelePublicIndexSession | None = None
        self._clinical_files: dict[Path, tuple[int, int, int, int, int]] = {}

    def _assert_active(self) -> None:
        if self._owner_pid != os.getpid() or self._stack is None or self._population is None or self._clinical is None:
            raise RuntimeError("population/clinical resource session is not active")

    def _assert_clinical_stable(self) -> None:
        try:
            for path, expected in self._clinical_files.items():
                if _identity(path.stat()) != expected:
                    raise CommonAlleleIndexError("common-allele public index validation failed")
        except OSError:
            raise CommonAlleleIndexError("common-allele public index validation failed") from None

    def __enter__(self) -> "PopulationClinicalResourceSession":
        if self._owner_pid != os.getpid():
            raise RuntimeError("resource session must be constructed in its owner process")
        if self._stack is not None:
            raise RuntimeError("population/clinical resource session is already active")
        stack = ExitStack()
        try:
            self._clinical_files = {
                path: _identity(path.stat())
                for path in (self._verified.projection_path, self._verified.sqlite_path)
            }
            self._population = stack.enter_context(PopulationResourceSession(
                **self._population_arguments
            ))
            # Authenticate contracts even if every occurrence takes a native bypass.
            self._population.query(())
            self._clinical = stack.enter_context(CommonAllelePublicIndexSession(
                self._verified, limits=self._clinical_limits,
            ))
            self._stack = stack
            self._assert_clinical_stable()
            return self
        except (Exception, KeyboardInterrupt) as error:
            primary = static_failure(
                error,
                "resource_session",
                operation="resource_session_open",
            )
            try:
                stack.close()
            except (Exception, KeyboardInterrupt) as cleanup_error:
                primary.add_secondary(static_failure(
                    cleanup_error,
                    "resource_session",
                    operation="resource_session_close",
                ))
            self._population = None
            self._clinical = None
            self._clinical_files = {}
            raise primary from None

    def query_population(
        self, requests: tuple[PopulationRequest, ...]
    ) -> PopulationAnnotationResult:
        self._assert_active()
        self._assert_clinical_stable()
        result = self._population.query(requests)
        self._assert_clinical_stable()
        return result

    def query_population_with_noncontribution(self, requests, *, _case_cache=None):
        """Explicit supplemental proof query under the same pinned session."""
        self._assert_active()
        self._assert_clinical_stable()
        result = self._population.query_with_noncontribution(requests, _case_cache=_case_cache)
        self._assert_clinical_stable()
        return result

    def _population_screening_facts(self, requests, *, _case_cache=None):
        """Borrow compact facts from the same strict population adapter."""
        self._assert_active()
        self._assert_clinical_stable()
        result = self._population._screening_facts(requests, _case_cache=_case_cache)
        self._assert_clinical_stable()
        return result

    def query_clinical(
        self, queries: tuple[PublicIndexQuery, ...]
    ) -> tuple[PublicIndexQueryResult, ...]:
        self._assert_active()
        budget = QueryBudget(self._clinical_limits, max_row_bytes=MAX_ROW_BYTES, origin="clinical_index")
        if (type(queries) is not tuple
                or any(type(query) is not PublicIndexQuery for query in queries)):
            raise CommonAlleleIndexError("common-allele public index validation failed")
        self._assert_clinical_stable()
        results = self._clinical.query_many(queries, budget=budget)
        self._assert_clinical_stable()
        return results

    def __exit__(self, *args: object) -> None:
        self._assert_active()
        stack = self._stack
        try:
            self._assert_clinical_stable()
        finally:
            self._stack = None
            self._population = None
            self._clinical = None
            self._clinical_files = {}
            stack.__exit__(*args)
