"""Cross-source active-ingredient and product identity resolution."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Iterable


SOURCES = {
    "chembl",
    "drugs-at-fda",
    "fda-label",
    "fda-product-status",
    "gtopdb",
    "lincs",
    "pubchem",
    "purple-book",
    "rxnorm",
}
CONCEPT_TYPES = {
    "active_ingredient",
    "chemical_class",
    "combination",
    "non_medication",
    "product",
    "salt",
    "stereoisomer",
}
STEREOCHEMISTRY_STATES = {"defined", "not_applicable", "unspecified"}


class DrugIdentityError(ValueError):
    """Raised when source identity assertions are contradictory or malformed."""


def _stable_id(prefix: str, value: str) -> str:
    return f"{prefix}:" + hashlib.sha256(value.encode()).hexdigest()[:24]


@dataclass(frozen=True, slots=True)
class DrugIdentityAssertion:
    source: str
    source_id: str
    preferred_name: str
    concept_type: str
    exact_substance_key: str | None = None
    active_moiety_keys: tuple[str, ...] = ()
    ingredient_keys: tuple[str, ...] = ()
    candidate_ingredient_sets: tuple[tuple[str, ...], ...] = ()
    stereochemistry: str = "not_applicable"
    regulator_product_id: str | None = None

    def __post_init__(self) -> None:
        if (
            self.source not in SOURCES
            or not self.source_id
            or not self.preferred_name
            or self.concept_type not in CONCEPT_TYPES
            or self.stereochemistry not in STEREOCHEMISTRY_STATES
            or len(self.active_moiety_keys) != len(set(self.active_moiety_keys))
            or len(self.ingredient_keys) != len(set(self.ingredient_keys))
            or any(not value for value in self.active_moiety_keys + self.ingredient_keys)
        ):
            raise DrugIdentityError("drug identity assertion is incomplete or invalid")
        if self.concept_type in {"active_ingredient", "salt", "stereoisomer"}:
            if self.exact_substance_key is None:
                raise DrugIdentityError("chemical concepts require an exact substance key")
        if self.concept_type == "salt" and len(self.active_moiety_keys) != 1:
            raise DrugIdentityError("a salt requires exactly one active-moiety identity")
        if self.concept_type == "stereoisomer" and self.stereochemistry != "defined":
            raise DrugIdentityError("a stereoisomer requires defined stereochemistry")
        if self.concept_type in {"combination", "product"} and not (
            self.ingredient_keys or self.candidate_ingredient_sets
        ):
            raise DrugIdentityError("combination and product records require ingredients")
        if self.concept_type == "product" and self.regulator_product_id is None:
            raise DrugIdentityError("product records require a regulator product identity")
        if self.concept_type in {"chemical_class", "non_medication"} and any(
            (
                self.exact_substance_key,
                self.active_moiety_keys,
                self.ingredient_keys,
                self.candidate_ingredient_sets,
            )
        ):
            raise DrugIdentityError(
                "a class or non-medication cannot masquerade as an ingredient"
            )

    @property
    def identifier(self) -> str:
        return f"{self.source}:{self.source_id}"


@dataclass(frozen=True, slots=True)
class DrugIdentityResolution:
    query_identifier: str
    state: str
    concept_type: str | None
    exact_substance_id: str | None
    active_ingredient_ids: tuple[str, ...]
    source_identifiers: tuple[str, ...]
    product_identifiers: tuple[str, ...]
    candidate_ingredient_sets: tuple[tuple[str, ...], ...]
    ambiguity_reasons: tuple[str, ...]


def _ingredient_ids(keys: Iterable[str]) -> tuple[str, ...]:
    return tuple(sorted({_stable_id("active-ingredient", key) for key in keys}))


def resolve_drug_identity(
    identifier: str,
    assertions: Iterable[DrugIdentityAssertion],
) -> DrugIdentityResolution:
    records = tuple(assertions)
    matches = tuple(record for record in records if record.identifier == identifier)
    if not matches:
        return DrugIdentityResolution(
            query_identifier=identifier,
            state="unresolved",
            concept_type=None,
            exact_substance_id=None,
            active_ingredient_ids=(),
            source_identifiers=(),
            product_identifiers=(),
            candidate_ingredient_sets=(),
            ambiguity_reasons=("identifier_not_found",),
        )
    signatures = {
        (
            record.concept_type,
            record.exact_substance_key,
            record.active_moiety_keys,
            record.ingredient_keys,
            record.candidate_ingredient_sets,
            record.stereochemistry,
            record.regulator_product_id,
        )
        for record in matches
    }
    if len(signatures) != 1:
        return DrugIdentityResolution(
            query_identifier=identifier,
            state="ambiguous",
            concept_type=None,
            exact_substance_id=None,
            active_ingredient_ids=(),
            source_identifiers=tuple(sorted({item.identifier for item in matches})),
            product_identifiers=(),
            candidate_ingredient_sets=(),
            ambiguity_reasons=("conflicting_source_assertions",),
        )
    record = matches[0]
    if record.concept_type == "chemical_class":
        return DrugIdentityResolution(
            query_identifier=identifier,
            state="class_only",
            concept_type=record.concept_type,
            exact_substance_id=None,
            active_ingredient_ids=(),
            source_identifiers=(identifier,),
            product_identifiers=(),
            candidate_ingredient_sets=(),
            ambiguity_reasons=("class_is_not_an_active_ingredient",),
        )
    if record.concept_type == "non_medication":
        return DrugIdentityResolution(
            query_identifier=identifier,
            state="non_medication",
            concept_type=record.concept_type,
            exact_substance_id=None,
            active_ingredient_ids=(),
            source_identifiers=(identifier,),
            product_identifiers=(),
            candidate_ingredient_sets=(),
            ambiguity_reasons=("concept_is_not_a_medication",),
        )
    if record.candidate_ingredient_sets:
        candidates = tuple(
            sorted({_ingredient_ids(candidate) for candidate in record.candidate_ingredient_sets})
        )
        return DrugIdentityResolution(
            query_identifier=identifier,
            state="ambiguous",
            concept_type=record.concept_type,
            exact_substance_id=(
                _stable_id("substance", record.exact_substance_key)
                if record.exact_substance_key
                else None
            ),
            active_ingredient_ids=(),
            source_identifiers=(identifier,),
            product_identifiers=(
                (f"{record.source}:{record.regulator_product_id}",)
                if record.regulator_product_id
                else ()
            ),
            candidate_ingredient_sets=candidates,
            ambiguity_reasons=("multiple_ingredient_interpretations",),
        )
    if record.concept_type == "active_ingredient":
        ingredient_keys = (record.exact_substance_key,)
    elif record.concept_type in {"salt", "stereoisomer"}:
        ingredient_keys = record.active_moiety_keys or (record.exact_substance_key,)
    else:
        ingredient_keys = record.ingredient_keys
    active_ids = _ingredient_ids(key for key in ingredient_keys if key is not None)
    equivalent_sources = tuple(
        sorted(
            {
                item.identifier
                for item in records
                if record.exact_substance_key is not None
                and item.exact_substance_key == record.exact_substance_key
            }
            | {identifier}
        )
    )
    products = tuple(
        sorted(
            {
                f"{item.source}:{item.regulator_product_id}"
                for item in records
                if item.regulator_product_id is not None
                and item.ingredient_keys == tuple(ingredient_keys)
            }
        )
    )
    state = {
        "active_ingredient": "resolved_active_ingredient",
        "combination": "resolved_combination",
        "product": "resolved_product",
        "salt": "resolved_salt",
        "stereoisomer": "resolved_stereoisomer",
    }[record.concept_type]
    return DrugIdentityResolution(
        query_identifier=identifier,
        state=state,
        concept_type=record.concept_type,
        exact_substance_id=(
            _stable_id("substance", record.exact_substance_key)
            if record.exact_substance_key
            else None
        ),
        active_ingredient_ids=active_ids,
        source_identifiers=equivalent_sources,
        product_identifiers=products,
        candidate_ingredient_sets=(),
        ambiguity_reasons=(),
    )
