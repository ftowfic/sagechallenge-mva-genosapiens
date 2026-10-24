"""Finite acceptance metadata backed by complete normalized scientific roots."""

from contextlib import contextmanager
from pathlib import Path

from genosapiens.wgs_scientific_products import open_scientific_products, product_descriptors
from genosapiens.wgs_validation_projection import NormalizedValidationProjection


@contextmanager
def open_validation_projection(products, audit, *, audit_directory, scratch_path,
                               epcr_assignment, population_selector=None, scientific_reader=None):
    """Reuse an operation's verified reader, or own one bounded disk replay."""
    arguments = dict(audit_directory=audit_directory, epcr_assignment=epcr_assignment,
        population_selector=population_selector)
    if scientific_reader is not None:
        with NormalizedValidationProjection(scientific_reader, products, audit, **arguments) as projection:
            yield projection
    else:
        with open_scientific_products(products.scientific, product_descriptors(products),
                directory=products.scientific_directory,
                expected_parent_bindings=dict(products.scientific.parent_bindings),
                scratch_directory=Path(scratch_path).parent) as reader:
            with NormalizedValidationProjection(reader, products, audit, **arguments) as projection:
                yield projection


def bounded_validation_metrics(products, audit, **arguments):
    """Keep selected groups and counts; bind complete recovery to saved roots."""
    with open_validation_projection(products, audit, **arguments) as projection:
        return projection.metrics()
