"""reconcILS: gene-tree/species-tree reconciliation with ILS."""
from __future__ import annotations

__all__ = ["reconcils", "reconcILS_validate", "SpeciesCopyMap"]


def __getattr__(name: str):
    if name == "reconcils":
        from .reconcILS import reconcils

        return reconcils
    if name == "reconcILS_validate":
        from .reconcILS_validate import reconcILS_validate

        return reconcILS_validate
    if name == "SpeciesCopyMap":
        from .copy_labels import SpeciesCopyMap

        return SpeciesCopyMap
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
