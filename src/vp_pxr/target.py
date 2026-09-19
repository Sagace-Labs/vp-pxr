"""The endpoints this pathway predicts.

Frozen so the data sources are auditable.

Both assays were verified live against PubChem on 2026-09-06, and screen the
same library:

    AID 1346982  Human pregnane X receptor (PXR) small molecule agonists,
                 qHTS assay                                9667 substances
    AID 1346977  the same, cell viability counter screen   9667 substances
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["CYTOTOX", "TARGET", "TARGETS", "Endpoint", "all_names", "get"]


@dataclass(frozen=True)
class Endpoint:
    """One molecular initiating event, identified by the assay that reads it out."""

    name: str
    pathway: str
    mie: str
    pubchem_aid: int
    assay_name: str


TARGET = Endpoint(
    name="PXR",
    pathway="pregnane X receptor / xenobiotic induction",
    mie=(
        "binding of the pregnane X receptor ligand-binding pocket, which induces "
        "CYP3A4 and the transporters that accompany it"
    ),
    pubchem_aid=1346982,
    assay_name=(
        "Human pregnane X receptor (PXR) small molecule agonists, qHTS assay"
    ),
)

CYTOTOX = Endpoint(
    name="VIABILITY",
    pathway="pregnane X receptor / xenobiotic induction",
    mie=(
        "loss of cell viability, which registers on the reporter readout as a "
        "consequence"
    ),
    pubchem_aid=1346977,
    assay_name=(
        "Human pregnane X receptor (PXR) small molecule agonists, qHTS cell "
        "viability counter screen"
    ),
)

TARGETS: dict[str, Endpoint] = {"PXR": TARGET, "VIABILITY": CYTOTOX}


def get(name: str) -> Endpoint:
    try:
        return TARGETS[name.upper()]
    except KeyError:
        raise KeyError(f"unknown target {name!r}; known: {sorted(TARGETS)}") from None


def all_names() -> list[str]:
    return sorted(TARGETS)
