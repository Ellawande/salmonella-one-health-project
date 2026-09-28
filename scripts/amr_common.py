"""Shared AMR classification logic, imported by scripts 02 and 06.

INTRINSIC_MARKERS lists CARD genes treated as chromosomal (intrinsic /
housekeeping) rather than horizontally acquired; is_acquired() is its accessor.

classify_pattern() turns the phylogenetic clumping test into the reporting
vocabulary used across the project.  Both scripts import from here so they can
never disagree about how a gene was classified -- the drift between them is
what made results/amr_inherited_vs_acquired.csv unreproducible.
"""
from __future__ import annotations

# --- reporting vocabulary -------------------------------------------------
PATTERN_SINGLE  = "single isolate (uninformative)"
PATTERN_NO_TEST = "invariant (no test)"
PATTERN_CLUMPED = "clumped (inherited)"
PATTERN_SCATTER = "scattered (repeated acquisition)"

CLUMPED_ALPHA = 0.05


INTRINSIC_MARKERS = (
    # --- aminoglycoside acetyltransferases that are chromosomal in enterobacteria
    "aac(6')-iaa", "aac(6')-iay", "aac(6')-iz", "aac(6')-iic",
    # --- RND / MFS / ABC efflux pumps (chromosomal).  "mdsabc" was a composite
    #     of three separate CARD names and so could never match any real hit;
    #     it is replaced by mdsa/mdsb/mdsc.
    "mdtk", "mdfa", "acra", "acrb", "acrd", "emra", "emrb", "emrd", "emry",
    "mdsa", "mdsb", "mdsc", "mdtb", "mdtc", "mdth", "msba", "yoji", "tolc",
    # --- efflux / stress regulators (chromosomal).  "ramA" was mixed-case and
    #     the match was case-sensitive, so it never fired; matching is now
    #     case-insensitive on both sides.
    "mara", "rama", "soxs", "emrr", "h-ns", "gols", "cpxa", "cpxr", "crp",
    "hfq", "roba", "sdia", "gadx", "gadw", "kdpe", "phop", "phoq",
    # --- LPS / membrane modification and peptide resistance (chromosomal)
    "ampr", "amps", "lpxc", "pmrf", "pmre", "ugd", "arna", "arnb",
    "baca", "baea", "baer", "epta", "eptb",
    # NOTE: "mcr-9" is kept deliberately -- it is chromosomally located in the
    # Salmonella it was first described in, though it is mobilisable.
    # NOTE: oqxA/oqxB are NOT listed: OqxAB is both plasmid- and chromosome-borne
    # in enterobacteria, so it is left classified as acquired pending a call.
    "mcr-9",
)


def is_acquired(label: str) -> bool:
    """True if a CARD hit looks plasmid-/transposon-borne rather than intrinsic."""
    low = label.lower()
    return not any(m.lower() in low for m in INTRINSIC_MARKERS)   # case-insensitive on both sides


def classify_pattern(n_isolates: int, p_clumped,
                     alpha: float = CLUMPED_ALPHA) -> str:
    """Turn the clumping test into a reporting label.

    A gene seen in only one isolate carries no phylogenetic information, so it
    is labelled uninformative before any p-value is considered.  When the test
    could not run (p is NaN -- the gene is invariant, so there is nothing to
    permute) that is stated explicitly rather than forced into either the
    inherited or repeated-acquisition bucket.
    """
    if n_isolates == 1:
        return PATTERN_SINGLE
    if p_clumped is None or p_clumped != p_clumped:       # NaN, no pandas needed
        return PATTERN_NO_TEST
    return PATTERN_CLUMPED if p_clumped < alpha else PATTERN_SCATTER
