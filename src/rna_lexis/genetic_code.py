"""The genetic code, and the codon -> amino-acid projection.

Translating codons to amino acids is a many-to-one *projection* of one alphabet
onto another (:class:`rna_lexis.alphabet.Projection`): 64 codons become 20 amino
acids and the stop signal.  Analyzing the projected text asks about the protein;
analyzing the codon text asks about the protein *and* the choice among synonymous
codons.  Running both and comparing is the direct way to see whether synonymous
choice matters.

The code is a replaceable table: ``TABLES`` maps a name to a codon table, and the
standard code is the only one included so far (the others, such as the
vertebrate mitochondrial code, are a table away; see the wishlist in
``docs/scoping_codon_units.md``).
"""

from __future__ import annotations

import itertools
from typing import Mapping

from rna_lexis.alphabet import Alphabet, Projection, codon

# One-letter amino-acid codes, lower case, and the stop signal.
AMINO_ACIDS = "acdefghiklmnpqrstvwy"
STOP = "*"

# The standard code, codons in T, C, A, G order (NCBI table 1).
_STANDARD_STRING = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
STANDARD: Mapping[str, str] = {
    "".join(p).lower(): aa.lower()
    for p, aa in zip(itertools.product("TCAG", repeat=3), _STANDARD_STRING)
}

TABLES = {"standard": STANDARD}


def amino_acid() -> Alphabet:
    """The 20 amino acids and the stop signal, one character each.

    Used for protein text (one-letter codes) and as the target of the
    codon-to-amino-acid projection.
    """
    return Alphabet(name="amino acid", symbols=tuple(AMINO_ACIDS) + (STOP,), mode="chars",
                    unit_name="residue")


def codon_to_amino_acid(table: str = "standard") -> Projection:
    """The projection from the 64 codons to amino acids under a genetic code.

    Args:
        table: Name of a code in :data:`TABLES` (``"standard"``).
    """
    try:
        code = TABLES[table]
    except KeyError:
        raise ValueError(f"unknown genetic code {table!r}; available: {sorted(TABLES)}") from None
    return Projection(source=codon(), target=amino_acid(), mapping=dict(code),
                      name=f"codon -> amino acid ({table} code)")


def translate_codons(internal: str, table: str = "standard") -> str:
    """Internal codon string -> internal amino-acid string."""
    return codon_to_amino_acid(table).apply(internal)
