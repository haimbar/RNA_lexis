"""Alphabets: the declared set of symbols a sequence is analyzed over.

An :class:`Alphabet` is declared, not observed.  Every symbol in it counts
whether or not it occurs in a given text, so an English document with no
``z`` is still analyzed over 26 letters, and a transcript that never uses one
of the 64 codons is still analyzed over 64.

Each symbol (a *unit*) is given one internal character, so the analysis code
can work on one-character-per-unit strings and report positions and lengths in
units.  The mapping to and from the external labels (``"atg"``, ``"a"``,
``"ba"``) lives here; :mod:`rna_lexis.tokenizers` turns text into internal
strings and :mod:`rna_lexis.codec` converts between units and nucleotides.

Internal characters
-------------------
A symbol whose label is a single lower-case ASCII letter is its own internal
character, so a nucleotide or Latin-letter text reads exactly as before.  Any
other symbol (a codon, a digit, a space, a syllable name) gets a character
from the Unicode private use area, which has no upper or lower case and so is
unaffected by ``lower()``/``upper()``.  Characters are assigned in the order
the symbols are declared, so the same alphabet always maps the same way.

The segment boundary is the internal character ``"_"``, as elsewhere in the
package, and ``WILDCARD`` (the last private use character) means "any symbol" in
a pattern.  Neither is ever a symbol of any alphabet.
"""

from __future__ import annotations

import itertools
import json
import os
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Optional

BOUNDARY = "_"
"""Internal character that separates segments; never an alphabet symbol."""

WILDCARD = "\uf8ff"
"""Internal character that stands for "any symbol" in a pattern; never an alphabet symbol.

The nucleotide code writes wildcards as ``n``, ``.``, ``?`` or ``*``, but those are
ordinary symbols in other alphabets (``n`` is a Latin letter), so patterns over an
:class:`Alphabet` use this reserved private use character instead."""

_PRIVATE_USE_START = 0xE000
_PRIVATE_USE_END = 0xF8FE            # U+F8FF is the wildcard

TOKENIZATION_MODES = ("chars", "fixed", "words")


def _is_identity_label(label: str) -> bool:
    """True for a single lower-case ASCII letter, which is its own character."""
    return len(label) == 1 and "a" <= label <= "z"


@dataclass(frozen=True)
class Alphabet:
    """A declared set of symbols plus how text is read into it.

    Attributes:
        name:           Display name, used in titles and result files.
        symbols:        The labels of the main symbols, in declared order.
        extras:         Extra single-character symbols (a space, digits,
                        punctuation) that belong to the alphabet beside the
                        main symbols.  They count in the statistics and can
                        appear inside motifs.
        codes:          Optional nucleotide code for each main symbol (for
                        the codon alphabet the code is the label itself).
                        Needed only to convert to and from nucleotide text.
        mode:           Default tokenization: ``"chars"`` (one character per
                        symbol), ``"fixed"`` (``width`` characters per symbol)
                        or ``"words"`` (whitespace-separated labels).
        width:          Characters per symbol when ``mode`` is ``"fixed"``.
        case_sensitive: When False (default) labels must be lower-case and
                        input text is lower-cased before it is read.
        aliases:        Input read as something else, e.g. ``u`` -> ``t``.  In
                        ``"chars"`` and ``"words"`` mode the target must be a
                        symbol; in ``"fixed"`` mode an alias is one character
                        for another, applied before characters are grouped.
        defaults:       Analysis defaults that belong to this alphabet
                        (background order, minimum motif length, ...).  Kept
                        as given; the analysis code decides how to use them.
    """

    name: str
    symbols: tuple
    extras: tuple = ()
    codes: Mapping = field(default_factory=dict)
    mode: str = "chars"
    width: int = 1
    case_sensitive: bool = False
    aliases: Mapping = field(default_factory=dict)
    defaults: Mapping = field(default_factory=dict)

    def __post_init__(self):
        object.__setattr__(self, "symbols", tuple(self.symbols))
        object.__setattr__(self, "extras", tuple(self.extras))
        object.__setattr__(self, "codes", dict(self.codes))
        object.__setattr__(self, "aliases", dict(self.aliases))
        object.__setattr__(self, "defaults", dict(self.defaults))
        self._validate()
        self._assign_characters()

    # ------------------------------------------------------------------
    # construction helpers
    # ------------------------------------------------------------------
    def _validate(self):
        if not self.name:
            raise ValueError("an alphabet needs a name")
        if self.mode not in TOKENIZATION_MODES:
            raise ValueError(f"mode must be one of {TOKENIZATION_MODES}, got {self.mode!r}")
        if self.width < 1:
            raise ValueError("width must be at least 1")
        if not self.symbols:
            raise ValueError("an alphabet needs at least one symbol")

        everything = list(self.symbols) + list(self.extras)
        for label in everything:
            if not isinstance(label, str) or label == "":
                raise ValueError(f"symbols must be non-empty strings, got {label!r}")
            if not self.case_sensitive and label != label.lower():
                raise ValueError(
                    f"symbol {label!r} is not lower-case; use lower-case labels "
                    "or set case_sensitive=True")
        if len(set(everything)) != len(everything):
            dup = sorted({x for x in everything if everything.count(x) > 1})
            raise ValueError(f"duplicate symbols: {dup}")
        for label in self.extras:
            if len(label) != 1:
                raise ValueError(f"extra symbols must be single characters, got {label!r}")
        if self.mode == "chars":
            bad = [x for x in self.symbols if len(x) != 1]
            if bad:
                raise ValueError(f"mode 'chars' needs single-character symbols, got {bad[:3]}")
        if self.mode == "fixed":
            bad = [x for x in self.symbols if len(x) != self.width]
            if bad:
                raise ValueError(
                    f"mode 'fixed' with width {self.width} needs {self.width}-character "
                    f"symbols, got {bad[:3]}")
        if self.mode == "words":
            bad = [x for x in self.symbols if any(c.isspace() for c in x)]
            if bad:
                raise ValueError(f"mode 'words' needs symbols without whitespace, got {bad[:3]}")

        unknown_codes = [x for x in self.codes if x not in self.symbols]
        if unknown_codes:
            raise ValueError(f"codes given for symbols not in the alphabet: {unknown_codes[:3]}")
        if self.codes:
            values = list(self.codes.values())
            if len(set(values)) != len(values):
                raise ValueError("two symbols have the same nucleotide code")
            if len({len(v) for v in values}) != 1:
                raise ValueError("all nucleotide codes must have the same length")

        for src, dst in self.aliases.items():
            if self.mode == "fixed":
                # read character by character before the characters are grouped
                if len(src) != 1 or len(dst) != 1:
                    raise ValueError(
                        f"alias {src!r} -> {dst!r}: in mode 'fixed' an alias is one "
                        "character for another")
            elif dst not in everything:
                raise ValueError(f"alias {src!r} -> {dst!r}: {dst!r} is not a symbol")

    def _assign_characters(self):
        to_char = {}
        next_code = _PRIVATE_USE_START
        for label in itertools.chain(self.symbols, self.extras):
            if _is_identity_label(label):
                to_char[label] = label
            else:
                if next_code > _PRIVATE_USE_END:
                    raise ValueError(
                        "too many symbols for the private use area "
                        f"({_PRIVATE_USE_END - _PRIVATE_USE_START + 1} available)")
                to_char[label] = chr(next_code)
                next_code += 1
        object.__setattr__(self, "_to_char", to_char)
        object.__setattr__(self, "_to_label", {c: x for x, c in to_char.items()})

    # ------------------------------------------------------------------
    # the declared set
    # ------------------------------------------------------------------
    @property
    def labels(self) -> tuple:
        """All declared labels: the main symbols, then the extras."""
        return self.symbols + self.extras

    @property
    def size(self) -> int:
        """Number of declared symbols, extras included.

        This is the alphabet size every statistic should use, whether or not
        a given text contains all of them.
        """
        return len(self.symbols) + len(self.extras)

    @property
    def chars(self) -> str:
        """The internal characters of all declared symbols, in declared order."""
        return "".join(self._to_char[x] for x in self.labels)

    def __contains__(self, label: str) -> bool:
        return label in self._to_char

    def char(self, label: str) -> str:
        """Internal character of a label."""
        try:
            return self._to_char[label]
        except KeyError:
            raise KeyError(f"{label!r} is not a symbol of alphabet {self.name!r}") from None

    def label(self, char: str) -> str:
        """Label of an internal character."""
        try:
            return self._to_label[char]
        except KeyError:
            raise KeyError(f"{char!r} is not an internal character of alphabet {self.name!r}") from None

    # ------------------------------------------------------------------
    # internal string <-> labels
    # ------------------------------------------------------------------
    def encode(self, labels: Iterable[str]) -> str:
        """Labels -> internal string.  A label not in the alphabet raises KeyError."""
        return "".join(self.char(x) for x in labels)

    def decode(self, internal: str, sep: str = "", boundary: str = "|") -> str:
        """Internal string -> labels joined by ``sep``.

        A segment boundary (``"_"``) is shown as ``boundary`` so it cannot be
        mistaken for a symbol.
        """
        return sep.join(boundary if c == BOUNDARY else self.label(c) for c in internal)

    def decode_list(self, internal: str) -> list:
        """Internal string -> list of labels (a segment boundary stays ``"_"``)."""
        return [BOUNDARY if c == BOUNDARY else self.label(c) for c in internal]

    # ------------------------------------------------------------------
    # variations
    # ------------------------------------------------------------------
    def with_extras(self, extras: Iterable[str], name: Optional[str] = None) -> "Alphabet":
        """A copy with extra symbols (spaces, digits, ...) added."""
        merged = self.extras + tuple(x for x in dict.fromkeys(extras) if x not in self.extras)
        return self._replace(extras=merged,
                             name=name or (self.name if merged == self.extras else f"{self.name}+extras"))

    def _replace(self, **changes) -> "Alphabet":
        fields = dict(name=self.name, symbols=self.symbols, extras=self.extras,
                      codes=self.codes, mode=self.mode, width=self.width,
                      case_sensitive=self.case_sensitive, aliases=self.aliases,
                      defaults=self.defaults)
        fields.update(changes)
        return Alphabet(**fields)

    # ------------------------------------------------------------------
    # files
    # ------------------------------------------------------------------
    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "mode": self.mode,
            "width": self.width,
            "case_sensitive": self.case_sensitive,
            "symbols": [{"label": x, **({"code": self.codes[x]} if x in self.codes else {})}
                        for x in self.symbols],
            "extras": list(self.extras),
            "aliases": dict(self.aliases),
            "defaults": dict(self.defaults),
        }

    @classmethod
    def from_dict(cls, d: Mapping) -> "Alphabet":
        symbols, codes = [], {}
        for item in d["symbols"]:
            if isinstance(item, str):
                symbols.append(item)
            else:
                symbols.append(item["label"])
                if item.get("code"):
                    codes[item["label"]] = item["code"]
        return cls(name=d["name"], symbols=symbols, extras=d.get("extras", ()),
                   codes=codes, mode=d.get("mode", "chars"), width=d.get("width", 1),
                   case_sensitive=d.get("case_sensitive", False),
                   aliases=d.get("aliases", {}), defaults=d.get("defaults", {}))

    def save(self, path: str) -> str:
        """Write the alphabet as a JSON file and return the path."""
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2, ensure_ascii=False)
            fh.write("\n")
        return path

    @classmethod
    def load(cls, path: str) -> "Alphabet":
        """Read an alphabet written by :meth:`save`."""
        with open(path, "r", encoding="utf-8") as fh:
            return cls.from_dict(json.load(fh))

    @classmethod
    def from_codebook(cls, name: str, codebook: Mapping, extras: Iterable[str] = (),
                      mode: str = "chars", width: int = 1, **kwargs) -> "Alphabet":
        """Build an alphabet from a codebook ``{symbol: nucleotide code}``.

        Every symbol in the codebook is declared, whether or not a given text
        uses it.  See :func:`rna_lexis.codec.read_codebook_csv` for reading a
        codebook from a file.
        """
        return cls(name=name, symbols=tuple(codebook), extras=tuple(extras),
                   codes=dict(codebook), mode=mode, width=width, **kwargs)


# ----------------------------------------------------------------------
# built-in alphabets
# ----------------------------------------------------------------------
def nucleotide() -> Alphabet:
    """a, c, g, t; ``u`` is read as ``t``."""
    return Alphabet(name="nucleotide", symbols=tuple("acgt"), mode="chars",
                    aliases={"u": "t"})


def codon() -> Alphabet:
    """The 64 codons over a, c, g, t, read three characters at a time.

    The code of each codon is the codon itself, so conversion to nucleotide
    text is the identity.  ``u`` is read as ``t``.
    """
    codons = ["".join(p) for p in itertools.product("acgt", repeat=3)]
    return Alphabet(name="codon", symbols=tuple(codons), codes={x: x for x in codons},
                    mode="fixed", width=3, aliases={"u": "t"})


def latin() -> Alphabet:
    """The 26 lower-case Latin letters."""
    return Alphabet(name="Latin letters", symbols=tuple("abcdefghijklmnopqrstuvwxyz"),
                    mode="chars")


BUILTIN_ALPHABETS = {
    "nucleotide": nucleotide,
    "codon": codon,
    "latin": latin,
}


def get_alphabet(name_or_path: str) -> Alphabet:
    """A built-in alphabet by name, or an alphabet file by path."""
    key = name_or_path.strip().lower()
    if key in BUILTIN_ALPHABETS:
        return BUILTIN_ALPHABETS[key]()
    if os.path.isfile(name_or_path):
        return Alphabet.load(name_or_path)
    raise ValueError(
        f"{name_or_path!r} is neither a built-in alphabet "
        f"({', '.join(sorted(BUILTIN_ALPHABETS))}) nor an alphabet file")
