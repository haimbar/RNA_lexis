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

MAX_ALPHABET_FILE_BYTES = 5 * 2**20
"""An alphabet file larger than this is refused: no real alphabet is, and a huge JSON file would only use up memory."""


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
        strip_marks:    When True, accents, breathings and other combining marks are removed
                        from the text before it is read (after lower-casing, before the aliases),
                        so a polytonic Greek text reads over the 24 plain letters: ``ἄνδρα``
                        is ``ανδρα``.  Default False.
        aliases:        Input read as something else, e.g. ``u`` -> ``t``.  In
                        ``"chars"`` and ``"words"`` mode the target must be a
                        symbol; in ``"fixed"`` mode an alias is one character
                        for another, applied before characters are grouped.
        defaults:       Analysis defaults that belong to this alphabet
                        (background order, minimum motif length, ...).  Kept
                        as given; the analysis code decides how to use them.
        unit_name:      What one symbol is called in axis labels and titles
                        (``"nt"``, ``"codon"``, ``"letter"``, ...).
        display_sep:    What goes between symbols when a motif is shown.
                        ``None`` (default) means nothing for single-character
                        symbols and a space for longer ones, so codons read as
                        ``atg aaa ttt``.
    """

    name: str
    symbols: tuple
    extras: tuple = ()
    codes: Mapping = field(default_factory=dict)
    mode: str = "chars"
    width: int = 1
    case_sensitive: bool = False
    strip_marks: bool = False
    aliases: Mapping = field(default_factory=dict)
    defaults: Mapping = field(default_factory=dict)
    unit_name: str = "unit"
    display_sep: Optional[str] = None

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
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("an alphabet needs a name (a non-empty string)")
        if self.mode not in TOKENIZATION_MODES:
            raise ValueError(f"mode must be one of {TOKENIZATION_MODES}, got {self.mode!r}")
        if isinstance(self.width, bool) or not isinstance(self.width, int) or self.width < 1:
            raise ValueError(f"width must be a whole number, at least 1, got {self.width!r}")
        if not isinstance(self.case_sensitive, bool):
            raise ValueError(f"case_sensitive must be true or false, got {self.case_sensitive!r}")
        if not isinstance(self.strip_marks, bool):
            raise ValueError(f"strip_marks must be true or false, got {self.strip_marks!r}")
        if not self.symbols:
            raise ValueError("an alphabet needs at least one symbol")
        for key, value in self.aliases.items():
            if not isinstance(key, str) or not isinstance(value, str) or not key:
                raise ValueError(f"an alias maps one string to another, got {key!r}: {value!r}")
        for key, value in self.defaults.items():
            if key in _SETTING_KEYS and (isinstance(value, bool) or not isinstance(value, int) or value < 1):
                raise ValueError(f"default {key!r} must be a whole number, at least 1, got {value!r}")
        if not isinstance(self.unit_name, str):
            raise ValueError(f"unit_name must be a string, got {self.unit_name!r}")
        if self.display_sep is not None and not isinstance(self.display_sep, str):
            raise ValueError(f"display_sep must be a string or null, got {self.display_sep!r}")
        for label, code in self.codes.items():
            if not isinstance(code, str) or not code:
                raise ValueError(f"the nucleotide code of {label!r} must be a non-empty string, got {code!r}")

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
            bad = [x for x in self.extras if x.isspace()]
            if bad:
                raise ValueError(
                    "mode 'words' reads whitespace as the separator between symbols, so whitespace "
                    f"cannot also be a symbol (got {bad[:3]}); write the separator into the text as a symbol "
                    "such as '|' instead")

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
    # display
    # ------------------------------------------------------------------
    @property
    def sep(self) -> str:
        """What goes between symbols when a motif is shown."""
        if self.display_sep is not None:
            return self.display_sep
        return "" if self.mode == "chars" else " "

    @property
    def label_width(self) -> int:
        """Length of the longest symbol label."""
        return max(len(x) for x in self.labels)

    @property
    def units_label(self) -> str:
        """Plural of ``unit_name`` for axis labels (``nt`` stays ``nt``)."""
        return self.unit_name if self.unit_name == "nt" else f"{self.unit_name}s"

    def show(self, internal: str, mask=None, sep: Optional[str] = None, upper: bool = False,
             pad: bool = False, boundary: str = "|") -> str:
        """An internal string as text a person can read.

        Args:
            internal: Internal string (or a slice of one).
            mask:     Optional sequence of bools, one per symbol; the flagged
                      symbols are shown in upper case ("this one differs from
                      the reference").  For a case-sensitive alphabet, where
                      upper case would change the symbol, they are shown in
                      square brackets instead.
            sep:      Between symbols; default :attr:`sep`.
            upper:    Show every symbol in upper case (the convention for
                      motif columns in result files).
            pad:      Pad every label to the longest label so columns line up.
            boundary: How a segment boundary is shown.
        """
        sep = self.sep if sep is None else sep
        if mask is not None and len(mask) != len(internal):
            raise ValueError("mask must have one entry per symbol")
        width = self.label_width
        parts = []
        for i, c in enumerate(internal):
            if c == BOUNDARY:
                parts.append(boundary)
                continue
            label = self.label(c)
            flagged = upper or (mask is not None and bool(mask[i]))
            if flagged:
                label = f"[{label}]" if self.case_sensitive else label.upper()
            parts.append(label.ljust(width) if pad else label)
        return sep.join(parts)

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
                      defaults=self.defaults, unit_name=self.unit_name,
                      display_sep=self.display_sep)
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
            "strip_marks": self.strip_marks,
            "symbols": [{"label": x, **({"code": self.codes[x]} if x in self.codes else {})}
                        for x in self.symbols],
            "extras": list(self.extras),
            "aliases": dict(self.aliases),
            "defaults": dict(self.defaults),
            "unit_name": self.unit_name,
            "display_sep": self.display_sep,
        }

    @classmethod
    def from_dict(cls, d: Mapping) -> "Alphabet":
        """Build an alphabet from the dict that :meth:`to_dict` writes.

        A description that is not valid (a missing field, a field of the wrong kind) raises ``ValueError``
        saying what is wrong, whatever the file contained.
        """
        if not isinstance(d, Mapping):
            raise ValueError(f"an alphabet description must be a JSON object, got {type(d).__name__}")
        for key in ("name", "symbols"):
            if key not in d:
                raise ValueError(f"the alphabet description has no {key!r}")
        if not isinstance(d["symbols"], (list, tuple)):
            raise ValueError("'symbols' must be a list")
        for key, kind in (("extras", (list, tuple, str)), ("aliases", Mapping), ("defaults", Mapping)):
            if d.get(key) is not None and not isinstance(d[key], kind):
                raise ValueError(f"{key!r} must be a {'list' if kind != Mapping else 'JSON object'}")
        symbols, codes = [], {}
        for item in d["symbols"]:
            if isinstance(item, str):
                symbols.append(item)
            elif isinstance(item, Mapping) and isinstance(item.get("label"), str):
                symbols.append(item["label"])
                if item.get("code"):
                    codes[item["label"]] = item["code"]
            else:
                raise ValueError(f"each symbol must be a string or an object with a 'label', got {item!r}")
        defaults = d.get("defaults") or {}
        unit_name = d.get("unit_name", "unit")
        if not defaults:
            # a file (or a session) from before the defaults were calibrated: use the built-in's
            builtin = _BUILTIN_BY_NAME.get(d["name"])
            if builtin is not None:
                built = builtin()
                if symbols == list(built.symbols):
                    defaults = dict(built.defaults)
                    if unit_name == "unit":
                        unit_name = built.unit_name
        try:
            return cls(name=d["name"], symbols=symbols, extras=d.get("extras") or (),
                       codes=codes, mode=d.get("mode", "chars"), width=d.get("width", 1),
                       case_sensitive=d.get("case_sensitive", False),
                       strip_marks=d.get("strip_marks", False),
                       aliases=d.get("aliases") or {}, defaults=defaults,
                       unit_name=unit_name, display_sep=d.get("display_sep"))
        except (TypeError, AttributeError) as exc:           # a field of the wrong kind that no check above names
            raise ValueError(f"the alphabet description is not valid: {exc}") from exc

    def same_symbols(self, other: "Alphabet") -> bool:
        """True if ``other`` reads text the same way: same name, symbols, extras, mode, width, case
        and aliases.  Analysis defaults, the unit name and the display separator may differ, so a
        session saved before a default was calibrated still matches the current built-in alphabet."""
        keys = ("name", "symbols", "extras", "mode", "width", "case_sensitive", "strip_marks", "aliases")
        mine, theirs = self.to_dict(), other.to_dict()
        return all(mine[k] == theirs[k] for k in keys if k not in ("symbols",)) and \
            [s["label"] for s in mine["symbols"]] == [s["label"] for s in theirs["symbols"]]

    def save(self, path: str) -> str:
        """Write the alphabet as a JSON file and return the path."""
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(self.to_dict(), fh, indent=2, ensure_ascii=False)
            fh.write("\n")
        return path

    @classmethod
    def load(cls, path: str) -> "Alphabet":
        """Read an alphabet written by :meth:`save`.

        Raises ``ValueError`` for a file that is too large to be an alphabet (over 5 MB), is not JSON, or
        does not describe an alphabet; the message names the file.
        """
        try:
            if os.path.getsize(path) > MAX_ALPHABET_FILE_BYTES:
                raise ValueError(f"{path}: too large to be an alphabet file "
                                 f"(over {MAX_ALPHABET_FILE_BYTES // 2**20} MB)")
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
        except (json.JSONDecodeError, UnicodeDecodeError, RecursionError) as exc:
            raise ValueError(f"{path}: not a JSON file ({exc})") from exc
        try:
            return cls.from_dict(data)
        except ValueError as exc:
            raise ValueError(f"{path}: {exc}") from exc

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
# projections: many-to-one maps between alphabets
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class Projection:
    """A many-to-one map from one alphabet onto another.

    Translating codons to amino acids is one; reading a methylated cytosine
    (``m``) as an ordinary one (``c``) is another.  Analyzing the projected text
    asks the question with some distinctions removed.

    Attributes:
        source:  The alphabet of the text to project.
        target:  The alphabet of the result.
        mapping: ``{source label: target label}``.  Every main symbol of the
                 source must be mapped; an extra symbol that is not mapped is
                 kept if the target declares it too, and is an error otherwise.
        name:    Display name.
    """

    source: Alphabet
    target: Alphabet
    mapping: Mapping
    name: str = "projection"

    def __post_init__(self):
        table = {}
        for label in self.source.labels:
            if label in self.mapping:
                dst = self.mapping[label]
            elif label in self.source.extras and label in self.target:
                dst = label
            else:
                raise ValueError(f"symbol {label!r} of {self.source.name!r} is not mapped")
            if dst not in self.target:
                raise ValueError(f"{label!r} maps to {dst!r}, which is not a symbol of {self.target.name!r}")
            table[self.source.char(label)] = self.target.char(dst)
        table[BOUNDARY] = BOUNDARY
        object.__setattr__(self, "table", table)

    def apply(self, internal: str) -> str:
        """Project an internal string of the source alphabet onto the target alphabet."""
        try:
            return "".join(self.table[c] for c in internal)
        except KeyError as exc:
            raise ValueError(f"{exc.args[0]!r} is not an internal character of {self.source.name!r}") from None

    def groups(self) -> dict:
        """``{target label: [source labels that map to it]}``."""
        out: dict = {}
        for src, dst in self.mapping.items():
            out.setdefault(dst, []).append(src)
        return out


# ----------------------------------------------------------------------
# built-in alphabets
# ----------------------------------------------------------------------
def nucleotide() -> Alphabet:
    """a, c, g, t; ``u`` is read as ``t``."""
    return Alphabet(name="nucleotide", symbols=tuple("acgt"), mode="chars",
                    aliases={"u": "t"}, unit_name="nt", defaults=dict(NUCLEOTIDE_SETTINGS))


def codon() -> Alphabet:
    """The 64 codons over a, c, g, t, read three characters at a time.

    The code of each codon is the codon itself, so conversion to nucleotide
    text is the identity.  ``u`` is read as ``t``.
    """
    codons = ["".join(p) for p in itertools.product("acgt", repeat=3)]
    return Alphabet(name="codon", symbols=tuple(codons), codes={x: x for x in codons},
                    mode="fixed", width=3, aliases={"u": "t"}, unit_name="codon",
                    defaults={"minxmlen": 5, "mincorelen": 3, "candidate_min_len": 3, "mincount": 2})


def latin() -> Alphabet:
    """The 26 lower-case Latin letters."""
    return Alphabet(name="Latin letters", symbols=tuple("abcdefghijklmnopqrstuvwxyz"),
                    mode="chars", unit_name="letter",
                    defaults={"minxmlen": 6, "mincorelen": 4, "candidate_min_len": 4, "mincount": 2})


_BUILTIN_BY_NAME = {
    "nucleotide": nucleotide,
    "codon": codon,
    "Latin letters": latin,
}


# Settings of the nucleotide pipeline, as they have always been.
NUCLEOTIDE_SETTINGS = {"minxmlen": 7, "mincorelen": 6, "candidate_min_len": 5, "mincount": 2}
_SETTING_KEYS = tuple(NUCLEOTIDE_SETTINGS)


def default_settings(alphabet=None) -> dict:
    """The discovery settings an alphabet starts with.

    ``minxmlen`` (shortest xmotif), ``mincorelen`` (shortest core, for the cores found from the
    xmotifs), ``candidate_min_len`` (shortest core candidate in *Rank core motifs*) and ``mincount``
    (fewest occurrences), all in units.  No alphabet is nucleotide text, which keeps the settings
    it has always had.  Built-in alphabets carry calibrated values (``calibration/README.md``); an
    alphabet file may give its own in ``defaults``; any other alphabet gets a rule: nucleotide
    settings for a small alphabet (8 symbols or fewer), and 5, 3, 3, 2 otherwise.  The shortest
    core to look for is the thing to choose: the shortest xmotif is that plus 2.
    """
    if alphabet is None:
        return dict(NUCLEOTIDE_SETTINGS)
    out = dict(NUCLEOTIDE_SETTINGS) if alphabet.size <= 8 else {
        "minxmlen": 5, "mincorelen": 3, "candidate_min_len": 3, "mincount": 2}
    out.update({k: v for k, v in alphabet.defaults.items() if k in _SETTING_KEYS})
    return out


def _amino_acid() -> Alphabet:
    from rna_lexis.genetic_code import amino_acid      # imported here: genetic_code imports this module
    return amino_acid()


BUILTIN_ALPHABETS = {
    "nucleotide": nucleotide,
    "codon": codon,
    "latin": latin,
    "amino-acid": _amino_acid,
}


_NAME_ALIASES = {"latin letters": "latin", "letters": "latin", "amino acid": "amino-acid",
                 "amino acids": "amino-acid", "aminoacid": "amino-acid", "amino_acid": "amino-acid",
                 "nucleotides": "nucleotide", "codons": "codon"}


def get_alphabet(name_or_path: str) -> Alphabet:
    """A built-in alphabet by name (also by its display name, e.g. ``"Latin letters"``), or an alphabet file by path."""
    if not isinstance(name_or_path, (str, os.PathLike)):
        raise ValueError(f"an alphabet is named by a string or the path of an alphabet file, got {name_or_path!r}")
    name_or_path = os.fspath(name_or_path)
    key = " ".join(name_or_path.strip().lower().split())
    key = _NAME_ALIASES.get(key, key)
    if key in BUILTIN_ALPHABETS:
        return BUILTIN_ALPHABETS[key]()
    path = os.path.expanduser(name_or_path.strip())
    if os.path.isfile(path):
        return Alphabet.load(path)
    raise ValueError(
        f"{name_or_path!r} is neither a built-in alphabet "
        f"({', '.join(sorted(BUILTIN_ALPHABETS))}) nor an alphabet file")
