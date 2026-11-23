# Robustness, performance and security review

A review of the codon-units branch before manual testing: how the code treats text that
is wrong in some way, how fast it is on long texts, and what it should not trust. What was
checked, what was changed, and what is left and why.

## How it was checked

- **Probes.** Hand-written calls with empty, blank, wrong-typed, mismatched and hostile input
  at every layer: the tokenizer, the alphabet and its files, the statistics, the session and
  FASTA readers, the command line.
- **A fuzz run of the menu.** Every interactive handler (about 20, in nucleotide, codon and
  Latin modes) was driven with 14 seeded random sequences of garbage answers: blank, words,
  `-1`, `nan`, `1e999`, `999999999999`, `3.5`, unmatched brackets, `../../etc/passwd`, a NUL
  byte, 3,000 characters, non-Latin text. The first run found 79 distinct places that raised
  or hung, and a later one a place that used all the memory; after the fixes, 22 remain, all discussed below.
- **Profiling** with `cProfile` on 196,000 codons (177 human coding sequences) and a 41,000-letter
  English text, before and after.
- **A search** of the source for `eval`, `exec`, `pickle`, `marshal`, `shell=True`, `os.system`,
  `yaml`, every `subprocess` and network call, and every place a user-supplied string reaches a
  file name, a web address, a regular expression or a spreadsheet.
- **Before and after.** Output of the nucleotide mode (79 menu, plot and statistics outputs),
  of the codon and Latin discovery and statistics, and of the batch CLI is identical to the
  previous commit under a fixed hash seed. Without a fixed seed two of those outputs differ
  between runs of the *same* code, because a list is built from a set; that was already so.

## 1. Unknown letters, gaps and nothing at all

| Case | Before | Now |
|---|---|---|
| A symbol not in the alphabet | skipped with a warning (as designed) | same; the warning is now also an *alphabet-mismatch* warning when half or more of the text is unknown |
| `None`, a number, `bytes`, a list containing `None` as text | `TypeError: 'int' object is not iterable` | `TypeError` saying text must be a string or a list of strings, and which segment |
| Empty or blank text, an empty list | reads as nothing, no error | same; every statistic returns an empty result instead of failing |
| Text with no symbol of the alphabet | empty result, a count in a warning | the same, plus the mismatch warning and a menu/CLI message that nothing could be read |
| `width=0`, a `width` the symbols do not have, a bad `frame`/`gap`/`mode` | `width=0` silently meant the default; a wrong width read nothing | `ValueError` naming the problem |
| A line break between two words when the alphabet declares a space | the words were glued (`to beor`) | one space |
| Whitespace as a symbol in `"words"` mode | accepted, never produced | refused: whitespace is the separator |
| A motif or text longer than the other, a one-letter text | handled | handled (tests) |
| A FASTA file that starts with a `;` comment | **header lines read as sequence** | headers and comments are skipped (`io.split_fasta`) |

The last row was a real defect, in nucleotide mode too, found because the first codon
test file began with a comment line.

## 2. Alphabet that does not fit the data

`AlphabetMismatchWarning`, raised by `tokenize` (so by the menu, the CLI and the library alike):

- half or more of the text is not in the alphabet (and at least 5 units, or nothing was read):
  English as nucleotides, a protein as codons, hyphenated syllables as words. The message gives the
  percentage, the most common unknown symbols and a hint for the mode (frame and width for codons;
  "separated by whitespace" for words);
- the text is at least 95% `a c g t u n`, at least 30 units, and the alphabet has 20 or more
  symbols and reads one character per unit: DNA as Latin letters or amino acids.

It is a warning, not an error: a deliberate choice is not blocked. It does not fire for IUPAC-style
alphabets (15 symbols), for DNA as codons, for English as Latin letters, or for short texts. The
nucleotide path of the batch CLI, which drops non-nucleotides silently by design, now says on stderr when
most letters were dropped and stops when none were (the result files are unchanged).

## 3. Robustness to what users type and what files contain

- **Prompts.** Twenty-one `int(...)`/`float(...)` conversions of typed text, in 12 handlers, could raise and dump the user
  to the main menu. A value that is blank gives the default, as before; one that is not a finite number
  of the right kind in range now says so and uses the default. Upper limits where a huge value is
  dangerous (the unit-profile flank, the coverage exponent).
- **Parameters of the statistics** are checked at the top of each function (`alpha`, `mutr`, counts,
  Markov order, pseudocount), so a wrong value fails with its name instead of deep inside or, worse,
  returning meaningless numbers.
- **Alphabet files and sessions** raise or print a message that names the file and the problem,
  whatever they contain (empty, binary, not JSON, a list, a number, a nesting that exhausts the stack,
  fields of the wrong kind, symbols outside the alphabet).
- **A damaged session no longer loops.** The menu reloaded a session whose alphabet could not be read on
  every pass of its main loop, printing "Unexpected error" for ever.
- **The batch CLI** reports a bad option, a missing file or a malformed motif list as one line starting
  `error:` in every mode, where it was a traceback in nucleotide mode. A `--motif-column` that does not
  exist used to fall back silently to the first column.
- **The main loop's catch-all** now separates "something you typed or a file cannot be used"
  (`That could not be done: ...`) from a real bug (`Unexpected error (ValueError): ...`).

What remains of the fuzz results: a path typed for an output file that cannot be written (a directory,
a missing folder, a name that is too long, a NUL byte) raises `OSError` or `ValueError` inside the
handler, and the main loop reports it; the handlers were called directly, without that loop. A search
with a single common letter on a long text (thousands of occurrences) is heavy on purpose and passes
the 8-second limit of the harness.

## 4. Performance

Same machine, same data, same results (see above):

| | before | after |
|---|--:|--:|
| `find_boundary`, 196,000 codons | 35 s | 0.8 s |
| `rank_core_candidates`, 196,000 codons, 300 xmotifs | 47 s | 1.1 s |
| `rank_core_candidates`, 41,000 letters, 889 xmotifs | 17 s | 0.7 s |
| `score_exact_motifs`, 20 motifs, 196,000 codons | 1.9 s | 0.08 s |
| `mutation_family_tests`, 588,000 nt, 2 motifs | 0.9 s | 0.18 s |

Where the time went, in order of how much it mattered:

1. `is_bounded` lower-cased the whole text on every call, and `find_boundary` called it tens of thousands
   of times with text it had already lower-cased (30 of 35 s).
2. `find_boundary` scanned the whole text at every word length, though a word of length `j` can only recur
   if its first `j-1` letters did. It now looks only where the shorter word recurred, and reuses the
   positions of the inner word that it already has.
3. `normalize_sequence` checked the entire text on every motif lookup (44 of 51 s in `rank_core_candidates`).
   It remembers the last few long texts, by identity, and uses set operations.
4. "Is this occurrence inside an xmotif" compared against every interval (422 million comparisons on the
   41,000-letter text); it is now one binary search.
5. The Hamming scan looped over windows in Python; it counts mismatches for all windows at once with NumPy.
6. `find_all_matches` lower-cased the text once per match.

Each replacement has a test that compares it with the plain definition it replaces on random inputs
(`tests/test_performance_paths.py`), and each test was checked to fail when the code is made wrong.

Not changed, and why: `find_longest_extensions` extends every pair of occurrences, so its time grows with the
square of the occurrences times the length of the extension (16 s for a 4-letter word occurring 238 times
in the doubled 41,000-letter text, which is an extreme repeat). Every pair is the definition of the result, so
the menu asks first when there are more than 300,000 pairs. `MarkovBackground` and the tokenizer are linear and
take under 0.5 s for 600,000 characters.

## 5. Security

**Looked for and not found:** `eval`, `exec`, `compile`, `pickle`, `marshal`, `shell=True`, `os.system`, YAML
loading, dynamic imports of user-named modules, `tempfile.mktemp`. Every `subprocess` call takes a list. The
only network calls are HTTPS to `rest.ensembl.org` and `api.genome.ucsc.edu`, with timeouts and the default
certificate checks.

**Fixed:**

| Risk | Fix |
|---|---|
| An identifier is put into a web address: `ENST1/../../x`, `ENST1?type=x`, a chromosome or genome with `;` could change which address is fetched | Identifiers are matched against what they can be (`ENST` digits and one version, `chr` and letters or digits, an assembly name, a cCRE accession); a range must be whole numbers, `0 <= start < end`, and at most 10 million bases |
| An answer from a server of any size is read into memory | Capped at 100 MB |
| A text cell such as `=HYPERLINK(...)` in a results table runs when the table is opened in Excel or LibreOffice (an alphabet file can name its symbols anything) | A leading apostrophe on any text cell starting with `=`, `+`, `-`, `@` that is not a number: the CSV writers and the summary table. Nucleotide output has only letters and is unchanged |
| A hostile or accidental huge or deeply nested JSON file used up memory or the stack | Alphabet files over 5 MB and sessions over 200 MB are refused before reading; `RecursionError` is caught |
| A session that is a list, or has text outside its own alphabet, crashed the menu | Checked on loading (`io.session_problem`) |
| A flank of 999999999999 in the unit profile allocated a table that size and the process was killed | Capped at 500 units |

**Left, with reasons:**

1. *Regular expressions typed as queries (nucleotide mode).* The search prompts treat what is typed as a
   regular expression, with `.` as the documented wildcard. A pattern with nested repeats such as
   `(a+)+b` on a long text can run for a very long time (catastrophic backtracking). The only person
   affected is the one typing it, in their own session, and `re` offers no time limit; a syntax error is
   now caught. Over an alphabet nothing but letters and `.` reaches the pattern. Closing this would mean
   replacing regular expressions with a purpose-built wildcard matcher, which changes nucleotide behavior,
   so it was not done.
2. *A session file names where results are written* (`dir`, `file_path`), and the menu changes to that
   folder. A session from someone else could point results at a folder where a file of the same name
   would be replaced. Names are fixed suffixes, nothing is run, and the user guide now says to open sessions
   only from people you trust. A confirmation when `dir` is not the folder the file is in would be a small
   addition.
3. *The summary table has a last-resort location in the shared temporary folder* with a predictable name
   (`init_summary`), used only when the session folder, the working directory and the home folder all
   fail. On a shared machine another user could place a link there first. Removing the last candidate
   would turn a rare success into a failure, so it was left.
4. *Dependencies were not audited* (no `pip-audit` or advisory database was available). `kaleido` is
   pinned below 1 on purpose (see `CLAUDE.md`); that series bundles an older browser engine, which here
   renders only figures the package has built itself.
5. *HTML figures contain labels taken from files* (FASTA headers, alphabet symbols). Plotly renders its text
   itself and allows only a few harmless tags, but this was not tried in a browser.
6. The helper that raises a LibreOffice window passes a file name to `wmctrl` and `xdotool` as an argument (a
   list, so no shell); a name beginning with `-` would be read as an option. No effect beyond that window.

## Re-running the checks

The unit tests (`python3 -m pytest -q`) include everything that can run without a terminal. The
menu fuzz harness and the golden-output comparisons were scripts kept outside the repository; the method
is above, and the harness is about 120 lines: replace `safe_input` with a function that returns a random
element of a list of bad answers and raises `EOFSignal` after 60 prompts, set a per-call alarm, call each
`*_input` handler, and report the exception type and the line.
