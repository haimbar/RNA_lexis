# Manual test plan: whole-unit analysis (branch `codon-units-alphabet`)

You test by hand; this plan says what to install, what to try, what you should see, and what to
keep so that I can reproduce anything you report. Tick the boxes as you go. Nothing here merges to
`main`; the branch stays as it is until you say you are done.

Commit under test: `d502fb2d38ca343c871653a8c74623c35db9135c` (change the hash below if we fix
something and you reinstall).

---

## 1. Install from the branch in a fresh virtual environment

A separate environment keeps this version away from any RNA_lexis you already use.

```bash
python3 -m venv ~/venvs/lexis-branch
source ~/venvs/lexis-branch/bin/activate
python -m pip install --upgrade pip
pip install "git+https://github.com/haimbar/RNA_lexis.git@d502fb2d38ca343c871653a8c74623c35db9135c"
```

Pinning the hash, not the branch name, means the install is exactly what I tested and can be
installed again later. After a fix, reinstall the new hash over it:

```bash
pip install --force-reinstall --no-deps "git+https://github.com/haimbar/RNA_lexis.git@<new hash>"
```

Check it:

```bash
python -c "import rna_lexis, sys; print(rna_lexis.__file__)"     # must be inside ~/venvs/lexis-branch
python -c "import json, importlib.metadata as m; print(json.loads(m.distribution('RNA_lexis').read_text('direct_url.json'))['vcs_info']['commit_id'])"
rna_lexis_stat_cli --help
python -c "import tkinter; print('tkinter ok')"            # needed for the file-picker windows
which weblogo                                               # nucleotide logo only; optional
```

- The banner in the menu shows the package version (0.2.7), which does not change between commits, so
  it cannot tell you which commit you have. Use the `commit_id` line above.
- The menu keeps its folder preference in `~/Library/Application Support/RNA_lexis/prefs.json`. A venv does
  not separate that file from other installs.
- Start the menu with `PYTHONHASHSEED=0 rna_lexis`. Some lists are built from sets, so their order can differ
  between runs; fixing the seed makes a run repeat exactly.

## 2. Set up the folders (outside the repository)

```
~/lexis-testing/
  inputs/        your original files, never edited; one list of SHA-256 sums (see section 3)
  alphabets/     alphabet JSON files you write
  runs/          one folder per test, named like  T12-antony-latin/
  notes.md       the log (template in section 6)
```

Nothing creates these folders for you. Make them once with
`mkdir -p ~/lexis-testing/{inputs,alphabets,runs}` (`controls/` is already there, with ready-made alphabet files),
and make each run folder before you start the test, for example
`mkdir ~/lexis-testing/runs/T16-synonymous && cp ~/lexis-testing/controls/planted/synonymous_demo.fa ~/lexis-testing/runs/T16-synonymous/`.
`script` will not create a missing folder for its log.

The menu writes its session file, tables and plots **next to the input and changes into that folder**. So
for each test, make a run folder and put a *copy* of the input in it before loading. Then everything one
test produced is in one place, and the originals stay clean.

**Careful with main menu item 8, *Clear workspace*.** It deletes every file in the working folder except
the loaded input and its session file, including your tables, plots and any other inputs there. It asks you
to type `YES`. Use it only in a run folder you are done with. It keeps the loaded input, its session file
and any other FASTA file (`.fasta`, `.fa`, `.fna`, `.fas`); check that in T73.

## 3. What to save so I can reproduce a result

| Save | How | Why |
|---|---|---|
| The exact environment | once per sitting: `python --version; pip freeze > runs/env.txt; sw_vers >> runs/env.txt` and the commit line from section 1 | same code, same libraries |
| The input files | keep copies under `inputs/`, and `shasum -a 256 inputs/* > inputs/SHA256.txt` | I can check I use the same bytes |
| Anything fetched from the web (ENST, genomic range, cCRE) | keep the session `.json` the menu writes; note the date | the databases change; the session holds the sequence you saw |
| **What you typed** | start the menu through `script`: `script -q runs/T12-antony-latin/terminal.log rna_lexis` | the log has every prompt, your answer and every warning, so I can replay it with the same answers. The file-picker windows are not in the log, so also write the file you picked in the notes |
| The session file | the `.json` the menu writes beside the input; do not rename it | holds the text as read, the alphabet, the reading frame, the cores and xmotifs, the settings |
| Alphabet files | in `alphabets/`, one file per alphabet, with a note of what it is for | needed to read the text |
| Settings changed | write them in the notes (*Show settings* before and after) | defaults depend on the alphabet |
| Outputs | leave generated tables and plots with their generated names (they contain the motif) | I compare numbers and look at figures |
| Batch CLI | the exact command lines, one per line, in `runs/<test>/commands.sh` | repeatable as is |
| Anything odd | a screenshot, and the terminal log | words fail where a picture does not |

You do not have to save anything for the passing checks except a tick. Save all of the above for
anything that looks wrong, surprising or interesting.

What I can do with it: I can read files anywhere under your home folder, install the same commit in a
scratch environment, load your session `.json` files, run your CLI lines, and replay the answers from a
`terminal.log` through the same handlers (the menu has no terminal in my tests, so I script the answers).
I cannot use the file-picker windows or look at a plot that is only on your screen, so save plots to files.

## 4. Inputs

**Yours**, and the reading to try first with each:

| Input | Read as | Also try |
|---|---|---|
| Genes: coding sequences (CDS) | codons, frame 0 | codons translated to amino acids; frames 1 and 2 (should give a different, mostly meaningless result) |
| Genes: full transcripts (ENST fetch) | nucleotides (automatic) | codons from a CDS only; transcripts include UTRs, so codon reading of a transcript is wrong unless you pick the frame |
| Mark Antony (English) | Latin letters | an alphabet file with a space declared (letters plus `" "`), to test word gaps |
| Whale and finch recordings | an alphabet file in `"words"` mode, one symbol per call or syllable | each bout or song as its own FASTA record, to see whether motifs span bouts |
| The other texts (Constitution, Genesis, Homer in Greek) | Latin letters; Greek letters from an alphabet file | compare two documents with the shared-motif diagram |

**Ready to use in the repository:** `calibration/data/human_cds.fasta` (177 human coding sequences, 588,000 nt; a
good size for timing, and its first line is a `;` comment, which is a test in itself), and the bundled NORAD
example (*Load example dataset*).

**A worked case that needs no file** (*Paste sequence*): a coding sequence in which the same seven-residue
peptide, Met-Lys-Val-Leu-Ala-Gly-Ile, appears four times, each time spelled with different synonymous codons
and surrounded by different filler.

```
attcatacggtgcgcctgagcaccacaaaccgcgatatgaaagttcttgccggtattttc
cagtcctacaatatggacgagcctcactatagtatgaaggtcctcgcgggcatctcgacg
caaatctttttgaactctggctcccaatccatgaaagtactagctggaatacacataagg
catcagggatgggtgtgctgctgtcctatgaaggtgctggcagggattacctgcgcgccc
gggcgagaaattagtattctgcac
```

Read as **codons** it contains no repeated word (the longest shared stretch is 2 codons). Read as **amino
acids** it contains `MKVLAGI` four times. If both readings agree, something is wrong.

**The controls pack** is in `~/lexis-testing/controls/` (made from your own Antony, whale, finch and gene files;
its `README.md` lists every file with what to expect, as observed on this branch). In short:
- `shuffled/`: the same material with the order destroyed. Expect almost no supported cores (Antony letters 0
  against 482 in the original; the whale shuffle gave 2; human genes 0).
- `planted/`: repeats whose number and place are known (`TRUTH.json`), a core split across two FASTA records,
  the synonymous-codon case as a file, and syllable files for an alphabet in words mode.
- `alphabets/`, `sessions/`, `hard_cases/`: good and deliberately bad alphabet files, session files and inputs
  (empty, binary, Latin-1, a `;` comment before the header, CRLF, `N` codons, a protein, 300 kb and 3 Mb of
  random nucleotides for timing).

## 5. The checks

Menu names are as on screen. "Expect" is what the branch is meant to do.

### A. Install and smoke (10 minutes)
- [ ] **T01** Section 1 checks all pass; `rna_lexis_stat_cli --help` lists four subcommands.
- [ ] **T02** `rna_lexis` starts, shows the input menu; *Load example dataset (NORAD)*, choose a save folder, reading *Automatic*. Expect: nucleotides, main menu shows the file name.
- [ ] **T03** Main menu *Plots* > *Core neighbors (detailed)* on NORAD, save as PNG. Expect a figure and no traceback.

### B. How the text is read (the reading menu after a new input)
- [ ] **T10** Automatic on DNA / RNA: nucleotides. On English or a protein: Latin letters (the message says so).
- [ ] **T11** *Codons*, frame 0, 1, 2 on one CDS. Expect: the number of codons read changes with the frame; a partial codon at the end is dropped with a warning; a codon with `N` is skipped as a whole (the frame does not shift).
- [ ] **T12** *Codons translated to amino acids*. Expect: stops shown as `*`; queries are typed as amino-acid letters; the settings line says the reading was translated.
- [ ] **T13** *Amino acids* on a protein FASTA; a letter such as `X` is skipped with a warning.
- [ ] **T14** *An alphabet from a file*: Greek letters (`~/Desktop/Homer`, with `controls/alphabets/greek.json`, which folds accents and the final sigma). Loading a text with no Latin letters must go straight to this option, with a note that only an alphabet from a file can read it; a `"words"` alphabet for syllables; one with `"extras": [" "]`. Expect: the declared size is used even for symbols that never occur.
- [ ] **T15** *Cancel* returns to the input menu.
- [ ] **T16** Synonymous-codon case (`controls/planted/synonymous_demo.fa`, or the pasted text above). As codons: no xmotif and no core, and *Find all matches* with the nucleotides `atgaaagttcttgccggtatt` finds 1. As amino acids: one xmotif and core, `MKVLAGI`, and *Find all matches* with `mkvlagi` finds 4.

### C. Wrong alphabet, unknown symbols, odd input
- [ ] **T20** Read English as *Codons*, and DNA as *Latin letters* or *Amino acids*. Expect a yellow warning that the alphabet may be wrong (the DNA one says it "looks like DNA or RNA"); the reading still goes ahead.
- [ ] **T21** Text with digits, punctuation, accents, emoji. Expect: skipped with a count and the commonest symbols named; no crash.
- [ ] **T22** Text with a declared space and hard line breaks: `to be` newline `or not` is read as `to be or not`.
- [ ] **T23** A FASTA file that begins with a `;` comment line, and a bare `>` header. Expect the letters of the headers are not in the sequence (compare the length with the file size).
- [ ] **T24** Empty file, binary file, Latin-1 file. Expect a message and a return to the menu, never a traceback or a hang.
- [ ] **T25** At every prompt: blank (default), `abc`, `-1`, `0`, `3.5`, `1e999`, `nan`, `999999999999`, a very long string, Ctrl+D. Expect: a red "`'abc' is not a whole number; using the default, 5`" or similar, never a traceback; Ctrl+D returns to the main menu.
- [ ] **T26** A pattern that is not valid at *Find all matches* in nucleotide mode (`[a`, `(`, `a{2,1}`): reported, asked again.
- [ ] **T27** A damaged session `.json` (delete a bracket; put a letter outside the alphabet in `txt`). Expect: "Cannot use the session ...: <reason>", and the menu goes on, not an endless loop.

### D. Settings per alphabet
- [ ] **T30** *Show settings* after each reading. Expect these starting values (min xmotif, min core, rank-core candidate, min occurrences): nucleotides 7, 6, 5, 2; codons and amino acids 5, 3, 3, 2; Latin letters 6, 4, 4, 2; another alphabet of more than 8 symbols 5, 3, 3, 2.
- [ ] **T31** *Change setting* then load text with a different alphabet. Expect: a setting you changed stays; one you did not follows the new alphabet.
- [ ] **T32** Load a new input (main menu 9) several times with different readings in one sitting. Nothing from the previous text should leak.

### E. Plots over each alphabet (nucleotides, codons, Latin letters)
For each item, run it on nucleotides first (it should behave as it always did), then on codons and letters.
- [ ] **T50** *Core neighbors (detailed)* and *(condensed)*: labels are codons or letters, positions are in units.
- [ ] **T51** *K-mers*: z-score, robust z-score, histogram, rank-frequency. A k so large that nothing repeats gives a message, not an error.
- [ ] **T52** *Logo (nucleotides) / unit profile*: nucleotides give the logo; other alphabets give the **unit profile**. Expect a **letter logo**: at each position next to the motif the letters that are clearly more common there than in the whole text, as tall as their enrichment, with a count above each stack (`3/53`); the motif as one grey block; vowels red and consonants blue (amino acids and codons colored by chemistry); the same letters listed with counts in the terminal. Choose style 2 for the older heatmap. Check a rare letter: it is small, and its count tells you how many occurrences it stands for. Try a flank of 0, 5, 50, 501 (the last is refused).
- [ ] **T53** *Coverage*, *Motif Match/Mutation*, *Self-similarity arc plot*, *Shared-motif diagram* (load a second text over the same alphabet; a session read over a different alphabet is refused), *Compare with another transcript*, *Compare motif coverage vs. homopolymer density*.
- [ ] **T54** *Compare motif coverage vs. GC-content* and the sequence logo on a non-nucleotide reading: the menu says it needs nucleotides and returns.
- [ ] **T55** Save each as PNG, PDF, SVG, HTML where offered; open the file.

### F. Sequence operations
- [ ] **T60** *Find all matches* (`.` is any one unit; for codons type nucleotides: `atgaaa.` is three codons), *Search with mutations*, *Motif extensions*, *Print core*, *Covered area*.
- [ ] **T61** *Motif extensions* and *Self-similarity arc plot* with a very common single letter (`e`, `a`) in a long text. Expect: it asks before comparing millions of pairs.
- [ ] **T62** *Motif spacing / periodicity test*, *Gapped motif search*, *Batch spacing test*. Anchors that are not in the alphabet are dropped with a warning.
- [ ] **T63** *Rank core motifs (Markov/FDR)* and *Mutation-family scoring*: open the CSV. Motifs are shown as labels (`ATG AAA TTT`). The output name is yours or the default; a path that cannot be written is reported.
- [ ] **T64** *Core neighbors (text export)*, *K-mer Markov analysis*, *Alignment score for two sequences* (nucleotides only: BLAST bit score and E-value are not offered for other alphabets).
- [ ] **T65** *Export hairpins* on nucleotides only; on others the menu says so.

### G. Sessions
- [ ] **T70** Quit and reload a saved session. Expect: same alphabet and reading, same cores, same tables; the main menu shows the reading.
- [ ] **T71** Open a session saved before this branch (any old `.json`). Expect: it loads as nucleotide text.
- [ ] **T72** *Summary statistics* twice: the second time it should not recompute (the table is only rewritten when the inputs change).
- [ ] **T73** *Clear workspace*, in a throwaway run folder that also holds a second FASTA file and a results table. Read the DELETED and KEPT lists before typing anything: the second FASTA file, the input and its session file must be under KEPT. Answer with just Enter first (nothing may be deleted), then `YES`.

### H. Batch CLI (nucleotides and each alphabet)
Put each line in `commands.sh`. Use `--input`, `--output`, and for codons `--alphabet codon --frame 0`.
- [ ] **T80** `score-exact`, `rank-cores`, `mutation-families`, `gapped-motif` on NORAD (no `--alphabet`): same as before this branch.
- [ ] **T81** The same with `--alphabet codon`, `--alphabet latin`, `--alphabet amino-acid`, `--alphabet codon --translate`, and an alphabet file path.
- [ ] **T82** `--on-unknown error` stops at the first unknown symbol with a message that names it.
- [ ] **T83** A missing file, bad JSON, a motif column that is not there, `--alpha 5`, `--frame 2` with latin: each is one line starting `error:`.
- [ ] **T84** Nucleotide mode on English text: a warning on stderr that most letters were dropped; on a file with no nucleotides, an error.

### I. Does it find real structure and not invent it? (the part that matters most)
- [ ] **T90** The files in `controls/shuffled/` (Antony, whale, finch, genes): how many cores are "statistically supported"? Expect few or none; the README gives what I got, for the originals and the shuffles.
- [ ] **T91** `controls/planted/planted_letters.txt` and `planted_codons.fasta`: the planted core is found, with the right count (5; and 3 exact plus 2 one-codon variants). `planted_split_across_records.fasta` must give 3, not 4.
- [ ] **T92** Mark Antony as letters: are the repeated phrases real phrases? Compare with what you know of the speech.
- [ ] **T93** A gene family or a gene with repeated domains, as codons and as amino acids: the amino-acid reading should find at least what the codon reading finds, and more when synonymous codons differ.
- [ ] **T94** Whales and finches: do the repeated units look like known phrase structure? Do any repeats span two bouts (put each bout in its own record and compare)?
- [ ] **T95** Rank-core output on a long real text: do many cores look "supported" because the text is highly repetitive (expected, since real text repeats), and is the ranking sensible?

### J. Speed and size
- [ ] **T100** `calibration/data/human_cds.fasta` as codons: *Load*, discovery, *Rank core motifs*. Expect seconds, not minutes. Note times.
- [ ] **T101** `controls/hard_cases/big_nt_300kb.fasta`, then `big_nt_3mb.fasta`: note load time, memory, and where it becomes impractical. Discovery on random nucleotides took 34 s at 300 kb and over 4 minutes at 3 Mb; real sequence repeats less than random text, so it should be faster.
- [ ] **T102** Ctrl+C during a long computation: returns cleanly or exits cleanly, no corrupted session file.

### K. Known gaps (do not report these as bugs)
Not built yet: several documents at once, with or without concatenation; finding similar documents in a
corpus; spelling variants such as *chuse* / *choose* (needs alignment-based matching); methylation as a fifth
letter (the projection `m` to `c` exists, class matching does not); the synonymous-codon enrichment test;
genetic codes other than the standard one in the menu; the relative-entropy logo for small alphabets; unit
distance matrices. In nucleotide mode, the batch CLI's default settings call cores "statistically supported" in
random text (documented in the README and CHANGELOG; unchanged on purpose).

## 6. Log entry (copy per observation into `notes.md`)

```
ID:            T52-03
Date / commit: 2026-.. / d502fb2
Input:         inputs/antony.txt  (sha256 ...)
Reading:       Latin letters           (menu option 5)
Path:          Plots > 4 (unit profile)
Typed:         motif "the", flank 5, mismatches 0, save PNG
Files:         runs/T52-antony/terminal.log, ..._profile_the.png
Expected:      ...
Saw:           ...
Severity:      wrong result / crash / confusing message / cosmetic / question
```

When you want me to look: tell me the run folders and the IDs. I will install the same commit, repeat
your steps from the logs and sessions, and tell you whether I get the same result before changing anything.
