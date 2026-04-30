# A simple menu for selecting RNAlang functions
import os
import re
import sys
import glob
import webbrowser
from os import getpid, chdir
from os.path import dirname
from re import sub
from math import ceil, floor, log10, log, exp

from rna_lexis.dialogs import openFile, openDir
from rna_lexis.algorithms import (
    count_kgrams, contains_only_rna, cover, find_boundary, cores,
    find_all_matches, print_core, find_with_mutations, extend_match_pair,
    find_longest_extensions, gen_hairpins,
    scramble_kmer_pvalues, markov_kmer_pvalues, decompose_motif,
)
from rna_lexis.alignment import gotoh_global, gotoh_local, print_alignment
from rna_lexis.plots import (
    plot_logo, plotzscore, plotkmerhist, plot_frequency_rank,
    plot_seq_nbrs, plot_coverage, plot_sequence_hits, plot_sequence_hits_detailed,
)
from rna_lexis.io import (
    read_text, save_session, load_session, init_summary, is_valid_session,
    open_file_with_default_software, _find_valid_sessions, fetch_enst_cdna,
    load_prefs, save_prefs,
)


class EOFSignal(Exception):
    """Exception raised when Ctrl+D (EOF) is detected to signal return to main menu."""
    pass


def _prompt_save(plotly=False):
    """Prompt for an output file name and format.

    Args:
        plotly: when True, include HTML as an extra format option.

    Returns:
        (filepath, scale) — filepath is '' if the user chose screen-only.
    """
    fn = safe_input(fmttxt(['Output file name',
                             '[leave blank to show on screen only]: '],
                            ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    scale = 1
    if fn:
        root = os.path.splitext(fn)[0]
        fmt_hint = ('[1: PNG standard (default), 2: PNG high-res 3x, 3: SVG, 4: HTML]'
                    if plotly else
                    '[1: PNG standard (default), 2: PNG high-res 3x, 3: SVG]')
        fmt = safe_input(fmttxt(['Output format', fmt_hint],
                                 ['bold', ''], ['yellow', 'cyan']) + ' ')
        if fmt == '2':
            fn, scale = f'{root}.png', 3
        elif fmt == '3':
            fn = f'{root}.svg'
        elif fmt == '4' and plotly:
            fn = f'{root}.html'
        else:
            fn = f'{root}.png'
        print(fmttxt([f"Saving as '{fn}'" + (' (high-res)' if scale > 1 else '')],
                     [''], ['cyan']))
    return fn, scale


def safe_input(prompt=''):
    """Wrapper for input() that catches Ctrl+D and raises EOFSignal.
    
    Args:
        prompt: The prompt string to display
        
    Returns:
        User input string
        
    Raises:
        EOFSignal: When user presses Ctrl+D
    """
    try:
        return input(prompt)
    except EOFError:
        print("\n(Returning to main menu...)")
        fn = globals().get('fn')
        if fn:
            print(f"Using file: {fn}")
        raise EOFSignal()


def fmttxt(txt, fnt, col, sep=' '):
    '''Format prompt text with ANSI codes.
    
    Args:
        txt: list of text parts
        fnt: list of formatting ('bold', 'dim', 'italic', 'underline') or empty strings
        col: list of colors ('black', 'red', 'green', 'yellow', 'blue', 'magenta', 'cyan', 'white') or empty strings
    
    Returns:
        Formatted string with text parts joined by spaces
    '''
    # Define mappings as dictionaries for O(1) lookup
    FORMATS = {
        'b': 1, 'bold': 1,
        'd': 2, 'dim': 2,
        'i': 3, 'italic': 3,
        'u': 4, 'underline': 4
    }
    
    COLORS = {
        'k': 30, 'black': 30,
        'r': 31, 'red': 31,
        'g': 32, 'green': 32,
        'y': 33, 'yellow': 33,
        'b': 34, 'blue': 34,
        'm': 35, 'magenta': 35,
        'c': 36, 'cyan': 36,
        'w': 37, 'white': 37
    }
    
    # Validate lengths
    if not (len(txt) == len(fnt) == len(col)):
        raise ValueError(f"Lists must have equal length: txt={len(txt)}, fnt={len(fnt)}, col={len(col)}")
    
    result = []
    for text, fmt, color in zip(txt, fnt, col):
        codes = []
        
        if fmt:
            key = fmt.lower()
            if key not in FORMATS:
                raise ValueError(f"Invalid format: {fmt}. Use: bold, dim, italic, underline")
            codes.append(str(FORMATS[key]))
        
        if color:
            key = color.lower()
            if key not in COLORS:
                raise ValueError(f"Invalid color: {color}")
            codes.append(str(COLORS[key]))
        
        if codes:
            text = f"\033[{';'.join(codes)}m{text}\033[0m"  # Add reset code
        
        result.append(text)
    
    return sep.join(result)


    
def print_hdr(fn, clr=True):
    """Print a session header line with the current file path.

    Args:
        fn:  File or session name to display after the 'File:' label.
        clr: When True (default), clear the terminal before printing.
             When False, only emit a blank line.
    """
    if clr == True:
        print("\033c", end="")
        # Package name banner: white text on light-blue background
        print("\033[104m\033[97m     RNA_lexis     \033[0m")
        print(fmttxt(["\nFile:\t", fn],['bold', ''], ['white', 'white']))
        print('\n')
    else:
        print('\n')

def choose_file(datadir=''):
    """Open a file-picker dialog and load the chosen file as a session dict.

    If the user selects a .json file it is loaded directly via load_session().
    Otherwise the file is read as plain text: non-alphabetic characters are
    stripped, RNA/DNA detection is applied (U→T), and a minimal session dict
    is returned.

    Args:
        datadir: Initial directory for the file-picker (default: system default).

    Returns:
        A session dict with keys ``file_path``, ``txt``, ``is_rna``, ``txtb``,
        and ``dir``.
    """
    file_path = openFile(initial_dir=datadir or None)
    if not file_path or not os.path.isfile(file_path):
        return None
    if file_path.endswith(".json"):
        return(load_session(file_path))
    else:
        # Detect FASTA format and strip header lines before building txt/txtb.
        gene_name = None
        with open(file_path, "r", encoding="utf-8", errors="ignore") as fh:
            raw_lines = fh.readlines()
        if raw_lines and raw_lines[0].startswith(">"):
            gene_name = raw_lines[0][1:].split()[0]
            raw_lines = [ln for ln in raw_lines if not ln.startswith(">")]
        raw = "".join(raw_lines)

        # Build txtb (spacing preserved) and txt (compact) from sequence only.
        txtb = raw.replace("\n", " ").lower()
        txtb = re.sub('[^a-z]+', ' ', txtb)
        txtb = re.sub(' +', ' ', txtb).strip()
        txt = sub(' +', '', txtb)

        if not txt:
            print(f"Error: '{file_path}' contains no alphabetic characters. "
                  "The file may use a non-Latin script or be empty.")
            return None
        is_rna = contains_only_rna(txt)
        if is_rna:
            txt = re.sub(r"[^A-Za-z]+", "", txt).lower().replace("u", "t")
    result = {'file_path': file_path, 'txt': txt,
              'is_rna': is_rna, 'txtb': txtb, 'dir': dirname(file_path)}
    if gene_name:
        result['gene_name'] = gene_name
    return result


def choose_enst(datadir=''):
    """Interactively prompt for an Ensembl transcript ID and load its cDNA.

    Before hitting the network, the function looks for a cached session JSON
    whose filename starts with the bare transcript ID in datadir.  If a cached
    file is found it is loaded and returned immediately.  Otherwise the cDNA
    sequence is fetched via the Ensembl REST API, the user is prompted for a
    save directory, and a minimal session dict is returned.

    Args:
        datadir: Directory to search for cached session files (default: '').

    Returns:
        A session dict with keys ``file_path``, ``txt``, ``txtb``, ``is_rna``,
        ``source`` (``'ensembl'``), ``enst``, and ``dir``.
    """
    while True:
        enst = safe_input(fmttxt(["Enter Ensembl transcript ID (ENST...): "],
                                 ['bold'], ['yellow']))
        enst = enst.strip()
        if not enst:
            print(fmttxt(["Please enter a transcript ID."], ['bold'], ['red']))
            continue
        base_id = enst.split('.')[0]  # strip version suffix (e.g. ENST00000123.5 -> ENST00000123)
        if datadir and os.path.isdir(datadir):
            matches = [os.path.join(root, f)
                       for root, _, files in os.walk(datadir)
                       for f in files
                       if f.startswith(base_id) and f.endswith('.json')]
            if matches:
                cached_path = matches[0]
                print(fmttxt([f"Loading cached session:", cached_path], ['bold', ''], ['green', 'cyan']))
                return load_session(cached_path)
        try:
            tid, label, seq = fetch_enst_cdna(enst)
            txt = sub('[^a-z]+', '', seq.lower())
            is_rna = contains_only_rna(txt)
            if is_rna:
                txt = re.sub(r"[^A-Za-z]+", "", txt).lower().replace("u", "t")
            pseudo_name = f"{tid}_{label}".replace(" ", "_")
            print(fmttxt([f"Fetched {tid} ({label}), length={len(txt)} nt"],
                         [''], ['green']))
            print(fmttxt([f"Choose the directory to save the data"],
                         [''], ['yellow']))
            path = openDir(initial_dir=datadir or None)
            print(path)
            return({
                'file_path': pseudo_name,
                'txt': txt,
                'txtb': txt,
                'is_rna': is_rna,
                'source': 'ensembl',
                'enst': tid,
                'dir': path
            })
        except Exception as exc:
            print(fmttxt([f"Failed to fetch ENST sequence: {exc}"], ['bold'], ['red']))


def load_from_paste(datadir=''):
    """Accept a pasted nucleotide (or text) sequence from the terminal.

    Reads lines until an empty line is entered, strips non-alphabetic
    characters, applies RNA/DNA detection (U→T), and prompts for a name
    and a save directory.

    Args:
        datadir: Initial directory for the save-directory picker (default: '').

    Returns:
        A session dict with keys ``file_path``, ``txt``, ``txtb``, ``is_rna``,
        ``source`` (``'paste'``), and ``dir``.

    Raises:
        ValueError: If no sequence characters are found in the pasted text.
    """
    print(fmttxt(["Paste sequence. Press Enter on an empty line to finish."], ['bold'], ['yellow']))
    first = input().strip()
    if not first:
        raise ValueError("No sequence pasted.")

    lines = [first]
    while True:
        line = input()
        if line.strip() == "":
            break
        lines.append(line.strip())

    txt = sub('[^a-zA-Z]+', '', "".join(lines)).lower()
    if not txt:
        raise ValueError("No sequence characters found in pasted text.")

    is_rna = contains_only_rna(txt)
    if is_rna:
        txt = re.sub(r"[^A-Za-z]+", "", txt).lower().replace("u", "t")
    print(fmttxt([f"Pasted sequence, length={len(txt)} nt"], [''], ['green']))

    default_name = f"PASTED_{len(txt)}nt"
    name = safe_input(fmttxt(["Enter a name for this sequence", f"[{default_name}]"],
                              ['bold', ''], ['yellow', 'cyan']) + " ")
    if name == '':
        name = default_name

    print(fmttxt(["Choose the directory to save the data"], [''], ['yellow']))
    path = openDir(initial_dir=datadir or None)

    return {
        'file_path': name,
        'txt': txt,
        'txtb': txt,
        'is_rna': is_rna,
        'source': 'paste',
        'dir': path
    }


def choose_input_source(datadir=''):
    """Show the input-source selector and delegate to the chosen loader.

    Presents a menu with four options: local file, Ensembl ENST fetch,
    paste, or Python prompt.  Returns -1 when the user selects the
    'Python prompt' option (the caller should handle this as a no-op or
    REPL mode).

    Args:
        datadir: Passed through to the chosen loader as the initial directory.

    Returns:
        A session dict from the selected loader, or -1 if the user chose
        'Python prompt'.
    """
    last_item = "Python prompt" if hasattr(sys, 'ps1') else "Quit"
    source_menu = ["Load from local file", "Fetch by Ensembl transcript ID (ENST)", "Paste sequence", last_item]
    val = show_menu("Input Source", "Choose Input", source_menu, clr=False)
    if val == len(source_menu):
        return -1
    if val == 1:
        return choose_file(datadir)
    if val == 2:
        return choose_enst(datadir)
    if val == 3:
        return load_from_paste(datadir)
    

# fmttxt([menu_ttl[menu_level]], ['bold'], ['yellow'])

def gen_menu(ttl, sttl, opts, clr='True', split=None):
    """Render a numbered menu to stdout.

    Clears the screen (via print_hdr) if clr is truthy, prints a decorated
    sub-title bar, and lists each option prefixed with its 1-based index.

    Args:
        ttl:   Top-level header / file name passed to print_hdr.
        sttl:  Sub-title shown centred inside a horizontal rule.
        opts:  List of option strings; each is printed as ``N. <option>``.
        clr:   Passed to print_hdr to control screen clearing (default True).
        split: If given, items 1..split are styled as primary (bold cyan) and
               items split+1 onward as utility (dim white).  Back/Quit are
               always bold red regardless of split.
    """
    print_hdr(ttl, clr)
    if not opts:
        return
    mxlen = max(len(s) for s in opts) + floor(log10(len(opts))) + 3
    if mxlen%2 == 1:
        mxlen += 1
    declen = max(ceil((mxlen-len(sttl)-2)/2), 0)
    hbarlen = 2*declen + 2 + len(sttl)
    sttl = fmttxt([sttl], ['bold'], ['green'])
    hd = '-'*declen + f" {sttl} " + '-'*declen
    print(hd)
    for i in range(len(opts)):
        if opts[i] == 'Quit':
            print(fmttxt([f"{i+1}. {opts[i]}"], ['bold'], ['red']))
        elif opts[i] == 'Back':
            print(fmttxt([f"{i+1}. {opts[i]}"], ['bold'], ['magenta']))
        elif split is not None and i >= split:
            print(fmttxt([f"{i+1}. {opts[i]}"], ['dim'], ['white']))
        elif split is not None:
            print(fmttxt([f"{i+1}. {opts[i]}"], ['bold'], ['cyan']))
        else:
            print(f"{i+1}. {opts[i]}")
    print('-'*hbarlen + '\n')


def show_menu(fn, ttl, submenu, clr: True, split=None):
    """Display a menu and block until the user enters a valid selection.

    Calls gen_menu() to render the numbered list, then loops reading input
    until an integer in [0, len(submenu)] is entered.  Entering 0 or the
    last option is typically treated as 'Back' / 'Quit' by callers.

    Args:
        fn:      File/session name forwarded to gen_menu as the header.
        ttl:     Sub-title of the menu.
        submenu: List of option strings.
        clr:     Passed to gen_menu to control screen clearing.
        split:   Forwarded to gen_menu; items before this index are styled as
                 primary, items from this index onward as utility.

    Returns:
        Integer in [0, len(submenu)] representing the user's choice.
        Returns -1 on EOFError (Ctrl+D).
    """
    gen_menu(fn, ttl, submenu, clr, split=split)
    while True:
        try:
            main_prompt = fmttxt(['Select a menu option', f'[1-{len(submenu)}]'], 
                         ['bold', ''], 
                         ['yellow', 'cyan'])
            slct = safe_input(main_prompt  + " ")
        except EOFError:
            print("\n")
            return -1
        try:
            vl = int(slct)
            if 0 <= vl <= len(submenu):
                return(int(slct))
            else:
                print(fmttxt([f"\nPlease enter a value between 1 and {len(submenu)}.\n"],
                      ['bold'], ['red']))
        except ValueError:
            print(fmttxt(["\nInvalid input! Please use numeric digits.\n"],['bold'], ['red']))


def neighbors_input(fn, txt, strs):
    """Interactive handler for the 'Core neighbors' plot.

    Prompts for a query sequence, neighbourhood half-width, the string list
    to use (cores or xmotifs), a plot title, an optional x-axis range, an
    output file, and whether to overlay hairpin regions.  Delegates to
    plot_seq_nbrs().

    Args:
        fn:   Current session file path (used for the default plot title and
              for locating a companion ``*_hairpins.csv`` file).
        txt:  Full source sequence.
        strs: Session strings dict with keys ``'corelist'`` and ``'xmotifs'``.
    """
    file_path = fn
    seq = safe_input(fmttxt(["Enter the sequence to analyze: "], ['bold'], ['yellow']))
    if len(find_all_matches(seq, txt)) < 1:
        print(fmttxt(["\nSequence not found."],['bold'], ['red']))
        safe_input(fmttxt(["Press any key to continue"], [''], ['white']))
        return None

    seq = seq.lower()
    wd_prompt = fmttxt(['Enter neighborhood width:', '[default: 40]'],
                         ['bold', ''],
                         ['yellow', 'cyan'])
    wd = safe_input(wd_prompt + " ")
    if wd == '':
        wd = 40
    else:
        wd = int(wd)
    wds_prompt = fmttxt(['Enter the strings to include in calculation ',
                        '[1: cores, or 2: xmotifs (default)]'],
                         ['bold', ''],
                         ['yellow', 'cyan'])
    wds = safe_input(wds_prompt + " ")
    if wds == '1':
        wds = strs["corelist"]
    else:
        wds = strs["xmotifs"]
    default_ttl = f"{os.path.basename(fn)} {seq}"
    ttl = safe_input(fmttxt(["Enter a title for the plot ", f"[Default: '{default_ttl}']"],
                       ['bold', ''], ['yellow', 'cyan']) + " ")
    if ttl == '':
        ttl = default_ttl
    fn, scale = _prompt_save(plotly=True)
    xrange = safe_input(fmttxt(["Enter the range to plot (leave blank for all, or use 'min, max')",
                    "[Default: '']"], ['bold', ''], ['yellow', 'cyan']) + " ")
    if xrange != '':
        try:
            parts = xrange.split(',')
            xrange = [int(parts[0].strip()), int(parts[1].strip())] if len(parts) == 2 else []
        except ValueError:
            xrange = []
    hairpins = []
    hairpins_csv = os.path.splitext(file_path)[0] + '_hairpins.csv'
    if os.path.isfile(hairpins_csv):
        ans = safe_input(fmttxt(['Show hairpin regions?', '[y/N]'],
                                ['bold', ''], ['yellow', 'cyan']) + ' ')
        if ans.strip().lower() == 'y':
            import csv
            with open(hairpins_csv, newline='') as _f:
                for row in csv.DictReader(_f):
                    hairpins.append({
                        'start':    int(row['start']),
                        'end':      int(row['end']),
                        'stem_seq': row.get('stem_seq', ''),
                        'loop_seq': row.get('loop_seq', ''),
                        'sequence': row.get('sequence', ''),
                    })
    plot_seq_nbrs(seq, wds, txt, sortby='CP', wd=wd, title=ttl, file=fn,
                  xrange=xrange, scale=scale, hairpins=hairpins)


def kmers_input(fn, txt, defvals):
    """Interactive handler for k-mer analysis plots.

    Prompts for a k-mer length and presents a sub-menu with four plot types:
    classical z-score, robust z-score, abundance histogram, and frequency-rank
    (Zipf) plot.  Each option further prompts for an output file path.

    Args:
        fn:      Current session file path (used as the menu header).
        txt:     Full source sequence.
        defvals: Global defaults dict (provides ``'clr'`` for menu rendering).
    """
    kmers_menu = ['Z-score', 'Robust Z-score', 'Abundance histogram',
                  'Rank-frequency', 'Return to plots']
    while True:
        k = safe_input(fmttxt(["Enter the length of the k-mers", "[Default: 5]"],
                     ['bold',''], ["yellow", "cyan"]) + " ")
        if k == '':
            k = 5
            break
        else:
            try:
                vl = int(k)
                if 1 <= vl <= len(txt):
                    k = vl
                    break
                else:
                    print(fmttxt([f"\nPlease enter an integer greater than 0.\n"],
                          ['bold'], ['red']))
            except ValueError:
                print(fmttxt(["\nInvalid input! Please enter an integer greater than 0.\n"],['bold'], ['red']))
    val = show_menu(fn, "K-mers plots", kmers_menu,
                            clr=defvals['clr'])
    if val == len(kmers_menu):
        return None
    kmers = count_kgrams(txt, k)
    if val == 1 or val == 2:
        zmin = safe_input(fmttxt(["Enter the z-score threshold", "[Default: 1.96]"],
                            ['bold',''], ["yellow", "cyan"]) + " ")
        if zmin == '':
            zmin = 1.96
        fn_out, scale = _prompt_save()
        if val == 1:
            plotzscore(kmers, float(zmin), robust=False, file=fn_out, scale=scale)
        else:
            plotzscore(kmers, float(zmin), robust=True, file=fn_out, scale=scale)
        return None
    if val == 3:
        fn_out, scale = _prompt_save()
        plotkmerhist(kmers, k, file=fn_out, scale=scale)
        return None
    if val == 4:
        fn_out, scale = _prompt_save()
        plot_frequency_rank(kmers, k, file=fn_out, scale=scale)


def logo_input(txt, strs):
    """Interactive handler for generating a WebLogo sequence logo.

    Prompts for a core/motif sequence, the number of flanking bases to
    include on each side, the maximum number of allowed mutations, and the
    output format (PDF or SVG).  Delegates to plot_logo().

    Args:
        txt:  Full source RNA/DNA sequence.
        strs: Session strings dict (currently unused, reserved for future
              use as a pre-filled sequence suggestion).
    """
    s0 = safe_input(fmttxt(["Enter the sequence to analyze: "], ['bold'], ['yellow']))
    s0 = sub(r'[^a-zA-Z]', '', s0).lower()
    if len(find_all_matches(s0, txt)) < 1:
        print(fmttxt(["\nSequence not found."],['bold'], ['red']))
        safe_input(fmttxt(["Press any key to continue"], [''], ['white']))
        return(None)
    while True:
        k = safe_input(fmttxt(["Enter the number of bases on both sides", "[Default: 0]"],
                     ['bold',''], ["yellow", "cyan"]) + " ")
        if k == '':
            k = 0
            break
        else:
            try:
                vl = int(k)
                if 0 <= vl <= min(len(txt), 100):
                    k = vl
                    break
                else:
                    print(fmttxt([f"\nPlease enter an integer greater than or equal to 0.\n"],
                          ['bold'], ['red']))
            except ValueError:
                print(fmttxt(["\nInvalid input! Please enter an integer greater than or equal to 0.\n"], ['bold'], ['red']))
    while True:
        nmut = safe_input(fmttxt(["Enter the maximum number of mutations in the sequence", "[Default: 0]"],
                     ['bold',''], ["yellow", "cyan"]) + " ")
        if nmut == '':
            nmut = 0
            break
        else:
            try:
                vl = int(nmut)
                if 0 <= vl <= round(len(s0)/6):
                    nmut = vl
                    break
                else:
                    print(fmttxt([f"\nPlease enter an integer greater than or equal to 0 and less than or equal to {round(len(s0)/6)}.\n"],
                          ['bold'], ['red']))
            except ValueError:
                print(fmttxt(["\nInvalid input! Please enter an integer greater than or equal to 0.\n"], ['bold'], ['red']))
    fmt = ""
    while fmt not in ["pdf", "svg"]:
        fmt = safe_input(fmttxt(["Enter the format", "[pdf (default), or svg]"],
                    ['bold',''], ["yellow", "cyan"]) + " ")
        if fmt == '':
            fmt = "pdf"
        fmt = fmt.lower()
    fn = s0 + "_logo"
    fn = safe_input(fmttxt([f"Enter the output file name: [{fn}]"],['bold'], ['yellow']))
    if fn == '':
        fn = s0 + "_logo"
    plot_logo(s0, k , txt, muts=nmut, outfile=fn, fmt=fmt)


def coverage_input(txt, strs):
    """Interactive handler for the weighted-coverage bar chart.

    Prompts for which string list to use (cores or xmotifs), the length
    exponent for the coverage score, and an output file path.  Computes
    cover() for each selected string and delegates to plot_coverage().

    Args:
        txt:  Full source sequence.
        strs: Session strings dict with keys ``'corelist'`` and ``'xmotifs'``.
    """
    wds = safe_input(fmttxt(["Enter the strings to include in calculation",
                        "[1: cores (default), 2: xmotifs]"], ['bold',''], ["yellow", "cyan"]) + " ")
    try:
        wds = strs["corelist"] if (wds == '' or int(wds) == 1) else strs["xmotifs"]
    except ValueError:
        wds = strs["corelist"]
    pwr = safe_input(fmttxt(["Enter the exponent a in length^a * occurrences", "[Default: 1.2]"],
                        ['bold',''], ["yellow", "cyan"]) + " ")
    if pwr == '':
        pwr = 1.2
    else:
        pwr = float(pwr)
    cvrs = dict()
    for i in range(len(wds)):
        cvrs[wds[i]] = cover(wds[i], txt, pwr=pwr)
    fn_out, scale = _prompt_save()
    plot_coverage(txt, cvrs, pwr, file=fn_out, scale=scale)


def find_match_input(txt, strs, minlen=4):
    """Interactive handler for regex-based sequence search.

    Prompts for a query string (supports ``.`` as a single-character wildcard),
    calls find_all_matches(), and prints each match position and matched string.

    Args:
        txt:    Full source sequence.
        strs:   Session strings dict (currently unused, present for interface
                consistency with other handlers).
        minlen: Minimum query length accepted (default 4).
    """
    while True:
        s0 = safe_input(fmttxt(["Enter the sequence to analyze (use . as a single letter wild card): "],
                      ['bold'], ["yellow"]) + " ")
        s0 = s0.lower()
        if len(s0) >= minlen:
            break
        else:
            print(fmttxt(["Sequence too short"], ['bold'], ['red']))
            safe_input("Press any key to continue.")
    ret = find_all_matches(s0, txt)
    retstr = find_all_matches(s0, txt, ret="str")
    print(f"\nFound {len(ret)} matches for {s0}:\n")
    for i in range(len(ret)):
        print(f"Match {i+1} at position {ret[i]}, match: {retstr[i]}")
    print('\n')
    safe_input("Press Enter to continue...")


def txt_coverage_input(txt, strs, minlen=4):
    """Interactive handler for displaying the weighted coverage of a single sequence.

    Prompts for a query string and a length exponent, then prints the
    coverage score (cover()) for that string in txt.

    Args:
        txt:    Full source sequence.
        strs:   Session strings dict (currently unused).
        minlen: Minimum query length accepted (default 4).
    """
    while True:
        s0 = safe_input(fmttxt(["Enter the sequence to analyze (use . as a single letter wild card): "],
                      ['bold'], ["yellow"]) + " ")
        s0 = s0.lower()
        if len(s0) >= minlen:
            break
        else:
            print(fmttxt(["Sequence too short"], ['bold'], ['red']))
            safe_input("Press any key to continue.")
    pwr = safe_input(fmttxt(["Enter the exponent a in length^a * occurrences", "[Default: 1]"],
                       ['bold',''], ["yellow", "cyan"]))
    if pwr == '':
        pwr = 1
    else:
        pwr = float(pwr)
    print(f"{s0}, length {len(s0)}, coverage {cover(s0, txt, pwr=pwr)}, out of {len(txt)}")
    safe_input("Press Enter to continue...")
 

def print_core_input(txt, strs, minlen=4):
    """Interactive handler for aligned xmotif display around a core.

    Prompts for a core sequence and prints all xmotifs that contain it,
    left-padded so the core aligns vertically, with occurrence counts
    prepended.  Delegates to print_core().

    Args:
        txt:    Full source sequence.
        strs:   Session strings dict; uses ``'xmotifs'`` as the candidate list.
        minlen: Minimum query length accepted (default 4).
    """
    while True:
        s0 = safe_input(fmttxt(["Enter the sequence to analyze: "],
                      ['bold'], ["yellow"]))
        s0 = s0.lower()
        if len(s0) >= minlen:
            break
        else:
            print(fmttxt(["Sequence too short"], ['bold'], ['red']))
            safe_input("Press any key to continue.")
    commoncore = print_core(txt, s0, strs["xmotifs"])
    if len(commoncore) > 1:
        print('\n'.join(commoncore) + '\n')
    safe_input("Press Enter to continue...")


def print_alignment_score(txt):
    """Interactive handler for pairwise sequence alignment and scoring.

    Prompts for two start positions and a common length within txt, then
    an alignment mode (local or global).  Runs the chosen Gotoh alignment,
    prints the formatted alignment, and reports the normalised score, bit
    score, and E-value (using standard nucleotide BLAST constants λ=1.28,
    K=0.46).

    Args:
        txt: Full source sequence.
    """
    n = len(txt)
    while True:
        p1 = safe_input(fmttxt(["Enter start position of first sequence", f"[0-{n-1}]"],
                                ['bold', ''], ['yellow', 'cyan']) + " ")
        try:
            p1 = int(p1)
            if 0 <= p1 < n:
                break
            print(fmttxt([f"Please enter a value between 0 and {n-1}."], ['bold'], ['red']))
        except ValueError:
            print(fmttxt(["Invalid input! Please enter an integer."], ['bold'], ['red']))

    while True:
        p2 = safe_input(fmttxt(["Enter start position of second sequence", f"[0-{n-1}]"],
                                ['bold', ''], ['yellow', 'cyan']) + " ")
        try:
            p2 = int(p2)
            if 0 <= p2 < n:
                break
            print(fmttxt([f"Please enter a value between 0 and {n-1}."], ['bold'], ['red']))
        except ValueError:
            print(fmttxt(["Invalid input! Please enter an integer."], ['bold'], ['red']))

    max_len = min(n - p1, n - p2)
    while True:
        ln = safe_input(fmttxt(["Enter sequence length", f"[1-{max_len}]"],
                                ['bold', ''], ['yellow', 'cyan']) + " ")
        try:
            ln = int(ln)
            if 1 <= ln <= max_len:
                break
            print(fmttxt([f"Please enter a value between 1 and {max_len}."], ['bold'], ['red']))
        except ValueError:
            print(fmttxt(["Invalid input! Please enter an integer."], ['bold'], ['red']))

    s1 = txt[p1:p1 + ln]
    s2 = txt[p2:p2 + ln]

    while True:
        mode = safe_input(fmttxt(["Alignment type", "[1: Local (default), 2: Global]"],
                                  ['bold', ''], ['yellow', 'cyan']) + " ")
        if mode in ('', '1'):
            res = gotoh_local(s1, s2)
            break
        elif mode == '2':
            res = gotoh_global(s1, s2)
            break
        print(fmttxt(["Please enter 1 or 2."], ['bold'], ['red']))

    print_alignment(res)

    # ── Normalised score (self-alignment) ────────────────────────────────
    # Self-alignment score = ln * match (match default = 2)
    self_score = ln * 2
    norm = res.score / self_score if self_score > 0 else 0.0
    print(fmttxt([f"Normalised score (self-alignment):  {norm:.3f}  "
                  f"[self-score = {self_score}]"], [''], ['cyan']))

    # ── Bit score and E-value ────────────────────────────────────────────
    # Constants for nucleotide scoring (λ=1.28, K=0.46; Altschul et al. 1990)
    # E-value uses the full transcript as the reference "database".
    _LAM, _K = 1.28, 0.46
    if res.score > 0:
        bit = (_LAM * res.score - log(_K)) / log(2)
        eval = _K * ln * len(txt) * exp(-_LAM * res.score)
        print(fmttxt([f"Bit score:                          {bit:.1f} bits"], [''], ['cyan']))
        print(fmttxt([f"E-value (vs transcript len {len(txt)}):  {eval:.2e}"], [''], ['cyan']))
    else:
        print(fmttxt(["Bit score / E-value: not applicable (score ≤ 0)"], [''], ['cyan']))

    safe_input(fmttxt(["Press Enter to return to the menu"], [''], ['white']))


def search_input(txt, mutr=1/6, M=4):
    """Prompt for a query string and find exact + approximate matches in txt.

    Mismatched characters are highlighted in red in the terminal output.
    """
    RED   = '\033[31m'
    RESET = '\033[0m'

    prompt = fmttxt(['Enter query string', '[blank to cancel]: '],
                    ['bold', ''], ['yellow', 'cyan'])
    seq = safe_input(prompt + ' ').strip().lower()
    if not seq:
        return

    mutr_str = safe_input(fmttxt(['Mutation rate: 1 per N letters', f'[default: {1/mutr:.1f}]: '],
                                   ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    if mutr_str:
        val = float(mutr_str)
        mutr = 1 / val if val > 0 else mutr

    exact_pos, approx, maxmut = find_with_mutations(seq, txt, mutr=mutr, M=M)
    L = len(seq)

    # ── header ──────────────────────────────────────────────────────────
    mut_note = (f'max {maxmut} mutation{"s" if maxmut != 1 else ""}'
                if maxmut else 'too short for mutation search')
    print(fmttxt([f'\n"{seq}"', f'L={L}, {mut_note}'], ['bold', ''], ['white', 'cyan']))

    # ── exact matches ────────────────────────────────────────────────────
    if exact_pos:
        print(fmttxt([f'  Exact ({len(exact_pos)}):'], ['bold'], ['green']))
        print(f'    positions: {exact_pos}')
    else:
        print(fmttxt(['  No exact matches'], [''], ['red']))

    # ── approximate matches ───────────────────────────────────────────────
    if maxmut > 0:
        if approx:
            print(fmttxt([f'  Approximate ≤{maxmut} mut ({len(approx)}):'], ['bold'], ['yellow']))
            for pos, matched, d in approx:
                highlighted = ''.join(
                    f'{RED}{mc.upper()}{RESET}' if mc != sc else mc
                    for sc, mc in zip(seq, matched)
                )
                print(f'    pos {pos:5d}  "{highlighted}"  dist={d}')
        else:
            print(fmttxt([f'  No approximate matches (≤{maxmut} mutations)'], [''], ['cyan']))

    print()
    safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))


def extend_match_input(txt, default_mutr=1/6):
    """Prompt for a seed sequence and mutation rate, then find and display the
    longest Hamming-bounded extension for every pair of exact seed occurrences."""
    RED   = '\033[31m'
    BOLD  = '\033[1m'
    RESET = '\033[0m'
    DIM   = '\033[2m'
    CYAN  = '\033[36m'

    # ── seed ────────────────────────────────────────────────────────────────
    seq = safe_input(fmttxt(['Enter seed sequence', '[blank to cancel]: '],
                             ['bold', ''], ['yellow', 'cyan']) + ' ').strip().lower()
    if not seq:
        return
    positions = find_all_matches(seq, txt, ret='pos')
    if len(positions) < 2:
        print(fmttxt([f'Seed "{seq}" found {len(positions)} time(s) — need ≥ 2.'],
                     ['bold'], ['red']))
        safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))
        return

    # ── mutation rate ────────────────────────────────────────────────────────
    mutr_str = safe_input(fmttxt(['Mutation rate: 1 per N letters',
                                   f'[default: {round(1/default_mutr)}]: '],
                                  ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    val = float(mutr_str) if mutr_str else round(1 / default_mutr)
    mutr = 1 / val if val > 0 else default_mutr

    # ── compute extensions ───────────────────────────────────────────────────
    print(fmttxt([f'Computing extensions for {len(positions)} occurrences...'],
                 [''], ['cyan']))
    results = find_longest_extensions(seq, txt, mutr=mutr)

    L = len(seq)
    max_N = 10       # show at most this many pairs
    print(fmttxt([f'\nFound {len(results)} pair(s). '
                  f'Showing top {min(max_N, len(results))}, '
                  f'sorted by total extension length.\n'],
                 ['bold'], ['white']))

    for idx, res in enumerate(results[:max_N]):
        p1      = res['pos1']
        p2      = res['pos2']
        left_e  = res['left_ext']
        right_e = res['right_ext']
        tlen    = res['total_len']
        h       = res['hamming']
        e1      = res['ext1']
        e2      = res['ext2']
        rate    = f"{h}/{tlen} ({100*h/tlen:.1f}%)" if tlen else '0/0'

        print(f"{BOLD}  Pair {idx+1}{RESET}  "
              f"pos1={CYAN}{p1}{RESET}  pos2={CYAN}{p2}{RESET}  "
              f"total_len={CYAN}{tlen}{RESET}  "
              f"(left={left_e}, seed={L}, right={right_e})  "
              f"hamming={CYAN}{rate}{RESET}")

        # Build display strings: lower-case extension, upper-case seed,
        # red upper-case for mismatches in extension parts.
        def _fmt_ext(ref, qry, left, right, sl):
            """Colour-code the extension around the seed (positions 0..left-1
            are the left flank, left..left+sl-1 are the seed, rest is right)."""
            out = []
            for i, (rc, qc) in enumerate(zip(ref, qry)):
                in_seed = (left <= i < left + sl)
                mismatch = (rc != qc)
                ch = qc.upper() if (in_seed or mismatch) else qc
                if mismatch:
                    ch = f"{RED}{ch}{RESET}"
                elif in_seed:
                    ch = f"{BOLD}{ch}{RESET}"
                out.append(ch)
            return ''.join(out)

        # Truncate very long sequences for display
        MAX_SHOW = 80
        if tlen <= MAX_SHOW:
            d1 = _fmt_ext(e1, e1, left_e, right_e, L)
            d2 = _fmt_ext(e1, e2, left_e, right_e, L)
            print(f"    [{DIM}{p1}{RESET}] {d1}")
            print(f"    [{DIM}{p2}{RESET}] {d2}")
        else:
            half = MAX_SHOW // 2
            head1 = _fmt_ext(e1[:half], e1[:half], left_e, right_e, L)
            head2 = _fmt_ext(e1[:half], e2[:half], left_e, right_e, L)
            print(f"    [{DIM}{p1}{RESET}] {head1}… ({tlen} nt total)")
            print(f"    [{DIM}{p2}{RESET}] {head2}…")
        print()

    # ── print one pair in full ────────────────────────────────────────────────
    shown = min(max_N, len(results))
    if shown > 0:
        pick_str = safe_input(
            fmttxt(['Print a pair in full', f'[1–{shown}, blank to skip]: '],
                   ['bold', ''], ['yellow', 'cyan']) + ' '
        ).strip()
        if pick_str.isdigit():
            pick = int(pick_str) - 1
            if 0 <= pick < shown:
                res    = results[pick]
                p1     = res['pos1']
                p2     = res['pos2']
                left_e = res['left_ext']
                e1     = res['ext1']
                e2     = res['ext2']
                tlen   = res['total_len']
                h      = res['hamming']
                print(f"\n{BOLD}  Pair {pick+1} — full sequences{RESET}  "
                      f"(total_len={CYAN}{tlen}{RESET}  hamming={CYAN}{h}{RESET})\n")
                print(f"  {DIM}pos1={p1}{RESET}")
                print(f"  {e1}\n")
                print(f"  {DIM}pos2={p2}{RESET}")
                print(f"  {e2}\n")
                d1 = _fmt_ext(e1, e1, left_e, right_e, L)
                d2 = _fmt_ext(e1, e2, left_e, right_e, L)
                print(f"  {DIM}(colour-coded, pos1 as reference){RESET}")
                print(f"  {d1}")
                print(f"  {d2}\n")

    # ── save all results to CSV ───────────────────────────────────────────────
    if results:
        import csv, re as _re
        slug = _re.sub(r'[^a-zA-Z0-9_-]', '_', seq)[:40]
        default_fn = f'extensions_{slug}.csv'
        fn_str = safe_input(
            fmttxt(['Save all pairs to CSV', f'[Enter for {default_fn}, Ctrl+D to skip]: '],
                   ['bold', ''], ['yellow', 'cyan']) + ' '
        )
        if fn_str is not None:
            fn_str = fn_str.strip()
            fn_csv = fn_str if fn_str else default_fn
            with open(fn_csv, 'w', newline='') as _f:
                w = csv.writer(_f)
                w.writerow(['rank', 'pos1', 'pos2', 'left_ext', 'right_ext',
                            'total_len', 'hamming', 'ext1', 'ext2'])
                for rank, res in enumerate(results, 1):
                    w.writerow([rank, res['pos1'], res['pos2'],
                                res['left_ext'], res['right_ext'],
                                res['total_len'], res['hamming'],
                                res['ext1'], res['ext2']])
            print(fmttxt([f'Saved {len(results)} pair(s) to:', fn_csv],
                         ['bold', ''], ['green', 'cyan']))
            open_file_with_default_software(fn_csv)

    safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))


def motif_extensions_input(txt):
    """Prompt for a motif, find matches (exact or with mutations), then extend
    each pair of exact occurrences and display lengths, Hamming distances, and
    colour-coded alignments."""
    RED   = '\033[31m'
    BOLD  = '\033[1m'
    RESET = '\033[0m'
    DIM   = '\033[2m'
    CYAN  = '\033[36m'
    GREEN = '\033[32m'

    # ── motif ───────────────────────────────────────────────────────────────
    motif = safe_input(fmttxt(['Enter motif', '[blank to cancel]: '],
                               ['bold', ''], ['yellow', 'cyan']) + ' ').strip().lower()
    if not motif:
        return

    # ── search method ────────────────────────────────────────────────────────
    method = safe_input(fmttxt(['Search method:', '[1] Exact  [2] With mutations  [default: 1]: '],
                                ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    use_mutations = (method == '2')

    if use_mutations:
        mutr_s_str = safe_input(fmttxt(['Search mutation rate: 1 per N letters', '[default: 6]: '],
                                        ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
        val_s = float(mutr_s_str) if mutr_s_str else 6
        mutr_search = 1 / val_s if val_s > 0 else 1 / 6
        exact_pos, approx, maxmut = find_with_mutations(motif, txt, mutr=mutr_search)
        n_exact = len(exact_pos)
        n_approx = len(approx)
        print(fmttxt([f'\nFound {n_exact} exact and {n_approx} approximate match(es) '
                      f'(max mismatches allowed: {maxmut}).'],
                     ['bold'], ['white']))
        MAX_LIST = 20
        for pos in exact_pos[:MAX_LIST]:
            print(f"  {GREEN}exact{RESET}   pos={CYAN}{pos}{RESET}  {txt[pos:pos+len(motif)]}")
        for pos, matched, dist in approx[:MAX_LIST]:
            print(f"  {RED}~dist={dist}{RESET}  pos={CYAN}{pos}{RESET}  {matched}")
        total_shown = min(MAX_LIST, n_exact + n_approx)
        if n_exact + n_approx > total_shown:
            print(f"  {DIM}... and {n_exact + n_approx - total_shown} more{RESET}")
    else:
        exact_pos = find_all_matches(motif, txt, ret='pos')
        n_exact = len(exact_pos)
        print(fmttxt([f'\nFound {n_exact} exact match(es).'], ['bold'], ['white']))
        MAX_LIST = 20
        for pos in exact_pos[:MAX_LIST]:
            print(f"  pos={CYAN}{pos}{RESET}  {txt[pos:pos+len(motif)]}")
        if n_exact > MAX_LIST:
            print(f"  {DIM}... and {n_exact - MAX_LIST} more{RESET}")

    print()

    # Build unified position list: (pos, seed_dist).
    # For approximate search include all found positions; for exact search, exact only.
    if use_mutations:
        all_positions = [(p, 0) for p in exact_pos] + [(p, d) for p, _, d in approx]
    else:
        all_positions = [(p, 0) for p in exact_pos]

    if len(all_positions) < 2:
        label = 'match(es)' if use_mutations else 'exact occurrence(s)'
        print(fmttxt([f'Need ≥ 2 {label} to compute extensions (found {len(all_positions)}).'],
                     ['bold'], ['red']))
        safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))
        return

    # ── extension mutation rate ───────────────────────────────────────────────
    mutr_e_str = safe_input(fmttxt(['Extension mutation rate: 1 per N letters', '[default: 6]: '],
                                    ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    val_e = float(mutr_e_str) if mutr_e_str else 6
    mutr_ext = 1 / val_e if val_e > 0 else 1 / 6

    # ── compute extensions ────────────────────────────────────────────────────
    # Extend every pair of positions.  For approximate matches, seed mismatches
    # (dist1 + dist2) are added to the Hamming total so the count is consistent
    # with what sequence_hits_input reports for the same motif.
    print(fmttxt([f'Computing extensions for {len(all_positions)} occurrence(s)...'], [''], ['cyan']))
    L = len(motif)
    results = []
    for i in range(len(all_positions)):
        for j in range(i + 1, len(all_positions)):
            p1, d1 = all_positions[i]
            p2, d2 = all_positions[j]
            l, r, e1, e2, flank_ham = extend_match_pair(motif, txt, p1, p2, mutr=mutr_ext)
            results.append({
                'pos1': p1, 'pos2': p2,
                'left_ext': l, 'right_ext': r,
                'total_len': L + l + r,
                'ext1': e1, 'ext2': e2,
                'hamming': d1 + d2 + flank_ham,
            })
    results.sort(key=lambda x: x['total_len'], reverse=True)

    max_N = 10
    print(fmttxt([f'\nFound {len(results)} pair(s). '
                  f'Showing top {min(max_N, len(results))}, '
                  f'sorted by total extension length.\n'],
                 ['bold'], ['white']))

    def _fmt_ext(ref, qry, left, sl):
        """Return a colour-coded string with the seed BOLD and mismatches RED.

        Args:
            ref: Reference sequence string.
            qry: Query sequence string (same length as ref).
            left: Number of left-extension characters before the seed starts.
            sl:  Length of the seed portion.
        """
        out = []
        for i, (rc, qc) in enumerate(zip(ref, qry)):
            in_seed  = (left <= i < left + sl)
            mismatch = (rc != qc)
            ch = qc.upper() if (in_seed or mismatch) else qc
            if mismatch:
                ch = f"{RED}{ch}{RESET}"
            elif in_seed:
                ch = f"{BOLD}{ch}{RESET}"
            out.append(ch)
        return ''.join(out)

    for idx, res in enumerate(results[:max_N]):
        p1      = res['pos1']
        p2      = res['pos2']
        left_e  = res['left_ext']
        right_e = res['right_ext']
        tlen    = res['total_len']
        h       = res['hamming']
        e1      = res['ext1']
        e2      = res['ext2']
        rate    = f"{h}/{tlen} ({100*h/tlen:.1f}%)" if tlen else '0/0'

        print(f"{BOLD}  Pair {idx+1}{RESET}  "
              f"pos1={CYAN}{p1}{RESET}  pos2={CYAN}{p2}{RESET}  "
              f"total_len={CYAN}{tlen}{RESET}  "
              f"(left={left_e}, seed={L}, right={right_e})  "
              f"hamming={CYAN}{rate}{RESET}")

        MAX_SHOW = 80
        if tlen <= MAX_SHOW:
            d1 = _fmt_ext(e1, e1, left_e, L)
            d2 = _fmt_ext(e1, e2, left_e, L)
            print(f"    [{DIM}{p1}{RESET}] {d1}")
            print(f"    [{DIM}{p2}{RESET}] {d2}")
        else:
            half = MAX_SHOW // 2
            head1 = _fmt_ext(e1[:half], e1[:half], left_e, L)
            head2 = _fmt_ext(e1[:half], e2[:half], left_e, L)
            print(f"    [{DIM}{p1}{RESET}] {head1}… ({tlen} nt total)")
            print(f"    [{DIM}{p2}{RESET}] {head2}…")
        print()

    # ── print one pair in full ────────────────────────────────────────────────
    shown = min(max_N, len(results))
    if shown > 0:
        pick_str = safe_input(
            fmttxt(['Print a pair in full', f'[1–{shown}, blank to skip]: '],
                   ['bold', ''], ['yellow', 'cyan']) + ' '
        ).strip()
        if pick_str.isdigit():
            pick = int(pick_str) - 1
            if 0 <= pick < shown:
                res    = results[pick]
                p1     = res['pos1']
                p2     = res['pos2']
                left_e = res['left_ext']
                e1     = res['ext1']
                e2     = res['ext2']
                tlen   = res['total_len']
                h      = res['hamming']
                print(f"\n{BOLD}  Pair {pick+1} — full sequences{RESET}  "
                      f"(total_len={CYAN}{tlen}{RESET}  hamming={CYAN}{h}{RESET})\n")
                d1 = _fmt_ext(e1, e1, left_e, L)
                d2 = _fmt_ext(e1, e2, left_e, L)
                print(f"  {DIM}pos1={p1}{RESET}")
                print(f"  {e1}\n")
                print(f"  {DIM}pos2={p2}{RESET}")
                print(f"  {e2}\n")
                print(f"  {DIM}(colour-coded, pos1 as reference){RESET}")
                print(f"  {d1}")
                print(f"  {d2}\n")

    # ── save all results to CSV ───────────────────────────────────────────────
    if results:
        import csv, re as _re
        slug = _re.sub(r'[^a-zA-Z0-9_-]', '_', motif)[:40]
        default_fn = f'extensions_{slug}.csv'
        fn_str = safe_input(
            fmttxt(['Save all pairs to CSV', f'[Enter for {default_fn}, Ctrl+D to skip]: '],
                   ['bold', ''], ['yellow', 'cyan']) + ' '
        )
        if fn_str is not None:
            fn_str = fn_str.strip()
            fn_csv = fn_str if fn_str else default_fn
            with open(fn_csv, 'w', newline='') as _f:
                w = csv.writer(_f)
                w.writerow(['rank', 'pos1', 'pos2', 'left_ext', 'right_ext',
                            'total_len', 'hamming', 'ext1', 'ext2'])
                for rank, res in enumerate(results, 1):
                    w.writerow([rank, res['pos1'], res['pos2'],
                                res['left_ext'], res['right_ext'],
                                res['total_len'], res['hamming'],
                                res['ext1'], res['ext2']])
            print(fmttxt([f'Saved {len(results)} pair(s) to:', fn_csv],
                         ['bold', ''], ['green', 'cyan']))
            open_file_with_default_software(fn_csv)

    safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))


def hairpins_input(txt, file_path, minSL=8, minLL=3, maxLL=40):
    """Prompt for hairpin parameters, run gen_hairpins, and export results to CSV."""
    import csv

    print(fmttxt(["\n--- Hairpin parameters ---"], ['bold'], ['cyan']))

    # --- minSL ---
    s = safe_input(fmttxt(['Min stem length (minSL)', f'[default: {minSL}]: '],
                           ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    if s:
        try:
            minSL = int(s)
        except ValueError:
            print(fmttxt(['Invalid value, using default.'], [''], ['red']))

    # --- minLL ---
    s = safe_input(fmttxt(['Min loop length (minLL)', f'[default: {minLL}]: '],
                           ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    if s:
        try:
            minLL = int(s)
        except ValueError:
            print(fmttxt(['Invalid value, using default.'], [''], ['red']))

    # --- maxLL ---
    s = safe_input(fmttxt(['Max loop length (maxLL)', f'[default: {maxLL}]: '],
                           ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    if s:
        try:
            maxLL = int(s)
        except ValueError:
            print(fmttxt(['Invalid value, using default.'], [''], ['red']))

    # --- output file ---
    default_out = os.path.splitext(file_path)[0] + '_hairpins.csv'
    fn = safe_input(fmttxt(['Output CSV file', f'[default: {os.path.basename(default_out)}]: '],
                            ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    if not fn:
        fn = default_out
    elif not os.path.isabs(fn):
        fn = os.path.join(os.path.dirname(file_path) or '.', fn)

    print(fmttxt(['Searching for hairpins...'], [''], ['cyan']))
    hairpins = gen_hairpins(txt, minSL=minSL, minLL=minLL, maxLL=maxLL)

    if not hairpins:
        print(fmttxt(['No hairpins found with the given parameters.'], ['bold'], ['red']))
        safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))
        return

    try:
        with open(fn, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['start', 'end', 'stem_len', 'loop_len', 'stem_seq', 'loop_seq', 'sequence'])
            for hp in hairpins:
                writer.writerow([hp.start, hp.end, hp.stem_len, hp.loop_len,
                                 hp.stem_seq, hp.loop_seq, txt[hp.start:hp.end]])
        print(fmttxt([f'Exported {len(hairpins)} hairpin(s) to:', fn],
                     ['bold', ''], ['green', 'cyan']))
        open_file_with_default_software(fn)
    except Exception as e:
        print(fmttxt([f'Error writing file: {e}'], ['bold'], ['red']))

    safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))


def markov_input(txt, file_path=''):
    """Prompt for k-mer length, then write a CSV of all observed k-mers with
    analytical Markov-model p-values — no shuffling required.

    For each k-mer the expected count under a (k-1)-th order Markov model is
    computed via the Prum/Schbath formula; a Poisson exact p-value tests
    whether the observed count is surprising given the (k-1)-mer frequencies.
    """
    k_str = safe_input(fmttxt(['k-mer length', '[default: 6]: '],
                               ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    k = int(k_str) if k_str.isdigit() and int(k_str) > 1 else 6

    base = os.path.splitext(file_path)[0] if file_path else 'sequence'
    default_csv = f'{base}_kmer{k}_markov.csv'
    fn_str = safe_input(fmttxt(['Output CSV file',
                                 f'[default: {os.path.basename(default_csv)}]: '],
                                ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    out_csv = fn_str if fn_str else default_csv

    print(fmttxt([f'Computing Markov p-values for k={k} …'], [''], ['cyan']))
    results = markov_kmer_pvalues(txt, k=k)

    import csv as _csv
    try:
        with open(out_csv, 'w', newline='') as fh:
            w = _csv.writer(fh)
            w.writerow([
                'kmer', 'real_count', 'expected_count',
                'pvalue_over', 'evalue_over', 'pvalue_over_bh',
                'pvalue_under', 'evalue_under', 'pvalue_under_bh',
                'direction',
            ])
            for row in results:
                w.writerow([
                    row['kmer'], row['real_count'],
                    f"{row['expected_count']:.3f}",
                    f"{row['pvalue_over']:.8f}",  f"{row['evalue_over']:.4f}",  f"{row['pvalue_over_bh']:.8f}",
                    f"{row['pvalue_under']:.8f}", f"{row['evalue_under']:.4f}", f"{row['pvalue_under_bh']:.8f}",
                    row['direction'],
                ])
        print(fmttxt([f'Saved {len(results)} k-mers to:', out_csv],
                     ['bold', ''], ['green', 'cyan']))
        open_file_with_default_software(out_csv)
    except Exception as e:
        print(fmttxt([f'Error writing file: {e}'], ['bold'], ['red']))

    safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))


def decompose_motif_input(txt, file_path=''):
    """Prompt for a motif, then write a CSV of hierarchical Markov decomposition.

    Every contiguous sub-k-mer of the motif is tested at each length k from
    len(motif) down to 2, revealing the shortest unit whose count is
    surprising beyond what the (k-1)-mer context already explains.
    """
    prompt = fmttxt(['Motif to decompose', '(e.g. ugcaug): '],
                    ['bold', ''], ['yellow', 'cyan'])
    motif = safe_input(prompt + ' ').strip().lower()
    if not motif:
        print(fmttxt(['No motif entered.'], ['bold'], ['red']))
        safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))
        return

    mink_str = safe_input(fmttxt(['Minimum k-mer length', '[default: 4]: '],
                                  ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    min_k = int(mink_str) if mink_str.isdigit() and int(mink_str) >= 2 else 4

    alpha_str = safe_input(fmttxt(['Significance threshold alpha', '[default: 0.05]: '],
                                   ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    try:
        alpha = float(alpha_str) if alpha_str else 0.05
    except ValueError:
        alpha = 0.05

    base = os.path.splitext(file_path)[0] if file_path else 'sequence'
    default_csv = f'{base}_decompose_{motif}.csv'
    fn_str = safe_input(fmttxt(['Output CSV file',
                                 f'[default: {os.path.basename(default_csv)}]: '],
                                ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    out_csv = fn_str if fn_str else default_csv

    print(fmttxt([f'Decomposing {motif!r} …'], [''], ['cyan']))
    results = decompose_motif(txt, motif, alpha=alpha, min_k=min_k)

    if not results:
        print(fmttxt(['Motif too short to decompose (need length >= 2).'], ['bold'], ['yellow']))
        safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))
        return

    import csv as _csv
    try:
        with open(out_csv, 'w', newline='') as fh:
            w = _csv.writer(fh)
            w.writerow([
                'level', 'kmer', 'real_count', 'expected_count',
                'pvalue_over', 'pvalue_under', 'pvalue_bh', 'direction', 'significant',
            ])
            for row in results:
                w.writerow([
                    row['level'], row['kmer'], row['real_count'],
                    f"{row['expected_count']:.3f}",
                    f"{row['pvalue_over']:.8f}",
                    f"{row['pvalue_under']:.8f}",
                    f"{row['pvalue_bh']:.8f}",
                    row['direction'], row['significant'],
                ])
        sig = [r for r in results if r['significant']]
        if sig:
            min_level = min(r['level'] for r in sig)
            min_kmers = [r['kmer'] for r in sig if r['level'] == min_level]
            print(fmttxt([f'Shortest significant unit: length {min_level}:',
                          ', '.join(min_kmers)],
                         ['bold', ''], ['green', 'cyan']))
        else:
            print(fmttxt(['No sub-unit reached significance at alpha =',
                          f'{alpha}'], ['bold', ''], ['yellow', 'cyan']))
        print(fmttxt([f'Saved {len(results)} entries to:', out_csv],
                     ['bold', ''], ['green', 'cyan']))
        open_file_with_default_software(out_csv)
    except Exception as e:
        print(fmttxt([f'Error writing file: {e}'], ['bold'], ['red']))

    safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))


def scramble_input(txt, file_path=''):
    """Prompt for k-mer length and shuffle parameters, then write a CSV of
    all observed k-mers with two-sided empirical p-values sorted ascending.

    Each shuffle preserves the nucleotide composition of the sequence.
    Two one-sided p-values are reported for each k-mer:
      pvalue_over  — fraction of shuffles with count >= real count
                     (small = over-represented)
      pvalue_under — fraction of shuffles with count <= real count
                     (small = under-represented)
    K-mers with no enrichment or depletion will have both p-values near 0.5.
    BH-adjusted p-values (FDR) and E-values (p * m, where m is the number of
    distinct k-mers tested) are provided separately for each direction.
    Rows are sorted by min(pvalue_over, pvalue_under).
    """
    k_str = safe_input(fmttxt(['k-mer length', '[default: 6]: '],
                               ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    k = int(k_str) if k_str.isdigit() and int(k_str) > 0 else 6

    n_str = safe_input(fmttxt(['Number of shuffles', '[default: 5000]: '],
                               ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    N = int(n_str) if n_str.isdigit() and int(n_str) > 0 else 5000

    seed_str = safe_input(fmttxt(['Random seed', '[default: 0]: '],
                                  ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    seed = int(seed_str) if seed_str.lstrip('-').isdigit() else 0

    base = os.path.splitext(file_path)[0] if file_path else 'sequence'
    default_csv = f'{base}_kmer{k}_pvalues.csv'
    fn_str = safe_input(fmttxt(['Output CSV file',
                                 f'[default: {os.path.basename(default_csv)}]: '],
                                ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    out_csv = fn_str if fn_str else default_csv

    print(fmttxt([f'Running {N} shuffles with k={k} …'], [''], ['cyan']))
    results = scramble_kmer_pvalues(txt, k=k, N=N, seed=seed)

    import csv as _csv
    try:
        with open(out_csv, 'w', newline='') as fh:
            w = _csv.writer(fh)
            w.writerow([
                'kmer', 'real_count', 'exceed_count', 'below_count',
                'pvalue_over', 'evalue_over', 'pvalue_over_bh',
                'pvalue_under', 'evalue_under', 'pvalue_under_bh',
                'direction',
            ])
            for row in results:
                w.writerow([
                    row['kmer'], row['real_count'],
                    row['exceed_count'], row['below_count'],
                    f"{row['pvalue_over']:.8f}",  f"{row['evalue_over']:.4f}",  f"{row['pvalue_over_bh']:.8f}",
                    f"{row['pvalue_under']:.8f}", f"{row['evalue_under']:.4f}", f"{row['pvalue_under_bh']:.8f}",
                    row['direction'],
                ])
        print(fmttxt([f'Saved {len(results)} k-mers to:', out_csv],
                     ['bold', ''], ['green', 'cyan']))
        open_file_with_default_software(out_csv)
    except Exception as e:
        print(fmttxt([f'Error writing file: {e}'], ['bold'], ['red']))

    safe_input(fmttxt(['Press Enter to continue'], [''], ['white']))


def sequence_hits_input(fn, txt):
    """Prompt for up to 3 sequences + mutation rate and plot their occurrences.

    For each sequence, finds all same-length exact and approximate (Hamming)
    matches in the text within the requested mutation rate.
    """
    seqs = []
    for i in range(3):
        prompt = fmttxt([f'Enter sequence {i+1}/3', '[blank to stop]: '],
                        ['bold', ''], ['yellow', 'cyan'])
        seq = safe_input(prompt + ' ').strip().lower()
        if not seq:
            break
        seqs.append(seq)

    if not seqs:
        return

    mutr_str = safe_input(fmttxt(['Mutation rate: 1 per N letters', '[default: 6]: '],
                                   ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
    mutr = 1 / (float(mutr_str) if mutr_str else 6)

    fn_out, scale = _prompt_save(plotly=True)

    cnd = safe_input(fmttxt(['Condense x-axis (skip large empty regions)?', '[y/N]: '],
                             ['bold', ''], ['yellow', 'cyan']) + ' ')
    condense = cnd.strip().lower() in ('y', 'yes')
    min_gap = 1000
    if condense:
        mg = safe_input(fmttxt(['Minimum gap size to compress', '[default: 1000]: '],
                                ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
        if mg:
            min_gap = int(mg)

    hairpins = []
    hairpins_csv = os.path.splitext(fn)[0] + '_hairpins.csv'
    if os.path.isfile(hairpins_csv):
        ans = safe_input(fmttxt(['Show hairpin regions?', '[y/N]'],
                                ['bold', ''], ['yellow', 'cyan']) + ' ')
        if ans.strip().lower() == 'y':
            import csv
            with open(hairpins_csv, newline='') as _f:
                for row in csv.DictReader(_f):
                    hairpins.append({
                        'start':    int(row['start']),
                        'end':      int(row['end']),
                        'stem_seq': row.get('stem_seq', ''),
                        'loop_seq': row.get('loop_seq', ''),
                        'sequence': row.get('sequence', ''),
                    })

    seq_data = []
    for seq in seqs:
        # M=len(seq) removes the hard cap so maxmut is determined purely by mutr
        exact_pos, approx, _ = find_with_mutations(seq, txt, mutr=mutr, M=len(seq))
        seq_data.append({'seq': seq, 'L': len(seq), 'exact': exact_pos, 'approx': approx})

    title = os.path.basename(fn)
    plot_sequence_hits(seq_data, len(txt), title=title, file=fn_out, scale=scale,
                       condense=condense, min_gap=min_gap, hairpins=hairpins)

    # ── Optional detailed text plot ───────────────────────────────────────
    det = safe_input(fmttxt(['Create detailed nucleotide-level plot?', '[y/N]: '],
                             ['bold', ''], ['yellow', 'cyan']) + ' ')
    if det.strip().lower() in ('y', 'yes'):
        rng_str = safe_input(fmttxt(['Position range',
                                      f'[start,end  default: 0,{len(txt)}]: '],
                                     ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
        if rng_str:
            try:
                parts     = rng_str.split(',')
                det_start = int(parts[0].strip())
                det_end   = int(parts[1].strip())
            except (IndexError, ValueError):
                det_start, det_end = 0, len(txt)
        else:
            det_start, det_end = 0, len(txt)
        if det_end - det_start > 500:
            print(fmttxt([f'Note: range {det_end - det_start} nt — '
                           'large ranges may render slowly in a browser.'],
                         [''], ['yellow']))
        det_file = safe_input(fmttxt(['Output HTML file',
                                       f'[default: {os.path.splitext(os.path.basename(fn))[0]}_detail.html]: '],
                                      ['bold', ''], ['yellow', 'cyan']) + ' ').strip()
        if not det_file:
            det_file = os.path.splitext(fn)[0] + '_detail.html'
        elif not det_file.lower().endswith('.html'):
            det_file += '.html'
        plot_sequence_hits_detailed(seq_data, txt, det_start, det_end,
                                    title=title, file=det_file)



def print_stats(txt, strs):
    '''Print summary statistics for the loaded text and wait for the user to press Enter.'''
    xmotifs  = strs['xmotifs']
    corelist = strs['corelist']
    print(fmttxt(["\n--- Summary Statistics ---"], ['bold'], ['cyan']))
    rows = [
        ('Text length',              len(txt)),
        ('Number of xmotifs',        len(xmotifs)),
        ('Longest xmotif',           max((len(x) for x in xmotifs),  default=0)),
        ('Shortest xmotif',          min((len(x) for x in xmotifs),  default=0)),
        ('Number of cores',          len(corelist)),
        ('Longest core',             max((len(x) for x in corelist), default=0)),
        ('Shortest core',            min((len(x) for x in corelist), default=0)),
    ]
    for lbl, val in rows:
        print(fmttxt([f"  {lbl}:", str(val)], ['bold', ''], ['white', 'cyan']))
    print()
    safe_input(fmttxt(["Press Enter to return to the menu"], [''], ['white']))


def print_settings(defvals, prefs=None):
    '''Print current default values and wait for the user to press Enter.'''
    print(fmttxt(["\n--- Current Settings ---"], ['bold'], ['cyan']))
    labels = {
        'minxmlen':   'Min xmotif length',
        'maxxmlen':   'Max xmotif length',
        'mincorelen': 'Min core length',
        'mincount':   'Min occurrences',
        'pwr':        'Coverage weight (pwr)',
        'clr':        'Clear screen',
        'datadir':    'Session data directory',
    }
    for key, lbl in labels.items():
        print(fmttxt([f"  {lbl}:", str(defvals[key])], ['bold', ''], ['white', 'cyan']))
    if prefs is not None:
        print(fmttxt(["\n--- Persistent Preferences ---"], ['bold'], ['cyan']))
        print(fmttxt(["  Default data directory:",
                      prefs.get('default_data_dir') or '(not set)'],
                     ['bold', ''], ['white', 'cyan']))
        print(fmttxt(["  Last used directory:",
                      prefs.get('last_used_dir') or '(not set)'],
                     ['bold', ''], ['white', 'cyan']))
    print()
    safe_input(fmttxt(["Press Enter to return to the menu"], [''], ['white']))


def get_choices(defvals, prefs=None):
    '''Prompt the user to modify the default values for global parameters.'''
    reload = False
    minxmlen = safe_input(fmttxt(['Minimum length for xmotifs', f'[{defvals["minxmlen"]}]: '],
                                  ['bold',''], ['yellow','cyan']))
    if minxmlen == '':
        minxmlen = defvals["minxmlen"]
    else:
        reload = True
    maxxmlen = safe_input(fmttxt(['Maximum length for xmotifs', f'[{defvals["maxxmlen"]}]: '],
                                  ['bold',''], ['yellow','cyan']))
    if maxxmlen == '':
        maxxmlen = defvals["maxxmlen"]
    else:
        reload = True
    mincorelen = safe_input(fmttxt(['Minimum length for cores', f'[{defvals["mincorelen"]}]: '],
                                  ['bold',''], ['yellow','cyan']))
    if mincorelen == '':
        mincorelen = defvals["mincorelen"]
    else:
        reload = True
    mincount = safe_input(fmttxt(['Minimum number of occurrences', f'[{defvals["mincount"]}]: '],
                                  ['bold',''], ['yellow','cyan']))
    if mincount == '':
        mincount = defvals["mincount"]
    else:
        reload = True
    pwr = safe_input(fmttxt(['Coverage weight parameter', f'[{defvals["pwr"]}]: '],
                                  ['bold',''], ['yellow','cyan']))
    if pwr == '':
        pwr = defvals["pwr"]
    else:
        reload = True
    clr = safe_input(fmttxt(['Clear screen after input (True/False)', f'[{defvals["clr"]}]: '],
                                  ['bold',''], ['yellow','cyan']))
    if clr == '':
        clr = defvals['clr']
    else:
        if clr.startswith('F'):
            clr = False
        else:
            clr = True
    datadir = safe_input(fmttxt(['Session data directory (Enter to browse, - to clear)',
                                  f'[{defvals["datadir"] or "not set"}]: '],
                                 ['bold', ''], ['yellow', 'cyan']))
    if datadir == '':
        datadir = defvals['datadir']
        if not datadir:
            print(fmttxt(["Browsing for data directory..."], [''], ['cyan']))
            datadir = openDir(initial_dir=None)
    elif datadir == '-':
        datadir = ''

    # Persistent default data directory (saved across sessions).
    cur_default = (prefs or {}).get('default_data_dir') or ''
    default_data_dir_input = safe_input(fmttxt(
        ['Default data directory (Enter to keep, - to clear, B to browse)',
         f'[{cur_default or "not set"}]: '],
        ['bold', ''], ['yellow', 'cyan']))
    if default_data_dir_input == '':
        default_data_dir = cur_default
    elif default_data_dir_input == '-':
        default_data_dir = ''
    elif default_data_dir_input.strip().upper() == 'B':
        print(fmttxt(["Browsing for default data directory..."], [''], ['cyan']))
        default_data_dir = openDir(initial_dir=cur_default or None) or cur_default
    else:
        default_data_dir = default_data_dir_input.strip()

    return({'minxmlen': int(minxmlen), 'maxxmlen': int(maxxmlen), 'mincorelen': int(mincorelen),
            'mincount': int(mincount), 'pwr': float(pwr), 'clr': bool(clr),
            'datadir': datadir, 'reload': reload, 'default_data_dir': default_data_dir})


def parsedata(txt, defvals):
    """Discover xmotifs and cores from a sequence using the current default settings.

    Runs find_boundary() to collect boundary-delimited recurring substrings,
    then cores() to extract the dominant cores.  Extended motifs that do not
    contain any core are appended as ``morecores`` and merged into the final
    core list, which is sorted by weighted coverage.

    If the longest xmotif found equals ``defvals['maxxmlen']``, the maximum is
    raised by 20 and the search is repeated until the longest result is strictly
    shorter than the ceiling.  ``defvals['maxxmlen']`` is updated in-place so
    the expanded value is reflected in the session settings.

    Args:
        txt:     Source sequence or text.
        defvals: Global defaults dict with keys ``'minxmlen'``, ``'maxxmlen'``,
                 ``'mincorelen'``, and ``'mincount'``.  ``'maxxmlen'`` may be
                 increased in-place if auto-expansion is triggered.

    Returns:
        Dict with keys ``'corelist'`` (combined cores + morecores, sorted by
        coverage) and ``'xmotifs'`` (all boundary motifs, sorted by coverage).
    """
    while True:
        xmotifs = find_boundary(txt, minlen=defvals['minxmlen'], maxlen=defvals['maxxmlen'],
                                at_least=defvals['mincount'])
        if not xmotifs:
            break
        longest = max(len(x) for x in xmotifs)
        if longest < defvals['maxxmlen']:
            break
        # The longest xmotif hits the ceiling — expand and re-run.
        defvals['maxxmlen'] += 20
        print(fmttxt([f"Longest xmotif ({longest} nt) equals the maximum length allowed.",
                      f" Increasing maximum to {defvals['maxxmlen']} and re-analysing…"],
                     ['bold', ''], ['yellow', 'cyan']))
    xmotifs.sort(key=lambda s: cover(s, txt), reverse=True)
    # Find the cores from the xmotifs and save in a list called corelist:
    corelist = cores(txt, xmotifs, minclen=defvals['mincorelen'])
    # sort them in order of area covered
    corelist.sort(key=lambda s: cover(s, txt), reverse=True)
    # If an xmotif did not contain a core, it will not be in corelist. 
    # we find those, and put them in morecores, and combine them into
    # core_xm
    morecores = []
    for i in range(len(xmotifs)):
        incore = False
        for j in range(len(corelist)):
            if corelist[j] in xmotifs[i]:
                incore = True
                break
        if not incore:
            morecores.append(xmotifs[i])

    morecores.sort(key=lambda s: cover(s, txt), reverse=True)
    core_xm = corelist + morecores
    core_xm.sort(key=lambda s: cover(s, txt), reverse=True)

    return({'corelist': core_xm, 'xmotifs': xmotifs})


def menus():
    """Run the interactive RNA_explore menu loop.

    Initialises default settings and the three-level menu structure (Main,
    Plots, Sequence operations).  On startup, scans the current working
    directory (or a user-chosen directory) for valid session JSON files and
    offers to load one automatically.

    The loop handles:
    * Session loading (file, Ensembl ENST, paste) and saving.
    * All plot actions (neighbours, k-mers, logo, coverage, sequence hits).
    * All sequence operations (search, alignment, hairpin export, extensions).
    * Settings changes with optional re-parsing of xmotifs/cores.
    * Ctrl+C exits cleanly; Ctrl+D resets to the main menu.
    """
    pid = getpid()
    prefs = load_prefs()
    _initial_dir = prefs.get('last_used_dir') or prefs.get('default_data_dir') or ''
    defvals = {'minxmlen': 7, 'maxxmlen': 60, 'mincorelen': 6, 'mincount': 2, 'pwr': 1.2,
               'clr': True, 'reload': False, 'datadir': _initial_dir}
    menu_ttl = ["Main Menu", "Plots", "Sequence operations"]
    submenu = [["Plots", "Sequence operations", "Open Core file", "Summary statistics",
               "Show settings", "Change setting", "Open User Guide", "Load new input", "Quit"],
               ["Core neighbors", "K-mers", "Logo", "Coverage", "Motif Match/Mutation", "Back"],
               ["Find all matches", "Search with mutations", "Motif extensions", "Print core", "Export hairpins to CSV", "Extend match pair", "Alignment score for two sequences", "K-mer scramble analysis", "Covered area", "Decompose motif", "Back"]]
    if not defvals['datadir']:
        cwd = os.getcwd()
        valid = _find_valid_sessions(cwd)
        if valid:
            defvals['datadir'] = cwd
        else:
            print(fmttxt(["No data directory set. Please select one."], ['bold'], ['yellow']))
            chosen = openDir(initial_dir=None)
            if chosen:
                defvals['datadir'] = chosen
                valid = _find_valid_sessions(chosen)
        if valid:
            if len(valid) == 1:
                globals()['fn'] = valid[0][1]
            else:
                dir_display = os.path.basename(defvals['datadir']) or defvals['datadir']
                idx = show_menu(dir_display, "Select a session to load",
                                [name for name, _ in valid], clr=False)
                if idx > 0:
                    globals()['fn'] = valid[idx - 1][1]
    menu_level = 0
    
    while True:
        try:
            while True:
                try:
                    txt
                except NameError:
                    if 'fn' in globals().keys():
                        data = load_session(globals()['fn'])
                        if data is None:
                            # Session file missing or corrupt — clear fn and
                            # fall through to the input-source menu.
                            del globals()['fn']
                            continue
                        file_path = data['file_path']
                        txt = data['txt']
                        txtb = data['txtb']
                        is_rna = data['is_rna']
                        if data.get("corelist"):
                            strs = {'corelist': data["corelist"], 'xmotifs': data["xmotifs"]}
                        else:
                            strs = parsedata(txt, defvals)
                    # If no text file was loaded, prompt the user to select a file
                    else:
                        data = choose_input_source(defvals['datadir'])
                        if data == -1:
                            return None
                        if data is None:
                            continue  # dialog cancelled – show the menu again
                        file_path = data['file_path']
                        txt = data['txt']
                        txtb = data['txtb']
                        is_rna = data['is_rna']
                        if data.get("corelist"):
                            strs = {'corelist': data["corelist"], 'xmotifs': data["xmotifs"]}
                        else:
                            strs = parsedata(txt, defvals)
                    workdir = data['dir']
                    if workdir and os.path.isdir(workdir):
                        chdir(workdir)
                        defvals['datadir'] = workdir
                    if defvals['datadir']:
                        prefs['last_used_dir'] = defvals['datadir']
                        save_prefs(prefs)
                    globals()['fn'] = file_path + ".json"
                    stats = {
                        'txt_len':      len(txt),
                        'xm_count':     len(strs['xmotifs']),
                        'xm_max_len':   max((len(x) for x in strs['xmotifs']), default=0),
                        'core_count':   len(strs['corelist']),
                        'core_max_len': max((len(x) for x in strs['corelist']), default=0),
                    }
                    save_session(file_path, data={'file_path': file_path,
                                                  'txt': txt,
                                                  'txtb': txtb,
                                                  'is_rna': is_rna,
                                                  'corelist': strs['corelist'],
                                                  'xmotifs': strs['xmotifs'],
                                                  'dir': workdir,
                                                  'stats': stats})
                    init_summary(file_path, strs['xmotifs'], strs['corelist'], txt)
                else:
                    val = show_menu(file_path, menu_ttl[menu_level], submenu[menu_level],
                                    clr=defvals['clr'],
                                    split={0: 4, 1: 5, 2: 10}.get(menu_level))
                    # Top level menu:
                    if menu_level == 0:
                        # Quit:
                        if val == len(submenu[menu_level]) or val <= 0:
                            return
                        if val == 1: # Change to plots menu
                            menu_level = 1
                            continue
                        if val == 2: # Change to sequence menu
                            menu_level = 2
                            continue
                        # show the core file (CSV)
                        if val == 3:
                            out_csv = f'{file_path}_init.csv'
                            open_file_with_default_software(out_csv)
                        if val == 4:
                            print_stats(txt, strs)
                            continue
                        if val == 5:
                            print_settings(defvals, prefs)
                            continue
                        if val == 6:
                            defvals = get_choices(defvals, prefs)
                            new_default = defvals.pop('default_data_dir', None)
                            if new_default is not None and new_default != prefs.get('default_data_dir'):
                                prefs['default_data_dir'] = new_default
                                save_prefs(prefs)
                            if defvals['reload']:
                                # recalculate xmotifs, cores, update JSON and CSV
                                strs = parsedata(txt, defvals)
                                stats = {
                                    'txt_len':      len(txt),
                                    'xm_count':     len(strs['xmotifs']),
                                    'xm_max_len':   max((len(x) for x in strs['xmotifs']), default=0),
                                    'core_count':   len(strs['corelist']),
                                    'core_max_len': max((len(x) for x in strs['corelist']), default=0),
                                }
                                save_session(file_path, data={'file_path': file_path,
                                                              'txt': txt,
                                                              'txtb': txtb,
                                                              'is_rna': is_rna,
                                                              'corelist': strs['corelist'],
                                                              'xmotifs': strs['xmotifs'],
                                                              'dir': workdir,
                                                              'stats': stats})
                                out_csv = f'{file_path}_init.csv'
                                if os.path.isfile(out_csv):
                                    os.remove(out_csv)
                                init_summary(file_path, strs['xmotifs'], strs['corelist'], txt)
                            continue
                        # Open the user guide in the default browser:
                        if val == 7:
                            _guide_md = os.path.abspath(os.path.join(
                                os.path.dirname(__file__), '..', '..', 'docs', 'user_guide.md'))
                            if os.path.isfile(_guide_md):
                                import markdown as _md
                                with open(_guide_md, 'r') as _f:
                                    _body = _md.markdown(_f.read(), extensions=['tables', 'fenced_code'])
                                _guide_html = os.path.splitext(_guide_md)[0] + '.html'
                                with open(_guide_html, 'w') as _f:
                                    _f.write(f'''<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>RNA_lexis User Guide</title>
<style>
body{{font-family:sans-serif;max-width:960px;margin:2em auto;padding:0 1.5em;line-height:1.6;color:#222;}}
h1,h2,h3,h4{{color:#1a4a7a;}}
h2{{border-bottom:1px solid #ccc;padding-bottom:0.2em;}}
table{{border-collapse:collapse;margin:1em 0;width:100%;}}
th,td{{border:1px solid #bbb;padding:6px 12px;text-align:left;}}
th{{background:#f0f4f8;}}
code{{background:#f4f4f4;padding:2px 5px;border-radius:3px;font-size:0.92em;}}
pre{{background:#f4f4f4;padding:1em;border-radius:5px;overflow-x:auto;}}
pre code{{background:none;padding:0;}}
hr{{border:none;border-top:1px solid #ddd;margin:2em 0;}}
p img{{display:block;margin:0 auto;}}
blockquote{{border-left:4px solid #ccc;margin:1em 0;padding:0.5em 1em;color:#555;}}
</style></head><body>{_body}</body></html>''')
                                webbrowser.open(f'file://{_guide_html}')
                            else:
                                print(fmttxt(["User guide not found:", _guide_md], ['bold', ''], ['red', 'white']))
                        # Load another file, and show the main menu again:
                        if val == 8:
                            del txt
                            del globals()['fn']
                            continue
                    # Plots:
                    if menu_level == 1:
                        match val:
                            case 0:
                                menu_level = 0
                            case 1:
                                neighbors_input(file_path, txt, strs)
                            case 2:
                                kmers_input(file_path, txt, defvals)
                            case 3:
                                logo_input(txt, strs)
                            case 4:
                                coverage_input(txt, strs)
                            case 5:
                                sequence_hits_input(file_path, txt)
                            case 6:
                                menu_level = 0 # Go to the main menu

                    # Sequences:
                    if menu_level == 2:
                        match val:
                            case 0:
                                menu_level = 0
                            case 1:
                                find_match_input(txt, strs)
                            case 2:
                                search_input(txt)
                            case 3:
                                motif_extensions_input(txt)
                            case 4:
                                print_core_input(txt, strs)
                            case 5:
                                hairpins_input(txt, file_path)
                            case 6:
                                extend_match_input(txt)
                            case 7:
                                print_alignment_score(txt)
                            case 8:
                                scramble_input(txt, file_path)
                            case 9:
                                txt_coverage_input(txt, strs)
                            case 10:
                                decompose_motif_input(txt, file_path)
                            case 11:
                                menu_level = 0 # Go to the main menu
        
        except EOFSignal:
            # Ctrl+D pressed: reset to main menu
            menu_level = 0


def main():
    """CLI entrypoint for the interactive RNAlang menu."""
    # Set dark-navy background for contrast with menu colours (OSC 11).
    # Terminals that do not support OSC sequences ignore this silently.
    print('\033]11;#0d1b2a\007', end='', flush=True)
    try:
        menus()
    except KeyboardInterrupt:
        print("\nExiting.")
    finally:
        print('\033]111;\007', end='', flush=True)  # restore user's default background


if __name__ == "__main__":
    main()
