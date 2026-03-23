"""Pure sequence-analysis algorithms — no I/O, no plotting."""

import re
import math
import itertools
import statistics
import numpy as np
from collections import Counter, defaultdict
from typing import List, NamedTuple


def contains_only_rna(s, allowed_chars='atucg'):
    """Return True if every character in s belongs to the allowed nucleotide set.

    Args:
        s:             Input string (case-insensitive; converted to lower-case
                       internally).
        allowed_chars: Characters considered valid nucleotides.  Defaults to
                       ``'atucg'`` so that both DNA (T) and RNA (U) are accepted.

    Returns:
        True when s contains only characters from allowed_chars, False otherwise.
    """
    if not s:
        return False
    s = s.lower()
    return all(char in allowed_chars for char in s)


def count_kgrams(s: str, k: int, skip='_', at_least=2, return_sorted=True, rev=False) -> dict:
    """Count all k-length substrings (k-grams) in s.

    Args:
        s:            Input string.
        k:            Substring length to count.
        skip:         Separator character; any k-gram that contains this
                      character is excluded from the count (default ``'_'``).
        at_least:     Minimum occurrence count for a k-gram to be included in
                      the result (default 2).
        return_sorted: When True (default), the returned dict is sorted by
                       count in descending order (most frequent first).
        rev:          When True, sort in ascending order instead of descending.
                      Has no effect when return_sorted is False.

    Returns:
        Dict mapping k-gram string → count, filtered to entries with count
        >= at_least and optionally sorted by count.
    """
    counter = Counter(
        s[j:j+k] for j in range(len(s) - k + 1)
        if skip not in s[j:j+k]
    )
    kgrams = {key: val for key, val in counter.items() if val >= at_least}
    if return_sorted:
        kgrams = {k: v for k, v in sorted(kgrams.items(), key=lambda item: item[1], reverse=rev)}
    return(kgrams)


def find_all_matches(pattern: str, txt: str, ret='pos', group=0) -> list:
    """Find all (possibly overlapping) occurrences of a regex pattern in txt.

    The search is case-insensitive: both pattern and txt are lower-cased before
    matching.  Because the position advances by one after each match, overlapping
    occurrences are all reported.

    Args:
        pattern: Regular expression pattern to search for.
        txt:     Source text to search in.
        ret:     ``'pos'`` (default) to return a list of start positions;
                 ``'str'`` to return the list of matched strings.
        group:   Capture group index to use when ret='str' (default 0 = whole
                 match).

    Returns:
        List of int start positions (ret='pos') or list of matched strings
        (ret='str').
    """
    if not pattern or not txt:
        return []
    pat = re.compile(pattern.lower())
    pos = 0
    out = []
    while m := pat.search(txt.lower(), pos):
        pos = m.start() + 1
        if ret == 'pos':
            out.append(m.start())
        else:
            out.append(m[group])
    return out


def stitch(txt: str, txtsp: str, maxlen = 2, min_prob = 0.5, min_odr = 2) -> str:
    """Merge consecutive short words in a space-separated string when evidence
    supports treating them as a single token.

    Words in txtsp are separated by ``'_'``.  For each short interior word
    (length <= maxlen), the function estimates whether it is more likely to
    combine with the word on its left or on its right by comparing empirical
    co-occurrence probabilities against txt.  If one combination is
    strongly preferred (probability > min_prob and the two options differ), the
    short word is absorbed into its preferred neighbour.

    Args:
        txt:      The full original text used as a reference corpus.
        txtsp:    A ``'_'``-separated string of candidate tokens.
        maxlen:   Maximum length of a word to be considered for merging
                  (default 2).
        min_prob: Minimum co-occurrence probability required for a merge to
                  occur (default 0.5).
        min_odr:  Reserved parameter (currently unused).

    Returns:
        Updated ``'_'``-separated string with short words merged where
        appropriate.
    """
    txtsp = txtsp.rstrip('_').lower()
    arr = txtsp.split('_')
    for i in range(1, len(arr)-1):
        if len(arr[i]) > maxlen:
            continue
        if min(len(arr[i-1]), len(arr[i+1])) <= maxlen:
            continue
        p1 = (2+len(find_all_matches(arr[i-1]+arr[i], txt)))/ \
            (4+len(find_all_matches(arr[i-1], txt, ret='str')))
        p2 = (2+len(find_all_matches(arr[i]+arr[i+1], txt)))/ \
            (4+len(find_all_matches(arr[i+1], txt, ret='str')))
        if abs(p1-p2) < 1e-2:
            continue
        if max(p1, p2) <= min_prob:
            continue
        if min(p1, p2) > 0.4:
            continue
        if p1 > p2:
            arr[i-1] = arr[i-1] + arr[i]
            arr[i] = ''
        else:
            arr[i+1] = arr[i] + arr[i+1]
            arr[i] = ''
    return('_'.join(arr))


def find_max_cover(txt, words, mincover=50):
    """Return the word from words that achieves the highest weighted coverage
    in txt, or an empty string if no word exceeds mincover.

    Coverage is computed by cover() which weights both the number of
    non-overlapping occurrences and the word length.

    Args:
        txt:       Source text to search in.
        words:     Candidate words (list of strings).
        mincover:  Minimum coverage score a word must achieve to be returned
                   (default 50).

    Returns:
        The highest-scoring word string, or ``''`` if no word meets the
        mincover threshold.
    """
    area = dict()
    maxval = 0
    maxstr = ''
    for i in range(len(words)):
        #area[words[i]] = len(find_all_matches(words[i], txt))*len(words[i])
        area[words[i]] = cover(words[i], txt, pwr=1.2)
        if area[words[i]] > maxval:
            maxstr = words[i]
            maxval = area[words[i]]
    if maxval < mincover:
        return('')
    return(maxstr)


def split_words(words, maxstr, txt, minlen=5):
    """Split extended motifs (xmotifs) around a core sequence.

    For each word in words that contains maxstr as a substring, the prefix
    and suffix flanking maxstr are extracted; those that are long enough
    (>= minlen) are kept and further expanded to their natural boundary via
    expand_to_boundary.  Words that do not contain maxstr are kept unchanged.

    Args:
        words:  List of extended motif strings.
        maxstr: The current dominant core sequence to split around.
        txt:    Full source text, passed to expand_to_boundary.
        minlen: Minimum length of a flanking fragment to retain (default 5).

    Returns:
        Set of surviving/derived word strings after the split.
    """
    newwords = set()
    for i in range(len(words)):
        if words[i] == maxstr:
            continue
        if maxstr in words[i]:
            prefpost = words[i].split(maxstr)
            if len(prefpost[0]) >= minlen:
                newwords.add(expand_to_boundary(prefpost[0], txt))
            if len(prefpost[1]) >= minlen:
                newwords.add(expand_to_boundary(prefpost[1], txt))
        else:
            newwords.add(words[i])
    return(newwords)


def cover(s, txt, pwr=1.1):
    """Compute the weighted coverage score of string s in txt.

    The score rewards both frequency and length: it is proportional to
    len(s)**pwr * (1 + number_of_non_overlapping_occurrences).  Strings
    with fewer than 2 non-overlapping occurrences score 0.

    Args:
        s:   Query string.
        txt: Source text.
        pwr: Length exponent; values > 1 favour longer strings (default 1.1).

    Returns:
        Float coverage score (0 when s appears fewer than twice).
    """
    found = find_all_matches(s, txt)
    if len(found) < 2:
        return(0)
    ls = len(s)
    cnt = 0
    s0 = [0]*len(txt)
    s0[found[0]:found[0]+ls] = [1]*ls
    for i in range(1, len(found)):
        if max(s0[found[i]:found[i]+ls]) == 1:
            continue # no overlap
        s0[found[i]:found[i]+ls] = [1]*ls
        cnt += 1
    return(ls**pwr * (1+cnt))

    #return(len(find_all_matches(s, txt))*len(s))


def cores(txt, xmotifs, minclen=5):
    """Iteratively extract core sequences from a list of extended motifs.

    At each iteration the xmotif with the highest weighted coverage (see
    cover()) is selected as the next core.  Its occurrences are masked in a
    working copy of txt (replaced with ``'_'``), and the remaining xmotifs
    are split around the chosen core so that only the flanking sub-sequences
    survive for the next round.  Iteration stops when no candidate achieves
    the minimum coverage threshold.

    Args:
        txt:      Source text.
        xmotifs:  List of extended motif strings (candidates for cores).
        minclen:  Minimum length for a flanking fragment to survive a split
                  (default 5).

    Returns:
        List of core strings ordered by the iteration in which they were
        selected (most dominant first).
    """
    xmotifs.sort(key=lambda s: len(s), reverse=False)
    txt2 = txt
    corelist = []
    maxstr = '-'
    mncvr = max(math.ceil(len(txt)**(1/3)), 13)
    while maxstr != '':
        maxstr = find_max_cover(txt2, xmotifs, mincover = mncvr)
        if maxstr == '':
            break
        corelist.append(maxstr)
        txt2 = txt2.replace(maxstr, '_')
        xmotifs = list(split_words(xmotifs, maxstr, txt, minclen))
    return(corelist)


def print_core(txt, core, xm):
    """Return an aligned, count-annotated list of xmotifs that contain core.

    All xmotifs in xm that contain core as a substring are padded on the left
    so that the core aligns vertically across all entries.  Each entry is
    prefixed with its occurrence count in txt.

    Args:
        txt:  Source text (used to count occurrences of each xmotif).
        core: Core sequence to align around (case-insensitive).
        xm:   List of extended motif strings to filter and align.

    Returns:
        Sorted list of strings, each of the form ``' <count> <padded_xmotif>'``.
    """
    core = core.lower()
    txt = txt.lower()
    xmmatch = dict()
    maxLeft = 0
    for i in range(len(xm)):
        xm[i] = xm[i].lower()
        if core in xm[i]:
            xmmatch[xm[i]] = find_all_matches(core, xm[i])[0]
            maxLeft = max(maxLeft, xmmatch[xm[i]])
    xmpadded = []
    for i in range(len(xmmatch)):
        adjkey = list(xmmatch.keys())[i]
        adjkey = adjkey.rjust(maxLeft - xmmatch[adjkey] + len(adjkey), ' ')
        xmpadded.append(adjkey)
    if len(xmpadded) > 1:
        for k in range(len(xmpadded)):
            cnt = len(find_all_matches(xmpadded[k].rstrip(' ').lstrip(' '), txt))
            xmpadded[k] = f'{cnt:6d} {xmpadded[k]}'
    return(sorted(xmpadded))


def is_bounded(s, txt, positions=None):
    """Return True if s is a maximal word — bounded on both left and right.

    A string is considered bounded when:
    1. Its set of positions in txt equals the set of positions of its interior
       (s[1:-1]) when the surrounding characters are included, i.e. no
       position of s[1:-1] is flanked by the same pair of characters as
       always occurs around s.
    2. The characters immediately to the left of all occurrences of s are not
       all identical (the left boundary is varied).
    3. The characters immediately to the right of all occurrences of s are not
       all identical (the right boundary is varied).

    Args:
        s:         Query string (case-insensitive).
        txt:       Source text (case-insensitive).
        positions: Pre-computed list of start positions of s in txt.  If None,
                   positions are computed internally (useful for caching).

    Returns:
        True when s is bounded on both sides, False otherwise.
    """
    s = s.lower()
    txt = txt.lower()

    # Get (or reuse) positions of s — str.find is faster than regex for literal strings
    if positions is None:
        positions = []
        start = 0
        while (idx := txt.find(s, start)) != -1:
            positions.append(idx)
            start = idx + 1
    if not positions:
        return False

    S0 = set(positions)
    n = len(txt)
    slen = len(s)
    s1 = s[1:-1]
    s1len = len(s1)

    # S1: set of positions one-before each occurrence of s1 (equivalent to '.{s1}.' match starts)
    # str.find replaces regex since s1 is a literal string
    S1 = set()
    start = 0
    while (idx := txt.find(s1, start)) != -1:
        if idx > 0 and idx + s1len < n:
            S1.add(idx - 1)
        start = idx + 1
    if S0 != S1:
        return False

    # SL / SR: derived directly from known positions — no text scan needed
    left_chars = {txt[pos - 1] for pos in positions if pos > 0}
    if len(left_chars) == 1:
        return False
    right_chars = {txt[pos + slen] for pos in positions if pos + slen < n}
    if len(right_chars) == 1:
        return False

    return True


def find_boundary(txt, minlen=7, maxlen=40, at_least=3):
    """Find all boundary-delimited recurring substrings in txt.

    For each length j from minlen to maxlen, collects all substrings of
    length j that appear at least at_least times and pass the is_bounded()
    test.  The scan stops early at the first length where no substring meets
    the at_least threshold.

    Args:
        txt:      Source text.
        minlen:   Minimum substring length to consider (default 7).
        maxlen:   Maximum substring length to consider (default 40).
        at_least: Minimum occurrence count required (default 3).

    Returns:
        List of bounded substring strings (may contain duplicates across
        lengths — use set() to deduplicate if needed).
    """
    txt_lower = txt.lower()
    n = len(txt_lower)
    words = []
    for j in range(minlen, maxlen+1):
        # Single pass: collect positions per k-gram (replaces count_kgrams + find_all_matches)
        pos_dict = {}
        for i in range(n - j + 1):
            sg = txt_lower[i:i+j]
            if '_' in sg:
                continue
            if sg not in pos_dict:
                pos_dict[sg] = []
            pos_dict[sg].append(i)
        qualifying = {wd: pos for wd, pos in pos_dict.items() if len(pos) >= at_least}
        if not qualifying:
            print('Maximum word length', j-1)
            break
        for wd, positions in qualifying.items():
            if is_bounded(wd, txt_lower, positions=positions):
                words.append(wd)
    return(words)


def expand_to_boundary(wd, txt, dir='b'):
    """Extend wd left and/or right until it becomes a bounded word in txt.

    Starting from a seed string (typically a core), characters are appended
    one at a time in the chosen direction(s) as long as all occurrences share
    the same flanking character, growing the word until the boundary condition
    is met.

    Args:
        wd:  Seed string to expand (case-insensitive).
        txt: Source text.
        dir: Direction of expansion — ``'b'`` (default) expands both sides,
             ``'L'`` expands left only, ``'R'`` expands right only.

    Returns:
        The expanded string, or None if wd occurs fewer than twice in txt.
    """
    wd = wd.lower()
    txt_lower = txt.lower()

    # Find all positions once (avoids repeated regex searches)
    positions = find_all_matches(wd, txt_lower, ret='pos')

    if len(positions) <= 1:
        return None

    # Check if already bounded using cached positions
    if _is_bounded_fast(wd, txt_lower, positions):
        return wd

    # Expand left: extend as long as the left-boundary character is unique
    if dir != 'R':
        wd, positions = _expand_left(wd, txt_lower, positions)

    # Expand right: extend as long as the right-boundary character is unique
    if dir != 'L':
        wd, positions = _expand_right(wd, txt_lower, positions)

    return wd


def _is_bounded_fast(wd, txt, positions):
    """Return True if wd is bounded on both left and right given pre-computed positions.

    Faster than is_bounded() because it skips the S0/S1 equality test and
    relies entirely on boundary-character diversity.  Used internally by
    expand_to_boundary after every expansion step.
    """
    wd_len = len(wd)
    left_chars = set()
    right_chars = set()

    # Collect unique boundary characters from all positions
    for pos in positions:
        if pos > 0:
            left_chars.add(txt[pos - 1])
        if pos + wd_len < len(txt):
            right_chars.add(txt[pos + wd_len])

    # Bounded if we have multiple unique boundary characters on both sides
    return len(left_chars) > 1 and len(right_chars) > 1


def _expand_left(wd, txt, positions):
    """Greedily extend wd one character to the left as long as all occurrences
    share the same left-flanking character and at least two copies exist.

    Returns:
        Tuple (expanded_wd, updated_positions) after all possible left
        extensions.
    """
    while True:
        left_extensions = {}

        # Collect left-extended versions for each position
        for pos in positions:
            if pos > 0:
                extended = txt[pos - 1] + wd
                left_extensions[extended] = left_extensions.get(extended, 0) + 1

        # Stop if there isn't exactly one unique extension or if count < 2
        if len(left_extensions) != 1:
            break

        extended, count = next(iter(left_extensions.items()))
        if count < 2:
            break

        wd = extended
        # Adjust positions: move left by 1
        positions = [p - 1 for p in positions]

    return wd, positions


def _expand_right(wd, txt, positions):
    """Greedily extend wd one character to the right as long as all occurrences
    share the same right-flanking character and at least two copies exist.

    Returns:
        Tuple (expanded_wd, positions) — positions are unchanged because
        right extension does not affect start coordinates.
    """
    while True:
        right_extensions = {}
        wd_len = len(wd)

        # Collect right-extended versions for each position
        for pos in positions:
            if pos + wd_len < len(txt):
                extended = wd + txt[pos + wd_len]
                right_extensions[extended] = right_extensions.get(extended, 0) + 1

        # Stop if there isn't exactly one unique extension or if count < 2
        if len(right_extensions) != 1:
            break

        extended, count = next(iter(right_extensions.items()))
        if count < 2:
            break

        wd = extended
        # Positions stay the same when extending right

    return wd, positions


def add_mut(s0, idx):
    """Upper-case the character at position idx in s0 to mark it as a mutation site.

    Used internally by allow_mutation() to flag positions that will later be
    converted to ``[^x]`` wildcard groups in a regex.

    Args:
        s0:  Input string (lower-case).
        idx: 0-based index of the position to mark.

    Returns:
        Copy of s0 with s0[idx] replaced by its upper-case equivalent.
    """
#    if not (0 <= idx < len(s0)):
#        raise ValueError("Invalid start or end index for string replacement.")
    return s0[:idx] + s0[idx].upper() + s0[idx+1:]


def allow_mutation(s0, n=1):
    """Build a regex that matches s0 with up to n single-character substitutions.

    Each combination of n positions in s0 is turned into a pattern where each
    chosen position is replaced by ``[^x]`` (any character except the original),
    and all resulting patterns are joined with ``|``.  The number of mutations is
    capped at floor(len(s0)/6).  Strings of length <= 6 are returned unchanged.

    Args:
        s0: Query string (lower-case).
        n:  Maximum number of substitution positions allowed (default 1).

    Returns:
        A regex string (alternation of all n-mutation variants), or s0 itself
        when it is too short to permit any mutation.
    """
    s0 = s0.lower()
    if len(s0) <= 6:
        return(s0)
    if n > len(s0)/6:
        n = math.floor(len(s0)/6)
    s1 = []
    combs = list(itertools.combinations(range(len(s0)), n))
    for comb in combs:
        stmp = s0
        for k in range(len(comb)):
            stmp = add_mut(stmp, comb[k])
        for j in range(len(stmp)-1, -1, -1):
            x = stmp[j]
            if x.upper() == x:
                x = f'[^{stmp[j].lower()}]'
            stmp = stmp[:j] + x + stmp[j+1:]
        s1.append(stmp)
    return('|'.join(s1))


def core_nbrs(core, corelist, txt, title="", wnd=75, rev=False):
    """Compute a Jaccard-style proximity score between core and every other
    sequence in corelist.

    For each candidate in corelist a window of ±wnd nucleotides around every
    occurrence of core and of the candidate is inspected.  The score for a
    candidate is the mean over all windows of min(cnt_core, cnt_cand) /
    max(cnt_core, cnt_cand), where the counts are the number of occurrences
    of each sequence inside the window.

    Args:
        core:     Query core sequence.
        corelist: List of sequences to compare against core.
        txt:      Source text.
        title:    Unused (reserved for plot title if this is ever plotted).
        wnd:      Half-window size in nucleotides (default 75).
        rev:      When True, sort the result in ascending score order
                  (default False = descending).

    Returns:
        Dict mapping each corelist entry to its mean proximity score, filtered
        to entries with score > 0 and sorted by score.
    """
    corelist = list(set(corelist).difference(set(core)))
    pos_1 = find_all_matches(core, txt)
    cntr = pos_1
    L = len(core)
    for i in range(len(cntr)):
        cntr[i] += math.floor(L/2)
    Dist = dict()
    for cn in range(len(corelist)):
        allcntr = cntr.copy()
        pos_2 = find_all_matches(corelist[cn], txt)
        L2 = len(corelist[cn])
        for j in range(len(pos_2)):
            if pos_2[j] + math.floor(L2/2) in allcntr:
                continue
            allcntr.append(pos_2[j] + math.floor(L2/2))
        cnt_1 = [0]*len(allcntr)
        cnt_2 = [0]*len(allcntr)
        J = [0]*len(allcntr)
        for k in range(len(allcntr)):
            # for each window, count how many copies of the two core are in it
            for i in range(len(pos_1)):
                if abs(pos_1[i] - allcntr[k]) <= wnd + math.floor(L/2):
                    cnt_1[k] += 1
            for i in range(len(pos_2)):
                if abs(pos_2[i] - allcntr[k]) <= wnd + math.floor(L2/2):
                    cnt_2[k] += 1
            denom = max(cnt_1[k], cnt_2[k])
            J[k] = min(cnt_1[k], cnt_2[k]) / denom if denom > 0 else 0.0
        Dist[corelist[cn]] = statistics.mean(J)
    return({k: v for k, v in sorted(Dist.items(), key=lambda item: item[1], reverse=rev) if v > 0})


def cond_prob_core(core, words, txt, wnd=75, rev=False):
    """Estimate the conditional probability P(word | core) for each word.

    For each occurrence of core in txt, the function counts how many times
    each word in words appears within a window of ±wnd nucleotides.  The
    conditional probability is defined as
    count_co-occurrences / count_core, capped at 1.0.

    Args:
        core:  Query core sequence.
        words: List of candidate sequences whose conditional probability is
               to be estimated.
        txt:   Source text.
        wnd:   Half-window size in nucleotides (default 75).
        rev:   When True, sort the result in ascending probability order
               (default False = descending).

    Returns:
        Dict mapping each word to its estimated conditional probability,
        filtered to entries with probability > 0 and sorted by probability.
    """
    pos_1 = find_all_matches(core, txt)
    cnt = [0]*len(words)
    pr = dict()
    for i in range(len(pos_1)):
        L = max(0, pos_1[i] - wnd)
        R = min(len(txt)+1, pos_1[i] + len(core) + wnd)
        for cn in range(len(words)):
            pos_2 = find_all_matches(words[cn], txt)
            for k in range(len(pos_2)):
                if L <= pos_2[k] + len(words[cn])/2 <= R:
                    cnt[cn] += 1
    for cn in range(len(words)):
        if len(pos_1) == 0:
            pr[words[cn]] = 0
        else:
            pr[words[cn]] = min(1.0, cnt[cn]/len(pos_1))
    return({k: v for k, v in sorted(pr.items(), key=lambda item: item[1], reverse=rev) if v > 0})


def extendRNA(rna, xtype=0):
    """Extend the RNA to include also the complementary, reverse, and
    reverse/complementary versions of a transcript.
    0: no extension
    1: add reverse
    2: add complementary
    3: add reverse complementary
    4: all versions
    """
    if xtype == 0:
        return(rna)
    trantab = rna.maketrans("ATGCatgc", "TACGtacg")
    rnarev = rna[::-1]
    rnacom = rna.translate(trantab)
    rnarevcom = rnacom[::-1]
    if xtype  == 1:
        return(rnarev)
    if xtype == 2:
        return(rnacom)
    if xtype == 3:
        return(rnarevcom)
    return(rna+rnacom+rnarev+rnarevcom)


class Hairpin(NamedTuple):
    """One hairpin structure found in an RNA/DNA sequence.

    Coordinates are 0-based indices into the original sequence:
      - stem 1 (5' arm) : rna[start : start+stem_len]
      - loop             : rna[start+stem_len : start+stem_len+loop_len]
      - stem 2 (3' arm) : rna[start+stem_len+loop_len : end]
    """
    start:    int   # first position of stem 1
    end:      int   # one-past-last position of stem 2
    stem_len: int   # number of paired bases in each arm
    loop_len: int   # number of unpaired bases in the loop
    stem_seq: str   # stem 1 sequence (5' arm)
    loop_seq: str   # loop sequence


def gen_hairpins(rna, minSL=8, minLL=3, maxLL=40) -> List[Hairpin]:
    """Find all hairpin structures in an RNA/DNA sequence.

    Args:
        rna:   nucleotide sequence (string, lowercase or uppercase)
        minSL: minimum stem length in base pairs (default 8)
        minLL: minimum loop length in nucleotides  (default 3)
        maxLL: maximum loop length in nucleotides  (default 40)

    Returns:
        List of Hairpin named tuples sorted by start position.
        Each hairpin covers rna[hp.start:hp.end] and carries the
        stem and loop sequences for direct downstream use.
        Returns an empty list when no hairpins are found.
    """
    Lrna = len(rna)
    revrna = extendRNA(rna, xtype=3)   # reverse complement

    # --- seed dictionaries: position -> minSL-mer ---
    seqdict = {f'{i:0>6}': rna[i:i+minSL]
               for i in range(Lrna - 2*minSL - minLL)}
    revdict  = {f'{i:0>6}': revrna[i:i+minSL]
                for i in range(Lrna - 2*minSL - minLL)
                if revrna[i:i+minSL] in seqdict.values()}
    if not revdict:
        return []

    seqdict = {f'{i:0>6}': rna[i:i+minSL]
               for i in range(Lrna - 2*minSL - minLL)
               if rna[i:i+minSL] in revdict.values()}
    if not seqdict:
        return []

    seq_to_keys = defaultdict(list)
    for k, v in seqdict.items():
        seq_to_keys[v].append(k)
    rev_to_keys = defaultdict(list)
    for k, v in revdict.items():
        rev_to_keys[v].append(k)

    hairpins: List[Hairpin] = []
    occupied = [False] * Lrna   # prevents reporting overlapping hairpins twice

    for seed, rev_positions in rev_to_keys.items():
        for pr in rev_positions:
            for px in seq_to_keys[seed]:
                leng = Lrna - int(pr) - int(px)
                if leng < 2*minSL + minLL:
                    continue
                if rna[int(px):int(px)+leng] == revrna[int(pr):int(pr)+leng]:
                    continue    # perfect self-reverse-complement — skip

                # extend the stem as far as bases match
                stemlen = 0
                for ps in range(leng):
                    if rna[int(px)+ps] == revrna[int(pr)+ps]:
                        stemlen += 1
                    else:
                        break

                looplen = leng - 2*stemlen
                if stemlen > minSL and minLL <= looplen <= maxLL:
                    start = int(px)
                    end   = start + leng
                    if not any(occupied[start:end]):
                        hairpins.append(Hairpin(
                            start    = start,
                            end      = end,
                            stem_len = stemlen,
                            loop_len = looplen,
                            stem_seq = rna[start : start + stemlen],
                            loop_seq = rna[start + stemlen : start + stemlen + looplen],
                        ))
                        occupied[start:end] = [True] * leng

    hairpins.sort(key=lambda hp: hp.start)
    return hairpins


def confmat(txtb, txt2):
    """Compute a confusion matrix comparing word boundaries in two texts.

    txtb is the ground-truth text where spaces mark true word boundaries; txt2
    is the text produced by the boundary-finding algorithm.  The two texts
    contain the same non-space characters in the same order; only the space
    positions differ.

    The confusion matrix cells are:
    * **tp** — boundaries present in both txtb and txt2 at the same position.
    * **fp** — boundaries in txt2 not present in txtb (false alarms).
    * **fn** — boundaries in txtb missing from txt2 (missed boundaries).
    * **tn** — non-boundary positions correctly left unmarked (derived as
               total_non_space - tp - fp - fn, returned as the fourth value).

    Args:
        txtb: Ground-truth text with spaces as boundary markers.
        txt2: Algorithm output text with spaces as boundary markers.

    Returns:
        Tuple (tp, fp, fn, tn).
    """
    tp = 0
    fp = 0
    fn = 0
    tot = len(txtb.replace(' ', ''))
    txtb = txtb.rstrip(' ').lstrip(' ')
    txt2 = txt2.rstrip(' ').lstrip(' ')
    while max(len(txtb), len(txt2)) > 0:
        idx1 = txtb.find(" ")
        idx2 = txt2.find(" ")
        if idx1 == idx2:
            if idx1 < 0:
                break
            tp += 1
            txtb = txtb[idx1+1:]
            txt2 = txt2[idx2+1:]
            continue
        if idx1 < idx2: # the boundary finding algorithm missed one
            if idx1 >= 0:
                fn += 1
                txt2 = txt2[idx1:]
            else:
                txt2 = ''
                fp += len(find_all_matches(' ', txtb))
            txtb = txtb[idx1+1:]
        else: # the boundary finding algorithm misidentified one
            if idx2 >= 0:
                fp += 1
                txtb = txtb[idx2:]
            else:
                txtb = ''
                fn += len(find_all_matches(' ', txt2))
            txt2 = txt2[idx2+1:]
    return(tp, fp, fn, tot-tp-fp-fn)


def zscore(x, robust=True):
    """Standardise an array of values to z-scores.

    Args:
        x:      Numeric array-like of values.
        robust: When True (default) use median and IQR (Q75 - Q25) as the
                location and scale estimates, giving a robust z-score
                resistant to outliers.  When False, use mean and standard
                deviation.

    Returns:
        Array of z-scores with the same shape as x.
    """
    if robust:
        m = statistics.median(x)
        s = np.percentile(x, 75) - np.percentile(x, 25)
    else:
        m = statistics.mean(x)
        s = np.std(x)
    if s == 0:
        return np.zeros_like(np.asarray(x), dtype=float)
    z = (x - m)/s
    return(z)


def find_with_mutations(seq, txt, mutr=1/6, M=4):
    """Find exact and approximate (Hamming) occurrences of seq in txt.

    Returns:
        exact_pos  - list of start positions where seq matches exactly
        approx     - list of (pos, matched_str, distance) for 0 < d <= maxmut
        maxmut     - maximum mutations allowed (0 when seq is too short)
    """
    seq = seq.lower()
    txt_lower = txt.lower()
    L = len(seq)
    exact_pos = find_all_matches(seq, txt_lower, ret='pos')
    maxmut = 0
    approx = []
    if L >= 1 / mutr:
        maxmut = min(M, math.floor(L * mutr))
        n_txt = len(txt_lower)
        for j in range(n_txt - L + 1):
            d = 0
            for k in range(L):
                if txt_lower[j + k] != seq[k]:
                    d += 1
                    if d > maxmut:
                        break
            if 0 < d <= maxmut:
                approx.append((j, txt_lower[j:j + L], d))
    return exact_pos, approx, maxmut


def extend_match_pair(s, txt, pos1, pos2, mutr=1/6):
    """Extend two exact occurrences of s (at pos1 and pos2 in txt) left and right,
    keeping the total Hamming distance ≤ floor(total_length * mutr).

    The right direction is extended first, consuming part of the budget; the left
    direction uses the remaining budget.  Within each direction the extension
    continues past individual mismatches as long as there is still hope of ending
    within budget (X-drop style stopping condition).

    Args:
        s     - seed string (must occur exactly at pos1 and pos2)
        txt   - full source text
        pos1  - 0-based start of first occurrence
        pos2  - 0-based start of second occurrence
        mutr  - maximum mismatch rate (e.g. 1/6 means 1 mismatch per 6 nt)

    Returns:
        (left_ext, right_ext, ext1, ext2, total_ham)
        left_ext, right_ext  — characters added on each side of the seed
        ext1, ext2           — the full extended substrings from txt
        total_ham            — Hamming distance of ext1 vs ext2
    """
    L = len(s)
    n = len(txt)

    # ── extend right ────────────────────────────────────────────────────────
    max_r = min(n - pos1 - L, n - pos2 - L)
    # absolute ceiling on mismatches (used as early-abort guard)
    max_budget_r = math.floor((L + max_r) * mutr)
    best_r, best_r_ham = 0, 0
    cur_ham = 0
    for k in range(max_r):
        cur_ham += (txt[pos1 + L + k] != txt[pos2 + L + k])
        if cur_ham > max_budget_r:
            break                             # no recovery possible
        ext = k + 1
        if cur_ham <= math.floor((L + ext) * mutr):
            best_r = ext
            best_r_ham = cur_ham

    # ── extend left (budget shared with right) ──────────────────────────────
    max_l = min(pos1, pos2)
    max_budget_l = math.floor((L + best_r + max_l) * mutr) - best_r_ham
    best_l, best_l_ham = 0, 0
    cur_ham = 0
    for k in range(max_l):
        cur_ham += (txt[pos1 - 1 - k] != txt[pos2 - 1 - k])
        if cur_ham > max_budget_l:
            break
        ext = k + 1
        if cur_ham <= math.floor((L + best_r + ext) * mutr) - best_r_ham:
            best_l = ext
            best_l_ham = cur_ham

    total_ham = best_r_ham + best_l_ham
    ext1 = txt[pos1 - best_l: pos1 + L + best_r]
    ext2 = txt[pos2 - best_l: pos2 + L + best_r]
    return best_l, best_r, ext1, ext2, total_ham


def find_longest_extensions(s, txt, mutr=1/6, outfile=''):
    """Find all pairs of exact occurrences of s in txt and compute the longest
    Hamming-bounded extension for each pair (see extend_match_pair).

    Args:
        s       - seed string
        txt     - source text
        mutr    - mismatch rate threshold
        outfile - if non-empty, write results to this CSV path; if the string
                  is exactly True, auto-generate a name from s (default: '')

    Returns:
        List of dicts sorted by total_len descending.  Each dict contains:
            pos1, pos2   — 0-based start positions of the two seed copies
            left_ext     — nt added to the left
            right_ext    — nt added to the right
            total_len    — left_ext + len(s) + right_ext
            ext1, ext2   — the extended substrings
            hamming      — Hamming distance of ext1 vs ext2
    """
    import csv as _csv
    import re as _re

    positions = find_all_matches(s, txt, ret='pos')
    L = len(s)
    results = []
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            p1, p2 = positions[i], positions[j]
            l, r, e1, e2, h = extend_match_pair(s, txt, p1, p2, mutr=mutr)
            results.append({
                'pos1': p1, 'pos2': p2,
                'left_ext': l, 'right_ext': r,
                'total_len': L + l + r,
                'ext1': e1, 'ext2': e2,
                'hamming': h,
            })
    results.sort(key=lambda x: x['total_len'], reverse=True)

    if outfile:
        if outfile is True:
            slug = _re.sub(r'[^a-zA-Z0-9_-]', '_', s)[:40]
            outfile = f'extensions_{slug}.csv'
        with open(outfile, 'w', newline='') as _f:
            w = _csv.writer(_f)
            w.writerow(['rank', 'pos1', 'pos2', 'left_ext', 'right_ext',
                        'total_len', 'hamming', 'ext1', 'ext2'])
            for rank, res in enumerate(results, 1):
                w.writerow([rank, res['pos1'], res['pos2'],
                            res['left_ext'], res['right_ext'],
                            res['total_len'], res['hamming'],
                            res['ext1'], res['ext2']])

    return results


# ── k-mer scramble / null-distribution analysis ──────────────────────────────

def scramble_kmer_pvalues(txt: str, k: int, N: int, seed: int = 0) -> list:
    """Compute empirical p-values for every k-mer observed in txt.

    Shuffles txt N times (preserving nucleotide composition) and counts each
    k-mer in every shuffled copy.  For each k-mer the p-value is the fraction
    of shuffles in which that k-mer's count was >= its real count in txt.

    Args:
        txt  - source sequence (lower-cased internally)
        k    - k-mer length
        N    - number of shuffles
        seed - random seed for reproducibility

    Returns:
        List of dicts sorted by p-value ascending (most over-represented
        first).  Each dict has: kmer, real_count, exceed_count, pvalue.
    """
    import random as _random
    _random.seed(seed)
    txt = txt.lower()
    chars = list(txt)
    n = len(txt)

    real_kmers = Counter(txt[i: i + k] for i in range(n - k + 1))

    # exceed[kmer] = number of shuffles where kmer count >= real count
    exceed = Counter()
    for _ in range(N):
        _random.shuffle(chars)
        s = ''.join(chars)
        shuffled = Counter(s[i: i + k] for i in range(n - k + 1))
        for kmer, real_count in real_kmers.items():
            if shuffled.get(kmer, 0) >= real_count:
                exceed[kmer] += 1

    results = [
        {'kmer': kmer, 'real_count': count,
         'exceed_count': exceed[kmer], 'pvalue': exceed[kmer] / N}
        for kmer, count in real_kmers.items()
    ]
    results.sort(key=lambda x: x['pvalue'])
    return results
