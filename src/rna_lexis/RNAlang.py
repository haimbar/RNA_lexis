"""Backward-compatibility shim for rna_lexis.RNAlang.

All symbols that previously lived in this module have been reorganised into
focused sub-modules.  This file re-exports everything so that existing code
using ``from rna_lexis.RNAlang import ...`` or ``import rna_lexis.RNAlang
as RNA`` continues to work unchanged.

New code should import directly from the appropriate sub-module:
  * rna_lexis.algorithms  — pure computation (no I/O, no plotting)
  * rna_lexis.alignment   — Gotoh global/local alignment
  * rna_lexis.plots       — matplotlib / Plotly visualisation
  * rna_lexis.io          — file, session, and network I/O
  * rna_lexis.dialogs     — tkinter file/directory chooser wrappers
"""

import re  # keep accessible as RNAlang.re for legacy callers

from rna_lexis.algorithms import (
    contains_only_rna,
    count_kgrams,
    find_all_matches,
    stitch,
    find_max_cover,
    split_words,
    cover,
    cores,
    print_core,
    is_bounded,
    find_boundary,
    expand_to_boundary,
    _is_bounded_fast,
    _expand_left,
    _expand_right,
    add_mut,
    allow_mutation,
    core_nbrs,
    cond_prob_core,
    extendRNA,
    Hairpin,
    gen_hairpins,
    confmat,
    zscore,
    find_with_mutations,
    extend_match_pair,
    find_longest_extensions,
)

from rna_lexis.alignment import (
    NEG_INF,
    AlignmentResult,
    score_sub,
    make_markers,
    print_alignment,
    gotoh_global,
    gotoh_local,
)

from rna_lexis.plots import (
    plot_seq_nbrs,
    plot_nbrs_condensed,
    export_nbrs_condensed,
    plot_logo,
    plotzscore,
    plotkmerhist,
    plot_frequency_rank,
    plot_coverage,
    _HIT_PALETTE,
    _HIT_PALETTE_RGBA,
    _blend_over_white,
    _compress_positions,
    plot_sequence_hits,
    plot_sequence_hits_detailed,
)

from rna_lexis.io import (
    _SESSION_REQUIRED_KEYS,
    SessionData,
    save_session,
    load_session,
    is_valid_session,
    open_file_with_default_software,
    _close_in_excel,
    read_text,
    open_pdf,
    init_summary,
)
