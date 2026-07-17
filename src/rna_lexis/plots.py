"""Visualisation functions (matplotlib and Plotly)."""

import csv
import os
import sys
import math
import shutil
import subprocess
import numpy as np
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeoutError
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.patches import PathPatch
import plotly.graph_objs as go
from collections import Counter, defaultdict

from rna_lexis.algorithms import (
    cond_prob_core, core_nbrs, allow_mutation, find_all_matches,
    zscore, contains_only_rna, hop_distance,
)
from rna_lexis.io import open_file_with_default_software, open_pdf


def plot_seq_nbrs(s0, s, txt, sortby='CP', wd=20, title='', file='',
                  xrange=[], scale=1, hairpins=[], min_occ=2):
    """Plot the genomic positions of a core sequence together with its
    neighbouring cores and single-nucleotide mutants.

    Produces an interactive Plotly figure with three horizontal bands:

    * **Row 0** – occurrences of the query s0 (light-blue rectangles; label
      in dark blue).
    * **Rows 1..n** – related cores from s that co-occur near s0.  Labels and
      rectangle borders are dark purple for non-overlapping neighbours and
      forest green for cores that overlap or contain s0.  Neighbours with
      fewer than ``min_occ`` co-occurrences (globally, across the full
      sequence) are omitted.
    * **Rows -1..-m** – single-nucleotide mutants of s0 (crimson labels;
      always shown regardless of ``min_occ``).

    A salmon-coloured window of width ``wd`` flanking each s0 occurrence is
    highlighted.  If hairpin regions are supplied, orange bands are drawn
    below row 0.  The horizontal legend uses a 13 pt font and 260 px entry
    width to prevent column overlap.

    Args:
        s0:       Query core sequence.
        s:        List of all cores to search for neighbours in.
        txt:      Full source sequence.
        sortby:   Metric used to rank neighbouring cores — ``'CP'`` or
                  ``'cp'`` for conditional probability (default); ``'IOU'``
                  or ``'iou'`` for Jaccard/intersection-over-union.
        wd:       Half-window size (in nt) around each s0 hit used to
                  identify co-occurring neighbours (default 20).
        title:    Plot title string.
        file:     Output path for saving (.png/.svg/.html).  When empty the
                  figure is only displayed (default '').
        xrange:   Optional [start, end] limits for the x-axis.  If empty or
                  invalid, the full sequence range is shown.
        scale:    Pixel scale for raster export (default 1; use 3 for
                  high-resolution PNG).
        hairpins: List of hairpin dicts (keys: start, end, stem_seq,
                  loop_seq) to draw as orange bands beneath row 0.
        min_occ:  Minimum number of global occurrences in the full sequence
                  required to include a neighbour.  Mutations are always kept
                  regardless of this threshold (default 2).
    """
    L = len(txt)
    lc = len(s0)
    maxlenstr = 0
    txtl = txt.lower()  # compute once; reused everywhere below

    if sortby.upper() == 'CP':
        s1 = cond_prob_core(s0, s, txt, wnd=wd, rev=True)
    else:
        s1 = core_nbrs(s0, s, txt, wnd=wd, rev=True)
    if s0 in s1:
        del s1[s0]
    related = list(s1.keys())

    if min_occ > 1:
        related = [seq for seq in related
                   if len(find_all_matches(seq, txtl, ret='pos')) >= min_occ]
    n = len(related)

    # discover mutations, then cache their positions once for both sorting and rendering
    mutset = {s0}
    if lc >= 6:
        for _ in range(math.floor(lc / 6)):
            prev = set(mutset)
            for mj in prev:
                mutset |= set(Counter(
                    find_all_matches(allow_mutation(mj), txtl, ret='str')
                ).keys())
            if mutset == prev:
                break
    arrmut_raw = list(mutset - {s0})
    mut_positions = {m: find_all_matches(m, txtl, ret='pos') for m in arrmut_raw}
    arrmut = sorted(arrmut_raw, key=lambda m: len(mut_positions[m]), reverse=True)
    nm = len(arrmut)

    plotxrange = [0, L]
    if len(xrange) == 2 and 0 <= xrange[0] < xrange[1] <= L:
        plotxrange = xrange

    arr0 = [p for p in find_all_matches(s0, txtl, ret='pos')
            if plotxrange[0] <= p < plotxrange[1]]

    # cache neighbour positions filtered to the visible range
    nbr_positions = {seq: [p for p in find_all_matches(seq, txtl, ret='pos')
                            if plotxrange[0] <= p < plotxrange[1]]
                     for seq in related}

    layout = go.Layout(
        showlegend=True, hovermode='closest',
        hoverlabel=dict(bgcolor='white', font_size=16, font_family='Rockwell'),
        xaxis=dict(range=[plotxrange[0], plotxrange[1] * 1.07], autorange=False),
        yaxis_range=[-nm - 0.5, n + 2], title=title,
    )
    fig = go.Figure(layout=layout)

    shapes = []
    traces = []
    annotations = []

    # s0 occurrences
    for p0 in arr0:
        shapes.append(dict(type='rect', x0=p0, x1=p0 + lc, y0=0, y1=0.9,
                           fillcolor='rgba(173,216,230,0.5)',
                           line=dict(color='blue', width=1)))
        shapes.append(dict(type='rect', xref='x', yref='paper',
                           x0=p0 - wd - lc, x1=p0 + wd + lc, y0=0, y1=1,
                           fillcolor='LightSalmon', opacity=0.3,
                           layer='below', line=dict(width=0)))
        traces.append(go.Scatter(
            x=[p0 + lc / 2], y=[0.5],
            hovertemplate=f'{s0} {p0}', showlegend=False, mode='markers',
            marker=dict(size=6, color='darkblue', symbol='x')))

    # hairpins
    for hp in hairpins:
        shapes.append(dict(type='rect', x0=hp['start'], x1=hp['end'],
                           y0=-0.15, y1=0,
                           fillcolor='rgba(255,165,0,0.6)',
                           line=dict(color='darkorange', width=1)))
        mid = (hp['start'] + hp['end']) / 2
        hover_hp = (f"pos: {hp['start']}–{hp['end']}<br>"
                    f"stem: {hp['stem_seq']}<br>loop: {hp['loop_seq']}")
        traces.append(go.Scatter(
            x=[mid], y=[-0.075], mode='markers',
            marker=dict(size=max(6, (hp['end'] - hp['start']) // 4),
                        color='darkorange', opacity=0),
            customdata=[hover_hp],
            hovertemplate='%{customdata}<extra>Hairpin</extra>',
            showlegend=False))

    annotations.append(dict(x=plotxrange[1] * 1.001, y=0.35, text=s0,
                             showarrow=False, xanchor='left',
                             font=dict(family='courier', size=16, color='darkblue')))

    # related cores
    for i, seq in enumerate(related):
        is_cont = (seq in s0) or (s0 in seq)
        col  = 'rgb(0,140,50)'      if is_cont else 'rgb(140,45,140)'
        colt = 'rgba(0,255,0,0.3)'  if is_cont else 'rgba(218,112,214,0.3)'
        arr = nbr_positions[seq]
        for p in arr:
            mrk, sz, colm = 'square', 6, 'blue'
            if not is_cont:
                for p0 in arr0:
                    if (p0 < p) and (p < p0 + lc < p + len(seq)):
                        traces.append(go.Scatter(
                            x=[p + len(seq) / 2], y=[i + 1.5],
                            hovertemplate=f'{seq} {p}',
                            showlegend=False, mode='markers',
                            marker=dict(size=8, color='green', symbol='triangle-right')))
                        mrk, sz, colm = 'triangle-right', 12, 'green'
                    if (p < p0) and (p0 < p + len(seq) < p0 + lc):
                        traces.append(go.Scatter(
                            x=[p + len(seq) / 2], y=[i + 1.5],
                            hovertemplate=f'{seq} {p}',
                            showlegend=False, mode='markers',
                            marker=dict(size=6, color='green', symbol='triangle-left')))
                        mrk, sz, colm = 'triangle-left', 12, 'green'
            shapes.append(dict(type='rect', x0=p, x1=p + len(seq),
                               y0=i + 1, y1=i + 1.9,
                               fillcolor=colt, line=dict(color=col, width=1)))
            traces.append(go.Scatter(
                x=[p + len(seq) / 2], y=[i + 1.5],
                hovertemplate=f'{seq} {p}',
                showlegend=False, mode='markers',
                marker=dict(size=sz, color=colm, symbol=mrk)))
        if len(seq) > maxlenstr:
            maxlenstr = len(seq)
        annotations.append(dict(x=plotxrange[1] * 1.001, y=i + 1.35, text=seq,
                                 showarrow=False, xanchor='left',
                                 font=dict(family='courier', size=16, color=col)))

    # mutations — one shape+trace per occurrence, display string computed once per mutation
    for j, mut in enumerate(arrmut):
        display = ''.join(
            c.upper() if (idx < lc and c != s0[idx]) else c
            for idx, c in enumerate(mut)
        )
        arr = [p for p in mut_positions[mut] if plotxrange[0] <= p < plotxrange[1]]
        for p in arr:
            mrk, sz, colm = 'square', 6, 'blue'
            for p0 in arr0:
                if (p0 < p) and (p < p0 + lc < p + len(mut)):
                    mrk, sz, colm = 'triangle-right', 12, 'darkgreen'
                if (p < p0) and (p0 < p + len(mut) < p0 + lc):
                    mrk, sz, colm = 'triangle-left', 12, 'darkgreen'
            shapes.append(dict(type='rect', x0=p, x1=p + len(mut),
                               y0=-1 - j, y1=-j - 0.1,
                               fillcolor='rgba(255,255,255,0.1)',
                               line=dict(color='red', width=1)))
            traces.append(go.Scatter(
                x=[p + len(mut) / 2], y=[-0.5 - j],
                hovertemplate=f'{display} {p}',
                showlegend=False, mode='markers',
                marker=dict(size=sz, color=colm, symbol=mrk)))
        if len(mut) > maxlenstr:
            maxlenstr = len(mut)
        annotations.append(dict(x=plotxrange[1] * 1.001, y=-0.55 - j, text=display,
                                 showarrow=False, xanchor='left',
                                 font=dict(family='courier', size=16, color='crimson')))

    def _leg(name, color, symbol, size=8):
        return go.Scatter(x=[None], y=[None], mode='markers', name=name,
                          showlegend=True,
                          marker=dict(size=size, color=color, symbol=symbol))

    traces += [
        _leg(s0,                                       'darkblue',         'x',              8),
        _leg('Neighboring core',                       'rgb(218,112,214)', 'square',         8),
        _leg('Complete overlap with s0',               'rgb(0,200,0)',     'square',         8),
    ]
    if nm > 0:
        traces.append(_leg('Mutation', 'red', 'square', 8))
    traces += [
        _leg('Partial overlap: core ends after s0',    'green', 'triangle-right', 10),
        _leg('Partial overlap: core starts before s0', 'green', 'triangle-left',  10),
    ]
    if hairpins:
        traces.append(go.Scatter(x=[None], y=[None], mode='lines',
                                  name='Hairpin region', showlegend=True,
                                  line=dict(color='darkorange', width=6)))

    avg_char_factor = 0.58
    right_padding_px = 20
    left_margin_px = 60
    plot_area_width_px = 700
    estimated_text_width_px = int(maxlenstr * 16 * avg_char_factor) + right_padding_px
    fig_width = left_margin_px + plot_area_width_px + estimated_text_width_px

    tick_vals = [0.45] + [i + 1.45 for i in range(n)] + [-0.55 - j for j in range(nm)]
    tick_text = ['0'] + [str(i + 1) for i in range(n)] + [str(j + 1) for j in range(nm)]

    fig.update_layout(
        shapes=shapes,
        annotations=annotations,
        width=fig_width,
        margin=dict(l=left_margin_px, r=estimated_text_width_px, t=40, b=80),
        legend=dict(orientation='h', yanchor='top', y=-0.05, xanchor='left', x=0,
                    font=dict(size=13), entrywidth=260),
        yaxis=dict(tickvals=tick_vals, ticktext=tick_text,
                   zeroline=True, zerolinewidth=2, zerolinecolor='lightgray'),
    )
    fig.add_traces(traces)
    fig.show()

    if file != '':
        _save_fig(fig, file, scale)


def _squish_lanes(intervals, gap=1):
    """Greedy interval-scheduling lane assignment.

    Args:
        intervals: iterable of (start, end, *rest) tuples, sorted by start.
        gap:       minimum nucleotide space required between consecutive
                   intervals assigned to the same lane.

    Returns:
        List of (start, end, *rest, lane) tuples (lane is 1-based).
    """
    lane_ends = []   # end position of the last interval placed in each lane
    result = []
    for item in intervals:
        start, end = item[0], item[1]
        placed = False
        for k, last_end in enumerate(lane_ends):
            if start >= last_end + gap:
                lane_ends[k] = end
                result.append((*item, k + 1))
                placed = True
                break
        if not placed:
            lane_ends.append(end)
            result.append((*item, len(lane_ends)))
    return result


_EXPORT_TIMEOUT = 45  # seconds before assuming Kaleido is hung


def _save_fig(fig, file, scale):
    """Write a Plotly figure to *file* (html or raster) and open it.

    PNG/SVG export runs with a timeout so a hung Kaleido process never freezes
    the terminal.  On failure an interactive HTML fallback is saved instead and
    a plain-text error report is written next to it.
    """
    try:
        low = file.lower()
        if low.endswith('.html') or low.endswith('.htm'):
            fig.write_html(file, auto_open=False)
            open_file_with_default_software(file)
            return
        root, ext = os.path.splitext(file)
        if not ext:
            file = f'{file}.png'
            root = file[:-4]
        with ThreadPoolExecutor(max_workers=1) as pool:
            fut = pool.submit(fig.write_image, file, scale=scale)
            try:
                fut.result(timeout=_EXPORT_TIMEOUT)
            except FuturesTimeoutError:
                msg = (
                    f'Plotly static image export timed out after {_EXPORT_TIMEOUT} s.\n'
                    f'This usually means Kaleido is not working correctly.\n'
                    f'Fix: pip install --upgrade "kaleido>=0.2,<1"'
                )
                print(msg)
                html_file = f'{root}_fallback.html'
                fig.write_html(html_file, auto_open=False)
                print(f'Saved interactive HTML fallback to: {html_file}')
                with open(f'{root}_export_message.txt', 'w', encoding='utf-8') as fh:
                    fh.write(msg + f'\nHTML fallback saved to: {html_file}\n')
                open_file_with_default_software(html_file)
                return
        open_file_with_default_software(file)
    except Exception as e:
        print(f'Could not save/open plot: {file}\n{e}')


def plot_nbrs_condensed(s0, s, txt, sortby='CP', wd=20, title='', file='',
                        xrange=[], scale=1, hairpins=[], min_occ=2):
    """Condensed neighbour plot: three bands — s0, all-neighbour density strip, mutations.

    Shows the query s0 at row 0, all neighbouring sequences collapsed into a
    single semi-transparent density strip at row 1 (green = overlapping with
    s0, orchid = non-overlapping; intensity reflects how many sequences share
    each position), and single-nucleotide mutants in squished lanes below row 0.
    No wd windows and no sequence labels are drawn inside the plot.

    Args:
        s0:       Query core sequence.
        s:        List of all cores/xmotifs to search for neighbours in.
        txt:      Full source sequence.
        sortby:   Ranking metric — ``'CP'`` (default) or ``'IOU'``.
        wd:       Half-window size in nt (default 20).
        title:    Plot title (optional).
        file:     Output path; ``_condensed`` is inserted before the extension.
        xrange:   Optional [start, end] x-axis limits.
        scale:    Pixel scale for raster export.
        hairpins: List of hairpin dicts (keys: start, end, stem_seq, loop_seq).
        min_occ:  Minimum number of global occurrences in the full sequence
                  required to include a neighbour (default 2).
    """
    L = len(txt)
    lc = len(s0)
    txtl = txt.lower()  # compute once; reused everywhere below
    plotxrange = [0, L]
    if len(xrange) == 2 and 0 <= xrange[0] < xrange[1] <= L:
        plotxrange = xrange

    if sortby.upper() == 'CP':
        s1 = cond_prob_core(s0, s, txt, wnd=wd, rev=True)
    else:
        s1 = core_nbrs(s0, s, txt, wnd=wd, rev=True)
    if s0 in s1:
        del s1[s0]
    related = list(s1.keys())
    if min_occ > 1:
        related = [seq for seq in related
                   if len(find_all_matches(seq, txtl, ret='pos')) >= min_occ]
    n_related = len(related)
    density_row = 1

    # discover mutations, then cache positions once for sorting and squish-building
    mutset = {s0}
    if lc >= 6:
        for _ in range(math.floor(lc / 6)):
            prev = set(mutset)
            for mj in prev:
                mutset |= set(Counter(
                    find_all_matches(allow_mutation(mj), txtl, ret='str')
                ).keys())
            if mutset == prev:
                break
    arrmut_raw = list(mutset - {s0})
    mut_positions = {m: find_all_matches(m, txtl, ret='pos') for m in arrmut_raw}
    arrmut = sorted(arrmut_raw, key=lambda m: len(mut_positions[m]), reverse=True)

    arr0 = [p for p in find_all_matches(s0, txtl, ret='pos')
            if plotxrange[0] <= p < plotxrange[1]]

    mutation_row = -1

    # cache related sequence positions filtered to the visible range
    related_positions = {seq: [p for p in find_all_matches(seq, txtl, ret='pos')
                                if plotxrange[0] <= p < plotxrange[1]]
                         for seq in related}

    layout = go.Layout(
        showlegend=True, hovermode='closest',
        hoverlabel=dict(bgcolor='white', font_size=16, font_family='Rockwell'),
        xaxis=dict(range=[plotxrange[0], plotxrange[1] * 1.07], autorange=False),
        yaxis_range=[(-1.3 if arrmut else -0.3), density_row + 1.1],
        title=(f'{title} (condensed)') if title else 'Core neighbours (condensed)',
    )
    fig = go.Figure(layout=layout)

    shapes = []
    traces = []

    # s0 occurrences
    for p0 in arr0:
        shapes.append(dict(type='rect', x0=p0, x1=p0 + lc, y0=0, y1=0.9,
                           fillcolor='rgba(173,216,230,0.5)',
                           line=dict(color='blue', width=1)))
        traces.append(go.Scatter(
            x=[p0 + lc / 2], y=[0.5],
            hovertemplate=f'{s0} {p0}', showlegend=False, mode='markers',
            marker=dict(size=6, color='darkblue', symbol='x')))

    # hairpins
    for hp in hairpins:
        shapes.append(dict(type='rect', x0=hp['start'], x1=hp['end'],
                           y0=-0.15, y1=0,
                           fillcolor='rgba(255,165,0,0.6)',
                           line=dict(color='darkorange', width=1)))
        mid = (hp['start'] + hp['end']) / 2
        hover_hp = (f"pos: {hp['start']}–{hp['end']}<br>"
                    f"stem: {hp['stem_seq']}<br>loop: {hp['loop_seq']}")
        traces.append(go.Scatter(
            x=[mid], y=[-0.075], mode='markers',
            marker=dict(size=max(6, (hp['end'] - hp['start']) // 4),
                        color='darkorange', opacity=0),
            customdata=[hover_hp],
            hovertemplate='%{customdata}<extra>Hairpin</extra>',
            showlegend=False))

    # density strip — all related sequences at row 1, using cached positions
    if n_related > 0:
        # alpha: floor at 0.15 so sparse positions stay visible; cap at 0.5 so
        # dense regions don't saturate; n_related fully-overlapping → ~90% opacity
        alpha = max(0.15, min(0.5, 1 - 0.1 ** (1 / n_related)))
        for seq in related:
            is_cont = (seq in s0) or (s0 in seq)
            a_fill   = f'{alpha:.3f}'
            a_border = f'{min(1.0, alpha * 2):.3f}'
            col_fill   = (f'rgba(0,200,0,{a_fill})'   if is_cont
                          else f'rgba(218,112,214,{a_fill})')
            col_border = (f'rgba(0,150,0,{a_border})' if is_cont
                          else f'rgba(180,80,180,{a_border})')
            occ = related_positions[seq]
            for p in occ:
                shapes.append(dict(type='rect', x0=p, x1=p + len(seq),
                                   y0=density_row, y1=density_row + 0.9,
                                   fillcolor=col_fill,
                                   line=dict(color=col_border, width=0.5)))
            if occ:
                traces.append(go.Scatter(
                    x=[p + len(seq) / 2 for p in occ],
                    y=[density_row + 0.5] * len(occ),
                    hovertemplate=f'{seq} %{{x:.0f}}<extra></extra>',
                    showlegend=False, mode='markers',
                    marker=dict(size=6, color='rgba(0,0,0,0)')))

    # dotted separator between s0 row and density strip
    shapes.append(dict(type='line', xref='paper', yref='y',
                       x0=0, x1=1, y0=density_row, y1=density_row,
                       line=dict(width=1, dash='dot', color='gray')))

    # mutations as a density strip — always a single band at mutation_row
    if arrmut:
        alpha_mut = max(0.15, min(0.5, 1 - 0.1 ** (1 / len(arrmut))))
        a_fill_mut   = f'{alpha_mut:.3f}'
        a_border_mut = f'{min(1.0, alpha_mut * 2):.3f}'
        for mut in arrmut:
            display = ''.join(
                c.upper() if (i < lc and c != s0[i]) else c
                for i, c in enumerate(mut)
            )
            occ = [p for p in mut_positions[mut] if plotxrange[0] <= p < plotxrange[1]]
            for p in occ:
                shapes.append(dict(type='rect', x0=p, x1=p + len(mut),
                                   y0=mutation_row, y1=mutation_row + 0.9,
                                   fillcolor=f'rgba(255,80,80,{a_fill_mut})',
                                   line=dict(color=f'rgba(200,0,0,{a_border_mut})', width=0.5)))
            if occ:
                traces.append(go.Scatter(
                    x=[p + len(mut) / 2 for p in occ],
                    y=[mutation_row + 0.5] * len(occ),
                    hovertemplate=f'{display} %{{x:.0f}}<extra></extra>',
                    showlegend=False, mode='markers',
                    marker=dict(size=6, color='rgba(0,0,0,0)')))

    tick_vals = ([0.45, density_row + 0.45]
                 + ([mutation_row + 0.45] if arrmut else []))
    tick_text = (['s0', f'neighbors ({n_related})']
                 + (['mutations'] if arrmut else []))

    def _leg(name, color, symbol, size=8):
        return go.Scatter(x=[None], y=[None], mode='markers', name=name,
                          showlegend=True,
                          marker=dict(size=size, color=color, symbol=symbol))

    traces += [
        _leg(s0,                                    'darkblue',         'x',      8),
        _leg('Neighbors (non-overlapping with s0)', 'rgb(140,45,140)',  'square', 8),
        _leg('Neighbors (overlapping with s0)',     'rgb(0,140,50)',    'square', 8),
    ]
    if arrmut:
        traces.append(_leg('Mutation', 'crimson', 'square', 8))
    if hairpins:
        traces.append(go.Scatter(x=[None], y=[None], mode='lines',
                                  name='Hairpin region', showlegend=True,
                                  line=dict(color='darkorange', width=6)))

    fig.update_layout(
        shapes=shapes,
        width=800,
        margin=dict(l=60, r=40, t=30, b=45),
        legend=dict(orientation='h', yanchor='top', y=-0.05, xanchor='left', x=0,
                    font=dict(size=13), entrywidth=260),
        yaxis=dict(tickvals=tick_vals, ticktext=tick_text,
                   zeroline=True, zerolinewidth=2, zerolinecolor='lightgray'),
    )
    fig.add_traces(traces)
    fig.show()

    if file:
        root, ext = os.path.splitext(file)
        _save_fig(fig, f'{root}_condensed{ext or ".png"}', scale)


def export_nbrs_condensed(s0, s, txt, sortby='CP', wd=20, xrange=[], file=''):
    """Export context windows around s0 occurrences and its neighbours as CSV.

    Uses the same logic as the condensed plot: finds neighbours via CP/IOU,
    discovers single-nucleotide mutations, then pools every occurrence of s0,
    its mutations, and its neighbours.  Each occurrence is expanded by ±wd;
    overlapping expanded intervals are merged into contiguous regions.
    Columns: start, end, seq — where seq = txt[start:end+1].

    Args:
        s0:     Query core sequence.
        s:      List of all cores/xmotifs to search for neighbours in.
        txt:    Full source sequence.
        sortby: Ranking metric — 'CP' (default) or 'IOU'.
        wd:     Half-window size in nt (default 20).
        xrange: Optional [start, end] positional limits.
        file:   Output CSV path. Always prints to stdout; also saves if given.

    Returns:
        Path of the written CSV, or None if no file path was given.
    """
    L = len(txt)
    lc = len(s0)
    txtl = txt.lower()
    plotxrange = [0, L]
    if len(xrange) == 2 and 0 <= xrange[0] < xrange[1] <= L:
        plotxrange = xrange

    # neighbours — identical computation to the condensed plot
    if sortby.upper() == 'CP':
        s1 = cond_prob_core(s0, s, txt, wnd=wd, rev=True)
    else:
        s1 = core_nbrs(s0, s, txt, wnd=wd, rev=True)
    if s0 in s1:
        del s1[s0]
    related = list(s1.keys())

    # mutations
    mutset = {s0}
    if lc >= 6:
        for _ in range(math.floor(lc / 6)):
            prev = set(mutset)
            for mj in prev:
                mutset |= set(Counter(
                    find_all_matches(allow_mutation(mj), txtl, ret='str')
                ).keys())
            if mutset == prev:
                break
    arrmut = list(mutset - {s0})

    # build ±wd expanded intervals for every occurrence of s0, mutations, and neighbours
    intervals = []
    for seq in [s0] + arrmut + related:
        L_seq = len(seq)
        for p in find_all_matches(seq, txtl, ret='pos'):
            if plotxrange[0] <= p < plotxrange[1]:
                intervals.append((max(0, p - wd), min(L - 1, p + L_seq + wd - 1)))

    intervals.sort()

    # merge overlapping intervals; end = end of the last sequence in each cluster
    merged = []
    for ws, we in intervals:
        if merged and ws <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], we)
        else:
            merged.append([ws, we])

    header = ['start', 'end', 'seq']
    rows = [[ws, we, txt[ws:we + 1]] for ws, we in merged]

    print(','.join(header))
    for row in rows:
        print(','.join(str(x) for x in row))

    if not file:
        return None

    root, ext = os.path.splitext(file)
    out = file if ext.lower() == '.csv' else f'{root}.csv'
    with open(out, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)
    print(f"Saved to '{out}'")
    open_file_with_default_software(out)
    return out


def plot_logo(wd, k, txt, outfile, fmt='pdf', muts=0):
    """Generate a sequence logo for a core and its flanking k nucleotides.

    Finds all occurrences of wd in txt (including up to muts mismatches if
    muts > 0), writes them as a FASTA file, and calls the WebLogo command-line
    tool to produce a logo plot.  The output is opened automatically.

    Args:
        wd:      Core/motif sequence to build the logo for.
        k:       Number of flanking nucleotides included on each side.
        txt:     Source RNA/DNA sequence (non-RNA sequences are rejected).
        outfile: Base output path (without extension); the logo is saved as
                 ``<outfile>.<fmt>``.
        fmt:     Output format — ``'pdf'`` (default) or ``'svg'`` (requires
                 pdf2svg in PATH).
        muts:    Number of allowed substitutions per occurrence (default 0).
    """
    if not contains_only_rna(txt):
        print('Logo plots are only available for RNA data')
        return(None)
    wd0 = wd
    pad = "".join(['.']*k)
    if muts > 0:
        wd = wd + '|' + allow_mutation(wd, n=muts)
    if k > 0:
        wd = wd.replace('|', f'{pad}|{pad}')
        wd = pad + wd + pad

    found = find_all_matches(f'{wd}', txt, ret='str')
    if len(found) == 0:
        print("No matching sequences found for logo plot. Try a different core/motif.")
        return(None)
    faseq = '\n'.join(found)
    weblogoexe = shutil.which("weblogo") or shutil.which("weblogo.exe")
    if weblogoexe is None and sys.platform == "win32":
        # In a venv on Windows, sys.executable is usually ...\venv\Scripts\python.exe
        py_dir = os.path.dirname(sys.executable)
        candidates = [
            os.path.join(py_dir, "weblogo.exe"),
            os.path.join(py_dir, "weblogo"),
            os.path.join(py_dir, "weblogo-script.py"),
            os.path.join(py_dir, "Scripts", "weblogo.exe"),
            os.path.join(py_dir, "Scripts", "weblogo"),
            os.path.join(py_dir, "Scripts", "weblogo-script.py"),
            os.path.join(os.getcwd(), "venv", "Scripts", "weblogo.exe"),
            os.path.join(os.getcwd(), ".venv", "Scripts", "weblogo.exe"),
        ]
        for candidate in candidates:
            if os.path.isfile(candidate):
                weblogoexe = candidate
                break
    if weblogoexe is None:
        print("Can't find the 'weblogo' executable, so logo plot was skipped.")
        print("Install inside your venv with:")
        print(r"  .\venv\Scripts\python.exe -m pip install weblogo")
        print("Then run again.")
        return(None)
    weblogoexe = weblogoexe.strip()
    with open(f"{outfile}.fa", "w", encoding="utf-8") as f:
        f.write(faseq)
    ttl = f'{wd0}, {len(found)} sequences'
    if fmt.lower() == "svg" and shutil.which("pdf2svg") is None:
        print("SVG output requires 'pdf2svg', which was not found in PATH.")
        print("Falling back to PDF output.")
        fmt = "pdf"
    cmd = [weblogoexe, "-U",  "probability", "-c", "classic", "-t", ttl, "-F", fmt, "-f", f"{outfile}.fa", "-o", f"{outfile}.{fmt}"]
    if weblogoexe.lower().endswith(".py"):
        cmd = [sys.executable] + cmd
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print("WebLogo failed to generate the plot.")
        if e.stderr:
            print(e.stderr.strip())
        elif e.stdout:
            print(e.stdout.strip())
        return(None)

    out_plot = f'{outfile}.{fmt}'
    if not os.path.isfile(out_plot) or os.path.getsize(out_plot) == 0:
        print(f"Logo output was not created: {out_plot}")
        return(None)
    open_pdf(out_plot)


def plotzscore(cnt, zmin=1.96, robust=True, file='', scale=1):
    """Plot normalised k-mer frequencies as z-scores with outlier labels.

    Each k-mer's frequency is normalised by the total count and then
    z-scored.  Outliers beyond ±zmin are labelled with their k-mer string and
    connected to their data points with leader lines.

    Args:
        cnt:    Dict mapping k-mer string to raw count (e.g. from count_kgrams).
        zmin:   Z-score threshold for labelling outliers (default 1.96 ≈ 5 %
                two-tailed).
        robust: Passed to zscore(); use robust (median/IQR) standardisation
                when True (default).
        file:   Output path for saving the figure.  When empty the figure is
                only displayed.
        scale:  DPI multiplier for raster export (default 1).
    """
    cnt1 = {k: v for k, v in sorted(cnt.items(),
                                    key=lambda item: item[1], reverse=True)}
    x0 = np.array(list(cnt1.values()))
    x0 = x0/np.sum(x0)
    z0 = zscore(x0, robust)

    _, ax = plt.subplots()
    ax.scatter(range(len(z0)), z0, s=0.75, c='red')
    ax.hlines(0, xmin=0, xmax=len(z0))
    ax.hlines([zmin, -zmin], xmin=0, xmax=len(z0), linestyles='dashed')

    keys = list(cnt1.keys())
    n = len(z0)
    y_range = float(np.max(z0) - np.min(z0)) or 1.0
    x_step = n * 0.05        # primary: shift each lower-ranked label further right
    y_nudge = y_range * 0.02 # secondary: tiny vertical nudge only when y values are identical

    # Positive outliers: sort descending — highest z (rank 0) stays at actual position
    pos = sorted([(i, float(z0[i])) for i in range(n) if z0[i] > zmin], key=lambda t: -t[1])
    # Negative outliers: sort ascending — most negative (rank 0) stays at actual position
    neg = sorted([(i, float(z0[i])) for i in range(n) if z0[i] < -zmin], key=lambda t: t[1])

    max_tx = float(n)
    for outliers, direction in ((pos, 1), (neg, -1)):
        last_ty = None
        for rank, (xi, yi) in enumerate(outliers):
            tx = xi + rank * x_step
            ty = yi
            if last_ty is not None and abs(ty - last_ty) < y_nudge:
                ty = last_ty + direction * y_nudge
            arrow = dict(arrowstyle='-', color='gray', lw=0.5) if rank > 0 else None
            ax.annotate(keys[xi], xy=(xi, yi), xytext=(tx, ty),
                        fontsize=8, arrowprops=arrow)
            print(keys[xi], yi)
            last_ty = ty
            max_tx = max(max_tx, tx)

    ax.set_xlim(left=-0.5, right=max_tx + n * 0.05)
    if file:
        try:
            plt.savefig(file, bbox_inches='tight', dpi=150 * scale)
            open_file_with_default_software(file)
        except Exception as e:
            print(f'Could not save plot: {file}\n{e}')
    plt.show()


def plotkmerhist(cnt, k, file='', scale=1):
    """Plot a histogram of k-mer abundance (frequency-of-frequency).

    The x-axis shows each distinct k-mer count value; the y-axis shows
    the number of distinct k-mers that occur exactly that many times.

    Args:
        cnt:   Dict mapping k-mer string to its count (e.g. from count_kgrams).
        k:     k-mer length (used only in the plot title).
        file:  Output path for saving the figure.  When empty the figure is
               only displayed.
        scale: DPI multiplier for raster export (default 1).
    """
    vals = Counter(cnt.values())
    plt.bar(vals.keys(), vals.values(), color='skyblue')
    plt.xlabel('k‑mer count')
    plt.ylabel('number of distinct k‑mers')
    plt.title(f'K-mer abundance (k={k})')
    if file:
        try:
            plt.savefig(file, bbox_inches='tight', dpi=150 * scale)
            open_file_with_default_software(file)
        except Exception as e:
            print(f'Could not save plot: {file}\n{e}')
    plt.show()


def plot_frequency_rank(freq_dict, k, file='', scale=1):
    """Plot frequency vs. rank on a log-log scale (Zipf's law diagnostic).

    k-mers are ranked from most to least frequent, and the resulting
    frequency-rank curve is plotted on logarithmic axes.  A power-law
    distribution (Zipf) appears as a straight line on this plot.

    Args:
        freq_dict: Dict mapping k-mer string to count.
        k:         k-mer length (used only in the plot title).
        file:      Output path for saving the figure.  When empty the figure
                   is only displayed.
        scale:     DPI multiplier for raster export (default 1).
    """
    sorted_items = sorted(freq_dict.items(), key=lambda item: item[1], reverse=True)
    frequencies = [item[1] for item in sorted_items]
    ranks = range(1, len(frequencies) + 1)
    plt.figure(figsize=(10, 6))
    plt.plot(ranks, frequencies, marker='o', linestyle='-')
    plt.xlabel('Rank')
    plt.ylabel('Frequency')
    plt.title(f"Frequency vs. Rank Plot (Zipf's Law), k={k}")
    plt.grid(True, which='both', linestyle='--')
    plt.xscale('log')
    plt.yscale('log')
    if file:
        try:
            plt.savefig(file, bbox_inches='tight', dpi=150 * scale)
            open_file_with_default_software(file)
        except Exception as e:
            print(f'Could not save plot: {file}\n{e}')
    plt.show()


def plot_coverage(txt, cvrs, pwr, file='', scale=1):
    """Plot weighted coverage scores for a collection of sequences as a bar chart.

    Sequences are sorted by ascending coverage so the highest-scoring
    sequence is the rightmost bar.

    Args:
        txt:   Source text (used only for its length in the plot title).
        cvrs:  Dict mapping sequence string to its coverage score
               (e.g. as returned by cover()).
        pwr:   Length exponent used when computing cvrs (shown in the title).
        file:  Output path for saving the figure.  When empty the figure is
               only displayed.
        scale: DPI multiplier for raster export (default 1).
    """
    cvrs = {k: v for k, v in sorted(cvrs.items(), key=lambda item: item[1])}
    plt.figure(figsize=(8, 5))
    plt.bar(cvrs.keys(), cvrs.values(), color='skyblue')
    plt.xlabel('Sequence')
    plt.xticks(rotation='vertical')
    plt.ylabel('Weighted coverage')
    plt.title(f'Text length={len(txt)}, weight={pwr}')
    if file:
        try:
            plt.savefig(file, bbox_inches='tight', dpi=150 * scale)
            open_file_with_default_software(file)
        except Exception as e:
            print(f'Could not save plot: {file}\n{e}')
    plt.show()


# alpha=0.55 lets overlapping sequences from different queries show through.
# All boxes use the same fill; approximate matches are distinguished by red
# mutation lines drawn inside the rectangle instead of a different opacity.
_HIT_PALETTE = [
    ('#1f77b4', 'rgba(31,119,180,0.55)'),  # blue
    ('#ff7f0e', 'rgba(255,127,14,0.55)'),  # orange
    ('#2ca02c', 'rgba(44,160,44,0.55)'),   # green
]

# Raw RGBA components for the detailed HTML view (alpha compositing over white)
_HIT_PALETTE_RGBA = [
    (31,  119, 180, 0.45),   # blue
    (255, 127,  14, 0.45),   # orange
    ( 44, 160,  44, 0.45),   # green
]


def _blend_over_white(indices):
    """Alpha-composite palette colors (by index) over white; return 'rgb(r,g,b)'."""
    r, g, b = 255.0, 255.0, 255.0
    for i in indices:
        cr, cg, cb, ca = _HIT_PALETTE_RGBA[i % len(_HIT_PALETTE_RGBA)]
        r = cr * ca + r * (1 - ca)
        g = cg * ca + g * (1 - ca)
        b = cb * ca + b * (1 - ca)
    return f'rgb({int(round(r))},{int(round(g))},{int(round(b))})'


def _compress_positions(all_intervals, txt_len, padding=50, stub_width=80,
                        min_gap=1000):
    """Compute a coordinate mapping that squashes large gaps between hit clusters.

    Intervals within `padding` of each other are merged into one visible segment.
    Only gaps whose actual size exceeds `min_gap` are compressed; smaller gaps
    are preserved at their true width. Each compressed gap is replaced by a
    fixed `stub_width` on the axis.

    Returns:
        transform(x)        maps original position to compressed coordinate
        tick_vals, tick_text original positions shown as x-axis ticks
        gap_xs              compressed-x positions for vertical separator lines
        total               total width of compressed axis
    """
    if not all_intervals:
        raise ValueError("_compress_positions: all_intervals must be non-empty")
    padded = sorted((max(0, s - padding), min(txt_len, e + padding))
                    for s, e in all_intervals)
    merged = [list(padded[0])]
    for s, e in padded[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])

    # Build (orig_s, orig_e, comp_s, comp_e) for each visible segment;
    # only compress gaps that exceed min_gap.
    segments, comp, gap_xs = [], 0, []
    for k, (s, e) in enumerate(merged):
        if k > 0:
            gap_size = s - merged[k - 1][1]
            if gap_size > min_gap:
                gap_xs.append(comp + stub_width / 2)
                comp += stub_width
            else:
                comp += gap_size   # keep the gap at its true width
        segments.append((s, e, comp, comp + (e - s)))
        comp += (e - s)
    total = comp

    def transform(x):
        """Map an original sequence position x to its compressed axis coordinate."""
        for s, e, cs, ce in segments:
            if s <= x <= e:
                return cs + (x - s)
        # x falls between two segments
        for k in range(len(segments) - 1):
            s_next = segments[k + 1][0]
            e_prev = segments[k][1]
            if e_prev < x < s_next:
                gap_size = s_next - e_prev
                if gap_size > min_gap:
                    # compressed gap — return midpoint of stub
                    return segments[k][3] + stub_width / 2
                else:
                    # uncompressed gap — map linearly
                    return segments[k][3] + (x - e_prev)
        return segments[0][2] if x < segments[0][0] else segments[-1][3]

    # Ticks at segment boundaries (deduplicated)
    tick_vals, tick_text, seen = [], [], set()
    for s, e, cs, ce in segments:
        for orig, cv in [(s, cs), (e, ce)]:
            if cv not in seen:
                tick_vals.append(cv)
                tick_text.append(str(orig))
                seen.add(cv)

    return transform, tick_vals, tick_text, gap_xs, total


def plot_sequence_hits(seq_data, txt_len, title='', file='', scale=1,
                       condense=False, min_gap=1000, hairpins=[]):
    """Plot occurrences of up to 3 sequences as stacked colored lanes.

    Each sequence occupies its own horizontal lane (seq 1 at bottom, seq 3 at
    top). Exact matches are solid rectangles; approximate matches carry red
    vertical lines at each mismatch position.

    Args:
        seq_data  - list of dicts: {seq, L, exact: [pos, ...],
                                    approx: [(pos, str, dist), ...]}
        txt_len   - full length of the source text (sets x-axis range)
        title     - plot title
        file      - optional output path (.png/.svg/.html); '' = display only
        scale     - pixel scale for raster export (3 = high-res)
        condense  - if True, compress gaps larger than min_gap on the x-axis
        min_gap   - minimum gap size (in source positions) to compress
    """
    # ── coordinate transform ──────────────────────────────────────────────
    if condense:
        all_ivs = [(pos, pos + d['L'])
                   for d in seq_data for pos in d['exact']] + \
                  [(pos, pos + d['L'])
                   for d in seq_data for pos, _, _ in d['approx']]
        if all_ivs:
            xform, tick_vals, tick_text, gap_xs, total = \
                _compress_positions(all_ivs, txt_len, min_gap=min_gap)
            x_axis = dict(title='Position', range=[0, total],
                          tickvals=tick_vals, ticktext=tick_text)
        else:
            condense = False  # nothing to compress

    if not condense:
        xform = lambda x: x
        x_axis = dict(title='Position', range=[0, txt_len])
        gap_xs = []

    # ── figure ────────────────────────────────────────────────────────────
    n_seqs = len(seq_data)
    y_min = -0.3 if hairpins else 0
    height = 200 + 130 * n_seqs + (60 if hairpins else 0)
    fig = go.Figure(layout=go.Layout(
        title=title or 'Sequence occurrences',
        hovermode='closest',
        showlegend=True,
        xaxis=x_axis,
        yaxis=dict(visible=False, range=[y_min, n_seqs]),
        width=1000,
        height=height,
        margin=dict(l=60, r=40, t=60, b=160),
        legend=dict(orientation='h', yanchor='top', y=-0.75, xanchor='left', x=0),
    ))

    # vertical dashed lines at compressed-gap positions
    for gx in gap_xs:
        fig.add_shape(type='line', x0=gx, x1=gx, y0=0, y1=n_seqs,
                      line=dict(color='gray', width=1, dash='dash'))

    for i, d in enumerate(seq_data):
        line_col, fill = _HIT_PALETTE[i % len(_HIT_PALETTE)]
        seq = d['seq']
        L = d['L']
        lane_y0 = i + 0.05
        lane_y1 = i + 0.95
        lane_mid = i + 0.5

        # Exact matches
        for pos in d['exact']:
            x0c = xform(pos)
            x1c = xform(pos + L)
            fig.add_shape(type='rect', x0=x0c, x1=x1c, y0=lane_y0, y1=lane_y1,
                          fillcolor=fill, line=dict(color=line_col, width=1))

        # Approximate matches — same fill as exact; red lines mark each mismatch
        for pos, matched, _dist in d['approx']:
            x0c = xform(pos)
            x1c = xform(pos + L)
            span = x1c - x0c
            fig.add_shape(type='rect', x0=x0c, x1=x1c, y0=lane_y0, y1=lane_y1,
                          fillcolor=fill, line=dict(color=line_col, width=1))
            for k, (sc, mc) in enumerate(zip(seq, matched)):
                if sc != mc:
                    xm = x0c + span * (k + 0.5) / L
                    fig.add_shape(type='line', x0=xm, x1=xm,
                                  y0=lane_mid, y1=lane_y1,
                                  line=dict(color='red', width=2))

        # Invisible scatter for hover — exact
        if d['exact']:
            fig.add_trace(go.Scatter(
                x=[(xform(pos) + xform(pos + L)) / 2 for pos in d['exact']],
                y=[lane_mid] * len(d['exact']),
                mode='markers',
                marker=dict(size=max(6, L // 2), color=line_col, opacity=0),
                customdata=[f'{seq}<br>pos={pos} (exact)' for pos in d['exact']],
                hovertemplate='%{customdata}<extra></extra>',
                showlegend=False,
            ))

        # Invisible scatter for hover — approx
        if d['approx']:
            fig.add_trace(go.Scatter(
                x=[(xform(pos) + xform(pos + L)) / 2 for pos, _, _ in d['approx']],
                y=[lane_mid] * len(d['approx']),
                mode='markers',
                marker=dict(size=max(6, L // 2), color=line_col, opacity=0),
                customdata=[
                    ''.join(
                        f'<span style="color:red"><b>{mc.upper()}</b></span>'
                        if mc != sc else mc
                        for sc, mc in zip(seq, matched)
                    ) + f'<br>pos={pos}  dist={dist}'
                    for pos, matched, dist in d['approx']
                ],
                hovertemplate='%{customdata}<extra></extra>',
                showlegend=False,
            ))

        # Legend entry (colored square + counts)
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(size=14, color=line_col, symbol='square'),
            name=f'{seq}   exact={len(d["exact"])}  ≈{len(d["approx"])}',
        ))

    # Fixed legend entries explaining visual markers
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode='lines',
        line=dict(color='red', width=2),
        name='red line — mutation position',
    ))
    if gap_xs:
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='lines',
            line=dict(color='gray', width=1, dash='dash'),
            name='dashed line — compressed gap',
        ))

    for hp in hairpins:
        x0c = xform(hp['start'])
        x1c = xform(hp['end'])
        fig.add_shape(type='rect', x0=x0c, x1=x1c, y0=-0.25, y1=-0.02,
                      fillcolor='rgba(255,165,0,0.6)',
                      line=dict(color='darkorange', width=1))
        mid = (x0c + x1c) / 2
        hover = (f"pos: {hp['start']}–{hp['end']}<br>"
                 f"stem: {hp['stem_seq']}<br>"
                 f"loop: {hp['loop_seq']}")
        fig.add_trace(go.Scatter(
            x=[mid], y=[-0.135],
            mode='markers',
            marker=dict(size=max(6, (hp['end'] - hp['start']) // 4),
                        color='darkorange', opacity=0),
            customdata=[hover],
            hovertemplate='%{customdata}<extra>Hairpin</extra>',
            showlegend=False,
        ))
    if hairpins:
        fig.add_trace(go.Scatter(
            x=[None], y=[None], mode='lines',
            line=dict(color='darkorange', width=6),
            name='Hairpin region',
        ))

    fig.show()
    if file:
        _save_fig(fig, file, scale)


def plot_sequence_hits_detailed(seq_data, txt, start, end, title='', file=''):
    """Nucleotide-level view of sequence hits within a specified range.

    Each character in txt[start:end] is shown as monospace text. Hit regions
    are highlighted with colored backgrounds (same palette as
    plot_sequence_hits). Within approximate hits, mutation positions are shown
    in red uppercase; all other characters are black lowercase. Output is
    always an HTML file with plain styled text (no plotting library needed).

    Args:
        seq_data  - list of dicts: {seq, L, exact: [pos, ...],
                                    approx: [(pos, str, dist), ...]}
        txt       - the full source text (str)
        start/end - 0-based slice of txt to display
        title     - page title
        file      - output HTML path (required; opens after saving)
    """
    start = max(start, 0)
    end   = min(end, len(txt))
    n     = end - start
    region = txt[start:end].lower()

    # Per-character state: set of sequence indices covering it, and mutation flag
    chars  = list(region)
    bg     = [set() for _ in range(n)]   # set of seq_data indices
    is_mut = [False] * n                  # True → red uppercase

    for i, d in enumerate(seq_data):
        seq = d['seq']
        L   = d['L']

        for pos in d['exact']:
            for k in range(L):
                abs_k = pos + k - start
                if 0 <= abs_k < n:
                    bg[abs_k].add(i)

        for pos, matched, _dist in d['approx']:
            for k in range(L):
                abs_k = pos + k - start
                if 0 <= abs_k < n:
                    bg[abs_k].add(i)
                    if seq[k] != matched[k]:
                        chars[abs_k] = matched[k].upper()
                        is_mut[abs_k] = True

    # Build rows of LINE characters each, prefixed with the genomic position
    LINE = 60
    pw   = len(str(end - 1))
    rows = []
    for row_s in range(0, n, LINE):
        row_e = min(row_s + LINE, n)
        label = str(start + row_s).rjust(pw)
        spans = []
        for k in range(row_s, row_e):
            c        = chars[k]
            indices  = sorted(bg[k])
            mut_sty  = 'color:red;font-weight:bold' if is_mut[k] else ''

            if not indices:
                spans.append(f'<span style="{mut_sty}">{c}</span>' if mut_sty else c)
            else:
                bg_col = _blend_over_white(indices)
                sty = f'background:{bg_col}'
                if mut_sty:
                    sty += f';{mut_sty}'
                spans.append(f'<span style="{sty}">{c}</span>')
        rows.append(
            f'<span style="color:#888;user-select:none">{label} </span>'
            + ''.join(spans)
        )

    # Legend
    legend_items = []
    for i, d in enumerate(seq_data):
        line_col, fill = _HIT_PALETTE[i % len(_HIT_PALETTE)]
        swatch = (f'<span style="display:inline-block;width:1em;height:1em;'
                  f'background:{fill};border:1px solid {line_col};'
                  f'vertical-align:middle;margin-right:4px"></span>')
        legend_items.append(
            f'<li>{swatch}<code>{d["seq"]}</code> &nbsp; '
            f'exact: {len(d["exact"])}, &nbsp; approx: {len(d["approx"])}</li>'
        )
    legend_items.append(
        '<li><span style="color:red;font-weight:bold">RED UPPERCASE</span>'
        ' &mdash; mutation position</li>'
    )
    legend_html = '<ul style="font-size:0.9em;margin-top:1em">' + \
                  ''.join(legend_items) + '</ul>'

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title or 'Sequence detail'}</title>
<style>
  body  {{ font-family: sans-serif; padding: 1.5em; }}
  h2    {{ margin-bottom: 0.3em; }}
  pre   {{ font-family: 'Courier New', Courier, monospace;
           font-size: 13px; line-height: 1.6;
           white-space: pre-wrap; word-break: break-all;
           border: 1px solid #ddd; padding: 1em;
           background: #fafafa; border-radius: 4px; }}
</style>
</head>
<body>
<h2>{title or 'Sequence detail'}</h2>
<p style="color:#555;font-size:0.9em">
  Positions {start}&ndash;{end-1} &nbsp;|&nbsp; {n} nucleotides
</p>
{legend_html}
<pre>{'<br>'.join(rows)}</pre>
</body>
</html>
"""

    if not file:
        file = 'sequence_detail.html'
    try:
        with open(file, 'w', encoding='utf-8') as fh:
            fh.write(html)
        open_file_with_default_software(file)
    except Exception as e:
        print(f'Could not save detail view: {file}\n{e}')


def _semicircle(x0, x1, n=80):
    """Return (x, y) arrays tracing a semicircle spanning x0 to x1."""
    cx, r = (x0 + x1) / 2, (x1 - x0) / 2
    t = np.linspace(0, np.pi, n)
    return cx + r * np.cos(t), r * np.sin(t)


_HOP_ARC_COLORS = {1: "#2166ac", 2: "#74add1", 3: "#fdae61", 4: "#d73027"}
_HOP_ARC_DEFAULT_COLOR = "#888888"  # hop >= 5 (beyond the 4-bucket legend)


def plot_self_similarity_arcs(seed, results, positions, unit=None,
                               spacing_stats=None, arcs='consecutive',
                               label_threshold=15, file='', scale=1):
    """Semicircular arc diagram of pairwise Hamming-bounded extensions.

    Visualises the output of find_longest_extensions() for a single seed
    motif within one sequence: one semicircular arc per pair of exact seed
    occurrences, thickness proportional to extension length, color encoding
    hop distance (physical spacing bucketed by the estimated tandem-repeat
    unit length), and node size proportional to the longest extension
    reachable from that position. Intended for self-similarity analysis of
    tandem-repeat regions.

    Args:
        seed:      Seed motif string, shown in the title.
        results:   Output of find_longest_extensions(seed, txt, ...) — list
                   of dicts with keys pos1, pos2, total_len, hamming, for
                   *every* pair of occurrences (n occurrences -> n*(n-1)/2
                   pairs). Node sizes are always computed from the full set
                   (see arcs, below, for which pairs are actually drawn).
        positions: Sorted list of unique occurrence positions (the arc
                   diagram's x-axis nodes) — typically the pos1/pos2 values
                   appearing in results.
        unit:      Estimated repeat-unit length in nt, used to bucket hop
                   distance (see algorithms.hop_distance()). If None
                   (default), estimated as the median gap between
                   consecutive sorted positions — no fixed value is assumed,
                   since repeat-unit length is specific to whichever
                   motif/region is being analysed.
        spacing_stats: Optional result dict from
                   statistical.spacing_periodicity_test(positions, len(txt))
                   — the same test used by "Motif spacing / periodicity
                   test". Requires >= 3 occurrences to be meaningful (2
                   occurrences give a single gap, not a distribution to
                   test); when given, its period and p-values are shown on
                   the plot as an additional reference line. When None, that
                   line is omitted — the arc diagram itself doesn't require
                   regular spacing to be drawn or interpreted.
        arcs:      Which pairs to draw arcs for. 'consecutive' (default)
                   draws only adjacent-in-position pairs (n-1 arcs for n
                   occurrences, scaling linearly) -- for a seed with many
                   occurrences the full pairwise set (n*(n-1)/2 arcs) is
                   illegible, and the periodicity signal this plot is about
                   is fundamentally a consecutive-spacing story anyway.
                   'all' draws every pair in `results`, matching the
                   original (pre-decluttering) behavior.
        label_threshold: When the number of *drawn* arcs exceeds this, the
                   per-arc spacing/Hamming/identity text labels are omitted
                   (the arcs themselves still draw) -- labels overlap into
                   unreadable noise well before the arcs do. Default 15.
        file:      Output path for saving the figure (PNG/SVG/PDF). When
                   empty the figure is only displayed.
        scale:     DPI multiplier for raster export (default 1).
    """
    if not results or len(positions) < 2:
        print('Need at least 2 occurrences with computed extensions to plot.')
        return

    if unit is None:
        gaps = [b - a for a, b in zip(positions[:-1], positions[1:])]
        unit = float(np.median(gaps)) if gaps else 1.0

    if arcs == 'consecutive':
        sorted_pos = sorted(positions)
        consecutive_pairs = set(zip(sorted_pos[:-1], sorted_pos[1:]))
        drawn = [r for r in results if (r['pos1'], r['pos2']) in consecutive_pairs]
    else:
        drawn = results

    show_labels = len(drawn) <= label_threshold

    tlen_vals = [r['total_len'] for r in drawn]
    tlen_min, tlen_max = min(tlen_vals), max(tlen_vals)
    lw_scale = lambda t: 0.6 + 3.5 * (t - tlen_min) / max(tlen_max - tlen_min, 1)

    def hop_color(spacing):
        return _HOP_ARC_COLORS.get(hop_distance(spacing, unit), _HOP_ARC_DEFAULT_COLOR)

    # Node sizes reflect the longest extension reachable from that position
    # across *all* pairs, even when only a subset of arcs is drawn -- a
    # node's best match elsewhere is still meaningful context.
    max_len_map = defaultdict(int)
    for r in results:
        max_len_map[r['pos1']] = max(max_len_map[r['pos1']], r['total_len'])
        max_len_map[r['pos2']] = max(max_len_map[r['pos2']], r['total_len'])

    node_sizes = np.array([max_len_map[n] for n in positions], dtype=float)
    # np.ptp(arr), not arr.ptp() -- the ndarray method was removed in NumPy 2.0.
    node_r = 8 + 28 * (node_sizes - node_sizes.min()) / max(np.ptp(node_sizes), 1)

    fig, ax = plt.subplots(figsize=(11, 5.2))

    for r in sorted(drawn, key=lambda r: r['total_len'], reverse=True):
        p0, p1 = r['pos1'], r['pos2']
        spacing = p1 - p0
        identity = 100 * (r['total_len'] - r['hamming']) / r['total_len']
        xs, ys = _semicircle(p0, p1)
        ax.plot(xs, ys, color=hop_color(spacing), lw=lw_scale(r['total_len']),
                 alpha=0.80, zorder=2)
        if show_labels:
            mid_y = ys.max()
            fontsize = 7.5 if spacing <= unit * 1.1 else 6.5
            ax.text((p0 + p1) / 2, mid_y + 2,
                    f"{spacing} nt  h={r['hamming']}  ({identity:.1f}%)",
                    ha="center", va="bottom", fontsize=fontsize, color="#444")

    for n, nr in zip(positions, node_r):
        ax.scatter(n, 0, s=nr ** 2 * 0.3, color="#2a6099",
                   zorder=5, edgecolors="white", linewidths=1.2)
        ax.text(n, -18, str(n), ha="center", va="top", fontsize=8, color="#333")

    max_hop = max(hop_distance(r['pos2'] - r['pos1'], unit) for r in drawn)
    for hops in range(1, min(max_hop, 4) + 1):
        ax.plot([], [], color=_HOP_ARC_COLORS[hops], lw=2.5,
                label=f"{hops} hop{'s' if hops > 1 else ''} (~{round(hops * unit)} nt)")
    if max_hop > 4:
        ax.plot([], [], color=_HOP_ARC_DEFAULT_COLOR, lw=2.5, label="≥ 5 hops")
    # Node-size key: real circle markers at the same radii as the plotted
    # nodes (not line-thickness swatches — those track arc thickness, a
    # different quantity, and would misrepresent circle size here).
    node_min, node_max = node_sizes.min(), node_sizes.max()
    node_mid = (node_min + node_max) / 2
    for t in sorted({node_min, node_mid, node_max}):
        r_t = 8 + 28 * (t - node_min) / max(node_max - node_min, 1)
        ax.scatter([], [], s=r_t ** 2 * 0.3, color="#2a6099",
                   edgecolors="white", linewidths=0.8, label=f"extension len={int(t)}")
    # handleheight/labelspacing are sized for the largest possible node-size
    # marker (node_r maxes out at 36 regardless of the data — see the node_r
    # formula above), so rows never overlap however the data happens to look.
    # Placed outside the axes (right side) rather than in a corner: arcs can
    # reach high on either side depending on the data, so no inside corner
    # is reliably free of overlap for arbitrary input.
    ax.legend(title="Arc color = hop distance\nNode size = max extension length",
              fontsize=8, title_fontsize=8, framealpha=0.9,
              handleheight=3.0, labelspacing=1.2,
              loc="upper left", bbox_to_anchor=(1.01, 1.02), borderaxespad=0)

    ax.axhline(0, color="#cccccc", lw=0.8, zorder=1)
    ax.set_xlim(min(positions) - 150, max(positions) + 150)
    ax.set_ylim(-55, max((r['pos2'] - r['pos1']) / 2 for r in drawn) * 1.25)
    ax.set_xlabel("Position (nt)", fontsize=11)
    ax.set_ylabel("Arc height = spacing / 2 (nt)", fontsize=10)
    title = f"Pairwise Hamming extensions — seed: '{seed}'"
    if arcs == 'consecutive' and len(drawn) < len(results):
        title += f"  ({len(drawn)} of {len(results)} pairs shown — consecutive-neighbor only)"
    title_lines = [title]
    if spacing_stats is not None:
        def _p(v):
            return f"{v:.4f}" if v >= 0.0001 else "<0.0001"
        title_lines.append(
            f"Spacing test (n={spacing_stats['m']}): T={spacing_stats['period']} nt  ·  "
            f"cluster p={_p(spacing_stats['p_cluster'])}  ·  "
            f"Rayleigh p={_p(spacing_stats['p_rayleigh'])}"
        )
    # Split across two lines -- the combined string is long enough to
    # overflow the figure on both sides (truncated on the left, covered by
    # the legend on the right) if kept on one.
    title_lines.append("Node size = max extension length  ·  Arc thickness = extension length")
    if show_labels:
        title_lines.append("Arc color = hop distance  ·  Labels: spacing | h=Hamming | identity%")
    else:
        title_lines.append(f"Arc color = hop distance  ·  Labels hidden ({len(drawn)} arcs > {label_threshold})")
    ax.set_title("\n".join(title_lines), fontsize=10)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(labelsize=9)

    plt.tight_layout()
    if file:
        try:
            plt.savefig(file, bbox_inches='tight', dpi=150 * scale)
            open_file_with_default_software(file)
        except Exception as e:
            print(f'Could not save plot: {file}\n{e}')
    plt.show()


# ColorBrewer Set1 -- 6 perceptually distinct colors, one per shared motif.
# A dedicated palette rather than _HIT_PALETTE (only 3 entries, designed for
# exact/1mm/2mm hit categories, not enough to keep up to 6 simultaneous
# motifs visually distinct here).
_SHARED_MOTIF_COLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628"]


def _nice_tick_step(n):
    """Return a readable axis-tick spacing for a sequence of length n."""
    for step in (1000, 500, 200, 100, 50, 20, 10):
        if n / step >= 3:
            return step
    return max(1, n // 4)


def plot_shared_motif_diagram(txt_a, txt_b, label_a, label_b, shared_motifs,
                               title='', file='', scale=1):
    """Two-row arc diagram connecting shared exact-motif positions between
    two sequences (e.g. a transcript and a candidate regulatory element).

    Both sequences are drawn as horizontal bars normalised to unit width, so
    hits are readable regardless of the length difference; nucleotide
    position ticks below/above each bar give the true scale. For each motif
    in `shared_motifs`, a colored tick mark is drawn at every exact
    occurrence in each sequence, and a cubic Bezier arc connects every
    (occurrence in A, occurrence in B) pair. Arc alpha is reduced
    automatically when a motif has many occurrences, so dense connections
    stay readable instead of saturating into a solid block.

    Args:
        txt_a:         Sequence A (drawn as the top bar).
        txt_b:         Sequence B (drawn as the bottom bar).
        label_a:       Label for sequence A.
        label_b:       Label for sequence B.
        shared_motifs: Motif strings to connect (e.g. from
                       statistical.shared_exact_motifs()) -- up to the first
                       6 are drawn, each in its own color; motifs without an
                       exact hit in both sequences are skipped.
        title:         Plot title; defaults to
                       "Shared-motif connections: {label_a} <-> {label_b}".
        file:          Output path for saving the figure (PNG/SVG/PDF). When
                       empty the figure is only displayed.
        scale:         DPI multiplier for raster export (default 1).
    """
    len_a, len_b = len(txt_a), len(txt_b)

    selected = []
    for m in shared_motifs:
        pos_a = find_all_matches(m, txt_a, ret='pos')
        pos_b = find_all_matches(m, txt_b, ret='pos')
        if pos_a and pos_b:
            selected.append((m, pos_a, pos_b))
        if len(selected) == len(_SHARED_MOTIF_COLORS):
            break

    if not selected:
        print('No motifs with exact hits in both sequences -- nothing to plot.')
        return

    Y_A, Y_B = 1.10, 0.00
    BAR_H, TICK_H = 0.055, 0.055

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.set_xlim(-0.04, 1.08)
    ax.set_ylim(-0.50, 1.45)
    ax.axis('off')

    for y_ctr, label, seqlen in [(Y_A, label_a, len_a), (Y_B, label_b, len_b)]:
        ax.add_patch(plt.Rectangle((0, y_ctr - BAR_H / 2), 1, BAR_H,
                                    facecolor='#cccccc', edgecolor='#444444',
                                    lw=0.8, zorder=2))
        ax.text(-0.01, y_ctr, label, ha='right', va='center',
                fontsize=9, fontweight='bold')
        ax.text(1.01, y_ctr, f'{seqlen:,} nt', ha='left', va='center',
                fontsize=8, color='#555555')

    for y_ctr, seqlen, above in [(Y_A, len_a, True), (Y_B, len_b, False)]:
        step = _nice_tick_step(seqlen)
        for tick in range(0, seqlen + 1, step):
            x = tick / seqlen
            yb = y_ctr + BAR_H / 2 if above else y_ctr - BAR_H / 2
            ya = yb + 0.01 if above else yb - 0.01
            ylt = ya + 0.01 if above else ya - 0.01
            ax.plot([x, x], [yb, ya], color='#999999', lw=0.5, zorder=3)
            ax.text(x, ylt, str(tick), ha='center',
                    va='bottom' if above else 'top', fontsize=6, color='#666666')

    y_arc_top = Y_A - BAR_H / 2
    y_arc_bot = Y_B + BAR_H / 2
    mid_y = (y_arc_top + y_arc_bot) / 2

    legend_handles = []
    for i, (m, pos_a, pos_b) in enumerate(selected):
        color = _SHARED_MOTIF_COLORS[i]
        n_a, n_b = len(pos_a), len(pos_b)
        arc_alpha = min(0.55, max(0.10, 0.85 / math.sqrt(n_a * n_b)))

        for p in pos_a:
            x = p / len_a
            ax.plot([x, x], [Y_A + BAR_H / 2, Y_A + BAR_H / 2 + TICK_H],
                    color=color, lw=1.5, alpha=0.90, solid_capstyle='round', zorder=5)
        for p in pos_b:
            x = p / len_b
            ax.plot([x, x], [Y_B - BAR_H / 2 - TICK_H, Y_B - BAR_H / 2],
                    color=color, lw=2.5, alpha=0.90, solid_capstyle='round', zorder=5)

        for pa in pos_a:
            x0 = pa / len_a
            for pb in pos_b:
                x1 = pb / len_b
                verts = [(x0, y_arc_top), (x0, mid_y), (x1, mid_y), (x1, y_arc_bot)]
                codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
                ax.add_patch(PathPatch(Path(verts, codes), facecolor='none',
                                        edgecolor=color, alpha=arc_alpha,
                                        lw=0.9, zorder=1))

        legend_handles.append(mpatches.Patch(
            color=color, alpha=0.85,
            label=f'{m}  ({n_a}× {label_a}, {n_b}× {label_b})'))

    ax.legend(handles=legend_handles, loc='upper center',
              bbox_to_anchor=(0.5, -0.10), ncol=min(3, len(legend_handles)),
              frameon=True, fontsize=8, title='shared motif  (exact hits)',
              title_fontsize=8)

    ax.set_title(title or f'Shared-motif connections: {label_a}  ↔  {label_b}',
                 fontsize=11, pad=6)
    fig.text(0.5, 0.01,
              'Both sequences are scaled to unit width; x-axis labels show actual nucleotide positions.',
              ha='center', fontsize=7, color='#666666', style='italic')

    plt.tight_layout(rect=[0, 0.09, 1, 1])
    if file:
        try:
            plt.savefig(file, bbox_inches='tight', dpi=150 * scale)
            open_file_with_default_software(file)
        except Exception as e:
            print(f'Could not save plot: {file}\n{e}')
    plt.show()


_COVERAGE_COMPARISON_COLORS = ("#1f77b4", "#d62728")


def plot_coverage_comparison(curve_a, label_a, curve_b, label_b,
                              dual_axis=False, title='', file='', scale=1):
    """Overlay two sliding-window coverage-fraction curves for comparison.

    Both curves are expected to already be resampled to the same length
    over a common 0-100% position axis (algorithms.sliding_coverage_fraction())
    -- this keeps the function agnostic to what each curve actually measures,
    so the same plot works for motif-coverage-vs-motif-coverage (two
    sequences), or motif-coverage-vs-GC-content / homopolymer-density (one
    sequence, two different sliding-window statistics).

    Args:
        curve_a, curve_b: Equal-length float arrays, values in [0, 1].
        label_a, label_b: Legend labels for each curve.
        dual_axis:        When False (default), both curves share one 0-1
                          axis -- directly comparable in absolute terms,
                          which is the more honest view when that's what you
                          want, but visually flattens whichever curve has
                          the smaller natural range (e.g. GC content rarely
                          nears 100%, while motif coverage often does).
                          When True, curve_b gets its own y-axis
                          (matplotlib twinx), each independently scaled to
                          its own data range -- matches how the published
                          figure this was validated against (SNHG14 motif
                          coverage vs. GC%) presents it, at the cost of the
                          two curves no longer being visually comparable in
                          absolute terms.
        title:            Plot title; defaults to "{label_a} vs {label_b}".
        file:             Output path for saving the figure (PNG/SVG/PDF).
                          When empty the figure is only displayed.
        scale:            DPI multiplier for raster export (default 1).
    """
    x = np.linspace(0, 100, len(curve_a))
    color_a, color_b = _COVERAGE_COMPARISON_COLORS

    fig, ax = plt.subplots(figsize=(10, 5))

    if dual_axis:
        ax2 = ax.twinx()
        for curve, label, color, axis in ((curve_a, label_a, color_a, ax),
                                           (curve_b, label_b, color_b, ax2)):
            axis.plot(x, curve, color=color, lw=1.3, label=label)
            axis.axhline(np.mean(curve), color=color, lw=1, ls='--', alpha=0.8)
            axis.tick_params(axis='y', labelcolor=color, labelsize=9)
        ax.set_ylabel(label_a, fontsize=10, color=color_a)
        ax2.set_ylabel(label_b, fontsize=10, color=color_b)
        ax2.spines[["top"]].set_visible(False)
        lines_a, labels_a = ax.get_legend_handles_labels()
        lines_b, labels_b = ax2.get_legend_handles_labels()
        ax.legend(lines_a + lines_b, labels_a + labels_b,
                  loc='upper right', fontsize=9, framealpha=0.9)
    else:
        for curve, label, color in ((curve_a, label_a, color_a), (curve_b, label_b, color_b)):
            ax.fill_between(x, curve, alpha=0.30, color=color)
            ax.plot(x, curve, color=color, lw=1.3, label=label)
            ax.axhline(np.mean(curve), color=color, lw=1, ls='--', alpha=0.8)
        ax.set_ylim(0, 1)
        ax.set_ylabel('Coverage fraction (sliding window)', fontsize=10)
        ax.legend(loc='upper right', fontsize=9, framealpha=0.9)

    ax.set_xlim(0, 100)
    ax.set_xlabel('Position (% of sequence length)', fontsize=11)
    ax.set_title(title or f'{label_a}  vs  {label_b}\nDashed lines = mean coverage', fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=9)

    plt.tight_layout()
    if file:
        try:
            plt.savefig(file, bbox_inches='tight', dpi=150 * scale)
            open_file_with_default_software(file)
        except Exception as e:
            print(f'Could not save plot: {file}\n{e}')
    plt.show()


