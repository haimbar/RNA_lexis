"""Visualisation functions (matplotlib and Plotly)."""

import os
import sys
import math
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from collections import Counter

from rna_lexis.algorithms import (
    cond_prob_core, core_nbrs, allow_mutation, find_all_matches,
    zscore, contains_only_rna,
)
from rna_lexis.io import open_file_with_default_software, open_pdf


def plot_seq_nbrs(s0, s, txt, sortby='CP', wd=20, title='', file='',
                  xrange=[], scale=1, hairpins=[]):
    """Plot the genomic positions of a core sequence together with its
    neighbouring cores and single-nucleotide mutants.

    Produces an interactive Plotly figure with three horizontal bands:
    * **Row 0** – occurrences of the query s0 (blue rectangles).
    * **Rows 1..n** – related cores from s that co-occur near s0 (orchid/green
      rectangles; green when s0 and the neighbour overlap).
    * **Rows -1..-m** – single-nucleotide mutants of s0 (red rectangles).
    A salmon-coloured window of width wd flanking each s0 occurrence is
    highlighted.  If hairpin regions are supplied, orange bands are drawn
    below row 0.

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
    """
    # sortby should be CP, cp, IOU, or iou
    L = len(txt)
    lc = len(s0)
    maxlenstr = 0

    # find other cores (from s) within wd from s0
    if sortby.upper() == 'CP':
        s1 = cond_prob_core(s0, s, txt, wnd=wd, rev=True)
    else:
        s1 = core_nbrs(s0, s, txt, wnd=wd, rev=True)
    if s0 in list(s1.keys()):
        del s1[s0]
    n = len(s1) # number of related cores
    related = list(s1.keys())

    # search for mutations
    mutset = set({s0})
    arrmut = []
    if lc >= 6:
        for ms in range(math.floor(lc/6)):
            tmpset = mutset
            for mut_j in mutset:
                mutsrch = allow_mutation(mut_j)
                arrmuttmp = set(Counter(find_all_matches(mutsrch, txt.lower(), ret='str')).keys())
                mutset = mutset.union(arrmuttmp)
            if mutset == tmpset:
                break
    arrmut = list(mutset - {s0})
    mutcnt = dict()
    for j in range(len(arrmut)):
        arr = find_all_matches(arrmut[j], txt.lower(), ret='pos')
        mutcnt[arrmut[j]] = len(arr)
    {k: v for k, v in sorted(mutcnt.items(), key=lambda item: item[1])}
    arrmut = list({k: v for k, v in sorted(mutcnt.items(), key=lambda item: item[1], reverse=True)}.keys())

    # plot all the occurrences of s0:
    nm = len(arrmut)  # number of mutations; s0 is at y=0, neighbors above (+), mutations below (-)
    plotxrange = [0, L]
    if len(xrange) == 2:
        if 0 <= xrange[0] < xrange[1] <= L:
            plotxrange = xrange
    layout = go.Layout(showlegend=True, hovermode='closest', hoverlabel=dict(
        bgcolor="white", font_size=16, font_family="Rockwell"),
        xaxis=dict(range=[plotxrange[0], plotxrange[1]*1.07], autorange=False),
        yaxis_range=[-nm-0.5, n+2], title=title)
    fig = go.Figure(layout=layout)
    arr0 = [p for p in find_all_matches(s0, txt.lower(), ret='pos')
            if plotxrange[0] <= p < plotxrange[1]]
    for j in range(len(arr0)):
        fig.add_shape(type="rect", x0=arr0[j], x1=arr0[j]+lc, y0=0, y1=0.9,
                fillcolor='rgba(173, 216, 230, 0.5)', line=dict(color='blue', width=1))
        # highlight the regions next to s0 (width = wd)
        fig.add_vrect(
            x0=arr0[j]-wd-lc, x1=arr0[j]+wd+lc,
            fillcolor="LightSalmon", opacity=0.3,
            layer="below", line_width=0,
        )
        fig.add_trace(go.Scatter(x=[arr0[j]+lc/2], y=[0.5],
            hovertemplate = f'{s0} {arr0[j]}', showlegend=False,
                             mode='markers', marker=dict(size=6, color='darkblue', symbol='x')))
    # draw hairpin regions as orange horizontal bands alongside s0
    for hp in hairpins:
        fig.add_shape(type="rect",
            x0=hp['start'], x1=hp['end'], y0=-0.15, y1=0,
            fillcolor='rgba(255,165,0,0.6)', line=dict(color='darkorange', width=1))
        mid = (hp['start'] + hp['end']) / 2
        hover = (f"pos: {hp['start']}–{hp['end']}<br>"
                 f"stem: {hp['stem_seq']}<br>"
                 f"loop: {hp['loop_seq']}")
        fig.add_trace(go.Scatter(
            x=[mid], y=[-0.075],
            mode='markers',
            marker=dict(size=max(6, (hp['end'] - hp['start']) // 4),
                        color='darkorange', opacity=0),
            customdata=[hover],
            hovertemplate='%{customdata}<extra>Hairpin</extra>',
            showlegend=False,
        ))
    # add the legend on the right
    fig.add_annotation(x=plotxrange[1]*1.001, y=0.35, text=s0, showarrow=False,
        xanchor="left", font=dict(family="courier", size=14, color='blue'))

    # plot all the related cores:
    for i in range(len(s1)):
        col = 'rgb(218,112,214)'
        colt = 'rgba(218,112,214, 0.3)'
        if (related[i] in s0) or (s0 in related[i]):
            col = 'rgb(0,255,0)'
            colt = 'rgba(0,255,0, 0.3)'
        arr = [p for p in find_all_matches(related[i], txt.lower(), ret='pos')
               if plotxrange[0] <= p < plotxrange[1]]
        for j in range(len(arr)):
            mrk = 'square'
            sz = 6
            colm = 'blue'
            if col == 'rgb(218,112,214)':
                for i0 in range(len(arr0)):
                # if s0 starts before s1, and ends after s1 starts
                    if (arr0[i0] < arr[j]) and (arr[j] < arr0[i0] + lc < arr[j] + len(related[i])):
                        fig.add_trace(go.Scatter(x=[arr[j]+len(related[i])/2], y=[i+1.5],
                            hovertemplate = f'{related[i]} {arr[j]}',
                             showlegend=False, mode='markers',
                             marker=dict(size=8, color='green', symbol='triangle-right')))
                        mrk = 'triangle-right'
                        sz = 12
                        colm = 'green'
                # if s1 starts before s0, and ends after s0 starts
                    if (arr[j] < arr0[i0]) and (arr0[i0] <  arr[j] + len(related[i]) < arr0[i0] + lc):
                        fig.add_trace(go.Scatter(x=[arr[j]+len(related[i])/2], y=[i+1.5], hovertemplate = f'{related[i]} {arr[j]}',
                             showlegend=False,
                             mode='markers', marker=dict(size=6, color='green', symbol='triangle-left')))
                        mrk = 'triangle-left'
                        sz = 12
                        colm = 'green'
            fig.add_shape(type="rect", x0=arr[j], x1=arr[j]+len(related[i]), y0=i+1, y1=i+1.9,
                fillcolor=colt, line=dict(color=col, width=1))
            fig.add_trace(go.Scatter(x=[arr[j]+len(related[i])/2], y=[i+1.5], hovertemplate = f'{related[i]} {arr[j]}',
                             showlegend=False,
                             mode='markers', marker=dict(size=sz, color=colm, symbol=mrk)))
        if len(related[i]) > maxlenstr:
            maxlenstr = len(related[i])
        fig.add_annotation(
                    x=plotxrange[1]*1.001, y=i+1.35, text=related[i],
                    showarrow=False, xanchor="left",
                    font=dict(family="courier", size=14, color=col)
                )
    if (len(s1) == 0):
        i = 0

    # plot all the mutations:
    for j in range(len(arrmut)):
        arr = [p for p in find_all_matches(arrmut[j], txt.lower(), ret='pos')
               if plotxrange[0] <= p < plotxrange[1]]
        for idx in range(len(s0)):
            if arrmut[j][idx] != s0[idx]:
                arrmut[j] = arrmut[j][:idx] + arrmut[j][idx].upper() + arrmut[j][idx+1:]
            for k in range(len(arr)):
                mrk = 'square'
                sz = 6
                colm = 'blue'
                for i0 in range(len(arr0)):
                # if s0 starts before s1, and ends after s1 starts
                    if (arr0[i0] < arr[k]) and (arr[k] < arr0[i0] + lc < arr[k] + len(arrmut[j])):
                        mrk = 'triangle-right'
                        sz = 12
                        colm = 'darkgreen'
                # if s1 starts before s0, and ends after s0 starts
                    if (arr[k] < arr0[i0]) and (arr0[i0] <  arr[k] + len(arrmut[j]) < arr0[i0] + lc):
                        mrk = 'triangle-left'
                        sz = 12
                        colm = 'darkgreen'
                fig.add_shape(type="rect", x0=arr[k], x1=arr[k]+len(arrmut[j]), y0=-1-j, y1=-j-0.1,
                    fillcolor='rgba(255, 255, 255, 0.1)', line=dict(color='red', width=1))
                fig.add_trace(go.Scatter(x=[arr[k]+len(arrmut[j])/2], y=[-0.5-j], hovertemplate = f'{arrmut[j]} {arr[k]}',
                             showlegend=False,
                             mode='markers', marker=dict(size=sz, color=colm, symbol=mrk)))
        if len(arrmut[j]) > maxlenstr:
            maxlenstr = len(arrmut[j])
        fig.add_annotation(
                    x=plotxrange[1]*1.001, y=-0.55-j, text=arrmut[j],
                    showarrow=False, xanchor="left",
                    font=dict(family="courier", size=14, color='red')
                )

    # --- legend entries (dummy invisible traces) ---
    def _leg(name, color, symbol, size=8):
        """Return an invisible Scatter trace used solely as a legend entry."""
        return go.Scatter(x=[None], y=[None], mode='markers', name=name,
                          showlegend=True,
                          marker=dict(size=size, color=color, symbol=symbol))

    fig.add_trace(_leg(s0,                          'darkblue',          'x',              size=8))
    fig.add_trace(_leg('Neighboring core',           'rgb(218,112,214)',  'square',         size=8))
    fig.add_trace(_leg('Complete overlap with s0',          'rgb(0,200,0)',      'square',         size=8))
    if nm > 0:
        fig.add_trace(_leg('Mutation',               'red',              'square',         size=8))
    fig.add_trace(_leg('Partial overlap: core ends after s0',  'green',          'triangle-right', size=10))
    fig.add_trace(_leg('Partial overlap: core starts before s0', 'green',        'triangle-left',  size=10))
    if hairpins:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='lines', name='Hairpin region',
                                 showlegend=True, line=dict(color='darkorange', width=6)))

    avg_char_factor = 0.58  # adjust if your font is narrow/wide
    right_padding_px = 20   # extra breathing room
    left_margin_px = 60     # to fit y axis labels
    plot_area_width_px = 700  # desired interior plot area width (excluding margins)
    estimated_text_width_px = int(maxlenstr * 16 * avg_char_factor) + right_padding_px
    fig_width = left_margin_px + plot_area_width_px + estimated_text_width_px

    # custom y-axis ticks: 0 at s0, positive counts in both directions
    tick_vals  = [0.45] + [i+1.45 for i in range(n)] + [-0.55-j for j in range(nm)]
    tick_text  = ['0']  + [str(i+1) for i in range(n)] + [str(j+1) for j in range(nm)]
    fig.update_layout(
        width=fig_width, margin=dict(l=left_margin_px, r=estimated_text_width_px, t=40, b=80),
        legend=dict(orientation='h', yanchor='top', y=-0.05, xanchor='left', x=0),
        yaxis=dict(tickvals=tick_vals, ticktext=tick_text,
                   zeroline=True, zerolinewidth=2, zerolinecolor='lightgray'))
    fig.show()
    if file != '':
        out_file = file
        try:
            low = out_file.lower()
            if low.endswith(".html") or low.endswith(".htm"):
                fig.write_html(out_file, auto_open=False)
            else:
                root, ext = os.path.splitext(out_file)
                if ext == "":
                    out_file = f"{out_file}.png"
                fig.write_image(out_file, scale=scale)
            open_file_with_default_software(out_file)
        except Exception as e:
            print(f"Could not save/open output plot: {out_file}")
            print(str(e))


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
        try:
            low = file.lower()
            if low.endswith('.html') or low.endswith('.htm'):
                fig.write_html(file, auto_open=False)
            else:
                root, ext = os.path.splitext(file)
                if not ext:
                    file = f'{file}.png'
                fig.write_image(file, scale=scale)
            open_file_with_default_software(file)
        except Exception as e:
            print(f'Could not save plot: {file}\n{e}')


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


