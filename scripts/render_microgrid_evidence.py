#!/usr/bin/env python3
"""Generate a vector TikZ figure directly from the recorded representative run."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/microgrid"
with (OUT / "representative_trajectory.csv").open() as f:
    rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]


def coords(key):
    return " ".join(f"({r['hour_end']:.2f},{r[key]:.5f})" for r in rows)


def axes(ymax, ticks, ylabel, yscale):
    s = [rf"\begin{{tikzpicture}}[x=0.145cm,y={yscale}cm,font=\scriptsize]",
         rf"\draw[->] (0,0) -- (49,0) node[right] {{h}};",
         rf"\draw[->] (0,0) -- (0,{ymax+.08}) node[above] {{{ylabel}}};"]
    for t in [0, 12, 24, 36, 48]:
        s.append(rf"\draw[gray!30] ({t},0) -- ({t},{ymax}); \node[below] at ({t},0) {{{t}}};")
    for t in ticks:
        s.append(rf"\draw[gray!30] (0,{t}) -- (48,{t}); \node[left] at (0,{t}) {{{t:g}}};")
    return s


soc = axes(1., [.1, .3, .5, .7, .9], "State of charge", 3.)
soc += [rf"\draw[gray,dashed] plot coordinates {{{coords('min_soc')}}};",
        rf"\draw[gray,dashed] plot coordinates {{{coords('max_soc')}}};",
        rf"\draw[blue!70!black,thick] plot coordinates {{{coords('mean_soc')}}};",
        r"\draw[blue!70!black,thick] (6,-.18)--(10,-.18) node[right,black] {Mean};",
        r"\draw[gray,dashed] (25,-.18)--(29,-.18) node[right,black] {Min/max};",
        r"\end{tikzpicture}"]
power = axes(1.8, [0, .4, .8, 1.2, 1.6], "Mean power (MW)", 1.66667)
power += [rf"\draw[black,thick] plot coordinates {{{coords('mean_forecast_MW')}}};",
          rf"\draw[gray,dotted,thick] plot coordinates {{{coords('mean_cap_MW')}}};",
          rf"\draw[orange!80!black,dashed,thick] plot coordinates {{{coords('RN_mean_bid_MW')}}};",
          rf"\draw[blue!70!black,thick] plot coordinates {{{coords('retuned_mean_bid_MW')}}};",
          r"\draw[black,thick] (0,-.32)--(3,-.32) node[right] {Forecast};",
          r"\draw[gray,dotted,thick] (16,-.32)--(19,-.32) node[right,black] {Cap};",
          r"\draw[orange!80!black,dashed,thick] (28,-.32)--(31,-.32) node[right,black] {RN};",
          r"\draw[blue!70!black,thick] (37,-.32)--(40,-.32) node[right,black] {Retuned};",
          r"\end{tikzpicture}"]
content = [r"\begin{figure*}[t]", r"\centering", *soc, r"\hspace{0.7cm}", *power,
           r"\caption{Sequential microgrid operation for $n=50$, $s=1$, and replicate 0. Left: battery state of charge under the prescribed local controller. Right: mean conditional export forecast, availability limit, and bidding quantities. All bidding strategies use the same local energy-management policy and physical trajectory.}",
           r"\label{fig:microgrid_trajectory}", r"\end{figure*}"]
(OUT / "trajectory_plot.tex").write_text("\n".join(content)+"\n")
standalone = [r"\documentclass{article}", r"\usepackage[paperwidth=18.5cm,paperheight=5.1cm,margin=0.2cm]{geometry}",
              r"\usepackage[T1]{fontenc}",r"\renewcommand{\rmdefault}{ptm}",r"\usepackage{tikz}",r"\pagestyle{empty}",r"\begin{document}\noindent",
              *soc, r"\hspace{0.7cm}", *power, r"\end{document}"]
(OUT / "trajectory_figure.tex").write_text("\n".join(standalone)+"\n")
print(OUT / "trajectory_plot.tex")
