"""Plot finite-tree probabilities by positive Bernoulli convolution (no simulation)."""
from math import ceil, exp, floor, log, sqrt, erfc
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.7,
                     'mathtext.fontset': 'stix', 'pdf.fonttype': 42,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.65, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6})


def finite_curves(n):
    length = log(n)
    k_center = floor(length)
    max_k = ceil(2.8 * length)
    # The terminal bin holds the whole remaining tail.
    mass = np.zeros(max_k + 3)
    mass[0] = 1
    jmax = min(n, ceil(exp(2.55 * sqrt(length))))
    upper = np.zeros(jmax)
    for i in range(n, 1, -1):
        if i <= jmax:
            upper[i - 1] = mass[k_center + 1:].sum()
        p = 1.0 / i
        mass[-1] += p * mass[-2]
        mass[1:-1] = (1.0 - p) * mass[1:-1] + p * mass[:-2]
        mass[0] *= 1.0 - p
    upper[0] = mass[k_center + 1:].sum()
    tails = np.cumsum(mass[::-1])[::-1]
    assert abs(tails[0] - 1.0) < 2e-12
    assert np.all(np.diff(upper) <= 1e-13)
    ks = np.arange(1, max_k + 1)
    x = np.log(np.arange(1, jmax + 1)) / sqrt(length)
    return ks / length, mass[ks] / tails[ks], x, upper / tails[k_center]


def draw():
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.18))
    fig.subplots_adjust(left=.085, right=.985, top=.87, bottom=.24, wspace=.32)
    agrid = np.linspace(.02, 2.8, 501)
    ygrid = np.linspace(0, 2.55, 401)
    for s, color, marker, style in [(6, '#b95125', 'o', '-'), (12, '#176287', 's', '-.')]:
        n = ceil(exp(s))
        a, root, y, surv = finite_curves(n)
        axes[0].plot(a, root, color=color, lw=1.15, marker=marker, ms=2.2, ls=style,
                     label=f'n = {n:,}')
        axes[1].step(y, surv, where='post', color=color, lw=1.15, ls=style)
    axes[0].plot(agrid, np.maximum(0, 1 - 1 / agrid), color='#242424',
                 lw=1.1, ls=(0,(4,2)), label='Limit')
    axes[1].plot(ygrid, [erfc(y / sqrt(2)) for y in ygrid],
                 color='#242424', lw=1.1, ls=(0,(4,2)))
    axes[0].set(xlim=(0, 2.8), ylim=(-.015, .78), xlabel=r'Prescribed distance $k/\log n$',
                ylabel='Conditional root probability', title='(a)  The common ancestor is the root')
    axes[1].set(xlim=(0,2.55), ylim=(0,1.035), xlabel=r'Scaled birth rank $y$',
                ylabel='Conditional upper tail', title='(b)  The critical distance')
    for ax in axes:
        ax.title.set_fontsize(9.0)
        ax.grid(axis='y', color='#c8ced0', alpha=.65, lw=.45)
        ax.set_axisbelow(True)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False,
               bbox_to_anchor=(.5,-.018), handlelength=2.5)
    destination = ROOT / 'figures' / 'conditional_laws.pdf'
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, metadata={'Title':'Conditional common-ancestor probabilities',
                                        'CreationDate': None, 'ModDate': None})
    plt.close(fig)


if __name__ == '__main__':
    draw()
