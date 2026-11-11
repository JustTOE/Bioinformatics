import os
import sys
import math
import urllib.request
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from enzymes import enzymes

# SARS-CoV-2 (Wuhan-Hu-1) and SARS-CoV (Tor2)
genomes = ["NC_045512.2", "NC_004718.3"]
if len(sys.argv) == 3:
    genomes = sys.argv[1:3]

ladder = [10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250]

# two bands count as "the same" if their sizes are within 2%
TOLERANCE = 0.02


def get_genome(acc):
    filename = acc + ".fasta"
    if not os.path.exists(filename):
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
               "?db=nuccore&rettype=fasta&retmode=text&id=" + acc)
        data = urllib.request.urlopen(url).read().decode()
        with open(filename, "w") as f:
            f.write(data)
    with open(filename) as f:
        lines = f.read().strip().split("\n")
    return lines[0][1:], "".join(lines[1:]).upper()


def digest(seq):
    cuts = []
    for name, site, cleavage in enzymes:
        offset = cleavage.index("^")
        i = seq.find(site)
        while i != -1:
            cuts.append((i + offset, name))
            i = seq.find(site, i + 1)
    cuts.sort()

    # the genomes are linear, so n cuts -> n+1 fragments
    fragments = []
    prev = 0
    for pos, name in cuts:
        fragments.append(pos - prev)
        prev = pos
    fragments.append(len(seq) - prev)
    return cuts, sorted(fragments, reverse=True)


def distance(bp):
    # semi-log: distance drops linearly with log10(size)
    # 10 kb runs 1 cm, 250 bp runs 9 cm
    bp = min(max(bp, 100), 20000)
    slope = 8 / (math.log10(10000) - math.log10(250))
    return 1 + slope * (math.log10(10000) - math.log10(bp))


def compare(a, b):
    left = list(b)
    shared = []
    only_a = []
    for x in a:
        match = None
        for y in left:
            if abs(x - y) / max(x, y) <= TOLERANCE:
                match = y
                break
        if match is None:
            only_a.append(x)
        else:
            shared.append((x, match))
            left.remove(match)
    return shared, only_a, left


def draw_band(ax, x, y, color):
    ax.add_patch(patches.Rectangle((x - 0.32, y - 0.06), 0.64, 0.12,
                                   facecolor=color, alpha=0.25, edgecolor="none"))
    ax.add_patch(patches.Rectangle((x - 0.28, y - 0.03), 0.56, 0.06,
                                   facecolor=color, edgecolor="none"))


def plot_gel(names, frags, shared_sizes, out="gel_comparison.png"):
    fig, ax = plt.subplots(figsize=(8, 9), facecolor="#111")
    ax.set_facecolor("#1a1a1a")

    for x in range(3):
        ax.add_patch(patches.Rectangle((x + 0.65, 0.25), 0.7, 0.3, facecolor="#333"))

    for bp in ladder:
        draw_band(ax, 1, distance(bp), "#7ec8ff")
        ax.text(0.55, distance(bp), str(bp), color="#7ec8ff", fontsize=7,
                ha="right", va="center")

    for lane, fr in enumerate(frags):
        x = lane + 2
        last_label = -1
        for bp in fr:
            # orange = band found in both lanes, white = only in this lane
            color = "#ffb347" if bp in shared_sizes[lane] else "#ffffff"
            y = distance(bp)
            draw_band(ax, x, y, color)
            # push the label down a bit if it would sit on top of the previous one
            label_y = max(y, last_label + 0.2)
            ax.text(x + 0.38, label_y, str(bp), color=color, fontsize=7, va="center")
            last_label = label_y

    ax.set_xlim(0, 4)
    ax.set_ylim(10, 0)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["Ladder (bp)"] + names, color="white", fontsize=9)
    ax.set_ylabel("distance from well (cm)", color="white")
    ax.tick_params(colors="white")
    ax.set_title("BamHI + SmaI + PstI digest\n(orange = shared band, white = different)",
                 color="white")
    plt.tight_layout()
    plt.savefig(out, dpi=200, facecolor=fig.get_facecolor())
    plt.close()
    print("saved", out)


names = []
frags = []
for acc in genomes:
    header, seq = get_genome(acc)
    cuts, fr = digest(seq)
    names.append(acc)
    frags.append(fr)

    print(header)
    print("length:", len(seq), "bp")
    for name, site, cleavage in enzymes:
        n = len([c for c in cuts if c[1] == name])
        print("  %-6s %d site(s)" % (name, n))
    print("  fragments:", fr)
    print()

shared, only_a, only_b = compare(frags[0], frags[1])

print("bands in both lanes:", len(shared))
for x, y in shared:
    print("  %6d ~ %d" % (x, y))
print("only in", names[0] + ":", only_a)
print("only in", names[1] + ":", only_b)

similarity = 2 * len(shared) / (len(frags[0]) + len(frags[1]))
print("band similarity: %.0f%%" % (similarity * 100))
if not only_a and not only_b:
    print("-> the two lanes look the same")
else:
    print("-> the two lanes are different")

shared_sizes = [[x for x, y in shared], [y for x, y in shared]]
plot_gel(names, frags, shared_sizes)
