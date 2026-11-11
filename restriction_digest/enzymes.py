import sys

# name, recognition seq, cleavage site (^ = cut on top strand)
enzymes = [
    ("BamHI", "GGATCC", "G^GATCC"),
    ("SmaI",  "CCCGGG", "CCC^GGG"),
    ("PstI",  "CTGCAG", "CTGCA^G"),
]

if __name__ == "__main__":
    filename = "NM_000546.6.fasta"
    if len(sys.argv) > 1:
        filename = sys.argv[1]

    with open(filename) as f:
        lines = f.read().strip().split("\n")
    header = lines[0][1:]
    seq = "".join(lines[1:]).upper()

    print(header, "(" + str(len(seq)) + " nt)")
    print()
    print("%-7s %-12s %-10s %-8s %s" % ("Enzyme", "Recognition", "Cleavage", "Site", "Cut after"))

    for name, site, cleavage in enzymes:
        offset = cleavage.index("^")
        # the bottom strand is cut at the mirrored position
        if offset == len(site) - offset:
            ends = "blunt"
        elif offset < len(site) - offset:
            ends = "5' overhang"
        else:
            ends = "3' overhang"

        positions = []
        i = seq.find(site)
        while i != -1:
            positions.append(i)
            i = seq.find(site, i + 1)

        if not positions:
            print("%-7s %-12s %-10s %-8s %s" % (name, site, cleavage, "-", "-"))
        for p in positions:
            print("%-7s %-12s %-10s %-8d %d   (%s)" % (name, site, cleavage, p + 1, p + offset, ends))
