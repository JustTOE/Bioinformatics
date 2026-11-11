import sys
import urllib.request

# TP53 mRNA, human
accession = "NM_000546.6"
if len(sys.argv) > 1:
    accession = sys.argv[1]

url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       "?db=nuccore&rettype=fasta&retmode=text&id=" + accession)
data = urllib.request.urlopen(url).read().decode()

lines = data.strip().split("\n")
header = lines[0][1:]
seq = "".join(lines[1:]).upper()

if len(seq) < 1000 or len(seq) > 3000:
    print("length is", len(seq), "- need something between 1000 and 3000")
    sys.exit(1)

with open(accession + ".fasta", "w") as f:
    f.write(">" + header + "\n")
    for i in range(0, len(seq), 70):
        f.write(seq[i:i+70] + "\n")

gc = (seq.count("G") + seq.count("C")) / len(seq) * 100
print(header)
print("length:", len(seq))
print("GC: %.2f%%" % gc)
for base in "ACGT":
    print(base, seq.count(base))
