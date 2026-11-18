#!/usr/bin/env python3
"""Pipeline outputs (out/SCEN) vs GATK's (truth/SCEN), by file name.  usage: compare.py germline|somatic
BAM: records and header apart from @PG; VCF: records and header apart from command / source / date lines;
archives: members' contents (F1R2 histograms by column name); tables and metrics: lines apart from '#' comments (Picard's command line)."""
import glob, gzip, hashlib, os, subprocess, sys, tarfile
os.chdir(os.path.dirname(os.path.abspath(__file__)))
scen = sys.argv[1]
A, T = f"out/{scen}", f"truth/{scen}"
def bam(p, header):
    out = subprocess.run(["samtools", "view", "--no-PG"] + (["-H"] if header else []) + [p], capture_output=True, check=True).stdout
    if header:
        out = b"\n".join(l for l in out.split(b"\n") if not l.startswith(b"@PG"))
    return hashlib.md5(out).hexdigest(), out.count(b"\n")
def lines(p):
    with (gzip.open if p.endswith(".gz") else open)(p, "rt") as f:
        return f.read().splitlines()
SKIP = ("##GATKCommandLine", "##variant-rs", "##source=", "##fileDate", "##filtering_status")
def cmp(a, t):
    if a.endswith(".bam"):
        (ra, na), (rt, nt) = bam(a, False), bam(t, False)
        return (f"identical ({na} records)" if ra == rt else f"DIFFERENT records ({na} vs {nt})") + ("" if bam(a, True)[0] == bam(t, True)[0] else "; header differs")
    if a.endswith(".tar.gz"):
        # F1R2 histograms: compared column by column (by name): GATK's column order follows Java identity hash codes
        # (a HashMap keyed by enum pairs) and differs between GATK runs
        def member(name, data):
            if not name.endswith("_histogram"):
                return data
            L = data.decode().split("\n")
            k = next(i for i, l in enumerate(L) if l.startswith("depth") or l.startswith("BIN") or l.startswith("set_size"))
            rows = [l.split("\t") for l in L[k:] if l]
            return (tuple(L[:k]), sorted(zip(*rows)))
        m = lambda p: {os.path.basename(x.name): member(x.name, tarfile.open(p).extractfile(x).read()) for x in tarfile.open(p).getmembers() if x.isfile()}
        return "identical contents" if m(a) == m(t) else "DIFFERENT contents"
    if ".vcf" in a:
        la, lt = [[l for l in lines(p) if not l.startswith(SKIP)] for p in (a, t)]
        n = sum(1 for l in la if not l.startswith("#"))
        return f"identical ({n} records)" if la == lt else f"DIFFERENT ({n} vs {sum(1 for l in lt if not l.startswith('#'))} records)"
    la, lt = [[l for l in lines(p) if not l.startswith("#")] for p in (a, t)]
    return "identical" if la == lt else "DIFFERENT"
rows = []
for a in sorted(glob.glob(f"{A}/**/*", recursive=True)):
    if not os.path.isfile(a) or "/lanes/" in a or a.endswith((".bai", ".tbi")):
        continue
    t = os.path.join(T, os.path.basename(a))
    rows.append((os.path.relpath(a, A), cmp(a, t) if os.path.exists(t) else "no GATK file"))
w = max(len(r[0]) for r in rows)
for r in rows:
    print(f"  {r[0]:{w}s}  {r[1]}")
bad = [r for r in rows if not r[1].startswith("identical")]
print(f"{len(rows)} files, {len(bad)} different or missing")
