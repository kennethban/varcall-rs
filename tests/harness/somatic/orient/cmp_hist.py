#!/usr/bin/env python3
"""Compare two F1R2 metrics files (histograms) as label -> column maps; prints the first differences."""
import sys
def read(p):
    L = open(p).read().split('\n')
    i = next(k for k, l in enumerate(L) if l.startswith('## HISTOGRAM'))
    labels = L[i + 1].split('\t')[1:]
    cols = {l: [] for l in labels}
    for l in L[i + 2:]:
        if not l: break
        f = l.split('\t')
        for lab, v in zip(labels, f[1:]): cols[lab].append(v)
    return L[:i], labels, cols
ha, la, a = read(sys.argv[1]); hb, lb, b = read(sys.argv[2])
same_order = la == lb
diff = [k for k in set(a) | set(b) if a.get(k) != b.get(k)]
print(("identical" if ha == hb and not diff else "DIFFER") + ("" if same_order else " (column order differs)"), len(diff), "columns differ", sorted(diff)[:5])
