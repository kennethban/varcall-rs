#!/usr/bin/env python3
"""Prototype of the GenomicsDB 1.5.5 query/merge as seen by GATK GenotypeGVCFs (see genomicsdb_merge.md).

Rebuilds received/<ws>.vcf (the VariantContexts the GenotypeGVCFs walker receives from -V gendb://<ws>, written by
htsjdk) from the per-sample input GVCFs listed in the workspace's sample map, then compares record by record.

usage: merge_proto.py <ws> [--map rb.map] [--interval chr8:1-15000000] [--out proto/<ws>.vcf] [--no-compare]
"""
import argparse, bisect, gzip, heapq, json, math, os, struct, sys
from collections import Counter, OrderedDict

# workspaces, sample maps, received/ and proto/: the harness data folder $VARCALL_DATA/genomicsdb
HERE = os.path.join(os.environ.get("VARCALL_DATA", "."), "genomicsdb")
REF_FA = os.environ.get("REF_FA") or os.path.join(os.environ.get("VARCALL_DATA", "."), "bqsr/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa")
NON_REF = "<NON_REF>"
SHARD = 100_000            # VariantWalkerBase.DEFAULT_DRIVING_VARIANTS_LOOKAHEAD_BASES (shard size and cache lookahead)
MAX_ALTS = 50              # GenomicsDBArgumentCollection default (htsjdk GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED)
MAX_GENOTYPE_COUNT = 1024  # GenotypeCalculationArgumentCollection default
WS_DEFAULTS = {"ws_rb": ("rb.map", "chr8:1-15000000"), "ws_raw": ("raw.map", "chr8:1-15000000"),
               "ws_rb_b1": ("rb.map", "chr8:1-15000000"), "ws11": ("rb11.map", "chr8:1-5000000")}
# GenomicsDB known fields with a default INFO combine operation (known_field_info.cc)
KNOWN_DEFAULT_OP = {"BaseQRankSum": "median", "ClippingRankSum": "median", "MQRankSum": "median",
                    "ReadPosRankSum": "median", "MQ": "median", "MQ0": "median", "ExcessHet": "median",
                    "RAW_MQ": "sum", "DP": "DP"}


# ----------------------------------------------------------------------------------------------- numeric helpers
def f32(x):
    return struct.unpack("f", struct.pack("f", x))[0]


def kputd(d):
    """htslib kstring.c kputd(): how htslib prints a float INFO/FORMAT value (float32 promoted to double)."""
    if d == 0:
        return "-0" if math.copysign(1.0, d) < 0 else "0"
    if math.isnan(d):
        return "nan"
    neg = d < 0
    if neg:
        d = -d
    if not (0.0001 <= d <= 999999):
        return ("-" if neg else "") + ("%g" % d)
    i = int(d * 10000000000)
    if d < .0001: i += 0
    elif d < 0.001: i += 5
    elif d < 0.01: i += 50
    elif d < 0.1: i += 500
    elif d < 1: i += 5000
    elif d < 10: i += 50000
    elif d < 100: i += 500000
    elif d < 1000: i += 5000000
    elif d < 10000: i += 50000000
    elif d < 100000: i += 500000000
    else: i += 5000000000
    digits = str(i)
    p = len(digits)
    if p <= 10:  # d < 1: keep 6 significant digits
        s = "0." + "0" * (10 - p) + digits[:6]
    else:
        s = (digits[:p - 10] + "." + digits[p - 10:])[:7]
        if s[6:7] == ".":
            s = s[:6]
    if '.' in s:  # cull trailing zeros, but an all-zero fraction keeps one zero ("2.0")
        s = s.rstrip('0')
        if s.endswith('.'):
            s += '0'

    return ("-" if neg else "") + s


def fmt3(x):
    return "%.3f" % x


# ----------------------------------------------------------------------------------------------- reference
class Fasta:
    def __init__(self, path):
        self.f = open(path, "rb")
        self.idx = {}
        for line in open(path + ".fai"):
            name, ln, off, lb, lbytes = line.split("\t")[:5]
            self.idx[name] = (int(ln), int(off), int(lb), int(lbytes))
        self.cache = {}

    def base(self, contig, pos1):
        ln, off, lb, lbytes = self.idx[contig]
        z = pos1 - 1
        self.f.seek(off + (z // lb) * lbytes + z % lb)
        return chr(self.f.read(1)[0])


# ----------------------------------------------------------------------------------------------- vid map
class Field:
    def __init__(self, d):
        self.name = d["name"]
        self.types = d.get("type", [])
        self.classes = d.get("vcf_field_class", [])
        self.lengths = [l.get("variable_length_descriptor") for l in d.get("length", [])]
        self.delims = d.get("vcf_delimiter")
        op = d.get("VCF_field_combine_operation")
        if op is None:
            op = KNOWN_DEFAULT_OP.get(self.name, "unknown")
        self.op = op
        self.remap_missing_with_non_ref = not d.get("disable_remap_missing_with_non_ref", False)
        self.two_d = len(self.lengths) == 2
        self.is_float = self.types and self.types[0].lower() == "float"


def load_vid(ws):
    v = json.load(open(os.path.join(HERE, ws, "vidmap.json")))
    return OrderedDict((f["name"], Field(f)) for f in v["fields"])


# ----------------------------------------------------------------------------------------------- input GVCFs
class Rec:
    __slots__ = ("begin", "end", "ref", "alts", "info", "fmt", "row")

    def __init__(self, begin, end, ref, alts, info, fmt, row):
        self.begin, self.end, self.ref, self.alts, self.info, self.fmt, self.row = begin, end, ref, alts, info, fmt, row


def read_gvcf(path, contig, lo, hi, row):
    recs = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            t = line.rstrip("\n").split("\t")
            if t[0] != contig:
                continue
            pos = int(t[1])
            if pos > hi:
                break
            ref = t[3]
            alts = t[4].split(",") if t[4] != "." else []
            info = {}
            if t[7] != ".":
                for kv in t[7].split(";"):
                    k, _, v = kv.partition("=")
                    info[k] = v
            end = int(info["END"]) if "END" in info else pos + len(ref) - 1
            if end < lo:
                continue
            fmt = dict(zip(t[8].split(":"), t[9].split(":"))) if len(t) > 9 else {}
            recs.append(Rec(pos, end, ref, alts, info, fmt, row))
    return recs


def parse_num(s, is_float):
    if s == "." or s == "":
        return None
    if is_float:
        try:
            return f32(float(s))
        except ValueError:
            return None
    try:
        return int(s)
    except ValueError:
        return int(float(s))


def parse_vec(s, is_float):
    if s is None or s == ".":
        return None
    return [parse_num(x, is_float) for x in s.split(",")]


def parse_2d(s, fld):
    """2D allele-specific field ('|' between alleles, ',' inside). Returns list (per allele) of element lists."""
    if s is None:
        return None
    out = []
    for ent in s.split("|"):
        if ent == "":
            out.append([])
        else:
            out.append([parse_num(x, fld.is_float) for x in ent.split(",")])
    return out


def parse_hist(s):
    """histogram 2D field (float bin, int count) pairs per allele."""
    if s is None:
        return None
    out = []
    for ent in s.split("|"):
        pairs = []
        if ent != "":
            el = ent.split(",")
            for i in range(0, len(el), 2):
                try:
                    b = f32(float(el[i]))
                except ValueError:
                    b = None
                c = None
                if i + 1 < len(el):
                    try:
                        c = int(el[i + 1])
                    except ValueError:
                        c = None
                if b is not None and math.isnan(b):
                    b = None
                pairs.append((b, c))
        out.append(pairs)
    return out


def parse_gt(s):
    if s is None or s == "" or s == ".":
        return None if s is None else [(None, False)]
    out = []
    i = 0
    phased = False
    cur = ""
    for ch in s:
        if ch in "/|":
            out.append((None if cur == "." else int(cur), phased))
            phased = ch == "|"
            cur = ""
        else:
            cur += ch
    out.append((None if cur == "." else int(cur), phased))
    return out


def is_symbolic(a):
    return a == NON_REF or a == "*" or (len(a) > 0 and ((a[0] == "<" and a[-1] == ">") or "[" in a or "]" in a))


def is_deletion(ref, alt):
    return len(ref) > 1 and (alt == "*" or (not is_symbolic(alt) and len(alt) < len(ref)))


def is_mnv(ref, alt):
    return (len(ref) > 1 and not is_symbolic(alt) and len(alt) == len(ref)
            and sum(1 for a, b in zip(ref, alt) if a != b) > 1)


def gt_index(alleles):
    """VCF genotype index of an allele combination (any ploidy)."""
    a = sorted(alleles)
    idx = 0
    for k, x in enumerate(a):
        # number of combinations with ploidy k+1 and max allele < x
        idx += math.comb(x + k, k + 1)
    return idx


def num_genotypes(nalts, ploidy):
    return math.comb(nalts + 1 + ploidy - 1, ploidy)


# ----------------------------------------------------------------------------------------------- the merger
class Call:
    """A VariantCall: a stored cell filled into a row (fresh state each time the cell is read by a query)."""
    __slots__ = ("rec", "row", "begin", "end", "ref", "alts", "deletion_or_mnv", "f", "info_invalid")

    def __init__(self, rec, M):
        self.rec, self.row, self.begin, self.end = rec, rec.row, rec.begin, rec.end
        self.ref = rec.ref
        self.alts = list(rec.alts)
        self.deletion_or_mnv = any(is_deletion(rec.ref, a) for a in rec.alts) or any(is_mnv(rec.ref, a) for a in rec.alts)
        self.info_invalid = False
        f = {}
        fm = rec.fmt
        f["GT"] = parse_gt(fm.get("GT"))
        for k in ("AD", "PL", "SB"):
            if k in M.vid:
                f[k] = parse_vec(fm.get(k), False)
        for k in ("GQ", "MIN_DP", "PS"):
            if k in M.vid:
                v = fm.get(k)
                f[k] = None if v is None else parse_num(v, False)
        f["DP_FORMAT"] = None if fm.get("DP") is None else parse_num(fm["DP"], False)
        for k in ("PGT", "PID"):
            if k in M.vid:
                v = fm.get(k)
                f[k] = None if (v is None or v == ".") else v
        info = rec.info
        f["DP"] = None if info.get("DP") is None else parse_num(info["DP"], False)
        for fld in M.info_fields:
            v = info.get(fld.name)
            if fld.op == "histogram_sum":
                f[fld.name] = parse_hist(v)
            elif fld.two_d:
                f[fld.name] = parse_2d(v, fld)
            elif fld.op in ("element_wise_sum",):
                f[fld.name] = parse_vec(v, fld.is_float)
            else:
                f[fld.name] = parse_num(v, fld.is_float) if v is not None else None
        self.f = f


class Merger:
    def __init__(self, ws, samples, fasta, contig):
        self.ws = ws
        self.vid = load_vid(ws)
        self.samples = samples
        self.n = len(samples)
        self.fasta = fasta
        self.contig = contig
        # INFO fields handled by BroadCombinedGVCFOperator (DP handled separately)
        self.info_fields = [f for f in self.vid.values() if "INFO" in f.classes and f.name not in ("END", "DP")
                            and f.op in ("sum", "median", "mean", "element_wise_sum", "histogram_sum", "concatenate")]
        self.has_fmt = {k: (k in self.vid) for k in ("AD", "PL", "SB", "GQ", "MIN_DP", "PS", "PGT", "PID", "GT")}

    # --- handle_deletions (BroadCombinedGVCFOperator::handle_deletions)
    def handle_deletions(self, s, calls):
        for c in calls:
            if c is None or not c.deletion_or_mnv or not (s > c.begin):
                continue
            c.info_invalid = True  # INFO fields of m_INFO_fields_vec (all non-histogram combine ops, not DP)
            alts = c.alts
            if alts[0] == "*" and (len(alts) == 1 or (len(alts) == 2 and alts[1] == NON_REF)):
                continue
            f = c.f
            gt = f["GT"]
            ploidy = len(gt) if gt else 0
            pl = f.get("PL")
            lowest = -1
            lowest_pl = 1 << 62
            non_ref_idx = None
            for i, a in enumerate(alts):
                ai = i + 1
                if is_deletion(c.ref, a) or is_mnv(c.ref, a):
                    if lowest < 0:
                        lowest = ai
                    if pl is not None:
                        gi = gt_index([ai] * ploidy)
                        if gi < len(pl) and pl[gi] is not None and pl[gi] < lowest_pl:
                            lowest_pl = pl[gi]
                            lowest = ai
                elif a == NON_REF:
                    non_ref_idx = ai
            has_nr = non_ref_idx is not None
            lut = {0: 0, lowest: 1}   # input -> reduced
            inv = {0: 0, 1: lowest}   # reduced -> input
            if has_nr:
                lut[non_ref_idx] = 2
                inv[2] = non_ref_idx
            nred = 3 if has_nr else 2
            # remap allele-dependent fields (1D: AD (R), PL (G); 2D: AS_* R fields incl. histograms)
            if f.get("AD") is not None:
                f["AD"] = remap_alleles(f["AD"], inv, nred, non_ref_idx if has_nr else None, False)
            if f.get("PL") is not None:
                f["PL"] = remap_genotypes(f["PL"], inv, nred, non_ref_idx if has_nr else None, ploidy)
            for fld in self.info_fields:
                if fld.two_d and f.get(fld.name) is not None:
                    f[fld.name] = remap_2d(f[fld.name], inv, nred, non_ref_idx if has_nr else None, fld.lengths[0] == "A")
            if gt:
                f["GT"] = [(a if a is None else lut.get(a, 2 if has_nr else None), ph) for a, ph in gt]
            c.ref = "N"
            c.alts = ["*", NON_REF] if has_nr else ["*"]

    def operate(self, s, e, calls):
        self.handle_deletions(s, calls)
        active = [c for c in calls if c is not None]
        # merge_reference_allele
        mref = "N"
        for c in active:
            if c.begin < s:
                continue
            r = c.ref
            if len(r) > len(mref):
                mref = r if mref[0] == "N" else mref + r[len(mref):]
            elif mref[0] == "N" and r[0] != "N":
                mref = r
        # merge_alt_alleles
        seen = {}
        malts = []
        luts = {}
        nonref_in = {}
        nr_exists = False
        for c in active:
            lut = {0: 0}
            suffix = mref[len(c.ref):] if len(c.ref) < len(mref) else ""
            for i, a in enumerate(c.alts):
                ai = i + 1
                if a == NON_REF:
                    nonref_in[c.row] = ai
                    nr_exists = True
                    continue
                if suffix and not is_symbolic(a):
                    a = a + suffix
                if a not in seen:
                    malts.append(a)
                    seen[a] = len(malts)
                lut[ai] = seen[a]
            luts[c.row] = lut
        if nr_exists:
            malts.append(NON_REF)
            nri = len(malts)
            for c in active:
                if c.row in nonref_in:
                    luts[c.row][nonref_in[c.row]] = nri
        nall = len(malts) + 1
        ref_block_only = len(mref) == 1 and malts == [NON_REF]
        remap = not ref_block_only
        too_many = len(malts) > MAX_ALTS
        # reference base for 'N'
        if mref == "N":
            b = self.fasta.base(self.contig, s)
            mref = b.upper() if b in "ACGTacgt" else "N"
        # inverse LUTs: merged -> input
        inv = {}
        for c in active:
            d = {}
            for k, v in luts[c.row].items():  # LUT merged->input: last added pair wins
                d[v] = k
            inv[c.row] = d
        # remapped per-call data
        R = {}
        for c in active:
            f = c.f
            r = {}
            gt = f["GT"]
            nr_in = nonref_in.get(c.row)
            if remap and gt is not None:
                lut = luts[c.row]
                r["GT"] = [(a if a is None else lut.get(a, (nall - 1) if (nr_exists) else None), ph) for a, ph in gt]
            else:
                r["GT"] = gt
            ploidy = len(gt) if gt else 0
            for k in ("AD", "PL"):
                v = f.get(k)
                if v is None or not remap:
                    r[k] = v
                    continue
                if k == "PL":
                    if too_many:
                        r[k] = None
                        continue
                    if num_genotypes(nall - 1, ploidy) > MAX_GENOTYPE_COUNT:
                        r[k] = None
                        continue
                    r[k] = remap_genotypes(v, inv[c.row], nall, nr_in if nr_exists else None, ploidy)
                else:
                    r[k] = remap_alleles(v, inv[c.row], nall, nr_in if nr_exists else None, False)
            R[c.row] = r
        # ---------------- INFO
        info = {}
        if e > s:
            info["END"] = str(e)
        for fld in self.info_fields:
            name = fld.name
            vals = []
            for c in active:
                if c.info_invalid:
                    continue
                v = c.f.get(name)
                if v is None:
                    continue
                vals.append((c, v))
            if not vals:
                continue
            if fld.op == "sum":
                good = [v for c, v in vals if v is not None]
                if good:
                    info[name] = str(sum(good)) if not fld.is_float else kputd(f32(sum(good)))
            elif fld.op == "median":
                good = sorted(v for c, v in vals if v is not None)
                if good:
                    info[name] = kputd(good[len(good) // 2])
            elif fld.op == "element_wise_sum" and not fld.two_d:
                res = []
                for c, v in vals:
                    for i, x in enumerate(v):
                        if i >= len(res):
                            res.append(x)
                        elif x is not None:
                            res[i] = x if res[i] is None else res[i] + x
                info[name] = ",".join("." if x is None else (kputd(f32(x)) if fld.is_float else str(x)) for x in res)
            elif fld.op == "element_wise_sum":
                res = []
                any_valid = False
                for c, v in vals:
                    if remap:
                        v = remap_2d(v, inv[c.row], nall, nonref_in.get(c.row) if (nr_exists and fld.remap_missing_with_non_ref) else None,
                                     fld.lengths[0] == "A")
                    while len(res) < len(v):
                        res.append([])
                    for i, ent in enumerate(v):
                        tgt = res[i]
                        while len(tgt) < len(ent):
                            tgt.append(None)
                        for j, x in enumerate(ent):
                            if x is None:
                                continue
                            any_valid = True
                            tgt[j] = x if tgt[j] is None else (f32(tgt[j] + x) if fld.is_float else tgt[j] + x)
                if any_valid:
                    info[name] = "|".join(",".join("" if x is None else (fmt3(x) if fld.is_float else str(x)) for x in ent) for ent in res)
            elif fld.op == "histogram_sum":
                res = []
                for c, v in vals:
                    if remap:
                        v = remap_2d(v, inv[c.row], nall, nonref_in.get(c.row) if (nr_exists and fld.remap_missing_with_non_ref) else None, False)
                    while len(res) < len(v):
                        res.append({})
                    for i, ent in enumerate(v):
                        h = res[i]
                        for b, cnt in ent:
                            if b is None or cnt is None:
                                continue
                            h[b] = h.get(b, 0) + cnt
                info[name] = "|".join(",".join("%s,%d" % (fmt3(b), h[b]) for b in sorted(h)) for h in res)
        # ---------------- FORMAT
        n = self.n
        bycall = [None] * n
        for c in active:
            bycall[c.row] = c
        gq_zero = set()
        pl_zero = set()
        fmt = {}  # name -> list per sample (None = missing)
        # GT
        if any(bycall[i] is not None and R[i]["GT"] is not None for i in range(n)):
            fmt["GT"] = [R[i]["GT"] if bycall[i] is not None else None for i in range(n)]
        for k in ("AD", "GQ", "MIN_DP", "PGT", "PID", "PL", "PS", "SB"):
            if not self.has_fmt.get(k):
                continue
            if k == "PL" and too_many:
                continue
            col = []
            for i in range(n):
                c = bycall[i]
                if c is None:
                    col.append(None)
                elif k in ("AD", "PL"):
                    col.append(R[i][k])
                else:
                    col.append(c.f.get(k))
            if all(v is None for v in col):
                continue
            fmt[k] = col
            if k == "GQ":
                gq_zero = {i for i, v in enumerate(col) if v == 0}
            if k == "PL":
                pl_zero = {i for i, v in enumerate(col) if v is not None and len(v) > 0 and v[0] == 0}
        if gq_zero and pl_zero and "GT" in fmt:
            gts = fmt["GT"]
            for i in gq_zero & pl_zero:
                if gts[i] is not None:
                    gts[i] = [(None, False)] * len(gts[i])
        # DP
        dpf = [bycall[i].f["DP_FORMAT"] if bycall[i] is not None else None for i in range(n)]
        dpi = [bycall[i].f["DP"] if bycall[i] is not None else None for i in range(n)]
        mdp = [bycall[i].f.get("MIN_DP") if bycall[i] is not None else None for i in range(n)] if "MIN_DP" in fmt else [None] * n
        if any(v is not None for v in dpf) or any(v is not None for v in dpi):
            tot = 0
            for i in range(n):
                v = dpi[i]
                if v is None:
                    v = mdp[i] if mdp[i] is not None else dpf[i]
                if v is not None:
                    tot += v
            if any(v is not None for v in dpf):
                fmt["DP"] = dpf
            if tot > 0 and not ref_block_only:
                info["DP"] = str(tot)
        return self.htsjdk_line(s, mref, malts, info, fmt)

    def htsjdk_line(self, s, ref, alts, info, fmt):
        """Format as htsjdk VCFCodec(parse) + VCFWriter(encode) would: sorted INFO keys, GT first then sorted keys."""
        infos = ";".join("%s=%s" % (k, info[k]) for k in sorted(info)) if info else "."
        keys = sorted(k for k in fmt if k != "GT")
        allkeys = ["GT"] + keys
        cols = []
        for i in range(self.n):
            gt = fmt.get("GT", [None] * self.n)[i]
            if gt is None:
                gts = "."
            else:
                ph = any(p for a, p in gt)
                gts = ("|" if ph else "/").join("." if a is None else str(a) for a, p in gt)
            vals = []
            for k in keys:
                v = fmt[k][i]
                if v is None:
                    vals.append(".")
                elif isinstance(v, list):
                    vals.append(",".join("." if x is None else str(x) for x in v))
                else:
                    vals.append(str(v))
            while vals and all(ch in ".," for ch in vals[-1]):
                vals.pop()
            cols.append(":".join([gts] + vals))
        return "\t".join([self.contig, str(s), ".", ref, ",".join(alts) if alts else ".", ".", ".", infos,
                          ":".join(allkeys)] + cols)


def remap_alleles(data, inv, nall, input_nr, alt_only):
    out = []
    rng = range(1, nall) if alt_only else range(nall)
    for j in rng:
        ij = inv.get(j)
        if ij is None:
            if input_nr is None:
                out.append(None)
                continue
            ij = input_nr
        ix = ij - 1 if alt_only else ij
        out.append(data[ix] if ix < len(data) else None)
    return out


def remap_genotypes(data, inv, nall, input_nr, ploidy):
    if ploidy == 0:
        return None
    if ploidy == 1:
        out = []
        for j in range(nall):
            ij = inv.get(j, input_nr)
            out.append(None if ij is None or ij >= len(data) else data[ij])
        return out
    if ploidy == 2:
        out = [None] * (nall * (nall + 1) // 2)
        for j in range(nall):
            ij = inv.get(j, input_nr)
            for k in range(j, nall):
                gi = k * (k + 1) // 2 + j
                if ij is None:
                    out[gi] = None
                    continue
                ik = inv.get(k, input_nr)
                if ik is None:
                    out[gi] = None
                    continue
                a, b = sorted((ij, ik))
                igi = b * (b + 1) // 2 + a
                out[gi] = data[igi] if igi < len(data) else None
        return out
    raise NotImplementedError("ploidy %d" % ploidy)


def remap_2d(data, inv, nall, input_nr, alt_only):
    out = []
    rng = range(1, nall) if alt_only else range(nall)
    for j in rng:
        ij = inv.get(j)
        if ij is None:
            if input_nr is None:
                out.append([])
                continue
            ij = input_nr
        ix = ij - 1 if alt_only else ij
        out.append(data[ix] if ix < len(data) else [])
    return out


# ----------------------------------------------------------------------------------------------- the scan
def scan_window(M, recs_by_row, cells, cell_begins, a, b, out):
    """VariantQueryProcessor::scan_and_operate for one query column interval [a, b]."""
    n = M.n
    calls = [None] * n
    # gt_get_column(a): per row, the first cell (begin cell or END copy) with coordinate >= a
    for r in range(n):
        recs = recs_by_row[r]
        begins = recs.begins
        i = bisect.bisect_left(begins, a)
        best = None  # (coord, kind, rec)  kind 0=begin cell, 1=END copy
        if i < len(recs.list):
            best = (recs.list[i].begin, 0, recs.list[i])
        j = i - 1
        k = 0
        while j >= 0 and k < 64:
            ce = recs.copy_end[j]
            if ce is not None and ce >= a:
                if best is None or ce < best[0]:
                    best = (ce, 1, recs.list[j])
            j -= 1
            k += 1
        if best is None:
            continue
        coord, kind, rc = best
        if kind == 1:   # END copy: the call spans [begin, truncated END copy column]
            c = Call(rc, M)
            c.end = coord
            calls[r] = c
        elif coord == a:
            calls[r] = Call(rc, M)
    cur = a if any(c is not None for c in calls) else -1
    ndel = sum(1 for c in calls if c is not None and c.deletion_or_mnv)

    def ranges(cur, nxt):
        nonlocal ndel
        while cur < nxt:
            act = [c for c in calls if c is not None]
            if not act:
                break
            top = min(c.end for c in act)
            m = top if top < nxt - 1 else nxt - 1
            if ndel:
                m = cur
            out.append(M.operate(cur, m, calls))
            for r in range(n):
                c = calls[r]
                if c is not None and c.end == m:
                    if c.deletion_or_mnv:
                        ndel -= 1
                    calls[r] = None
            cur = m + 1
        return cur

    i0 = bisect.bisect_right(cell_begins, a)
    for idx in range(i0, len(cells)):
        rc = cells[idx]
        cb = rc.begin
        if cb > b:
            break
        if cur < 0:
            cur = cb
        if cb != cur:
            cur = ranges(cur, cb)
            cur = cb
        r = rc.row
        old = calls[r]
        if old is not None and old.end >= cb:
            if not old.deletion_or_mnv and not (len(old.rec.ref) == 1 and old.rec.alts == [NON_REF]):
                raise RuntimeError("Unhandled overlapping variants at %d row %d" % (cb, r))
            if old.deletion_or_mnv:
                ndel -= 1
        c = Call(rc, M)
        calls[r] = c
        if c.deletion_or_mnv:
            ndel += 1
    if cur >= 0:
        ranges(cur, b + 1)


class RowRecs:
    """One sample's stored cells. The loader (load_operators.cc LoaderArrayWriter) keeps each record's begin cell with
    its original END, but when the next record of the same row begins at or before that END it moves the END copy
    cell to next.begin-1 (and drops it if that equals the record's begin)."""
    def __init__(self, lst):
        self.list = lst
        self.begins = [r.begin for r in lst]
        self.copy_end = []
        for k, r in enumerate(lst):
            ce = r.end
            if k + 1 < len(lst) and lst[k + 1].begin <= ce:
                ce = lst[k + 1].begin - 1
            self.copy_end.append(ce if ce > r.begin else None)


def build(ws, mapfile, interval, out_path):
    contig, rng = interval.split(":")
    S, E = (int(x) for x in rng.replace(",", "").split("-"))
    callset = json.load(open(os.path.join(HERE, ws, "callset.json")))["callsets"]
    order = sorted(callset, key=lambda c: int(c["row_idx"]))
    samples = [c["sample_name"] for c in order]
    paths = dict(l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, mapfile)) if l.strip())
    arrays = [d for d in os.listdir(os.path.join(HERE, ws)) if d.count("$") == 2 and d.startswith(contig + "$")]
    arr_lo, arr_hi = (int(x) for x in arrays[0].split("$")[1:])
    fasta = Fasta(REF_FA)
    M = Merger(ws, samples, fasta, contig)
    rows = []
    for r, s in enumerate(samples):
        recs = read_gvcf(paths[s], contig, arr_lo, arr_hi, r)
        rows.append(RowRecs(recs))
        print("read %s: %d records" % (s, len(recs)), file=sys.stderr)
    cells = sorted((rc for rr in rows for rc in rr.list), key=lambda x: (x.begin, x.row))
    cell_begins = [c.begin for c in cells]
    out = []
    # GATK VariantLocusWalker: shards of SHARD bases from the interval start; FeatureCache refills on a miss with
    # [shardStart, shardEnd + SHARD]; GenomicsDBFeatureIterator clips the query to the array bounds.
    a = S
    with open(out_path, "w") as fo:
        while a <= E:
            shard_end = min(a + SHARD - 1, E)
            qa, qb = max(a, arr_lo), min(shard_end + SHARD, arr_hi)
            out.clear()
            if qa <= qb:
                scan_window(M, rows, cells, cell_begins, qa, qb, out)
            for line in out:
                pos = int(line.split("\t", 2)[1])
                if pos <= E:
                    fo.write(line + "\n")
            # next cache miss: the shard starting 2*SHARD after this window's start
            a += 2 * SHARD
    return samples


# ----------------------------------------------------------------------------------------------- comparison
def records(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line[0] != "#":
                yield line.rstrip("\n")


def compare(proto, received, nshow=8):
    A = {}
    stats = Counter()
    shown = 0
    it_p = records(proto)
    it_r = records(received)
    P = {}
    # both sorted by POS (single contig)
    p = next(it_p, None)
    r = next(it_r, None)
    diffs = []
    while p is not None or r is not None:
        pp = int(p.split("\t", 2)[1]) if p is not None else None
        rp = int(r.split("\t", 2)[1]) if r is not None else None
        if r is None or (p is not None and pp < rp):
            stats["only_proto"] += 1
            if shown < nshow:
                diffs.append("ONLY_PROTO " + p[:300]); shown += 1
            p = next(it_p, None)
            continue
        if p is None or rp < pp:
            stats["only_received"] += 1
            if shown < nshow:
                diffs.append("ONLY_RECV  " + r[:300]); shown += 1
            r = next(it_r, None)
            continue
        stats["both"] += 1
        if p == r:
            stats["line_match"] += 1
        else:
            tp, tr = p.split("\t"), r.split("\t")
            names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            for k in range(max(len(tp), len(tr))):
                a_ = tp[k] if k < len(tp) else None
                b_ = tr[k] if k < len(tr) else None
                if a_ != b_:
                    nm = names[k] if k < 9 else "SAMPLE"
                    stats["diff_" + nm] += 1
                    if nm == "INFO":
                        ia = dict(x.split("=", 1) if "=" in x else (x, "") for x in a_.split(";"))
                        ib = dict(x.split("=", 1) if "=" in x else (x, "") for x in b_.split(";"))
                        for key in set(ia) | set(ib):
                            if ia.get(key) != ib.get(key):
                                stats["diff_INFO_" + key] += 1
            if shown < nshow:
                diffs.append("PROTO " + p[:600] + "\nRECV  " + r[:600])
                shown += 1
        p = next(it_p, None)
        r = next(it_r, None)
    return stats, diffs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ws")
    ap.add_argument("--map")
    ap.add_argument("--interval")
    ap.add_argument("--out")
    ap.add_argument("--no-compare", action="store_true")
    ap.add_argument("--show", type=int, default=8)
    ap.add_argument("--queries", help="file of 'begin end' lines: emulate GenomicsDBFeatureReader.query() for each")
    ap.add_argument("--queries-received", help="harness/GdbQuery output for the same queries, to compare")
    a = ap.parse_args()
    mp, iv = WS_DEFAULTS.get(a.ws, (None, None))
    mp = a.map or mp
    iv = a.interval or iv
    if a.queries:
        out = a.out or os.path.join(HERE, "proto", a.ws + ".queries.vcf")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        run_queries(a.ws, mp, iv.split(":")[0], a.queries, out)
        if a.queries_received:
            compare_queries(out, a.queries_received, a.show)
        return
    out = a.out or os.path.join(HERE, "proto", a.ws + ".vcf")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    build(a.ws, mp, iv, out)
    if not a.no_compare:
        stats, diffs = compare(out, os.path.join(HERE, "received", a.ws + ".vcf"), a.show)
        for d in diffs:
            print(d)
        both = stats["both"]
        print("%s: received=%d proto=%d common=%d whole-line match=%d (%.4f%%)" % (
            a.ws, both + stats["only_received"], both + stats["only_proto"], both, stats["line_match"],
            100.0 * stats["line_match"] / max(1, both + stats["only_received"])))
        for k in sorted(stats):
            if k.startswith("diff_") or k.startswith("only"):
                print("  %s %d" % (k, stats[k]))




# ----------------------------------------------------------------------------------------------- arbitrary queries
def run_queries(ws, mapfile, contig, qfile, out_path):
    """Emulate GenomicsDBFeatureReader.query(contig, b, e) for each 'b e' line of qfile (see harness/GdbQuery.java)."""
    callset = json.load(open(os.path.join(HERE, ws, "callset.json")))["callsets"]
    samples = [c["sample_name"] for c in sorted(callset, key=lambda c: int(c["row_idx"]))]
    paths = dict(l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, mapfile)) if l.strip())
    arrays = [d for d in os.listdir(os.path.join(HERE, ws)) if d.count("$") == 2 and d.startswith(contig + "$")]
    arr_lo, arr_hi = (int(x) for x in arrays[0].split("$")[1:])
    M = Merger(ws, samples, Fasta(REF_FA), contig)
    rows = [RowRecs(read_gvcf(paths[s], contig, arr_lo, arr_hi, r)) for r, s in enumerate(samples)]
    cells = sorted((rc for rr in rows for rc in rr.list), key=lambda x: (x.begin, x.row))
    cell_begins = [c.begin for c in cells]
    with open(out_path, "w") as fo:
        for line in open(qfile):
            if not line.strip():
                continue
            b, e = (int(x) for x in line.split())
            fo.write("#QUERY\t%d\t%d\n" % (b, e))
            out = []
            qa, qb = max(b, arr_lo), min(e, arr_hi)
            if qa <= qb:
                scan_window(M, rows, cells, cell_begins, qa, qb, out)
            for l in out:
                fo.write(l + "\n")


def compare_queries(proto, received, nshow=8):
    def blocks(path):
        cur, lines = None, []
        for l in open(path):
            l = l.rstrip("\n")
            if l.startswith("#QUERY"):
                if cur is not None:
                    yield cur, lines
                cur, lines = l, []
            elif not l.startswith("#"):
                lines.append(l)
        if cur is not None:
            yield cur, lines
    def norm(line):
        """GdbQuery writes undecoded genotypes verbatim (htsjdk lazy fast path): re-order as GT + sorted keys and
        drop trailing missing values, like a decoded write."""
        t = line.split("\t")
        if len(t) < 10:
            return line
        keys = t[8].split(":")
        order = ["GT"] + sorted(k for k in keys if k != "GT")
        present = [k for k in order if any(dict(zip(keys, c.split(":"))).get(k, ".").strip(".,") != "" for c in t[9:])]
        if "GT" not in present:
            present = ["GT"] + present
        cols = []
        for c in t[9:]:
            d = dict(zip(keys, c.split(":")))
            v = [d.get(k, ".") for k in present]
            while len(v) > 1 and v[-1].strip(".,") == "":
                v.pop()
            cols.append(":".join(v))
        return "\t".join(t[:8] + [":".join(present)] + cols)
    st = Counter()
    shown = 0
    for (qp, lp), (qr, lr) in zip(blocks(proto), blocks(received)):
        lr = [norm(x) for x in lr]
        assert qp == qr, (qp, qr)
        st["queries"] += 1
        if lp == lr:
            st["query_match"] += 1
        else:
            if shown < nshow:
                shown += 1
                print(qp)
                for x, y in zip(lp + [""] * (len(lr) - len(lp)), lr + [""] * (len(lp) - len(lr))):
                    if x != y:
                        print("  PROTO " + x[:500] + "\n  RECV  " + y[:500])
                        break
        st["records"] += len(lr)
        st["record_match"] += sum(1 for x, y in zip(lp, lr) if x == y) if len(lp) == len(lr) else 0
    print(dict(st))




if __name__ == "__main__":
    main()
