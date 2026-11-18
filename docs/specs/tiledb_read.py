#!/usr/bin/env python3
"""Stdlib-only reader for GenomicsDB 1.5.5 workspaces written by GATK GenomicsDBImport.

Prototype for genomicsdb_tiledb.md (same directory).  Only zlib is needed: every tile written by
GenomicsDBImport is compressed with TILEDB_GZIP (=1), which is a *zlib* stream
(deflateInit, 0x78 header), optionally preceded by a delta-encode filter (TILEDB_DELTA_ENCODE
= 0x10) for the coordinates and for var-length offsets.  The book-keeping file is a real
gzip (RFC 1952) stream.

Usage:
  tiledb_read.py dump WORKSPACE [--array NAME] [--limit N]
      one line per cell: row <TAB> column <TAB> contig:pos <TAB> attr=value ...
  tiledb_read.py validate WORKSPACE SAMPLE_MAP
      checks every cell against the sample GVCFs (see genomicsdb_tiledb.md section 8)
  tiledb_read.py info WORKSPACE
      prints schema / fragments / book-keeping summary
"""
import gzip
import json
import os
import struct
import sys
import zlib

# ----------------------------------------------------------------------------------------
# TileDB constants (dependencies/TileDB/core/include/c_api/tiledb_constants.h)
# ----------------------------------------------------------------------------------------
T_INT32, T_INT64, T_FLOAT32, T_FLOAT64, T_CHAR, T_INT8, T_UINT8, T_INT16, T_UINT16, T_UINT32, T_UINT64 = range(11)
TYPE_FMT = {T_INT32: 'i', T_INT64: 'q', T_FLOAT32: 'f', T_FLOAT64: 'd', T_CHAR: 'c', T_INT8: 'b',
            T_UINT8: 'B', T_INT16: 'h', T_UINT16: 'H', T_UINT32: 'I', T_UINT64: 'Q'}
TYPE_SIZE = {k: struct.calcsize(v) for k, v in TYPE_FMT.items()}
VAR_NUM = 0x7FFFFFFF            # TILEDB_VAR_NUM (INT_MAX) as cell_val_num
COMPRESS_MASK = 0x0F            # low nibble: compressor
PRE_COMPRESS_MASK = 0x30        # TILEDB_DELTA_ENCODE (1<<4), TILEDB_BIT_SHUFFLE (2<<4)
NO_COMPRESSION, GZIP = 0, 1
DELTA_ENCODE, BIT_SHUFFLE = 0x10, 0x20

# null / missing sentinels (bit patterns)
TILEDB_NULL_INT32 = 0x7FFFFFFF      # TILEDB_EMPTY_INT32 = INT32_MAX
TILEDB_NULL_FLOAT_BITS = 0x7F7FFFFF  # TILEDB_EMPTY_FLOAT32 = FLT_MAX
BCF_INT32_MISSING = -0x80000000      # bcf_int32_missing
BCF_INT32_VECTOR_END = -0x7FFFFFFF   # bcf_int32_vector_end
BCF_FLOAT_MISSING_BITS = 0x7F800001
BCF_FLOAT_VECTOR_END_BITS = 0x7F800002


# ----------------------------------------------------------------------------------------
# Array schema  (__array_schema.tdb, uncompressed, little endian)
# ----------------------------------------------------------------------------------------
class Attribute:
    def __init__(self, idx, name):
        self.idx, self.name = idx, name
        self.type = self.cell_val_num = self.compression = self.level = None
        self.off_compression = self.off_level = None

    @property
    def var(self):
        return self.cell_val_num == VAR_NUM

    @property
    def type_size(self):
        return TYPE_SIZE[self.type]

    @property
    def cell_size(self):  # bytes per cell in <attr>.tdb
        return 8 if self.var else self.type_size * self.cell_val_num


class Schema:
    def __init__(self, path):
        b = open(path, 'rb').read()
        o = 0

        def rd(fmt):
            nonlocal o
            v = struct.unpack_from('<' + fmt, b, o)
            o += struct.calcsize('<' + fmt)
            return v if len(v) > 1 else v[0]

        def rstr():
            nonlocal o
            n = rd('i')
            s = b[o:o + n].decode()
            o += n
            return s

        tag = rd('I')
        # version tag present iff 0xFFFFFFFF - tag <= 0xF ; GenomicsDB 1.5.5 writes 0xFFFFFFFD (version 2)
        if 0xFFFFFFFF - tag <= 0xF:
            self.version = 0xFFFFFFFF - tag
            self.workspace = rstr()
        else:  # pre-versioned schema (no tag, no workspace field) - never written by GenomicsDB 1.5.5
            raise ValueError('unversioned array schema not supported by this prototype')
        self.array_name = rstr()
        self.dense = rd('B') != 0
        self.tile_order, self.cell_order = rd('b'), rd('b')
        self.capacity = rd('q')
        n = rd('i')
        self.attrs = [Attribute(i, rstr()) for i in range(n)]
        nd = rd('i')
        self.dims = [rstr() for _ in range(nd)]
        dsz = rd('i')
        self.domain_bytes = b[o:o + dsz]
        o += dsz
        tsz = rd('i')
        self.tile_extents_bytes = b[o:o + tsz]
        o += tsz
        types = [rd('b') for _ in range(n + 1)]
        cvn = [rd('i') for _ in range(n)]
        comp = [rd('B') for _ in range(n + 1)]
        lvl = [rd('b') for _ in range(n + 1)]
        ocomp = [rd('B') for _ in range(n)]
        olvl = [rd('b') for _ in range(n)]
        assert o == len(b), (o, len(b))
        for i, a in enumerate(self.attrs):
            a.type, a.cell_val_num, a.compression, a.level = types[i], cvn[i], comp[i], lvl[i]
            a.off_compression, a.off_level = ocomp[i], olvl[i]
        self.coords = Attribute(n, '__coords')
        self.coords.type, self.coords.cell_val_num = types[n], nd
        self.coords.compression, self.coords.level = comp[n], lvl[n]
        assert self.coords.type == T_INT64 and nd == 2 and not self.dense
        self.domain = struct.unpack('<%dq' % (dsz // 8), self.domain_bytes)
        self.by_name = {a.name: a for a in self.attrs}


# ----------------------------------------------------------------------------------------
# Tile codec
# ----------------------------------------------------------------------------------------
def delta_decode(raw, fmt, stride):
    vals = list(struct.unpack('<%d%s' % (len(raw) // struct.calcsize(fmt), fmt), raw))
    for j in range(stride, len(vals)):
        vals[j] += vals[j - stride]
    if fmt == 'Q':
        vals = [v & 0xFFFFFFFFFFFFFFFF for v in vals]
    return vals


def decompress(data, compression, expected_size):
    c = compression & COMPRESS_MASK
    if c == NO_COMPRESSION:
        out = data
    elif c == GZIP:
        out = zlib.decompress(data)   # zlib framing (deflateInit), not gzip framing
    else:
        raise NotImplementedError('compression %d (zstd=2, lz4=3, blosc*=4..9, rle=10) not needed by GATK' % c)
    assert len(out) == expected_size, (len(out), expected_size)
    if compression & PRE_COMPRESS_MASK == BIT_SHUFFLE:
        raise NotImplementedError('bit shuffle filter (only with enable_bit_shuffle_gt)')
    return out


# ----------------------------------------------------------------------------------------
# Fragment + book-keeping
# ----------------------------------------------------------------------------------------
class Fragment:
    def __init__(self, schema, path):
        self.schema, self.path = schema, path
        self.name = os.path.basename(path)
        # name = "__" + uuid(36 chars) + pthread_id(decimal) + "_" + ms_since_epoch
        self.timestamp = int(self.name[2:].split('_', 1)[1])
        n = len(schema.attrs)
        b = gzip.open(os.path.join(path, '__book_keeping.tdb.gz'), 'rb').read()
        o = 0

        def rd(fmt, cnt=1):
            nonlocal o
            v = struct.unpack_from('<%d%s' % (cnt, fmt), b, o)
            o += struct.calcsize('<%d%s' % (cnt, fmt))
            return v

        ned_size = rd('Q')[0]
        self.non_empty_domain = rd('q', ned_size // 8) if ned_size else None
        nm = rd('q')[0]
        flat = rd('q', 4 * nm)
        self.mbrs = [flat[4 * i:4 * i + 4] for i in range(nm)]           # row_lo,row_hi,col_lo,col_hi
        nb = rd('q')[0]
        flat = rd('q', 4 * nb)
        self.bounding = [flat[4 * i:4 * i + 4] for i in range(nb)]       # first(row,col), last(row,col)
        self.tile_offsets = []
        for _ in range(n + 1):
            k = rd('q')[0]
            self.tile_offsets.append(list(rd('q', k)) if k else [])
        self.tile_var_offsets = []
        for _ in range(n):
            k = rd('q')[0]
            self.tile_var_offsets.append(list(rd('q', k)) if k else [])
        self.tile_var_sizes = []
        for _ in range(n):
            k = rd('q')[0]
            self.tile_var_sizes.append(list(rd('Q', k)) if k else [])
        self.last_tile_cell_num = rd('q')[0]
        assert o == len(b), (o, len(b))
        self.tile_num = nm
        self.cell_num = (nm - 1) * schema.capacity + self.last_tile_cell_num if nm else 0
        self._files = {}

    def _file(self, fname):
        f = self._files.get(fname)
        if f is None:
            f = self._files[fname] = open(os.path.join(self.path, fname), 'rb').read()
        return f

    def tile_cells(self, t):
        return self.last_tile_cell_num if t == self.tile_num - 1 else self.schema.capacity

    @staticmethod
    def _chunk(data, offsets, t):
        end = offsets[t + 1] if t + 1 < len(offsets) else len(data)
        return data[offsets[t]:end]

    def coords(self, t):
        a = self.schema.coords
        raw = decompress(self._chunk(self._file('__coords.tdb'), self.tile_offsets[a.idx], t),
                         a.compression, self.tile_cells(t) * 16)
        if a.compression & PRE_COMPRESS_MASK == DELTA_ENCODE:
            v = delta_decode(raw, 'q', 2)                 # stride = dim_num, per-dimension deltas
        else:
            v = list(struct.unpack('<%dq' % (len(raw) // 8), raw))
        return [(v[2 * i], v[2 * i + 1]) for i in range(len(v) // 2)]

    def attr_tile(self, a, t):
        """Fixed attr: list of per-cell byte strings. Var attr: list of per-cell byte strings."""
        ncell = self.tile_cells(t)
        data = self._file(a.name + '.tdb')
        if not a.var:
            raw = decompress(self._chunk(data, self.tile_offsets[a.idx], t), a.compression, ncell * a.cell_size)
            cs = a.cell_size
            return [raw[i * cs:(i + 1) * cs] for i in range(ncell)]
        # offsets tile: uint64 absolute offsets into the fragment-wide var stream
        raw = decompress(self._chunk(data, self.tile_offsets[a.idx], t), a.off_compression, ncell * 8)
        if a.off_compression & PRE_COMPRESS_MASK == DELTA_ENCODE:
            offs = delta_decode(raw, 'Q', 1)
        else:
            offs = list(struct.unpack('<%dQ' % ncell, raw))
        vsize = self.tile_var_sizes[a.idx][t]
        if vsize:
            var = decompress(self._chunk(self._file(a.name + '_var.tdb'), self.tile_var_offsets[a.idx], t),
                             a.compression, vsize)
        else:
            var = b''
        base = offs[0]
        ends = offs[1:] + [base + vsize]
        return [var[s - base:e - base] for s, e in zip(offs, ends)]


def list_fragments(schema, array_dir):
    frs = []
    for d in os.listdir(array_dir):
        p = os.path.join(array_dir, d)
        # a directory is a (visible, committed) fragment iff it contains __tiledb_fragment.tdb
        if os.path.isdir(p) and os.path.isfile(os.path.join(p, '__tiledb_fragment.tdb')):
            frs.append(p)
    # TileDB sorts by the integer after the first '_' that follows the leading "__" (= ms timestamp)
    frs.sort(key=lambda p: int(os.path.basename(p)[2:].split('_', 1)[1]))
    return [Fragment(schema, p) for p in frs]


def iter_fragment_cells(fr):
    """Yields (row, col, {attr_name: raw_bytes}) in on-disk order (column-major: col, then row)."""
    attrs = fr.schema.attrs
    for t in range(fr.tile_num):
        co = fr.coords(t)
        cols = [fr.attr_tile(a, t) for a in attrs]
        for i, (r, c) in enumerate(co):
            yield r, c, {a.name: cols[k][i] for k, a in enumerate(attrs)}


def iter_array_cells(schema, frags):
    """Union of fragments in column-major order; for equal (row,col) the newest fragment wins."""
    if len(frags) == 1:
        yield from iter_fragment_cells(frags[0])
        return
    import heapq
    its = [iter_fragment_cells(f) for f in frags]
    heap = []
    for fid, it in enumerate(its):
        x = next(it, None)
        if x:
            heap.append(((x[1], x[0], -fid), x, fid))
    heapq.heapify(heap)
    last = None
    while heap:
        key, x, fid = heapq.heappop(heap)
        if last is not None and (key[0], key[1]) == last:
            pass                      # older duplicate - dropped
        else:
            last = (key[0], key[1])
            yield x
        y = next(its[fid], None)
        if y:
            heapq.heappush(heap, ((y[1], y[0], -fid), y, fid))


# ----------------------------------------------------------------------------------------
# Workspace + GenomicsDB semantics (vid map)
# ----------------------------------------------------------------------------------------
class Workspace:
    def __init__(self, ws):
        self.ws = ws
        self.vid = json.load(open(os.path.join(ws, 'vidmap.json')))
        cs = json.load(open(os.path.join(ws, 'callset.json')))
        self.row_to_sample = {int(c['row_idx']): c['sample_name'] for c in cs['callsets']}
        self.contigs = sorted(((int(c['tiledb_column_offset']), c['name'], int(c['length']))
                               for c in self.vid['contigs']))
        self.contig_offset = {n: o for o, n, _ in self.contigs}
        self.fields = self.vid['fields']
        self.field_by_name = {f['name']: f for f in self.fields}
        # GenomicsDB "global field index" (the FILTER attribute stores these).  VidMapper walks the
        # vidmap "fields" list in order; right after a field that is both INFO and FORMAT it inserts
        # "<name>_FORMAT", and after a field whose type has >1 tuple elements it inserts
        # "<name>_tuple_element_<i>" for each element (vid_mapper_pb.cc + VidMapper::flatten_field).
        self.global_fields = []
        for f in self.fields:
            self.global_fields.append(f['name'])
            cls = f.get('vcf_field_class', [])
            both = 'INFO' in cls and 'FORMAT' in cls
            if both:
                self.global_fields.append(f['name'] + '_FORMAT')
            nt = len(f.get('type', []))
            if nt > 1:
                for base in ([f['name'], f['name'] + '_FORMAT'] if both else [f['name']]):
                    self.global_fields += ['%s_tuple_element_%d' % (base, i) for i in range(nt)]
        self.field_idx = {n: i for i, n in enumerate(self.global_fields)}
        self.arrays = sorted(d for d in os.listdir(ws)
                             if os.path.isfile(os.path.join(ws, d, '__array_schema.tdb')))

    def locate(self, col):
        import bisect
        i = bisect.bisect_right([o for o, _, _ in self.contigs], col) - 1
        o, n, _ = self.contigs[i]
        return n, col - o + 1

    def open_array(self, name):
        d = os.path.join(self.ws, name)
        schema = Schema(os.path.join(d, '__array_schema.tdb'))
        return schema, list_fragments(schema, d)

    def field_for_attr(self, aname):
        """vidmap field, tuple element index (or None), is_format."""
        if aname.endswith('_FORMAT') and aname[:-7] in self.field_by_name:
            return self.field_by_name[aname[:-7]], None, True
        if '_tuple_element_' in aname:
            base, k = aname.rsplit('_tuple_element_', 1)
            return self.field_by_name[base], int(k), False
        f = self.field_by_name.get(aname)
        return f, None, (f is not None and 'FORMAT' in f.get('vcf_field_class', []) and
                         'INFO' not in f.get('vcf_field_class', []))


VID_TYPE = {'Integer': 'i', 'int': 'i', 'Float': 'f', 'float': 'f', 'String': 'c', 'char': 'c',
            'Flag': 'c', 'flag': 'c'}


def is_multid(field):
    return field is not None and len(field.get('length', [])) > 1


def decode_multid(raw, elem):
    """2-D GenomicsDB multi-D vector (e.g. AS_RAW_MQ 'R','var').  Returns list of lists.
    layout: u64 S | S bytes data (concatenated inner vectors) | u64 N | u64 offsets[N+1] (relative to data)
    float elements are returned as uint32 bit patterns."""
    if not raw:
        return None
    S = struct.unpack_from('<Q', raw, 0)[0]
    data = raw[8:8 + S]
    N = struct.unpack_from('<Q', raw, 8 + S)[0]
    offs = struct.unpack_from('<%dQ' % (N + 1), raw, 16 + S)
    assert len(raw) == 16 + S + 8 * (N + 1) and offs[-1] == S
    out = []
    for k in range(N):
        seg = data[offs[k]:offs[k + 1]]
        out.append(list(struct.unpack('<%d%s' % (len(seg) // 4, 'i' if elem == 'i' else 'I'), seg)))
    return out


def fmt_int(v):
    if v == TILEDB_NULL_INT32:
        return 'NULL'
    if v == BCF_INT32_MISSING:
        return '.'
    if v == BCF_INT32_VECTOR_END:
        return 'VEND'
    return str(v)


def f32_repr(bits):
    """float32 given as its uint32 bit pattern (Python floats would quieten the bcf sNaN sentinels)."""
    v = struct.unpack('<f', struct.pack('<I', bits))[0]
    if bits == TILEDB_NULL_FLOAT_BITS:
        return 'NULL'
    if bits == BCF_FLOAT_MISSING_BITS:
        return '.'
    if bits == BCF_FLOAT_VECTOR_END_BITS:
        return 'VEND'
    for p in range(1, 10):
        s = '%.*g' % (p, v)
        if struct.unpack('<f', struct.pack('<f', float(s)))[0] == v:
            return s
    return repr(v)


def format_value(ws, a, raw):
    """Human-readable decoding of one attribute value (see spec section 7)."""
    if a.type == T_CHAR:
        f, k, _ = ws.field_for_attr(a.name)
        if is_multid(f):
            elem = VID_TYPE[f['type'][k if k is not None else 0]]
            v = decode_multid(raw, elem)
            if v is None:
                return ''
            conv = fmt_int if elem == 'i' else f32_repr
            return '|'.join(','.join(conv(x) for x in vec) for vec in v)
        return raw.decode('latin-1')
    fmt = 'I' if a.type == T_FLOAT32 else TYPE_FMT[a.type]
    vals = struct.unpack('<%d%s' % (len(raw) // a.type_size, fmt), raw)
    if a.type == T_FLOAT32:
        return ','.join(f32_repr(x) for x in vals)
    if a.type == T_INT32:
        if a.name == 'FILTER':
            names = ws.global_fields
            return ','.join(names[x] if 0 <= x < len(names) else str(x) for x in vals)
        return ','.join(fmt_int(x) for x in vals)
    return ','.join(str(x) for x in vals)


def cmd_dump(ws_path, array=None, limit=None):
    ws = Workspace(ws_path)
    out = sys.stdout
    n = 0
    for arr in ws.arrays:
        if array and arr != array:
            continue
        schema, frags = ws.open_array(arr)
        for row, col, vals in iter_array_cells(schema, frags):
            ctg, pos = ws.locate(col)
            parts = ['%d\t%d\t%s:%d' % (row, col, ctg, pos)]
            for a in schema.attrs:
                raw = vals[a.name]
                if a.name == 'END':
                    e = struct.unpack('<q', raw)[0]
                    parts.append('END=%d(%s:%d)%s' % (e, ctg, ws.locate(e)[1], '' if e >= col else '[endcopy]'))
                elif a.var and not raw:
                    continue                       # empty var-length value = field absent
                else:
                    parts.append('%s=%s' % (a.name, format_value(ws, a, raw)))
            out.write('\t'.join(parts) + '\n')
            n += 1
            if limit and n >= limit:
                return


def cmd_info(ws_path):
    ws = Workspace(ws_path)
    print('arrays:', ws.arrays)
    for arr in ws.arrays:
        schema, frags = ws.open_array(arr)
        print('array %s: schema v%d name=%r capacity=%d cell_order=%d tile_order=%d dims=%s domain=%s' % (
            arr, schema.version, schema.array_name, schema.capacity, schema.cell_order, schema.tile_order,
            schema.dims, schema.domain))
        for a in schema.attrs + [schema.coords]:
            print('  attr %2d %-40s type=%d cvn=%s comp=0x%02x lvl=%d off_comp=%s' % (
                a.idx, a.name, a.type, 'var' if a.var else a.cell_val_num, a.compression, a.level,
                '0x%02x' % a.off_compression if a.off_compression is not None else '-'))
        for f in frags:
            print('  fragment %s ts=%d tiles=%d cells=%d last_tile=%d ned=%s' % (
                f.name, f.timestamp, f.tile_num, f.cell_num, f.last_tile_cell_num, f.non_empty_domain))
            print('    first MBR', f.mbrs[0], 'bounding', f.bounding[0], 'last MBR', f.mbrs[-1])


# ----------------------------------------------------------------------------------------
# Validation against the source GVCFs
# ----------------------------------------------------------------------------------------
def htslib_var_type(ref, alt):
    """Subset of htslib bcf_set_variant_type sufficient for GenomicsDB deletion handling.
    Returns 'NON_REF', 'OVERLAP' ('*'), 'SNP', 'MNP', 'INDEL', 'REF', 'OTHER'."""
    if alt == '*':
        return 'OVERLAP'
    if len(ref) == 1 and len(alt) == 1:
        if alt == '.' or ref.upper() == alt.upper() or alt == 'X':
            return 'REF'
        return 'SNP'
    if alt.startswith('<'):
        if alt in ('<X>',):
            return 'REF'
        if alt in ('<*>', '<NON_REF>'):
            return 'NON_REF'
        return 'OTHER'
    if alt[0] in '[]' or alt[-1] in '[]' or '[' in alt or ']' in alt:
        return 'OTHER'
    r = ref.upper()
    a = alt.upper()
    i = 0
    while i < len(r) and i < len(a) and r[i] == a[i]:
        i += 1
    if i < len(a) and i == len(r):
        return 'INDEL'
    if i < len(r) and i == len(a):
        return 'INDEL'
    if i == len(r) and i == len(a):
        return 'REF'
    re_, ae = len(r) - 1, len(a) - 1
    while re_ > i and ae > i and r[re_] == a[ae]:
        re_ -= 1
        ae -= 1
    if ae == i:
        if re_ == i:
            return 'SNP'
        return 'INDEL' if r[re_] == a[ae] else 'OTHER'
    if re_ == i:
        return 'INDEL' if r[re_] == a[ae] else 'OTHER'
    return 'MNP' if (re_ - i) == (ae - i) and (re_ - i) > 1 else 'OTHER'


def gvcf_end(pos, ref, alts, info):
    """1-based END GenomicsDB stores (treat_deletions_as_intervals=true, GATK default)."""
    if 'END' in info:
        return int(info['END'])
    for alt in alts:
        t = htslib_var_type(ref, alt)
        if t in ('INDEL', 'OVERLAP') and len(ref) > len(alt):
            return pos + len(ref) - 1
    return pos


def f32(x):
    return struct.unpack('<f', struct.pack('<f', float(x)))[0]


def pack_ints(vals):
    return struct.pack('<%di' % len(vals), *vals)


def pack_floats(vals):
    return b''.join(struct.pack('<I', BCF_FLOAT_MISSING_BITS) if v is None else struct.pack('<f', v)
                    for v in vals)


def parse_num(tok, elem):
    if tok == '.' or tok == '':
        return BCF_INT32_MISSING if elem == 'i' else None
    if elem == 'i':
        return int(tok)
    return f32(tok)


def strtol0(tok):
    """C strtoll(..., base 0) on the full token (as used by the multi-D parser)."""
    t = tok.strip()
    neg = t.startswith('-')
    s = t.lstrip('+-')
    try:
        if s.lower().startswith('0x'):
            v = int(s, 16)
        elif len(s) > 1 and s.startswith('0'):
            v = int(s, 8)
        else:
            v = int(s)
    except ValueError:
        v = 0
    return -v if neg else v


def expected_multid(s, field):
    """Bytes for each tuple element of a 2-D multi-D vector field given its VCF string."""
    types = [VID_TYPE[t] for t in field['type']]
    nt = len(types)
    d_outer, d_inner = field.get('vcf_delimiter', ['|', ','])[:2]
    outs = []
    groups = s.split(d_outer)
    per_elem = [[] for _ in range(nt)]          # per tuple element: list of inner vectors
    for g in groups:
        toks = g.split(d_inner)
        vecs = [[] for _ in range(nt)]
        for j, tok in enumerate(toks):
            k = j % nt
            if types[k] == 'i':
                v = BCF_INT32_MISSING if tok == '' or tok.lower() == 'nan' else strtol0(tok)
                v = ((v + 2 ** 31) % 2 ** 32) - 2 ** 31
            else:
                v = None if tok == '' or tok.lower() == 'nan' else f32(tok)
            vecs[k].append(v)
        m = max(len(v) for v in vecs)
        for k in range(nt):
            while len(vecs[k]) < m:
                vecs[k].append(BCF_INT32_MISSING if types[k] == 'i' else None)
            per_elem[k].append(vecs[k])
    for k in range(nt):
        datas = [pack_ints(v) if types[k] == 'i' else pack_floats(v) for v in per_elem[k]]
        offs = [0]
        for d in datas:
            offs.append(offs[-1] + len(d))
        data = b''.join(datas)
        outs.append(struct.pack('<Q', len(data)) + data + struct.pack('<Q', len(datas)) +
                    struct.pack('<%dQ' % len(offs), *offs))
    return outs


def expected_attrs(ws, schema, rec):
    """Expected raw bytes of every attribute (except END) for a parsed GVCF record."""
    chrom, pos, vid, ref, alts, qual, filt, info, fmt = rec
    exp = {}
    exp['REF'] = ref.encode()
    exp['ALT'] = '|'.join('&' if a == '<NON_REF>' else a for a in alts).encode()
    if 'ID' in schema.by_name:
        exp['ID'] = b'' if vid == '.' else vid.encode()
    exp['QUAL'] = struct.pack('<I', TILEDB_NULL_FLOAT_BITS) if qual == '.' else struct.pack('<f', f32(qual))
    flts = [] if filt == '.' else filt.split(';')
    exp['FILTER'] = pack_ints([ws.field_idx.get(x, -1) for x in flts])   # PASS (not in vidmap) -> -1
    for a in schema.attrs:
        if a.name in exp or a.name == 'END':
            continue
        f, k, is_fmt = ws.field_for_attr(a.name)
        src = fmt if is_fmt else info
        key = f['name']
        val = src.get(key)
        if is_multid(f):
            if val is None or val == '.':
                exp[a.name] = b''
            else:
                exp[a.name] = expected_multid(val, f)[k if k is not None else 0]
            continue
        if key == 'GT':
            if val is None or val == '.':
                # '.' alone in BCF is a single bcf_gt_missing (0) element, not bcf_int32_missing
                exp[a.name] = pack_ints([-1]) if val == '.' else b''
                continue
            import re
            al = re.split(r'[/|]', val)
            seps = re.findall(r'[/|]', val)
            out = [-1 if al[0] == '.' else int(al[0])]
            for s_, x in zip(seps, al[1:]):
                out += [1 if s_ == '|' else 0, -1 if x == '.' else int(x)]
            exp[a.name] = pack_ints(out)
            continue
        elem = VID_TYPE[f['type'][0]]
        if a.type == T_CHAR:
            exp[a.name] = b'' if val is None or val == '.' else val.encode()
            continue
        if val is True:  # flag
            val = '1'
        missing = val is None or val == '.'
        if a.var:
            if missing:
                exp[a.name] = b''
            else:
                vals = [parse_num(t, elem) for t in val.split(',')]
                exp[a.name] = pack_ints(vals) if elem == 'i' else pack_floats(vals)
        else:
            n = a.cell_val_num
            if missing:
                exp[a.name] = (struct.pack('<i', TILEDB_NULL_INT32) if elem == 'i'
                               else struct.pack('<I', TILEDB_NULL_FLOAT_BITS)) * n
            else:
                vals = [parse_num(t, elem) for t in val.split(',')]
                exp[a.name] = pack_ints(vals) if elem == 'i' else pack_floats(vals)
    return exp


def read_gvcf(path, contig, start, end):
    """Records overlapping [start,end] on contig (1-based, inclusive)."""
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            p = line.rstrip('\n').split('\t')
            if p[0] != contig:
                continue
            pos = int(p[1])
            if pos > end:
                continue
            info = {}
            if p[7] != '.':
                for kv in p[7].split(';'):
                    if '=' in kv:
                        k, v = kv.split('=', 1)
                        info[k] = v
                    else:
                        info[kv] = True
            fmt = dict(zip(p[8].split(':'), p[9].split(':'))) if len(p) > 9 else {}
            rec = (p[0], pos, p[2], p[3], p[4].split(',') if p[4] != '.' else [], p[5], p[6], info, fmt)
            e = gvcf_end(pos, p[3], rec[4], info)
            if e < start:
                continue
            yield rec, e


def cmd_validate(ws_path, sample_map):
    ws = Workspace(ws_path)
    smap = dict(l.rstrip('\n').split('\t') for l in open(sample_map) if l.strip())
    sample_to_row = {s: r for r, s in ws.row_to_sample.items()}
    total = {'cells': 0, 'begin_cells': 0, 'end_copy_cells': 0, 'records': 0, 'attr_values': 0,
             'mismatch': 0, 'missing_cells': 0, 'extra_cells': 0, 'truncated_end_copies': 0,
             'deletion_intervals': 0}
    for arr in ws.arrays:
        contig, s, e = arr.split('$')
        s, e = int(s), int(e)
        off = ws.contig_offset[contig]
        schema, frags = ws.open_array(arr)
        # expected cells: (row, col) -> (END value, record)
        expected = {}
        for sample, path in smap.items():
            row = sample_to_row[sample]
            pending = None          # [copy_col, orig_begin_col, rec]
            for rec, end1 in read_gvcf(path, contig, s, e):
                total['records'] += 1
                b = off + rec[1] - 1
                ec = off + end1 - 1
                if 'END' not in rec[7] and ec > b:
                    total['deletion_intervals'] += 1
                if pending is not None:
                    if pending[0] >= b:                  # overlap: END copy truncated to b-1
                        total['truncated_end_copies'] += 1
                        pending[0] = b - 1
                        if pending[0] != pending[1]:
                            expected[(row, pending[0])] = (pending[1], pending[2], 'endcopy')
                    else:
                        expected[(row, pending[0])] = (pending[1], pending[2], 'endcopy')
                    pending = None
                assert (row, b) not in expected, ('duplicate begin', sample, rec[1])
                expected[(row, b)] = (ec, rec, 'begin')
                if ec > b:
                    pending = [ec, b, rec]
            if pending is not None:
                expected[(row, pending[0])] = (pending[1], pending[2], 'endcopy')
        exp_cache = {}
        nmis_print = 0
        seen = set()
        for row, col, vals in iter_array_cells(schema, frags):
            total['cells'] += 1
            key = (row, col)
            seen.add(key)
            if key not in expected:
                total['extra_cells'] += 1
                if nmis_print < 20:
                    print('EXTRA cell', row, col, ws.locate(col), file=sys.stderr)
                    nmis_print += 1
                continue
            end_exp, rec, kind = expected[key]
            total['begin_cells' if kind == 'begin' else 'end_copy_cells'] += 1
            got_end = struct.unpack('<q', vals['END'])[0]
            bad = []
            if got_end != end_exp:
                bad.append(('END', got_end, end_exp))
            rid = id(rec)
            ex = exp_cache.get(rid)
            if ex is None:
                ex = expected_attrs(ws, schema, rec)
                if kind == 'begin' and end_exp > col:
                    exp_cache[rid] = ex         # reused by its END copy
            else:
                del exp_cache[rid]
            for name, want in ex.items():
                total['attr_values'] += 1
                if vals[name] != want:
                    bad.append((name, vals[name], want))
            if bad:
                total['mismatch'] += 1
                if nmis_print < 20:
                    print('MISMATCH', kind, row, col, ws.locate(col), rec[1:7], file=sys.stderr)
                    for b_ in bad:
                        print('   ', b_, file=sys.stderr)
                    nmis_print += 1
        for key in expected:
            if key not in seen:
                total['missing_cells'] += 1
                if nmis_print < 40:
                    print('MISSING cell', key, ws.locate(key[1]), expected[key][2], expected[key][1][1:6],
                          file=sys.stderr)
                    nmis_print += 1
        print('array %s: %d fragments' % (arr, len(frags)))
    print(' '.join('%s=%d' % kv for kv in total.items()))
    return total['mismatch'] == 0 and total['missing_cells'] == 0 and total['extra_cells'] == 0


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    cmd, ws = argv[1], argv[2]
    if cmd == 'dump':
        array = limit = None
        rest = argv[3:]
        while rest:
            if rest[0] == '--array':
                array, rest = rest[1], rest[2:]
            elif rest[0] == '--limit':
                limit, rest = int(rest[1]), rest[2:]
            else:
                raise SystemExit('unknown option ' + rest[0])
        try:
            cmd_dump(ws, array, limit)
        except BrokenPipeError:
            pass
        return 0
    if cmd == 'info':
        cmd_info(ws)
        return 0
    if cmd == 'validate':
        return 0 if cmd_validate(ws, argv[3]) else 1
    print(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main(sys.argv))
