# GenomicsDB 1.5.5 workspace — on-disk format (storage layer)

Scope: everything a reader needs to enumerate the **cells** of a workspace written by GATK 4.7
`GenomicsDBImport` (GenomicsDB 1.5.5, bundled TileDB fork commit 86e7245) and decode every
attribute value.  Query/merge semantics (how GenotypeGVCFs combines cells into variants) are out
of scope.

Prototype: `tiledb_read.py` (same directory, Python 3 stdlib only).  `dump` prints every cell,
`info` prints schema/fragment metadata, `validate` re-derives every cell from the source GVCFs
(section 9 gives the results: 0 mismatches over 13.4 M cells and 532 M attribute values).

Source references are relative to `genomicsdb-src/` (`T/` = `dependencies/TileDB/core/`).

All integers are **little endian**.  `i32/i64/u64/f32` are the obvious fixed-width types.
`size_t` and `off_t` are 8 bytes (x86-64 Linux).

---------------------------------------------------------------------------------------------

## 1. Workspace layout

```
ws/
  __tiledb_workspace.tdb          empty file; marks a TileDB workspace
  callset.json                    sample -> row
  vidmap.json                     fields + contigs (column offsets)
  vcfheader.vcf                   merged VCF header of all inputs (text)
  <array>/                        one TileDB array per imported interval
    __array_schema.tdb            binary schema (section 2)
    .__consolidation_lock         empty file (created with the array; flock target)
    genomicsdb_meta_dir/
      genomicsdb_column_bounds.json
      genomicsdb_meta_<uuid>.json   (one or more)
    __<uuid><tid>_<ms>/           fragment directories (section 3)
```

### 1.1 Arrays and intervals
* GATK passes one *partition* per (merged) `-L` interval, with contig positions, and asks for
  `generate_array_name_from_partition_bounds`.  The Java importer names the array
  `String.format("%s$%d$%d", contig, start, end)` (1-based inclusive interval), e.g. `chr8$1$15000000`
  (`src/main/java/org/genomicsdb/Constants.java`, `GenomicsDBImporter.updateConfigPartitionsAndLbUb`).
  If the partition was given by TileDB columns instead (`--merge-contigs-into-num-partitions`),
  the array is named by its partition rank (`"0"`, `"1"`, ...); then use
  `genomicsdb_column_bounds.json` for its range.
* A directory is an array iff it contains `__array_schema.tdb`.  Several arrays in one workspace
  are completely independent (a GVCF record overlapping two intervals is stored in both).
* `genomicsdb_column_bounds.json`: `{"min_column": C0, "max_column": C1}` = the partition
  in global columns (`offset+start-1`, `offset+end-1`).  ws_rb: `1391350276 .. 1406350275`.
  NB: cells may lie outside these bounds (section 5.4).
* `genomicsdb_meta_<uuid>.json`: `{"lb_row_idx": L, "max_valid_row_idx_in_array": M}`.  One file is
  written per import batch that extends the bounds; readers take `min(L)` and `max(M)` over all
  files (`VariantArrayInfo::read_row_bounds_from_metadata`).  `--consolidate` replaces them with a
  single file.  ws_rb_b1 has three files (M = 0, 1, 2), ws_rb one (M = 2).

### 1.2 callset.json
```json
{"callsets": [{"sample_name": "HG002", "row_idx": "0", "idx_in_file": "0", "stream_name": "HG002_stream"}, ...]}
```
`row_idx` (a decimal **string**) is the TileDB row coordinate.  GATK assigns rows in sorted
sample-name order (Java `TreeMap`, i.e. UTF-16 code-unit order = byte order for ASCII): ws_rb has
HG002=0, HG002pcr=1, NA12878=2 although rb.map lists NA12878 second.  Incremental imports continue
at `max_valid_row_idx_in_array + 1`.

### 1.3 vidmap.json
```json
{"fields":  [{"name": "AD", "type": ["Integer"], "vcf_field_class": ["FORMAT"],
              "length": [{"variable_length_descriptor": "R"}]}, ...],
 "contigs": [{"name": "chr1", "length": "248956422", "tiledb_column_offset": "0"}, ...]}
```
* **Global column** of a 1-based position: `column = tiledb_column_offset(contig) + POS - 1`.
  Offsets are the cumulative contig lengths in header order (chr8 = 1391350276 in GRCh38).
  To map back: last contig with `offset <= column`, `POS = column - offset + 1`.
* Field entries: `type` list (one entry per tuple element; `Integer|int`→i32, `Float|float`→f32,
  `String|char`→char), `length` list (one descriptor per dimension: `"1"`, `"2"`.., `"A"`, `"R"`,
  `"G"`, `"var"`, `"PP"` (GT: phased ploidy)), optional `vcf_delimiter` (one per dimension, e.g.
  `["|", ","]`), `VCF_field_combine_operation`, `disable_remap_missing_with_non_ref`.
  `ID` has no `vcf_field_class`; `LowQual` is a `FILTER` field; `DP` is `["FORMAT","INFO"]`.
  GATK lists: `ID`, FILTER fields, FORMAT fields (sorted), INFO fields (sorted).
* **Global field index** (used by the FILTER attribute): walk `fields` in order, assigning
  consecutive indices; immediately after a field that is both INFO and FORMAT insert
  `<name>_FORMAT`; immediately after a field with k>1 tuple types insert
  `<name>_tuple_element_0..k-1` (`vid_mapper_pb.cc` + `VidMapper::flatten_field`).  ws_rb:
  `0 ID, 1 LowQual, 2 AD, 3 DP, 4 DP_FORMAT, 5 GQ, ... 23 AS_RAW_BaseQRankSum,
  24 .._tuple_element_0, 25 .._tuple_element_1, 26 AS_RAW_MQ, ...`.

---------------------------------------------------------------------------------------------

## 2. `__array_schema.tdb` (uncompressed binary; `T/src/array/array_schema.cc` serialize/deserialize)

```
u32  version_tag            0xFFFFFFFD  (= 0xFFFFFFFF - 2 -> schema version 2)
                            [if 0xFFFFFFFF - tag > 15 the file is a pre-version schema with no
                             tag and no workspace field; never produced by 1.5.5]
i32  workspace_len, bytes   always 0 / empty
i32  array_name_len, bytes  "<workspace path as given to GATK>/<array>", e.g. "ws_rb/chr8$1$15000000"
                            (informational only — use the directory name)
u8   dense                  0 (sparse)
i8   tile_order             0 (ROW_MAJOR; irrelevant, no tile extents)
i8   cell_order             1 (COL_MAJOR)
i64  capacity               1000 (cells per tile)
i32  attribute_num  N
     N x { i32 len, bytes name }
i32  dim_num  D = 2
     D x { i32 len, bytes name }        "samples", "position"
i32  domain_size = 2*D*8 = 32
     i64 x 4                            [0, INT64_MAX-1, 0, INT64_MAX-1]  (lo,hi per dim)
i32  tile_extents_size = 0             (no extents)
i8   type[N+1]                          attribute types, last = coordinates type (1 = INT64)
i32  cell_val_num[N]                    1, 2, 3, 4 ... or 0x7FFFFFFF (TILEDB_VAR_NUM)
u8   compression[N+1]                   last = coordinates
i8   compression_level[N+1]
u8   offsets_compression[N]
i8   offsets_compression_level[N]
(end of file — the reader asserts offset == file size)
```
Types: 0 INT32, 1 INT64, 2 FLOAT32, 3 FLOAT64, 4 CHAR, 5 INT8, 6 UINT8, 7 INT16, 8 UINT16,
9 UINT32, 10 UINT64.  Compression byte = compressor (low nibble) | pre-filter:
0 none, 1 GZIP, 2 ZSTD, 3 LZ4, 4..9 BLOSC*, 10 RLE; `0x10` DELTA_ENCODE, `0x20` BIT_SHUFFLE.

Values written by GenomicsDBImport (defaults; GATK exposes no option to change them):
* every attribute: compression `0x01` (GZIP), level `-1` (Z_DEFAULT_COMPRESSION);
* var-length attributes: offsets_compression `0x11` (GZIP + DELTA_ENCODE), level −1;
  fixed-length attributes: offsets_compression 0, level 0 (unused);
* coordinates: compression `0x11`, level −1.

ws_rb byte example (first 0x50 bytes):
```
fd ff ff ff | 00 00 00 00 | 15 00 00 00 "ws_rb/chr8$1$15000000" | 00 | 00 | 01 |
e8 03 00 00 00 00 00 00 | 34 00 00 00 (52 attrs) | 03 00 00 00 "END" | 03 00 00 00 "REF" |
03 00 00 00 "ALT" | 02 00 00 00 "ID" | 04 00 00 00 "QUAL" ...
```

### 2.1 Attribute list (built by `VidMapper::build_tiledb_array_schema`)
In this order:
1. `END` INT64 ×1
2. `REF` CHAR var, `ALT` CHAR var, `ID` CHAR var (only if vidmap has field `ID` — GATK always
   writes it), `QUAL` FLOAT32 ×1, `FILTER` INT32 var
3. INFO fields: every global field (1.3 order) with the INFO class, except `END` and except
   multi-tuple parents (their `_tuple_element_i` children are used instead)
4. FORMAT fields, same rule; a field that is also INFO appears under `<name>_FORMAT` (e.g. `DP_FORMAT`)

Attribute type/cell_val_num: fixed-length descriptor (`"1"`, `"2"`, `"3"`, `"4"`) → that many
elements of the field type; anything else (`A`,`R`,`G`,`var`,`PP`) → `TILEDB_VAR_NUM`.
**Multi-dimensional fields** (two `length` descriptors, e.g. `AS_RAW_MQ` = `R`,`var` with
delimiters `|`,`,`) are stored as **CHAR var** (an opaque byte blob, section 7.6) — also each
`_tuple_element_i`.  `GT` is INT32 var.  `END` is INT64.

ws_rb (52 attributes): `END REF ALT ID QUAL FILTER DP AC AF AN AS_BaseQRankSum AS_FS
AS_InbreedingCoeff AS_MQ AS_MQRankSum AS_QD AS_QUALapprox AS_RAW_BaseQRankSum_tuple_element_0
AS_RAW_BaseQRankSum_tuple_element_1 AS_RAW_MQ AS_RAW_MQRankSum_tuple_element_0/1
AS_RAW_ReadPosRankSum_tuple_element_0/1 AS_ReadPosRankSum AS_SB_TABLE AS_SOR AS_VarDP BaseQRankSum
ExcessHet FS InbreedingCoeff MQ MQRankSum MQ_DP QD QUALapprox RAW_GT_COUNT(i32×3)
RAW_MQandDP(i32×2) ReadPosRankSum SOR VarDP AD DP_FORMAT GQ GT MIN_DP PGT PID PL PS SB(i32×4)`.
ws_raw has 36 (different header → different vidmap; e.g. `MLEAC`, `MLEAF`, no `AS_QUALapprox`).
**The schema is per workspace — always drive decoding from the schema + vidmap.**

---------------------------------------------------------------------------------------------

## 3. Fragments

### 3.1 Naming, visibility, order
* Name: `"__" + uuid(36 chars, lowercase with '-') + pthread_id(decimal u64) + "_" + ms_since_epoch`
  (`Array::new_fragment_name`), e.g.
  `__af2b4013-a03d-43d6-a80a-d2ae8a63108c128536578537152_1794728249843`.
  While being written the directory is `.__...`; on finalize it is renamed to `__...` and the
  empty marker file `__tiledb_fragment.tdb` is created.
* **A sub-directory of the array is a live fragment iff it contains `__tiledb_fragment.tdb`**
  (`is_fragment`).  Consolidation deletes the marker of old fragments before deleting their
  directories, so a marker-less leftover directory must be ignored.
* Order: ascending by the integer after the first `_` following the leading `__` (the ms
  timestamp; uuids contain no `_`) (`StorageManager::sort_fragment_names`).  Later = newer.
* A fragment is sparse iff `__coords.tdb` exists (always, here).

### 3.2 Files in a fragment
For every attribute `A` (all of them, even if never written; empty files are 0 bytes):
* fixed-size `A`: `A.tdb` = concatenation of compressed tiles of values
  (cell_val_num × type size bytes per cell);
* var-size `A`: `A.tdb` = compressed tiles of **u64 offsets** (one per cell), and
  `A_var.tdb` = compressed tiles of the concatenated values.  Fixed attributes have no `_var` file.
* `__coords.tdb`: compressed tiles of `(i64 row, i64 column)` per cell.
* `__book_keeping.tdb.gz`: metadata (3.3).
* `__tiledb_fragment.tdb`: 0-byte marker.

Tiling: cells are stored in **column-major cell order = sorted by (column, row)**; tile *t*
holds cells `[t*1000, t*1000+1000)`; every tile has `capacity` = 1000 cells except the last,
which has `last_tile_cell_num` cells.  All attributes use the same tiling.

### 3.3 `__book_keeping.tdb.gz` (`T/src/fragment/book_keeping.cc`)
A standard **gzip** member (deflateInit2 with windowBits|16; Python `gzip`, Rust
`flate2::read::GzDecoder`).  Decompressed content:
```
u64 non_empty_domain_size        32
i64 x 4                          non-empty domain [row_lo,row_hi,col_lo,col_hi]
                                 (GenomicsDB writes with no subarray, so this is the full
                                  domain 0..INT64_MAX-1 — useless for pruning; use MBRs)
i64 mbr_num  = T (number of tiles)
T x i64[4]                       MBR per tile: row_lo, row_hi, col_lo, col_hi
i64 bounding_coords_num = T
T x i64[4]                       first cell (row,col), last cell (row,col) of each tile
for a in 0..N (N attributes + coordinates as index N):
    i64 k; k x i64                 tile_offsets[a]: byte offset of each compressed tile in a.tdb
                                   (for var attrs: of the offsets tile)
for a in 0..N-1:
    i64 k; k x i64                 tile_var_offsets[a]: byte offset in a_var.tdb (k = 0 for fixed)
for a in 0..N-1:
    i64 k; k x u64                 tile_var_sizes[a]: UNCOMPRESSED size of each var tile
i64 last_tile_cell_num           cells in the last tile (= capacity if the count is a multiple)
```
Coordinates attribute index = N (after all attributes).  ws_rb: T = 1360 tiles, 1 359 702 cells,
last tile 702 cells; MBR[0] = (0, 2, 1391350276, 1391426101), bounding[0] =
(0, 1391350276) .. (1, 1391426101).

### 3.4 Tile framing / codecs (`T/src/fragment/write_state.cc`, `read_state.cc`, `codec*.cc`)
* Compressed size of tile *t* = `offsets[t+1] - offsets[t]`, or `file_size - offsets[t]` for the
  last tile.  No per-tile header.
* Uncompressed size: fixed attr `cells(t) * cell_size`; offsets tile `cells(t) * 8`; coords
  `cells(t) * 16`; var tile `tile_var_sizes[a][t]`.
* GZIP (1) = a **zlib** stream (`deflateInit`, 2-byte header `78 9c`, adler32 trailer) —
  Python `zlib.decompress`, Rust `flate2::read::ZlibDecoder`.  No other codec is used.
* DELTA_ENCODE (0x10) is applied before compression; to decode, after inflating:
  `v[i] += v[i - stride]` for i ≥ stride, per tile (restarts at 0 each tile).
  Coordinates: element type i64, stride = dim_num = 2 (rows and columns delta'd separately).
  Var offsets: element type u64, stride 1.
* A var tile whose uncompressed size is 0 writes nothing (`tile_var_offsets[t] ==
  tile_var_offsets[t+1]`, size 0).  E.g. `AC_var.tdb` is 0 bytes in ws_rb (no record has AC).
* **Var offsets are absolute** positions in the fragment-wide uncompressed var stream (they keep
  growing across tiles and across write batches).  Within tile *t* with offsets `o[0..n)`:
  cell *i* = `var_tile[o[i]-o[0] .. (o[i+1] if i+1<n else o[0]+tile_var_sizes[t]) - o[0]]`.
  (`o[0]` of tile t equals the sum of `tile_var_sizes` of tiles < t.)

Worked bytes (ws_rb, tile 0):
```
__coords.tdb  tile0: 1091 bytes "78 9c cd 9b bd 6e ..." -> 16000 bytes
  raw (delta) i64: 0 1391350276 | 1 0 | 1 0 | ...        decoded: (0,1391350276) (1,1391350276)
  (2,1391350276) (0,1391410431) (0,1391410432) (0,1391410445) ...
END.tdb tile0 -> i64: 1391410431 1391410620 1391410827 1391350276 ...
ALT.tdb tile0 (offsets) raw u64 deltas 0 1 1 1 1 1 1 1 1 1 3 3 -> 0 1 2 3 4 5 6 7 8 9 12 15
ALT_var.tdb tile0 (1098 bytes): "&&&&&&&&&G|&G|&&&&&T|&T|&..."; tile1 offsets start at 1098.
```

### 3.5 Reading a fragment (pseudo-code)
```
for t in 0..T:
  n = (t == T-1) ? last_tile_cell_num : capacity
  coords = delta_decode_i64_stride2(inflate(coords.tdb[slice t]))
  for each attribute a:
     fixed: vals = inflate(a.tdb[slice t])                  # n * cell_size bytes
     var:   offs = delta_decode_u64(inflate(a.tdb[slice t])) # n u64
            data = tile_var_sizes[a][t] ? inflate(a_var.tdb[var slice t]) : ""
```
MBRs/bounding coords allow skipping tiles by column range without inflating anything.

---------------------------------------------------------------------------------------------

## 4. Several fragments; consolidation

* One fragment per import batch (and per incremental import).  With `--batch-size B` and S samples,
  batch *j* contains rows `[jB, jB+B)`; each fragment holds the complete cells of its rows
  (ws_rb_b1: 3 fragments, rows 0, 1, 2 respectively, 425 174 / 514 361 / 420 167 cells).
* Read = **union** of all live fragments, merged in (column, row) order
  (`ArrayReadState::sort_fragment_cell_ranges`).  If two fragments contain the **same (row,
  column)**, only the cell of the **newest** fragment (larger timestamp → later in the sorted list)
  is returned (`SmallerPQFragmentCellRange`: "latest fragment wins"; older ranges are trimmed).
  There is no tombstone/deletion marker in this TileDB version.  GATK batches use disjoint rows,
  so duplicates do not occur in normal imports.
* Within one fragment duplicates cannot occur (the loader rejects overlapping begins for a row).
* `--consolidate` (`tiledb_array_consolidate`, `T/src/array/array.cc`): if the array has ≤ 1
  fragment it is a no-op (ws_rb, ws_raw, ws11: one batch → one fragment, untouched).  Otherwise
  all fragments are read through the merge above and re-written as one new fragment (normal
  sorted write, same capacity/codecs); old markers then directories are deleted.  The meta json
  files are merged into one.  Verified: a scratch workspace imported with `--batch-size 1
  --consolidate` produced a single fragment whose every `.tdb` file and decompressed book-keeping
  is **byte-identical** to ws_rb's (single-batch) fragment.

---------------------------------------------------------------------------------------------

## 5. GVCF record → cell: overview

GATK streams each sample's `VariantContext`s (htsjdk → BCF2 in memory) to GenomicsDB; the C++
loader (`src/main/cpp/src/vcf/vcf2binary.cc`) converts each BCF record into one *cell*
(`convert_VCF_to_binary_for_callset`), and `LoaderArrayWriter` (`src/loader/load_operators.cc`)
writes cells in (column,row) order, adding END-copy cells.  Every record of every sample that
GATK's reader returns for the interval becomes exactly one **begin cell**; records spanning more
than one position additionally produce at most one **END-copy cell**.  Nothing else is written.

### 5.1 Begin cell
* `row` = callset row of the sample; `column` = `offset + POS - 1`.
* `END` (i64, **0-based global column**, not a VCF position):
  * INFO `END` present: `offset + END - 1`;
  * otherwise, *deletion as interval* (`treat_deletions_as_intervals`, proto default **true**,
    GATK does not change it): scan ALT alleles 1..n in order; the first allele whose htslib
    variant type is INDEL or spanning deletion (`*`) **and** `strlen(REF) > strlen(ALT)` sets
    `END = column + len(REF) - 1`;
  * otherwise `END = column`.
  htslib's type: symbolic `<...>` alleles are never INDEL; after stripping the common prefix,
  a pure deletion/insertion is INDEL; e.g. `GAAAA→G` is INDEL, so `chr8:131675 GAAAA G,<NON_REF>`
  gets END = column+4.  (All 22 102 deletion records in the test GVCFs agree with the simple
  rule "some non-symbolic ALT shorter than REF"; see section 10 for the untested corner cases.)
* All other attributes: section 7.

### 5.2 END-copy cell
If `END > column`, a byte-identical copy of the cell is written **at column = END** with the
`END` attribute set to the **begin column** (so `END < column` identifies an END copy;
`END == column` is a single-position record; `END > column` a begin cell of an interval).
All other attributes are identical to the begin cell.

ws_rb, rows 0..2 at chr8:1 (ref blocks), dump excerpt:
```
row col         END                          REF ALT
0   1391350276  1391410431 (chr8:60156)      N   &       begin, END=60156
1   1391350276  1391410620 (chr8:60345)      N   &
2   1391350276  1391410827 (chr8:60552)      N   &
0   1391410431  1391350276 (chr8:1) [copy]   N   &       END copy of row 0 at 60156
0   1391410432  1391410445 (chr8:60170)      A   &
```

### 5.3 Overlap truncation of END copies (only deletions overlap in practice)
For each row the writer remembers the last END.  When the next record of the same row begins
at column `b <= pending END copy column` (e.g. records inside a deletion's span), the pending END
copy is moved to column `b - 1` (its `END` attribute stays = original begin) and written only if
`b - 1 != original begin`.  The begin cell is **not** modified (keeps the full END).
```
pending = None                                  # per row
for rec in records of the row (file order):
    b = begin column, e = END column (5.1)
    if pending:
        if pending.col >= b:  pending.col = b-1; emit copy if pending.col != pending.begin
        else:                 emit copy at pending.col
        pending = None
    emit begin cell (b, END=e)
    if e > b: pending = (col=e, begin=b)
emit pending if any
```
Example ws_raw row 0: `131675 GAAAA→G` (END 131679) followed by `131677 AAAAG→A,*,...`:
begin cell at 131675 (END=131679), copy truncated to **131676** (END attr 131675), begin at
131677 (END 131681), its copy at 131681, then the ref block at 131682.  If the next record had
started at 131676 the copy would have been dropped entirely.
Counts: ws_rb 540, ws_raw 1 597, ws11 612 truncated copies.
(Two records of one row with the same begin, or a begin before a previous begin, make the
loader throw — they never reach disk.)

### 5.4 Interval boundaries
* GATK feeds every record **overlapping** the interval.  A record beginning inside but ending
  after the interval end is stored unchanged, including its END copy **beyond `max_column`**
  (ws11: `chr8:4999977 END=5000244` → copy at column of 5000244 in array `chr8$1$5000000`;
  last MBR col_hi = chr8:5000255).
* From code (not exercised — all test intervals start at 1): a record beginning **before** the
  interval start but overlapping it is kept (the latest such record per row) and written at
  its original begin column, i.e. **before `min_column`** (`handle_intervals_spanning_partition_begin`).
* Readers must therefore not assume `min_column <= column <= max_column`.

---------------------------------------------------------------------------------------------

## 6. Coordinates summary
* dims: `samples` (row, i64), `position` (column, i64); domain `[0, 2^63-2]` both.
* Cell order inside a fragment and in the merged view: column ascending, then row ascending.

---------------------------------------------------------------------------------------------

## 7. Attribute value encodings

Values below are what lands in the attribute buffers (`BufferVariantCell::set_cell` strips the
loader's per-field length prefix; var attributes store only the elements, element count =
`bytes / element_size`).  "Missing" means: field absent from the record (or FORMAT key absent
for the sample), or the BCF value is a single missing value (VCF `.`).

### 7.1 Sentinels
| meaning | int32 | float32 bits | where |
|---|---|---|---|
| TileDB null (whole fixed field missing) | `0x7FFFFFFF` (INT32_MAX) | `0x7F7FFFFF` (FLT_MAX) | fixed-length attrs, QUAL |
| BCF missing element (`.` inside a vector) | `0x80000000` | `0x7F800001` | elements of vectors, multi-D |
| BCF vector end | `0x80000001` | `0x7F800002` | not stored (loader stops at it) |
| missing var-length field | — empty (0 bytes) — | | var attrs |
NB the float sentinels are NaN patterns (0x7F800001 is a signalling NaN): compare bits, and do
not round-trip them through a host float type that may quieten them.

### 7.2 Mandatory attributes
* `REF`: the REF string bytes (no terminator).
* `ALT`: ALT alleles joined with `|`; `<NON_REF>` is replaced by `&`
  (`TILEDB_NON_REF_VARIANT_REPRESENTATION`).  Other alleles (incl. `*`) verbatim.  ALT `.` →
  empty.  e.g. `G,<NON_REF>` → `47 7c 26` = `G|&`; `A,*,AAAAGAAAG,<NON_REF>` → `A|*|AAAAGAAAG|&`.
* `ID`: the ID string if it is not empty and not `.`; else empty.
* `QUAL`: f32; missing → FLT_MAX (`ff ff 7f 7f`).  e.g. 101.13 → `8f 42 ca 42`.
* `FILTER`: i32 var — one **global field index** (1.3) per filter in the record; `.` → empty.
  `PASS` is not in the vidmap, so it would be stored as `-1` (from code; not in test data: every
  test record has FILTER `.`).  `LowQual` → 1 in these vidmaps.

### 7.3 Ordinary INFO / FORMAT fields (1-D)
Element type from vidmap (`Integer`→i32, `Float`→f32, `String`→char).
* fixed length k (cell_val_num = k): if missing → k × TileDB null; otherwise the k values as
  stored in BCF (an inner `.` → BCF missing `0x80000000` / `0x7F800001`).
* var length (`A`, `R`, `G`, `var`): missing → empty; otherwise all elements (inner `.` → BCF
  missing).  Strings (`PGT`, `PID`): the bytes up to the first NUL, no terminator.
* Values are the BCF values htsjdk encoded: floats are `(float)Double.parseDouble(text)`
  (Python: `struct.pack('<f', float(text))`), ints are exact.
* FORMAT `DP` → attribute `DP_FORMAT`; INFO `DP` → `DP`.
* Flag fields (none in GATK GVCFs): char, value 1 when present (from code, untested).

### 7.4 GT (`length "PP"` → phased ploidy)
i32 var, `2*ploidy - 1` elements: `allele0, phase1, allele1, phase2, allele2, ...` where
allele = `bcf_gt_allele` = VCF allele index, **−1 for `.`**, and phase_k = 1 if the separator
before allele k is `|` else 0.  `0/0` → `00000000 00000000 00000000` = [0,0,0];
`1|1` → [1,1,1]; `0/1` → [0,0,1]; `./.` → [-1,0,-1]; `.|.` → [-1,1,-1] (both occur in ws_raw
and were verified); haploid `1` → [1] (from code).  GT absent → empty.

### 7.5 INFO "sum" fields with multi-sample inputs
`divide_up_among_samples` only applies when one input VCF has > 1 sample; GATK streams one
sample per input, so values are stored verbatim.

### 7.6 Multi-dimensional fields (`length` has 2 descriptors; stored as CHAR var)
GATK declares `AS_QUALapprox`, `AS_VarDP`, `AS_SB_TABLE` (int), `AS_RAW_MQ` (float) and
`AS_RAW_BaseQRankSum`, `AS_RAW_MQRankSum`, `AS_RAW_ReadPosRankSum` (tuple `[float, int]`) as
`R`,`var` with delimiters `|` (outer) and `,` (inner).  The VCF string is parsed by
`GenomicsDBMultiDVectorField::parse_and_store_numeric`:
* split on `|` into outer entries (one per allele incl. REF);
* split each entry on `,`; token *j* belongs to tuple element `j % k` (k = number of types;
  k = 1 for single-type fields);
* token → value: empty or `NaN` (case-insensitive) → BCF missing (`0x80000000` /
  `0x7F800001`); int: `strtoll(token, base 0)` truncated to i32 (so `0x..` hex, leading-0 octal
  — never occurs); float: `strtof`;
* within one outer entry, all tuple elements are padded with BCF missing to the same length
  (max over elements);
* each tuple element is serialized independently and stored in its own attribute
  (`<field>` for k = 1, `<field>_tuple_element_<i>` for k > 1):
```
u64 S                       total bytes of all inner vectors
S bytes                     inner vectors concatenated (i32 or f32 elements)
u64 N                       number of outer entries
u64 offsets[N+1]            byte offsets of each inner vector relative to the data start
                            (offsets[0] = 0, offsets[N] = S)
total = 16 + S + 8*(N+1)
```
* missing/absent field (or `.`) → empty attribute.
Examples (ws_rb row 0, chr8:60346; ws_raw row 0, chr8:131675):
```
AS_QUALapprox "|115|0" -> 0c00000000000000 | 00000080 73000000 00000000 | 0300000000000000 |
                          0000000000000000 0400000000000000 0800000000000000 0c00000000000000
                          (= [.] [115] [0])
AS_RAW_MQ "0.00|2322.00|0.00" -> S=12, data 00000000 00201145 00000000, N=3, offs 0,4,8,12
AS_RAW_BaseQRankSum "||" -> element_0: S=12 data 0100807f x3 (float missing), N=3
                            element_1: S=12 data 00000080 x3 (int missing; padding), N=3
AS_RAW_BaseQRankSum "|-1.1,1|NaN" -> element_0 [.] [-1.1] [.]   (01 00 80 7f | cd cc 8c bf | 01 00 80 7f)
                                     element_1 [.] [1]    [.]   (00 00 00 80 | 01 00 00 00 | 00 00 00 80)
AS_SB_TABLE "0,0|4,0|0,0" -> S=24, data 0,0,4,0,0,0 (i32), N=3, offs 0,8,16,24
```

### 7.7 Complete worked cell (ws_rb, row 0, column 1391410621 = chr8:60346, tile 0 cell 9)
GVCF: `A  G,<NON_REF>  101.13  .  AS_QUALapprox=|115|0;...;DP=4;QUALapprox=115;RAW_GT_COUNT=0,0,1;
RAW_MQandDP=2322,4;VarDP=4  GT:AD:DP:GQ:PL:SB  1/1:0,4,0:4:12:115,12,0,115,12,115:0,0,4,0`
```
END        bd 3d ef 52 00 00 00 00      = 1391410621 (single position)
REF        41                          "A"
ALT        47 7c 26                    "G|&"
ID         (empty)                     QUAL 8f 42 ca 42 = 101.13     FILTER (empty)
DP         04 00 00 00                 MQ ff ff 7f 7f (null)          PS/MIN_DP ff ff ff 7f (null)
RAW_MQandDP 12 09 00 00 04 00 00 00    = 2322,4
AD         00000000 04000000 00000000  GT 01000000 00000000 01000000 = [1,0,1] (1/1)
PL         115,12,0,115,12,115         SB 0,0,4,0        PGT (empty)
```

---------------------------------------------------------------------------------------------

## 8. Minimal reader algorithm (for the Rust port)
1. Parse `callset.json` (row→sample), `vidmap.json` (contigs→offset; fields→types/length/global
   index), list arrays (dirs with `__array_schema.tdb`).
2. Per array: parse schema (section 2); list live fragments (marker file), sort by timestamp;
   gunzip + parse each book-keeping file.
3. For a column range [c0,c1]: per fragment, select tiles with `mbr.col_hi >= c0 && mbr.col_lo
   <= c1` (remember END copies may be needed for intervals that begin before c0 — query
   semantics belong to the other spec); inflate coords and the needed attributes; filter cells.
4. k-way merge fragment streams by (column, row); on equal (row, column) keep the newest fragment.
5. Decode attributes per section 7.  Classify each cell by `END` vs column (begin / single /
   END copy).

Rust crates: `flate2` (zlib + gzip), `serde_json`.  Python: `zlib`, `gzip`, `json`, `struct`
(no third-party codec needed).

---------------------------------------------------------------------------------------------

## 9. Validation (tiledb_read.py validate)

For each workspace the prototype reads every cell (multi-fragment merge included), rebuilds the
expected cell set from the sample GVCFs with sections 5–7 (begin cells, END copies incl.
truncation, deletion intervals), and compares the coordinates, `END`, and the exact bytes of every
other attribute.

| workspace | fragments | cells | begin | END copies | attr values compared | mismatches | missing | extra |
|---|---|---|---|---|---|---|---|---|
| ws_rb    | 1 | 1 359 702 |   916 176 |   443 526 |  69 344 802 | 0 | 0 | 0 |
| ws_rb_b1 | 3 | 1 359 702 |   916 176 |   443 526 |  69 344 802 | 0 | 0 | 0 |
| ws11     | 1 | 1 814 409 | 1 240 357 |   574 052 |  92 534 859 | 0 | 0 | 0 |
| ws_raw   | 1 | 8 854 450 | 5 720 915 | 3 133 535 | 301 051 300 | 0 | 0 | 0 |
| scratch `--batch-size 1 --consolidate` | 1 | 1 359 702 | same as ws_rb, files byte-identical to ws_rb | | | 0 | 0 | 0 |

Deletion intervals (no INFO END, END from REF length): ws_rb 10 226, ws_raw 11 876, ws11 11 089.
Runtime: ~55 s (ws_rb) … 216 s (ws_raw) in CPython.

---------------------------------------------------------------------------------------------

## 10. Uncertainties / paths not exercised by the test data
1. **htslib variant-type corner cases** for END of deletions: GenomicsDB uses its htslib fork
   (`bcf_get_variant_type(...) == VCF_INDEL || == VCF_SPANNING_DELETION`).  The submodule source
   is not in the local checkout; `tiledb_read.py:htslib_var_type` re-implements upstream htslib's
   classification.  In the test data every record with a non-symbolic ALT shorter than REF is a
   plain INDEL, so complex alleles (e.g. `ATG→CT`, classified OTHER by htslib → END = POS) and a
   `*`-only allele with a long REF are not verified.
2. Records starting before the interval start (5.4) — from code only.
3. Duplicate (row, column) across fragments (newest wins) — from code only (GATK never does this).
4. FILTER `PASS` → −1 and non-`.` filters — from code only (all test FILTERs are `.`).
5. Flag fields, multi-sample input VCFs, BCF vector-end padding — not produced by GATK GVCFs.
6. Non-default import options (`disable_delta_encode_*`, `enable_bit_shuffle_gt`,
   `enable_lz4_compression_gt`, other compression types) are not settable from GATK; the
   prototype raises on anything other than zlib(+delta).
