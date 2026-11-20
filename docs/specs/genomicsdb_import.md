# GenomicsDBImport (GATK 4.7.0.0, GenomicsDB 1.5.5): workspace metadata, intervals, import flow

Scope: everything `GenomicsDBImport` writes or decides **besides** the TileDB storage format and the
record→cell derivation, which are in `genomicsdb_tiledb.md` (sections 1, 5 and 7 there are assumed).
Covered: `vidmap.json`, `callset.json`, `vcfheader.vcf`, the per-array files, how `-L` becomes arrays,
batching, what records reach the C++ loader, the two reader paths (default and
`--bypass-feature-reader`), and the checks the tool runs, with their error messages.

Source references:
* `G/` = `gatk-src/src/main/java/org/broadinstitute/hellbender/`
* `J/` = `genomicsdb-src/src/main/java/org/genomicsdb/`
* `C/` = `genomicsdb-src/src/main/cpp/`
* `T/` = `genomicsdb-src/dependencies/TileDB/core/`
* `H/` = `htsjdk-src/htsjdk/` (5.0.0)

Every claim carries one of these tags:
* **[V]**: checked against the GATK-written workspaces in `gatk/genomicsdb/`: ws_rb, ws_raw, ws_rb_b1,
  ws11, ws_wg and ws_mid (built by `tests/harness/genomicsdb/truth.sh` and `truth2.sh`), plus
  ws_rb_bypass and ws_raw_bypass. The last two are the rb.map and raw.map imports re-run with
  `--bypass-feature-reader --genomicsdb-shared-posixfs-optimizations`; their commands are in their
  `vcfheader.vcf`.
* **[S]**: taken from source only.

Byte-exact checks: a Python re-implementation of sections 2, 3 and 4 (algorithm in 2.7 and 4.1)
reproduces `vidmap.json`, `callset.json` and `vcfheader.vcf` **byte for byte** for ws_rb, ws_raw,
ws_rb_b1, ws11, ws_rb_bypass and ws_raw_bypass. For `vcfheader.vcf`, the `##GATKCommandLine` line is
copied from the target file. ws_wg and ws_mid cannot be rebuilt because their input GVCF is gone.

---------------------------------------------------------------------------------------------

## 1. What is written, and in which order [S, order of files V by timestamps]

`GenomicsDBImport.onTraversalStart` → `GenomicsDBImporter(ImportConfig)` → `executeImport`:

1. Workspace directory + empty `__tiledb_workspace.tdb`
   (`G/tools/genomicsdb/GenomicsDBImport.java:1038-1062` → `T/src/c_api/tiledb_utils.cc:103-136`).
2. `vidmap.json`, first version, built from the merged header (`J/importer/GenomicsDBImporter.java:254-259`).
3. `vcfheader.vcf` (`:262-264`).
4. `vidmap.json` again, after GATK's overrides (`G/.../GATKGenomicsDBUtils.java:54-103` →
   `GenomicsDBImporter.updateProtobufVidMapping:301-310`). Only this final version matters (section 2).
5. Only with `--merge-contigs-into-num-partitions`: `vidmap.json` a third time, with the contig list
   regenerated (2.6).
6. For each batch, then each array (partition):
   * create the array if it does not exist: directory, `__array_schema.tdb`, `.__consolidation_lock`;
   * `genomicsdb_meta_dir/genomicsdb_column_bounds.json` (first batch only);
   * `genomicsdb_meta_<uuid>.json` (when the row bounds grow);
   * one fragment;
   * on the last batch with `--consolidate`: consolidation (5.4).
7. `callset.json`, **last**, only after every batch succeeded (`GenomicsDBImporter.postImportFileCleanup:629-646`).
   The harness's `[ -f callset.json ]` success test relies on this.

All files and directories created by GenomicsDB/TileDB have mode **0700** (`-rwx------`, `drwx------`) [V].
The cause: `open(..., S_IRWXU)` and `mkdir(dir, S_IRWXU)` (`T/src/storage/storage_posixfs.cc:197,431,460`;
`T/src/storage/storage_manager.cc:489,1550,1707`). The JSON files and `vcfheader.vcf` are written through
`GenomicsDBUtils.writeToFile` → TileDB, so they are 0700 too [V].

---------------------------------------------------------------------------------------------

## 2. `vidmap.json`

### 2.1 Which header the field list comes from [S; V for the sample-name-map case]
`G/.../GenomicsDBImport.java:557-641`:

* **`--sample-name-map`** (WARP, all test workspaces): the header of the **first sample in sorted
  sample-name order** (`SampleNameMap` is a `TreeMap`; `:600-605`). Other inputs' headers are never
  merged. A field that exists only in another sample's header is not in the vidmap and is silently not
  stored.
  * Read by htsjdk's `VCFCodec`, which applies `VCFStandardHeaderLines.repairStandardHeaderLines`
    (`H/variant/vcf/VCFStandardHeaderLines.java:62-83`): a standard INFO/FORMAT line with the wrong
    Number or Type is replaced by the standard one. `variant-rs/src/vcfheader.rs:repair_standard_line`
    already implements this.
  * Line order = **input file order** (`getMetaDataInInputOrder()`, `H/variant/vcf/VCFHeader.java:475-494`:
    a `##fileformat` line, then the lines in file order, exact duplicates once).
  * Then GATK appends `##source=GenomicsDBImport` and the `##GATKCommandLine` line (4.2). They are not
    fields.
* **`-V` files**: `VCFUtils.smartMergeHeaders(headers, true)` (`H/variant/vcf/VCFUtils.java:48-127`).
  It walks the headers in `-V` order, each in **sorted** order (`getMetaDataInSortedOrder`), and keys
  each line by `key` or `key-ID`. The first occurrence of a key wins, with three exceptions:
  * compound lines whose Number differs → the kept line is changed **in place** to `Number=.`
    (vidmap length `"VAR"`, 2.2);
  * an Integer line meeting a Float line → the Float line is kept;
  * any other type conflict → `IllegalStateException("Incompatible header types, collision between these two types: ...")`.

  A second plain line with the same key (e.g. a second `##source=`) is dropped with a warning. Each
  `-V` header must have exactly one sample (2.6 / 6.4).
* **`--header <vcf>`**: that file's header, in input order (`:584-587`, `:606-608`).

The vidmap contigs and the sequence dictionary used to parse `-L` both come from this header
(`getBestAvailableSequenceDictionary` override, `:1127-1134`). No reference is needed.

### 2.2 Field entries (`J/importer/extensions/VidMapExtensions.java:49-184`) [V]
`fields[0]` is always `{"name":"ID","type":["char"],"length":[{"variable_length_descriptor":"var"}]}`
(no `vcf_field_class`). Then each header line, in the order of 2.1, is handled as follows. All other
lines (ALT, GVCFBlock, GATKCommandLine, source, fileformat, ...) are ignored.

| header line | name | type | length (`variable_length_descriptor`) | vcf_field_class |
|---|---|---|---|---|
| `##FILTER=<ID=X>` | X | `["Integer"]` | `"1"` | `["FILTER"]` |
| `##FORMAT=<ID=GT,...>` | GT | `["int"]` | `"PP"` | `["FORMAT"]` |
| `##FORMAT` Type=String | X | `["String"]` | `"var"` (whatever Number says) | `["FORMAT"]` |
| `##FORMAT` other types | X | `[Type]` | `L(Number)` | `["FORMAT"]` |
| `##INFO` in the AS set below | X | see below | `"R"`, `"var"` | `["INFO"]` |
| `##INFO` Type=String | X | `["String"]` | `"var"` | `["INFO"]` |
| `##INFO` other types | X | `[Type]` | `L(Number)` | `["INFO"]` |
| `##contig=<ID=c,length=n>` | goes to `contigs` (2.5) | | | |

* `Type` is htsjdk's `VCFHeaderLineType.toString()`: `Integer`, `Float`, `String`, `Flag`, `Character`.
* `L(Number)` (`getVcfHeaderLength:198-240`) maps `A`/`R`/`G` → `"A"`/`"R"`/`"G"`, `.` → **`"VAR"`**
  (upper case), and an integer n → `"n"`. Exception: `Number=0` with `Type=Flag` → `"1"`.
  * The C++ side upper-cases length descriptors (`C/src/utils/vid_mapper.cc:750-770`), so `"var"` and
    `"VAR"` mean the same. They are still different bytes in the file.
  * No test header has `Number=.` on a non-String field, so `"VAR"` is [S].
* **`Type=Character`** is written to the vidmap as `"Character"`. The C++ type table has no such name
  (`vid_mapper.cc:54-102`), so the import aborts with `VERIFY_OR_THROW ... "Unhandled field type"`
  (`vid_mapper.cc:782-787`) [S]. Flag → `"Flag"` (char, BCF_HT_FLAG) [S].
* Allele-specific INFO sets (`J/importer/Constants.java:13-25`), applied by **name**, whatever the
  header says:
  * `AS_RAW_BaseQRankSum`, `AS_RAW_MQRankSum`, `AS_RAW_ReadPosRankSum`: type `["float","int"]`,
    `VCF_field_combine_operation` `"histogram_sum"`;
  * `AS_RAW_MQ`: type `["float"]`, `"element_wise_sum"`;
  * `AS_SB_TABLE`, `AS_QUALapprox`, `AS_VarDP`: type `["int"]`, `"element_wise_sum"`.

  All seven get `length` `[{"R"},{"var"}]` and `vcf_delimiter` `["|", ","]`.
* **Same name in several classes** (`checkForDuplicateVcfFieldNames:242-300`). Example: FORMAT DP and
  INFO DP. When the second line with that name is met, the classes already seen for the name are
  checked in the fixed order FILTER, INFO, FORMAT, and the first one found decides:
  * same type string and same length string → **no new entry**; the class is appended to the existing
    entry's `vcf_field_class`. ws_rb: the FORMAT DP line comes first in the input, so DP is
    `["FORMAT","INFO"]` [V];
  * otherwise → a new entry named `<ID>_<CLASS>` (e.g. `AF_FORMAT`) with `"vcf_name": "<ID>"`.

  Comparison uses the strings computed above (FORMAT GT is `int`/`PP`; String is `var`). A FILTER line
  counts as `Integer`/`1`, so `##FILTER=<ID=X>` would merge with an `Integer, Number=1` INFO X.
  A second line with the same name **and** class is skipped, with a message on stderr.
* `##FILTER=<ID=PASS>` is not special: if the input declares it, `PASS` becomes a field [S]. No GATK
  GVCF declares it.

### 2.3 GATK's overrides (`G/.../GATKGenomicsDBUtils.java:54-103`) [V for the fields present]
Applied by **field name**, only if the field exists, in this order:

| field | change |
|---|---|
| `RAW_MQandDP` | `VCF_field_combine_operation` = `"element_wise_sum"` |
| `AS_RAW_MQ` | rebuilt as AS: type `["float"]`, `"element_wise_sum"` (same as 2.2) |
| `MQ_DP`, `QUALapprox`, `VarDP`, `TREE_SCORE`, `ASSEMBLED_HAPS`, `FILTERED_HAPS` | `"sum"` |
| `RAW_GT_COUNT`, `HAPCOMP`, `HAPDOM` | `"element_wise_sum"` |
| `AS_QUALapprox` | rebuilt as AS: type `["int"]`, `"element_wise_sum"` |
| `AS_RAW_MQ`, `AS_SB_TABLE` | `"disable_remap_missing_with_non_ref": true` |

How an "AS rebuild" works (`updateAlleleSpecificINFOFieldCombineOperation:232-279`):
* it clears `length`, `vcf_delimiter` and `type`;
* it adds length `R`,`var`, delimiters `|`,`,` and the type(s);
* it sets the combine operation (`histogram_sum` → `[float,int]`; `element_wise_float_sum` →
  `[float]` + `element_wise_sum`; otherwise `[int]` + `element_wise_sum`);
* `vcf_field_class` is untouched, so a field that is not INFO keeps its class.

Unlike 2.2, these overrides also apply when the field is not INFO. `vcf_type` is never set.

### 2.4 Worked example: every entry of ws_rb (45) and ws_raw (28) [V]
ws_rb, from HG002.rb.g.vcf.gz's header (which is already in htsjdk sorted order: ALT, FILTER, FORMAT,
GATKCommandLine, GVCFBlock, INFO, contig, source):

```
#   name                    type          class           length      extra
0   ID                      char          -               var
1   LowQual                 Integer       FILTER          1
2   AD                      Integer       FORMAT          R
3   DP                      Integer       FORMAT,INFO     1           (INFO DP line merged, 2.2)
4   GQ                      Integer       FORMAT          1
5   GT                      int           FORMAT          PP          (header: Number=1,String)
6   MIN_DP                  Integer       FORMAT          1
7   PGT, 8 PID              String        FORMAT          var         (Number=1,String)
9   PL                      Integer       FORMAT          G
10  PS                      Integer       FORMAT          1
11  SB                      Integer       FORMAT          4
12  AC                      Integer       INFO            A
13  AF                      Float         INFO            A
14  AN                      Integer       INFO            1
15-20 AS_BaseQRankSum AS_FS AS_InbreedingCoeff AS_MQ AS_MQRankSum AS_QD   Float INFO A
21  AS_QUALapprox           int           INFO            R,var       |,  element_wise_sum   (header String)
22  AS_RAW_BaseQRankSum     float,int     INFO            R,var       |,  histogram_sum
23  AS_RAW_MQ               float         INFO            R,var       |,  element_wise_sum + disable_remap
24  AS_RAW_MQRankSum        float,int     INFO            R,var       |,  histogram_sum
25  AS_RAW_ReadPosRankSum   float,int     INFO            R,var       |,  histogram_sum
26  AS_ReadPosRankSum       Float         INFO            A
27  AS_SB_TABLE             int           INFO            R,var       |,  element_wise_sum + disable_remap
28  AS_SOR                  Float         INFO            A
29  AS_VarDP                int           INFO            R,var       |,  element_wise_sum
30  BaseQRankSum            Float         INFO            1
31  END                     Integer       INFO            1
32-34 ExcessHet FS InbreedingCoeff        Float INFO 1
35  MQ, 36 MQRankSum        Float         INFO            1
37  MQ_DP                   Integer       INFO            1           sum
38  QD                      Float         INFO            1
39  QUALapprox              Integer       INFO            1           sum
40  RAW_GT_COUNT            Integer       INFO            3           element_wise_sum
41  RAW_MQandDP             Integer       INFO            2           element_wise_sum
42  ReadPosRankSum, 43 SOR  Float         INFO            1
44  VarDP                   Integer       INFO            1           sum
```

The order is the input order: FILTER, then FORMAT, then INFO. Inside each class it looks sorted
because htsjdk writes headers sorted by line text, e.g. `MQRankSum` < `MQ_DP` since `R` (0x52) <
`_` (0x5F). These numbers are the vidmap list positions, not the global field indices of
`genomicsdb_tiledb.md` 1.3: those insert `DP_FORMAT` and the `_tuple_element_k` entries.

* ws_raw (raw HaplotypeCaller GVCF) has 28 entries:
  * `ID LowQual AD DP GQ GT MIN_DP PGT PID PL PS SB`;
  * `AS_InbreedingCoeff AS_QD AS_RAW_BaseQRankSum AS_RAW_MQ AS_RAW_MQRankSum AS_RAW_ReadPosRankSum
    AS_SB_TABLE BaseQRankSum END ExcessHet InbreedingCoeff MLEAC MLEAF MQRankSum RAW_MQandDP
    ReadPosRankSum`.

  Every entry ws_raw shares with ws_rb is identical; MLEAC (Integer, A) and MLEAF (Float, A) are its own.
* ws_wg and ws_mid (another reblocked GVCF) have the ws_rb list minus the five `AS_RAW_*`/`AS_SB_TABLE`
  fields. ws11 and ws_rb_b1 are byte-identical to ws_rb.

### 2.5 Contigs [V]
* One entry per `##contig` line, in the order of 2.1:
  `{"name": ID, "length": length, "tiledb_column_offset": sum of the previous lengths}`.
  The first offset is 0. Lengths and offsets are int64 (strings in JSON).
* ws_rb: 195 contigs `chr1 (0) ... chr8 (1391350276) ... chrEBV (3099750718)`. ws_wg: chr20 2713028904,
  chr21 2777473071, chr22 2824183054.
* A contig line without `length` gets the htsjdk default length (`SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH` = 0)
  [S, unverified].
* C++ checks when the vidmap is loaded (`C/src/utils/vid_mapper_pb.cc:174-250`):
  * duplicate contig name → `ProtoBufBasedVidMapperException : Duplicate contigs found: <name>`;
  * overlapping column ranges → `... Overlapping contigs found`;
  * duplicate field name → `... Duplicate fields exist in vid map`.

### 2.6 `--merge-contigs-into-num-partitions N` (N > 0, not incremental) [S]
`GenomicsDBImport.traverse:875-879` → `GenomicsDBImporter.coalesceContigsIntoNumPartitions:913-965`.

1. Every interval must be a whole contig: start 1 and end ≥ length. Otherwise:
   `UserException("--merge-contigs-into-num-partitions requires that entire contigs be specified for input intervals. Input interval contained: Contig:%s, Start:%d, End:%d while reference contig was: Contig:%s, Start:1, End:%d")`.
2. **The vidmap contigs are cut down to the `-L` contigs** (in vidmap order) and the offsets are
   **recomputed from 0** (`regenerateVidUsingContigList:877-906`), but only if that changes the number
   of contigs. `vidmap.json` is rewritten. All column coordinates change accordingly.
3. Partitions are cut by TileDB column:
   * `genomeLength = lastOffset + lastLength + 1` and `goal = genomeLength / N` (integer division);
   * walk the contigs adding lengths; when the running sum ≥ goal, close the partition at that
     contig's last column;
   * the last partition ends at `genomeLength - 1`, one column past the last contig.

   Example, `-L chr20 -L chr21 -L chr22`, N=2: offsets 0 / 64444167 / 111154150; partitions
   `[0, 111154149]` and `[111154150, 161972618]`.
4. Arrays are named by rank, `"0"`, `"1"`, ... (5.1). With `--genomicsdb-update-workspace-path` the
   option is ignored, with a warning.

### 2.7 Exact JSON text (protobuf-java 3.7.1 `JsonFormat.printer().preservingProtoFieldNames()`) [V]
The GATK jar ships protobuf-java 3.7.1 and `JsonFormat$PrinterImpl$GsonHolder` (a default `Gson`,
i.e. HTML-safe escaping).

* **Field order** = declaration order in the .proto file (`genomicsdb-src/src/resources/*.proto`):
  * `GenomicsDBFieldInfo`: `name, type, vcf_field_class, vcf_type, length, vcf_delimiter,
    VCF_field_combine_operation, vcf_name, disable_remap_missing_with_non_ref`;
  * `Chromosome`: `name, length, tiledb_column_offset`;
  * `VidMappingPB`: `fields, contigs`;
  * `SampleIDToTileDBIDMap`: `sample_name, row_idx, idx_in_file, stream_name`;
  * `CallsetMappingPB`: `callsets`.
* **What is printed:** only present fields (proto2) and non-empty repeated fields.
* **Value forms:**
  * int64 → quoted decimal (`"length": "248956422"`); int32 → bare number;
  * bool → `true`;
  * strings get Gson HTML-safe escaping: `"`→`\"`, `\`→`\\`, `\t \b \n \r \f`; other chars < 0x20,
    and `< > & = '` U+2028 U+2029 → `\u00XX` (lower-case hex); everything else is raw UTF-8.

  The escaping is [S]: no test string contains such characters. It matters for sample names and
  paths, e.g. `O'Brien` → `O'Brien`.
* **Layout** (2-space steps, `, ` between array items, no newline at end of file):

```python
def pmsg(m, ind=0):                      # m: ordered (key, value) pairs, absent/empty ones removed
    if not m: return "{\n" + " "*ind + "}"
    pad = " " * (ind + 2)
    return "{\n" + ",\n".join(pad + jstr(k) + ": " + pval(v, ind + 2) for k, v in m) + "\n" + " "*ind + "}"
def pval(v, ind):
    if list:    return "[" + ", ".join(pval(x, ind) for x in v) + "]"
    if message: return pmsg(v, ind)      # nested message: indent of the *field* + 2
    ...
```

which gives

```
{
  "fields": [{
    "name": "ID",
    "type": ["char"],
    "length": [{
      "variable_length_descriptor": "var"
    }]
  }, {
  ...
    "VCF_field_combine_operation": "element_wise_sum",
    "disable_remap_missing_with_non_ref": true
  }, ...],
  "contigs": [{
    "name": "chr1",
    "length": "248956422",
    "tiledb_column_offset": "0"
  }, {
  ...
  }]
}
```

---------------------------------------------------------------------------------------------

## 3. `callset.json` [V]

Built by `CallSetMapExtensions.generateSortedCallSetMap` (`J/importer/extensions/CallSetMapExtensions.java:68-167`),
written at the very end (section 1, step 7), same printer as 2.7.

* **Rows:** samples in `SampleNameMap` order, i.e. a Java `TreeMap` of the sample names
  (`String.compareTo`, UTF-16 code units). `row_idx` = 0, 1, 2, ... (`lbRowIdx + i`, 0 for a new
  workspace). The order of the map file does not matter: rb.map lists HG002, NA12878, HG002pcr, and
  the rows are HG002=0, HG002pcr=1, NA12878=2.
* **`idx_in_file`** is always `"0"`.
* **`stream_name`:**
  * default path: `<sample>_stream`;
  * `--bypass-feature-reader`: the VCF location exactly as given in the map, or the `-V` string (`URI.toString()`)
    (ws_rb_bypass: `/home/kban/.../HG002.rb.g.vcf.gz`).

```
{
  "callsets": [{
    "sample_name": "HG002",
    "row_idx": "0",
    "idx_in_file": "0",
    "stream_name": "HG002_stream"
  }, {
  ...
  }]
}
```

* **Batches** (`--batch-size B`) do not change the file. There is one callset map over all samples,
  and batch k imports rows `[k·B, k·B+B-1]` of it. ws_rb_b1 (B=1) has the same `callset.json` as ws_rb.
* **Incremental import** (`--genomicsdb-update-workspace-path`) [S]:
  * the new samples (sorted) get rows starting at `1 + max_valid_row_idx_in_array` of the first
    partition's array (`GenomicsDBImporter.java:187-215`);
  * they are appended after the existing entries (`mergeCallsetsForIncrementalImport:249-264`);
  * a sample name already present →
    `GenomicsDBException("Duplicate sample name found: <s>. Sample was originally in <uri>")`;
  * `callset.json.fragmentlist` (the existing fragment names, `\n`-joined) is written first and
    deleted at the end;
  * `vidmap.json` and `vcfheader.vcf` are not rewritten.

---------------------------------------------------------------------------------------------

## 4. `vcfheader.vcf`

### 4.1 Content and text [V]
`JsonFileExtensions.writeVcfHeaderFile` (`J/importer/extensions/JsonFileExtensions.java:126-140`):
`new VCFHeader(mergedHeaderLines)` (no samples) is written with
`new VariantContextWriterBuilder().clearOptions().setOutputStream(..)` → `VCFWriter.writeHeader`
(`H/variant/variantcontext/writer/VCFWriter.java:200-240`).

* `mergedHeaderLines` = the lines of 2.1 (first sample's header for the sample-name map, in input order,
  `##fileformat` dropped) plus GATK's default tool lines (`G/.../GATKTool.java:985-993` →
  `G/utils/variant/GATKVariantContextUtils.java:70-82`):
  * `##source=GenomicsDBImport`;
  * `##GATKCommandLine=<ID=GenomicsDBImport,CommandLine="...",Version="4.7.0.0",Date="...">`.

  Both are absent with `--add-output-vcf-command-line false`.
* Output:
  * `##fileformat=VCFv4.2`;
  * then each line once, **sorted** by its text after `##`. This is Java `String.compareTo`, except
    that two `contig` lines compare by their position in the header, so the contig block keeps header
    order and sits where `contig=<` sorts (after `INFO`, before `source`);
  * then `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n` with **no FORMAT/sample columns**;
  * the file ends with `\n`.
* `variant-rs/src/vcfheader.rs:write_header(input_meta, extra, &[])` produces exactly this.
* Each line is printed with htsjdk's `toString`. Lines that htsjdk itself wrote round-trip byte for
  byte. Other inputs may be normalized:
  * compound lines are rebuilt as `ID,Number,Type,Description[,Source][,Version]`, so other keys are
    dropped;
  * values with `,` or space, and `Description`/`Source`/`Version`, are quoted (`H/variant/vcf/VCFHeaderLine.java:145-171`).
* Diff of ws_rb/vcfheader.vcf against HG002.rb.g.vcf.gz's header: the `GenomicsDBImport` command-line
  line is added (sorted between `##FORMAT=<ID=SB` and `##GATKCommandLine=<ID=HaplotypeCaller`),
  `##source=GenomicsDBImport` is added before `##source=HaplotypeCaller`, and the `FORMAT\tHG002`
  columns are gone. Rebuilt byte-exactly for the six workspaces listed at the top.
* With `-V`, the lines are `smartMergeHeaders`' output (2.1): duplicate `##source=` lines are
  collapsed, and contigs from all inputs are kept in first-seen order [S].

### 4.2 The `##GATKCommandLine` line [V format; content is GATK-specific]
```
##GATKCommandLine=<ID=GenomicsDBImport,CommandLine="GenomicsDBImport <args>",Version="4.7.0.0",Date="November 20, 2026 at 9:37:05 PM SGT">
```
`<args>` is Barclay's `getCommandLine()`, read off the three logs:
* first every **explicitly given** argument, then every other argument with a non-null value
  (defaults), each as `--long-name value` (booleans `true`/`false`, the empty string as nothing,
  hence `--gcs-project-for-requester-pays  --disable...`);
* within each group, in field declaration order: GenomicsDBImport's own fields
  (`GenomicsDBImport.java:228-393`) first, then GATKTool's (`--intervals`, `--interval-set-rule`, ...);
* `-L` appears as `--intervals <string>`; a flag given bare shows as `--consolidate true`.

The ws_rb line in full is in `ws_rb/vcfheader.vcf:14`. The key order is ID, CommandLine, Version, Date
(a `HashMap`'s iteration order). The project convention (`vcfheader.rs:command_line`) writes our own
command line in this position, which is the only intended difference.

---------------------------------------------------------------------------------------------

## 5. Arrays (one per partition)

### 5.1 From `-L` to partitions and array names [V for the examples]
1. **Parse** `-L` with GATK's standard interval collection against the merged header's dictionary
   (2.1): `--interval-set-rule UNION`, `--interval-padding 0`, `--interval-merging-rule ALL`,
   `-XL`. The result is sorted in dictionary order, and overlapping or **abutting** intervals are
   merged (`G/.../GenomicsDBImport.java:1077-1110`; `variant-rs/src/intervals.rs` implements this).
   * A whole contig `-L chr20` becomes `chr20:1-64444167`.
   * No `-L` (and not incremental) → `UserException("No intervals specified")`.
   * More than 100 intervals without `--merge-input-intervals` → a warning only.
2. **`--merge-input-intervals`** → `IntervalUtils.getSpanningIntervals` (`G/utils/IntervalUtils.java:1326-1331`):
   one interval per contig, from the smallest start to the largest end, sorted by dictionary. Data in
   the gaps is imported too.
3. Each interval becomes one partition, `begin = {contig, start}` and `end = {contig, end}`, with
   `generate_array_name_from_partition_bounds` (`GenomicsDBImport.java:747-787`).
4. **Array name** = `String.format("%s$%d$%d", contig, start, end)`, the 1-based inclusive interval
   after steps 1-2 (`GenomicsDBImporter.updateConfigPartitionsAndLbUb:734-754`). With
   `--merge-contigs-into-num-partitions`: `"0"`, `"1"`, ... by rank (2.6).

| command | arrays |
|---|---|
| `-L chr8:1-15000000` | `chr8$1$15000000` (ws_rb) |
| `-L chr20 -L chr21 -L chr22` | `chr20$1$64444167`, `chr21$1$46709983`, `chr22$1$50818468` (ws_wg) |
| `-L chr21:20000000-30000000` | `chr21$20000000$30000000` (ws_mid) |
| `-L chr1:100-200 -L chr1:500-600` | `chr1$100$200`, `chr1$500$600`; with `--merge-input-intervals`: `chr1$100$600` [S] |
| `-L chr1:100-200 -L chr1:201-300` | `chr1$100$300` (abutting, merged by rule ALL) [S] |

Column range of a partition (`C/src/config/genomicsdb_config_base.cc:63-83`):
`[offset+start-1, offset+end-1]`.
* An end beyond the contig length is clamped, with a warning. GATK's interval validation prevents
  this case.
* If partitions overlapped, the earlier one's end would be cut to the next one's begin - 1
  (`C/src/config/pb_config.cc:173-182`). This cannot happen after merge rule ALL.

### 5.2 `genomicsdb_meta_dir/genomicsdb_column_bounds.json` [V]
Written once, when the array is first opened for writing, if absent (`C/src/genomicsdb/variant_storage_manager.cc:852-877`).
The format is rapidjson `PrettyWriter`: 4-space indent, `"key": value`, int64 numbers, no newline at end.

```
{
    "min_column": 1391350276,
    "max_column": 1406350275
}
```

| workspace | min_column | max_column |
|---|---|---|
| ws_rb | 1391350276 | 1406350275 (= chr8 offset + 15000000 - 1) |
| ws_mid | 2797473070 | 2807473070 (chr21 offset 2777473071 + 20000000 - 1, + 30000000 - 1) |
| ws_wg chr20 | 2713028904 | 2777473070 |

### 5.3 `genomicsdb_meta_dir/genomicsdb_meta_<uuid>.json` [V]
Written when the array is opened for a batch (`C/src/loader/load_operators.cc:150-151` →
`variant_storage_manager.cc:424-466`).

* The values passed in are `lb = lb_callset_row_idx` of the batch and
  `max = min(ub_callset_row_idx, max row in callset map)`, where `ub = lb + batchSize - 1`.
* A new file is written only if `max` exceeds the stored maximum, or no file exists yet, or `lb` is
  below the stored minimum. Its content is then the **running** min(lb) and max(max).
* `<uuid>` is a random v4 UUID, lower case (`uuid_generate` + `uuid_unparse`).
* Same format as 5.2:

```
{
    "lb_row_idx": 0,
    "max_valid_row_idx_in_array": 2
}
```

* ws_rb (one batch of 3): one file, 0/2.
* ws_rb_b1 (batch size 1, 3 batches, no consolidation): three files, 0/0, 0/1, 0/2. Every batch
  raised the maximum, and `lb` stays 0 because later batches' lb (1, 2) is not below the stored 0.
* ws11: 0/10. ws_wg and ws_mid: 0/0 per array.

### 5.4 Fragments and `--consolidate` [V]
* One fragment per batch per array (ws_rb_b1 has 3 fragments). Fragment naming and the merged view
  are in `genomicsdb_tiledb.md` 3-4.
* `--consolidate`: only after the **last** batch (`GenomicsDBImporter.java:659-663`), per array
  (`variant_storage_manager.cc:468-520`):
  1. read every `genomicsdb_meta*.json`, delete them, and write **one** new file with the overall
     min/max (new uuid);
  2. `tiledb_array_consolidate`, with buffer size = segment size.
* With a single fragment the fragment is left as is. ws_rb, ws_raw and ws11 keep the fragment whose
  timestamp is the start of the write (15:37:29.8), while consolidation ran at 15:37:41.
* The merged cells of ws_rb (batch 50, consolidated) and ws_rb_b1 (batch 1, 3 fragments) are
  identical (all 1 359 702 dumped cells).

### 5.5 Other files [V]
* `.__consolidation_lock`: empty, mode 0700, created with the array (`T/src/storage/storage_manager.cc:486-492`).
  It is present in every array, including ws_rb_bypass, which ran with shared-posixfs on.
* `__tiledb_workspace.tdb`: empty, workspace root.
* `__array_schema.tdb`: see `genomicsdb_tiledb.md` 2.
* An array directory contains nothing else.

### 5.6 What the options change

| option | effect on workspace content |
|---|---|
| `--batch-size B` (0 = all) | one fragment and possibly one meta file per batch (5.3, 5.4); merged cells unchanged [V] |
| `--consolidate` | one fragment per array (if >1 batch), one meta file; cells unchanged [V] |
| `--reader-threads n` | none. Used only with exactly 1 interval: readers are opened in parallel and their query started early (`InitializedQueryWrapper`, `GenomicsDBImport.java:701-718,1141-1186`); with more intervals there is a warning and serial opening [V log of ws_wg] |
| `--merge-input-intervals` | array names/ranges (5.1) and which records are fed |
| `--genomicsdb-shared-posixfs-optimizations` | none: TileDB disables file locking and keeps write handles open [V ws_rb_bypass] |
| `--bypass-feature-reader` | `callset.json` stream names (3); reader path (6.3); cells identical on all test data [V] |
| `--genomicsdb-segment-size` (1048576) | none: size of TileDB write buffers and the consolidation buffer; tiles are cut every 1000 cells regardless [S] |
| `--genomicsdb-vcf-buffer-size` (16384) | none, unless too small: `size_per_column_partition = vcfBufferSize × min(B, nSamples)`; Java stream buffer = that / readers in batch (grows for a single oversized record); C++ per-sample cell buffer = that / samples. A record whose cells do not fit an empty buffer → fatal `Buffer(size_per_column_partition) did not have space to hold a line fully - increase buffer size` (`C/src/loader/tiledb_loader_file_base.cc:446`) [S] |
| `--overwrite-existing-genomicsdb-workspace` | an existing **workspace** (has `__tiledb_workspace.tdb`) is deleted recursively and recreated; a non-workspace directory still fails (6.4) [S] |
| `--max-num-intervals-to-import-in-parallel n` | none: arrays of a batch are imported by n threads (fragment timestamps differ) |
| `--merge-contigs-into-num-partitions N` | vidmap contigs and offsets, array names and ranges (2.6) |
| `--validate-sample-name-map` | none in this path: it is stored in ImportConfig but not used by `GenomicsDBImporter` [S] |
| `--add-output-vcf-command-line false` | no `##source`/`##GATKCommandLine` lines in vcfheader.vcf |

Fixed by GenomicsDB, not settable from GATK (`GenomicsDBImporter.addExplicitValuesToImportConfiguration:372-383`):
`treat_deletions_as_intervals=true`, `compress_tiledb_array=true`, `num_cells_per_tile=1000`,
`row_based_partitioning=false`, `produce_tiledb_array=true`. GATK sets
`fail_if_updating = !incremental` (only on the first batch) and `consolidate`/`shared-posixfs` as given.

---------------------------------------------------------------------------------------------

## 6. Import flow

### 6.1 Loops (`GenomicsDBImporter.executeImport:606-624`, `iterateOverSamplesInBatches:648-680`) [S, V by logs]
```
B' = (batch-size <= 0) ? nSamples : min(batch-size, nSamples)
for lb in 0, B', 2B', ... < nSamples:                       # rows lb..lb+B'-1 (sorted sample order)
    for each partition (array) in parallel (max-num-intervals-to-import-in-parallel threads):
        open readers for the batch's samples (GATK FeatureReader, or none with --bypass)
        C++ loader: create array if needed, write column bounds / meta (5.2, 5.3)
        read each sample's records for the partition, convert to cells, merge in (column,row) order,
        write one fragment; on the last batch with --consolidate: consolidate
    "Done importing batch k/K"
write callset.json
```
Readers are reopened for every partition and batch.

### 6.2 Default path: htsjdk reads, then **VCF text** goes to the C++ loader [S, V by equal cells]
* Per sample: `AbstractFeatureReader.getFeatureReader(uri, index, new VCFCodec(), requireIndex=true)`
  (`GenomicsDBImport.java:971-1028`), then `query(contig, start, end)` for the partition
  (`GenomicsDBImporter.java:151`).
* The query returns every record that **overlaps** `[start, end]`, i.e. `rec.start <= end && rec.end >= start`.
  `rec.end` is INFO `END` if present, else `POS + len(REF) - 1` (`H/variant/vcf/AbstractVCFCodec.java:438-447`).
* **GATK passes `passAsVcf = true`** (`ImportConfig(..., validate, true, ...)`, `GenomicsDBImport.java:837`), so
  `addSortedVariantContextIterator` uses `OutputType.VCF_STREAM` (`GenomicsDBImporter.java:170-171`).
  This corrects `genomicsdb_tiledb.md` 5, which says "BCF2 in memory".
* The record is re-serialized as **VCF text** by htsjdk's `VCFWriter` (`setOutputVCFStream`, no
  options except no indexing, `J/importer/GenomicsDBImporterStreamWrapper.java:66-90`). The C++ side
  parses that text with htslib (`bcf_hdr_deserialize` / `bcf_deserialize` with `is_bcf=0`,
  `C/src/vcf/vcf2binary.cc:76-101`).
* The stream header is **that sample's own header**: its htsjdk-repaired lines, sorted, `VCFv4.2`, with
  that file's sample column. It is not the merged header.
* For each vidmap field missing from a stream/file header, a definition is added
  (`VCFAdapter::add_field_to_hdr_if_missing`, `C/src/vcf/vcf_adapter.cc:65-140`). For multi-D/tuple
  fields the definition is forced to `Number=1,Type=String`.

How `VCFEncoder.write` (`H/variant/vcf/VCFEncoder.java:104-195`) rewrites each field:

| column | written as |
|---|---|
| CHROM, POS, ID | verbatim (ID `.` stays `.`) |
| REF | **upper-cased** by the codec (`AbstractVCFCodec.java:428`) |
| ALT | display strings; non-symbolic alleles **upper-cased** (`SimpleAllele.java:60-75`); none → `.` |
| QUAL | `.` if missing, else `String.format("%.2f", q)` with a trailing `.00` removed (`:214-220`): **rounded to 2 decimals** |
| FILTER | `.`, `PASS`, or the filters **sorted** and `;`-joined; a filter missing from the header → error (6.4) |
| INFO | keys sorted (`TreeMap`); values are the text from the file (comma lists re-joined); a Flag → bare key; `Flag=0` → dropped by the codec; a non-Flag key without `=` → `KEY=.`; `.` if empty |
| FORMAT + sample | the **raw text from the file**, passed through undecoded (`LazyGenotypesContext` with String data, `:165-169`) |

For GATK-written GVCFs (QUAL already 2 decimals, upper case, ≤1 filter) the text is the original text.
This is why `genomicsdb_tiledb.md` 9's validation against the source GVCFs had 0 mismatches.

### 6.3 `--bypass-feature-reader`: GenomicsDB reads the files with htslib [S, V by equal cells]
* GATK opens no readers (`sampleToReaderMapCreator = null`, `GenomicsDBImport.java:838`). The callset
  stream names are the file paths (3), and the C++ `VCFReader` opens each file itself
  (`C/src/vcf/vcf2binary.cc:130-237`).
  * It reads the header with `bcf_hdr_read` (the file's own header, **not** repaired by htsjdk, plus
    the vid fields added as in 6.2).
  * It builds a synced reader whose regions are all vidmap contigs present in the file.
  * It `bcf_sr_seek`s to `(contig, start-1)` and reads records until
    `offset + pos > partition end` (`seek_and_fetch_position:574-666`).
  * The sample at `idx_in_file` 0 is imported, i.e. the first sample column, whatever its name.
* GATK still reads the first sample's header (2.1), and the sample-name map requires `<path>` to
  end in `.vcf.gz` (case-insensitive) with a readable `<path>.tbi` (`GATKGenomicsDBUtils.java:347-358`):
  `UserException("Input variant files must be block compressed vcfs when using bypass-feature-reader, but <p> does not end with the standard file extension .vcf.gz")`.
  An explicit index column in the map →
  `UserException("Indices were specified for some VCFs in the sample name map file, but --bypass-feature-reader was also specified. ...")`.
* **Differences from 6.2** [S]. None of them occurs in the test data, and cells were identical for
  ws_rb vs ws_rb_bypass (1 359 702 cells) and ws_raw vs ws_raw_bypass (8 854 450 cells) [V].
  * No MNP check (6.4).
  * QUAL is not rounded to 2 decimals.
  * REF/ALT are not upper-cased.
  * FILTER is kept in file order.
  * INFO/FORMAT keys missing from the file header are accepted by htslib (a dummy definition and a
    warning) instead of an htsjdk error. They are still not stored unless in the vidmap.
  * Standard-header repair does not happen. If a file declares a vidmap field with another type,
    `bcf_get_*_values` returns -2 and the field is stored as missing (`vcf2binary.cc:722`). The same
    applies in 6.2 when a non-first sample's header disagrees with the vidmap.
  * Two samples mapped to the same file →
    `ProtoBufBasedVidMapperException : Attempting to import a sample from file/stream <f> multiple times under aliases '<a>' and '<b>' ...`
    (`vid_mapper_pb.cc:108-115`). In 6.2 the stream names differ, so this is allowed.
* With `--avoid-nio` (requires `--bypass-feature-reader` and `--header`), GATK does not open the inputs at all.

### 6.4 Validation and error messages

**Workspace** (`GenomicsDBImport.java:1038-1062`, `T/src/c_api/tiledb_utils.cc:103-136`) [S]
* Path is a file →
  `"Error creating GenomicsDB workspace: <ws> already exists and is not a directory"`.
* Existing workspace without overwrite → `"... <ws> already exists"`.
* Existing directory that is not a workspace (even an empty one, since `mkdir` fails), or a workspace
  inside another workspace/array → `"Error creating GenomicsDB workspace: <ws>"`. The tool docs say
  "empty or non-existent", but the code requires non-existent.
* Incremental import with a missing workspace →
  `"We require an existing valid workspace when incremental import is set"`.

**Sample-name map** (`G/.../SampleNameMap.java:107-159`) [S]
* Lines are split on tab with `split("\\t",-1)`, so a trailing tab counts as a field.
* Not 2 or 3 fields → `BadInput("Sample name map file must have 2 or 3 fields per line in the format:\nSample\tFile\nor:\nSample\tFile\tIndex\nbut found line: \"<l>\" with <n> fields")`.
* Empty file → `"At least 1 sample is required but none were found in the sample mapping file"`.
* Sample name empty, or with leading/trailing whitespace, or blank file → `"Sample name map file must have lines in the format: ..."`.
* The file and index are `trim()`med. A blank index → `"Found a line in the sample name map file with an empty or all-whitespace value for the index: ..."`.
* Duplicate sample → `"Found two mappings for the same sample: <s>\n<path>\n<previous>"`.
* Not a URI → `"Malformed URI: ..."`.
* The name in the map is **not** checked against the file's sample (renaming is allowed).

**`-V`** (`GenomicsDBImport.java:559-583`, `SampleNameMap.addSample:191-211`) [S]
* The sample name is the header's first sample.
* A header with ≠1 sample →
  `"Input GVCF: <p> was expected to contain a single sample but actually contained <n> samples."`.
* Duplicate sample → `"Duplicate sample: <s>. Sample was found in both <p2> and <p1>."`.
* `-V` and `--sample-name-map` are mutually exclusive. One of them is required:
  `"One of --variant or --sample-name-map must be specified"`.

**Header** [S]
* No contig lines → `"The merged vcf header has no sequence dictionary. Please provide a header that contains a sequence dictionary."`.
* The `-L` contigs must be in that dictionary (GATK interval parsing errors).

**Records, default path** (`GenomicsDBImport.java:984-1023`) [S]
* `GATKVariantContextUtils.isUnmixedMnpIgnoringNonRef` (`G/utils/variant/GATKVariantContextUtils.java:2175-2186`)
  is true when len(REF) ≥ 2 and every allele is either `<NON_REF>` or has len(REF) bases. Then:
  `BadInput("GenomicsDBImport does not support GVCFs with MNPs. MNP found at <ctg>:<pos> in VCF <abs path>")`.
* INFO/FORMAT/FILTER key not in that file's header → `IllegalStateException("Key <k> found in VariantContext field <INFO|FORMAT|FILTER> at <ctg>:<pos> but this key isn't defined in the VCFHeader.  We require all VCFs to have complete VCF headers by default.")`.
  This comes from `VCFEncoder.java:222-227`, and an htsjdk lazy FORMAT passthrough skips the FORMAT check.
* VCF syntax errors are codec `TribbleException`s, e.g. `"The VCF specification does not allow for whitespace in the INFO field..."`.
* Reader creation failure → `"Failed to create reader from <uri> because of the following error:\n\t<msg>"`.

**C++ loader** (both paths) [S, except the first]
* Two records of one sample with the same or a decreasing start (so unsorted or duplicate-position
  input) → `LoadOperatorException("ERROR: two cells in incorrect order found\nPrevious cell: <row>, <begin>, <end>\nNew cell: <row>, <begin>")`
  (`C/src/loader/load_operators.cc:247-252`; see `genomicsdb_tiledb.md` 5.3).
* An input that is not sorted by contig/position cannot be tabix-indexed, so it fails earlier.
* A cell outside the partition → `LoadOperatorException("Found cell [ <row>, [ <b>, <e> ] ] that does not belong to TileDB/GenomicsDB partition with row_bounds [ .. ] column_bounds [ .. ]")`
  (`load_operators.cc:75-98`). It fires if a cell begins after the partition end, or if the **first**
  cell (column order) ends before the partition start. That needs a record overlapping the start only
  through len(REF) while its GenomicsDB END is shorter, e.g. a complex non-deletion allele (10.2).
* A record on a contig not in the vidmap → `VERIFY_OR_THROW(global_contig_idx >= 0)` (`vcf2binary.cc:563-572`).
* The array already exists while `fail_if_updating` → `"Array <ws>/<a> exists and flag \"fail_if_updating\" is set to true in the loader JSON configuration"`.
* Buffer too small → 5.6.

---------------------------------------------------------------------------------------------

## 7. Things that affect cell content that `genomicsdb_tiledb.md` does not cover

1. **Which fields are stored** = the vidmap (2.1/2.2), i.e. only fields declared in the first
   (sorted) sample's header, with the types and lengths derived there. GATK's AS overrides (2.3)
   decide how `AS_*` strings are parsed (tiledb 7.6).
2. **Rows** = sorted sample order (3), independent of map order and batches.
3. **Records fed per partition** [V]: every record overlapping `[start, end]` by htsjdk/tabix
   (6.2/6.3). Records that begin before `start` are handled by
   `handle_intervals_spanning_partition_begin` (`load_operators.cc:35-72`):
   * per row, only the **last** such record (in file order) whose GenomicsDB END ≥ partition start is
     kept;
   * a later record ending before the start cancels it;
   * the kept record is written at its original begin column (< `min_column`), with its END copy.

   ws_mid shows this: the first cell is row 0, column 2797472380 = chr21:19999310 (a ref block, END
   chr21:20000907), stored below `min_column` 2797473070, and its END copy is at chr21:20000907.
   Records beginning in the interval and ending after it are kept whole (tiledb 5.4, ws11).
4. **Default-path text normalization** (6.2): QUAL rounded to `%.2f`, REF/ALT upper-cased, filters
   sorted. These are no-ops for GATK GVCFs but would differ from a direct htslib parse of the file.
   A port that parses the GVCF text itself must apply them to match GATK's default path.
   `--bypass-feature-reader` skips them.
5. **Partition bounds** from the vidmap offsets. With `--merge-contigs-into-num-partitions` these are
   *recomputed* offsets (2.6), so the same position gets a different column than without the option.
6. **Batching and consolidation** do not change the merged cells (5.4, [V]).

---------------------------------------------------------------------------------------------

## 8. Checklist for a byte-identical port (metadata)

* `__tiledb_workspace.tdb`: empty.
* `vidmap.json`: 2.2 + 2.3 + 2.5 printed per 2.7, no newline at end.
* `callset.json`: section 3 printed per 2.7, no newline at end.
* `vcfheader.vcf`: section 4 (`vcfheader::write_header`, no samples).
* Per array:
  * directory name per 5.1;
  * `__array_schema.tdb` (tiledb 2);
  * empty `.__consolidation_lock`;
  * `genomicsdb_meta_dir/genomicsdb_column_bounds.json` (5.2);
  * one `genomicsdb_meta_<uuid4>.json` per row-bound increase, or one after consolidation (5.3);
  * fragments (tiledb 3).
* Modes 0700 if you want `ls -l` identical. Write `callset.json` last.

---------------------------------------------------------------------------------------------

## 9. Verification log

| check | result |
|---|---|
| vidmap.json rebuilt from first sorted sample's header (2.2-2.3, 2.5, 2.7) | byte-identical: ws_rb, ws_raw, ws_rb_b1, ws11, ws_rb_bypass, ws_raw_bypass |
| callset.json rebuilt (3), incl. bypass stream names | byte-identical, same six |
| vcfheader.vcf rebuilt (4.1), GATKCommandLine line copied | byte-identical, same six (ws_wg/ws_mid input missing) |
| column bounds / meta files / lock / workspace file / modes | as described, all workspaces |
| ws_rb vs ws_rb_bypass, ws_raw vs ws_raw_bypass, ws_rb vs ws_rb_b1 (`tiledb_read.py dump`) | identical dumps (1 359 702 / 8 854 450 / 1 359 702 cells) |
| ws_mid record before the interval start | stored at chr21:19999310 with END copy (7.3) |

---------------------------------------------------------------------------------------------

## 10. Open questions / not verified

1. **First sample by sorted order vs map order** for the merged header: in all test maps the first
   line is also the first sorted name, so the data cannot tell them apart. The source (`TreeMap`
   iterator, `GenomicsDBImport.java:600`) says sorted.
2. **Complex alleles overlapping the interval start only through REF length** (e.g. `ATG→CT` starting
   1 bp before the start) can make the first cell end before the partition. That throws if it is the
   first cell in column order, and is otherwise dropped. Not exercised.
3. **Gson escaping** (2.7) is inferred from the jar contents (default `Gson`, no
   `disableHtmlEscaping`). No test string contains an escaped character.
4. **`Number=.` → `"VAR"`, `Character` type, Flag, FILTER/INFO name collisions, `-V` path,
   `--header`, incremental import, `--merge-contigs-into-num-partitions`**: source only.
5. **Existing empty directory as workspace path**: the code (TileDB `create_dir` fails if the
   directory exists) contradicts the tool doc ("empty or non-existent"). Not run.
6. **Segment size and fragment bytes**: the claim that tiles are cut every 1000 cells regardless of
   write-call boundaries is from TileDB's sorted-write design and was not checked byte for byte.
   Cells are unaffected either way.
