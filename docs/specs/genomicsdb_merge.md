# GenomicsDB query/merge semantics as seen by GATK GenotypeGVCFs (`-V gendb://ws`)

GATK 4.7 + GenomicsDB 1.5.5. This document says which VariantContexts (VCs) the GenotypeGVCFs walker receives,
given each sample's stored GVCF records. What is stored (cells, END copies, field values) is specified separately;
here a sample's stored data is taken to be its GVCF records as GenomicsDBImport read them, with the vidmap fields
parsed to their vidmap types (Float values are float32).

Verified by `merge_proto.py`, which rebuilds the received stream from the input GVCFs: 100% whole-line match on
all four workspaces (§12).

Sources: `GenotypeGVCFs.getGenomicsDBOptions`, `GATKGenomicsDBUtils.createExportConfiguration`,
`FeatureDataSource`/`FeatureCache`/`VariantLocusWalker`, and GenomicsDB `GenomicsDBFeatureReader/Iterator`,
`query_variants.cc` (scan_and_operate, handle_gvcf_ranges, gt_get_column), `broad_combined_gvcf.cc`,
`variant_operations.cc`, `variant_field_handler.cc` and `known_field_info.cc`.

## 0. Pipeline in one paragraph

GenotypeGVCFs opens a `GenomicsDBFeatureReader` with `VCFCodec`. The default is VCFCodec, not BCF2Codec:
`--genomicsdb-use-bcf-codec` is false. Each reader query runs GenomicsDB's `BroadCombinedGVCFOperator` over the
array. The operator writes VCF **text**: htslib formats the record, then GenomicsDB streams it through JNI.
htsjdk's `VCFCodec` parses that text into VCs. Therefore:

* Every INFO value is a `String`, or an `ArrayList<String>` when the text contains `,`.
* FORMAT `GT`/`AD`/`DP`/`GQ`/`PL` are typed. They are decoded into a `FastGenotype` with alleles, `int[]` and
  `int`.
* Every other FORMAT value (`SB`, `MIN_DP`, `PGT`, `PID`, `PS`) is a `String` extended attribute. For example,
  `SB` is the String `"0,0,4,0"`, not a list.

Genotypes are a `LazyGenotypesContext`, so this matches a multi-sample GVCF file read with VCFCodec. The
`received/<ws>.types.txt` file has the full class census.

## 1. Export configuration (what GenotypeGVCFs asks for)

| setting | value |
|---|---|
| produce_GT_field | **true** (GenotypeGVCFs forces `callGenotypes=true`) |
| produce_GT_with_min_PL_value_for_spanning_deletions | false |
| produce_FILTER_field | unset, so **false** |
| sites_only_query | false |
| max_diploid_alt_alleles_that_can_be_genotyped | `--genomicsdb-max-alternate-alleles`, default **50** |
| max_genotype_count | `--max-genotype-count`, default **1024** |
| query attributes | none given, so every vidmap field is queried |
| arrays | `generate_array_name_from_partition_bounds`: the workspace directories `contig$begin$end` |
| vid/callset/header | `vidmap.json`, `callset.json`, `vcfheader.vcf` of the workspace |

* **Header.** The header GenotypeGVCFs sees is `vcfheader.vcf` plus one added line,
  `##FILTER=<ID=PASS,...>`. It is the merged header written by GenomicsDBImport, so INFO/FORMAT `Type`s come from
  the input GVCF headers. For example, `AS_RAW_MQ`, `AS_SB_TABLE`, `AS_QUALapprox`, `AS_VarDP` and
  `AS_RAW_*RankSum` are `Number=1,Type=String`; `RAW_MQandDP` is Integer/2; `RAW_GT_COUNT` is Integer/3.
* **Sample order.** Samples appear in `callset.json` `row_idx` order. GenomicsDBImport assigns rows in sorted
  sample-name order.

## 2. Query windows: what the walker actually queries

GenotypeGVCFs is a `VariantLocusWalker`. Its driving `FeatureDataSource` has cache lookahead
`L = 100_000`.

**Traversal intervals.** With `--merge-input-intervals`, the -L intervals collapse to one spanning interval per
contig (`IntervalUtils.getSpanningIntervals`). Each traversal interval `[S,E]` is split into shards of `L` bases
starting at `S`: `[S, S+L-1]`, `[S+L, S+2L-1]`, … The last shard is clipped at `E`.

**Cache windows.** For each shard the walker first calls `query(shard)`:

* On a cache miss, `FeatureCache` refills with `reader.query(contig, shardStart, shardEnd + L)`.
* The next shard is inside that window, so it is a hit. The shard after that misses again.

The windows are therefore `[S + 2kL, min(S+2kL+L-1, E) + L]`: `[1,200000]`, `[200001,400000]`, … for `S=1`.
A window is refilled even when it contains no records.

**Array clipping.** `GenomicsDBFeatureIterator` splits a query over the arrays (directories) overlapping it,
sorted by `ChrArrayFolderComparator`. It clips the query to each array, `[max(qb, arrBegin), min(qe, arrEnd)]`,
and concatenates the results.

Each array query is an independent GenomicsDB scan (§3). Its first record can start exactly at the window start.
A block or deletion that began earlier is **clipped**: it reappears as a new record starting at the window start.
The last record is clipped at the window end; for example, `chr8:198163 END=200000` is followed by
`chr8:200001 END=200336`.

**Which VCs a locus sees.** At locus `X`, the walker sees the VCs from the window containing `X`'s shard that
overlap `X`. Overlap means `vc.start <= X <= vc.end`, with `vc.end` from END or `start+len(REF)-1` as VCFCodec
computes it.

Within one window, records never start at the same position. The only overlaps come from a record at a
deletion's start. That record has no END and its REF is the full deleted sequence, so its VCFCodec end reaches
past the following 1-bp records (§3). A deletion record from the previous window is **not** visible after the
window boundary.

**Received stream.** The stream is the concatenation of the per-window GenomicsDB outputs, keeping records with
`start <= E`. Records cannot start before `S`. `received/<ws>.vcf` holds exactly this stream. The harness writes
each VC object once, and no object is first received after its start (`late 0` in `received/<ws>.events.txt`).

## 3. Segmentation of one query `[a,b]` (`scan_and_operate`)

### Per-sample calls

State: one *call* slot per sample (row). Each call has `begin`, `end` and fields.

* A stored record's `end` is its INFO `END` if present, otherwise `POS + len(REF) - 1`.
* A record **contains a deletion** if `len(REF) > 1` and some ALT is `*`, or is a non-symbolic allele shorter
  than REF.
* A record **contains an MNV** if `len(REF) > 1` and some non-symbolic ALT of equal length differs from REF at
  more than one base.
* *Symbolic* means `<NON_REF>`, `*`, `<...>`, or an allele containing `[` or `]`.
* A **reference block** has `len(REF) == 1` and ALT exactly `[<NON_REF>]`.

The deletion/MNV flag and the ref-block flag are computed once, when the record is filled into the slot.

### Start of the query (`gt_get_column(a)`)

**Stored cells, as written by the loader (`load_operators.cc LoaderArrayWriter`).** Each record has:

* a *begin cell* at `POS`, which keeps the record's original END;
* an *END copy* at `copyEnd`, only when `copyEnd > POS`.

`copyEnd` is `end`, except when the sample's **next** record begins at or before `end`. In that case the loader
truncates `copyEnd` to `next.POS - 1`, and drops the copy if that equals `POS`.

Example: a 50-bp deletion at 149869 (end 149918) followed by a record at 149918 has its END copy at 149917.

**Per sample.** Find the first of its cells at coordinate `>= a`, in column order.

* If that cell is a begin cell at exactly `a`, the sample starts with that record as its call.
* If it is an END copy, the sample starts with that record as its call, ending at `copyEnd`. This is a record
  that spans `a`; when its END copy was truncated, it is a deletion whose remainder was never replaced.
* Otherwise the sample has no call at `a`.

This gives the same result as the forward-scan replacement rule below, so window starts behave consistently.

If any sample has a call, `cur = a`. Otherwise `cur` is the first forward cell's begin.

### Forward scan

Begin cells with `a < POS <= b` are visited in (POS, row) order. END copies are skipped. For each cell:

1. If `POS != cur`, flush ranges up to `POS - 1` (below), then set `cur = POS`.
2. **Replace** the sample's slot with this record. An overlapping older call is dropped completely; its
   remainder is never resumed.
   * Example: HC writes a deletion at 266482 and a `*`/variant record at 266484 for the same sample. From 266485
     until the deletion's end, that sample has **no call** and its genotype is `.` (§7).
   * If the dropped call was neither a deletion/MNV nor a reference block, GenomicsDB aborts with
     "Unhandled overlapping variants".

At the end of the query, flush ranges up to `b`.

### Flush ranges (`handle_gvcf_ranges(cur, next)`)

While some sample has an active call and `cur < next`:

1. Set `m = min(min active end, next-1)`.
2. If any **active** call contains a deletion or MNV, set `m = cur`. This single-position stepping applies to all
   samples, including at the deletion's own start position.
3. Emit one record for `[cur, m]` (§4–§8).
4. Drop the calls whose `end == m`.
5. Set `cur = m+1`.

If no sample has an active call, nothing is emitted for that gap.

### Consequences

* A record starts at every position where any sample's record begins.
* A record also starts the day after any active call ends (`end+1`) and at every window start.
* Each record spans to the next such breakpoint minus 1 (`END`), unless a deletion/MNV is active. In that case
  every position is its own 1-bp record with no END.

## 4. Spanning deletions (`handle_deletions`, run first for each emitted record)

This applies to every active call that contains a deletion or MNV and has `call.begin < s`, where `s` is the
record start:

1. **Invalidate the call's INFO fields.** This covers every vidmap INFO field with a sum, median, mean,
   element_wise_sum, concatenate or histogram_sum operation, *including* the AS histograms. It does not cover
   `DP`: the INFO DP of a spanning deletion still counts (§6).
2. If the call's ALTs are already `[*]` or `[*, <NON_REF>]`, stop here. This covers a call reduced at an earlier
   position, or an input record that was already like that.
3. Otherwise **reduce** the call, in place for the rest of the call's life in this query:
   * **Pick the `*` allele.** It is the first deletion/MNV ALT. If PL exists, it is instead the deletion/MNV ALT
     whose homozygous genotype (`[i]*ploidy`) has the strictly lowest PL; ties keep the first.
   * **Reduced alleles.** Reduced allele 0 is REF. Reduced allele 1 is the chosen ALT. Reduced allele 2 is
     `<NON_REF>`, if the call had one.
   * **Remap fields.** AD (R), PL (G) and every 2-D allele-specific INFO field are remapped to the reduced
     alleles (§5.3, `has_NON_REF` fill). All other ALTs, such as SNP/insertion alleles of the same record, are
     dropped, and their slots are filled from `<NON_REF>`.
   * **Remap GT.** Each allele maps through the reduced LUT; an unmapped allele becomes 2 if the call has
     `<NON_REF>`, otherwise missing.
   * **Rewrite REF/ALT.** REF becomes `"N"`. ALT becomes `["*", "<NON_REF>"]`, or `["*"]` without
     `<NON_REF>`.

## 5. Alleles

### 5.1 REF (`merge_reference_allele`)

Start with `mref = "N"`. Then, for every active call with `begin == s` (calls that began earlier are ignored), in
row order:

* If its REF is longer than `mref`, `mref` becomes that REF when `mref` starts with `N`. Otherwise its extra
  suffix is appended to `mref`.
* Else, if `mref` starts with `N` and this REF does not, `mref` becomes this REF.

If the final `mref` is exactly `"N"`, it is replaced by the reference base at `s`. The base is upper-cased if it is
one of `ACGTacgt`; any other base stays `N`. So:

* Blocks clipped at a window start or split mid-block get the reference base.
* `*` records inside deletions get the reference base.
* A stored REF of `N`, such as the `chr8:1` block, stays `N`.

### 5.2 ALT (`merge_alt_alleles`)

Visit **all** active calls in row order. This includes calls that began before `s`: spanning ref blocks contribute
`<NON_REF>`, reduced deletions contribute `*` and `<NON_REF>`, and other spanning records contribute their original
ALTs.

For each ALT of a call:

* `<NON_REF>` is noted and skipped.
* Any other ALT is padded first: if `len(callREF) < len(mref)` and the ALT is not symbolic, the suffix
  `mref[len(callREF):]` is appended.
* A padded ALT not seen before is appended to the merged list. Equal strings are deduplicated.
* The call's LUT records input index → merged index.

After all calls, `<NON_REF>` is appended last if any call had it.

A record is **reference-block-only** if `len(mref) == 1` and the merged ALTs are exactly `[<NON_REF>]`. Such a
record gets no remapping and no INFO DP.

### 5.3 Remapping (only when the record is not reference-block-only)

The LUT is input → merged. Its inverse, merged → input, is used for remapping; for duplicates, the last pair added
wins. In the rules below, "fill from `<NON_REF>`" means: when an allele has no input mapping and the merged record
has `<NON_REF>`, use the sample's own `<NON_REF>` entry; if the sample has no `<NON_REF>`, the value is missing.

* **GT.** Each allele maps through the LUT. An unmapped allele becomes the merged `<NON_REF>` index if the merged
  record has `<NON_REF>`, otherwise missing. Missing stays missing, and phase flags are kept.
  Ploidy is the GT length.
* **AD (R).** For each merged allele, take the input value, filled from `<NON_REF>`. An out-of-range index gives a
  missing value.
* **PL (G).** Uses ploidy 1, ploidy 2 (`bcf_alleles2gt`) or general ploidy. For each merged genotype, map both
  alleles, each filled from `<NON_REF>`; if either is unmapped the value is missing.
  * If there are more than `max_alt_alleles` merged ALTs (50 by default, counting `<NON_REF>`), PL is dropped for
    every sample. This applies to all G-length fields.
  * If `C(nalts+ploidy, ploidy) > max_genotype_count` (1024 by default), PL becomes missing for that sample only.
* **2-D allele-specific INFO fields (R).** Each merged allele gets the sample's entry for the mapped input allele.
  Missing entries are filled from `<NON_REF>` only if the field does **not** have
  `disable_remap_missing_with_non_ref`. That flag is set on `AS_RAW_MQ` and `AS_SB_TABLE`; for those, an
  unmapped allele gives an empty entry. `AS_QUALapprox`, `AS_VarDP` and the AS_RAW_*RankSum histograms are filled.
* **SB.** Number=4, so it is not remapped.

## 6. INFO

QUAL is always missing (`.`). FILTER is never produced: `.`, unfiltered. ID is the `;`-join of the distinct IDs of
the active calls; the order comes from an `unordered_set` and is not deterministic for more than one ID. With no IDs
it is `.`. HC GVCFs carry no IDs.

`END=e` is present iff `e > s`.

Only calls with a valid value contribute; invalidated spanning-deletion calls (§4) do not. The combine operation is
taken from the vidmap. A field without one gets the known-field default:

* **median:** BaseQRankSum, ClippingRankSum, MQRankSum, ReadPosRankSum, MQ, MQ0, ExcessHet.
* **sum:** RAW_MQ.
* **DP rule:** DP.

Every other INFO field has no operation and is **dropped**. In these workspaces that includes AC, AF, AN, the
`AS_*` Number=A Floats, FS, SOR, QD, InbreedingCoeff, MLEAC and MLEAF.

| operation (fields here) | result | text |
|---|---|---|
| `sum` (QUALapprox, VarDP, MQ_DP, all Integer) | Σ of first value over valid calls, as int64 | decimal integer |
| `median` (BaseQRankSum, MQRankSum, ReadPosRankSum, ExcessHet, MQ) | values sorted; element `n//2`, the upper median, no averaging | htslib `kputd` of the float32 (below) |
| `element_wise_sum`, 1-D (RAW_MQandDP, RAW_GT_COUNT) | element-wise Σ, as int64. Length is the max over calls; elements missing everywhere print `.` | `a,b[,c]`, so htsjdk gives `List<String>` |
| `element_wise_sum`, 2-D (AS_QUALapprox int, AS_VarDP int, AS_SB_TABLE int, AS_RAW_MQ float) | after §5.3 remap, per allele entry, element-wise Σ. Floats sum in **float32** in row order. Present iff some element is valid | entries joined by `|`, elements by `,`. Valid floats are `%.3f` (`0.000|5852.000|0.000`), ints plain, missing elements and empty entries are empty strings (`|247|0`, `0,0|7,1|0,0`) |
| `histogram_sum` (AS_RAW_BaseQRankSum, AS_RAW_MQRankSum, AS_RAW_ReadPosRankSum; (float bin, int count) pairs) | after remap, per allele entry, a map bin→Σcount keyed by float32 bin. Pairs with a missing or NaN bin or count are ignored. Present iff some call has the field, even with no pairs (`||`) | per entry, pairs sorted by bin ascending, printed `%.3f,count`. Entries joined by `|`, e.g. `|-1.900,1,-0.600,1|`, `||`, `|||` |
| DP (INFO) | per sample *i* over **all** samples: INFO DP if valid, else MIN_DP if valid, else FORMAT DP if valid; Σ, int64. Output iff Σ > 0 **and** the record is not reference-block-only. Spanning-deletion calls keep their INFO DP | integer |

Input texts are parsed into stored types first:

* `0.00` becomes float32 0, and the histogram bin `0.1` becomes float32 0.1, printed `0.100`.
* `NaN` histogram entries (from raw GVCFs, e.g. `|0.1,1|NaN`) contribute nothing.

**htslib `kputd`, used for median floats.** Let `d` be the float32 promoted to double.

* 0 prints `0`; negative zero prints `-0`.
* For `|d|` outside `[1e-4, 999999]`, the output is `%g`.
* Otherwise take `i = (uint64)(d*1e10)`, then add the rounding term: 5, 50, …, 5e9 for decades
  `[1e-4,1e-3)` … `[1e5,1e6)`.
  * For `d < 1`: `0.` + zero padding + the first 6 digits.
  * Otherwise keep 6 significant digits.
* Trailing zeros are culled, **but an all-zero fraction keeps one zero**. So `2` prints `2.0`, `-1.0` stays `-1.0`
  and `1.383` stays `1.383`, while zero prints `0`.

The prototype's `kputd()` is the reference. This difference matters only because the text becomes the INFO String
GATK sees.

Key order in the text is the operator's own. htsjdk keeps a map and re-sorts on writing, so order does not matter.

## 7. FORMAT (per record; samples in row order)

Each field is emitted iff at least one active call has a valid value. A sample with no call, or without the field,
is missing.

| field | value |
|---|---|
| GT | Remapped GT (§5.3), or the reduced GT for spanning deletions (§4). It is printed htslib-style: the first allele is unphased and each later allele is phased if its stored phase flag is set (`0|1`, `1|1`, `2|1`). |
| GT, no call | A sample with **no call** prints `.`, a single no-call allele, which htsjdk reads as a haploid no-call. This differs from `./.`. |
| GQ | Stored value. |
| AD, PL | Remapped. Missing prints `.`. Ref blocks without AD/PL print `.`. Raw-GVCF ref blocks have PL `0,x,y`, remapped like any sample, e.g. `0,3,42,3,42,42`. |
| SB, PGT, PID, PS, MIN_DP | Stored values, not remapped. |
| DP | FORMAT DP. Emitted if any sample has it. |

**GQ0/PL0 rule** (GATK PR 8715, only with produce_GT_field). For every sample whose `GQ == 0` **and** whose first
PL value `== 0`, GT becomes all-missing at the same ploidy (`./.`). This happens in raw-GVCF ref blocks: a block
with `GQ 0` and `PL 0,0,0` gives `./.:0:0:0:0,0,0`. Reblocked ref blocks have no PL, so the rule never fires for
them; a reblocked ref block is a called `0/0`.

**After htsjdk parsing (what GATK holds):**

* Genotype keys are written as GT, then the rest sorted alphabetically.
* Missing values (`.`) are not set on the genotype.
* Trailing missing fields disappear when htsjdk re-writes. In memory they are simply absent.

Example: a ref-block sample at a variant record is `0/0:.:2:0`, with GT/AD/DP/GQ and no AD/PL/SB.

## 8. Worked example (ws_rb chr8:282168)

| sample | its record at 282168 |
|---|---|
| HG002 | `C T,<NON_REF>` |
| HG002pcr | `C *,T,<NON_REF>` |
| NA12878 | a deletion that started earlier, reduced to `N *,<NON_REF>` |

The merged record:

* **REF/ALT:** `C` / `T,*,<NON_REF>`.
* **INFO sums:** `QUALapprox = 1294+894` (NA12878 invalidated) and `AS_QUALapprox=|2168|0|0`.
* **INFO DP:** `DP=89`, which is 29 + 20 + NA12878's INFO DP 40.
* **NA12878 genotype:** `0|2:28,0,12,0:…`, with PL `420,504,1606,0,1102,1066,504,1606,1102,1606`.
* **Next record:** 282169 starts the 1-bp `*,<NON_REF>` records. HG002pcr is `.` there because its overlapping
  282168 record replaced its own deletion.

## 9. Things that do not happen (checked)

* No QUAL, no FILTER, no ID in HC/Reblock data, and no AC/AF/AN/QD/FS/… INFO.
* No record is received twice (windows are disjoint).
* Batch layout does not matter: `ws_rb_b1` (batch size 1, 3 fragments, no consolidation) is byte-identical to
  `ws_rb`.

## 10. Not exercised by the test data (implemented from source only)

The following have no example in these workspaces:

* more than 50 merged ALTs;
* more than 1024 genotypes;
* non-diploid samples;
* MNV-only spanning records;
* multiple IDs;
* `concatenate`/`mean`/`move_to_FORMAT` operations;
* multiple arrays per contig.

The prototype does not implement general-ploidy PL remapping. The windowed queries in §11 do exercise window
starts inside deletions, overlaps and blocks.

## 11. Arbitrary-window check

`harness/GdbQuery.java` runs `GenomicsDBFeatureReader.query(contig,b,e)` for chosen intervals, the same call
FeatureDataSource makes. The intervals are in `harness/queries/*.q250`. They start inside deletions, at deletion
ends, at/after same-sample overlapping records, and at random positions. `merge_proto.py --queries` emulates the
same intervals.

| ws | queries matching exactly | records |
|---|---|---|
| ws_rb | 250/250 | 4,878/4,878 |
| ws_raw | 250/250 | 15,951/15,951 |

Two of these windows exposed the END-copy truncation rule of §3. Before it was added, the prototype disagreed
with GenomicsDB on window starts inside overlapped deletions.

Rerun with:

```
java -cp gatk.jar:harness/build ...GdbQuery $PWD/ws_raw REF chr8 harness/queries/ws_raw.q250 out.vcf
merge_proto.py ws_raw --queries harness/queries/ws_raw.q250 --queries-received out.vcf
```

## 12. Match statistics (merge_proto.py vs received/<ws>.vcf, whole line)

| ws | received records | whole-line match |
|---|---|---|
| ws_rb | 652,056 | 652,056 (100%) |
| ws_rb_b1 | 652,056 | 652,056 (100%) |
| ws_raw | 4,598,569 | 4,598,569 (100%) |
| ws11 | 242,544 | 242,544 (100%) |
