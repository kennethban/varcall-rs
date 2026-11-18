# GenotypeGVCFs (GATK 4.7) annotation layer — exact specification

Scope: `GenotypeGVCFs -V file.g.vcf` with by-variant traversal (no `--include-non-variant-sites`, no
`--force-output-intervals`, no `--keep-combined-raw-annotations`, no `--keep-specific-combined`, no somatic
input), annotation sets *default* (`StandardAnnotation` group) and *WARP* (`-G StandardAnnotation -G
AS_StandardAnnotation -D dbsnp`). Every input GVCF record is processed alone. The allele / genotype remapping of
the merger and `GenotypingEngine.calculateGenotypes` are assumed known and are only referenced.

Sources: `tools/walkers/{GenotypeGVCFs,GenotypeGVCFsEngine,ReferenceConfidenceVariantContextMerger}.java`,
`tools/walkers/annotator/**`, `cmdline/GATKPlugin/GATKAnnotationPluginDescriptor.java`, `utils/{Histogram,Utils,
FisherExactTest,GenotypeUtils}.java`, htsjdk 5.0.0 `VCFEncoder`, `VariantContextUtils`, `VariantContext`.

Everything in this document was checked against the four truth files (section 9); every formula listed as
"validated" reproduced 100 % of the records.

---------------------------------------------------------------------------------------------------------------

## 0. Java number formatting primitives (used everywhere below)

* `F(x, n)` = `String.format("%.<n>f", x)` (JDK 17, default locale = en/US here). Algorithm: take the *shortest
  round-trip decimal digits* of the double (`FloatingDecimal`, i.e. what `Double.toString` prints), then round
  those decimal digits **HALF_UP** at `n` places. So `F(1.005,2)="1.01"`, `F(2.675,2)="2.68"`, `F(0.125,2)="0.13"`,
  `F(-12.345,2)="-12.35"` (C/Rust would give 1.00 / 2.67 / 0.12 / -12.34). Negative zero after rounding keeps
  the sign: `F(-0.0001,2)="-0.00"`. `F(NaN,n)="NaN"`, `F(+inf,n)="Infinity"`, `F(-inf,n)="-Infinity"`.
  (`hc_genotype::javamath::format_fixed` already implements this.)
* `E3(x)` = `String.format(Locale.US,"%.3e",x)`: same rule — shortest digits, HALF_UP to 4 significant digits —
  mantissa `d.ddd`, then `e`, sign, at least two exponent digits: `-0.431 -> "-4.310e-01"`, `-1.514 ->
  "-1.514e+00"`, `-12.345 -> "-1.235e+01"`, `-10.125 -> "-1.013e+01"`.
  **Warning**: `hc_genotype::javamath::format_sci3` uses Rust's `{:.3e}` (exact-binary, half-even) and is wrong
  for 5+ significant digit ties (e.g. `-10.125` gives `-1.012e+01`). In GenotypeGVCFs this path is hit by every
  negative BaseQRankSum/MQRankSum/ReadPosRankSum (section 3.4), so fix it to the HALF_UP-on-shortest-digits rule.
  With HC inputs (`%.3f` z-scores, |z| < 10 → ≤ 4 significant digits) no rounding actually occurs, but
  |z| ≥ 10 inputs (5 significant digits) would expose it.
* `VCFDOUBLE(d)` = htsjdk `VCFEncoder.formatVCFDouble(double d)`, applied by the writer to every **`Double`**
  attribute value (INFO and FORMAT, also to each element of a `List<Double>` / `double[]`):
  ```
  if d < 1:
      if d < 0.01:            # includes ALL negative values, 0.0 and -0.0
          if |d| >= 1e-20: E3(d)
          else:            "0.00"     # 0.0 and -0.0 -> "0.00"
      else:                F(d,3)     # 0.01 <= d < 1
  else:                    F(d,2)     # d >= 1 (also +Infinity -> "Infinity")
  ```
  NaN: `NaN < 1` is false → `F(NaN,2)` = `"NaN"`.
* Values that are Java `String` are written verbatim; `Integer`/`Long` via `toString`; `Boolean true` (flags)
  as the bare key; `List`/array elements joined with `,`; `null` → `.`.

Which representation each output key has (String vs Integer vs Double) is therefore part of the spec and is
given for every key below.

---------------------------------------------------------------------------------------------------------------

## 1. Annotation set and its order

`GATKAnnotationPluginDescriptor.getResolvedInstances()` puts the instances in a `TreeSet` ordered by
`getClass().getSimpleName()` (Java `String.compareTo`, i.e. ASCII: `A..Z` < `_` < `a..z`). That order is the
order of `infoAnnotations`, which matters for (a) the RNG draw order, (b) which annotation's `combineRawData`
consumes a shared raw key, (c) `putAll` override order.

Group membership is discovered from the interfaces *directly declared* on each class (plus their
super-interfaces), not from superclasses.

**StandardAnnotation** (12 classes; GenotypeGVCFs' tool-default group, so present in both runs):

| class | kind | output key(s) | reducible (raw keys) |
|---|---|---|---|
| BaseQualityRankSumTest | INFO | BaseQRankSum | no |
| ChromosomeCounts | INFO | AN, AC, AF | no |
| Coverage | INFO | DP | no |
| DepthPerAlleleBySample | FORMAT | AD | no |
| ExcessHet | INFO | ExcessHet | **no** (its raw methods are commented-out overrides) |
| FisherStrand | INFO | FS | no |
| InbreedingCoeff | INFO | InbreedingCoeff | no |
| MappingQualityRankSumTest | INFO | MQRankSum | no |
| QualByDepth | INFO | QD | no |
| RMSMappingQuality | INFO | MQ | yes: `RAW_MQandDP` |
| ReadPosRankSumTest | INFO | ReadPosRankSum | no |
| StrandOddsRatio | INFO | SOR | no |

**AS_StandardAnnotation** (8 classes, WARP only):

| class | kind | output key | reducible (raw keys, primary first) |
|---|---|---|---|
| AS_BaseQualityRankSumTest | INFO | AS_BaseQRankSum | `AS_RAW_BaseQRankSum` |
| AS_FisherStrand | INFO | AS_FS | `AS_SB_TABLE` |
| AS_InbreedingCoeff | INFO | AS_InbreedingCoeff | no |
| AS_MappingQualityRankSumTest | INFO | AS_MQRankSum | `AS_RAW_MQRankSum` |
| AS_QualByDepth | INFO | AS_QD | `AS_QUALapprox`, `AS_QUAL`, `AS_VarDP` |
| AS_RMSMappingQuality | INFO | AS_MQ | `AS_RAW_MQ` |
| AS_ReadPosRankSumTest | INFO | AS_ReadPosRankSum | `AS_RAW_ReadPosRankSum` |
| AS_StrandOddsRatio | INFO | AS_SOR | `AS_SB_TABLE` |

StrandBiasBySample (FORMAT SB) is in neither group (it is StandardMutectAnnotation / HC-only) → `keepSB = false`.

Resolved order, default run (11 INFO + 1 FORMAT):
`BaseQualityRankSumTest, ChromosomeCounts, Coverage, [DepthPerAlleleBySample], ExcessHet, FisherStrand,
InbreedingCoeff, MappingQualityRankSumTest, QualByDepth, RMSMappingQuality, ReadPosRankSumTest, StrandOddsRatio`.

Resolved order, WARP run (19 INFO + 1 FORMAT):
`AS_BaseQualityRankSumTest, AS_FisherStrand, AS_InbreedingCoeff, AS_MappingQualityRankSumTest, AS_QualByDepth,
AS_RMSMappingQuality, AS_ReadPosRankSumTest, AS_StrandOddsRatio, BaseQualityRankSumTest, ChromosomeCounts,
Coverage, [DepthPerAlleleBySample], ExcessHet, FisherStrand, InbreedingCoeff, MappingQualityRankSumTest,
QualByDepth, RMSMappingQuality, ReadPosRankSumTest, StrandOddsRatio`
(`AS_RMSMappingQuality` < `AS_ReadPos…` because `'M'` < `'e'`; likewise `RMSMappingQuality` < `ReadPos…`).

`reducibleKeys` (= `isRequestedReducibleRawKey`) = union of `getRawKeyNames()` of the reducible INFO annotations:
* default: `{RAW_MQandDP}`
* WARP: `{AS_RAW_BaseQRankSum, AS_SB_TABLE, AS_RAW_MQRankSum, AS_QUALapprox, AS_QUAL, AS_VarDP, AS_RAW_MQ,
  AS_RAW_ReadPosRankSum, RAW_MQandDP}`

Other engine settings: `useRaw = false` (so `annotate()`, never `annotateRawData()`), `keepCombined = false`,
`rawAnnotationsToKeep = []`. `doAlleleSpecificCalcs` of the genotyping engine = "any INFO annotation is an
`AlleleSpecificAnnotation`" → **true in WARP, false in default**: only WARP gets `AS_QUAL` from genotyping.

---------------------------------------------------------------------------------------------------------------

## 2. Merger: attributes of the merged VC (`ReferenceConfidenceVariantContextMerger.merge`)

Input: one GVCF record `vc` at `loc`, `removeNonRefSymbolicAllele = true`. Let `remapped` = the record's alleles
after remapping (ref, alts, **including `<NON_REF>`** at its position), and `allelesList` = merged output alleles
(ref + real alts, `<NON_REF>` removed; `*` kept if the record has it).

### 2.1 addReferenceConfidenceAttributes (builds `annotationMap: LinkedHashMap<String, List<?>>`)

For each `(key, value)` of `vc.getAttributes()` (htsjdk lazily-decoded: a `String` if the INFO value has no
comma, else a `List<String>` of the comma-split pieces; a flag is `Boolean true`):

1. skip keys `TLOD`, `SQ`, `POPAF`, `MMQ`, `MBQ`, `MPOS`, `MFRL`.
2. if `key ∈ reducibleKeys`: `raw = join(",", getAttributeAsList(key).map(toString))` — for a single record this
   reproduces the original INFO text exactly (e.g. `AS_SB_TABLE=0,2|0,4|0,0` → `"0,2|0,4|0,0"`,
   `AS_RAW_BaseQRankSum=|-0.5,1|NaN` → `"|-0.5,1|NaN"`, `RAW_MQandDP=2322,4` → `"2322,4"`). Append
   `AlleleSpecificAnnotationData(remapped, raw)` to `annotationMap[key]`.
3. else (non-raw): create `annotationMap[key] = []` (even if nothing is added later), then with
   `s = value.toString()`:
   * if `s` has no `,`: add `parse(key, s)`;
   * else split `s` on regex `\[|,|\]`, drop empty tokens, add `parse(key, tok)` for each.
   * `parse` uses the **input header** type of `key`: `Integer` → `Integer.parseInt(tok)`; `Float` →
     `Double.parseDouble(tok)` (so `"0"` becomes `Double 0.0`, `"NaN"` becomes NaN); any other type
     (String, Flag, Character) → `NumberFormatException`. No header line: `Double` if `tok` contains `.`, else
     `Integer.parseInt`.
   * A `NumberFormatException` aborts the rest of this key (one-shot warning); values already added stay.
     Consequences: String-typed keys (`AS_RAW_*`, `AS_SB_TABLE`, `AS_QUALapprox`, `AS_VarDP` when not requested)
     and flags end with an empty list → dropped. For a multi-valued **Integer** key the list's `toString` is
     `"[a, b, c]"`, the tokens after the first carry a leading space, `Integer.parseInt(" b")` throws → only the
     first value survives (irrelevant in practice: MLEAC/AC are removed below, RAW_GT_COUNT is removed in §4).
     `Double.parseDouble` tolerates the space, so multi-valued Float keys keep all values.

### 2.2 combineAnnotations(allelesList, annotationMap) → `combined: Map<String,Object>`

```
for annot in infoAnnotations (resolved order) if annot is ReducibleAnnotation:
    for rawKey in annot.getRawKeyNames():            # primary first, then secondaries
        if annotationMap.containsKey(rawKey):
            m = annot.combineRawData(allelesList, annotationMap[rawKey])
            if m != null: combined.putAll(m)
            annotationMap.keySet().removeAll(annot.getRawKeyNames())   # ALL raw keys of this annotation
```
So yes: **after combine, `annotationMap` no longer contains any raw key of a requested reducible annotation**
(and a raw key shared by two annotations, `AS_SB_TABLE`, is combined only by the first one, AS_FisherStrand;
AS_StrandOddsRatio finds it gone and does nothing — harmless, both finalize from the same combined string).
Note the AS_QD quirk: whichever of `AS_QUALapprox`, `AS_QUAL`, `AS_VarDP` is found first is combined (always as
integers, output under the primary key `AS_QUALapprox`) and the other two are discarded. For ReblockGVCF input
(`AS_QUALapprox` and `AS_VarDP` both present) `AS_VarDP` is therefore **dropped** and never reaches AS_QD.

The combined values (all `String`) per annotation are given in section 3.

### 2.3 mergeAttributes

```
attributes = LinkedHashMap()
attributes.putAll(combined)
for key in annotationMap.keys():             # only non-raw keys remain
    values = annotationMap[key]
    if values non-empty:
        attributes[key] = values[0] if len==1 else sortedAscending(values)[len/2]   # Utils.getMedianValue: upper median, no averaging
if depth > 0: attributes["DP"] = String.valueOf(depth)          # String!
remove AC, AF, AN, MLEAC, MLEAF, END, ECNT(EVENT_COUNT_IN_HAPLOTYPE), EVENT_COUNT_IN_REGION
```
`depth = calculateVCDepth(vc)` = `vc.getAttributeAsInt("DP",0)` if the record has INFO DP, else
Σ over genotypes of (`MIN_DP` if present else FORMAT `DP` else 0).
The merged VC's attributes are then stored as `new TreeMap<>(attributes)`; ID = input IDs (non-`.`) joined by `,`,
else `.`; QUAL missing.

So the value types after merge:
* `DP`: String (e.g. `"4"`); if depth == 0 the median `Integer` of the input DP stays (such records are never
  output: `regenotypeVC` requires `DP > 0`).
* every other numeric non-raw key (`BaseQRankSum`, `MQRankSum`, `ReadPosRankSum`, `ExcessHet`, `InbreedingCoeff`,
  `QUALapprox`, `VarDP`, `RAW_GT_COUNT`, `MQ`, `FS`, `SOR`, `QD`, … if present in the input): a single `Double`
  (Float header) or `Integer` (Integer header) = the value itself for a single-valued key, or the upper median of
  the comma-separated values.
* raw keys: the combined String (section 3).

---------------------------------------------------------------------------------------------------------------

## 3. Per-annotation behaviour

Notation: `vc` in `finalizeRawData(vc, originalVC)` is `withGenotypingAnnotations` (regenotyped VC: alleles
possibly subset by genotyping, genotypes with subset AD/PL, **attributes = merged attributes + MLEAC, MLEAF
(+ AS_QUAL in WARP)**), `originalVC` is the merged VC (all `allelesList`). Finalizers read only from `vc` /
`originalVC`, never from the map being built, so finalization order affects only `putAll` overrides and the RNG.
`annotate(ref, vc, null)` is called on the trimmed, no-call-adjusted VC (§5) and its returned map is `putAll`ed
over the existing attributes; an empty map leaves any existing value (e.g. the merged median) **unchanged**.

### 3.1 RMSMappingQuality → MQ (both runs)

* combineRawData: for each input datum parse `raw.trim()`, remove `[`/`]`, split on `", *"`; must be exactly 2
  tokens → `Long` sumSq, `Long` dp; sum over inputs; result `RAW_MQandDP = String.format("%d,%d", sumSq, dp)`.
* finalizeRawData:
  * if `vc` has `RAW_MQandDP`: parse again →
    `MQ = F(Math.sqrt(sumSq / (double) dp), 2)` (String). `dp == 0` → `"Infinity"` or `"NaN"` (0/0).
  * else if `vc` has deprecated `RAW_MQ`: `UserException` unless `--allow-old-rms-mapping-quality-annotation-data`
    (then the old formula; out of scope).
  * else return empty map (MQ untouched).
  * `RAW_MQandDP` is then removed (raw key, keepCombined=false).
* annotate (likelihoods == null) → empty map.
* **MQ already present in the input**: it is a non-raw Float key → merged `Double` median. If `RAW_MQandDP` is
  also present, finalize overwrites it with the String above. If not, it survives and is written as
  `VCFDOUBLE(inputMQ)` (e.g. `60.00`, `0.500`, `5.000e-03`).
* Present iff input has `RAW_MQandDP` (or an input `MQ`).
* Validated: `MQ == F(sqrt(sumSq/dp),2)` for 48067/48067 (warp), 47574/47574 (noas), 48080/48080 (warp_rb).

### 3.2 AS_RMSMappingQuality → AS_MQ (WARP)

* combineRawData: `combined = AlleleSpecificAnnotationData(allelesList)` (every allele → null). For each input
  datum: split raw on `\|` (Java `split`, trailing empty strings removed); token `i` belongs to `remapped[i]`;
  empty token → null, else `Double.parseDouble`. Then for each allele `a` of `allelesList`: if input value
  non-null: `combined[a] = (combined[a] ?: 0) + value` (double add). `<NON_REF>`'s value is dropped (not in
  `allelesList`). Output `AS_RAW_MQ = join("|", for a in allelesList: F(combined[a] ?: 0.0, 2))`.
  **The raw sums are rounded to 2 decimals here** and re-parsed at finalize.
* finalizeRawData: parse `vc[AS_RAW_MQ]` against `originalVC.getAlleles()` (same split rules) → `perAllele`.
  `ADsum[alt]` = Σ over **all** genotypes of `vc` that have AD of `AD[i]` (i = alt index in `vc`; no GT or depth
  filter, unlike AS_QD). For each alt of `vc` in order:
  `perAllele[alt] == null ? "." : F(Math.sqrt(perAllele[alt] / (double) ADsum[alt]), 2)`, joined by `,`.
  0/0 → `"NaN"`, x/0 → `"Infinity"`. Result String `AS_MQ`; (also re-emits `AS_RAW_MQ`, removed right after).
* annotate (likelihoods null) → empty.
* Present iff input has `AS_RAW_MQ`.
* Validated 47960/47960 records (warp, no trimming) and 48065/48065 (warp_rb).

### 3.3 AS_RankSumTest family → AS_BaseQRankSum, AS_MQRankSum, AS_ReadPosRankSum (WARP)

Raw keys `AS_RAW_BaseQRankSum`, `AS_RAW_MQRankSum`, `AS_RAW_ReadPosRankSum`. HC writes them as
`|v,c,v,c…|v,c…|NaN` (ref slot empty, one `|`-field per remapped allele incl. `<NON_REF>`; an allele with no
data is empty or `NaN`).

`Histogram` (binSize = 0.1, `BIN_EPSILON = 0.01`): bin key of value `d` =
`Math.round(Math.floor((d + 0.01*0.1) / 0.1))` (long → int); a histogram is a map key → count.

* parseRawDataString(data): strip one leading `[` and trailing char if the string starts with `[`; split on
  `\|` (trailing empties removed); token `i` → allele `data.alleles[i]`; split token on `,`; for `j = 0,2,4,…`:
  if `tok[j]` non-empty: `v = Double.parseDouble(tok[j].trim())`; if `v` is not NaN and `tok[j+1]` non-empty:
  `add(v, Integer.parseInt(tok[j+1].trim()))` (count < 1 → exception). Every allele gets a (possibly empty)
  histogram.
* combineRawData: `combined` = empty histogram for every allele of `allelesList`; for each input datum parse
  (against `remapped`) and add each histogram of an allele present in `allelesList` into combined
  (count-wise merge of keys). Output String
  `for i, a in enumerate(allelesList): if a is not ref: s += "|" (i != 0 always) ; if hist non-empty: s += hist.toString()`
  where `toString` = keys ascending, `F(key*0.1, 1) + "," + count`, joined by `,`. E.g. `["C","T"]` with an empty
  alt histogram → `"|"`.
* finalizeRawData: parse `vc[raw]` against `originalVC.getAlleles()`; for every non-ref allele (a HashMap; all
  alleles have a histogram) `median()`:
  ```
  n = Σcounts; if n == 0 → null
  medianIndex = (n+1)/2 (int); counter = 0; first = null
  for key ascending:
      counter += count[key]
      if counter > medianIndex: return first==null ? key*0.1 : (first+key)/2.0*0.1
      if counter == medianIndex:
          if n odd: return key*0.1
          else first = (double) key
  return null
  ```
  (computed in this exact double arithmetic: `key * binSize` with `binSize = 0.1`, e.g. `3*0.1 =
  0.30000000000000004`). Output for each alt of `vc` in order: `median == null ? "." : F(median, 3)`, joined by
  `,` → e.g. `.`, `-0.500`, `0.000,.`. The shortcut "no results → return empty map" only triggers when
  `originalVC` has no alt allele, so the key is present for **every** record whose input has the raw key.
* annotate (likelihoods null) → empty.
* Validated 47960/47960 (warp) and 48065/48065 (warp_rb), all three keys.

### 3.4 BaseQRankSum, MQRankSum, ReadPosRankSum (non-AS RankSumTest; both runs)

`annotate(…, likelihoods = null)`: refQuals/altQuals stay empty → empty map. Hence the output value is **the
merged non-raw value: `Double.parseDouble(input text)` written with `VCFDOUBLE`**:
* `1.111 → "1.11"`, `2.000 → "2.00"`, `0.842 → "0.842"`, `0.000 → "0.00"`, `-0.000 → "0.00"`,
  `-0.431 → "-4.310e-01"`, `-1.514 → "-1.514e+00"`, `0.005 → "5.000e-03"`.
* Present iff present in the input record. (Validated: 32061/32107/32069 present and 16006/15960/15998 absent,
  all matching, warp; same for noas and warp_rb.)

### 3.5 AS_StrandBiasTest → AS_FS, AS_SOR (WARP)

Raw key `AS_SB_TABLE` = `refF,refR|alt1F,alt1R|…|nonrefF,nonrefR`.

* parseRawDataString(data): `AnnotationUtils.getAlleleLengthListOfStringFromRawData`: if it starts with `[`
  strip first/last char and all whitespace; `split("\\|", -1)` (empties kept). The number of fields must equal
  the number of alleles of `data` (else `IllegalStateException`). Each field split on `,`, non-empty tokens
  `Integer.parseInt(trim)` → list (normally `[fwd, rev]`, empty list for an empty field).
* combineRawData (done once, by AS_FisherStrand): `combined` over `allelesList`; for each input datum, for each
  allele `a` of `allelesList` with a non-null list in the datum: if `combined[a]` non-null add element-wise
  (`[0]`, `[1]`), else `combined[a] = copy` (empty list stays empty). Output
  `AS_SB_TABLE = join("|", for a in allelesList: join(",", combined[a] ?: [0,0]))`.
* finalizeRawData (both AS_FS and AS_SOR, each from `vc[AS_SB_TABLE]`, parsed against `originalVC` alleles):
  * AS_FS: for each non-ref allele `a`: alt list empty → null; else table `[[ref0, ref1], [alt0, alt1]]`
    (ref list = the ref allele's list), value = `phredScaleErrorRate(max(pValueForContingencyTable(table),
    1e-320))` (FisherStrand's p-value incl. normalization, §3.6).
  * AS_SOR: for each allele (ref included, ignored on output): empty → null; else
    `StrandOddsRatio.calculateSOR([[ref0, ref1],[alt0, alt1]])` (§3.7).
  * Output for each alt of `vc` in order: `null → "."`, else `F(value, 3)`; joined by `,`. String.
  * The raw key is re-emitted and removed.
* annotate: `StrandBiasTest.annotate` → if any genotype has `SB`, `calculateAnnotationFromGTfield`, which the AS
  classes override to return an **empty map** (AS values from finalize stay). Without SB and likelihoods null →
  empty.
* Present iff input has `AS_SB_TABLE`. AS_SOR validated 47960/47960 (warp), 48065/48065 (warp_rb).

### 3.6 FisherStrand → FS (both runs, from FORMAT SB)

`StrandBiasTest.annotate(vc)`: if `!vc.isVariant()` → empty. If **any** genotype has an `SB` attribute →
`calculateAnnotationFromGTfield(genotypes)`; else (likelihoods null) → empty.

* `getTableFromSamples(genotypes, minCount = 2)`: `sb = [0,0,0,0]`, `foundData = false`; for each genotype with
  `SB`: `foundData = true`; `data` = its 4 ints (`SB` may be a String `"a,b,c,d"`, an `ArrayList<String>` — after
  `subsetAlleleSpecificFormatFields` it is one — or `ArrayList<Integer>`); if `data[0]+data[1]+data[2]+data[3] >
  minCount` add element-wise. Return `foundData ? [[sb0,sb1],[sb2,sb3]] : null` (null → FS absent).
  SB is per sample `refFwd, refRev, altFwd, altRev` with all alts pooled; it is never re-subset when genotyping
  drops alleles.
* `pValueForContingencyTable(t)`: `sum = t00+t01+t10+t11` (addExact); if `sum > 400` (`TARGET_TABLE_SIZE*2`):
  `norm = sum / 200.0`, each cell `(int)(cell / norm)` (truncation). Then `FisherExactTest.twoSidedPValue`
  (R algorithm: `m=t00+t01, n=t10+t11, k=t00+t10, lo=max(0,k-n), hi=min(k,m)`; if `hi-lo+1 <= 1` → 1.0; else
  `logd[x] = HypergeometricDistribution(m+n, m, k).logProbability(x)` (Apache commons-math 3),
  `threshold = logd[t00] * (1 - 1e-6)`, `p = 10^log10SumLog10({logd*log10(e) : logd <= threshold})`, `min(p,1)`).
* `FS = F(phredScaleErrorRate(max(p, 1e-320)), 3)` where `phredScaleErrorRate(e) = |−10 · max(log10(e),
  log10(Double.MIN_VALUE))|` (so `-0.0` → `0.0`). String. (`hc_genotype::annotate::fisher_strand` implements all
  of this.)
* Present iff some genotype has SB. With SB sum ≤ 2 the table is all zeros → `FS=0.000`.

### 3.7 StrandOddsRatio → SOR (both runs, from FORMAT SB)

Same dispatch as FS, `getTableFromSamples(genotypes, minCount = 0)` (a sample is added iff its SB sum > 0).
`calculateSOR(t)`: `t00..t11 += 1.0`;
`ratio = (t00/t01)*(t11/t10) + (t01/t00)*(t10/t11)`; `refRatio = min(t00,t01)/max(t00,t01)`;
`altRatio = min(t10,t11)/max(t10,t11)`; `SOR = ln(ratio) + ln(refRatio) − ln(altRatio)` (`Math.log`);
output `F(SOR, 3)`, String. Null table → absent. Validated 48067/48067 (warp), 47574 (noas), 48080 (warp_rb).

### 3.8 QualByDepth → QD (both runs)

```
annotate(vc):
    if !(vc.hasLog10PError() || vc has QUALapprox): return {}
    if no genotypes: return {}
    depth = getDepth(genotypes, likelihoods=null)
    if depth == 0: return {}
    qual = -10.0 * vc.getLog10PError()          # the exact QUAL double, not the printed QUAL
    QD = fixTooHighQD(qual / depth)              # int depth promoted to double
    return {QD: F(QD, 2)}                         # String

getDepth: depth = 0; ADrestricted = 0
    for g in genotypes:
        if !(g.isHet() || g.isHomVar()): continue           # het includes 1/2; no-calls/hom-ref skipped
        if g.hasAD():
            tot = (int) Σ AD
            if tot != 0:
                if tot - AD[0] > 1: ADrestricted += tot
                depth += tot
                continue
        if g.hasDP(): depth += g.getDP()                    # (likelihoods null branch)
    return ADrestricted > 0 ? ADrestricted : depth

fixTooHighQD(QD): QD < 35 ? QD : 30 + nextGaussian() * 3          # NaN also takes the jitter branch
```
`nextGaussian` is `java.util.Random.nextGaussian()` (legacy polar Box–Muller with the cached second value,
`StrictMath.log/sqrt`) on GATK's single global generator `Utils.getRandomGenerator()` = `new Random(47382911L)`,
never reseeded in GenotypeGVCFs. First draws: `-1.5453553803199442, -0.42311771903704026, 0.3245691215126965,
-0.9206444338433659, -0.6004014126294749, -1.6650408935335503` → `25.36, 28.73, 30.97, 27.24, 28.20, 25.00`.
(`hc_genotype::javamath`'s Java Random has `next_gaussian`.)

### 3.9 AS_QualByDepth → AS_QD (WARP)

* combineRawData (on whichever of its raw keys is present first, see §2.2): `combined` over `allelesList`; for
  each datum split on `\|` (trailing empties removed); token `i` → `remapped[i]`, empty or `.` → null else
  `Integer.parseInt`; integer-add into combined; output `AS_QUALapprox = join("|", String.format("%d", v ?: 0))`.
  (Irrelevant for the output: finalize uses AS_QUAL first.)
* finalizeRawData(vc, originalVC):
  ```
  if !(vc has AS_QUAL) and !(vc has AS_QUALapprox): return null
  if no genotypes: return null
  standardDepth = originalVC has AS_VarDP ? split(AS_VarDP,"|").map(parseInt)       # never true here (§2.2)
                                          : getAlleleDepths(vc.genotypes)
  if standardDepth == null: return null          # no genotype has AD
  quals = parseQualList(vc)                       # AS_QUAL preferred
  for i in 0..len(quals)-1:
      AS_QD_i = fixTooHighQD( quals[i] / ((double) standardDepth[i+1] + (double) standardDepth[0]) )
  return {AS_QD: join(",", F(AS_QD_i, 2)), (AS_QUALapprox: join(",", quals) if vc had it)}
  ```
  * `getAlleleDepths`: length = AD length of the first genotype that has AD; all 0; for each genotype that is
    het or hom-var **and** has AD with `ΣAD − AD[0] > 1`: add its AD element-wise. (So a lone het with alt
    depth 1 gives all-zero depths.)
  * `parseQualList`: `AS_QUAL` (a `List<Integer>` put by the genotyper; must have `nAlleles−1` entries, else
    `IllegalStateException`) → `Integer.parseInt(toString)`; else `AS_QUALapprox` split `\|` with `-1`, must have
    `nAlleles` fields, alt fields parsed (empty → 0).
  * Denominator 0 → `q/0.0 = Infinity` (or NaN for 0/0) → **jitter** (both fail `< 35`).
  * `AS_QUAL` from `GenotypingEngine.composeCallAttributes` (only when doAlleleSpecificCalcs): if
    `AFresult.getAllelesUsedInGenotyping().size() > 2` (more than one alt *went into* the AF calculation, even if
    some were dropped from the output): for each non-ref allele of the output alleles,
    `(int) Math.round(-10 * AFresult.getLog10PosteriorOfAlleleAbsent(a))`; else one value
    `(int) Math.round(-10 * AFresult.log10ProbOnlyRefAlleleExists())` (≈ `round(QUAL)`; validated: for every
    biallelic-input record, non-jittered `AS_QD == F(round(QUAL)/(AD0+AD1), 2)`, 32443/32443).
* annotate (called again in annotateContext): `finalizeRawData(vc, vc)` — by then AS_QUAL / AS_QUALapprox have
  been removed (§4) → null → empty map, **no RNG draw**.
* Present whenever WARP runs genotyping (AS_QUAL always exists), even if the input has no AS data; absent only
  if no genotype has AD.

### 3.10 RNG draw order (validated on all four truth files)

Per output record, in processing order: first `finalizeAnnotations` (AS_QD, one potential draw per alt in alt
order), then `annotateContext` (QD, one potential draw). Nothing else in GenotypeGVCFs consumes the generator.
Simulating this reproduced every jittered QD of default (14466 draws), noas, warp (29271 draws incl. AS_QD) and
warp_rb with **zero mismatches**. Caveat: a record whose AS_QD was finalized but that is then not annotated/
written (properly polymorphic after genotyping but `!isPolymorphicInSamples()` or `DP ≤ 0`) still consumes its
AS_QD draws; none occurred in the truth data.

### 3.11 ChromosomeCounts → AN, AC, AF (both runs)

`VariantContextUtils.calculateChromosomeCounts(vc, new LinkedHashMap, removeStale=true, founders = ∅)` on the
no-call-adjusted VC (§5.1):
* `AN = Σ over non-filtered genotypes of #called (non-`.`) alleles` → `Integer`.
* If `AN == 0`: nothing is returned (no AN/AC/AF keys produced; existing ones would stay — they were removed by
  the merger anyway).
* Else `AN` (Integer); for each alt: `AC = Σ g.countAllele(alt)` (`Integer`, or `List<Integer>` if >1 alt),
  `AF = (double) AC / (double) AN` (`Double` or `List<Double>`).
* Written: `AN=2`, `AC=1` / `AC=1,1`, `AF=VCFDOUBLE(x)` → `0.500`, `1.00`, `0.250`, `0.333`, `0.00`.

### 3.12 Coverage → DP (both runs)

`annotate` with likelihoods null → empty. Output DP = the merger's String `depth` (input INFO DP, else Σ
MIN_DP/DP of the genotypes). Validated 100 %.

### 3.13 ExcessHet (both runs)

Not reducible; `RAW_GT_COUNT` is never finalized (just removed, §4). `annotate(vc)`:
* genotypes = all (no pedigree); if empty or `!vc.isVariant()` → empty.
* `GenotypeUtils.computeDiploidGenotypeCounts(vc, genotypes, round=true)`:
  for each g: skip unless `(hasPL && ploidy==2) || (isCalled && ploidy==2 && (hasPL || hasGQ))`;
  if no PL and hom-ref → refs += 1; else `L = normalizeFromLog10ToLinearSpace(PL/-10 vector)`; if >1 alt:
  `maxInd = argmax L`, if its allele pair has no ref allele → homs += 1, continue; else choose the het index
  `GLindex(0,a)` with the largest L over alts (start with index 1), its `GLindex(a,a)`; `L' = normalizeSumToOne
  ([L[0], L[het], L[hom]])`; then `refs += fastRound(L'0)`; `if L'0 != L'1: hets += fastRound(L'1)`;
  `if L'2 != L'1: homs += fastRound(L'2)`.
* `sampleCount = #genotypes isCalled && diploid && (hasPL || hasGQ)`; if `< 1` → empty.
* `p = exactTest(hets, refs, homs)` (Wigginton exact test, as in HC's `exact_test`); `p < 10e-60` →
  `PHRED_SCALED_MIN_P_VALUE = 160.0`; else `-10*log10(p) + 0.0`. Output `F(v, 4)` String.
* Single sample: always `0.0000`. If annotate returns empty, an input ExcessHet (merged `Double`) would survive
  as `VCFDOUBLE` (e.g. `0.00`) — does not happen for called sites.

### 3.14 InbreedingCoeff (both runs) — absent for < 10 samples

`annotate`: `genotypes.size() < 10 || !vc.isVariant()` → empty. Counts as ExcessHet (rounded);
`sampleCount = #(isCalledAndDiploidWithPLorGQ || diploidWithPL)`; `< 10` → empty.
`p = (2·refs + hets) / (2·(refs+hets+homs))`, `q = 1−p`, `F = 1 − hets / (2·p·q·sampleCount)`; output `F(F, 4)`.
If empty, an input InbreedingCoeff (median Double) would pass through via `VCFDOUBLE`.

### 3.15 AS_InbreedingCoeff (WARP) — absent for < 10 samples

`HeterozygosityCalculator(vc)`: sampleCount = #genotypes `isCalled && ploidy 2 && (hasPL || hasGQ)`; `< 10` →
empty. Otherwise (no rounding) for each counted g and each alt `a` (altIndex k = 1..): `L` = linear-normalized PLs
(or approximated from GQ: PLs `0 / GQ / round(30/-10/log10(0.5))*GQ = 10*GQ`); for `i` in `0..nAlleles-1`: if
`i == k`: `alleleCount[a] += 2·L[GL(k,k)]`; else `h = L[GL(min(i,k), max(i,k))]`: `hetCount[a] += h`,
`alleleCount[a] += h`, `alleleCount[ref] += h`; after the `i` loop `alleleCount[ref] += 2·L[0]` (note: added once
per alt). Per alt: biallelic → `refAF = refAC/(altAC+refAC)`, `altAF = 1−refAF`,
`F = 1 − het/(2·refAF·altAF·sampleCount)`; multi-allelic → `altAF = altAC/AN` (`AN = getCalledChrCount()`),
`F = 1 − het/(2·(1−altAF)·altAF·sampleCount)`. Output `join(",", F(F_a, 4))`.

### 3.16 DepthPerAlleleBySample (FORMAT AD, both runs)

`annotate(…, likelihoods = null)` returns immediately: AD is whatever genotyping produced (subset AD), plus the
cleanup fill-in of §5.3. No other FORMAT annotation runs.

### 3.17 dbSNP (`-D`; WARP): ID and DB (`VariantOverlapAnnotator`)

After all INFO annotations, on the annotated VC (trimmed alleles):
* `annotateRsID`: candidates = dbSNP records from the FeatureContext (the input record's interval) **whose start
  == vc start**; skip filtered ones; a candidate matches if any of its biallelic `Event`s
  (`splitVariantContextToEvents(…, trim=true, SET_TO_NO_CALL_NO_ANNOTATIONS, true)`) equals one of the vc's.
  `rsID = join(";", IDs of matching candidates in file order)`. If non-null: vc without ID → ID = rsID; else if
  `!vc.getID().contains(rsID)` (substring test) → `ID + ";" + rsID`.
* `annotateOverlaps`: same matching → `DB = Boolean.TRUE` (written as the bare flag `DB`).
* Same algorithm as the HC port (`variant-rs/src/dbsnp.rs`). Validated: DB present ⇔ ID ≠ `.` for all records.

---------------------------------------------------------------------------------------------------------------

## 4. finalizeAnnotations(withGenotypingAnnotations, originalVC)

```
map = LinkedHashMap(vc.attributes)          # merged attrs + MLEAC, MLEAF (+ AS_QUAL)
saved = {}                                   # rawAnnotationsToKeep is empty
for annot in infoAnnotations (resolved order) if Reducible:
    m = annot.finalizeRawData(vc, originalVC)
    if m != null: map.putAll(m)
    if !keepCombined: map.removeAll(annot.getRawKeyNames())
if !keepCombined: map.remove(AS_QUAL, QUALapprox, VarDP, RAW_GT_COUNT)
map.putAll(saved)
return VariantContextBuilder(vc).attributes(map)
```
Reducible finalizers that run: default → RMSMappingQuality only. WARP → AS_BaseQRankSum, AS_FS, AS_MQRankSum,
AS_QD, AS_MQ, AS_ReadPosRankSum, AS_SOR, MQ (in that order; AS_QD's jitter draws happen here).

Keys written / removed:

| step | adds | removes |
|---|---|---|
| AS_BaseQRankSum | AS_BaseQRankSum, (AS_RAW_BaseQRankSum) | AS_RAW_BaseQRankSum |
| AS_FS | AS_FS, (AS_SB_TABLE) | AS_SB_TABLE |
| AS_MQRankSum | AS_MQRankSum, (AS_RAW_MQRankSum) | AS_RAW_MQRankSum |
| AS_QD | AS_QD, (AS_QUALapprox if present) | AS_QUALapprox, AS_QUAL, AS_VarDP |
| AS_MQ | AS_MQ, (AS_RAW_MQ) | AS_RAW_MQ |
| AS_ReadPosRankSum | AS_ReadPosRankSum, (AS_RAW_ReadPosRankSum) | AS_RAW_ReadPosRankSum |
| AS_SOR | AS_SOR, (AS_SB_TABLE) | AS_SB_TABLE |
| MQ | MQ | RAW_MQandDP |
| manual | — | AS_QUAL, QUALapprox, VarDP, RAW_GT_COUNT |

With `keepCombined = true` (not our case) the raw keys stay (with their re-emitted, allele-subset values), the
four manual removals are skipped, and AS_QUAL / AS_QUALapprox header lines are added.

Not touched: deprecated `RAW_MQ`, `MQ_DP` (not raw keys of RMSMappingQuality) — a numeric input key of that
name would pass through as its median.

---------------------------------------------------------------------------------------------------------------

## 5. GenotypeGVCFsEngine steps around the annotations

Order inside `regenotypeVC` for a variant site with `DP > 0` (merged DP attribute parsed as int):
1. `regenotyped = calculateGenotypes(originalVC)`; null → no output.
2. `isProperlyPolymorphic(regenotyped)` (has an alt; if biallelic the alt is not `*`/`<*:DEL>` and the VC is not
   symbolic; if the first alt is `*` and the second `<NON_REF>` → false) else no output.
3. `addGenotypingAnnotations(originalVC.attributes, regenotyped)`: attrs = copy of merged attributes, then
   `MLEAC`, `MLEAF` (from the regenotyped VC; `List<Integer>` / `List<Double>`), `NDA` if present, `AS_QUAL` if
   present; all other regenotyped attributes discarded.
4. `finalizeAnnotations` (§4).
5. `reverseTrimAlleles` (§5.5).
6. `subsetAlleleSpecificFormatFields` (§5.4) with `relevantIndices[i] = originalVC.alleles.indexOf(regenotyped
   .alleles[i])`.
7. if `result.isPolymorphicInSamples()` (some genotype carries a non-ref called allele) and `DP > 0`:
   `assignNoCallsAnnotationExcludedGenotypes` (§5.1) → `annotateContext(vc, features, ref, null, all)` (§6) →
   `cleanupGenotypeAnnotations(annotated, createRefGTs=false, keepSB=false)` (§5.3). Else no output.
8. `apply`: skip if the result is spanning-deletion-only (single alt `*`) — cannot happen after step 2.

### 5.1 excludeFromAnnotations / assignNoCallsAnnotationExcludedGenotypes

`excludeFromAnnotations(g) = (g.isHomRef() || g.isNoCall()) && (!g.hasDP() || g.getDP() == 0) && g.hasGQ() &&
g.getGQ() == 0`. Applied twice: in the merger to each input genotype (before the genotype call), and again in
step 7 to the regenotyped genotypes: such a genotype gets `ploidy × NO_CALL` alleles (all other fields kept).
The change is **permanent** — the no-called genotypes are what `annotateContext` sees (AN, ExcessHet, …) and what
is written.

### 5.2 Which genotype fields exist going in

Genotype FORMAT data entering annotation: GT, AD (subset), DP, GQ, PL (subset) from genotyping, plus extended
attributes carried from the input: `MIN_DP`, `PGT`, `PID`, `PS`, `SB` (as in the input genotype).

### 5.3 cleanupGenotypeAnnotations(vc, createRefGTs, keepSB) — per genotype

```
attrs = copy(g.extendedAttributes)
depth = g.hasDP() ? g.DP : 0
if g has MIN_DP: depth = parseInt(MIN_DP); builder.DP(depth); attrs.remove(MIN_DP)
if !keepSB: attrs.remove(SB)
if g.isHomVar() && g has PGT: attrs[PGT] = "1|1"
if !g.hasAD() && vc.isVariant() && depth > 0: AD = int[nAlleles]; AD[0] = depth
if createRefGTs:                        # only --include-non-variant-sites monomorphic path
    if depth > 0 && g.hasGQ():
        alleles = GQ > 0 ? ploidy×REF : ploidy×NO_CALL
        noGQ(); attrs[RGQ] = old GQ
    else:
        alleles = ploidy×NO_CALL; noGQ(); noDP()
    noPL()
builder.noAttributes().attributes(attrs)
```
In our path (`createRefGTs = false`, `keepSB = false`): SB and MIN_DP disappear (MIN_DP → DP), PGT `0|1`/`1|0` →
`1|1` for hom-var genotypes (both `1/1` and `1|1`), missing AD filled as `[depth,0,…]`; RGQ is never written.
`isHomVar` = all alleles called and the same non-ref allele.

FORMAT column = `GT` first, then the remaining keys present in any genotype **sorted** (Java String order):
`AD, DP, GQ, PGT, PID, PL, PS` (htsjdk `calcVCFGenotypeKeys`); per sample trailing missing values are dropped.
Observed outputs: `GT:AD:DP:GQ:PL` and `GT:AD:DP:GQ:PGT:PID:PL:PS`.

### 5.4 subsetAlleleSpecificFormatFields(outputHeader, genotypes, relevantIndices)

For each genotype and each **extended** attribute key: `h = outputHeader.getFormatHeaderLine(key)` (a key
missing from the header → NPE). If `h` count type is fixed INTEGER and count == 1 → value unchanged; otherwise
`generateAnnotationValueVector(h.countType, attributeToList(value), relevantIndices)`:
`A` → `[orig[relevantIndices[i]-1] for i in 1..]`, `R` → `[orig[relevantIndices[i]] for i in 0..]` (index past
the end → filler null), anything else (G, `.`, fixed count ≠ 1) → the list itself. `attributeToList(String)` =
`new ArrayList(split(","))`. For HC GVCF inputs this only turns `SB` (Number=4) into an `ArrayList<String>` (same
text; removed later anyway). Inline AD/PL are not touched here.

### 5.5 GATKVariantContextUtils.reverseTrimAlleles(vc) = trimAlleles(vc, false, true)

```
if nAlleles <= 1 or any allele (other than '*') has length 1: return vc unchanged
seqs = bases of alleles that are neither symbolic nor '*'   (ranges [0,len) each)
# AlignmentUtils.normalizeAlleles(seqs, ranges, maxShift=0, trim=true):
endTrim = 0; minSize = min len
while minSize > 0 and all seqs share the last base of their range: shrink every range end by 1; minSize--; endTrim++
startTrim = 0
while minSize > 0 and all seqs share the first base of their range: advance every range start by 1; minSize--; startTrim++
emptyAllele = any range now empty
endClip = (emptyAllele && startTrim == 0) ? endTrim - 1 : endTrim
if endClip == 0: return vc unchanged
every non-symbolic, non-'*' allele: bases[0 .. len-endClip) (ref flag kept); '*'/symbolic unchanged
genotype alleles remapped through the same map; start unchanged; stop = start + len(newRef) - 1
```
Attributes are not touched (AS values were computed per pre-trim allele and keep their alt order). Examples
(truth): input `CAAA → C,CA`, genotyping keeps `CA` → `CAAA,CA` → `CAA,C`; `CAGAGAGAGAG → C,CAGAG,…` keeping
`C` and the long one → no trim (a length-1 allele exists); `CA → *,AA,C,CAA` keeping `*,CAA` → `CA,*,CAA` →
`C,*,CA`.

---------------------------------------------------------------------------------------------------------------

## 6. annotateContext(vc, features, ref, likelihoods = null, predicate = all)

1. Genotype annotations: for each genotype, each FORMAT annotation (`DepthPerAlleleBySample`, no-op with null
   likelihoods) → genotypes rebuilt but unchanged.
2. `infoAnnotMap = LinkedHashMap(vc.attributes)`; (no `-E` expressions); for each INFO annotation in resolved
   order: `m = annot.annotate(ref, vcWithGenotypes, null)`; if `m != null`: `putAll(m)`.
   Effective producers with null likelihoods: AS_InbreedingCoeff (≥10 samples), ChromosomeCounts (AN/AC/AF),
   ExcessHet, FisherStrand (FS from SB), InbreedingCoeff (≥10 samples), QualByDepth (QD, RNG), StrandOddsRatio
   (SOR from SB). Every other annotation returns an empty map (AS_* → empty; AS_QD → empty because AS_QUAL is
   gone; RankSumTests, RMSMappingQuality, AS_RMSMappingQuality, Coverage → empty), so the finalized / merged
   values stay.
3. dbSNP: `annotateOverlaps(annotateRsID(vc))` (§3.17).

---------------------------------------------------------------------------------------------------------------

## 7. Final INFO content per output record

htsjdk writes INFO from a `TreeMap<String,String>` → keys in Java `String.compareTo` order (ASCII: uppercase <
`_` < lowercase; a prefix sorts first). Every attribute is written (no header check failure possible here);
flags with `true` as bare key. No key is missing from the header (all are in the input header or added in §8).

Default run (`StandardAnnotation`), HC GVCF input — keys and value sources:

| key | present when | value / type |
|---|---|---|
| AC | AN > 0 | Integer / List<Integer> |
| AF | AN > 0 | VCFDOUBLE each |
| AN | AN > 0 | Integer |
| BaseQRankSum | input has it | VCFDOUBLE(input) |
| DP | always | String(depth) |
| ExcessHet | ≥1 called diploid sample with PL/GQ | `F(·,4)` |
| FS | some genotype has SB | `F(·,3)` |
| InbreedingCoeff | ≥10 samples | `F(·,4)` |
| MLEAC / MLEAF | always (genotyper) | Integer list / VCFDOUBLE list |
| MQ | input RAW_MQandDP | `F(sqrt(s/d),2)` |
| MQRankSum | input has it | VCFDOUBLE(input) |
| QD | QD depth > 0 | `F(·,2)` with jitter |
| ReadPosRankSum | input has it | VCFDOUBLE(input) |
| SOR | some genotype has SB | `F(·,3)` |
| any other numeric input key not listed in §2.3/§4 removals | input has it | Integer.toString / VCFDOUBLE of the median |

(AS_RAW_*, AS_SB_TABLE: dropped by the merger as non-numeric; AS_QUALapprox, AS_VarDP likewise; QUALapprox,
VarDP, RAW_GT_COUNT: removed in §4; MLEAC/MLEAF/AC/AF/AN/END: removed by the merger and recomputed.)

WARP adds: `AS_BaseQRankSum, AS_FS, AS_MQ, AS_MQRankSum, AS_QD, AS_ReadPosRankSum, AS_SOR` (Strings, one value
per output alt, `,`-joined), `AS_InbreedingCoeff` (≥10 samples), `DB`; ID from dbSNP. Truth key order example:
`AC;AF;AN;AS_BaseQRankSum;AS_FS;AS_MQ;AS_MQRankSum;AS_QD;AS_ReadPosRankSum;AS_SOR;BaseQRankSum;DB;DP;ExcessHet;
FS;MLEAC;MLEAF;MQ;MQRankSum;QD;ReadPosRankSum;SOR`.

ReblockGVCF input (warp_rb): `AS_QUALapprox` (AS_QD raw, combined then removed), `AS_VarDP` (dropped at combine),
`QUALapprox`, `VarDP`, `RAW_GT_COUNT` (Integer medians, removed in §4) — output key set identical to warp.

---------------------------------------------------------------------------------------------------------------

## 8. Output header (`setupVCFWriter`)

```
lines = LinkedHashSet(inputHeader.getMetaDataInInputOrder())
lines += defaultToolVCFHeaderLines            # ##GATKCommandLine=<ID=GenotypeGVCFs,...> and ##source=GenotypeGVCFs
lines.removeIf(key startsWith "GVCFBlock")    # all ##GVCFBlock… lines
lines += annotationEngine.getVCFAnnotationDescriptions(false):
          INFO annotations in resolved order: reducible → getDescriptions() only (keepCombined=false, so no raw
          lines); others getDescriptions(); then FORMAT annotations (AD); then one INFO line per overlap name
          (DB → VCFStandardHeaderLines "dbSNP Membership" when -D given)
lines += genotypingEngine.getAppropriateVCFInfoHeaders()   # NDA only with --annotate-with-num-discovered-alleles
lines += INFO MLEAC, INFO MLEAF, FORMAT RGQ, INFO DP (standard)
if keepCombined: INFO AS_QUAL, INFO AS_QUALapprox
if dbsnp: INFO DB (standard)
lines += FILTER LowQual
header = VCFHeader(lines, TreeSet(sampleNames))
```
Lines equal to an input line collapse (LinkedHashSet / TreeSet equality on the full text); a line with the same
ID but different text would appear twice. The writer emits `##fileformat=VCFv4.2` first, then
`getMetaDataInSortedOrder()` = TreeSet by `toString()` (contig lines among themselves by contig index) — i.e.
`##ALT`, `##FILTER`, `##FORMAT` (by ID), `##GATKCommandLine` (GenotypeGVCFs < HaplotypeCaller < ReblockGVCF),
`##INFO` (by ID), `##contig`, `##source`. Input lines are kept even if their key is never written (`AS_RAW_*`,
`AS_SB_TABLE`, `END`, `MIN_DP`, `SB`, `RAW_MQandDP`, `NON_REF` ALT, …).

Lines added relative to the HC GVCF header (taken verbatim from warp.vcf.gz): FORMAT RGQ; INFO AC, AF, AN, FS,
MQ, QD, SOR; WARP also INFO AS_BaseQRankSum, AS_FS, AS_MQ, AS_MQRankSum, AS_ReadPosRankSum, AS_SOR, DB (AS_QD,
AS_InbreedingCoeff, BaseQRankSum, DP, ExcessHet, InbreedingCoeff, MLEAC, MLEAF, MQRankSum, ReadPosRankSum,
AD, LowQual already existed in the input with identical text). For the reblocked input, only RGQ, DB, MLEAC,
MLEAF are new (the rest are already in the ReblockGVCF header). Copy the exact header-line texts from the truth
headers; they are `GATKVCFHeaderLines` / `VCFStandardHeaderLines` constants.

---------------------------------------------------------------------------------------------------------------

## 9. Validation (hand-worked records and whole-file checks)

Whole-file checks (Python re-implementation of the formulas above with Java-style HALF_UP formatting, against
the truth): MQ, DP, BaseQRankSum/MQRankSum/ReadPosRankSum (presence + VCFDOUBLE), SOR, AS_MQ, AS_BaseQRankSum,
AS_MQRankSum, AS_ReadPosRankSum, AS_SOR: **100 %** on warp, warp_rb, noas (AS keys on warp/warp_rb for every
untrimmed record). QD (non-jitter, within the 2-decimal QUAL rounding of the check): 100 %. AS_QD biallelic
(`round(QUAL)` as AS_QUAL): 32443/32443. RNG sequence (AS_QD alts then QD, seed 47382911): 0 mismatches over all
four files.

### chr8:60346 (warp input → warp output)
Input `A G,<NON_REF> 101.13 … AS_RAW_MQ=0.00|2322.00|0.00; AS_SB_TABLE=0,0|4,0|0,0; AS_RAW_*RankSum=||;
DP=4; ExcessHet=0.0000; RAW_MQandDP=2322,4  GT:AD:DP:GQ:PL:SB 1/1:0,4,0:4:12:…:0,0,4,0`.
* MQ = sqrt(2322/4)=24.0935 → `24.09`; AS_MQ = sqrt(2322.00/ADalt 4) → `24.09`.
* QD: hom-var, AD total 4, alt 4 > 1 → depth 4; QUAL 101.14/4 = 25.285 → `25.29`.
  AS_QD: AS_QUAL = round(101.14…) = 101; 101/(4+0) = 25.25 → `25.25`.
* SOR (SB 0,0,4,0): t=(1,1,5,1): ratio = 1·(1/5)+1·5 = 5.2; ln 5.2 + ln 1 − ln(1/5) = 1.6487+1.6094 = 3.2581 →
  `3.258`; AS_SOR from `0,0|4,0` identical `3.258`. FS: SB sum 4 > 2, table [[0,0],[4,0]], support size 1 →
  p = 1 → `0.000`; AS_FS `0.000`.
* AS_*RankSum raw `||` → combined `|` → empty histogram → median null → `.`.
* AC=2 AF=1.00 AN=2, DP=4 (String from INFO DP), ExcessHet=0.0000, DB + ID rs375050473.
Output `AC=2;AF=1.00;AN=2;AS_BaseQRankSum=.;AS_FS=0.000;AS_MQ=24.09;AS_MQRankSum=.;AS_QD=25.25;
AS_ReadPosRankSum=.;AS_SOR=3.258;DB;DP=4;ExcessHet=0.0000;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=24.09;QD=25.29;
SOR=3.258` ✓.

### chr8:60422
`RAW_MQandDP=2412,4` → MQ = sqrt(603) = 24.556 → `24.56`; `AS_RAW_MQ=0.00|1683.00|0.00`, AD 0,3 →
AS_MQ = sqrt(561) = 23.685 → `23.69` (MQ uses the raw read count 4, AS_MQ the AD count 3).
QD = 71.84/3 = 23.95 (AD-restricted depth 3, not DP 4); AS_QD = 72/3 = `24.00`. SB 0,0,1,2: t=(1,1,2,3):
ratio = 3/2+2/3 = 2.1667, ln = 0.7732; altRatio 2/3 → +0.4055 → `1.179`; FS: sum 3 > 2, [[0,0],[1,2]] → `0.000`. ✓

### chr8:61922 (first jitters)
`RAW_MQandDP=2312,2`, AD 0,2, SB 0,0,0,2, PGT 0|1, GT 1|1.
* AS_QD = 78/2 = 39 ≥ 35 → 30 + 3·g1 = 30 + 3·(−1.54536) = 25.364 → `25.36` (draw #1);
  QD = 78.32/2 = 39.16 → 30 + 3·g2 = 28.730 → `28.73` (draw #2). In the default run (no AS_QD) QD takes draw #1
  → `QD=25.36` ✓ (default.vcf.gz).
* SOR: t=(1,1,1,3): ratio = 3 + 1/3 = 3.333, ln = 1.2040; altRatio 1/3 → +1.0986 → `2.303`. FS: SB sum 2 is not
  > 2 → zero table → `0.000`. AS_SOR from `0,0|0,0` → ln 2 = `0.693` (differs from SOR!).
* FORMAT: SB removed, PGT `0|1` → `1|1` (hom-var): `1|1:0,2:2:6:1|1:61922_T_C:90,6,0:61922` ✓.
Following records 61935 / 61965 consume draws 3,4 / 5,6: AS_QD 30.97, QD 27.24 / AS_QD 28.20, QD 25.00 ✓.

### chr8:61986 (rank sums)
Input `AS_RAW_BaseQRankSum=|-0.5,1|NaN; AS_RAW_MQRankSum=|1.1,1|NaN; AS_RAW_ReadPosRankSum=|0.8,1|NaN;
BaseQRankSum=-0.431; MQRankSum=1.111; ReadPosRankSum=0.842; AS_RAW_MQ=1685.00|5445.00|0.00; RAW_MQandDP=7130,6;
AD 2,4,0; SB 0,2,0,4`.
* AS_BaseQRankSum: bin floor((−0.5+0.001)/0.1) = −5, n = 1 → −5·0.1 = −0.5 → `-0.500`; AS_MQRankSum: bin
  floor(11.01) = 11 → 1.1000000000000001 → `1.100`; AS_ReadPosRankSum → `0.800`. `NaN` fields (NON_REF) ignored.
* BaseQRankSum Double −0.431 → `-4.310e-01`; MQRankSum 1.111 → `1.11`; ReadPosRankSum 0.842 → `0.842`.
* MQ = sqrt(7130/6) = 34.472 → `34.47`; AS_MQ = sqrt(5445/4) = 36.895 → `36.90`.
* QD = 126.64/6 = `21.11`; AS_QD = 127/(4+2) = `21.17`. SOR: t=(1,3,1,5): ratio 5/3+3/5 = 2.2667 →
  0.8183 − 1.0986 + 1.6094 = `1.329`; AS_SOR same. FS/AS_FS: k = 0 → support 1 → `0.000`. ✓

### chr8:73072 (non-trivial FS)
SB 0,9,4,6 = AS_SB_TABLE `0,9|4,6` → table [[0,9],[4,6]], m=9, n=10, k=4, support 0..4; p = Σ P(x) with
P(x) ≤ P(0)·(1+1e-6) = 0.086687 → FS = −10·log10 p = 10.6204 → `10.620`, AS_FS `10.620`; SOR `0.019`.
MQ sqrt(17776/19) = `30.59`; AS_MQ sqrt(12124/10) = `34.82`; QD 260.64/19 = `13.72`; AS_QD 261/19 = `13.74`;
ReadPosRankSum −1.514 → `-1.514e+00`; MQRankSum 2.000 → `2.00`. ✓

### chr8:131677 (multi-allelic with `*`)
Input `AAAAG A,*,AAAAGAAAG,<NON_REF>; AS_RAW_MQ=0.00|12577.00|8146.00|576.00|0.00;
AS_SB_TABLE=0,0|5,1|3,2|0,1|0,0; RAW_MQandDP=21299,12; AD 0,6,5,1,0; SB 0,0,8,4`; genotyping keeps `A,*`, AD 0,6,5.
No trim (`A` has length 1). AS_MQ = sqrt(12577/6), sqrt(8146/5) = `45.78,40.36`; AS_SOR: A (5,1): ratio 2/6+6/2
= 3.333 → 1.2040+1.0986 = `2.303`; `*` (3,2): 3/4+4/3 = 2.0833 → 0.7340+0.2877 = `1.022`; AS_FS `0.000,0.000`;
AS rank sums `.,.`. MQ sqrt(21299/12) = `42.13`; QD 182.87/(0+6+5) = `16.62`; AS_QD `30.50,26.80` (AS_QUAL
183/6, 134/5 — per-allele posteriors because 3 alts entered genotyping); SOR (SB 0,0,8,4): t=(1,1,9,5) → 5/9+9/5
= 2.3556 → 0.8568+0.5878 = `1.445`; FS: m = 0 → support 1 → `0.000`; AC=1,1 AF=0.500,0.500; DB rs367782712. ✓

### chr8:288148 (genotyping drops an allele, reverse trim)
Input `CAAA C,CA,<NON_REF>`, AD 0,2,8,0, GT 1/2 → output `CAA C`, GT 1/1, AD 0,8. AS values keyed on the
pre-trim allele `CA` (index 2): AS_MQ sqrt(28800/8) = `60.00`; AS_SOR from (6,2): ratio 3/7+7/3 = 2.7619 →
1.0159+0.8473 = `1.863`; AS_QD 144/8 = `18.00`; QD 280.20/8 = 35.03 ≥ 35 → jitter `32.83`; SOR (SB 0,0,7,3)
`1.609`; MQ sqrt(57600/16) `60.00`. ✓

### chr8:1417628 (zero AD)
Output alt kept with AD 0,0, FORMAT DP 0, `RAW_MQandDP=113229,35`. QD: AD total 0 → falls to FORMAT DP 0 →
depth 0 → **QD absent**; AS_MQ = sqrt(0/0) = `NaN`; AS_QD = q/(0+0) = Infinity → jitter `32.56`; MQ = sqrt(113229
/35) = `56.88`; DP=35 (INFO DP). ✓

---------------------------------------------------------------------------------------------------------------

## 10. Uncertainties / not covered

* `--include-non-variant-sites` path (monomorphic sites: cleanup with `createRefGTs=true` → RGQ, no PL; annotate
  with RankSumTest / RMSMappingQuality / AS_RMSMappingQuality skipped) is described only in §5.3; not validated.
* Multi-sample (CombineGVCFs) inputs: InbreedingCoeff / AS_InbreedingCoeff / ExcessHet formulas are from the
  source only (truth files are single-sample, where they are absent / always 0.0000). Median rule for non-raw
  keys only matters for >1 input record — never in GenotypeGVCFs file input (one record at a time), except the
  "median of the comma-separated values of one key".
* FisherStrand p-value depends on Apache commons-math `HypergeometricDistribution.logProbability`
  (saddle-point); reuse `hc_genotype`'s port, which was validated for HaplotypeCaller.
* `%.3e` HALF_UP behaviour for 5+ significant digit inputs (|rank sum| ≥ 10) is from a JDK 17 jshell test, not
  from the truth (no such value in the data).
* JDK 17 `FloatingDecimal` is not always the *shortest* repr (pre-JDK19 `Double.toString` occasionally emits one
  extra digit); `F()` ties could in principle differ from a Ryu-based shortest-digit implementation in rare cases.
