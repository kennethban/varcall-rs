#!/usr/bin/env nextflow
// FASTQ -> germline (HaplotypeCaller) and somatic (Mutect2) short variants with variant-rs, GATK Best Practices order.
// Every variant-rs step gives GATK 4.7's output; test/truth.sh runs the same chain with GATK and test/compare.sh checks.
//
//   per sample   align-rs mem: bwa-mem2 2.2.1's alignments on the GPU,  (a BAM + index per lane)
//                sorted BAM + index written directly (or bwa-mem2 | samtools sort)
//   per sample   MarkDuplicates over all lanes -> BaseRecalibrator (whole genome, all threads)
//   germline     HaplotypeCaller with the recal table applied on the fly (no recalibrated BAM) -> hard filters
//   somatic      ApplyBQSR -> Mutect2 (tumor + matched normal, or tumor only) with F1R2 counts,
//                GetPileupSummaries -> CalculateContamination, LearnReadOrientationModel -> FilterMutectCalls
//
// Heavy steps take the whole machine one at a time (multi-threaded, GPU PairHMM for the callers) rather than
// being scattered into per-interval tasks. With several GPUs, the GPU steps run side by side, one per GPU.

def sh(x) { "'${x}'" }
// GPUs: those listed in CUDA_VISIBLE_DEVICES when the pipeline is launched (e.g. CUDA_VISIBLE_DEVICES=0,1,2,3).
// Nextflow then gives each GPU task (accelerator 1) one of them, as the task's CUDA_VISIBLE_DEVICES, and runs at most
// one GPU task per GPU at a time. Unset: no GPU assignment, every task sees all GPUs and uses the first.
def gpu_ids() { (System.getenv('CUDA_VISIBLE_DEVICES') ?: '').tokenize(',') }
def gpu_task(uses_gpu) { uses_gpu && gpu_ids() ? 1 : null }
// align-rs tasks per sample: one per GPU (each with a share of the lanes; a lane's BAM is the same whichever task
// aligns it), one without GPUs or with bwa-mem2
def align_groups() { params.aligner == 'align-rs' ? Math.max(1, gpu_ids().size()) : 1 }
// The variant-rs binary's build time, written into every variant-rs command: Nextflow's task hash covers a script's text
// and the params it references, not the binary it calls (nor what a function returns), so this is what makes
// -resume rerun the variant-rs steps after a rebuild while the alignments stay cached.
params.variant_rs_build = file(params.variant_rs).lastModified()
// likewise for align-rs: a rebuild realigns (its output is bwa-mem2's either way, but a changed binary is retested)
params.align_rs_build = params.aligner == 'align-rs' ? file(params.align_rs).lastModified() : 0
// -L as the tasks get it: an interval file given relative to the launch directory made absolute (tasks run in their
// work directories); a region such as chr20:1-1000 as it is (file() would read "chr20:" as a URI scheme)
def is_ivfile(x) { x ==~ /.*\.(bed|interval_list|intervals|list)$/ }
params.intervals_arg = params.intervals ? (is_ivfile(params.intervals) ? file(params.intervals).toAbsolutePath().toString() : params.intervals) : null
// HaplotypeCaller per interval (each its own GATK run, as sarek's scatter: no waiting on one random generator across
// the genome) or one run over all intervals. auto: per interval for a short list (the 356 whole-genome calling
// regions: 672 -> 410 s), one run for a long one (295k exome targets: per target 317 s, one run 30 s; exomes have
// few positions where GATK's downsampling draws, so one run does not wait).
def count_intervals(f) {
    file(f).readLines().count { l -> l.trim() && !(['@', '#', 'track', 'browser'].any { p -> l.startsWith(p) }) }
}
params.hc_scatter = !params.intervals_arg ? false :
    params.scatter_intervals != 'auto' ? (params.scatter_intervals as boolean) :
    is_ivfile(params.intervals_arg) && count_intervals(params.intervals_arg) <= 1000
// (What the scripts use are params, not function results: Nextflow's task hash covers a script's text and the params
// it references, so -resume reruns a step when they change.)

process ALIGN {
    tag "${meta.sample} ${lanes[0]}${lanes.size() > 1 ? " +${lanes.size() - 1} lanes" : ''}"
    label 'heavy'
    accelerator gpu_task(params.aligner == 'align-rs')
    publishDir "${params.outdir}/lanes", mode: params.publish_dir_mode, enabled: params.save_lane_bams

    input:
    // a sample's lanes (samplesheet order), the number of FASTQs of each (1 single-end, 2 paired) and the FASTQs in
    // that order (each in its own directory: lanes may share file names)
    tuple val(meta), val(lanes), val(n_fq), path(reads, stageAs: 'in?/*')

    output:
    tuple val(meta), path("${meta.sample}.*.bam"), path("${meta.sample}.*.bam.bai")

    script:
    def fq = [reads].flatten()
    def first = 0
    def rows = lanes.withIndex().collect { lane, i ->
        def r = fq[first..<(first + n_fq[i])]
        first += n_fq[i]
        [bam: "${meta.sample}.${lane}.bam", rg: "@RG\\tID:${meta.sample}.${lane}\\tSM:${meta.sample}\\tPL:ILLUMINA\\tLB:${meta.sample}\\tPU:${lane}", reads: r]
    }
    def TAB = '\\t'
    def NL = '\\n'
    // align-rs (../../align-rs): bwa-mem2 mem's records exactly, seeding / extension / mate rescue on the GPU (CPU when
    // there is none), written as the coordinate-sorted BAM and .bai that samtools sort / index would give (records,
    // index bytes). One process per sample (--lanes: out.bam, @RG, FASTQs per line): the index loads once, and each
    // lane's final BAM pass (sort, write, index) overlaps the next lane's alignment; every lane's BAM is the one
    // aligning it alone gives. A 30x lane: 357 s vs 549 + 17 s. Level 1: the lane BAM is read once, by MarkDuplicates.
    def tsv = rows.collect { r -> "printf '${(['%s'] * (2 + r.reads.size())).join(TAB)}${NL}' ${r.bam} '${r.rg}' ${r.reads.join(' ')} >> lanes.tsv" }
    // samtools sort: 4 GB in all (16 x 256M). bwa-mem2 holds ~21 GB (the 17 GB index): with sort's former 8 GB the
    // two filled this 30 GB machine and each sort spill swapped the index out (a 30x lane 801 s -> 549 s with 4 GB,
    // records identical).
    def bwa = rows.collect { r ->
        "bwa-mem2 mem -t ${task.cpus} -K 100000000 -Y -R '${r.rg}' ${sh(params.fasta)} ${r.reads.join(' ')} | samtools sort -@ 16 -m 256M -l 1 -T sort.tmp -o ${r.bam} -\n" +
        "    samtools index ${r.bam}"
    }
    if (params.aligner == 'align-rs')
    """
    : > lanes.tsv
    ${tsv.join('\n    ')}
    ALIGN_RS_BUILD=${params.align_rs_build} ${sh(params.align_rs)} mem -t ${task.cpus} -K 100000000 -Y --gpu -l 1 --lanes lanes.tsv ${sh(params.fasta)}
    """
    else
    // pipefail: a failed bwa-mem2 must fail the task, not leave the sorted part of the lane as its BAM
    """
    set -o pipefail
    ${bwa.join('\n    ')}
    """
}

process MARKDUPLICATES {
    tag "${meta.sample}"
    label 'heavy'
    publishDir path: { "${params.outdir}/preprocessing/${meta.sample}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("${meta.sample}.md.bam"), path("${meta.sample}.md*.bai"), emit: bam
    path "${meta.sample}.md.metrics", emit: metrics

    script:
    def ins = bams.collect { "-I ${it}" }.join(' ')
    """
    VARIANT_RS_BUILD=${params.variant_rs_build} ${params.variant_rs} markdup ${ins} -O ${meta.sample}.md.bam -M ${meta.sample}.md.metrics \\
        --CREATE_INDEX true --tmp-dir . --threads ${task.cpus}
    """
}

process BASERECALIBRATOR {
    tag "${meta.sample}"
    label 'heavy'
    publishDir path: { "${params.outdir}/preprocessing/${meta.sample}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}.recal.table")

    script:
    def ks = params.known_sites.tokenize(',').collect { "--known-sites ${sh(it)}" }.join(' ')
    """
    VARIANT_RS_BUILD=${params.variant_rs_build} ${params.variant_rs} bqsr -I ${bam} -R ${sh(params.fasta)} ${ks} -O ${meta.sample}.recal.table --threads ${task.cpus}
    """
}

process HAPLOTYPECALLER {
    tag "${meta.sample}"
    label 'heavy'
    accelerator gpu_task(params.pairhmm != 'cpu')
    publishDir path: { "${params.outdir}/germline/${meta.sample}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai), path(table)

    output:
    tuple val(meta), path("${meta.sample}.{vcf,g.vcf}.gz"), path("${meta.sample}.{vcf,g.vcf}.gz.tbi")

    script:
    def ext = params.gvcf ? 'g.vcf.gz' : 'vcf.gz'
    // scatter_intervals: each interval is its own GATK run (as sarek's HaplotypeCaller scatter), in one process
    def opt = (params.intervals_arg ? "-L ${sh(params.intervals_arg)} " + (params.hc_scatter ? '--scatter-intervals ' : '') : '') +
        (params.dbsnp ? "--dbsnp ${sh(params.dbsnp)} " : '') + (params.gvcf ? '-ERC GVCF ' : '') +
        (params.hc_workers ? "--contig-workers ${params.hc_workers} " : '')
    """
    VARIANT_RS_BUILD=${params.variant_rs_build} ${params.variant_rs} -I ${bam} --bqsr-recal-file ${table} -R ${sh(params.fasta)} ${opt}-O ${meta.sample}.${ext} \\
        --pairhmm ${params.pairhmm} --threads ${task.cpus}
    """
}

process HARDFILTER {
    tag "${meta.sample}"
    publishDir path: { "${params.outdir}/germline/${meta.sample}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.sample}.filtered.vcf.gz"), path("${meta.sample}.filtered.vcf.gz.tbi")

    script:
    """
    VARIANT_RS_BUILD=${params.variant_rs_build} ${params.variant_rs} filter -V ${vcf} -O ${meta.sample}.filtered.vcf.gz
    """
}

process APPLYBQSR {
    tag "${meta.sample}"
    label 'heavy'
    publishDir path: { "${params.outdir}/preprocessing/${meta.sample}" }, mode: params.publish_dir_mode, enabled: params.save_recal_bam

    input:
    tuple val(meta), path(bam), path(bai), path(table)

    output:
    tuple val(meta), path("${meta.sample}.recal.bam"), path("${meta.sample}.recal*.bai")

    script:
    """
    VARIANT_RS_BUILD=${params.variant_rs_build} ${params.variant_rs} applybqsr -I ${bam} -R ${sh(params.fasta)} --bqsr-recal-file ${table} -O ${meta.sample}.recal.bam \\
        --threads ${task.cpus}
    """
}

process MUTECT2 {
    tag "${meta.id}"
    label 'heavy'
    accelerator gpu_task(params.pairhmm != 'cpu')
    publishDir path: { "${params.outdir}/somatic/${meta.id}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bams), path(bais)

    output:
    tuple val(meta), path("${meta.id}.mutect2.vcf.gz"), path("${meta.id}.mutect2.vcf.gz.tbi"), path("${meta.id}.mutect2.vcf.gz.stats"), emit: vcf
    tuple val(meta), path("${meta.id}.f1r2.tar.gz"), emit: f1r2

    script:
    def ins = bams.collect { "-I ${it}" }.join(' ')
    // -normal: the normal's SM as its BAM has it (a samplesheet sample name need not be the read groups' SM: sarek
    // writes patient_sample, the SEQC2 BAMs WES_FD_N_1)
    def normal = meta.normal ? "-normal \"\$(samtools view -H ${meta.normal}.recal.bam | grep -m1 '^@RG' | tr '\\t' '\\n' | sed -n 's/^SM://p')\" " : ''
    def opt = normal + (params.intervals_arg ? "-L ${sh(params.intervals_arg)} " : '') +
        (params.germline_resource ? "--germline-resource ${sh(params.germline_resource)} " : '') + (params.pon ? "-pon ${sh(params.pon)} " : '')
    """
    VARIANT_RS_THREADS=${task.cpus} VARIANT_RS_BUILD=${params.variant_rs_build} ${params.variant_rs} mutect2 -R ${sh(params.fasta)} ${ins} ${opt}--f1r2-tar-gz ${meta.id}.f1r2.tar.gz \\
        -O ${meta.id}.mutect2.vcf.gz --pairhmm ${params.pairhmm}
    """
}

process GETPILEUPSUMMARIES {
    tag "${meta.sample}"
    label 'heavy'
    publishDir "${params.outdir}/somatic/pileups", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}.pileups.table")

    script:
    """
    VARIANT_RS_THREADS=${task.cpus} VARIANT_RS_BUILD=${params.variant_rs_build} ${params.variant_rs} getpileupsummaries -I ${bam} -V ${sh(params.pileup_sites)} \\
        -L ${sh(params.intervals_arg && params.pileup_on_intervals ? params.intervals_arg : params.pileup_sites)} -O ${meta.sample}.pileups.table
    """
}

process CALCULATECONTAMINATION {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/somatic/${meta.id}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(tumor), path(normal)

    output:
    tuple val(meta), path("${meta.id}.contamination.table"), path("${meta.id}.segments.table")

    script:
    def matched = normal ? "-matched ${normal}" : ''
    """
    VARIANT_RS_BUILD=${params.variant_rs_build} ${params.variant_rs} calculatecontamination -I ${tumor} ${matched} -O ${meta.id}.contamination.table \\
        --tumor-segmentation ${meta.id}.segments.table
    """
}

process LEARNREADORIENTATIONMODEL {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/somatic/${meta.id}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(f1r2)

    output:
    tuple val(meta), path("${meta.id}.priors.tar.gz")

    script:
    """
    VARIANT_RS_BUILD=${params.variant_rs_build} ${params.variant_rs} learnreadorientationmodel -I ${f1r2} -O ${meta.id}.priors.tar.gz
    """
}

process FILTERMUTECTCALLS {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/somatic/${meta.id}" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats), path(contamination), path(segments), path(priors)

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi")

    script:
    """
    VARIANT_RS_BUILD=${params.variant_rs_build} ${params.variant_rs} filtermutectcalls -R ${sh(params.fasta)} -V ${vcf} --contamination-table ${contamination} \\
        --tumor-segmentation ${segments} --ob-priors ${priors} -O ${meta.id}.filtered.vcf.gz
    """
}

workflow {
    main:
    def tools = params.tools.tokenize(',')
    if (!(params.aligner in ['align-rs', 'bwa-mem2'])) {
        error "aligner: align-rs or bwa-mem2, not ${params.aligner}"
    }
    if (gpu_ids()) {
        if (gpu_ids().unique(false).size() != gpu_ids().size()) {
            error "CUDA_VISIBLE_DEVICES=${System.getenv('CUDA_VISIBLE_DEVICES')}: each GPU once"
        }
        log.info "GPUs: ${gpu_ids().join(', ')} (one GPU task per GPU)"
    } else if (System.getenv('CUDA_VISIBLE_DEVICES') == null) {
        def n = 0
        try { n = ['nvidia-smi', '-L'].execute().text.readLines().count { it.startsWith('GPU ') } } catch (Exception e) { }
        if (n > 1) {
            log.warn "${n} GPUs found, but only the first is used: launch with CUDA_VISIBLE_DEVICES=${(0..<n).join(',')} to use them all"
        }
    }
    // samplesheet: patient,sample,status,lane,fastq_1,fastq_2 (status 0 normal, 1 tumor; as nf-core/sarek; fastq_2 empty:
    // single-end)
    // --step markduplicates: start from coordinate-sorted BAMs instead (columns patient,sample,status,bam[,bai]; one
    // row per BAM, several per sample are merged by MarkDuplicates; a BAM without an index is indexed in memory)
    def rows = channel.fromPath(params.input).splitCsv(header: true)
    def lane_bams
    if (params.step == 'markduplicates') {
        lane_bams = rows.map { r ->
            def bai = r.bai ?: ([r.bam + '.bai', r.bam.replaceAll(/\.bam$/, '.bai')].find { file(it).exists() })
            tuple([patient: r.patient, sample: r.sample, status: (r.status ?: '0') as int], file(r.bam), bai ? file(bai) : null)
        }
    } else {
        // one ALIGN task per sample and GPU (align_groups), the lanes dealt out in samplesheet order
        def samples = rows
            .map { r -> tuple([patient: r.patient, sample: r.sample, status: (r.status ?: '0') as int], r.lane, r.fastq_2 ? [file(r.fastq_1), file(r.fastq_2)] : [file(r.fastq_1)]) }
            .groupTuple()
            .flatMap { meta, lanes, reads ->
                // a lane's BAM, its @RG ID and PU are named by sample and lane: two rows with the same pair would
                // write one file twice
                def dup = lanes.countBy { it }.findAll { it.value > 1 }.keySet()
                if (dup) {
                    error "samplesheet: sample ${meta.sample} has lane ${dup.join(', ')} on more than one row (sample + lane name each lane's BAM and read group)"
                }
                def n = Math.min(align_groups(), lanes.size())
                (0..<n).collect { g ->
                    def ix = (0..<lanes.size()).findAll { it % n == g }
                    tuple(meta, ix.collect { lanes[it] }, ix.collect { reads[it].size() }, ix.collect { reads[it] }.flatten())
                }
            }
        lane_bams = ALIGN(samples).flatMap { meta, bams, bais ->
            def bl = [bams].flatten()
            def il = [bais].flatten()
            bl.collect { b -> tuple(meta, b, il.find { it.name == b.name + '.bai' }) }
        }
    }
    def aligned = lane_bams
        .map { meta, bam, bai -> tuple(meta.sample, meta, bam, bai) }
        .groupTuple()
        .map { sample, metas, bams, bais -> tuple(metas[0], bams.sort { it.name }, bais.findAll { it != null }.sort { it.name }) }
    def md = MARKDUPLICATES(aligned).bam
    def tables = BASERECALIBRATOR(md)
    def md_table = md.join(tables, by: 0)

    def hc_vcf = channel.empty()
    def m2_vcf = channel.empty()
    if (tools.contains('haplotypecaller')) {
        // germline calling on normal samples (status 0), as sarek
        hc_vcf = HAPLOTYPECALLER(md_table.filter { meta, bam, bai, table -> meta.status == 0 })
        if (!params.gvcf) {
            HARDFILTER(hc_vcf)
        }
    }
    if (tools.contains('mutect2')) {
        def recal = APPLYBQSR(md_table.filter { meta, bam, bai, table -> meta.status == 1 || params.mutect2_normals })
        def by_patient = recal.map { meta, bam, bai -> tuple(meta.patient, meta, bam, bai) }.groupTuple()
        // one call per tumor: with the patient's normal when there is one, else tumor only
        def pairs = by_patient.flatMap { patient, metas, bams, bais ->
            def n = (0..<metas.size()).find { metas[it].status == 0 }
            (0..<metas.size()).findAll { metas[it].status == 1 }.collect { t ->
                def id = n == null ? metas[t].sample : "${metas[t].sample}_vs_${metas[n].sample}"
                def meta = [id: id, patient: patient, tumor: metas[t].sample, normal: n == null ? null : metas[n].sample]
                n == null ? tuple(meta, [bams[t]], [bais[t]]) : tuple(meta, [bams[t], bams[n]], [bais[t], bais[n]])
            }
        }
        def m2 = MUTECT2(pairs)
        def pileups = GETPILEUPSUMMARIES(recal)
        def pile_by_sample = pileups.map { meta, table -> tuple(meta.sample, table) }
        // combine, not join: a normal can be shared by several tumors
        def cc_in = pairs.map { meta, bams, bais -> tuple(meta.tumor, meta) }
            .combine(pile_by_sample, by: 0)
            .map { tumor, meta, ttable -> tuple(meta.normal ?: '-', meta, ttable) }
            .combine(pile_by_sample.mix(channel.of(tuple('-', []))), by: 0)
            .map { normal, meta, ttable, ntable -> tuple(meta, ttable, ntable) }
        def contamination = CALCULATECONTAMINATION(cc_in)
        def priors = LEARNREADORIENTATIONMODEL(m2.f1r2)
        m2_vcf = FILTERMUTECTCALLS(m2.vcf.join(contamination, by: 0).join(priors, by: 0))
    }
}
