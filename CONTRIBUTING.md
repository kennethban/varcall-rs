# Contributing

Contributions are welcome: speed-ups, new GATK arguments or tools, other GPUs, bug reports with an input that
shows the difference.

## The one rule

**Outputs stay identical to the reference tool**: bwa-mem2 2.2.1 for align-rs, GATK 4.7.0.0 / Picard 3.5.0 for
variant-rs. A change that is faster but moves one QUAL by 0.01, or swaps two records with equal keys, is a bug.
If you cannot make something exact, variant-rs refuses that case with the reason, as it does today for arguments
it does not reproduce.

## Reporting a difference

If an output of align-rs or variant-rs differs from bwa-mem2 2.2.1's or GATK 4.7.0.0's for the same command, please
open an issue with the **Output differs** form. Include:

* both command lines, and the first differing records;
* an input that shows it, as small as you can make it: a region of a BAM (`samtools view -b in.bam chr1:1000-2000`),
  a few FASTQ records, or a public data set and the region;
* the varcall-rs commit, and whether the CPU path differs too (align-rs without `--gpu`, variant-rs `--pairhmm cpu`).

Other bugs (crashes, errors, a refused command that should run) go in a **Bug report**. Security problems: see
[SECURITY.md](SECURITY.md).

## Before opening a pull request

1. Build with `environment.yml` (see [docs/building.md](docs/building.md)) and run:
   ```bash
   cargo test --release --workspace
   tests/selftest.sh
   ```
   Both must pass. `tests/selftest.sh` exits with the number of differences. CI runs both (CPU only) on every
   pull request; with a GPU, run them locally too.
2. If you changed a path that the synthetic data does not exercise, run the harness area that does
   (`tests/harness/<area>/cmp*.sh`; data: [tests/harness/README.md](tests/harness/README.md)). Or describe the
   check you ran. Say which in the PR.
3. **Bug fixes** come with a test that fails before the fix: a unit test, a case in `tests/synth`, or a harness
   script with its truth command.
4. **Optimizations** come with before / after timings: `ALIGN_RS_PROFILE=1`, `VARIANT_RS_TRACE_STAGES=1`, or the
   tool's own `time:` line. Where the fast path relies on an invariant, add a check switch that runs the plain
   version beside it ([docs/optimizations.md](docs/optimizations.md) has examples).
5. **New GATK behaviour:** name the GATK / htsjdk / Picard source it follows (class and method). Map the
   argument in `variant-rs/src/gatk_cli.rs` (`variant-rs/gatkargs/` lists every GATK argument of each tool). Produce truth with the
   GATK 4.7.0.0 release jar.

## Code

* Rust 2021, the workspace's release profile. The code is not formatted with rustfmt (long lines), so don't
  reformat files you aren't otherwise changing; keep diffs to what the change needs.
* Comments say *why*, and name the reference code being reproduced (`mem_matesw`,
  `ReadThreadingAssembler.runLocalAssembly`, ...). Quirks of the reference that are reproduced on purpose are
  marked as such.
* No new runtime dependencies without discussion. The CPU path must keep building without CUDA.
* Environment variables use the `ALIGN_RS_` / `VARIANT_RS_` prefixes.

## Where to start

* [docs/architecture.md](docs/architecture.md) maps the stages to files.
* The crate READMEs in `variant-rs/crates/` and [docs/notes/](docs/notes/) record how each GATK behaviour was
  matched and what was measured.
* Open items: testing on GPUs other than sm_120 and on more CPUs, more GATK arguments, multi-GPU.

## License

By contributing you agree that your contribution is licensed under the MIT License of this repository.
