# Contributing

Contributions are welcome, including optimizations, support for further GATK arguments or tools, support for other
GPUs, and bug reports with an input that reproduces the problem.

## Requirement: identical output

Every output must remain identical to that of the reference tool: bwa-mem2 2.2.1 for align-rs, and GATK 4.7.0.0
with Picard 3.5.0 for variant-rs. A change that improves speed but alters a single QUAL value by 0.01, or exchanges
two records with equal keys, is treated as a bug. Where a case cannot be reproduced exactly, variant-rs should
refuse it with an explanation, as it currently does for the arguments that it does not reproduce.

## Reporting a difference

If an output of align-rs or variant-rs differs from bwa-mem2 2.2.1's or GATK 4.7.0.0's for the same command, please
open an issue with the **Output differs** form. Include:

* both command lines, and the first differing records;
* an input that shows it, as small as you can make it: a region of a BAM (`samtools view -b in.bam chr1:1000-2000`),
  a few FASTQ records, or a public data set and the region;
* the varcall-rs commit, and whether the CPU path differs too (align-rs without `--gpu`, variant-rs `--pairhmm cpu`).

Other defects (crashes, errors, or a refused command that should run) should be filed with the **Bug report**
form. Security problems: see
[SECURITY.md](SECURITY.md).

## Before opening a pull request

1. Build with `environment.yml` (see [docs/building.md](docs/building.md)) and run:
   ```bash
   cargo test --release --workspace
   tests/selftest.sh
   ```
   Both must pass. `tests/selftest.sh` exits with the number of differences. CI runs both (CPU only) on every
   pull request; with a GPU, run them locally too.
2. If the change affects a code path that the synthetic data do not exercise, run the corresponding area of the
   harness (`tests/harness/<area>/cmp*.sh`; data: [tests/harness/README.md](tests/harness/README.md)) or describe
   the check that was run instead, and state this in the pull request.
3. **Bug fixes** come with a test that fails before the fix: a unit test, a case in `tests/synth`, or a harness
   script with its truth command.
4. **Optimizations** come with before / after timings: `ALIGN_RS_PROFILE=1`, `VARIANT_RS_TRACE_STAGES=1`, or the
   tool's own `time:` line. Where the fast path relies on an invariant, add a check switch that runs the plain
   version beside it ([docs/optimizations.md](docs/optimizations.md) has examples).
5. **New GATK behaviour:** name the GATK / htsjdk / Picard source it follows (class and method). Map the
   argument in `variant-rs/src/gatk_cli.rs` (`variant-rs/gatkargs/` lists every GATK argument of each tool), then
   regenerate [docs/gatk-arguments.md](docs/gatk-arguments.md) with `python3 variant-rs/gatkargs/doc.py` (CI checks
   it is current). Produce truth with the GATK 4.7.0.0 release jar.

## Code

* The code uses Rust 2021 and the release profile of the workspace. It is not formatted with rustfmt (lines are
  long), so files should not be reformatted, and a diff should contain only what the change requires.
* Comments explain the reason for the code and name the reference code being reproduced (`mem_matesw`,
  `ReadThreadingAssembler.runLocalAssembly`, ...). Quirks of the reference that are reproduced on purpose are
  marked as such.
* No new runtime dependencies without discussion. The CPU path must keep building without CUDA.
* Environment variables use the `ALIGN_RS_` / `VARIANT_RS_` prefixes.

## Where to start

* [docs/architecture.md](docs/architecture.md) maps the stages to files.
* The crate READMEs in `variant-rs/crates/` and [docs/notes/](docs/notes/) record how each GATK behaviour was
  matched and what was measured.
* Open items: testing on GPUs other than sm_120 and sm_89 and on more CPUs, more GATK arguments, one align-rs
  process on several GPUs (the pipeline spreads a sample's lanes over them instead). The speed plan and what was
  tried: [PLAN.md](PLAN.md).

## License

By contributing you agree that your contribution is licensed under the MIT License of this repository.
