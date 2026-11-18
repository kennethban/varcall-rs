## What and why

<!-- What does this change, and why? Link the issue it fixes. -->

## Checks

- [ ] `cargo test --release --workspace` and `tests/selftest.sh` pass (CI runs both on the CPU; with a GPU, run them locally too)
- [ ] Bug fix: a test that fails before the fix (unit test, `tests/synth` case, or harness script)
- [ ] Optimization: before / after timings, and the outputs it was compared on
- [ ] New GATK behaviour: the GATK / Picard / htsjdk source it follows (class and method), truth from the GATK 4.7.0.0 release jar
- [ ] Harness areas run, if the change touches paths the synthetic data does not reach:

<!-- See CONTRIBUTING.md. -->
