#!/usr/bin/env bash
# Writes the test samplesheets and params files (*.csv, *.yaml) from their templates (*.in), with @DATA@ replaced by
# $VARCALL_DATA (the data directory: tests/harness/README.md). Run before nextflow run ... -params-file test/X.yaml.
cd "$(dirname "$0")" || exit 1
: "${VARCALL_DATA:?set VARCALL_DATA to the data directory}"
for f in *.csv.in *.yaml.in; do sed "s|@DATA@|$VARCALL_DATA|g" $f > ${f%.in}; done
