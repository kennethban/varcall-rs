#!/usr/bin/env bash
# Runs the pipeline (main.nf) with Nextflow, choosing its GPUs when CUDA_VISIBLE_DEVICES is not set: the GPUs no
# process is using (another job's, another user's), up to VARCALL_GPUS if set. CUDA_VISIBLE_DEVICES, when set (by
# you, or by a scheduler such as SLURM), is used as it is. The arguments go to `nextflow run main.nf`.
#
#   pipeline/run.sh -params-file my.yaml [--cpus N] [-resume] ...
#   VARCALL_GPUS=4 pipeline/run.sh -params-file my.yaml      # at most 4 GPUs
#   CUDA_VISIBLE_DEVICES=2,3 pipeline/run.sh -params-file my.yaml
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)

if [ -z "${CUDA_VISIBLE_DEVICES+set}" ] && command -v nvidia-smi > /dev/null; then
    max=${VARCALL_GPUS:-1000}
    # GPUs with a compute process on them (another job, another user's) are left alone
    busy=$(nvidia-smi --query-compute-apps=gpu_uuid --format=csv,noheader 2> /dev/null | sort -u || true)
    free=()
    while IFS=', ' read -r idx uuid; do
        [ -n "$idx" ] || continue
        if [ -n "$busy" ] && grep -qxF "$uuid" <<< "$busy"; then
            continue
        fi
        free+=("$idx")
    done < <(nvidia-smi --query-gpu=index,uuid --format=csv,noheader 2> /dev/null || true)
    if [ ${#free[@]} -gt 0 ]; then
        pick=("${free[@]:0:max}")
        CUDA_VISIBLE_DEVICES=$(IFS=,; echo "${pick[*]}")
        export CUDA_VISIBLE_DEVICES
        echo "[run.sh] GPUs ${CUDA_VISIBLE_DEVICES}: the idle ones${VARCALL_GPUS:+, at most VARCALL_GPUS=$VARCALL_GPUS} (CUDA_VISIBLE_DEVICES chooses them yourself)" >&2
    else
        echo "[run.sh] no idle GPU: every step uses the first one (set CUDA_VISIBLE_DEVICES to choose)" >&2
    fi
fi

exec nextflow run "$here/main.nf" "$@"
