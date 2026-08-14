#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
  echo "Usage: $0 <output-dir> [sample-json] [era]" >&2
  exit 1
fi

OUTPUT_DIR="$1"
SAMPLE_JSON="${2:-source_cleanup/json/samples/SR_medium_muon/RunIII2024Summer24NanoAODv15_corrected_temp.json}"
ERA="${3:-RunIII2024Summer24NanoAODv15}"

CHANNELS=(
  TTHH_DL
  TTHH_SL
  ttbarDL
  ttbarSL
  TTBB_DL
  TTBB_SL
  TTHBB
  TTW
  TTZ_high
  TTZ_low
  TTHnonBB
  TTTT
)

mkdir -p "${OUTPUT_DIR}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 "${SCRIPT_DIR}/extractSignalFromSampleJson.py" \
  --sampleJson "${SAMPLE_JSON}" \
  --output-dir "${OUTPUT_DIR}" \
  --era "${ERA}" \
  --channels "${CHANNELS[@]}"
