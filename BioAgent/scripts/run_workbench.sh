#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"

HOST="${BIOAGENT_WORKBENCH_HOST:-127.0.0.1}"
PORT="${BIOAGENT_WORKBENCH_PORT:-8013}"

python -m uvicorn bioagent.webapp.app:app --host "$HOST" --port "$PORT"
