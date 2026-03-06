#!/usr/bin/env bash
# 示例：sh submit_job.sh task_1.sh "test"
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  sh submit_job.sh <task_sh> <job_name> [output_json]

Defaults:
  template: job.json
  output_json: task.json

Env vars:
  JOB_TEMPLATE   Path to template JSON (default: job.json)
  BOHR_SUBMIT_CMD Override submit command (default: "bohr job submit")
  NO_SUBMIT      If set to 1, only generate JSON and skip submit
USAGE
}

# 传入参数少于两个，打印说明书
if [ $# -lt 2 ]; then
  usage
  exit 1
fi

task_sh="$1"
job_name="$2"
output_json="${3:-}"

template_json="${JOB_TEMPLATE:-job.json}"
submit_cmd="${BOHR_SUBMIT_CMD:-bohr job submit}"

if [ ! -f "$task_sh" ]; then
  echo "Task file not found: $task_sh" >&2
  exit 1
fi

if [ ! -f "$template_json" ]; then
  echo "Template JSON not found: $template_json" >&2
  exit 1
fi

task_dir="$(cd "$(dirname "$task_sh")" && pwd)" # print the directory of task.sh
task_base="$(basename "$task_sh")"

if [ -z "$output_json" ]; then
  output_json="task.json"
  # task_stem="${task_base%.sh}"
  # output_json="${task_stem}.json"
fi

python3 - "$template_json" "$task_dir" "$task_base" "$job_name" "$output_json" <<'PY'
import json
import os
import sys
from pathlib import Path

template_path = Path(sys.argv[1])
task_dir = Path(sys.argv[2])
task_base = sys.argv[3]
job_name = sys.argv[4]
out_path = Path(sys.argv[5])
if not out_path.is_absolute():
    out_path = task_dir / out_path

with template_path.open('r', encoding='utf-8') as f:
    data = json.load(f)

# Update fields based on task file and job name
# Keep other fields unchanged from the template.
data['job_name'] = job_name
# Keep command relative to input_directory
# We set input_directory to the task directory (relative to cwd when possible).
# Use a relative path if possible to keep JSON portable.
try:
    # Force input_directory to be the task's directory relative to repo root,
    # defaulting to "jobs" for this project.
    rel_task_dir = task_dir.relative_to(Path.cwd())
    input_dir = str(rel_task_dir) if str(rel_task_dir) else 'jobs'
except Exception:
    input_dir = 'jobs'

data['input_directory'] = input_dir
# Command should be executed inside input_directory
# so we only need the script basename.
data['command'] = f"sh {task_base}"

out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open('w', encoding='utf-8') as f:
    json.dump(data, f, indent=4, ensure_ascii=False)
    f.write("\n")

print(f"Wrote {out_path}")
PY

if [ "${NO_SUBMIT:-0}" = "1" ]; then
  echo "NO_SUBMIT=1 set, skipping submit."
  exit 0
fi

# Default bohr usage in your environment:
# bohr job submit -i <job.json> -p <package_path>
if [ -n "${BOHR_SUBMIT_CMD:-}" ]; then
  $submit_cmd "$output_json"
else
  $submit_cmd -i "$output_json" -p "$task_dir"
fi
