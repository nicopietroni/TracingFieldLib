#!/usr/bin/env bash
# Batch symmetry detection of a folder of meshes (.obj and .stl) with sym_align.
#
# usage: ./batch_symmetry.sh <input_folder> [sym_align options...]
#   e.g. ./batch_symmetry.sh remeshed -checkedge 0.5
#
#  - runs (cores - 2) jobs in parallel (override with JOBS=n)
#  - symmetric meshes: the aligned copy (symmetry plane moved onto YZ) goes
#    to <input_folder>/symm/<name>.obj
#  - non symmetric meshes: a copy of the original goes to
#    <input_folder>/no_symm/
#  - the input folder is not modified; meshes already in symm/ or no_symm/
#    are skipped, so re-running continues from where it stopped
#  - errors, crashes, timeouts (TIMEOUT_S, default 600) and meshes above
#    the memory limit (MAXMEM_GB, default 16) go to sym_fail.txt, the logs
#    to <input_folder>/sym_logs/
#  - sym_fail.txt and symmetry.csv (written by sym_align) are in the folder
#    of the executable
#
# The executable is expected next to this script (override with
# SYM_ALIGN=/path/to/sym_align).

set -u

IN="${1:?usage: $0 <input_folder> [sym_align options...]}"
shift
OPTS="$*"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EXE="${SYM_ALIGN:-$SCRIPT_DIR/sym_align}"
if [ ! -x "$EXE" ]; then
  echo "Error: executable not found: $EXE (set SYM_ALIGN)" >&2
  exit 1
fi
EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"
EXE_DIR="$(dirname "$EXE")"

if [ ! -d "$IN" ]; then
  echo "Error: not a folder: $IN" >&2
  exit 1
fi
IN="$(cd "$IN" && pwd)"

# ---- parallel jobs: all cores but 2 (at least 1)
if [ -z "${JOBS:-}" ]; then
  NCPU="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)"
  JOBS=$((NCPU - 2))
fi
[ "$JOBS" -lt 1 ] && JOBS=1

TIMEOUT_S="${TIMEOUT_S:-600}"
MAXMEM_GB="${MAXMEM_GB:-16}"
POLL_S="${POLL_S:-0.2}"
MEM_POLL_S="${MEM_POLL_S:-5}"

SYMDIR="$IN/symm"
NOSYMDIR="$IN/no_symm"
LOGDIR="$IN/sym_logs"
WORKROOT="$IN/.symwork"
FAILFILE="$EXE_DIR/sym_fail.txt"
CSVFILE="$EXE_DIR/symmetry.csv"
mkdir -p "$SYMDIR" "$NOSYMDIR" "$LOGDIR"
rm -rf "$WORKROOT" # leftovers of an interrupted run
mkdir -p "$WORKROOT"

# header written here once, so parallel jobs never race on it
# (must match the header written by main.cpp of sym_align)
if [ ! -f "$CSVFILE" ]; then
  echo "mesh,symmetric,outliers,max_err,sample_score,nx,ny,nz,offset" >"$CSVFILE"
fi
touch "$FAILFILE"

# ---- memory of a process in MB (0 if it is gone)
mem_mb() {
  local M=""
  if [ "$(uname)" = "Darwin" ]; then
    # footprint as in Activity Monitor (includes compressed memory)
    M="$(top -l 1 -pid "$1" -stats mem 2>/dev/null | tail -1 |
      awk '{v=$1; gsub(/[+-]/,"",v); u=substr(v,length(v)); n=substr(v,1,length(v)-1)+0;
            if(u=="G")n*=1024; else if(u=="K")n/=1024; else if(u=="B")n/=1048576;
            else if(u!="M")n=0; printf "%d", n}')"
  else
    M="$(ps -o rss= -p "$1" 2>/dev/null | awk '{printf "%d", $1/1024}')"
  fi
  echo "${M:-0}"
}

# ---- progress: one line per mesh in PROGFILE (appends are atomic), prints
# "[done/total] " with the elapsed time from the start
progress() {
  echo x >>"$PROGFILE"
  local K E
  K="$(wc -l <"$PROGFILE" | tr -d ' ')"
  E=$(($(date +%s) - START_T))
  printf '[%s/%s %dm%02ds] ' "$K" "$NUM_TOTAL" $((E / 60)) $((E % 60))
}

# ---- one mesh
process_one() {
  local SRC="$1"
  local BASE NAME WORK MESH LOG OUT RC KILLED
  BASE="$(basename "$SRC")"
  NAME="${BASE%.*}"

  # already classified in a previous run
  if [ -e "$SYMDIR/$NAME.obj" ] || [ -e "$NOSYMDIR/$BASE" ]; then
    progress >/dev/null
    return
  fi

  WORK="$(mktemp -d "$WORKROOT/job.XXXXXX")"
  MESH="$WORK/$BASE"
  LOG="$LOGDIR/$NAME.log"
  # linked in the job folder: sym_align writes <name>_sym.obj there
  ln -s "$SRC" "$MESH"

  # $OPTS unquoted on purpose: one word per option
  # shellcheck disable=SC2086
  "$EXE" "$MESH" $OPTS >"$LOG" 2>&1 &
  local PID=$! T0=$SECONDS MEM_MB TMEM=$SECONDS
  KILLED=""
  # SPEED: a mesh takes a few seconds, so the end is checked every 0.2 s;
  # the memory (top is slow on macOS) only every MEM_POLL_S seconds
  while kill -0 "$PID" 2>/dev/null; do
    if [ $((SECONDS - TMEM)) -ge "$MEM_POLL_S" ]; then
      TMEM=$SECONDS
      MEM_MB="$(mem_mb "$PID")"
      if [ "$MEM_MB" -gt $((MAXMEM_GB * 1024)) ]; then
        KILLED="memory limit (${MEM_MB} MB > ${MAXMEM_GB} GB)"
      fi
    fi
    if [ $((SECONDS - T0)) -gt "$TIMEOUT_S" ]; then
      KILLED="timeout (${TIMEOUT_S} s)"
    fi
    if [ -n "$KILLED" ]; then
      kill -9 "$PID" 2>/dev/null
      break
    fi
    sleep "$POLL_S"
  done
  wait "$PID" 2>/dev/null
  RC=$?
  OUT="$WORK/${NAME}_sym.obj"

  if [ -n "$KILLED" ]; then
    printf '%s\t%s\n' "$BASE" "$KILLED" >>"$FAILFILE"
    echo "$(progress)FAIL     $BASE ($KILLED)"
  elif [ "$RC" -eq 0 ] && [ -s "$OUT" ]; then
    # symmetric: the aligned mesh
    mv "$OUT" "$SYMDIR/$NAME.obj"
    rm -f "$LOG"
    echo "$(progress)SYMM     $BASE"
  elif [ "$RC" -eq 1 ]; then
    # not symmetric: a copy of the original
    cp "$SRC" "$NOSYMDIR/$BASE"
    rm -f "$LOG"
    echo "$(progress)NO_SYMM  $BASE"
  else
    local REASON
    if [ "$RC" -gt 128 ]; then
      REASON="crash (signal $((RC - 128)))"
    else
      REASON="exit code $RC: $(grep -m1 -i 'error' "$LOG" | tr '\t' ' ')"
    fi
    printf '%s\t%s\n' "$BASE" "$REASON" >>"$FAILFILE"
    echo "$(progress)FAIL     $BASE ($REASON)"
  fi
  rm -rf "$WORK"
}
export -f process_one mem_mb progress
PROGFILE="$WORKROOT/progress"
START_T=$(date +%s)
export PROGFILE START_T
export EXE OPTS TIMEOUT_S MAXMEM_GB POLL_S MEM_POLL_S SYMDIR NOSYMDIR LOGDIR WORKROOT FAILFILE

TOTAL="$(find "$IN" -maxdepth 1 -type f \( -iname '*.obj' -o -iname '*.stl' \) | wc -l | tr -d ' ')"
NUM_TOTAL=$TOTAL
export NUM_TOTAL
echo "Meshes: $TOTAL   jobs: $JOBS   options: ${OPTS:-none}   timeout: ${TIMEOUT_S}s   memory limit: ${MAXMEM_GB} GB"
echo "Output: $SYMDIR, $NOSYMDIR"
echo "Reports: $FAILFILE, $CSVFILE"

find "$IN" -maxdepth 1 -type f \( -iname '*.obj' -o -iname '*.stl' \) -print0 |
  xargs -0 -n 1 -P "$JOBS" bash -c 'process_one "$1"' _

rm -rf "$WORKROOT"
echo "Finished: $(ls "$SYMDIR" | wc -l | tr -d ' ') symmetric, $(ls "$NOSYMDIR" | wc -l | tr -d ' ') not symmetric, $(wc -l <"$FAILFILE" | tr -d ' ') lines in sym_fail.txt"
