#!/usr/bin/env bash
# Re-runs sym_align on the meshes of a no_symm folder: the ones that are now
# symmetric are aligned and moved to the symm folder.
#
# usage: ./recheck_symmetry.sh <no_symm_folder> <symm_folder> [sym_align options...]
#   e.g. ./recheck_symmetry.sh remeshed_no_symm remeshed_symm -checkedge 0.5
#
#  - symmetric: the aligned mesh goes to <symm_folder>/<name>.obj and the
#    file is removed from <no_symm_folder>
#  - still not symmetric: left where it is
#  - errors, crashes, timeouts (TIMEOUT_S, default 600) and meshes above the
#    memory limit (MAXMEM_GB, default 16): left where they are, listed in
#    sym_fail.txt (folder of the executable), log in <no_symm_folder>/sym_logs/
#  - runs (cores - 2) jobs in parallel (override with JOBS=n)
#
# The executable is expected next to this script (override with
# SYM_ALIGN=/path/to/sym_align).

set -u

NOSYMDIR="${1:?usage: $0 <no_symm_folder> <symm_folder> [sym_align options...]}"
SYMDIR="${2:?usage: $0 <no_symm_folder> <symm_folder> [sym_align options...]}"
shift 2
OPTS="$*"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
EXE="${SYM_ALIGN:-$SCRIPT_DIR/sym_align}"
if [ ! -x "$EXE" ]; then
  echo "Error: executable not found: $EXE (set SYM_ALIGN)" >&2
  exit 1
fi
EXE="$(cd "$(dirname "$EXE")" && pwd)/$(basename "$EXE")"
EXE_DIR="$(dirname "$EXE")"

for D in "$NOSYMDIR" "$SYMDIR"; do
  if [ ! -d "$D" ]; then
    echo "Error: not a folder: $D" >&2
    exit 1
  fi
done
NOSYMDIR="$(cd "$NOSYMDIR" && pwd)"
SYMDIR="$(cd "$SYMDIR" && pwd)"

if [ -z "${JOBS:-}" ]; then
  NCPU="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 2)"
  JOBS=$((NCPU - 2))
fi
[ "$JOBS" -lt 1 ] && JOBS=1

TIMEOUT_S="${TIMEOUT_S:-600}"
MAXMEM_GB="${MAXMEM_GB:-16}"
POLL_S="${POLL_S:-0.2}"
MEM_POLL_S="${MEM_POLL_S:-5}"

LOGDIR="$NOSYMDIR/sym_logs"
WORKROOT="$NOSYMDIR/.symwork"
FAILFILE="$EXE_DIR/sym_fail.txt"
CSVFILE="$EXE_DIR/symmetry.csv"
mkdir -p "$LOGDIR"
rm -rf "$WORKROOT" # leftovers of an interrupted run
mkdir -p "$WORKROOT"
if [ ! -f "$CSVFILE" ]; then
  echo "mesh,symmetric,outliers,max_err,sample_score,nx,ny,nz,offset" >"$CSVFILE"
fi
touch "$FAILFILE"

# ---- memory of a process in MB (0 if it is gone)
mem_mb() {
  local M=""
  if [ "$(uname)" = "Darwin" ]; then
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
  WORK="$(mktemp -d "$WORKROOT/job.XXXXXX")"
  MESH="$WORK/$BASE"
  LOG="$LOGDIR/$NAME.log"
  ln -s "$SRC" "$MESH"

  # shellcheck disable=SC2086
  "$EXE" "$MESH" $OPTS >"$LOG" 2>&1 &
  local PID=$! T0=$SECONDS MEM_MB TMEM=$SECONDS
  KILLED=""
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
    # now symmetric: aligned mesh to symm, original out of no_symm
    mv "$OUT" "$SYMDIR/$NAME.obj"
    rm -f "$SRC" "$LOG"
    echo "$(progress)MOVED    $BASE"
  elif [ "$RC" -eq 1 ]; then
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

N0="$(find "$NOSYMDIR" -maxdepth 1 -type f \( -iname '*.obj' -o -iname '*.stl' \) | wc -l | tr -d ' ')"
NUM_TOTAL=$N0
export NUM_TOTAL
echo "Meshes to recheck: $N0   jobs: $JOBS   options: ${OPTS:-none}"

find "$NOSYMDIR" -maxdepth 1 -type f \( -iname '*.obj' -o -iname '*.stl' \) -print0 |
  xargs -0 -n 1 -P "$JOBS" bash -c 'process_one "$1"' _

rm -rf "$WORKROOT"
N1="$(find "$NOSYMDIR" -maxdepth 1 -type f \( -iname '*.obj' -o -iname '*.stl' \) | wc -l | tr -d ' ')"
echo "Finished: $((N0 - N1)) moved to $SYMDIR, $N1 left in $NOSYMDIR"
