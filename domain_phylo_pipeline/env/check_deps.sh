#!/usr/bin/env bash
# Check that the required software and Python packages are available on Ubuntu,
# and print their versions.
set -uo pipefail

ok=0
check() {
    local cmd="$1" name="$2"
    if command -v "$cmd" >/dev/null 2>&1; then
        printf "  [OK]   %-10s %s\n" "$name" "$(command -v "$cmd")"
    else
        printf "  [MISS] %-10s not installed / not active\n" "$name"
        ok=1
    fi
}

echo "== Executables =="
check hmmsearch HMMER
check hmmscan   HMMER
check hmmfetch  HMMER
check hmmpress  HMMER
check mafft     MAFFT
check trimal    trimAl
if command -v iqtree2 >/dev/null 2>&1; then
    check iqtree2 IQ-TREE
elif command -v iqtree >/dev/null 2>&1; then
    check iqtree IQ-TREE
else
    printf "  [MISS] %-10s not installed / not active\n" "IQ-TREE"
    ok=1
fi
check cd-hit    CD-HIT
check python3   Python

echo
echo "== Versions =="
hmmsearch -h 2>&1 | head -n 2 || true
mafft --version 2>&1 || true
trimal --version 2>&1 | head -n 1 || true
if command -v iqtree2 >/dev/null 2>&1; then
    iqtree2 --version 2>&1 | head -n 1 || true
elif command -v iqtree >/dev/null 2>&1; then
    iqtree --version 2>&1 | head -n 1 || true
fi
cd-hit 2>&1 | head -n 2 || true

echo
echo "== Python packages =="
python3 - <<'PY' || ok=1
mods = ["Bio", "numpy", "scipy", "pandas"]
for m in mods:
    try:
        mod = __import__(m)
        print(f"  [OK]   {m:8s} {getattr(mod, '__version__', '?')}")
    except Exception as exc:
        print(f"  [MISS] {m:8s} {exc}")
PY

exit "$ok"
