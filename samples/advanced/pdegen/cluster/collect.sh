#!/bin/bash
# Harvest results/*.txt into one table, sorted by mode, size and node count.
#
#   ./collect.sh              markdown table on stdout
#   ./collect.sh --tsv        tab separated, for gnuplot or a spreadsheet
#
# The comparison that matters is time per iteration at equal node count: CG
# against SSTEPCG. Where s-step crosses below CG is the point of the method.

fmt=${1:---md}
rows=$(
for f in results/*.txt; do
  [ -e "$f" ] || continue
  b=$(basename "$f" .txt)
  method=${b%%_*}
  nodes=$(sed -E 's/.*_n([0-9]+)_.*/\1/'  <<<"$b")
  ranks=$(sed -E 's/.*_p([0-9]+)_.*/\1/'  <<<"$b")
  mode=$(sed -E 's/.*_p[0-9]+_([a-z]+).*/\1/' <<<"$b")
  size=$(sed -E 's/.*_d([0-9]+)$/\1/;t;s/.*/-/' <<<"$b")

  get () { grep -m1 "^$1" "$f" | sed -E 's/.*:[[:space:]]*//' | tr -d ' '; }
  n=$(get "Linear system size")
  it=$(get "Iterations to convergence")
  ts=$(get "Time to solve system")
  ti=$(get "Time per iteration")
  er=$(get "Relative error estimate on exit")
  s=$(get "step size s")
  [ -n "$it" ] || continue
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$mode" "${size:--}" "$nodes" "$ranks" "$method" "${s:--}" "$it" "$ts" "$ti" "$er"
done | sort -t$'\t' -k1,1 -k2,2n -k3,3n -k5,5
)

[ -z "$rows" ] && { echo "no results found under results/"; exit 1; }

if [ "$fmt" = "--tsv" ]; then
  printf "mode\tsize\tnodes\tranks\tmethod\ts\titers\tt_solve\tt_iter\trel_err\n"
  echo "$rows"
else
  printf "| %-8s | %-9s | %5s | %5s | %-7s | %2s | %6s | %11s | %11s | %11s |\n" \
    mode size nodes ranks method s iters t_solve t_iter rel_err
  printf "|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|\n" \
    ---------- ----------- ------- ------- --------- ---- -------- ------------- ------------- -------------
  echo "$rows" | while IFS=$'\t' read -r mode size nodes ranks method s it ts ti er; do
    printf "| %-8s | %-9s | %5s | %5s | %-7s | %2s | %6s | %11s | %11s | %11s |\n" \
      "$mode" "$size" "$nodes" "$ranks" "$method" "$s" "$it" "$ts" "$ti" "$er"
  done
fi
