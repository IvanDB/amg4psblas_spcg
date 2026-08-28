#!/bin/bash
# Harvest results/*.txt into one table, sorted by mode, size and node count.
#
#   ./collect.sh              markdown table, repetitions aggregated
#   ./collect.sh --tsv        tab separated, for gnuplot or a spreadsheet
#   ./collect.sh --raw        one row per run, nothing aggregated
#
# Past 8 nodes a single run varies by 15-20% with network contention, so the
# minimum over repetitions is the honest figure: it is the run least polluted
# by someone else's traffic. The spread is reported next to it.
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
  size=$(sed -E 's/.*_d([0-9]+)(_.*)?$/\1/;t;s/.*/-/' <<<"$b")
  swp=$(sed -E 's/.*_j([0-9]+)(_.*)?$/\1/;t;s/.*/-/' <<<"$b")

  get () { grep -m1 "^$1" "$f" | sed -E 's/.*:[[:space:]]*//' | tr -d ' '; }
  n=$(get "Linear system size")
  it=$(get "Iterations to convergence")
  ts=$(get "Time to solve system")
  ti=$(get "Time per iteration")
  er=$(get "Relative error estimate on exit")
  s=$(get "step size s")
  [ -n "$it" ] || continue
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$mode" "${size:--}" "${swp:--}" "$nodes" "$ranks" "$method" "${s:--}" "$it" "$ts" "$ti" "$er"
done | sort -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4n -k6,6
)

[ -z "$rows" ] && { echo "no results found under results/"; exit 1; }

# Aggregate repetitions of the same configuration: keep the minimum time per
# iteration and report how far the slowest run was above it.
if [ "$fmt" != "--raw" ]; then
  rows=$(echo "$rows" | awk -F'\t' '
    { k = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8
      n[k]++
      if (!(k in mn) || $10+0 < mn[k]) { mn[k] = $10+0; ts[k] = $9; er[k] = $11 }
      if (!(k in mx) || $10+0 > mx[k]) mx[k] = $10+0 }
    END { for (k in n)
            printf "%s\t%d\t%s\t%.5E\t%.1f%%\t%s\n",
                   k, n[k], ts[k], mn[k], (mx[k]/mn[k]-1)*100, er[k] }' \
    | sort -t$'\t' -k1,1 -k2,2n -k3,3n -k4,4n -k6,6)
fi

if [ "$fmt" = "--tsv" ]; then
  if [ "$fmt" = "--raw" ]; then
    printf "mode\tsize\tnodes\tranks\tmethod\ts\titers\tt_solve\tt_iter\trel_err\n"
  else
    printf "mode\tsize\tsweeps\tnodes\tranks\tmethod\ts\titers\truns\tt_solve\tt_iter_min\tspread\trel_err\n"
  fi
  echo "$rows"
elif [ "$fmt" = "--raw" ]; then
  printf "| %-8s | %-6s | %5s | %5s | %-7s | %2s | %6s | %11s | %11s | %11s |\n" \
    mode size nodes ranks method s iters t_solve t_iter rel_err
  printf "|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|\n" \
    ---------- -------- ------- ------- --------- ---- -------- ------------- ------------- -------------
  echo "$rows" | while IFS=$'\t' read -r mode size nodes ranks method s it ts ti er; do
    printf "| %-8s | %-6s | %5s | %5s | %-7s | %2s | %6s | %11s | %11s | %11s |\n" \
      "$mode" "$size" "$nodes" "$ranks" "$method" "$s" "$it" "$ts" "$ti" "$er"
  done
else
  printf "| %-8s | %-6s | %5s | %5s | %-7s | %2s | %6s | %4s | %11s | %7s | %11s |\n" \
    mode size nodes ranks method s iters runs t_iter_min spread rel_err
  printf "|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|%s|\n" \
    ---------- -------- ------- ------- --------- ---- -------- ------ ------------- --------- -------------
  echo "$rows" | while IFS=$'\t' read -r mode size nodes ranks method s it n ts ti sp er; do
    printf "| %-8s | %-6s | %5s | %5s | %-7s | %2s | %6s | %4s | %11s | %7s | %11s |\n" \
      "$mode" "$size" "$nodes" "$ranks" "$method" "$s" "$it" "$n" "$ti" "$sp" "$er"
  done
fi
