#!/bin/bash
# Pair CG against s-step at every configuration and derive what a scaling
# study is actually read from. collect.sh lists runs; this interprets them.
#
#   ./scaling.sh [mode]        mode defaults to converge
#
# Columns:
#   t_cg, t_ss   best time per iteration over repetitions, seconds
#   ss/cg        the headline number. Below 1 means s-step wins.
#   eff          t(first point) / t(this point), per method.
#                Weak scaling ideal is 1.00 flat: the local problem is fixed,
#                so time should not move. Strong scaling ideal is the node
#                ratio, printed alongside.
#
# Weak and strong are told apart by whether IDIM moves with the node count.

want=${1:-converge}

# Locate the columns by NAME rather than by position: an older collect.sh
# emits ten fields instead of twelve, and a fixed index silently reads the
# residual where the timing should be.
tsv=$(bash "$(dirname "$0")/collect.sh" --tsv) || exit 1
hdr=$(head -1 <<<"$tsv")
col () { awk -F'\t' -v n="$1" '{for(i=1;i<=NF;i++) if($i==n){print i; exit}}' <<<"$hdr"; }
c_mode=$(col mode); c_size=$(col size); c_nodes=$(col nodes)
c_ranks=$(col ranks); c_meth=$(col method)
c_t=$(col t_iter_min); [ -z "$c_t" ] && c_t=$(col t_iter)
if [ -z "$c_t" ] || [ -z "$c_nodes" ]; then
  echo "collect.sh gave columns this script does not recognise:"; echo "  $hdr"
  echo "The two are out of step; update both from the repository."; exit 1
fi

rows=$(tail -n +2 <<<"$tsv" \
       | awk -F'\t' -v w="$want" -v m=$c_mode -v n=$c_nodes -v r=$c_ranks \
             -v z=$c_size -v e=$c_meth -v t=$c_t \
             '$m==w {print $n"\t"$r"\t"$z"\t"$e"\t"$t}')
[ -z "$rows" ] && { echo "no results for mode '$want'"; exit 1; }

# Two campaigns in results/ get silently merged: this script keys on the node
# count, so a run at 64M and a weak run at the same node count overwrite each
# other and produce a row that is half one and half the other. Refuse instead.
mixed=$(cut -f3 <<<"$rows" | sort -u | tr '\n' ' ')
if grep -q -- '-' <<<"$mixed" && grep -qE '[0-9]' <<<"$mixed"; then
  echo "ERROR: results/ holds runs with a fixed size and runs with IDIM set:"
  echo "         sizes present: $mixed"
  echo "       These are different campaigns and cannot go in one table."
  echo "       Move the ones without _d in the name aside, for example:"
  echo "         mkdir -p results/old"
  echo "         for f in results/*.txt; do case \"\$f\" in *_d*) ;; *) mv \"\$f\" results/old/ ;; esac; done"
  exit 1
fi

echo "$rows" | sort -k1,1n -k4,4 | awk -F'\t' '
  { nodes=$1; ranks=$2; size=$3; meth=$4; t=$5+0
    key=nodes; N[key]=nodes; R[key]=ranks; S[key]=size
    if (meth=="cg") cg[key]=t; else ss[key]=t }
  END {
    n=0; for (k in N) ord[++n]=k+0
    for (i=1;i<n;i++) for (j=i+1;j<=n;j++) if (ord[i]>ord[j]) { t=ord[i]; ord[i]=ord[j]; ord[j]=t }

    sizes=0; first=ord[1]
    for (i=1;i<=n;i++) if (S[ord[i]] != S[first]) sizes=1
    kind = sizes ? "WEAK" : "STRONG"

    printf "%s scaling, mode %s\n\n", kind, "'"$want"'"
    printf "| %5s | %5s | %6s | %10s | %11s | %11s | %6s | %7s | %7s |\n",
           "nodes","ranks","IDIM","unk/rank","t_cg","t_ss","ss/cg","eff_cg","eff_ss"
    printf "|%s|%s|%s|%s|%s|%s|%s|%s|%s|\n",
           "-------","-------","--------","------------","-------------","-------------","--------","---------","---------"
    for (i=1;i<=n;i++) {
      k=ord[i]
      d = (S[k]=="-") ? 400 : S[k]+0
      unk = d*d*d/R[k]
      if (i==1) { b_cg=cg[k]; b_ss=ss[k] }
      printf "| %5d | %5d | %6s | %10d | %11.5E | %11.5E | %6.3f | %7.2f | %7.2f |\n",
             N[k], R[k], S[k], unk, cg[k], ss[k], ss[k]/cg[k], b_cg/cg[k], b_ss/ss[k]
    }
    printf "\nideal eff: "
    if (kind=="WEAK") printf "1.00 flat, the local problem does not change\n"
    else { printf "the node ratio ->"; for (i=1;i<=n;i++) printf " %dx", N[ord[i]]/N[first]; printf "\n" }
  }'
