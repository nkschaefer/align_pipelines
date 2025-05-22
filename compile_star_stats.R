#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "USAGE: summarize_star_stats.sh <output_dir>"
    exit 1
fi

dirn=$1

find $dirn -mindepth 2 -maxdepth 2 -type f -name "Summary.csv" | head -1 | xargs cat | cut -d',' -f1 | sed 's/ /_/g' | tr '\n' '\t' | sed "s/$/\n/" | sed "s/^/libname\t/" 

find $dirn -mindepth 2 -maxdepth 2 -type f -name "Summary.csv" | while read fn; do
    libn="${fn%/Summary.csv}"
    libn="${libn##*/}"
    cat $fn | cut -d',' -f2 | tr '\n' '\t' | sed 's/$/\n/' | sed "s/^/${libn}\t/"
done
