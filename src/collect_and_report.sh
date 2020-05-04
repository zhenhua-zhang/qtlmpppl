#!/bin/bash
while read -r line; do
    zgrep "${line}" qtls.tsv.gz | gzip > ../per_traits/${line}.qtls.tsv
done < $(cut -f2 qtls.tsv.gz)
