#!/bin/bash

# Extract gene names from gene annotations used by the 10x Genomics Cell Ranger for initial processing

gtffile="data/genome/M.musculus/refdata-gex-mm10-2020-A/genes/genes.gtf"
tsvfile="data/gene_names.tsv"

(
  echo -e "gene_id\\tgene_name"
  grep gene_name $gtffile | sed -e "s/.*gene_id\([^;]*\);.*gene_name\([^;]*\);.*/\1\\t\2/" | tr -d ' "' | sort -u
) > $tsvfile
