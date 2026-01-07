#! /bin/sh

cellranger  count --id cellranger_barcode13 --transcriptome /home/jj2024/00.db/01.hg38/10.singlecell_reference/refdata-gex-GRCh38-2024-A  --feature-ref feature_reference_13.csv  --libraries libraries.csv  --create-bam false
