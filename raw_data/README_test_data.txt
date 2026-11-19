Bulk RNA-seq test data (mouse, GRCm39 / Ensembl 113)

Public data from DDBJ BioProject PRJDB17560 (downloaded from ENA):
"Sex differences in the function of hepatic VLDLR during protein deficiency in mice".
Mouse liver, Illumina NovaSeq 6000, paired-end. Each run is downsampled to 500,000
read pairs with seqtk sample -s11 (same seed for R1 and R2). Rebuild with:
    bash raw_data/make_testdata.sh

Files are renamed to SLX-11111.<Barcode>.AAA111.r_{1,2}.fq.gz. Pool, Barcode and index
Sequence in inputs/metadata.csv are placeholders, and Batch (1/2) is ARTIFICIAL (added only
to exercise the numeric PCA/correlation plots). Sex, diet and genotype are the real labels.

Barcode  Run        Sex  Condition 1 (diet)  Condition 2 (genotype)
UDI001   DRR530875  F    control (20% casein)  WT
UDI002   DRR530877  F    low_protein (2%)      WT
UDI003   DRR530879  F    control               Vldlr_KO
UDI004   DRR530881  F    low_protein           Vldlr_KO
UDI005   DRR530874  M    control               WT
UDI006   DRR530876  M    low_protein           WT
UDI007   DRR530878  M    control               Vldlr_KO
UDI008   DRR530880  M    low_protein           Vldlr_KO
