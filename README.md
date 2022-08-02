# bootRanges.paper
Code analysis for bootranges

**script:**
1. *liver*: Derive overlap rate of  caQTLs detected in human liver tissue-SNPs associated with total cholesterol: 
   * *GRanges_construct.R*(Preprocess) -> *segBootranges.R*(perform bootRanges) -> *var.R*(Variance of the rate of overlaps and $z$ score) -> *fit_splines.R*(conditional density plot) -> *plot.R*
   * *bootstrap.R*: perform conventional shuffling.
2. *macrophage.R*: Derive overlap counts of differential expression genes and differential accessibility peaks.
And fitted peneralized splines to optimize DEG logFC thresholds. 
3. *multiomics.R*: modeling Single Cell Multiome ATAC + Gene Expression assay by 10x Genomics
4. *bootranges_vs_GSC.R*: Comparing bootRanges and GSC running time on genome-scale statistics 
