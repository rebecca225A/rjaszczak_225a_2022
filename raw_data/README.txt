Rebecca (I) generated a downsampled, demultiplexed Seurat object for BMS 225A.
R code used to generate these data can be provided, but here are the methods:

- Downloaded processed seurat object from: 
https://figshare.com/articles/dataset/Processed_naive_T_cell_single-cell_RNA-seq_Seurat_object/11886891
- Assigned "Sample" as active identity
- Extracted a list of "high quality" cell IDs that passed singlet checks and qc filtering,
grouped by "Sample" identity.
- Set seed and took 10% of each cell list.
 
I also wanted to include some low quality cells to practice filtering based on qc metrics.
- Downloaded raw cellranger counts from: 
https://figshare.com/articles/dataset/naive_T_cell_single-cell_RNA-seq_raw_counts_and_annotation/11894637
- Demultiplexed ALL Tcells
- Excluded doublets, but left in cells which had lower read counts and high mitochondrial c
ontent
- Generated a list of cells NOT found in the high quality cell list
- Sampled 10% of the low quality cell list.
- 10% low quality cells + 10% high quality cells lists were combined.
- The DropletUtils package used to write out barcodes.tsv, genes.tsv, matrix.mtx 

Thus our class dataset looks like a 'standard' 10X library that students might 
encounter in the future, but is a managable size for laptops with less RAM.

While multiplexing is valuable, I wanted to introduce analysis of single cell RNA seq
data without that extra variable of complexity.
