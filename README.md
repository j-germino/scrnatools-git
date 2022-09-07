
# scrnatools

Helper methods for scRNAseq processing using scanpy + scVI

# To install from github:
```
git clone https://github.com/ntranoslab/scrnatools-git.git
cd sc-rna-tools-git
pip install -e scrnatools # allows live editing of source files
```

# To install from PyPI:
```
pip install scrnatools
```

# Contents:

## Plotting methods (rna.pl)
**gene_umap** - Plots gene expression on a UMAP with quantile thresholds for the colorbar mapping<br/>
**gene_tsne** - Plots gene expression on a TSNE with quantile thresholds for the colorbar mapping<br/>
**gene_density plot** - Plots the desity of a gene's expression on a UMAP/TSNE<br/>
**qc_plotting** - general QC plots showing read/cell, num genes/cell, percent mito reads, etc.<br/>

## QC methods (rna.qc)
**scrublet** - Uses scrublet to filter doublets from data based on a doublet score threshold<br/>
**filter_cells** - filter cells from data based on total counts, num genes, and percent mito reads thresholds<br/>

## Tools (rna.tl)
**cell_type_similarity** - Calculates cosine similarity per cell to samples in a reference dataset <br/>
**create_cell_type_signature** - Creates pseudobulk cell type signatures from single cell data <br/>
**get_immgen_similarity_signatures** - Downloads and loads immgen cell type signatures (REQUIRES WGET) <br/>
**cluster_de** - Uses scVI's differential expression method to identify marker genes for each cluster based on bayes
factor, mean lfc, and non zero expression proportion thresholds<br/>
**get_expression_matrix** - Gets a cell x gene matrix with expression values from the .X, .raw attributes or a layer of
an AnnData object<br/>
**log_density_ratio** - Calculates the relative log density of condition vs control samples in a dataset<br/>
<br/>
Detailed examples of each method and a sample analysis pipeline to be run in google colab is contained in
examples/scrnatools_example.ipynb. <br/><br/>
For a complete explanation of the scRNAseq pipelines used with these methods refer to the scanpy and scvi-tools documentation
