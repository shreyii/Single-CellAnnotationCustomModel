# Custom Reference Model for Single-Cell RNA seq Annotation 

There are broadly two ways to annotate cell types for single-cell rna seq data, using manual annotation or reference based annotation.
For reference based annotation, selecting the relevant dataset for the study and tissue sample in question is crucial to avoid incorrect mapping. It also ensures that key cell types and cell populations of interest are accurately represented. 

This code contains creating reference models to annotate cell types using CellTypist to train Tabula Muris dataset, which is the selected reference dataset to annotate our query data. 

The model is trained for specific tissues and as a whole containing all tissues. 

The ref data and query data are in anndata formats. 

##### References: 
1. https://github.com/czbiohub-sf/tabula-muris

2. https://github.com/Teichlab/celltypist
