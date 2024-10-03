# Custom Reference Model for Single-Cell RNA seq Annotation 

There are broadly two ways to annotate cell types for single-cell rna seq data, using manual annotation or reference based annotation.
For reference based annotation, selecting the relevant dataset for the study and tissue sample in question is crucial to avoid incorrect mapping. It also ensures that key cell types and cell populations of interest are accurately represented. 

This code contains creating reference models to annotate cell types using CellTypist to train Tabula Muris dataset which is the selected reference dataset to annotate our query data for mouse. 

The model is trained for specific tissues and as a whole containing all tissues. 

The ref data and query data are in anndata formats. 

## References 
1. The Tabula Muris Consortium., Overall coordination., Logistical coordination. et al. Single-cell transcriptomics of 20 mouse organs creates a Tabula Muris. Nature 562, 367–372 (2018). https://doi.org/10.1038/s41586-018-0590-4

2. Domínguez Conde, C., Xu, C., Jarvis, L. B., Rainbow, D. B., Wells, S. B., Gomes, T., Howlett, S. K., Suchanek, O., Polanski, K., King, H. W., Mamanova, L., Huang, N., Szabo, P. A., Richardson, L., Bolt, L., Fasouli, E. S., Mahbubani, K. T., Prete, M., Tuck, L., Richoz, N., Tuong, Z. K., Campos, L., Mousa, H. S., Needham, E. J., Pritchard, S., Li, T., Elmentaite, R., Park, J., Rahmani, E., Chen, D., Menon, D. K., Bayraktar, O. A., James, L. K., Meyer, K. B., Yosef, N., Clatworthy, M. R., Sims, P. A., Farber, D. L., Saeb-Parsy, K., Jones, J. L., & Teichmann, S. A. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science, 376(6594). https://doi.org/10.1126/science.abl5197
