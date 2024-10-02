# Creating Custom ref model using CellTypist for single-cell annotation based on Tabula Muris dataset

import os
import scanpy as sc
import matplotlib
import scvi
import celltypist
from celltypist import models

# Enabling `force_update = True` will overwrite existing (old) models.
models.download_models(force_update = True)

# Path of models 
models.models_path

# Pre existing models 
models.models_description() 

# Read the reference dataset - Tabula Muris
# Define the directory where your files are saved
save_dir = "c:/Users/rshreya/Downloads/"

# Construct the full path to reference 
tm_droplet_path = os.path.join(save_dir, "TM_droplet.h5ad")

# Read the datasets with a backup URL in case the files are not found locally
tm_droplet = sc.read(
    tm_droplet_path,
    backup_url="https://figshare.com/ndownloader/files/23938934",
)

# Add condition or method 
tm_droplet.obs["tech"] = "10x"

# Save a copy of counts
tm_droplet.layers["counts"] = tm_droplet.X.copy()

# Preprocess and normalize counts of ref data
sc.pp.normalize_total(tm_droplet, target_sum=1e4)
sc.pp.log1p(tm_droplet) # important to log normalize 
tm_droplet.raw = tm_droplet  # keep full dimension safe
sc.pp.highly_variable_genes(
    tm_droplet,
    flavor="cell_ranger", #for 10x data, flavor is cell ranger  
    n_top_genes=5000,
    layer="counts",
    batch_key="tech",
    subset=True,
)

# Structure of reference dataset 
str(tm_droplet)

# Train entire ref model using celltypist.train
new_model = celltypist.train(tm_droplet[:, tm_droplet.var.highly_variable], labels = tm_droplet.obs['cell_ontology_class'], check_expression = False, use_SGD = True, mini_batch = True)

# Reassigning 
Celltypist_TM_Ref_Model = new_model

# Save model locally 
Celltypist_TM_Ref_Model.write('C:/Users/rshreya/Downloads/Celltypist_TM_Ref_Model.pkl')

# Save model under CellTypist models  
Celltypist_TM_Ref_Model.write(f'{models.models_path}/Celltypist_TM_Ref_Model.pkl')

# Load ref model 
model = models.Model.load(model = 'Celltypist_TM_Ref_Model.pkl')

model.cell_types
model.features

# Ref model specific to tissue origin - accurate for cell type prediction of specific tissue  
"""**Heart**"""
tm_droplet_heart = tm_droplet[
    (tm_droplet.obs.tissue == "Heart_and_Aorta")   #Change reference tissue 
    & (~tm_droplet.obs.cell_ontology_class.isna())
].copy()

tm_droplet_heart.obs["tech"] = "10x"
tm_droplet_heart.layers["counts"] = tm_droplet_heart.X.copy()

sc.pp.normalize_total(tm_droplet_heart, target_sum=1e4)
sc.pp.log1p(tm_droplet_heart)
tm_droplet_heart.raw = tm_droplet_heart  # keep full dimension safe

sc.pp.highly_variable_genes(
    tm_droplet_heart,
    flavor="cell_ranger", #flavor is cell ranger for 10x data 
    n_top_genes=5000,
    layer="counts",
    batch_key="tech",
    subset=True,
)

Celltypist_TM_Ref_Model_Heart = celltypist.train(tm_droplet_heart[:, tm_droplet_heart.var.highly_variable], labels = tm_droplet_heart.obs['cell_ontology_class'], check_expression = False, use_SGD = True, mini_batch = True)

# Save model locally 
Celltypist_TM_Ref_Model_Heart.write('C:/Users/rshreya/Downloads/Celltypist_TM_Ref_Model_Heart.pkl')

# Save model under celltypist models for direct use
Celltypist_TM_Ref_Model_Heart.write(f'{models.models_path}/Celltypist_TM_Ref_Model_Heart.pkl')

model_heart = models.Model.load(model = 'Celltypist_TM_Ref_Model_Heart.pkl')

model_heart.cell_types
model_heart.features

"""**Spleen**"""

tm_droplet_spleen = tm_droplet[
    (tm_droplet.obs.tissue == "Spleen")   #Change reference tissue 
    & (~tm_droplet.obs.cell_ontology_class.isna())
].copy()

tm_droplet_spleen.obs["tech"] = "10x"

tm_droplet_spleen.layers["counts"] = tm_droplet_spleen.X.copy()

sc.pp.normalize_total(tm_droplet_spleen, target_sum=1e4)
sc.pp.log1p(tm_droplet_spleen)
tm_droplet_spleen.raw = tm_droplet_spleen  # keep full dimension safe

sc.pp.highly_variable_genes(
    tm_droplet_spleen,
    flavor="cell_ranger", #flavor is cell ranger for 10x data 
    n_top_genes=5000,
    layer="counts",
    batch_key="tech",
    subset=True,
)

new_model_spleen = celltypist.train(tm_droplet_spleen[:, tm_droplet_spleen.var.highly_variable], labels = tm_droplet_spleen.obs['cell_ontology_class'], check_expression = False, use_SGD = True, mini_batch = True)

Celltypist_TM_Ref_Model_Spleen = new_model_spleen

#Write out the model.
Celltypist_TM_Ref_Model_Spleen.write('C:/Users/rshreya/Downloads/Celltypist_TM_Ref_Model_Spleen.pkl')
Celltypist_TM_Ref_Model_Spleen.write(f'{models.models_path}/Celltypist_TM_Ref_Model_Spleen.pkl')

model_spleen = models.Model.load(model = 'Celltypist_TM_Ref_Model_Spleen.pkl')

model_spleen.cell_types
model_spleen.features


"""**BM**"""

tm_droplet_bm = tm_droplet[
    (tm_droplet.obs.tissue == "Marrow")   #Change reference tissue 
    & (~tm_droplet.obs.cell_ontology_class.isna())
].copy()


tm_droplet_bm.obs["tech"] = "10x"

tm_droplet_bm.layers["counts"] = tm_droplet_bm.X.copy()

sc.pp.normalize_total(tm_droplet_bm, target_sum=1e4)
sc.pp.log1p(tm_droplet_bm)
tm_droplet_bm.raw = tm_droplet_bm  # keep full dimension safe

sc.pp.highly_variable_genes(
    tm_droplet_bm,
    flavor="cell_ranger", #flavor is cell ranger for 10x data 
    n_top_genes=5000,
    layer="counts",
    batch_key="tech",
    subset=True,
)

Celltypist_TM_Ref_Model_BM = celltypist.train(tm_droplet_bm[:, tm_droplet_bm.var.highly_variable], labels = tm_droplet_bm.obs['cell_ontology_class'], check_expression = False, use_SGD = True, mini_batch = True)

Celltypist_TM_Ref_Model_BM.write('C:/Users/rshreya/Downloads/Celltypist_TM_Ref_Model_BM.pkl')
Celltypist_TM_Ref_Model_BM.write(f'{models.models_path}/Celltypist_TM_Ref_Model_BM.pkl')

model_bm = models.Model.load(model = 'Celltypist_TM_Ref_Model_BM.pkl')

model_bm.cell_types
model_bm.features


