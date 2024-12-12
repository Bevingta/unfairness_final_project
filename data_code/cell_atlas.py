# 2. Option A: Download from Cell Atlas dataset
import scanpy as sc

# Download a specific dataset (e.g., Paul et al. data)
adata = sc.datasets.paul15()

# Save the dataset locally if needed
adata.write('paul15_data.h5ad')

import urllib.request
import os

# Create a directory for the data
os.makedirs('data', exist_ok=True)

# Download 10X PBMC dataset
url = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5"
urllib.request.urlretrieve(url, "data/pbmc_1k_v3.h5")

# Load the downloaded 10X data
adata = sc.read_10x_h5("data/pbmc_1k_v3.h5")

# Save in h5ad format
adata.write('data/pbmc_1k_v3.h5ad')
