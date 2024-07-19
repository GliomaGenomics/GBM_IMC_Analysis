from astir.data import from_csv_yaml
import os
import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Setting the directories
expression_mat_path = 'data/downstream/astir_exprs_matrix.csv'
yaml_marker_path = 'data/downstream/astir_marker.yml'
out_dir = 'outputs/cell_phenotyping/Astir'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Check to see if the expression matrix loads correctly
pd.read_csv(expression_mat_path, index_col=0).head()

# Creating an astir object using the expression matrix and the markers
ast = from_csv_yaml(expression_mat_path, marker_yaml=yaml_marker_path)
print(ast)

# Model Parameters ------------------------------------------------

# Create batch size proportional to the number of cells
N = ast.get_type_dataset().get_exprs_df().shape[0]
batch_size = int(N/100)
# Number of training epochs
max_epochs = 1000
# Set learning rate
learning_rate = 1e-3
# Set initial epochs
initial_epochs = 10

# --------------------------------------------------------------------

# Fitting the cell types
ast.fit_type(max_epochs = max_epochs,
             batch_size = batch_size,
             learning_rate = learning_rate,
             n_init_epochs = initial_epochs)

# Fitting the cell states
ast.fit_state(max_epochs = max_epochs,
             batch_size = batch_size,
             learning_rate = learning_rate,
             n_init_epochs = initial_epochs)


# plotting the learning rates for the cell_type and states 
plt.figure(figsize=(10,7))
plt.plot(np.arange(len(ast.get_type_losses())), ast.get_type_losses())
plt.ylabel("Loss")
plt.xlabel("Epoch")
plt.title("Cell Type Astir Learnging rate")
plt.savefig(os.path.join(out_dir,'celltype_learning_rate.png'), dpi=300)


plt.figure(figsize=(10,7))
plt.plot(np.arange(len(ast.get_state_losses())), ast.get_state_losses())
plt.ylabel("Loss")
plt.xlabel("Epoch")
plt.title("Cell State Astir Learnging rate")
plt.savefig(os.path.join(out_dir,'cellstate_learning_rate.png'), dpi=300)

# Astir automatically creates a design matrix based on on additional 
# covariates in the data such as batches. This is done using additional columns already 
# present in the input data. In the example data of this vignette we have included 
# a 'batch' column for the purposes of illustration.

# We can get the cell type assignments in one of two ways:

# get_celltypes(): this returns the most likely cell type or classifies a 
# cell as unknown if no cell type has a probability above 0.7. 
# (This threshold can be altered by the user with the threshold argument)

# get_celltype_probabilities(): this returns the probabilty of each cell being assigned to any 
# given cell type.

ast.get_celltype_probabilities()
assignments = ast.get_celltype_probabilities()
assignments

# Assignment counts
ast.get_celltypes().value_counts()

ast.get_celltypes(threshold=0.5).value_counts()

# Assignment probabilities
ast.get_celltype_probabilities()

ast.diagnostics_celltype()

ast.type_to_csv(os.path.join(out_dir, 'cell-types.csv'))
ast.state_to_csv(os.path.join(out_dir, 'cell-states.csv'))
