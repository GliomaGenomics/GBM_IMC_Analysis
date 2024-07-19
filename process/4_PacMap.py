import os
import pacmap
import numpy as np

data_dir = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/cellpose/data/downstream/dim_reductions/'

input_file = os.path.join(data_dir, 'harmony_embedding_matrix.csv')
output_file = os.path.join(data_dir, 'PacMAP_harmony_embeddings.csv')

# load in the processed single cell expression matrix
# you can change it with any dataset that is in the ndarray format, with the shape (N, D)
# where N is the number of samples and D is the dimension of each sample

data = np.loadtxt(fname=input_file, delimiter=',', skiprows=1)

# calculate the number of neighbors for PacMAP
if np.shape(data)[0] < 10000:
    neighbor_select = 10
else:
    neighbor_select = round(10 + 15 * (np.log10(np.shape(data)[0])-4))

# initializing the pacmap instance
# Setting n_neighbors to "None" leads to a default choice shown below in "parameter" section
embedding = pacmap.PaCMAP(n_components=2, n_neighbors=neighbor_select) 

# fit the data (The index of transformed data corresponds to the index of the original data)
data_transformed = embedding.fit_transform(data, init="pca")

# save the pacmap reduced dims
np.savetxt(output_file, data_transformed, delimiter=',')