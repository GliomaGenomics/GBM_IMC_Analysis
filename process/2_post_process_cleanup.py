import os as os
from cellpose import io
import shutil as sh
import matplotlib.pyplot as plt

# Postprocessing directories
workdir = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/Data/cellpose/drive-download-20230529T071000Z-001/'
full_ROI_combined_dir = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/Data/cellpose/data/combined_2_channel_full_ROIs/'
png_mask_dir = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/Data/cellpose/outputs/masks/png/'
tif_mask_dir = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/Data/cellpose/outputs/masks/tif/'
outlines_dir = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/Data/cellpose/outputs/outlines/'
cp_output_dir = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/Data/cellpose/outputs/QC/'


seg_output_files = os.listdir(workdir)

for filename in seg_output_files:
    if not filename.startswith('.'):  # Exclude hidden files
        
        if filename.endswith('_combined.tiff'):
           sh.move(os.path.join(workdir, filename), os.path.join(full_ROI_combined_dir, filename))
        
        elif filename.endswith('_masks.png'):
           sh.move(os.path.join(workdir, filename), os.path.join(png_mask_dir, filename))
        
        elif filename.endswith('_masks.tif'):
            sh.move(os.path.join(workdir, filename), os.path.join(tif_mask_dir, filename))
        
        elif filename.endswith('cp_outlines.txt'):
            sh.move(os.path.join(workdir, filename), os.path.join(outlines_dir, filename))
        
        elif filename.endswith('cp_output.png'):
            sh.move(os.path.join(workdir, filename), os.path.join(cp_output_dir, filename))
        
        else:
            print('No files to move')



plot_out_dir = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/Data/cellpose/outputs/test/'

if not os.path.isdir(plot_out_dir):
    os.mkdir(plot_out_dir)


qc_plots = os.listdir(cp_output_dir)

plt.figure(figsize=(20,10), dpi=300)
plt.imshow(io.imread(os.path.join(cp_output_dir, qc_plots[0])))
plt.axis('off')
plt.savefig(os.path.join(plot_out_dir, qc_plots[0]))
