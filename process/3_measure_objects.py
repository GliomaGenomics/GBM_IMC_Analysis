import steinbock.measurement.intensities as intensity
import steinbock.measurement.regionprops as regions
import steinbock.io as sio
import os
import numpy as np
import pandas as pd
import steinbock.utils as sutils

# Set the directories
outpath = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/Data/cellpose/outputs/'

intensity_out = os.path.join(outpath, 'intensities')
regions_out = os.path.join(outpath, 'regionprops')

if not os.path.exists(intensity_out):
    os.makedirs(intensity_out)

if not os.path.exists(regions_out):
    os.makedirs(regions_out)


# Marker Panel
panelpath = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/Data/cellpose/data/panel/panel.csv'
channel_col = 'name'
panel_file = pd.read_csv(panelpath)

if channel_col in panel_file.columns:
    channel_names = panel_file[channel_col].tolist()
    print(f"\nFound {len(channel_names)} channel names")
else:
    print(f"channel_col is not present in the panel_file")


# Segmentation Masks
mask_path = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/Data/cellpose/outputs/masks/tif/'
mask_suffix = '_combined_cp_masks.tif'
masks = np.sort(os.listdir(mask_path))
masks = [i for i in masks if mask_suffix in i]

# Full-stack Images
image_path = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/Data/cellpose/data/fullstacks/'
image_suffix = '_fullstack.tiff'
images = np.sort(os.listdir(image_path))
images = [i for i in images if image_suffix in i]

image_id = [i.split(image_suffix)[0] for i in images] 
mask_id = [i.split(mask_suffix)[0] for i in masks]

if image_id == mask_id:
    print(f"\nprocessing {len(image_id)} images ...")

    for img_file, mask_file, file_id in zip(images, masks, image_id):

        print(f"\nprocessing {file_id} ...")
        outfile = f"{file_id}.csv"

        img = sio.read_image(os.path.join(image_path, img_file))
        mask = sio.read_mask(os.path.join(mask_path, mask_file))

        print('\tCalculating intensities ...')
        marker_exp = intensity.measure_intensites(img=img, mask=mask, channel_names=channel_names,
                                                intensity_aggregation=intensity.IntensityAggregation.MEAN)
        
        print('\tCalculating region props ...')
        region_measures = ['area','centroid', 'eccentricity', 'axis_major_length','axis_minor_length']
        region_props = regions.measure_regionprops(img=img, mask=mask, skimage_regionprops=region_measures)
    
        marker_exp.to_csv(os.path.join(intensity_out, outfile), )
        region_props.to_csv(os.path.join(regions_out, outfile))

else:
    print("Images and masks have a different order and/or lengths")


