# This script takes the raw multi-stack .tiff files when have been converted from the propriety 
# .mcd files and combines the channels to generate a single RGB file per ROI and also generates
# crops of a defined size in order to be used for labelling cells for cellpose  

import os as os
import tifffile as tiff
import random
import pandas as pd
import numpy as np
import argparse
import math


# Preprocessing arguments
workdir = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/cellpose/data/'
panelpath = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/cellpose/data/panel/panel.csv'
imgdir = '/Users/shoaibaliajaib/Library/CloudStorage/OneDrive-UniversityofLeeds/Hyperion_Analysis/cellpose/data/fullstacks/'
imgsuffix = '_fullstack.tiff'
channel_nuc = 'nuclear'
# channel_mem = 'membrane'
channel_mem = 'cytoplasm'

# renaming all images in the fullstacks dir
files = os.listdir(imgdir)

for filename in files:
    if not filename.startswith('.'):  # Exclude hidden files
        old_path = os.path.join(imgdir, filename)
        new_filename = filename.replace('.tiff', imgsuffix)  # Modify the filename as desired
        new_path = os.path.join(imgdir, new_filename)
        os.rename(old_path, new_path)


def combine_channels(images_path, panel_path, channel_type, out_path, imgsuffix, imagej=True):
    images = []
    for fn in os.listdir(images_path):
        if fn.endswith(imgsuffix):
            images.append(fn)
    print('process images: {}'.format(images))
    print('process images number: {}'.format(len(images)))
    panel = pd.read_csv(panel_path)

    for i in range(len(images)):
        img = tiff.imread(os.path.join(images_path, images[i]))
        if img.shape[0] > 3:
            result = np.zeros((len(channel_type), img.shape[1], img.shape[2]))
            # ["nuclear","membrane_cytoplasm"]
            for j in range(len(channel_type)):
                channels = list(panel[panel[channel_type[j]] == 1].index)
                result[j, :, :] = np.sum(img[channels, :, :], axis=0)
            # print(result.shape)
            tiff.imwrite(os.path.join(out_path, images[i].replace(
                imgsuffix.split('.tiff')[0], '_combined')), result, imagej=imagej)


def combine2Img(images_path, combined_path, panel_info, imgsuffix, channel_nuc='nuclear', channel_mem='membrane_cytoplasm'):

    if not os.path.exists(combined_path):
        os.makedirs(combined_path)

    combine_channels(images_path=images_path,
                     panel_path=panel_info,
                     channel_type=[channel_nuc, channel_mem],
                     out_path=combined_path,
                     imgsuffix=imgsuffix,
                     imagej=False)


def cropImg(crop_combined_path, combined_path, crop_size = 100, n_image_crops = 3):

    if not os.path.exists(crop_combined_path):
        os.makedirs(crop_combined_path)
    
    files = [file for file in os.listdir(combined_path) if not file.startswith('.')]
    
    for file in files:
        img = tiff.imread(os.path.join(combined_path, file))
        
        num_channels = img.shape[0]  # Assumes the channels are the first dimension
        
        x_pixels = img.shape[1]
        y_pixels = img.shape[2]

        if (y_pixels > crop_size) and (x_pixels > crop_size):

            max_x = img.shape[2] - crop_size
            max_y = img.shape[1] - crop_size

            for crop in range(n_image_crops):
                
                # Generate random coordinates within the allowable range
                x = random.randint(0, max_x)
                y = random.randint(0, max_y)
                
                # Extract the crop from each channel
                cropped_image = img[:, y:y+crop_size, x:x+crop_size]
                
                # Save the crop as a new TIFF image
                tiff.imwrite(os.path.join(crop_combined_path,file.replace('combined.tiff', f"combined_{crop_size}x{crop_size}_{crop+1}.tiff")),
                        data = cropped_image)


combine2Img(images_path=imgdir, 
            combined_path= os.path.join(workdir, 'combined_nuclear_cytoplasm'), 
            panel_info=panelpath, imgsuffix=imgsuffix, 
            channel_nuc=channel_nuc, channel_mem=channel_mem)

cropImg(os.path.join(workdir, 'combined_nucleus_membrane__crops'),
        os.path.join(workdir, 'combined_nucleus_membrane_images'), 
        crop_size = 200, n_image_crops= 1)



# def mergeImg(prop_path, combined_path, out_path):

#     files = os.listdir(prop_path)
#     files = [i for i in files if '_pred_Probabilities' in i]

#     if not os.path.exists(out_path):
#         os.makedirs(out_path)

#     sampledic = {}

#     for each in files:
#         sample = each.split('_combined_cropx')[0]
#         if sample not in sampledic.keys():
#             sampledic[sample] = [each]
#         else:
#             sampledic[sample].append(each)

#     for each in sampledic.keys():
#         rois = sampledic[each]
#         imgs = {}
#         for roi in rois:
#             x = roi.split('_pred_Probabilities')[0].split('_crop')[1][1:]
#             y = roi.split('_pred_Probabilities')[0].split('_crop')[2][1:]
#             imgs[str(x)+'_'+str(y)] = tiff.imread(os.path.join(prop_path, roi))
#         original = tiff.imread(os.path.join(
#             combined_path, each+'_combined.tiff'))

#         z, y, x = original.shape
#         ori_prop = np.zeros((512*math.ceil(x/512), 512*math.ceil(x/512), 3))

#         for i in range(0, x, 512):
#             for j in range(0, y, 512):
#                 pos = str(i)+'_'+str(j)
#                 ori_prop[i:i+512, j:j+512] = imgs[pos]

#         tiff.imwrite('{}/{}_pred_Probabilities.tiff'.format(out_path,
#                      each), ori_prop.astype('uint16'))








# if __name__ == '__main__':
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--process', type=str, help='pre or post')
#     parser.add_argument('--workdir', type=str,
#                         default="/workspace/data/predict_BRCA2/", help='the path of the data folder')
#     parser.add_argument('--propn', type=str, default="BRCA1_threshold-99.7_withAugnoise-0.5",
#                         help='the name of predicted probability in the data folder')
#     parser.add_argument('--panel', type=str, default='/workspace/data/panel/BRCA2.csv',
#                         help='the path of the panel file => same ordered with img channel: at least need contain two columns("nuclear","membrane_cytoplasm")')
#     parser.add_argument('--channel_n', type=str, default='nuclear',
#                         help='the column name for nuclear channel in panel file')
#     parser.add_argument('--channel_m', type=str, default='membrane_cytoplasm',
#                         help='the column name for membrane channel in panel file')
#     parser.add_argument('--imgdir', type=str, default="/workspace/data/predict_BRCA2/fullstacks/",
#                         help='the path of the img data folder')
#     parser.add_argument('--imgsuf', type=str,
#                         default='_fullstack.tiff', help='images suffix')
#     args = parser.parse_args()

#     workdir = args.workdir
#     process = args.process
#     panelpath = args.panel
#     imgsuffix = args.imgsuf
#     imgdir = args.imgdir
#     prop_path_name = args.propn
#     channel_mem = args.channel_m
#     channel_nuc = args.channel_n

#     if process == 'pre':
#         combine2Img(imgdir, os.path.join(workdir, 'combined2_image'), panelpath,
#                     imgsuffix, channel_nuc=channel_nuc, channel_mem=channel_mem)
#         cropImg(os.path.join(workdir, 'crop512_combined2_image'),
#                 os.path.join(workdir, 'combined2_image'))
#     elif process == 'post':
#         prop_path = os.path.join(
#             workdir, 'predictionProbability', prop_path_name)
#         combined_path = os.path.join(workdir, 'combined2_image')
#         out_path = os.path.join(prop_path, 'ori_prop')
#         mergeImg(prop_path, combined_path, out_path)
