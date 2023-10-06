import numpy as np
import os
import sys
sys.path.append('/home/pablo/Desktop/PhD/projects/pyjiyama/src')
from pyjiyama import (square_stack4D, centroid_correction_3d_based_on_mid_plane, generate_fijiyama_file_system, 
                     generate_fijiyama_stacks, openfiji, remove_dir, create_transformations_folders, move_transformation,
                     create_dir, correct_path, read_img_with_resolution, get_file_names)
import os

# home = os.path.expanduser("~")
# path_data = home + "/Desktop/PhD/projects/Data/blastocysts/Lana/20230607_CAG_H2B_GFP_16_cells/stack_2_channel_0_obj_bottom/crop/volumes/"
# embcode = "20230607_CAG_H2B_GFP_16_cells_stack2_part2"

# files = get_file_names(path_data)

# for file in files: 
#     if ".tif" not in file:
#         files.remove(file) 
        
# embcode = "20230607_CAG_H2B_GFP_16_cells_stack2_part2"

# # sort filename 
# current_order = []
# for filename in files:
#     if ".tif" not in filename: continue
#     idx = filename.find(".tif")
#     filecode=int(filename[:idx])
#     current_order.append(filecode)
    
# idxs_sort = np.argsort(current_order) 
# files = [files[idx] for idx in idxs_sort]

# _IMGS = []
# for filename in files[50:]:
#     _IMG, xyres, zres = read_img_with_resolution(path_data + filename)
#     _IMGS.append(_IMG[0])

# _IMGS = np.array(_IMGS)

# IMGS = square_stack4D(_IMGS)
# del _IMGS

# # Combine channels into single stack
# t, z, x, y = np.where(IMGS > 255)
# IMGS[t, z, x, y] = 255
# IMGS = IMGS.astype("uint8")

# ### PREPROCESSING ###
# # Run centroid correction prior to Fijiyama registration to improve performance
# IMGS_corrected, extra_IMGS_corrected = centroid_correction_3d_based_on_mid_plane(IMGS, extra_IMGS=[])
# del IMGS
# IMGS_corrected = IMGS_corrected.astype("uint8")
# extra_IMGS_corrected = np.array(extra_IMGS_corrected).astype("uint8")

# ### FJIYAMA REGISTRATION ###
# # Create Fijiyama file system (input and output folders)
# (
#     path_registered,
#     path_output,
#     path_movies_reg,
# ) = generate_fijiyama_file_system(
#     path_data, "movies_registered", embcode
# )

# # Save registration stacks into input folder
# generate_fijiyama_stacks(
#     path_registered, IMGS_corrected, xyres, zres, file_format="t%d.tif"
# )

# # Open Imagej to run fijiyama registration
# path_to_fiji = "/opt/Fiji.app/ImageJ-linux64"
# openfiji(path_to_fiji=path_to_fiji)

# # Remove stacks used for registration
# remove_dir(path_registered)

# # Create transformation folders
# (
#     path_trans_emb_global,
#     path_trans_emb_steps,
# ) = create_transformations_folders(
#     path_movies_reg, embcode, trans_folder="transformations"
# )

# # Move transformations
# move_transformation(
#     path_output, path_trans_emb_global, path_trans_emb_steps
# )

# # Remove Fijiyama output folder
# remove_dir(path_output)

# ### APPLY TRANSFORMATIONS ###

# path_movies_reg_embcode = create_dir(
#     path_movies_reg, embcode, return_path=True, rem=True
# )
# IMGS_chs = np.array([IMGS_corrected])

# # Expand channels
# registered_IMGS_chs = np.zeros_like(np.array(IMGS_chs))

# pth_beanshell = "/opt/Fiji.app/beanshell/bsh-2.0b4.jar"
# for ch, IMGS_ch in enumerate(IMGS_chs):
#     path_movies_reg_embcode_ch = create_dir(
#         path_movies_reg_embcode, "%d" % ch, return_path=True, rem=True
#     )
#     path_movies_reg_embcode_ch_reg = create_dir(
#         path_movies_reg_embcode, "registered_%d" % ch, return_path=True, rem=True
#     )

#     generate_fijiyama_stacks(
#         path_movies_reg_embcode_ch,
#         IMGS_ch,
#         xyres,
#         zres,
#         file_format="t%d.tif",
#         rem=True,
#     )

#     # Define where you have the beanshell class to be called from beanshell
#     text_to_write = "\n".join(
#         [
#             pth_beanshell,
#             path_trans_emb_global,
#             path_movies_reg_embcode_ch,
#             correct_path(path_movies_reg_embcode_ch_reg),
#         ]
#     )
#     # Save path information in a text file to be open in beanshell.
#     temporal_file = correct_path(home) + "tmp.txt"
#     with open(temporal_file, "w") as the_file:
#         the_file.write(text_to_write)

#     # Run Beanshell script
#     pth_beanshell_script = (
#         correct_path(
#             "/home/pablo/Desktop/PhD/projects/pyjiyama/src/pyjiyama/"
#         )
#         + "utils/apply_transformation.bsh"
#     )
#     import subprocess

#     subprocess.run(
#         [path_to_fiji, "--headless", pth_beanshell_script]
#     )  # , stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

#     remove_dir(path_movies_reg_embcode_ch)
#     # Remove path file
#     os.remove(temporal_file)

#     tfiles = os.listdir(path_movies_reg_embcode_ch_reg)
#     for _t, tfile in enumerate(tfiles):
#         registered_IMGS_chs[ch]
#         t = int(tfile.split(".")[0].split("t")[-1]) - 1
#         IMG_t, xyres, zres = read_img_with_resolution(
#             correct_path(path_movies_reg_embcode_ch_reg) + tfile,
#             channel=None,
#         )
#         registered_IMGS_chs[ch][t] = IMG_t

# registered_IMGS_chs = registered_IMGS_chs.astype("uint8")


# # Reorder stack from CTZXY to the imagej structure TZCXY
# sh = registered_IMGS_chs.shape
# final_registered_IMGS_chs = np.zeros((sh[1], sh[2], sh[0], sh[3], sh[4])).astype(
#     "uint8"
# )
# for ch in range(sh[0]):
#     for t in range(sh[1]):
#         for z in range(sh[2]):
#             final_registered_IMGS_chs[t, z, ch] = registered_IMGS_chs[ch, t, z]

# # Save final stack and remove unncessary intermediate files
# fullpath = path_movies_reg_embcode + ".tif"
# mdata = {"axes": "TZCYX", "spacing": zres, "unit": "um"}
# import tifffile

# tifffile.imwrite(
#     fullpath,
#     final_registered_IMGS_chs,
#     imagej=True,
#     resolution=(1 / xyres, 1 / xyres),
#     metadata=mdata,
# )
# remove_dir(path_movies_reg_embcode)
