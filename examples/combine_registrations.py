import numpy as np
import os
from pyjiyama import (square_stack4D, centroid_correction_3d_based_on_mid_plane, generate_fijiyama_file_system, 
                     generate_fijiyama_stacks, openfiji, remove_dir, create_transformations_folders, move_transformation,
                     create_dir, correct_path, read_img_with_resolution, get_file_names)

home = os.path.expanduser("~")
path_data = home + "/Desktop/PhD/projects/Data/blastocysts/Lana/20230607_CAG_H2B_GFP_16_cells/stack_2_channel_0_obj_bottom/crop/"
embcode = "20230607_CAG_H2B_GFP_16_cells_stack2"

files = get_file_names(path_data)
for file in reversed(files):
    if ".tif" not in file:
        files.remove(file) 
        continue
    
    if "part" not in file:
        files.remove(file)
        
# sort filename 
current_order = []
for filename in files:
    idx = filename.find(".tif")
    filecode=int(filename[idx-1:idx])
    current_order.append(filecode)
        
idxs_sort = np.argsort(current_order) 
files = [files[idx] for idx in idxs_sort]

_IMGS = []
for filename in files:
    _IMG, xyres, zres = read_img_with_resolution(path_data + filename)
    _IMGS.extend(_IMG)


merged_IMGS = np.array(_IMGS)
del _IMGS

IMGS_corrected = merged_IMGS[49:51]
extra_IMGS_corrected = merged_IMGS[51:]

(
    path_registered,
    path_output,
    path_movies_reg,
) = generate_fijiyama_file_system(
    path_data, "movies_registered", embcode
)

generate_fijiyama_stacks(
    path_registered, IMGS_corrected, xyres, zres, file_format="t%d.tif"
)

path_to_fiji = "/opt/Fiji.app/ImageJ-linux64"
openfiji(path_to_fiji=path_to_fiji)

remove_dir(path_registered)

# Create transformation folders
(
    path_trans_emb_global,
    path_trans_emb_steps,
) = create_transformations_folders(
    path_movies_reg, embcode, trans_folder="transformations"
)

# Move transformations
move_transformation(
    path_output, path_trans_emb_global, path_trans_emb_steps
)


remove_dir(path_output)

path_movies_reg_embcode = create_dir(
    path_movies_reg, embcode, return_path=True, rem=True
)

import shutil
files_trans = get_file_names(path_trans_emb_global)
file_trans = list(filter(lambda x: "2" in x, files_trans))[0]

for i in range(extra_IMGS_corrected.shape[1]):
    file_trans_new = file_trans.replace("2.txt", str(2+i+1)+".txt")
    shutil.copy(correct_path(path_trans_emb_global)+file_trans,correct_path(path_trans_emb_global)+file_trans_new)


IMGS_chs = np.array([np.concatenate((IMGS_corrected, extra_IMGS_corrected), axis=0)])
registered_IMGS_chs = np.zeros_like(np.array(IMGS_chs))

pth_beanshell = "/opt/Fiji.app/beanshell/bsh-2.0b4.jar"

from pyjiyama import __file__ as pyjiyama_path

for ch, IMGS_ch in enumerate(IMGS_chs):
    path_movies_reg_embcode_ch = create_dir(
        path_movies_reg_embcode, "%d" % ch, return_path=True, rem=True
    )
    path_movies_reg_embcode_ch_reg = create_dir(
        path_movies_reg_embcode, "registered_%d" % ch, return_path=True, rem=True
    )

    generate_fijiyama_stacks(
        path_movies_reg_embcode_ch,
        IMGS_ch,
        xyres,
        zres,
        file_format="t%d.tif",
        rem=True,
    )

    # Define where you have the beanshell class to be called from beanshell
    text_to_write = "\n".join(
        [
            pth_beanshell,
            path_trans_emb_global,
            path_movies_reg_embcode_ch,
            correct_path(path_movies_reg_embcode_ch_reg),
        ]
    )
    # Save path information in a text file to be open in beanshell.
    temporal_file = correct_path(home) + "tmp.txt"
    with open(temporal_file, "w") as the_file:
        the_file.write(text_to_write)

    # Run Beanshell script
    idx = pyjiyama_path.rfind('/')
    pth_beanshell_script = (
        correct_path(
            pyjiyama_path[:idx]
        )
        + "utils/apply_transformation.bsh"
    )
    import subprocess

    subprocess.run(
        [path_to_fiji, "--headless", pth_beanshell_script]
    )  # , stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    remove_dir(path_movies_reg_embcode_ch)
    # Remove path file
    os.remove(temporal_file)

    tfiles = os.listdir(path_movies_reg_embcode_ch_reg)
    for _t, tfile in enumerate(tfiles):
        registered_IMGS_chs[ch]
        t = int(tfile.split(".")[0].split("t")[-1]) - 1
        IMG_t, xyres, zres = read_img_with_resolution(
            correct_path(path_movies_reg_embcode_ch_reg) + tfile,
            channel=None,
        )
        registered_IMGS_chs[ch][t] = IMG_t

registered_IMGS_chs = registered_IMGS_chs.astype("uint8")

merged_IMGS[49:] = registered_IMGS_chs[0]
registered_IMGS_chs = np.array([merged_IMGS])

sh = registered_IMGS_chs.shape
final_registered_IMGS_chs = np.zeros((sh[1], sh[2], sh[0], sh[3], sh[4])).astype(
    "uint8"
)

for ch in range(sh[0]):
    for t in range(sh[1]):
        for z in range(sh[2]):
            final_registered_IMGS_chs[t, z, ch] = registered_IMGS_chs[ch, t, z]

fullpath = path_movies_reg_embcode + ".tif"

mdata = {"axes": "TZCYX", "spacing": zres, "unit": "um"}
import tifffile

tifffile.imwrite(
    fullpath,
    final_registered_IMGS_chs,
    imagej=True,
    resolution=(1 / xyres, 1 / xyres),
    metadata=mdata,
)
remove_dir(path_movies_reg_embcode)