# Image Registration with ImageJ Fijiyama

Helper python package for 3D image registration using the Image J package Fijiyama (https://imagej.net/plugins/fijiyama, R. Fernandez and C. Moisy, (2020) Fijiyama: a registration tool for 3D multimodal time-lapse imaging. *Bioinformatics*)

Feel free to distribute and use for whatever non-profit use. 

## Installation
`$ python -m pip install 'pyjiyama @ git+https://github.com/dsb-lab/pyjiyama'`

## Usage

For a ready to use example check `examples/example_registration.ipynb`

Let's start by importing the packages


```python
import numpy as np
import os
import subprocess
import sys
from pyjiyama import (square_stack4D, centroid_correction_3d_based_on_mid_plane, generate_fijiyama_file_system, 
                     generate_fijiyama_stacks, openfiji, remove_dir, create_transformations_folders, move_transformation,
                     create_dir, correct_path, read_img_with_resolution, get_file_names)
```

Now define your path to the data in the following way and the code you want to give to the registered series.


```python
home = os.path.expanduser("~")
path_data = "/path/to/data/"
embcode = "my_embryo_name"
```
#### DISCLAIMER: These first cells are only for the case in which you have your data as individual time files.

Load the file list in that folder


```python
files = get_file_names(path_data)
```

I like to make sure there is no file names in the list other that tif files in the following manner


```python
files = get_file_names(path_data)
for file in files: 
    if ".tif" not in file:
        files.remove(file) 
```

Usually these files are not properly ordered so that's the next step. 


```python
# sort filename 
current_order = []
for filename in files:
    idx = filename.find(".tif")
    filecode=int(filename[:idx])
    current_order.append(filecode)
    
idxs_sort = np.argsort(current_order) 
files = [files[idx] for idx in idxs_sort]
```

Now create 4D stack of the files you want to register


```python
_IMGS = []
for filename in files:
    _IMG, xyres, zres = read_img_with_resolution(path_data + filename)
    _IMGS.append(_IMG[0])

_IMGS = np.array(_IMGS)
```
#### DISCLAIMER: In case you have your stack as a 4D stack already, uncomment the first line of the next cell and start directly here.

Square stacks so that `x` and `y` dimensions are equal and remove non-squared stack to free some memory


```python
# _IMGS, xyres, zres = read_img_with_resolution(path_data + filename)
IMGS = square_stack4D(_IMGS)
del _IMGS
```

Make sure that there are no pixels over `255` (above saturation value, this can happen, for example, if you combine two channels for the registration) and convert to `uint8` for memory saving


```python
t, z, x, y = np.where(IMGS > 255)
IMGS[t, z, x, y] = 255
IMGS = IMGS.astype("uint8")
```

Now we have the data ready to start the registration process. 

### Centroid correction

Before passing the stack to fijiyama I like to ran a centroid correction in case there is a lot of traslation on the `xy` plane to reduce possible errors in fijiyama and improve fijiyama performance. If you have more than one channel pass them as a list of 4D arrays to the `extra_IMGS` keyword argument. If you do not do this, when you apply the fijiyama registration to the extra channels is going to look displaced since the fijiyama registration is going to be performed on the centroid corrected stacks.


```python
IMGS_corrected, extra_IMGS_corrected = centroid_correction_3d_based_on_mid_plane(IMGS, extra_IMGS=[])
del IMGS
IMGS_corrected = IMGS_corrected.astype("uint8")
extra_IMGS_corrected = np.array(extra_IMGS_corrected).astype("uint8")
```

### Fijiyama registration

Now we create the file system that Fijiyama requires. What this does is create a folder names `movies_registered` (or whatever name you want to give it). Inside this, it adds an output folder and an input folder where the split images are going to be saved


```python
(
    path_registered,
    path_output,
    path_movies_reg,
) = generate_fijiyama_file_system(
    path_data, "movies_registered", embcode
)
```

Now we save the stakcs for the individual times in the input folder. The format I use for the file names is `t{Time}.tif`


```python
generate_fijiyama_stacks(
    path_registered, IMGS_corrected, xyres, zres, file_format="t%d.tif"
)
```

Now we are ready to open Fiji and run Fijiyama as described in the instructions from Fijiyama plugin site (https://imagej.net/plugins/fijiyama). Put your path to ImageJ on `path_to_fiji`.


```python
path_to_fiji = "/opt/Fiji.app/ImageJ-linux64"
openfiji(path_to_fiji=path_to_fiji)
```

Now registration should be finished and is time to organize the remove temporal files created and organize the results. First, let't remove the temporal input stacks that we created. 


```python
remove_dir(path_registered)
```

Create a folder to store registration transformation files and move them from the Fijiyama created folder to our newly created folder 


```python
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
```

Now we remove the Fijiyama output folder since the only thing we care about now is the transformations files and we have already moved them from the Fijiyama output folder. If you still would like to keep this folder skip the following cell.


```python
remove_dir(path_output)

```

### Apply transformations

Is time to apply the obtained transformations to the original stack or to any other. For example, if we had an stack with more than one channel, we can now apply the transformations to every channel. 

We start by creating a temporal folder to save the 3D stacks to which we are going to apply the transformations. We also create the 5D stack with the shape `(C, T, Z, X, Y)` which we want to register.


```python
path_movies_reg_embcode = create_dir(
    path_movies_reg, embcode, return_path=True, rem=True
)
IMGS_chs = np.array([IMGS_corrected, *extra_IMGS_corrected])
```

Preallocate the array where the registered images will be allocated


```python
registered_IMGS_chs = np.zeros_like(np.array(IMGS_chs))
```

We can now loop for each channel and apply the corresponding transformations to each of the stacks using a beanshell script. However, before doing that make sure yout path to beanshell is correct. In my case is the following (If yours is different modify that line).


```python
pth_beanshell = "/opt/Fiji.app/beanshell/bsh-2.0b4.jar"
```

We can now apply the transformations now. The first line in the cell finds your pyjiyama installation to locate the beanshell script used for applying the transformations. 


```python
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

```

    Java HotSpot(TM) 64-Bit Server VM warning: ignoring option PermSize=128m; support was removed in 8.0
    Java HotSpot(TM) 64-Bit Server VM warning: Using incremental CMS is deprecated and will likely be removed in a future release


Reorder the stack from `(C, T, Z, X, Y)` to `(T, Z, C, X, Y)`


```python
sh = registered_IMGS_chs.shape
final_registered_IMGS_chs = np.zeros((sh[1], sh[2], sh[0], sh[3], sh[4])).astype(
    "uint8"
)
for ch in range(sh[0]):
    for t in range(sh[1]):
        for z in range(sh[2]):
            final_registered_IMGS_chs[t, z, ch] = registered_IMGS_chs[ch, t, z]
```

### Save results

Save the registered 5D hyperstack as a tif file. 


```python
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
```
