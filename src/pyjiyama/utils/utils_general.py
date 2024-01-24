import os
import shutil

import numpy as np
from tifffile import TiffFile


def get_file_names(path_data):
    files = os.listdir(path_data)
    return files

def get_ordered_tif_files(path_data):
    total_files = get_file_names(correct_path(path_data))
    for i in range(len(total_files)).__reversed__(): 
        file = total_files[i]
        if ".tif" not in file:
            _ = total_files.pop(i) 
            
    # # sort filename 
    # current_order = []
    # for filename in total_files:
    #     if ".tif" not in filename: continue
    #     idx = filename.find(".tif")
    #     filecode=int(filename[idx-4:idx])
    #     current_order.append(filecode)
        
    # idxs_sort = np.argsort(currenWt_order) 
    # total_files = [total_files[idx] for idx in idxs_sort]

    # import os
    # offset = 0
    # for i, file in enumerate(total_files):
    #     os.rename(correct_path(path_data+embcode)+file, correct_path(path_data+embcode)+"{}.tif".format(i+offset))

    # sort filename 

    for i in range(len(total_files)).__reversed__(): 
        file = total_files[i]
        if ".tif" not in file:
            total_files.pop(i) 
            
    current_order = []
    for filename in total_files:
        if ".tif" not in filename: continue
        idx = filename.find(".tif")
        filecode=int(filename[:idx])
        current_order.append(filecode)
        
    idxs_sort = np.argsort(current_order) 
    total_files = [total_files[idx] for idx in idxs_sort]
    return total_files

def get_file_embcode(path_data, f, returnfiles=False):
    """
    Parameters
    ----------
    path_data : str
        The path to the directory containing emb
    f : str or int
        if str returns path_data/emb
        if int returns the emb element in path_data

    Returns
    -------
    file, name
        full file path and file name.
    """
    files = os.listdir(path_data)

    fid = -1
    if isinstance(f, str):
        for i, file in enumerate(files):
            if f in file:
                fid = i

        if fid == -1:
            raise Exception("given file name extract is not present in any file name")
    else:
        fid = f

    if fid > len(files):
        raise Exception("given file index is greater than number of files")

    file = files[fid]
    name = file.split(".")[0]
    if returnfiles:
        return file, name, files
    return file, name


def read_img_with_resolution(path_to_file, channel=None, stack=True):
    """
    Parameters
    ----------
    path_to_file : str
        The path to the tif file.
    channel : int or None
        if None assumes the tif file contains only one channel
        if int selects that channel from the tif

    Returns
    -------
    IMGS, xyres, zres
        4D numpy array with shape (t, z, x, y), x and y resolution and z resolution
    """
    with TiffFile(path_to_file) as tif:
        preIMGS = tif.asarray()
        shapeimg = preIMGS.shape
        if stack:
            if channel == None:
                if len(shapeimg) == 3:
                    IMGS = np.array([tif.asarray()])
                else:
                    IMGS = np.array(tif.asarray())
            else:
                if len(shapeimg) == 4:
                    IMGS = np.array([tif.asarray()[:, channel, :, :]])
                else:
                    IMGS = np.array(tif.asarray()[:, :, channel, :, :])
        else:
            if channel == None:
                if len(shapeimg) == 2:
                    IMGS = np.array([tif.asarray()])
                else:
                    IMGS = np.array(tif.asarray())
            else:
                if len(shapeimg) == 3:
                    IMGS = np.array([tif.asarray()[channel, :, :]])
                else:
                    IMGS = np.array(tif.asarray()[:, channel, :, :])

        if len(IMGS.shape) == 3:
            IMGS = np.array([IMGS])
        imagej_metadata = tif.imagej_metadata
        tags = tif.pages[0].tags
        # parse X, Y resolution
        try:
            npix, unit = tags["XResolution"].value
            xres = unit / npix
        except KeyError:
            xres = None

        try:
            npix, unit = tags["YResolution"].value
            yres = unit / npix
        except KeyError:
            yres = None

        try:
            zres = imagej_metadata["spacing"]
        except KeyError:
            zres = None

        if xres == yres:
            xyres = xres
        else:
            xyres = (xres, yres)
    return IMGS, xyres, zres

def get_max_xy_z_resolutions(path_data, files, file_step=30):
    maxsize_xy=0
    maxsize_z =0
    for file in files[::file_step]:
        _IMG, xyres, zres = read_img_with_resolution(correct_path(path_data) + file)
        max_dim = np.max(_IMG.shape[-2:])
        zdim = _IMG.shape[1]
        maxsize_xy = max(maxsize_xy, max_dim)
        maxsize_z = max(maxsize_z, zdim)

    return maxsize_xy, maxsize_z


def remove_dir(path, dir=""):
    try:
        shutil.rmtree(path + dir)
    except FileNotFoundError:
        return


def create_dir(path, dir="", rem=False, return_path=False):
    if dir != "":
        path = correct_path(path)
    try:
        os.mkdir(path + dir)
        if return_path:
            return path + dir
        else:
            return
    except FileExistsError:
        if rem:
            remove_dir(path + dir)
            create_dir(path, dir)
        else:
            pass

        if return_path:
            return path + dir
        else:
            return

    raise Exception("something is wrong with the dir creation")


def correct_path(path):
    if path[-1] != "/":
        path = path + "/"
    return path


def square_stack2D(img):
    x, y = img.shape
    if x == y:
        return img

    x, y = img.shape

    xydif = np.abs(x - y)
    crop = xydif / 2
    left = np.floor(crop).astype("int32")
    right = np.ceil(crop).astype("int32") * -1

    if x > y:
        new_img = img[left:right, :]
    else:
        new_img = img[:, left:right]

    return new_img


def square_stack3D(stack):
    slices = stack.shape[0]
    testimg = square_stack2D(stack[0])
    new_stack = np.zeros((slices, *testimg.shape), dtype="uint8")
    for z in range(slices):
        new_stack[z] = square_stack2D(stack[z])
    return new_stack


def square_stack4D(hyperstack):
    times = hyperstack.shape[0]
    teststack = square_stack3D(hyperstack[0])
    new_hyperstack = np.zeros((times, *teststack.shape), dtype="uint8")
    for t in range(times):
        new_hyperstack[t] = square_stack3D(hyperstack[t])
    return new_hyperstack


def pad_image_and_square_array2D(IMG, required_size=None, pad_val=0):
    """Squares and pads input image
    
    Parameters
    ----------
    IMG : ndarray
        2D ndarray with the unsquared, unpadded image
    required_size: None or int
        if None, pads and squares to max dim of IMG
        if int, pads both dimensions to that number

    Returns
    -------
    IMG_padded:  ndarray
        2d ndarray with the image squared and padded
    """
    
    # indexed to be used for image reconstruction
    sh = IMG.shape
    if required_size is None:
        required_size=max(*sh)
    ishdiff = required_size-sh[0]
    top_ishdiff = np.int32(np.ceil(ishdiff/2))
    bot_ishdiff = np.int32(np.floor(ishdiff/2))
    jshdiff = required_size-sh[1]
    lef_jshdiff = np.int32(np.ceil(jshdiff/2))
    rig_jshdiff = np.int32(np.floor(jshdiff/2))
    
    IMG_padded = np.ones((required_size, required_size))*pad_val
    IMG_padded[top_ishdiff:-bot_ishdiff, +lef_jshdiff:-rig_jshdiff] = IMG
    return IMG_padded


def pad_image_and_square_array3D(stack, required_size_xy=None, required_size_z=None, pad_val=0):
    slices = stack.shape[0]
    if required_size_z is None:
        new_slices=slices
    else:
        new_slices=required_size_z
    testimg = pad_image_and_square_array2D(stack[0], required_size=required_size_xy, pad_val=pad_val)
    new_stack = np.zeros((new_slices, *testimg.shape), dtype="uint8")
    offset = np.floor((new_slices-slices)/2).astype('int32')
    for z in range(slices):
        new_stack[z+offset] = pad_image_and_square_array2D(stack[z], required_size=required_size_xy,pad_val=pad_val)
    return new_stack


def pad_image_and_square_array4D(hyperstack, required_size_xy=None, required_size_z=None, pad_val=0):
    times = hyperstack.shape[0]
    teststack = pad_image_and_square_array3D(hyperstack[0], required_size_xy=required_size_xy, required_size_z=required_size_z, pad_val=pad_val)
    new_hyperstack = np.ones((times, *teststack.shape), dtype="uint8")*np.uint8(np.rint(pad_val))
    for t in range(times):
        new_hyperstack[t] = pad_image_and_square_array3D(hyperstack[t], required_size_xy=required_size_xy, required_size_z=required_size_z,pad_val=pad_val)
    return new_hyperstack
