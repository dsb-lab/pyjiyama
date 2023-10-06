import os
import shutil

import numpy as np
from tifffile import TiffFile


def get_file_names(path_data):
    files = os.listdir(path_data)
    return files


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
        return

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
