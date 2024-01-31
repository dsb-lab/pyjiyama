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
        if imagej_metadata is None: 
            xyres = 1
            zres = 1
        else:
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
            return correct_path(path + dir)
        else:
            return
    except FileExistsError:
        if rem:
            remove_dir(path + dir)
            create_dir(path, dir)
        else:
            pass

        if return_path:
            return correct_path(path + dir)
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
    pad_val: int
        value used for padding
        
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
    bot_ishdiff = -np.int32(np.floor(ishdiff/2))
    if bot_ishdiff==0:
        bot_ishdiff=None
    jshdiff = required_size-sh[1]
    lef_jshdiff = np.int32(np.ceil(jshdiff/2))
    rig_jshdiff = -np.int32(np.floor(jshdiff/2))
    if rig_jshdiff==0:
        rig_jshdiff=None
    IMG_padded = np.ones((required_size, required_size))*pad_val
    IMG_padded[top_ishdiff:bot_ishdiff, lef_jshdiff:rig_jshdiff] = IMG
    return IMG_padded


def pad_image_and_square_array3D(stack, required_size_xy=None, required_size_z=None, pad_val=0):
    """Squares and pads input image in 3D
    
    Parameters
    ----------
    stack : ndarray
        3D ndarray with the unsquared, unpadded image
    required_size_xy: None or int
        if None, pads and squares to max xy dim of IMG
        if int, pads both x and y dimensions to that number
    required_size_z: None or int
        if None, does not pad on z dim
        if int, pads stack whole xy planes until given value starting from the top
    pad_val: int
        value used for padding

    Returns
    -------
    new_stack:  ndarray
        3d ndarray with the image squared and padded
    """
    slices = stack.shape[0]
    if required_size_z is None:
        new_slices=slices
    else:
        new_slices=required_size_z
    pad_val = np.uint8(np.rint(pad_val))
    testimg = pad_image_and_square_array2D(stack[0], required_size=required_size_xy, pad_val=pad_val)
    new_stack = np.ones((new_slices, *testimg.shape), dtype="uint8")*pad_val
    offset = np.floor((new_slices-slices)/2).astype('int32')
    for z in range(slices):
        new_stack[z+offset] = pad_image_and_square_array2D(stack[z], required_size=required_size_xy, pad_val=pad_val)
    return new_stack



import matplotlib.pyplot as plt

from copy import copy, deepcopy

import numpy as np
from matplotlib.transforms import TransformedPatchPath
from matplotlib.widgets import LassoSelector, Slider


class Slider_t(Slider):
    def __init__(self, *args, **kwargs):
        Slider.__init__(self, *args, **kwargs)
        ax = kwargs["ax"]
        vmin = kwargs["valmin"]
        vmax = kwargs["valmax"]
        vstp = kwargs["valstep"]
        colr = kwargs["initcolor"]
        for v in range(vmin + vstp, vmax, vstp):
            vline = ax.axvline(
                v, 0, 1, color=colr, lw=1, clip_path=TransformedPatchPath(self.track)
            )

class Slider_z(Slider):
    def __init__(self, *args, **kwargs):
        Slider.__init__(self, *args, **kwargs)
        ax = kwargs["ax"]
        vmin = kwargs["valmin"]
        vmax = kwargs["valmax"]
        vstp = kwargs["valstep"]
        colr = kwargs["initcolor"]
        for v in range(vmin + vstp, vmax, vstp):
            vline = ax.axvline(
                v, 0, 1, color=colr, lw=1, clip_path=TransformedPatchPath(self.track)
            )

from matplotlib.widgets import RectangleSelector
class MyRectangleSelector(RectangleSelector):
    def __init__(self, *args, **kwargs):
        sh = kwargs.pop('shape')
        self.fig = kwargs.pop('fig')
        self.tmp_xlims = [0, sh[-1]]
        self.tmp_ylims = [0, sh[-2]]
        args = (args[0], self.selection_callback)
        RectangleSelector.__init__(self, *args, **kwargs)

    def selection_callback(self, eclick, erelease): 
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        self.tmp_xlims[0] = min(x1, x2)
        self.tmp_xlims[1] = max(x1, x2)
        self.tmp_ylims[0] = min(y1, y2)
        self.tmp_ylims[1] = max(y1, y2)
        
class PyjiyamaPlotter(object):
    def __init__(self, IMGS):
        
        fig, ax = plt.subplots(figsize=(10,10))
        self.fig = fig
        self.ax = ax
        
        self.IMGS = IMGS
        self.t = 1
        self.z = 1
        self.times = len(self.IMGS)
        self.slices = self.IMGS[0].shape[0]
        
        self.ctrl_press = self.fig.canvas.mpl_connect(
            "key_press_event", self.on_key_press
        )
        self.ctrl_release = self.fig.canvas.mpl_connect(
            "key_release_event", self.on_key_release
        )
        
        self.ctrl_is_held = False

    def on_key_press(self, event):
        if event.key == "control":
            self.ctrl_is_held = True
        elif event.key == 'enter':
            self.xlims = np.rint(self.RS.tmp_xlims).astype('int32')
            self.ylims = np.rint(self.RS.tmp_ylims).astype('int32')
            plt.close()

    def on_key_release(self, event):
        if event.key == "control":
            self.ctrl_is_held = False
    
    def update_slider_t(self, t):
        self.t = t
        self.replot_axis()
        self.update()

    def update_slider_z(self, z):
        self.z = z
        self.replot_axis()
        self.update()

    def time_scroll(self, event):
        if event.button == "up":
            self.t = self.t + 1
        elif event.button == "down":
            self.t = self.t - 1
        
        self.t = max(self.t, 1)
        self.t = min(self.t, self.times)
        self.set_val_t_slider(self.t)

    def slice_scroll(self, event):
        if event.button == "up":
            self.z = self.z - 1
        elif event.button == "down":
            self.z = self.z + 1

        self.z = max(self.z, 1)
        self.z = min(self.z, self.slices)
        self.set_val_z_slider(self.z)

    def onscroll(self, event):
        if self.ctrl_is_held:
            self.time_scroll(event)
        else:
            self.slice_scroll(event)

    def update(self):
        self.fig.canvas.draw_idle()


    def plot(
        self,
    ):
        
        # Make a horizontal slider to control the time.
        axslide = self.fig.add_axes([0.10, 0.01, 0.75, 0.03])

        sliderstr2 = "/%d)" % (self.times)
        valfmt = "%d" + sliderstr2

        time_slider = Slider_t(
            ax=axslide,
            label="time",
            initcolor="r",
            valmin=1,
            valmax=self.times,
            valinit=1,
            valstep=1,
            valfmt=valfmt,
            track_color=[0.8, 0.8, 0, 0.5],
            facecolor=[0.8, 0.8, 0, 1.0],
        )

        # Make a horizontal slider to control the zs.
        axslide = self.fig.add_axes([0.10, 0.04, 0.75, 0.03])
        sliderstr = "/%d" % (self.slices)
        zslide_val_fmt ="%d" + sliderstr
        z_slider = Slider_z(
            ax=axslide,
            label="z slice",
            initcolor="r",
            valmin=1,
            valmax=self.slices,
            valinit=1,
            valstep=1,
            valfmt=zslide_val_fmt,
            track_color=[0, 0.7, 0, 0.5],
            facecolor=[0, 0.7, 0, 1.0],
        )

        # Point to sliders
        time_slider.on_changed(self.update_slider_t)
        self.set_val_t_slider = time_slider.set_val

        z_slider.on_changed(self.update_slider_z)
        self.set_val_z_slider = z_slider.set_val
        
        scl = self.fig.canvas.mpl_connect("scroll_event", self.onscroll)

        self.ax.axis(False)
        _ = self.ax.set_xticks([])
        _ = self.ax.set_yticks([])
        
        t = 0
        z = 0
        img = self.IMGS[t][z, :, :]
        self.imshow = self.ax.imshow(img, vmin=0, vmax=255)

        self.zlimits=[0, self.slices-1]
        ZPicker(self.ax, self.fig.canvas, self.zpicker_callback)
        self.RS = MyRectangleSelector(
            self.ax,
            useblit=True,
            button=[3], 
            minspanx=5, minspany=5,
            spancoords='pixels',
            interactive=True,
            fig = self.fig,
            shape = img.shape
            )
        plt.subplots_adjust(bottom=0.075)
        plt.show()

    def replot_axis(self):
        img = self.IMGS[self.t-1][self.z-1, :, :]
        self.imshow.set_array(img)

    def zpicker_callback(self):
        if self.ctrl_is_held:
            self.zlimits[1] = self.z
        else:
            self.zlimits[0] = self.z
        
        print("current z limits", self.zlimits)
        
class ZPicker(object):
    def __init__(self, ax, canvas, callback):
        self.ax = ax
        self.cid = canvas.mpl_connect("button_press_event", self)
        self.canvas = canvas
        self.callback = callback
        
    def __call__(self, event):
        if event.dblclick == True:
            if event.button == 1:
                if event.inaxes == self.ax:
                    self.callback()

    def stopit(self):
        self.canvas.mpl_disconnect(self.cid)