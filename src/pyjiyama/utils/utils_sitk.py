
import numpy as np
import SimpleITK as sitk
from .utils_general import get_file_names

def to_affine_array_representation(transformation):
    pt0 = [0, 0, 0]
    pt0_tr = transformation.TransformPoint(pt0)
    pt1 = [1, 0, 0]
    pt1_tr = transformation.TransformPoint(pt1)
    pt2 = [0, 1, 0]
    pt2_tr = transformation.TransformPoint(pt2)
    pt3 = [0, 0, 1]
    pt3_tr = transformation.TransformPoint(pt3)

    return np.array([
        [pt1_tr[0] - pt0_tr[0], pt2_tr[0] - pt0_tr[0], pt3_tr[0] - pt0_tr[0], pt0_tr[0]],
        [pt1_tr[1] - pt0_tr[1], pt2_tr[1] - pt0_tr[1], pt3_tr[1] - pt0_tr[1], pt0_tr[1]],
        [pt1_tr[2] - pt0_tr[2], pt2_tr[2] - pt0_tr[2], pt3_tr[2] - pt0_tr[2], pt0_tr[2]],
        [0, 0, 0, 1]
    ])
    

def combine_transformations_with_preview(pre_transform_full_path, transforms_paths):
    tr = sitk.ReadTransform(pre_transform_full_path)
    transforms = [tr]
    filenames = get_file_names(transforms_paths)
    filenames = [fn for fn in filenames if ".txt" in fn]
    for time in range(1,len(filenames)+1):
        full_file_name = transforms_paths + 'transform_global_t{}.txt'.format(time)
        tr = sitk.ReadTransform(full_file_name)
        transforms.append(tr)
                
        global_transformation_composite = sitk.CompositeTransform([transforms[0], transforms[time]])    

        affine_array_representation = to_affine_array_representation(global_transformation_composite)

        rotation = [i for sublist in affine_array_representation[:3,:3] for i in sublist]
        translation = [i for i in affine_array_representation[:3,-1]]

        new_transform = sitk.AffineTransform(3)
        new_transform.SetMatrix(rotation)
        new_transform.SetTranslation(translation)

        new_transform.WriteTransform(full_file_name)