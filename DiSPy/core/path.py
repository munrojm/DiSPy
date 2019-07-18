import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.groups import SymmOp

from DiSPy.core.dg import DistortionGroup

from DiSPy.core.vecutils import closewrapped

# -- Path object and its attributes

class Path:
    def __init__(self, num_im: int):
        self._images = [None]*num_im
        self._distortion_group = [None]
        self._distortion_group_std = [None]
        self._img_sym_data = [None]*num_im

    @property
    def images(self):
        return self._images

    @property
    def distortion_group(self):
        return self._distortion_group

    @property
    def img_sym_data(self):
        return self._img_sym_data

    @images.setter
    def images(self, image_list):
        if len(image_list) != len(self._images):
            raise ValueError("Image list has wrong length.")
        elif not all([isinstance(i, Structure) for i in image_list]):
            raise ValueError(
                "Images must be pymaten structures.")
        else:
            self._images = image_list

    @distortion_group.setter
    def distortion_group(self, dg):
        if not isinstance(dg, DistortionGroup):
            raise ValueError(
                "Symmetry operations in group data must be instances of SymmOp.")
        else:
            self._distortion_group = dg

    @img_sym_data.setter
    def img_sym_data(self, img_sym_data_list):
        if len(img_sym_data_list) != len(self._images):
            raise ValueError("Symmetry data list has wrong length.")
        else:
            self._img_sym_data = img_sym_data_list


    def gen_atom_map(self, basis, vectol):

        images = self._images
        DG = self._distortion_group.matrices
        numIm = len(images)
        numAtoms = len(images[0].frac_coords)
        num_unstar = self._distortion_group.num_unstar

        a_map=np.zeros((len(DG),numIm,numAtoms,numAtoms))

        for i in range(len(DG)):
            for j in range(1,numIm-1):
                atoms1 = images[j].frac_coords
                num1 = images[j].species

                for k in range(0,numAtoms):

                    t_coord = np.dot(DG[i].rotation_matrix,atoms1[k])
                    t_coord = (t_coord + DG[i].translation_vector)%1.0

                    if i < num_unstar:
                       atoms2 = images[j].frac_coords
                       num2 = images[j].species
                       t_im = j
                    else:
                       atoms2 = images[numIm-1-j].frac_coords
                       num2 = images[numIm-1-j].species
                       t_im = numIm-1-j

                    for l in range(0,numAtoms):
                        if (closewrapped (t_coord,atoms2[l],vectol) and \
                            num1[k] == num2[l] and basis[k] == 1):
                            a_map[i,j,k,l] = 1

        return a_map
