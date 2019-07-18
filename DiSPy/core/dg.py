from DiSPy.core.vecutils import *

from pymatgen.symmetry.groups import SymmOp
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from dataclasses import dataclass
from typing import List

import spglib

import numpy as np


@dataclass
class DistortionGroup:

    matrices: List
    std_matrices: List
    trans_mat: List
    oshift: List
    label: str
    num_unstar: int
    img_sym_dat: List

    def __init__(self, images: List, symprec: float, angle_tolerance: float, gentol: float,
                 vectol: List[float], vectol2: List[float]):

        dataset = self.get_img_sym_data(
            images, symprec, angle_tolerance, gentol)

        self.img_sym_dat = dataset

        # -- Finds H, the symmetry operations common for each image
        # H[0] is the list of rotations, H[1] is the list of translations

        numIm = len(images)

        H = [[], []]
        rotH = dataset[0]['rotations']
        tranH = dataset[0]['translations']

        for i in range(0, len(rotH)):
            add = True
            for data in dataset:
                found = False
                for j in range(0, len(data['rotations'])):
                    if np.allclose(data['rotations'][j], rotH[i], atol=gentol, rtol=0.0) and closewrapped(data['translations'][j], tranH[i], vectol2):
                        found = True
                        break
                if not found:
                    add = False
                    break
            if add:
                H[0].append(rotH[i])
                H[1].append(tranH[i])

        # -- Finds A*, the starred symmetry operations that map all images to their opposite image
        # Astar[0] is the list of rotations, Astar[1] is the list of translations
        Astar = [[], []]
        rotA = dataset[int(numIm/2)]['rotations']
        tranA = dataset[int(numIm/2)]['translations']

        for i in range(0, len(rotA)):
            add = True
            for j in range(0, int(numIm/2)):
                positions = images[j].frac_coords
                positions = np.dot(positions, np.transpose(rotA[i]))
                positions = positions + findtranslation(images[int(
                    numIm/2)], rotA[i], tranA[i], gentol, vectol2, symprec, angle_tolerance)

                newimage = images[j].copy()
                new_species = newimage.species
                for k in range(len(new_species)):
                    newimage.replace(k, new_species[k], positions[k])

                if not atomsequal(newimage, images[numIm-1-j], vectol):
                    add = False
                    break
            if add:
                Astar[0].append(rotA[i])
                Astar[1].append(tranA[i])

        # -- Finds the general distortion group with a direct product of H and A*
        # DG[0] is the list of rotations, DG[1] is the list of translations
        # The first len(H[0]) operations are unstarred operations; the rest are starred.
        DG = []
        for i in range(0, len(H[0])):
            DG.append(SymmOp.from_rotation_and_translation(H[0][i], H[1][i]))

        for i in range(0, len(H[0])):
            for j in range(0, len(Astar[0])):
                add = True
                testrot = np.dot(H[0][i], Astar[0][j])
                testtran = standardize(
                    np.dot(Astar[1][j], np.transpose(H[0][i]))+H[1][i], gentol)
                for k in range(len(H[0]), len(DG)):
                    if np.allclose(testrot, DG[k].rotation_matrix, atol=gentol) and \
                            closewrapped(testtran, DG[k].translation_vector, vectol):
                        add = False
                        break
                if add:
                    DG.append(SymmOp.from_rotation_and_translation(
                        testrot, testtran))

        self.num_unstar = len(H[0])
        self.matrices = DG

    def __len__(self):
        return len(self.matrices)

    def __iter__(self):
        return self.matrices

    def obtain_tmat(self, images, vectol, symprec, angle_tolerance):
        # -- Generate lattice with DG in input basis

        numIm = len(images)
        DG = self.matrices

        cp_pos = images[0].frac_coords
        cp_lattice = images[int(numIm/2)].lattice.matrix
        cp_labels = images[0].species
        newimage = images[int(numIm/2)].copy()

        for i in [int(numIm/2), numIm-1]:

            cp_pos = np.append(
                cp_pos, images[i].frac_coords, axis=0)
            cp_labels = np.append(cp_labels, images[i].species)

            newimage.remove_sites(range(len(newimage.species)))

        new_pos = [cp_pos[0]]
        new_labels = [cp_labels[0]]

        for j in range(len(cp_pos[:, 0])):
            atom = cp_pos[j]
            for i in range(len(DG)):
                t_coord = np.dot(atom, np.transpose(DG[i].rotation_matrix))
                t_coord = (t_coord + DG[i].translation_vector) % 1.0

                if not any([closewrapped(t_coord, m, vectol) for m in new_pos]):
                    new_pos = np.append(new_pos, [t_coord], axis=0)
                    new_labels = np.append(new_labels, cp_labels[j])

        for k in range(len(new_labels)):
            newimage.append(new_labels[k], new_pos[k])

        dataset2 = SpacegroupAnalyzer(
            newimage, symprec=symprec, angle_tolerance=angle_tolerance).get_symmetry_dataset()
        return dataset2['transformation_matrix'], dataset2['origin_shift']

    @staticmethod
    def get_img_sym_data(images, symprec, angle_tolerance, gentol):
        # -- Gets symmetry data for each images

        numIm = len(images)
        dataset = []
        for image in images:
            dataset.append(SpacegroupAnalyzer(
                image, symprec=symprec, angle_tolerance=angle_tolerance).get_symmetry_dataset())

        for i in range(0, numIm):
            for j in range(0, len(dataset[i]['translations'])):
                dataset[i]['translations'][j] = standardize(
                    dataset[i]['translations'][j], gentol)

        return dataset

    def standardize(self, images, vectol, symprec, tmat, oshift):
        DG_std = []
        DG = self.matrices

        for i in range(len(DG)):

            std_mat = np.dot(
                np.dot(tmat, DG[i].rotation_matrix), np.linalg.inv(tmat))

            std_trans = np.dot(tmat, DG[i].translation_vector) + oshift - \
                np.dot(
                    np.dot(np.dot(tmat, DG[i].rotation_matrix), np.linalg.inv(tmat)), oshift)

            DG_std.append(
                SymmOp.from_rotation_and_translation(std_mat, std_trans))

        self.std_matrices = DG_std
        self.trans_mat = tmat
        self.oshift = oshift

    # -- Obtain isomorphic spacegroup of distortion group
    def get_iso(self, symprec):

        DG = self.matrices

        mat_list = [mat.rotation_matrix for mat in DG]
        vec_list = [mat.translation_vector for mat in DG]

        h_number = spglib.get_hall_number_from_symmetry(np.around(
            mat_list, decimals=4), np.around(vec_list, decimals=6), symprec = symprec)
        sg_type_data=spglib.get_spacegroup_type(h_number)
        iso_sg=sg_type_data['international_short']
        iso_sg_num=sg_type_data['number']

        return iso_sg, iso_sg_num
