import DiSPy.core.stokes as stokes

from pymatgen.symmetry.groups import SymmOp
from pymatgen.core.structure import Structure
from dataclasses import dataclass
from typing import List

import numpy as np

import pkg_resources


@dataclass(frozen=True)
class IrrepMat:

    stokes_number: int
    label: str
    k_params: List[float]
    dimension: int
    symmop: List
    matrix: List


@dataclass
class Irrep:

    irrep_matrices: List
    group_data: List
    label: str
    stokes_number: int
    dimension: int

    def __len__(self):
        return len(self.irrep_matrices)


class IrrepTools:
    def __init__(self):

        stokes.pir_data_read(pkg_resources.resource_filename("DiSPy", "PIR_data.txt"))

    def get_irrep_mat(self, stokes_number: int, dimension: int, symmop: SymmOp, k_params=(None, None, None)):
        """

        Args:
            stokes_number: Unique number of the irrep as defined by ISO-IR (Stoked et al.)
            dimension: Dimesion of the irrep
            symmop: Instance of SymmOp defining the matrix representation of
            an element in the group in the standard setting (has
            to have all integer elements)
            k_params: Arbitrary for high symmetry k-point lookup

        Returns: An instance of IrrepMat containing matrix and irrep information

        """

        # -- Irrep matrix
        matrix = stokes.pir_data_get_irmatrix(
            stokes_number, k_params, symmop.affine_matrix, np.zeros((dimension, dimension)), dimension
        )

        # -- Irrep label
        PIR = open(pkg_resources.resource_filename("DiSPy", "PIR_data.txt"), "r")
        for line in PIR:
            if str(stokes_number) in line:
                label = line[24:32].strip()
        PIR.close()

        return IrrepMat(
            stokes_number=stokes_number,
            dimension=dimension,
            k_params=k_params,
            symmop=symmop,
            matrix=matrix,
            label=label,
        )

    def get_irrep(self, group_data: List[SymmOp], stokes_number: int, dimension: int, k_params=(None, None, None)):
        """

        Args:
            stokes_number: Unique number of the irrep as defined by ISO-IR (Stoked et al.)
            dimension: Dimension of the irrep
            k_params: Arbitrary for high symmetry k-point lookup
            group_data: List of matrix vector pairs in the symmetry group
            as instances of SymmOp

        Returns: An instance of Irrep containing all matrix and irrep information

        """

        if not all([isinstance(i, SymmOp) for i in group_data]):
            raise ValueError("Matrix-vector pairs in group data must be instances of SymmOp.")

        irrep_matrices = []
        for op in group_data:
            mult = 1.0
            tvec = op.translation_vector

            while any(np.logical_and((mult * tvec) % 1.0 > 0.05, (mult * tvec) % 1.0 < 0.95)):
                mult += 1.0

            std_mat = op.affine_matrix * mult
            std_trans = mult * op.translation_vector
            std_mat[0:3, 3] = np.rint(std_trans)
            std_mat = SymmOp(std_mat.astype(int))

            irrep_matrices.append(
                self.get_irrep_mat(stokes_number=stokes_number, dimension=dimension, k_params=k_params, symmop=std_mat)
            )

        return Irrep(
            irrep_matrices=irrep_matrices,
            group_data=group_data,
            label=irrep_matrices[0].label,
            stokes_number=stokes_number,
            dimension=dimension,
        )

    @staticmethod
    def possible_irreps(iso_sg):
        """

        Args:
            iso_sg: Symbol of spacegroup isomorphic to distortion group

        Returns: Text listing possible irreps

        """

        PIR = open(pkg_resources.resource_filename("DiSPy", "PIR_data.txt"), "r")
        formatted_iso_sg = iso_sg.replace("*", "")
        out_str = ""
        for line in PIR:
            if formatted_iso_sg[0:1].upper() + formatted_iso_sg[1:].lower() in line:
                out_str = out_str + "Irrep #" + line[1:7].strip() + ": " + line[24:32].strip() + "\n"
        PIR.close()
        return out_str

    def projection_diag(
        self,
        images: List[Structure],
        symmop_list: List[SymmOp],
        irrep: Irrep,
        num_unstar: int,
        a_map: List,
        basis: List,
    ):
        """

        Args:
            images: List of pymatgen structures
            symmop_list: List of SymmOp instances used in projection
            irrep: Instance of Irrep for a given group
            num_unstar: Number of unstarred ops. in distortion group.
            This can be set to the length of the group is perojection using
            a regular group is required.
            a_map: Matrix providing atom mappings for each symmetry operation
            in symmop_list.
            basis: List of atoms to include in the basis.

        Returns: Vector after projection onto symmetry invariant subspace associated
        with the chosen irrep using first diagonal entry of irrep matrices.

        """

        numIm = len(images)
        numAtoms = len(images[0].frac_coords)
        irrep_matrices = irrep.irrep_matrices

        mode = np.zeros((numIm, numAtoms, 3))
        vec_list = np.zeros((1, numIm, numAtoms, 3))

        for i in range(1, numIm - 1):
            for j in range(0, numAtoms):
                if basis[j] == 1:
                    for m in range(0, 3):
                        temp_mode = np.zeros((1, numIm, numAtoms, 3))
                        for l in range(len(symmop_list)):

                            if l < num_unstar:
                                t_im = i
                            else:
                                t_im = numIm - 1 - i

                            temp_mode[0, t_im, np.nonzero(a_map[l, t_im, j, :])[0][0], :] += (  # type: ignore
                                irrep_matrices[l].matrix[0, 0] * symmop_list[l].rotation_matrix[:, m]
                            )

                        if not np.any(
                            [
                                np.sum(np.multiply(temp_mode[0, :, :, :], vec_list[t, :, :, :]))
                                for t in range(0, len(vec_list[:, 0, 0, 0]))
                            ]
                        ) and np.any(temp_mode):

                            vec_list = np.concatenate((vec_list, temp_mode))

        # if iv.r_co:
        #     for i in range(0, len(vec_list[:, 0, 0, 0])):
        #         mode += random.uniform(-1.0, 1.0)*vec_list[i, :, :, :]
        # else:

        mode = np.sum(vec_list, axis=0)  # type: ignore

        return mode

    def projection_odiag(
        self,
        images: List[Structure],
        symmop_list: List[SymmOp],
        num_unstar: int,
        irrep: Irrep,
        vec: List,
        index: int,
        a_map: List,
        basis: List,
    ):
        """

        Args:
            images: List of pymatgen structures
            symmop_list: List of SymmOp instances used in projection
            irrep: Instance of Irrep for a given group
            vec: List representing vector obtained from original projection with diagonal element
            index: Off diagonal matrix index
            num_unstar: Number of unstarred ops. in distortion group.
            This can be set to the length of the group is perojection using
            a regular group is required.
            a_map: Matrix providing atom mappings for each symmetry operation
            in symmop_list.
            basis: List of atoms to include in the basis.

        Returns: Vector after projection onto symmetry invariant subspace associated
        with the chosen irrep using off-diagonal elements of irrep matrices.

        """

        numIm = len(images)
        numAtoms = len(images[0].frac_coords)
        irrep_matrices = irrep.irrep_matrices

        mode = np.zeros((numIm, numAtoms, 3))
        vec_list = np.zeros((1, numIm, numAtoms, 3))

        for i in range(1, numIm - 1):
            for j in range(0, numAtoms):
                if basis[j] == 1 and any([vec[i, j, p] != 0 for p in range(0, 3)]):  # type: ignore
                    for m in range(0, 3):
                        temp_mode = np.zeros((1, numIm, numAtoms, 3))
                        for l in range(0, len(symmop_list)):

                            if l < num_unstar:
                                t_im = i
                            else:
                                t_im = numIm - 1 - i

                            n_coord = np.dot(symmop_list[l].rotation_matrix, vec[i, j, :])  # type: ignore

                            temp_mode[0, t_im, np.nonzero(a_map[l, t_im, j, :])[0][0], :] += (  # type: ignore
                                irrep_matrices[l].matrix[index, 0] * n_coord
                            )

                        vec_list = np.concatenate((vec_list, temp_mode))

        mode = np.sum(vec_list, axis=0)  # type: ignore

        return mode
