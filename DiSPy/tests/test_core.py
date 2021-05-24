import unittest

from DiSPy.core.irreps import IrrepTools

from pymatgen.symmetry.groups import SymmOp

import numpy as np


def test_get_irrep_mat():
    stokes_number = 7611
    dimension = 3
    symmop = SymmOp([[0, 3, 0, -1], [-3, 3, 0, 1], [0, 0, -3, 1], [0, 0, 0, 3]])

    irr = IrrepTools()
    irrep_mat = irr.get_irrep_mat(stokes_number=stokes_number, dimension=dimension, symmop=symmop)

    assert irrep_mat.matrix.tolist() == [[0, 0, -1], [-1, 0, 0], [0, 1, 0]]
    assert irrep_mat.label == "F1+"


def test_get_irrep():
    stokes_number = 2857
    dimension = 1

    matrices = [
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
        [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
        [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]],
        [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
        [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
        [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
        [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
        [[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]],
        [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]],
    ]

    vectors = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.5],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, -0.5, 0.0],
        [1.4999999959999997, -0.4999999999999998, 0.5],
        [1.4999999959999997, -0.5, 0.0],
        [0.5, -0.4999999999999998, 0.5],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 0.4999868850768938],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.4999868850768938],
        [1.4999998980000022, -0.5, 0.0],
        [0.5, -0.5, 0.4999868850768938],
        [0.5, -0.5, 0.0],
        [1.4999998980000022, -0.5, 0.4999868850768938],
    ]

    group_data = []
    for entry in range(len(matrices)):
        group_data.append(
            SymmOp.from_rotation_and_translation(rotation_matrix=matrices[entry], translation_vec=vectors[entry])
        )

    ir_mat_data = [
        [[1]],
        [[-1]],
        [[-1]],
        [[1]],
        [[-1]],
        [[1]],
        [[1]],
        [[-1]],
        [[1]],
        [[-1]],
        [[-1]],
        [[1]],
        [[-1]],
        [[1]],
        [[1]],
        [[-1]],
    ]

    irr = IrrepTools()
    irrep = irr.get_irrep(
        group_data=group_data, stokes_number=stokes_number, dimension=dimension, k_params=(None, None, None)
    )

    irmats = irrep.irrep_matrices
    for op_num in range(len(irmats)):
        assert irmats[op_num].matrix == ir_mat_data[op_num]

    irrep.label == "Y4+"


# if __name__ == "__main__":
#     unittest.main()
