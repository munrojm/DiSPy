import numpy as np
from DiSPy.core.irreps import IrrepTools

# -- Function to obtain the irrep. matrices for a DG and irrep.


def get_irrep_matrices(path, io):

    DG_std = path.distortion_group.std_matrices

    io.print("\n\nIrrep. #" + str(io.irr_num) + " chosen...")

    if io.irr_num != 0:
        irrep_tools = IrrepTools()
        irrep = irrep_tools.get_irrep(
            group_data=DG_std, stokes_number=io.irr_num, dimension=io.irr_dim, k_params=(None, None, None)
        )

    else:
        raise ValueError("No irrep chosen.")

    io.print("\n\n")

    io.print("------- Irrep. Matrices:\n")

    for j in range(0, len(irrep)):
        irrep_matrix = irrep.irrep_matrices[j].matrix
        io.print("Symmetry Element " + str(j + 1) + " : \n" + str(irrep_matrix) + "\n")

    return irrep
