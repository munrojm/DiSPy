import numpy as np
from DiSPy.core.irreps import IrrepTools

# -- Function to obtain the irrep. matrices for a DG and irrep.


def get_irrep_matrices(path, io):

    DG_std = path.distortion_group.std_matrices

    if io.irr_num != 0:
        irrep_tools = IrrepTools()
        irrep = irrep_tools.get_irrep(
            group_data=DG_std, stokes_number=io.irr_num, dimension=io.irr_dim, k_params=(None, None, None)
        )

    else:
        raise ValueError("No irrep chosen.")

    return irrep
