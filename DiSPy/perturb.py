import numpy as np
import random
from DiSPy.core.vecutils import closewrapped
from DiSPy.core.irreps import IrrepTools


###
# -- Main function to generate perturbation
###
def gen_perturb(path, irrep, io):

    images = path.images
    DG = path.distortion_group.matrices
    num_unstar = path.distortion_group.num_unstar
    numIm = len(images)

    irrep_tools = IrrepTools()

    # -- Generate starting basis
    atoms1 = images[0].frac_coords
    numAtoms = len(atoms1)
    basis = np.zeros(numAtoms)
    atoms2 = images[numIm - 1].frac_coords
    m_vec = np.dot(np.linalg.inv(images[0].lattice.matrix), [io.min_move, io.min_move, io.min_move])

    for i in range(0, numAtoms):
        if not closewrapped(atoms1[i, :], atoms2[i, :], m_vec):
            basis[i] = 1

    # -- Generate matrix showing atom mapping for each operations
    a_map = path.gen_atom_map(basis=basis, vectol=io.vectol)

    # -- Generate and apply modes for an irrep of arbitrary dimension
    pt = np.zeros((numIm, numAtoms, 3))

    pt_mode_init = irrep_tools.projection_diag(
        images=images, symmop_list=DG, irrep=irrep, num_unstar=num_unstar, a_map=a_map, basis=basis
    )
    pt_mode_init = pt_mode_init / np.linalg.norm(np.ndarray.flatten(pt_mode_init))

    pt += io.m_co[0] * pt_mode_init

    if io.irr_dim > 1:
        for i in range(1, io.irr_dim):
            pt_mode = irrep_tools.projection_odiag(
                images=images,
                symmop_list=DG,
                num_unstar=num_unstar,
                irrep=irrep,
                vec=pt_mode_init,
                index=i,
                a_map=a_map,
                basis=basis,
            )

            pt += io.m_co[i] * (pt_mode / np.linalg.norm(np.ndarray.flatten(pt_mode)))

    p_vec = [io.p_mag / np.linalg.norm(images[0].lattice.matrix[m, :]) for m in range(0, 3)]
    # print pt[3]/np.amax(abs(pt[3]))
    # print(np.amax(p_vec*(pt[:,:,:]/np.amax(abs(pt)))))
    images_alt = []

    for i in range(numIm):
        image_copy = images[i].copy()
        alt_species = image_copy.species

        perturbed_coords = image_copy.frac_coords + p_vec * (pt[i, :, :] / np.amax(abs(pt)))
        for j in range(numAtoms):
            image_copy.replace(j, alt_species[j], perturbed_coords[j])

        images_alt.append(image_copy)

    return images_alt, basis
