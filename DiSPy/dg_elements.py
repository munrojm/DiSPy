import numpy as np
import pkg_resources


from DiSPy.core.dg import DistortionGroup
from DiSPy.core.opid import OperationIdentification
from DiSPy.core.irreps import IrrepTools


# -- Get the elements of the distortion group
def get_DG(path, io):

    path.distortion_group = DistortionGroup.from_images(
        images=path.images,
        symprec=io.symprec,
        angle_tolerance=io.angtol,
        gentol=io.gentol,
        vectol=io.vectol,
        vectol2=io.vectol2,
    )

    path.img_sym_data = path.distortion_group.get_img_sym_data(
        images=path.images, symprec=io.symprec, angle_tolerance=io.angtol, gentol=io.gentol
    )


# -- Get elements of the distortion group in a standard setting.
def get_DG_std(path, io):

    DG = path.distortion_group
    dataset = path.img_sym_data
    try:
        num_unstar = DG.num_unstar
    except AttributeError:
        raise RuntimeError("Must obtain distortion group in crystal basis before standard basis.")

    if io.trnum == 0:
        io.tr_mat, io.oshift = path.distortion_group.obtain_tmat(
            images=path.images, vectol=io.vectol, symprec=io.symprec, angle_tolerance=io.angtol
        )
        io.tr_mat = np.around(io.tr_mat, decimals=5)
        io.oshift = np.around(io.oshift, decimals=5)
    elif io.trnum > 0:
        io.tr_mat = np.around(dataset[io.trnum - 1]["transformation_matrix"], decimals=5)
        io.oshift = np.around(dataset[io.trnum - 1]["origin_shift"], decimals=5)

    DG.standardize(images=path.images, vectol=io.vectol, symprec=io.symprec, tmat=io.tr_mat, oshift=io.oshift)

    path.distortion_group_std = DG


# Obtain the distortion group name given the elements of the standardized group
def get_dist_name(path, io):

    DG = path.distortion_group

    (iso_sg_name, iso_sg_num) = DG.get_iso(io.gentol)
    dat = _line_dat(iso_sg_num)
    DG_std = DG.std_matrices

    pos_print = []

    if dat[1] == "single" and DG.num_unstar != len(DG_std):
        _print_ds_name(["e"], iso_sg_name, io)
        return iso_sg_name
    else:
        for i in range(len(DG_std)):

            if i > DG.num_unstar - 1:
                op_iden = OperationIdentification(DG_std[i], io.gentol)
                idf = op_iden.info.split()

                try:
                    mat = idf[7] + " " + idf[8]
                except IndexError:
                    mat = idf[7]

                if "mirror" in idf:
                    about = _to_xyz(" ".join(idf[10:13]))
                elif "rotation" in idf or "rotoinversion" in idf:
                    about = _to_xyz(" ".join(idf[11:14]))
                else:
                    about = " "

                if str(_to_symbol(mat)) + str(about) + " " in " ".join(dat):
                    pos_print.append(dat[dat.index(str(_to_symbol(mat)) + str(about)) + 1])

        final_name = str(iso_sg_name)

        p_count = 0

        for i in sorted(list(set(pos_print))):
            if i == "e":
                final_name = final_name + "*"
            else:
                final_name = (
                    final_name[0 : 1 + (int(i) + int(p_count))] + "*" + final_name[1 + (int(i) + int(p_count)) :]
                )
                p_count += 1

        return final_name


# - Search SPG name dictionary
def _line_dat(spg_num):
    with open(pkg_resources.resource_filename("DiSPy", "SPG_dict.txt")) as symbol_dict:
        for line in symbol_dict:
            l_split = line.split()
            if l_split[0] == str(spg_num):
                return l_split


def _to_xyz(vec):
    n = {"0 0 1": "z", "0 1 0": "y", "1 0 0": "x", "1 1 1": "xyz", "1 1-bar 0": "xby", "1 1 0": "xy", "0 1-bar 0": "y"}

    return n[str(vec)]


def _to_symbol(sym):
    n = {
        "twofold rotation": "2",
        "threefold rotation": "3",
        "fourfold rotation": "4",
        "sixfold rotation": "6",
        "threefold rotoinversion": "-3",
        "fourfold rotoinversion": "-4",
        "sixfold rotoinversion": "-6",
        "mirror across": "m",
        "inversion with": "i",
        "inversion.": "i",
    }

    return n[str(sym)]
