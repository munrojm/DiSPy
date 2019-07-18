import spglib
import numpy as np
import pkg_resources


from DiSPy.core.dg import DistortionGroup
from DiSPy.core.opid import OperationIdentification
from DiSPy.core.irreps import IrrepTools

from pymatgen.symmetry.groups import SymmOp

# -- Get the elements of the distortion group
def get_DG(path, io):

    path.distortion_group = DistortionGroup(images=path.images, symprec=io.symprec, angle_tolerance=io.angtol,
                                            gentol=io.gentol, vectol=io.vectol, vectol2=io.vectol2)

    path.img_sym_data = path.distortion_group.get_img_sym_data(images=path.images, symprec=io.symprec,
                                                               angle_tolerance=io.angtol, gentol=io.gentol)

    io.print("\n-------\n-------\n------- Symmetry of images:\n-------\n-------\n")

    for image_num in range(io.numIm):

        io.print("Image " + str(image_num+1)+": " +
                 str(path.img_sym_data[image_num]['international'])+' ('+str(path.img_sym_data[image_num]['number'])+')')

    # Output information about operations
    io.print("\n-------\n------- Elements of distortion group in the basis of the inputted structures:\n-------\n")

    for i in range(len(path.distortion_group)):
        if i < path.distortion_group.num_unstar:
            io.print("Symmetry Element " + str(i+1) + " (Unstarred):")

        elif i >= path.distortion_group.num_unstar:
            io.print("Symmetry Element " + str(i+1) + " (Starred):")

        op_iden = OperationIdentification(
            path.distortion_group.matrices[i], io.gentol)
        io.print("Rotation:\n"+str(path.distortion_group.matrices[i].rotation_matrix)+"\nTranslation:\n"+str(
            path.distortion_group.matrices[i].translation_vector)+'\n'+op_iden.info+"\n")

# -- Get elements of the distortion group in a standard setting.
def get_DG_std(path, io):

    images = path.images
    DG = path.distortion_group
    dataset = path.img_sym_data
    num_unstar = DG.num_unstar

    # -- Obtain transformation to standard basis
    io.print(
        "\n-------\n------- Elements of distortion group in the standard basis:\n-------\n")

    if io.trnum < 0:
        io.print("Using user inputted transformation matrix and origin shift...")

    elif io.trnum == 0:
        io.tr_mat, io.oshift = path.distortion_group.obtain_tmat(
            images=path.images, vectol=io.vectol, symprec=io.symprec, angle_tolerance=io.angtol)
        io.tr_mat = np.around(io.tr_mat, decimals=5)
        io.oshift = np.around(io.oshift, decimals=5)
    else:
        io.tr_mat = np.around(
            dataset[io.trnum-1]['transformation_matrix'], decimals=5)
        io.oshift = np.around(dataset[io.trnum-1]['origin_shift'], decimals=5)

        io.print(
            "Using transformation matrix and origin shift from image "+str(io.trnum)+"...")

    io.print("Transformation matrix and origin shift...\n")
    io.print("Matrix:")
    io.print(str(io.tr_mat))

    io.print("Origin shift:")
    io.print(str(io.oshift))

    io.print("\n\n")

    DG.standardize(
        images=path.images, vectol=io.vectol, symprec=io.symprec, tmat=io.tr_mat, oshift=io.oshift)
    for i in range(len(DG)):
        if i < num_unstar:
            io.print("Symmetry Element " + str(i+1) + " (Unstarred):")

        elif i >= num_unstar:
            io.print("Symmetry Element " + str(i+1) + " (Starred):")

        op_iden = OperationIdentification(DG.std_matrices[i], io.gentol)

        io.print("Rotation:\n"+str(np.rint(
            DG.std_matrices[i].rotation_matrix))+"\nTranslation:\n"+ \
        str(DG.std_matrices[i].translation_vector)+'\n'+op_iden.info+'\n')


# Obtain the distortion group name given the elements of the standardized group
def get_dist_name(path, io):

    DG = path.distortion_group

    (iso_sg_name, iso_sg_num) = DG.get_iso(io.gentol)
    dat = _line_dat(iso_sg_num)
    DG_std = DG.std_matrices

    pos_print = []

    if dat[1] == "single" and DG.num_unstar != len(DG_std):
        _print_ds_name(["e"],iso_sg_name,io)
        return iso_sg_name
    else:
        for i in range(len(DG_std)):

            if i > DG.num_unstar-1:
                op_iden = OperationIdentification(DG_std[i], io.gentol)
                idf = op_iden.info.split()

                try:
                    mat = idf[7]+" "+idf[8]
                except:
                    mat = idf[7]

                if "mirror" in idf:
                    about = _to_xyz(" ".join(idf[10:13]))
                elif "rotation" in idf or "rotoinversion" in idf:
                    about = _to_xyz(" ".join(idf[11:14]))
                else:
                    about = " "

                if str(_to_symbol(mat))+str(about)+" " in " ".join(dat):
                    pos_print.append(dat[dat.index(str(_to_symbol(mat))+str(about))+1])

        _print_ds_name(pos_print,iso_sg_name,io)

        return iso_sg_name

def print_irreps(iso_sg_name,io):
        io.print ("\n-------\n------- Possible irreps of the distortion group:\n-------\n")
        out_text = IrrepTools.possible_irreps(iso_sg_name)

        io.print(out_text)

def _print_ds_name(pos_print,iso_sg_name,io):
    final_name = str(iso_sg_name)

    p_count = 0

    for i in sorted(list(set(pos_print))):
        if i == "e":
            final_name = final_name+"*"
        else:
            final_name = final_name[0:1+(int(i)+int(p_count))]+"*"+final_name[1+(int(i)+int(p_count)):]
            p_count += 1

    io.print ("\n-------\n------- Distortion group:\n-------\n")
    io.print (final_name)

# - Search SPG name dictionary
def _line_dat(spg_num):
    with open(pkg_resources.resource_filename('DiSPy', 'SPG_dict.txt')) as symbol_dict:
        for line in symbol_dict:
            l_split = line.split()
            if l_split[0] == str(spg_num):
                return l_split

def _to_xyz(vec):
    n = {"0 0 1":"z", \
         "0 1 0":"y", \
         "1 0 0":"x", \
         "1 1 1":"xyz", \
         "1 1-bar 0":"xby", \
         "1 1 0":"xy",   \
         "0 1-bar 0":"y"}

    return n[str(vec)]

def _to_symbol(sym):
    n = {"twofold rotation":"2", \
         "threefold rotation":"3", \
         "fourfold rotation":"4", \
         "sixfold rotation":"6", \
         "threefold rotoinversion":"-3", \
         "fourfold rotoinversion":"-4", \
         "sixfold rotoinversion":"-6", \
         "mirror across":"m",  \
         "inversion with":"i", \
         "inversion.":"i"}

    return n[str(sym)]



