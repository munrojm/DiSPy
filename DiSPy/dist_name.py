import numpy as np
import spglib
from DiSPy import op_id_utils
import pkg_resources

# Obtain the distortion group name given the elements of the group
def get_dist_name(path, iv):

    (iso_sg_name, iso_sg_num) = get_iso(path, iv)
    dat = line_dat(iso_sg_num)
    DG_std = path.get_DG_std()

    pos_print = []

    if dat[1] == "single" and path.numUstar != len(DG_std[0]):
        print_ds_name(["e"],iso_sg_name,iv)
        return iso_sg_name
    else:
        for i in range(0, len(DG_std[0])):

            if i > path.numUstar-1:
                id_line = op_id_utils.operationAttributes(DG_std[0][i][0:3,0:3],DG_std[1][i], iv.gentol)
                idf = id_line.split()

                try:
                    mat = idf[7]+" "+idf[8]
                except:
                    mat = idf[7]

                if "mirror" in idf:
                    about = to_xyz(" ".join(idf[10:13]))
                elif "rotation" in idf or "rotoinversion" in idf:
                    about = to_xyz(" ".join(idf[11:14]))
                else:
                    about = " "

                if str(to_symbol(mat))+str(about)+" " in " ".join(dat):
                    pos_print.append(dat[dat.index(str(to_symbol(mat))+str(about))+1])

        print_ds_name(pos_print,iso_sg_name,iv)

        return iso_sg_name

# - Output name
def print_ds_name(pos_print,iso_sg_name,iv):
    outputfile = open(iv.image_dir + "/../results/output.out", "a")
    final_name = str(iso_sg_name)

    p_count = 0

    for i in sorted(list(set(pos_print))):
        if i == "e":
            final_name = final_name+"*"
        else:
            final_name = final_name[0:1+(int(i)+int(p_count))]+"*"+final_name[1+(int(i)+int(p_count)):]
            p_count += 1

    print ("\n-------\n------- Distortion group:\n-------\n")
    print (final_name)

    outputfile.write("\n-------\n------- Distortion group:\n-------\n\n")
    outputfile.write(final_name)

# - Search SPG name dictionary
def line_dat(spg_num):
    with open(pkg_resources.resource_filename('DiSPy', 'SPG_dict.txt')) as symbol_dict:
        for line in symbol_dict:
            l_split = line.split()
            if l_split[0] == str(spg_num):
                return l_split

def to_xyz(vec):
    n = {"0 0 1":"z", \
         "0 1 0":"y", \
         "1 0 0":"x", \
         "1 1 1":"xyz", \
         "1 1-bar 0":"xby", \
         "1 1 0":"xy",   \
         "0 1-bar 0":"y"}

    return n[str(vec)]

def to_symbol(sym):
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

# -- Obtain isomorphic spacegroup of distortion group
def get_iso(path, iv):

    DG = path.get_DG()

    outputfile = open(iv.image_dir + "/../results/output.out", "a")

    h_number = spglib.get_hall_number_from_symmetry(np.around(DG[0][:],decimals=4),np.around(DG[1][:],decimals=6),symprec=iv.gentol)
    sg_type_data = spglib.get_spacegroup_type(h_number)
    iso_sg = sg_type_data['international_short']
    iso_sg_num = sg_type_data['number']


    return iso_sg, iso_sg_num
