import numpy as np
import spglib


# -- Obtain and print isomorphic spacegroup of distortion group
def get_iso(path, iv):

    DG = path.get_DG()

    outputfile = open(iv.image_dir + "/../results/output.out", "a")

    h_number = spglib.get_hall_number_from_symmetry(np.around(DG[0][:],decimals=4),np.around(DG[1][:],decimals=6),symprec=iv.gentol)
    sg_type_data = spglib.get_spacegroup_type(h_number)
    iso_sg = sg_type_data['international_short']

    print ("\n-------\n------- Isomorphic space group of distortion group:\n-------\n")
    print (iso_sg)

    outputfile.write("\n-------\n------- Isomorphic space group of distortion group:\n-------\n\n")
    outputfile.write(iso_sg)

    return iso_sg
