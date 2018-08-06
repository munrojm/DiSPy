import os
from ase.io import read, write
import time

def o_images(i_path,p_path,iv):

    images = i_path.get_images()
    images_alt = p_path.get_images()

    outputfile = open(iv.image_dir + "/../results/output.out", "a")

    # -- Ouput new structures to files
    os.mkdir(iv.image_dir + "/../results/output_structures")

    for i in range(0,iv.numIm):
        write(iv.image_dir + "/../results/output_structures/"+str(i),images_alt[i],format=iv.o_format)

    # -- Generates XYZ Movie Files
    xyz1 = open(iv.image_dir + "/../results/images1.xyz", "w")
    for image in images:
        write(iv.image_dir + "/../results/temp", image, format="xyz")
        tempfile = open(iv.image_dir + "/../results/temp", "r")
        for line in tempfile:
            xyz1.write(line)
        tempfile.close()
    xyz1.close()

    if iv.irr_num != 0:
        xyz2 = open(iv.image_dir + "/../results/images2.xyz", "w")
        for image in images_alt:
            write(iv.image_dir + "/../results/temp", image, format="xyz")
            tempfile = open(iv.image_dir + "/../results/temp", "r")
            for line in tempfile:
                xyz2.write(line)
            tempfile.close()
        xyz2.close()

    os.remove(iv.image_dir + "/../results/temp")

    print("\n\n==================================================")
    outputfile.write("\n\n\n==================================================")

    print('\n\nTask completed on ' + time.strftime("%c")+"\n\n")
    outputfile.write('\n\n\nTask completed on ' + time.strftime("%c")+"\n\n\n")


    print("If you found this program useful, please consider citing:\n")
    print("\"J.M. Munro et. al. Discovering minimum energy pathways via distortion symmetry groups. Phys. Rev. B (2018).\"\n")
    print("\"B.K. VanLeeuwen & V. Gopalan. The antisymmetry of distortions. Nat. Commun. (2015).\"\n")
 
    outputfile.write("If you found this program useful, please consider citing:\n\n")
    outputfile.write("\"J.M. Munro et. al. Discovering minimum energy pathways via distortion symmetry groups. Phys. Rev. B (2018).\"\n\n")
    outputfile.write("\"B.K. VanLeeuwen & V. Gopalan. The antisymmetry of distortions. Nat. Commun. (2015).\"\n\n")
 
