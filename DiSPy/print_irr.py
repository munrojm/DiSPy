import numpy as np
import spglib


# -- Obtain and print isomorphic spacegroup of distortion group
def get_iso(DG, iv):
    outputfile = open(iv.image_dir + "/../results/output.out", "a")

    print ("\n-------\n------- Possible irreps of the distortion group:\n-------\n")
    outputfile.write("\n-------\n------- Possible irreps of the distortion group:\n-------\n\n")



    PIR = open("PIR_data.txt", "r")
    for line in PIR:
        if iso_sg[0:1].upper() + iso_sg[1:].lower() in line:
            if printresults:
                print("Irrep #" + line[1:7].strip() + ": " + line[24:32].strip())
            outputfile.write("Irrep #" + line[1:7].strip() + ": " + line[24:32].strip()+'\n')
    PIR.close()
