import numpy as np
import DiSPy.irreps as irr


def print_irreps(iso_sg, iv):
    outputfile = open(iv.image_dir + "/../results/output.out", "a")

    print ("\n-------\n------- Possible irreps of the distortion group:\n-------\n")
    outputfile.write("\n\n-------\n------- Possible irreps of the distortion group:\n-------\n\n")

    PIR = open('PIR_data.txt', "r")
    for line in PIR:
        if iso_sg[0:1].upper() + iso_sg[1:].lower() in line:
            print("Irrep #" + line[1:7].strip() + ": " + line[24:32].strip())
            outputfile.write("Irrep #" + line[1:7].strip() + ": " + line[24:32].strip()+'\n')
    PIR.close()

def get_irreps(path, iv):

    DG_std = path.get_DG_std()

    outputfile = open(iv.image_dir + "/../results/output.out", "a")

    print ("\n\nIrrep. #"+str(iv.irr_num)+" chosen...")
    outputfile.write("\n\nIrrep. #"+str(iv.irr_num)+" chosen...")

    if iv.irr_num != 0:
        irreps = np.zeros((len(DG_std[0]),iv.irr_dim, iv.irr_dim))
        kvec_z = np.zeros(3)
        DG_std_aug = []

        # -- Obtain augmented integer matrices for elements of DG in the std. basis
        irr.pir_data_read()

        for i in range(0,len(DG_std[0])):
            mult = 1.0

            while( any(np.logical_and( (mult*DG_std[1][i])%1.0 > 0.05, (mult*DG_std[1][i])%1.0 < 0.95))):
                mult += 1.0

            std_trans = mult*DG_std[1][i]
            std_mat = DG_std[0][i]*mult
            std_mat[0:3,3] = np.rint(std_trans)
            std_mat = std_mat.astype(int)


            DG_std_aug.append(std_mat)

            # -- Obtain irrep matrices
            irreps[i,:,:]=irr.pir_data_get_irmatrix(iv.irr_num,kvec_z,std_mat,irreps[i,:,:],iv.irr_dim)


        print("\n\n")
        outputfile.write("\n\n")

        print ("------- Irrep. Matrices:\n")
        outputfile.write("------- Irrep. Matrices:\n\n")
        for j in range(0,len(irreps)):
            print ("Symmetry Element "+str(j+1)+" : \n"+str(irreps[j])+"\n")
            outputfile.write("Symmetry Element "+str(j+1)+" : \n"+str(irreps[j])+"\n\n")

    return irreps
