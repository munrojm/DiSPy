import numpy as np
import spglib
from DiSPy.vec_util import closewrapped
from DiSPy.vec_util import standardize
from DiSPy.op_id_utils import operationAttributes

# -- Function to obtain the transformation matrix and origin shift to std. setting
def obtain_tmat(DG,images,iv):
    # -- Generate lattice with DG in input basis

    cp_pos = images[0].get_scaled_positions()
    cp_lattice = images[int(iv.numIm/2)].get_cell()
    cp_labels = images[0].get_atomic_numbers()

    for i in [int(iv.numIm/2),iv.numIm-1]:

        cp_pos = np.append(cp_pos,images[i].get_scaled_positions(),axis=0)
        cp_labels = np.append(cp_labels,images[i].get_atomic_numbers())

    new_pos = [cp_pos[0]]
    new_labels = [cp_labels[0]]

    for j in range(0,len(cp_pos[:,0])):
        atom = cp_pos[j]
        for i in range(0, len(DG[0])):
            t_coord = np.dot(atom, np.transpose(DG[0][i]))
            t_coord = (t_coord + DG[1][i])%1.0

            if not any([closewrapped (t_coord,m,iv.vectol) for m in new_pos]):
                new_pos = np.append(new_pos,[t_coord],axis=0)
                new_labels = np.append(new_labels,cp_labels[j])


    dataset2 = spglib.get_symmetry_dataset((cp_lattice,new_pos,new_labels),symprec=iv.symprec)
    return (dataset2['transformation_matrix'], dataset2['origin_shift'])

def get_DG_std(path, iv):

    images = path.get_images()
    DG = path.get_DG()
    Dataset = path.get_img_sym_data()
    numUstar = path.numUstar

    outputfile = open(iv.image_dir + "/../results/output.out", "a")

    iv.numIm = len(images)

    # -- Obtain transformation to standard basis
    print("\n-------\n------- Elements of distortion group in the standard basis:\n-------\n")
    outputfile.write("\n-------\n------- Elements of distortion group in the standard basis:\n-------\n\n")

    if iv.trnum < 0:
        print("Using user inputted transformation matrix and origin shift...")
        outputfile.write("Using user inputted transformation matrix and origin shift...\n")

    elif iv.trnum == 0:
        iv.tr_mat, iv.oshift = obtain_tmat(DG,images,iv)
        iv.tr_mat = np.around(iv.tr_mat,decimals=5)
        iv.oshift = np.around(iv.oshift,decimals=5)
    else:
        iv.tr_mat = np.around(Dataset[iv.trnum-1]['transformation_matrix'],decimals=5)
        iv.oshift = np.around(Dataset[iv.trnum-1]['origin_shift'],decimals=5)

        print ("Using transformation matrix and origin shift from image "+str(iv.trnum)+"...")
        outputfile.write("Using transformation matrix and origin shift from image "+str(iv.trnum)+"...\n")


    print ("Transformation matrix and origin shift...\n")
    print ("Matrix:")
    outputfile.write("Transformation matrix and origin shift...\n\n")
    outputfile.write("Matrix:\n")
    print(iv.tr_mat)
    outputfile.write(str(iv.tr_mat)+"\n")

    print ("Origin shift:")
    outputfile.write("Origin shift:\n")
    print (iv.oshift)
    outputfile.write(str(iv.oshift)+"\n")



    print ("\n\n")
    outputfile.write("\n\n")

    DG_std = []
    DG_std.append([])
    DG_std.append([])

    for i in range(0,len(DG[0])):
        std_mat = np.zeros((4,4))
        std_trans = np.zeros(3)

        std_mat[0:3,0:3] = np.dot(np.dot(iv.tr_mat, DG[0][i]), np.linalg.inv(iv.tr_mat))
        std_mat[3,3] = 1


        std_trans = np.dot(iv.tr_mat, DG[1][i]) + iv.oshift - \
                     np.dot(np.dot(np.dot(iv.tr_mat, DG[0][i]),np.linalg.inv(iv.tr_mat)), iv.oshift)


        DG_std[0].append(std_mat)
        DG_std[1].append(std_trans)


        if i < numUstar:
            print("Symmetry Element "+ str(i+1) + " (Unstarred):")
            outputfile.write("Symmetry Element "+ str(i+1) + " (Unstarred):\n")
        elif i >= numUstar:
            print("Symmetry Element "+ str(i+1) + " (Starred):")
            outputfile.write("Symmetry Element "+ str(i+1) + " (Starred):\n")

        op_info = operationAttributes(np.rint(std_mat[0:3,0:3]),std_trans[:],iv.gentol)
        print("Rotation:")
        print(np.rint(std_mat[0:3,0:3]))
        print("Translation:")
        print(std_trans[:])
        print(op_info)
        print("")
        outputfile.write("Rotation:\n"+str(np.rint(std_mat[0:3,0:3]))+"\nTranslation:\n"+str(std_trans[:])+'\n'+op_info+'\n\n')

    outputfile.close()
    
    path.set_DG_std(DG_std)

