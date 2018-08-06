import spglib
import numpy as np
from DiSPy.vec_util import *
from DiSPy.op_id_utils import *

## -- Get the elements of the distortion group
def get_DG(path,iv):

    images = path.get_images()

    numIm = len(images)

    outputfile = open(iv.image_dir + "/../results/output.out", "a")

    print("\n-------\n-------\n------- Symmetry of images:\n-------\n-------\n")
    outputfile.write("\n-------\n-------\n------- Symmetry of images:\n-------\n-------\n\n")

    # -- Gets symmetry data for each image
    Dataset = []
    for image in images:
        image.wrap()
        print("Image "+ str(images.index(image)+1)+": "+str(spglib.get_spacegroup(image,symprec=iv.symprec,angle_tolerance=iv.angtol)))
 
        outputfile.write("Image "+ str(images.index(image)+1)+": "+str(spglib.get_spacegroup(image, symprec=iv.symprec, angle_tolerance=iv.angtol)+'\n'))

        Dataset.append(spglib.get_symmetry_dataset(image, symprec=iv.symprec, angle_tolerance=iv.angtol))

    for i in range(0, iv.numIm):
        for j in range(0, len(Dataset[i]['translations'])):
            Dataset[i]['translations'][j] = standardize(Dataset[i]['translations'][j],iv.gentol)

    


    # -- Finds H, the symmetry operations common for each image
    # H[0] is the list of rotations, H[1] is the list of translations
    H = []
    H.append([])
    H.append([])
    rotH = Dataset[0]['rotations']
    tranH = Dataset[0]['translations']

    for i in range (0,len(rotH)):
        add = True
        for Data in Dataset:
            found = False
            for j in range (0,len(Data['rotations'])):
                if np.allclose(Data['rotations'][j],rotH[i],atol=iv.gentol,rtol=0.0) and closewrapped(Data['translations'][j],tranH[i],iv.vectol2):
                    found = True
                    break
            if not found:
                add = False
                break
        if add:
            H[0].append(rotH[i])
            H[1].append(tranH[i])


    # -- Finds A*, the starred symmetry operations that map all images to their opposite image
    # Astar[0] is the list of rotations, Astar[1] is the list of translations
    Astar = []
    Astar.append([])
    Astar.append([])
    rotA = Dataset[int(iv.numIm/2)]['rotations']
    tranA = Dataset[int(iv.numIm/2)]['translations']

    for i in range (0,len(rotA)):
        add = True
        for j in range(0,int(iv.numIm/2)):
            positions = images[j].get_scaled_positions()
            positions = np.dot(positions, np.transpose(rotA[i]))
            positions = positions + findtranslation(images[int(iv.numIm/2)],rotA[i],tranA[i],iv.gentol,iv.vectol2,iv.symprec,iv.angtol)
            

            newimage = images[j].copy()
            newimage.set_scaled_positions(positions)
            newimage.wrap()

            if not atomsequal(newimage,images[iv.numIm-1-j],iv.vectol):
                add = False
                break
        if add:
            Astar[0].append(rotA[i])
            Astar[1].append(tranA[i])


    # -- Finds the general distortion group with a direct product of H and A*
    # DG[0] is the list of rotations, DG[1] is the list of translations
    # The first len(H[0]) operations are unstarred operations; the rest are starred.
    DG = []
    DG.append([])
    DG.append([])
    for i in range (0,len(H[0])):
        DG[0].append(H[0][i])
        DG[1].append(H[1][i])

    for i in range (0,len(H[0])):
        for j in range (0,len(Astar[0])):
            add = True
            testrot = np.dot(H[0][i],Astar[0][j])
            testtran = standardize(np.dot(Astar[1][j],np.transpose(H[0][i]))+H[1][i],iv.gentol)
            for k in range (len(H[0]),len(DG[0])):
                if np.allclose(testrot,DG[0][k],atol=iv.gentol) and closewrapped(testtran,DG[1][k],iv.vectol):
                    add = False
                    break
            if add:
                DG[0].append(testrot)
                DG[1].append(testtran)

    numUstar = len(H[0])

    # Output information about operations
    print("\n-------\n------- Elements of distortion group in the basis of the inputted structures:\n-------\n")
    outputfile.write("\n-------\n------- Elements of distortion group in the basis of the inputted structures:\n-------\n\n")

    for i in range(0,len(DG[0])):
        if i < numUstar:
            print("Symmetry Element "+ str(i+1) + " (Unstarred):")
            outputfile.write("Symmetry Element "+ str(i+1) + " (Unstarred):\n")
        elif i >= numUstar:
            print("Symmetry Element "+ str(i+1) + " (Starred):")
            outputfile.write("Symmetry Element "+ str(i+1) + " (Starred):\n")

        op_info = operationAttributes(DG[0][i],DG[1][i],iv.gentol)

        print("Rotation:")
        print(DG[0][i])
        print("Translation:")
        print(DG[1][i])
        print(op_info)
        print("")
        outputfile.write("Rotation:\n"+str(DG[0][i])+"\nTranslation:\n"+str(DG[1][i])+'\n'+op_info+"\n\n")

    outputfile.close()
    
    path.set_DG(DG)
    path.numUstar = numUstar
    path.set_img_sym_data(Dataset)

