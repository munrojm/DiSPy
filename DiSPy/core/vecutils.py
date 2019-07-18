import numpy as np
import spglib

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


# -- Checking whether two structure objects are equivalent in configuration
def atomsequal (atom1,atom2,tolerance):

    pos1 = atom1.frac_coords
    pos2 = atom2.frac_coords
    num1 = atom1.species
    num2 = atom2.species

    for i in range (0,len(num1)):
        found = False
        for j in range (0,len(num1)):
            if closewrapped(pos1[i],pos2[j],tolerance) and num1[i] == num2[j]:
                found = True
                break
        if not found:
            return False
    return True

# -- Checking if two positions are approximately equal given wrapping
def closewrapped (pos1,pos2,tolerance):
    if len(pos1) != len(pos2):
        return False
    for i in range (0,len(pos1)):
        if abs(pos1[i]-pos2[i]) > tolerance[i] and abs(pos1[i]-pos2[i]) < 1.0 - tolerance[i]:
            return False
    return True

# -- Finding the appropriate translation for an image
def findtranslation (image,rotation,translation,gentol,vectol2,symprec,angtol):
    final = translation
    #image.wrap()
    imageData = SpacegroupAnalyzer(image, symprec=symprec, angle_tolerance=angtol).get_symmetry_dataset()
    for i in range (0,len(imageData['rotations'])):
        if np.allclose(rotation,imageData['rotations'][i],atol=gentol) and closewrapped(translation,imageData['translations'][i],vectol2):
            final = imageData['translations'][i]
            break
    return final

def standardize(vector,gentol):
    for i in range(0, len(vector)):
        if abs(np.around(vector[i], 0) - vector[i]) < gentol:
            vector[i] = 0
        else:
            vector[i] = vector[i] - np.floor(vector[i])
    return vector
