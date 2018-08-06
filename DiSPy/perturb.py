import numpy as np
import random
from DiSPy.vec_util import closewrapped

# -- Function to obtain mapping matrix for atoms
def gen_atom_map(DG, images, numUstar, basis, iv):
    a_map=np.zeros((len(DG[0]),iv.numIm,iv.numAtoms,iv.numAtoms))

    for i in range(0,len(DG[0])):
        for j in range(1,iv.numIm-1):
            atoms1 = images[j].get_scaled_positions()
            num1 = images[j].get_atomic_numbers()

            for k in range(0,iv.numAtoms):

                t_coord = np.dot(DG[0][i],atoms1[k])
                t_coord = (t_coord + DG[1][i])%1.0

                if i < numUstar:
                   atoms2 = images[j].get_scaled_positions()
                   num2 = images[j].get_atomic_numbers()
                   t_im = j
                else:
                   atoms2 = images[iv.numIm-1-j].get_scaled_positions()
                   num2 = images[iv.numIm-1-j].get_atomic_numbers()
                   t_im = iv.numIm-1-j

                for l in range(0,iv.numAtoms):
                    if (closewrapped (t_coord,atoms2[l],iv.vectol) and \
                        num1[k] == num2[l] and basis[k] == 1):
                        a_map[i,j,k,l] = 1



    return a_map

# -- Function to generate mode for a given diagonal irrep matrix entry
def get_perturbation(images,DG,irreps,numUstar,a_map,basis,iv):

    mode = np.zeros((iv.numIm,iv.numAtoms,3))
    vec_list = np.zeros((1,iv.numIm,iv.numAtoms,3))

    for i in range(1,iv.numIm-1):
        atoms1 = images[i].get_scaled_positions()
        num1 = images[i].get_atomic_numbers()

        for j in range(0,iv.numAtoms):
            if basis[j] == 1:
                for m in range(0,3):
                    temp_mode = np.zeros((1,iv.numIm,iv.numAtoms,3))
                    for l in range(0,len(DG[0])):

                        #t_coord = np.dot(DG[0][l],atoms1[j])
                        #t_coord = (t_coord + DG[1][l])%1.0

                        if l < numUstar:
                           #atoms2 = images[i].get_scaled_positions()
                           #num2 = images[i].get_atomic_numbers()
                           t_im = i
                        else:
                           #atoms2 = images[numIm-1-i].get_scaled_positions()
                           #num2 = images[numIm-1-i].get_atomic_numbers()
                           t_im = iv.numIm-1-i

                        #for k in range(0,numAtoms):
                            #if (closewrapped (t_coord,atoms2[k],vectol) and \
                                #num1[j] == num2[k] and basis[j] == 1):

                        # was in commented for-loop above
                        temp_mode[0,t_im,np.nonzero(a_map[l,t_im,j,:])[0][0],:] += irreps[l,0,0]*DG[0][l][:,m]

                    if not np.any([np.sum(np.multiply(temp_mode[0,:,:,:], vec_list[t,:,:,:])) \
                    for t in range(0,len(vec_list[:,0,0,0]))])                       \
                    and np.any(temp_mode):

                        vec_list = np.concatenate((vec_list, temp_mode))


    if iv.r_co:
        for i in range(0,len(vec_list[:,0,0,0])):
            mode += random.uniform(-1.0, 1.0)*vec_list[i,:,:,:]
    else:
        mode = np.sum(vec_list,axis=0)

    return mode

# -- Function to generate mode for a given off-diagonal irrep matrix entry
def get_perturbation_odiag(DG,images,numUstar,irreps,mat_ent,vec,a_map,basis,iv):
    mode = np.zeros((iv.numIm,iv.numAtoms,3))
    vec_list = np.zeros((1,iv.numIm,iv.numAtoms,3))



    for i in range(1,iv.numIm-1):
        atoms1 = images[i].get_scaled_positions()
        num1 = images[i].get_atomic_numbers()

        for j in range(0,iv.numAtoms):
            if basis[j] == 1 and any([vec[i,j,p] != 0 for p in range(0,3)]):
                for m in range(0,3):
                    temp_mode = np.zeros((1,iv.numIm,iv.numAtoms,3))
                    for l in range(0,len(DG[0])):

                        #t_coord = np.dot(DG[0][l],atoms1[j])
                        #t_coord = (t_coord + DG[1][l])%1.0

                        if l < numUstar:
                           #atoms2 = images[i].get_scaled_positions()
                           #num2 = images[i].get_atomic_numbers()
                           t_im = i
                        else:
                           #atoms2 = images[numIm-1-i].get_scaled_positions()
                           #num2 = images[numIm-1-i].get_atomic_numbers()
                           t_im = iv.numIm-1-i

                        #for k in range(0,numAtoms):
                            #if (closewrapped (t_coord,atoms2[k],vectol) and \
                                #num1[j] == num2[k] and basis[j] == 1):


                        n_coord = np.dot(DG[0][l],vec[i,j,:])

                        temp_mode[0,t_im,np.nonzero(a_map[l,t_im,j,:])[0][0],:] += irreps[l,mat_ent,0]*n_coord



                    vec_list = np.concatenate((vec_list, temp_mode))


    mode = np.sum(vec_list,axis=0)


    return mode



###
### -- Main function to generate perturbation
###
def gen_perturb(path,irreps,iv):

    images = path.get_images()
    DG = path.get_DG()
    numUstar = path.numUstar

    outputfile = open(iv.image_dir + "/../results/output.out", "a")

    print ("\n\n\n\n"+("="*27+"\n")*2+"\nGenerating perturbations...\n\n"+("="*27+"\n")*2+"\n")
    outputfile.write("\n\n\n\n"+("="*27+"\n")*2+"\nGenerating perturbations...\n\n"+("="*27+"\n")*2+"\n")

    print ("***THE FOLLOWING DATA IS FOR THE PERTURBED IMAGES***\n\n")
    outputfile.write("***THE FOLLOWING DATA IS FOR THE PERTURBED IMAGES***\n\n")


    # -- Generate starting basis
    basis = np.zeros(iv.numAtoms)
    atoms1 = images[0].get_scaled_positions()
    atoms2 = images[iv.numIm-1].get_scaled_positions()
    m_vec = np.dot(np.linalg.inv(images[0].get_cell()),[iv.min_move,iv.min_move,iv.min_move])

    for i in range(0,iv.numAtoms):
        if not closewrapped (atoms1[i,:],atoms2[i,:],m_vec):
            basis[i] = 1

    print ("------- Atoms included in basis:")
    outputfile.write("------- Atoms included in basis:\n")

    symbols=images[0].get_chemical_symbols()
    for i in range(0,iv.numAtoms):
        if basis[i] == 0:
            print  (str(symbols[i])+" No")
            outputfile.write(str(symbols[i])+" No\n")
        else:
            print (str(symbols[i])+" Yes")
            outputfile.write(str(symbols[i])+" Yes\n")


    # -- Generate matrix showing atom mapping for each operations
    a_map = gen_atom_map(DG,images,numUstar,basis,iv)

    #for i in range(0,len(DG[0])):
     #  for j in range(1,iv.numIm-1):
      #     for k in range(0,iv.numAtoms):
       #        if sum(a_map[i][j][k][:]) == 0:
        #           print ("*************Atom matching error!!!!")


    # -- Generate and apply modes for an irrep of arbitrary dimension
    pt = np.zeros((iv.numIm,iv.numAtoms,3))


    pt_mode_init = get_perturbation(images,DG,irreps,numUstar,a_map,basis,iv)
    pt_mode_init = pt_mode_init/np.linalg.norm(np.ndarray.flatten(pt_mode_init))

    pt +=  iv.m_co[0]*pt_mode_init

    if iv.irr_dim > 1:
        for i in range(1,iv.irr_dim):
            pt_mode = get_perturbation_odiag(DG,images,numUstar,irreps,i,pt_mode_init,a_map,basis,iv)
            pt += iv.m_co[i]*(pt_mode/np.linalg.norm(np.ndarray.flatten(pt_mode)))


    p_vec = [iv.p_mag/np.linalg.norm(images[0].get_cell()[m,:]) for m in range(0,3)]
    #print pt[3]/np.amax(abs(pt[3]))
    #print(np.amax(p_vec*(pt[:,:,:]/np.amax(abs(pt)))))
    images_alt=[]

    for i in range(0,iv.numIm):
        images_alt.append(images[i].copy())
        images_alt[i].set_scaled_positions((images_alt[i].get_scaled_positions() + p_vec*(pt[i,:,:]/np.amax(abs(pt)))))


    return images_alt
