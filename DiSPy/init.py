import argparse
from ase.io import read, write
import os
import numpy as np
import time

class in_var():

    def __init__(self, input_dir):

        # -- Initialize parameters
        self.perturb = False
        self.interpolate = False
        self.numIm = 0
        self.symprec = 1e-3
        self.angtol = -1.0
        self.image_dir = ""
        self.gentol = 0.001
        self.vectol = [0.001,0.001,0.001]
        self.vectol2 = [-1,-1,-1]
        self.angstrom1 = False
        self.angstrom2 = False
        self.trnum = 0
        self.tr_mat = np.identity(3)
        self.oshift = np.zeros(3)
        self.irr_num = 0
        self.irr_dim = 1
        self.m_co = [1]
        self.min_move = 0.0
        self.r_co = False
        self.p_mag = 0.1
        self.o_format = 'vasp'
        self.i_format = ''
        self.super_dim = [1.0,1.0,1.0]
        self.temptol = []
        self.temptol2 = []


        # -- Finds specified input file
        # Format *(in any order):
        # perturb = true/false                      *(to perturb or not)
        # interpolate = true/false                  *(to linearly interpolate or not)
        # images = ##                               *(if interpolation is on, specify the number of images (odd positive integer))
        # symprec = ##                              *(symmetry tolerance)
        # angle_tolerance = ##                      *(angular tolerance)
        # general_tolerance = ##                    *(general tolerance)
        # vector_tolerance (f/a) = ##,##,##         *(vector tolerance, fractional/angstrom)
        # translation_tolerance (f/a) = ##,##,##    *(translation vector tolerance, fractional/angstrom)
        # image_dir = './dir/'                      *(image directory to be read)
        # trans_num = ##                            *(image to take transformation matrix and origin shift from (0 - use user inputted))
        # trans_mat = ##,##,##,##,##,##,##,##,##    *(entries of transformation matrix)
        # oshift = ##,##,##                         *(entries of origin shift [frac])
        # min_move = ##                             *(minimum movement in any dimension to have atom included in basis [A])
        # input_format = ''                         *(format of input files - see ase listing)
        # output_format = ''                        *(format of outputted files - see ase listing)
        # irr_num = ##                              *(irrep number from listings in 'PIR_data.txt' by Stokes et al.)
        # irr_dim = #                               *(dimension of irrep)
        # mode_types = #,#,...                      *(types of modes and combination order)
        # mode_coeff = ##,##,...                    *(coefficients of modes)
        # random = true/false                       *(turn on or off random coefficients for modes)
        # perturb_mag = ##                          *(maximum magnitude of perturbation [A])
        # super_dim = ##,##,##                      *(vector showing number of unit cells in supercell [frac])

        ###
        ## FUNCTION TO READ INPUT FILE AND SET GLOBAL VARIABLES
        ##


        # -- Reads input file and adjusts parameters accordingly
        f = open(input_dir)
        line = f.readline()

        while line != "":
            temp = line.strip().upper().replace(" ","")
            if temp[:12] == "INTERPOLATE=":
                if temp[12:] == "TRUE" or temp[12:] == "T":
                    self.interpolate = True
            elif temp[:8] == "PERTURB=":
                if temp[8:] == "TRUE" or temp[8:] == "T":
                    self.perturb = True
            elif temp[:7] == "IMAGES=":
                try:
                    self.numIm = int(temp[7:])
                except:
                    self.numIm = 0
            elif temp[:8] == "SYMPREC=":
                try:
                    self.symprec = float(temp[8:])
                except:
                    self.symprec = 1e-3
            elif temp[:16] == "ANGLE_TOLERANCE=":
                try:
                    self.angtol = float(temp[16:])
                except:
                    self.angtol = -1.0
            elif temp[:18] == "GENERAL_TOLERANCE=":
                try:
                    self.gentol = float(temp[18:])
                except:
                    self.gentol = 0.001
            elif temp[:20] == "VECTOR_TOLERANCE(F)=":
                angstrom1 = False
                try:
                    self.vectol = []
                    temp = temp[20:] + ","
                    index = temp.index(",")
                    while index != -1:
                        self.vectol.append(float(temp[:index]))
                        temp = temp[index+1:]
                        try:
                            index = temp.index(",")
                        except ValueError:
                            index = -1
                except:
                    self.vectol = [0.001,0.001,0.001]
            elif temp[:20] == "VECTOR_TOLERANCE(A)=":
                try:
                    self.angstrom1 = True
                    self.temptol = []
                    temp = temp[20:] + ","
                    index = temp.index(",")
                    while index != -1:
                        self.temptol.append(float(temp[:index]))
                        temp = temp[index+1:]
                        try:
                            index = temp.index(",")
                        except ValueError:
                            index = -1
                except:
                    self.vectol = [0.001,0.001,0.001]
                    self.angstrom = False
            elif temp[:25] == "TRANSLATION_TOLERANCE(F)=":
                self.angstrom2 = False
                try:
                    self.vectol2 = []
                    temp = temp[25:] + ","
                    index = temp.index(",")
                    while index != -1:
                        self.vectol2.append(float(temp[:index]))
                        temp = temp[index+1:]
                        try:
                            index = temp.index(",")
                        except ValueError:
                            index = -1
                except:
                    self.vectol2 = [-1,-1,-1]
            elif temp[:25] == "TRANSLATION_TOLERANCE(A)=":
                try:
                    self.angstrom2 = True
                    self.temptol2 = []
                    temp = temp[25:] + ","
                    index = temp.index(",")
                    while index != -1:
                        self.temptol2.append(float(temp[:index]))
                        temp = temp[index+1:]
                        try:
                            index = temp.index(",")
                        except ValueError:
                            index = -1
                except:
                    self.vectol2 = [-1,-1,-1]
                    self.angstrom2 = False
            elif temp[:10] == "IMAGE_DIR=":
                self.image_dir = line.strip().replace(" ","")[10:]
            elif temp[:10] == "TRANS_NUM=":
                try:
                    self.trnum = int(temp[10:])
                except:
                    self.trnum = 0
            elif temp[:10] == "TRANS_MAT=":
                self.tr_mat = []
                temp = temp[10:] + ","
                index = temp.index(",")
                while index != -1:
                    self.tr_mat.append(float(temp[:index]))
                    temp = temp[index+1:]
                    try:
                        index = temp.index(",")
                    except ValueError:
                        index = -1
                self.tr_mat = np.reshape(self.tr_mat, (3,3))

            elif temp[:7] == "OSHIFT=":
                self.oshift = []
                temp = temp[7:] + ","
                index = temp.index(",")
                while index != -1:
                    self.oshift.append(float(temp[:index]))
                    temp = temp[index+1:]
                    try:
                        index = temp.index(",")
                    except ValueError:
                        index = -1
            elif temp[:10] == "SUPER_DIM=":
                try:
                    self.super_dim = []
                    temp = temp[10:] + ","
                    index = temp.index(",")
                    while index != -1:
                        self.super_dim.append(float(temp[:index]))
                        temp = temp[index+1:]
                        try:
                            index = temp.index(",")
                        except ValueError:
                            index = -1
                except:
                    self.super_dim = [1.0,1.0,1.0]
            elif temp[:8] == "IRR_NUM=":
                try:
                    self.irr_num = int(temp[8:])
                except:
                    self.irr_num = 0
            elif temp[:8] == "IRR_DIM=":
                try:
                    self.irr_dim = int(temp[8:])
                except:
                    self.irr_dim = 1
            elif temp[:11] == "MODE_COEFF=":
                try:
                    self.m_co = []
                    temp = temp[11:] + ","
                    index = temp.index(",")
                    while index != -1:
                        self.m_co.append(float(temp[:index]))
                        temp = temp[index+1:]
                        try:
                            index = temp.index(",")
                        except ValueError:
                            index = -1
                except:
                    self.m_co = [1]
            elif temp[:9] == "MIN_MOVE=":
                try:
                    self.min_move = float(temp[9:])
                except:
                    self.min_move = 0
            elif temp[:7] == "RANDOM=":
                if temp[7:] == "TRUE" or temp[7:] == "T":
                    self.r_co = True
            elif temp[:12] == "PERTURB_MAG=":
                try:
                    self.p_mag = float(temp[12:])
                except:
                    self.p_mag = 0.1
            elif temp[:14] == "OUTPUT_FORMAT=":
                self.o_format = str.lower(temp[14:])
            elif temp[:13] == "INPUT_FORMAT=":
                self.i_format = str.lower(temp[13:])

            line = f.readline()

        if np.array_equal(self.vectol2,[-1,-1,-1]):
            self.vectol2 = self.vectol
            self.temptol2 = self.temptol
            self.angstrom2 = self.angstrom1

        # -- Quick fix to remove .DS_Store file from MacOS
        try:
            os.remove(self.image_dir+'/'+'.DS_Store')
        except:
            pass

        # -- Create and manage the results directory
        if os.path.isdir(self.image_dir + "/../results"):
            i = 1
            while os.path.isdir(self.image_dir + "/../results_old_v" + str(i)):
                i = i + 1
            os.rename(self.image_dir + "/../results", self.image_dir + "/../results_old_v" + str(i))
        os.mkdir(self.image_dir + "/../results")

        # -- Create output file
        outputfile = open(self.image_dir + "/../results/output.out", "w")
        outputfile.write("""  
  ____  _ ____  ____        
 |  _ \(_) ___||  _ \ _   _ 
 | | | | \___ \| |_) | | | |
 | |_| | |___) |  __/| |_| |
 |____/|_|____/|_|    \__, |
                      |___/ 
""")
        outputfile.write("\nDiSPy v0.1.0 run started on " + time.strftime("%c") + "\n\n")
        outputfile.write("This code uses the tabulation of irreducible representations \nof the crystallographic space groups by H.T. Stokes and coworkers:\n\n")
        outputfile.write("\"H.T. Stokes et. al., Acta Cryst. A. 69, 388-395 (2013).\"\n\n")
        outputfile.write("===================================================================\n")
        outputfile.close()

        print("""  
  ____  _ ____  ____        
 |  _ \(_) ___||  _ \ _   _ 
 | | | | \___ \| |_) | | | |
 | |_| | |___) |  __/| |_| |
 |____/|_|____/|_|    \__, |
                      |___/ 
""")
        print("\nDiSPy v0.1.0 run started on " + time.strftime("%c") + "\n")
        print("This code uses the tabulation of irreducible representations \nof the crystallographic space groups by H.T. Stokes and coworkers:\n")
        print("\"H.T. Stokes et. al., Acta Cryst. A. 69, 388-395 (2013).\"\n")
        print("===================================================================")

    def get_images(self):

        # -- Finds image files from image directory
        # Files should be numbered in numerical/alphebetic order; no other files should be in the directory.
        self.image_dir,_,files = next(os.walk(self.image_dir))

        files = sorted(files)

        # -- Reads image files into images and does linear interpolation if specified
        images = []

        # -- Interpolation
        if self.interpolate:
            begin = read(self.image_dir+'/'+files[0],format=self.i_format)
            end = read(self.image_dir+'/'+files[len(files)-1],format=self.i_format)
            begin.wrap()
            end.wrap()
            b_pos = begin.get_scaled_positions()
            e_pos = end.get_scaled_positions()
            diff = e_pos - b_pos
            for i in diff:
                for j in range(0,3):
                    if i[j] > 0.5:
                        i[j] = -1.0 + i[j]
                    if i[j] < -0.5:
                        i[j] = 1.0 + i[j]
            for i in range (0,self.numIm):
                image = begin
                image.set_scaled_positions(b_pos + (i*(diff/(self.numIm-1))))
                image.wrap()
                images.append(image.copy())

        else:
            for file in sorted(files):
                if file != ".DS_Store" and file != ".localized":
                    image = read(self.image_dir+'/'+file,format=self.i_format)
                    image.wrap()
                    images.append(image.copy())


        self.numIm = len(images)
        self.numAtoms = len(images[0].get_scaled_positions())


        if self.angstrom1:
            self.vectol = [self.temptol[m]/np.linalg.norm(images[0].get_cell()[m,:]) for m in range(0,3)]
        if self.angstrom2:
            self.vectol2 = [self.temptol2[m]/np.linalg.norm(images[0].get_cell()[m,:]) for m in range(0,3)]


        return images
