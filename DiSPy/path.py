import numpy as np


class make_path:
    def __init__(self,num_im):
        self.images = [None]*num_im
        self.DG = [None]
        self.DG_std = [None]
        self.numUstar = None
        self.img_sym_data = [None]*num_im


    def set_images(self,ase_images):
        if len(ase_images) != len(self.images):
            raise ValueError("Input list has wrong length or shape.")
        else:
            self.images = ase_images

    def get_images(self):
        return self.images[:]

    def set_DG(self,DG):
        if np.shape(DG[0][0]) != (3,3) or np.shape(DG[1][0]) != (3,):
            raise ValueError("Input array has wrong length or shape.")
        else:
            self.DG = DG

    def get_DG(self):
        return self.DG[:]

    def set_DG_std(self,DG_std):
        if np.shape(DG_std[0][0]) != (4,4) or np.shape(DG_std[1][0]) != (3,):
            raise ValueError("Input array has wrong length or shape.")
        else:
            self.DG_std = DG_std

    def get_DG_std(self):
        return self.DG_std[:]

    def set_img_sym_data(self,img_sym_data):
        if type(img_sym_data[0]) != dict or len(img_sym_data) != len(self.images):
            raise ValueError("Input list has wrong length or data types.")
        else:
            self.img_sym_data = img_sym_data

    def get_img_sym_data(self):
        return self.img_sym_data[:]
