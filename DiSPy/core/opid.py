import numpy as np


class OperationIdentification:

    def __init__(self, symmop, gentol):
        rot = symmop.rotation_matrix
        tran = symmop.translation_vector

        self._info = self.operationAttributes(rot, tran, gentol)

    @property
    def info(self):
        return self._info

    @staticmethod
    def _magic(mat, order, gentol):
        axismat = np.zeros((3, 3))
        for i in range(order):
            axismat += np.linalg.matrix_power(mat, i)
        vec = [axismat[0][0], axismat[1][0], axismat[2][0]]
        if np.allclose(vec, [0, 0, 0], atol=gentol):
            vec = [axismat[0][1], axismat[1][1], axismat[2][1]]
            if np.allclose(vec, [0, 0, 0], atol=gentol):
                vec = [axismat[0][2], axismat[1][2], axismat[2][2]]
        if abs(vec[0]) > gentol:
            vec = np.dot(vec, 1 / abs(vec[0]))
        else:
            if abs(vec[1]) > gentol:
                vec = np.dot(vec, 1 / abs(vec[1]))
            else:
                if abs(vec[2]) > gentol:
                    vec = np.dot(vec, 1 / abs(vec[2]))
        return vec, axismat

    @staticmethod
    def _vectorID(vec, gentol):
        isInt = False
        if abs(vec[0]) > gentol:
            vec = np.dot(vec, 1 / abs(vec[0]))
        else:
            if abs(vec[1]) > gentol:
                vec = np.dot(vec, 1 / abs(vec[1]))
            else:
                if abs(vec[2]) > gentol:
                    vec = np.dot(vec, 1 / abs(vec[2]))
        if np.allclose(vec, np.around(vec, decimals=0), atol=gentol):
            vec = np.around(vec, decimals=0)
            isInt = True
        elif np.allclose(np.dot(vec, 2), np.around(np.dot(vec, 2), decimals=0), atol=gentol):
            vec = np.around(np.dot(vec, 2))
            isInt = True
        elif np.allclose(np.dot(vec, 3), np.around(np.dot(vec, 3), decimals=0), atol=gentol):
            vec = np.around(np.dot(vec, 3))
            isInt = True
        elif np.allclose(np.dot(vec, 4), np.around(np.dot(vec, 4), decimals=0), atol=gentol):
            vec = np.around(np.dot(vec, 4))
            isInt = True
        else:
            vec = np.dot(
                vec, 1 / np.sqrt(np.square(vec[0]) + np.square(vec[1]) + np.square(vec[2])))
        vectext = ""
        if isInt:
            for element in vec:
                if element < 0:
                    vectext += str(abs(int(element))) + "-bar "
                else:
                    vectext += str(int(element)) + " "
        vectext = vectext.strip()
        return isInt, vec, vectext

    # -- Function to identify a symmetry operation

    def operationAttributes(self, rot, tran, gentol):
        det = np.linalg.det(rot)
        if abs(det - 1) < gentol:
            det = 1
        elif abs(det + 1) < gentol:
            det = -1
        trace = np.trace(rot)
        type = 3
        if abs(abs(trace) - 3) < gentol:
            type = 1
        elif abs(abs(trace) - 2) < gentol:
            type = 6
        elif abs(trace) < gentol:
            type = 3
        else:
            type += int(round(trace * det))
        order = type
        if type % 2 == 1 and abs(det + 1) < gentol:
            order *= 2
        axis = [0, 0, 0]
        axisInt = False
        if type != 1:
            axis, _ = self._magic(np.dot(rot, det), order, gentol)
            axisInt, axis, axistext = self._vectorID(axis, gentol)
        translate = False
        if not np.allclose(tran, [0, 0, 0], atol=gentol):
            translate = True
            _, magicmatrix = self._magic(rot, order, gentol)
            moved = np.dot(magicmatrix, tran)
        result = "The operation has been identified as "
        if type == 1 and not translate:
            if det == 1:
                result += "the identity."
            else:
                result += "an inversion."
        elif type == 1 and translate:
            if det == 1:
                result += "a translation of " + str(moved) + "."
            else:
                result += "an inversion with a translation of " + \
                    str(tran) + " but no intrinsic translation component."
        else:
            if type == 2 and det == -1:
                result += "a mirror across the "
            else:
                if type == 2:
                    result += "a twofold "
                elif type == 3:
                    result += "a threefold "
                elif type == 4:
                    result += "a fourfold "
                elif type == 6:
                    result += "a sixfold "
                if det == 1:
                    result += "rotation along the "
                if det == -1:
                    result += "rotoinversion along the "
            if axisInt:
                if type == 2 and det == -1:
                    result += axistext + " plane"
                else:
                    result += axistext + " axis"
            else:
                if type == 2 and det == -1:
                    result += str(axis) + " plane"
                else:
                    result += str(axis) + " axis"
            if translate:
                result += " with a translation of " + str(tran)
                if np.allclose(moved, [0, 0, 0], atol=gentol):
                    if type == 2 and det == -1:
                        result += " but no intrinsic glide component."
                    else:
                        result += " but no intrinsic screw component."
                else:
                    if type == 2 and det == -1:
                        result += " and an intrinsic glide component of "
                    else:
                        result += " and an intrinsic screw component of "
                    result += str(moved) + "."
            else:
                result += "."
        return result
