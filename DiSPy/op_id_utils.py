import numpy as np

def magic(mat, order, gentol):
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


def vectorID(vec, gentol):
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
        vec = np.dot(vec, 1 / np.sqrt(np.square(vec[0]) + np.square(vec[1]) + np.square(vec[2])))
    vectext = ""
    if isInt:
        for element in vec:
            if element < 0:
                vectext += str(abs(int(element))) + "-bar "
            else:
                vectext += str(int(element)) + " "
    vectext = vectext.strip()
    return isInt, vec, vectext


def operationAttributes(rot, tran, gentol):
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
        axis, _ = magic(np.dot(rot, det), order, gentol)
        axisInt, axis, axistext = vectorID(axis, gentol)
    translate = False
    if not np.allclose(tran, [0, 0, 0], atol=gentol):
        translate = True
        _, magicmatrix = magic(rot, order, gentol)
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
            result += "an inversion with a translation of " + str(tran) + " but no intrinsic translation component."
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


#-- Print operations in Bilbao format
def groupBB (path,iv):

    DG = path.get_DG()

    outputfile = open(iv.image_dir + "/../results/output.out", "a")
    num = 0
    group = ""
    s = ""
    rot = []
    tran = []
    for i in  range(0,len(DG[0])):
        for j in range(0,3):
            first = True
            if abs(DG[0][i][j][0]) > iv.gentol:
                s += printNumber(DG[0][i][j][0], True, 'x', first, iv)
                first = False
            if abs(DG[0][i][j][1]) > iv.gentol:
                s += printNumber(DG[0][i][j][1], True, 'y', first, iv)
                first = False
            if abs(DG[0][i][j][2]) > iv.gentol:
                s += printNumber(DG[0][i][j][2], True, 'z', first, iv)
                first = False
            if abs(DG[1][i][j]) > iv.gentol:
                s += printNumber(DG[1][i][j], False, 'a', first, iv)
            if j != 2:
                s += ","
        s += '%0A'
    print("\n-------\n------- Symmetry operations in Bilbao IDENTIFY GROUP format:\n-------\n")
    outputfile.write("\n-------\n------- Symmetry operations in Bilbao IDENTIFY GROUP format:\n-------\n\n")

    print(s.replace("%0A","\n"))
    outputfile.write(s.replace("%0A","\n"))

# -- Function for printing Bilbao input in a nice way
def printNumber (num,type,var,first,iv):
    if abs(num - 1) < iv.gentol:
        if type:
            if first:
                return var
            else:
                return "+" + var
        else:
            return "+1"
    elif abs(num + 1) < iv.gentol:
        if type:
            return "-" + var
        else:
            return " - 1"
    elif abs(num) < iv.gentol:
        if type:
            return ""
        else:
            return ""
    elif abs(num - 0.5) < iv.gentol:
        if type:
            if first:
                return "1/2" + var
            else:
                return "+1/2" + var
        else:
            return "+1/2"
    elif abs(num + 0.5) < iv.gentol:
        if type:
            return "-1/2" + var
        else:
            return "-1/2"
    elif abs(num - 1.0/4) < iv.gentol:
        if type:
            if first:
                return "1/4" + var
            else:
                return "+1/4" + var
        else:
            return "+1/4"
    elif abs(num + 1.0/4) < iv.gentol:
        if type:
            return "-1/4" + var
        else:
            return "-1/4"
    elif abs(num - 1.0/3) < iv.gentol:
        if type:
            if first:
                return "1/3" + var
            else:
                return "+1/3" + var
        else:
            return "+1/3"
    elif abs(num + 1.0/3) < iv.gentol:
        if type:
            return "-1/3" + var
        else:
            return "-1/3"
    elif abs(num - 2.0/3) < iv.gentol:
        if type:
            if first:
                return "2/3" + var
            else:
                return "+2/3" + var
        else:
            return "+2/3"
    elif abs(num + 2.0/3) < iv.gentol:
        if type:
            return "-2/3" + var
        else:
            return "-2/3"
    elif abs(num - 3.0/4) < iv.gentol:
        if type:
            if first:
                return "3/4" + var
            else:
                return "+3/4" + var
        else:
            return "+3/4"
    elif abs(num + 3.0/4) < iv.gentol:
        if type:
            return "-3/4" + var
        else:
            return "-3/4"
    elif abs(num - 1.0/6) < iv.gentol:
        if type:
            if first:
                return "1/6" + var
            else:
                return "+1/6" + var
        else:
            return "+1/6"
    elif abs(num + 1.0/6) < iv.gentol:
        if type:
            return "-1/6" + var
        else:
            return "-1/6"
    elif abs(num - 5.0/6) < iv.gentol:
        if type:
            if first:
                return "5/6" + var
            else:
                return "+5/6" + var
        else:
            return "+5/6"
    elif abs(num + 5.0/6) < iv.gentol:
        if type:
            return "-5/6" + var
        else:
            return "-5/6"
    else:
        if type:
            if num > 0 and not first:
                return "+" + str(num) + var
            else:
                return str(num) + var
        else:
            if num > 0:
                return "+" + str(num)
            else:
                return "-" + str(num)
