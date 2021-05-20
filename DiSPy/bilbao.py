import numpy as np

# -- Print operations in Bilbao format


def groupBB(path, io):

    DG = path.distortion_group.matrices

    s = ""
    for i in range(len(DG)):
        for j in range(0, 3):
            first = True
            if abs(DG[i].rotation_matrix[j][0]) > io.gentol:
                s += _printNumber(DG[i].rotation_matrix[j][0], True, "x", first, io.gentol)
                first = False
            if abs(DG[i].rotation_matrix[j][1]) > io.gentol:
                s += _printNumber(DG[i].rotation_matrix[j][1], True, "y", first, io.gentol)
                first = False
            if abs(DG[i].rotation_matrix[j][2]) > io.gentol:
                s += _printNumber(DG[i].rotation_matrix[j][2], True, "z", first, io.gentol)
                first = False
            if abs(DG[i].translation_vector[j]) > io.gentol:
                s += _printNumber(DG[i].translation_vector[j], False, "a", first, io.gentol)
            if j != 2:
                s += ","
        s += "%0A"
    io.print("\n-------\n------- Symmetry operations in Bilbao IDENTIFY GROUP format:\n-------\n")
    io.print(s.replace("%0A", "\n"))


# -- Function for printing Bilbao input in a nice way
def _printNumber(num, type, var, first, gentol):
    if abs(num - 1) < gentol:
        if type:
            if first:
                return var
            else:
                return "+" + var
        else:
            return "+1"
    elif abs(num + 1) < gentol:
        if type:
            return "-" + var
        else:
            return " - 1"
    elif abs(num) < gentol:
        if type:
            return ""
        else:
            return ""
    elif abs(num - 0.5) < gentol:
        if type:
            if first:
                return "1/2" + var
            else:
                return "+1/2" + var
        else:
            return "+1/2"
    elif abs(num + 0.5) < gentol:
        if type:
            return "-1/2" + var
        else:
            return "-1/2"
    elif abs(num - 1.0 / 4) < gentol:
        if type:
            if first:
                return "1/4" + var
            else:
                return "+1/4" + var
        else:
            return "+1/4"
    elif abs(num + 1.0 / 4) < gentol:
        if type:
            return "-1/4" + var
        else:
            return "-1/4"
    elif abs(num - 1.0 / 3) < gentol:
        if type:
            if first:
                return "1/3" + var
            else:
                return "+1/3" + var
        else:
            return "+1/3"
    elif abs(num + 1.0 / 3) < gentol:
        if type:
            return "-1/3" + var
        else:
            return "-1/3"
    elif abs(num - 2.0 / 3) < gentol:
        if type:
            if first:
                return "2/3" + var
            else:
                return "+2/3" + var
        else:
            return "+2/3"
    elif abs(num + 2.0 / 3) < gentol:
        if type:
            return "-2/3" + var
        else:
            return "-2/3"
    elif abs(num - 3.0 / 4) < gentol:
        if type:
            if first:
                return "3/4" + var
            else:
                return "+3/4" + var
        else:
            return "+3/4"
    elif abs(num + 3.0 / 4) < gentol:
        if type:
            return "-3/4" + var
        else:
            return "-3/4"
    elif abs(num - 1.0 / 6) < gentol:
        if type:
            if first:
                return "1/6" + var
            else:
                return "+1/6" + var
        else:
            return "+1/6"
    elif abs(num + 1.0 / 6) < gentol:
        if type:
            return "-1/6" + var
        else:
            return "-1/6"
    elif abs(num - 5.0 / 6) < gentol:
        if type:
            if first:
                return "5/6" + var
            else:
                return "+5/6" + var
        else:
            return "+5/6"
    elif abs(num + 5.0 / 6) < gentol:
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
