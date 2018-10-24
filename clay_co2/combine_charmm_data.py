import os


class Charmm(object):
    def __init__(self, filename):
        os.system('ln ~/apps/bin/par_all36_cgenff.prm charmm.prm')
        self.filename = filename
        self.get_name_from_data()
        os.system('rm charmm.prm')

    @staticmethod
    def sort_ele(toptype, *args):
        ele_proiri = {'C': '0', 'N': '1', 'O': '2', 'P': '3', 'S': '4', 'H': '5'}
        if toptype == "bond":
            if args[0] == args[1]:
                return args[0] + '-' + args[1]
            first = ele_proiri[args[0][0]] + args[0][1:]
            last = ele_proiri[args[1][0]] + args[1][1:]
            for ix in range(min(len(first), len(last))):
                if first[ix] < last[ix]:
                    return args[0] + '-' + args[1]
                elif first[ix] > last[ix]:
                    return args[1] + '-' + args[0]
            if len(first) < len(last):
                return args[0] + '-' + args[1]
            else:
                return args[1] + '-' + args[0]

        if toptype == "angle":
            if args[0] == args[2]:
                return args[0] + '-' + args[1] + '-' + args[2]
            first = ele_proiri[args[0][0]] + args[0][1:]
            last = ele_proiri[args[2][0]] + args[2][1:]
            for ix in range(min(len(first), len(last))):
                if first[ix] < last[ix]:
                    return args[0] + '-' + args[1] + '-' + args[2]
                elif first[ix] > last[ix]:
                    return args[2] + '-' + args[1] + '-' + args[0]
            if len(first) < len(last):
                return args[0] + '-' + args[1] + '-' + args[2]
            else:
                return args[2] + '-' + args[1] + '-' + args[0]

        if toptype == "dihedral":
            if args[1] != args[2]:
                first = ele_proiri[args[1][0]] + args[1][1:]
                last = ele_proiri[args[2][0]] + args[2][1:]
                for ix in range(min(len(first), len(last))):
                    if first[ix] < last[ix]:
                        return args[0] + '-' + args[1] + '-' + args[2] + '-' + args[3]
                    elif first[ix] > last[ix]:
                        return args[3] + '-' + args[2] + '-' + args[1] + '-' + args[0]
                if len(first) < len(last):
                    return args[0] + '-' + args[1] + '-' + args[2] + '-' + args[3]
                else:
                    return args[3] + '-' + args[2] + '-' + args[1] + '-' + args[0]
            else:
                if args[0] == args[3]:
                    return args[0] + '-' + args[1] + '-' + args[2] + '-' + args[3]
                else:
                    first = ele_proiri[args[0][0]] + args[0][1:]
                    last = ele_proiri[args[3][0]] + args[3][1:]
                    for ix in range(min(len(first), len(last))):
                        if first[ix] < last[ix]:
                            return args[0] + '-' + args[1] + '-' + args[2] + '-' + args[3]
                        elif first[ix] > last[ix]:
                            return args[3] + '-' + args[2] + '-' + args[1] + '-' + args[0]
                    if len(first) < len(last):
                        return args[0] + '-' + args[1] + '-' + args[2] + '-' + args[3]
                    else:
                        return args[3] + '-' + args[2] + '-' + args[1] + '-' + args[0]

    def get_name_from_data(self, workdir='./'):
        atom_in_data = []
        atom_in_str = []
        atom_tmp = []
        bond_tmp = []
        ang_tmp = []
        dih_tmp = []
        imp_tmp = []
        atommass = []
        atomcharge = []
        atomtype = []
        bondtype = []
        angtype = []
        dihtype = []
        imptype = []
        atompara = []
        bondpara = []
        angpara = []
        dihpara = []
        file2 = open(workdir + self.filename + ".data")
        file3 = open(workdir + self.filename + ".txt")
        if os.path.exists("charmm.prm"):
            filepara = open("charmm.prm")
            paralines = filepara.readlines()
            filepara.close()
        else:
            paralines = []
            print("can't find parameter file for charmm field")

        # get atom name from the str files
        for perline in file3.readlines():
            tt = perline.split()
            if len(tt) and tt[0] == "ATOM":
                atom_in_str.append(tt[2])
                atomcharge.append(tt[3])
        file3.close()

        # find the atom section
        while 1:
            tt = file2.readline()
            if len(tt.split()) and tt.split()[0] == "Atoms":
                break
        file2.readline()
        # get atom coordination and atom name in the original data file
        while 1:
            tt = file2.readline().split()
            if len(tt):
                atom_tmp.append([tt[4], tt[5], tt[6]])
                atom_in_data.append(tt[11])
            else:
                break
        for ix in range(len(atom_in_str)):
            if atom_in_str[ix] not in atomtype:
                atomtype.append(atom_in_str[ix])
        for ix in range(len(atomtype)):
            tt_test = atomtype[ix].split('-')
            whe_find = 0
            for perline in range(len(paralines)):
                tt = paralines[len(paralines) - perline - 1].split()
                # if tt[:len(tt_test)] == [ix for ix in tt_test] and tt[4] == '!':
                if tt[:len(tt_test)] == [ix for ix in tt_test]:
                    atompara.append(['%.10f' % (-float(tt[2])), '%.10f' % (float(tt[3]) * 2 ** (5 / 6))])
                    whe_find = 1
                    break
            if not whe_find:
                print("can't find the parameter for atom:", atomtype[ix])
        for ix in range(len(atomtype)):
            tt_test = atomtype[ix].split('-')
            whe_find = 0
            for perline in range(len(paralines)):
                tt = paralines[perline].split()
                if len(tt) > 2 and tt[0] == "MASS" and tt[2] == tt_test[0]:
                    atommass.append('%.6f' % (float(tt[3])))
                    whe_find = 1
                    break
            if not whe_find:
                print("can't find the mass for atom:", atomtype[ix])

        # get bond information
        file2.readline()
        file2.readline()
        while 1:
            tt = file2.readline().split()
            if len(tt):
                bond_tmp.append([tt[2], tt[3], self.sort_ele("bond", atom_in_str[int(tt[2]) - 1],
                                                             atom_in_str[int(tt[3]) - 1])])
            else:
                break
        for ix in range(len(bond_tmp)):
            if bond_tmp[ix][2] not in bondtype:
                bondtype.append(bond_tmp[ix][2])
        for ix in range(len(bondtype)):
            tt_test = bondtype[ix].split('-')
            whe_find = 0
            for perline in range(len(paralines)):
                tt = paralines[perline].split()
                if tt[:len(tt_test)] == [ix for ix in tt_test] and tt[4] == '!':
                    bondpara.append([bondtype[ix], tt[2], tt[3]])
                    whe_find = 1
                    break
            if not whe_find:
                print("can't find the parameter for bond:", bondtype[ix])

        # get angle information
        file2.readline()
        file2.readline()
        while 1:
            tt = file2.readline().split()
            if len(tt):
                ang_tmp.append([tt[2], tt[3], tt[4], self.sort_ele(
                    "angle", atom_in_str[int(tt[2]) - 1], atom_in_str[int(tt[3]) - 1], atom_in_str[int(tt[4]) - 1])])
            else:
                break
        for ix in range(len(ang_tmp)):
            if ang_tmp[ix][3] not in angtype:
                angtype.append(ang_tmp[ix][3])
        for ix in range(len(angtype)):
            whe_find = 0
            tt_test = angtype[ix].split('-')
            for perline in range(len(paralines)):
                tt = paralines[perline].split()
                if tt[:len(tt_test)] == [ix for ix in tt_test]:
                    whe_find = 1
                    if tt[5] == '!':
                        angpara.append(['harmonic', tt[3], tt[4]])
                    else:
                        angpara.append(['charmm', tt[3], tt[4], tt[5], tt[6]])
                    break
            if not whe_find:
                print("can't find the parameter for ang:", angtype[ix])

        # get dihedral information
        file2.readline()
        file2.readline()
        while 1:
            tt = file2.readline().split()
            if len(tt):
                dih_tmp.append([tt[2], tt[3], tt[4], tt[5],
                                self.sort_ele("dihedral", atom_in_str[int(tt[2]) - 1], atom_in_str[int(tt[3]) - 1],
                                              atom_in_str[int(tt[4]) - 1], atom_in_str[int(tt[5]) - 1])])
            else:
                break
        for ix in range(len(dih_tmp)):
            if dih_tmp[ix][4] not in dihtype:
                dihtype.append(dih_tmp[ix][4])
        for ix in range(len(dihtype)):
            whe_find = 0
            tt_test = dihtype[ix].split('-')
            for perline in range(len(paralines)):
                tt = paralines[perline].split()
                if tt[:len(tt_test)] == [ix for ix in tt_test]:
                    whe_find = 1
                    if tt[6] == '!':
                        dihpara.append(['harmonic', tt[4], tt[5]])
                    else:
                        dihpara.append(['charmm', tt[4], tt[5], tt[6], '0.00'])
                    break
            if not whe_find:
                print("can't find the parameter for dih:", dihtype[ix])

        file2.close()
        fileout = open(workdir + self.filename + 'tran.data', 'w')
        print("the monte", file=fileout)
        print(str(len(atom_tmp)) + " atoms", file=fileout)
        print(str(len(bond_tmp)) + " bonds", file=fileout)
        print(str(len(ang_tmp)) + " angles", file=fileout)
        print(str(len(dih_tmp)) + " dihedrals", file=fileout)
        print(str(len(imp_tmp)) + " impropers\n", file=fileout)
        print(str(len(atomtype)) + " atom types", file=fileout)
        print(str(len(bondtype)) + " bond types", file=fileout)
        print(str(len(angtype)) + " angle types", file=fileout)
        print(str(len(dihtype)) + " dihedral types", file=fileout)
        print(str(len(imptype)) + " improper types\n", file=fileout)
        print("\nMasses\n", file=fileout)
        for i in range(len(atomtype)):
            print("%d %s %s" % (i + 1, atommass[i], "# " + atomtype[i]), file=fileout)
        print("\nPair Coeffs\n", file=fileout)
        for i in range(len(atomtype)):
            print("%d %s %s" % (i + 1, ' '.join(atompara[i]), "# " + atomtype[i]), file=fileout)
        print("\nBond Coeffs\n", file=fileout)
        for i in range(len(bondtype)):
            print("%d %s %s %s" % (i + 1, bondpara[i][1], bondpara[i][2], "# " + bondtype[i]), file=fileout)
        print("\nAngle Coeffs\n", file=fileout)
        for i in range(len(angtype)):
            print("%d %s %s" % (i + 1, ' '.join(angpara[i]), "# " + angtype[i]), file=fileout)
        if len(dihtype):
            print("\nDihedral Coeffs\n", file=fileout)
            for i in range(len(dihtype)):
                print("%d %s %s" % (i + 1, ' '.join(dihpara[i]), "# " + dihtype[i]), file=fileout)
        if len(imptype):
            print("\nImproper Coeffs\n", file=fileout)
            for i in range(len(imptype)):
                print("%d %s" % (i + 1, "# " + imptype[i]), file=fileout)
        print("\nAtoms\n", file=fileout)
        for i in range(len(atom_tmp)):
            print("%8d %5d %4d %.6f %10.6f %10.6f %10.6f %-10s" %
                  (i + 1, 1, atomtype.index(atom_in_str[i]) + 1, float(atomcharge[i]), float(atom_tmp[i][0]),
                   float(atom_tmp[i][1]), float(atom_tmp[i][2]), "# " + atom_in_str[i]), file=fileout)
        print("\nBonds\n", file=fileout)
        for i in range(len(bond_tmp)):
            print("%4d %4d %8s %8s" % (i + 1, bondtype.index(bond_tmp[i][2]) + 1,
                                       bond_tmp[i][0], bond_tmp[i][1]), file=fileout)
        print("\nAngles\n", file=fileout)
        for i in range(len(ang_tmp)):
            print("%4d %4d %8s %8s %8s" % (i + 1, angtype.index(ang_tmp[i][3]) + 1,
                                           ang_tmp[i][0], ang_tmp[i][1], ang_tmp[i][2]), file=fileout)
        if len(dihtype):
            print("\nDihedrals\n", file=fileout)
            for i in range(len(dih_tmp)):
                print("%4d %4d %8s %8s %8s %8s" % (i + 1, dihtype.index(dih_tmp[i][4]) + 1, dih_tmp[i][0],
                                                   dih_tmp[i][1], dih_tmp[i][2], dih_tmp[i][3]), file=fileout)
        if len(imptype):
            print("\nImpropers\n", file=fileout)
            for i in range(len(imp_tmp)):
                print("%4d %4d %8s %8s %8s %8s" % (i + 1, imptype.index(imp_tmp[i][4]) + 1, imp_tmp[i][0],
                                                   imp_tmp[i][1], imp_tmp[i][2], imp_tmp[i][3]), file=fileout)
        fileout.close()


test1 = Charmm('ss')
