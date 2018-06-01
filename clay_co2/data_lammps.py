import math
import numpy as np
import os
import random
alsubt = 25
sisubt = 12
addzlen = 15.5
add_z1 = [["na", 23, 10, 37], ["h2o", 18, 1.5, 550]]
add_z2 = [["na", 23, 10, 37], ["h2o", 18, 1.5, 550]]
add_z3 = [["na", 23, 10, 37], ["co2", 44, 1.5, 100], ["h2o", 18, 1.5, 1000]]
#add_z4 = [["na", 23, 10, 10], ["h2o", 18, 1.5, 10]]
add_y = [["co2", 44, 1, 500], ["h2o", 18, 1.3, 5000]]  # name mol/g density number
atoms = []    # the list to store all the atoms in the form: ele_name x y z
box = []      # the box size read from the initial clay
ang = []     # the list involves  all the angle
bond = []     # the bond involves all the bond
dih = []
imp = []
atom_turn = []  # the list to store all the atom types
bond_turn = []  # the list to store all the bond types
ang_turn = []
dih_turn = []
imp_turn = []
list_dic = {"Si": 1, "Al": 2, "O": 3, "O1": 4, "O2": 5, "O3": 6,
            "O-H": 7, "H": 8, "eow": 9, "ehw": 10}  # the order to sort the atoms in initial clay
char_dic = {"st": 2.1000, "ao": 1.5750, "ob": -1.0500, "obts": -1.1688, "obos": -1.1808,
            "ohy": -0.9500, "ohs": -1.0808, "mgo": 1.3600, "at": 1.5750,
            "hoy": 0.4250, "eow": -0.8200, "ehw": 0.4100, "na": 1.0000,
            "ow": -0.8200, "hw": 0.4100, "emgo": 1.2292, "cl": -1.0000,
            "oc": -0.3256, "co": 0.6512}   # charge dictionary
mass = {"st": 28.08600, "at": 26.98150, "ao": 26.98150,  "mgo": 24.30500, "emgo": 24.30500, "ob": 15.99940,
        "obts": 15.99940, "obos": 15.99940, "ohy": 15.99940, "ohs": 15.99940, "hoy": 1.007940,
        "eow": 15.99940, "ehw": 1.007940, "ow": 15.99940, "hw": 1.007940, "na": 22.989768,
        "cl": 35.45000, "oc": 15.99940, "co": 12.00000}  # mass dictionary
pair_coeff = {"st": [0.0000018405, 3.30203000], "at": [0.0000018405, 3.30203000],
              "ao": [0.0000013298, 4.2712400], "mgo": [0.0000009030, 5.2643200],
              "ob": [0.1554000000, 3.1655400], "obos": [0.1554000000, 3.1655400],
              "obts": [0.1554000000, 3.1655400], "ohy": [0.1554000000, 3.1655400],
              "ohs": [0.1554000000, 3.1655400], "hoy": [0.0000000000, 0.0000000],
              "eow": [0.1553000000, 3.1660000], "ehw": [0.0000000000, 0.0000000],
              "ow": [0.1553000000, 3.1660000], "hw": [0.0000000000, 0.0000000],
              "na": [0.1301000000, 2.3500100], "emgo": [0.0000009030, 5.2643200],
              "cl": [0.1001000000, 4.3999620], "oc": [0.1645070000, 3.0280000],
              "co": [0.0559270000, 2.8000000]}   # pair_coeff dictionary
bond_coeff = {"hoy-ohy": ["morse", 132.2941, 2.1350, 0.9572], "hoy-ohs": ["morse", 132.2941, 2.1350, 0.9572],
              "eow-ehw": ["harmonic", 554.1300, 1.000], "hw-ow": ["harmonic", 450.2698, 1.000],
              "oc-co": ["harmonic", 2017.9254, 1.162]}
ang_coeff = {"hoy-ohy-ao": [15.00000, 110.10000], "hoy-ohs-mgo": [6.00000, 110.10000],
             "hoy-ohs-ao": [15.00000, 110.10000], "hoy-ohy-st": [3.00000, 110.10000],
             "hoy-ohs-at": [4.00000, 110.10000], "ehw-eow-ehw": [45.77000, 109.47000],
             "hw-ow-hw": [45.77000, 109.47000], "oc-co-oc": [108.00669, 180.00000]}
dih_coeff = {}
imp_coeff = {}
mol_bond = []
mol_ang = []
mol_dih = []
mol_imp = []
molnum = []
moltimes = []
molname = []
moltot = []
whemol = []


def sort_atom(target):
    global list_dic
    onewater = []
    for ix in range(len(target)):
        target[ix].append(list_dic[target[ix][0]])
    atoms_2 = sorted(target, key=lambda x: x[4])
    for ix in range(len(atoms_2)):
        atoms_2[ix].remove(atoms_2[ix][4])
    for ix in range(len(atoms_2)):
        if atoms_2[ix][0] == "eow":
            water_tmp = []
            tmp_1 = [atoms_2[ix][1], atoms_2[ix][2], atoms_2[ix][3]]
            water_tmp.append(atoms_2[ix])
            for jx in range(len(atoms_2)):
                if atoms_2[jx][0] == "ehw":
                    tmp_2 = [atoms_2[jx][1], atoms_2[jx][2], atoms_2[jx][3]]
                    if caldis(tmp_1, tmp_2, 1.5) == 0:
                        water_tmp.append(atoms_2[jx])
            onewater.append(water_tmp)
    xx = []
    for ix in range(len(atoms_2)):
        if atoms_2[len(atoms_2)-ix-1][0] == "ehw" or atoms_2[len(atoms_2)-ix-1][0] == "eow":
            xx.append(len(atoms_2)-ix-1)
    for ix in range(len(xx)):
        del atoms_2[xx[ix]]
    for ix in range(len(onewater)):
        atoms_2.append(onewater[ix][0])
        atoms_2.append(onewater[ix][1])
        atoms_2.append(onewater[ix][2])
    return atoms_2


def readmon():
    global atoms
    global box
    global molnum
    global moltimes
    global whemol
    global molname
    global mol_bond
    global mol_ang
    global mol_dih
    global mol_imp
    atom_mon = []
    mol_bond.append([])
    mol_ang.append([])
    mol_dih.append([])
    mol_imp.append([])
    mol_bond.append([])
    mol_ang.append([])
    mol_dih.append([])
    mol_imp.append([])
    xyzpath = "./edge_mon"  # the file to read clay
    file1 = open(xyzpath)
    box_tmp = file1.readline()
    box = [float(box_tmp.split()[0]), float(box_tmp.split()[1]), float(box_tmp.split()[2])]
    content = file1.readlines()
    file1.close()
    for ix in content:
        tmp_1 = ix.split()[0]
        tmp_2 = float(ix.split()[1])
        tmp_3 = float(ix.split()[2])
        tmp_4 = float(ix.split()[3])
        atom_mon.append([tmp_1, tmp_2, tmp_3, tmp_4])
    numwat = 0
    for ix in atom_mon:
        if ix[0] == "eow":
            numwat += 1
    molnum.append(len(content) - 3 * numwat)
    moltimes.append(1)
    whemol.append(0)
    molname.append("")
    molnum.append(3)
    moltimes.append(numwat)
    whemol.append(0)
    molname.append("")
    atoms.extend(sort_atom(atom_mon))


def caldis(atom1, atom2, dis):
    global box
    length = dis
    halfx = box[0]/2
    halfy = box[1]/2
    halfz = box[2]/2
    lenx = halfx - abs(halfx - abs(atom1[0]-atom2[0]))
    if lenx > length:
        return 1
    leny = halfy - abs(halfy - abs(atom1[1]-atom2[1]))
    if leny > length:
        return 1
    lenz = halfz - abs(halfz - abs(atom1[2]-atom2[2]))
    if lenz > length:
        return 1
    lenxyz = math.sqrt(lenx*lenx + leny*leny + lenz*lenz)
    if lenxyz < length:
        return 0


def readmon2(addz):
    global atoms
    global box
    global molnum
    global moltimes
    global whemol
    global mol_bond
    global mol_ang
    global mol_dih
    global mol_imp
    atom_mon = []
    mol_bond.append([])
    mol_ang.append([])
    mol_dih.append([])
    mol_imp.append([])
    mol_bond.append([])
    mol_ang.append([])
    mol_dih.append([])
    mol_imp.append([])
    xyzpath = "./edge_mon"  # the file to read clay
    arr_zz = np.array(np.array(atoms)[:, 3], dtype=float)
    zz_sort = np.sort(arr_zz)
    tmp_z1 = zz_sort[len(zz_sort)-1] - zz_sort[0] + 1.5
    tmp_z2 = addz
    tmp_z = max(tmp_z1, tmp_z2)
    print("addz =","%f" %(tmp_z))
    file1 = open(xyzpath)
    file1.readline()
    box[2] += tmp_z + 10
    content = file1.readlines()
    file1.close()
    for ix in content:
        tmp_1 = ix.split()[0]
        tmp_2 = float(ix.split()[1])
        tmp_3 = float(ix.split()[2])
        tmp_4 = float(ix.split()[3]) + tmp_z
        atom_mon.append([tmp_1, tmp_2, tmp_3, tmp_4])
    numwat = 0
    for ix in atom_mon:
        if ix[0] == "eow":
            numwat += 1
    molnum.append(len(content) - 3 * numwat)
    moltimes.append(1)
    whemol.append(0)
    molname.append(0)
    molnum.append(3)
    moltimes.append(numwat)
    whemol.append(0)
    molname.append(0)
    atoms.extend(sort_atom(atom_mon))


def findang():
    global atoms
    global ang
    global bond
    for ix in range(len(atoms)):
        if atoms[ix][0] in ["hoy"]:
            for jx in range(len(atoms)):
                if atoms[jx][0] in ["ohy", "ohs"]:
                    tmp_1 = [atoms[ix][1], atoms[ix][2], atoms[ix][3]]
                    tmp_2 = [atoms[jx][1], atoms[jx][2], atoms[jx][3]]
                    if caldis(tmp_1, tmp_2, 1.5) == 0:
                        bond.append([atoms[ix][0] + '-' + atoms[jx][0], ix + 1, jx + 1])
                        for kx in range(len(atoms)):
                            if atoms[kx][0] in ["mgo", "ao", "st", "at"]:
                                tmp_1 = [atoms[jx][1], atoms[jx][2], atoms[jx][3]]
                                tmp_2 = [atoms[kx][1], atoms[kx][2], atoms[kx][3]]
                                if caldis(tmp_1, tmp_2, 2.0) == 0:
                                    ang.append([atoms[ix][0]+'-'+atoms[jx][0]+'-'+atoms[kx][0],
                                                ix+1, jx+1, kx+1])


def tranohho():
    global atoms
    global ang
    global bond
    for ix in range(len(atoms)):
        if atoms[ix][0] == "H":
            atoms[ix][0] = "ho"
            for jx in range(len(atoms)):
                if atoms[jx][0] == "O" or atoms[jx][0] == "O-H":
                    tmp_1 = [atoms[ix][1], atoms[ix][2], atoms[ix][3]]
                    tmp_2 = [atoms[jx][1], atoms[jx][2], atoms[jx][3]]
                    if caldis(tmp_1, tmp_2, 1.5) == 0:
                        atoms[jx][0] = "oh"
                        # bond.append([atoms[ix][0] + '-' + atoms[jx][0], ix + 1, jx + 1])
                        # for kx in range(len(atoms)):
                        #     if atoms[kx][0] == "Al" or atoms[kx][0] == "Si":
                        #         tmp_1 = [atoms[jx][1], atoms[jx][2], atoms[jx][3]]
                        #         tmp_2 = [atoms[kx][1], atoms[kx][2], atoms[kx][3]]
                        #         if caldis(tmp_1, tmp_2, 2.0) == 0:
                        #             ang.append([atoms[ix][0]+'-'+atoms[jx][0]+'-'+atoms[kx][0],
                        #                         ix+1, jx+1, kx+1])


def findwatmon():
    global atoms
    global ang
    global bond
    for ix in range(len(atoms)):
        if atoms[ix][0] == "eow":
            tmp_1 = [atoms[ix][1], atoms[ix][2], atoms[ix][3]]
            tmp_3 = []
            for jx in range(len(atoms)):
                if atoms[jx][0] == "ehw":
                    tmp_2 = [atoms[jx][1], atoms[jx][2], atoms[jx][3]]
                    if caldis(tmp_1, tmp_2, 1.5) == 0:
                        tmp_3.append(jx)
                        bond.append([atoms[ix][0] + '-' + atoms[jx][0], ix + 1, jx + 1])
            ang.append([atoms[tmp_3[1]][0]+'-'+atoms[ix][0]+'-'+atoms[tmp_3[0]][0],
                        tmp_3[1]+1, ix+1, tmp_3[0]+1])


def subti(alfin, sifin):
    global atoms
    alindex = []
    siindex = []
    alsub = []
    sisub = []
    allal = []
    allsi = []
    siedge = []
    aledge = []
    for ix in range(len(atoms)):
        if atoms[ix][0] == "Al":
            allal.append(atoms[ix][2])
    for ix in range(len(atoms)):
        if atoms[ix][0] == "Si":
            allsi.append(atoms[ix][2])
    almin = sorted(allal)[0]
    almax = sorted(allal)[len(allal)-1]
    simin = sorted(allsi)[0]
    simax = sorted(allsi)[len(allsi)-1]
    for ix in range(len(atoms)):
        if atoms[ix][0] == "Al" and (atoms[ix][2] > (almax - 1) or atoms[ix][2] < (almin + 1)):
            aledge.append(ix)
        elif atoms[ix][0] == "Si" and (atoms[ix][2] > (simax - 1) or atoms[ix][2] < (simin + 1)):
            siedge.append(ix)
    for ix in range(len(atoms)):
        if atoms[ix][0] == "Al":
            atoms[ix][0] = "ao"
            alindex.append(ix)
        elif atoms[ix][0] == "Si":
            atoms[ix][0] = "st"
            siindex.append(ix)
        elif atoms[ix][0] == "O" or atoms[ix][0] == "O1" or atoms[ix][0] == "O2" or atoms[ix][0] == "O3":
            atoms[ix][0] = "ob"
    while len(alsub) < alfin:
        tmp = random.choice(alindex)
        judge = 1
        if tmp in alsub:
            continue
        if tmp in aledge:
            continue
        else:
            for jx in range(len(alsub)):
                atom_1 = [atoms[alsub[jx]][1], atoms[alsub[jx]][2], atoms[alsub[jx]][3]]
                atom_2 = [atoms[tmp][1], atoms[tmp][2], atoms[tmp][3]]
                if caldis(atom_1, atom_2, 3.5) == 0:
                    judge = 0
                    break
        if judge == 1:
            alsub.append(tmp)
            atoms[tmp][0] = "mgo"
    while len(sisub) < sifin:
        tmp = random.choice(siindex)
        judge = 1
        if tmp in sisub:
            continue
        if tmp in siedge:
            continue
        else:
            for jx in range(len(alsub)):
                atom_1 = [atoms[alsub[jx]][1], atoms[alsub[jx]][2], atoms[alsub[jx]][3]]
                atom_2 = [atoms[tmp][1], atoms[tmp][2], atoms[tmp][3]]
                if caldis(atom_1, atom_2, 3.5) == 0:
                    judge = 0
                    break
            for jx in range(len(sisub)):
                atom_1 = [atoms[sisub[jx]][1], atoms[sisub[jx]][2], atoms[sisub[jx]][3]]
                atom_2 = [atoms[tmp][1], atoms[tmp][2], atoms[tmp][3]]
                if caldis(atom_1, atom_2, 3.5) == 0:
                    judge = 0
                    break
        if judge == 1:
            sisub.append(tmp)
            atoms[tmp][0] = "at"
    for ix in alsub:
        for jx in range(len(atoms)):
            atom_1 = [atoms[ix][1], atoms[ix][2], atoms[ix][3]]
            atom_2 = [atoms[jx][1], atoms[jx][2], atoms[jx][3]]
            if caldis(atom_1, atom_2, 2.0) == 0 and atom_1 != atom_2:
                if atoms[jx][0] == "ob":
                    atoms[jx][0] = "obos"
                elif atoms[jx][0] == "oh":
                    atoms[jx][0] = "ohs"
    for ix in sisub:
        for jx in range(len(atoms)):
            atom_1 = [atoms[ix][1], atoms[ix][2], atoms[ix][3]]
            atom_2 = [atoms[jx][1], atoms[jx][2], atoms[jx][3]]
            if caldis(atom_1, atom_2, 2.0) == 0 and atom_1 != atom_2:
                if atoms[jx][0] == "ob":
                    atoms[jx][0] = "obts"
                elif atoms[jx][0] == "oh" or atoms[jx][0] == "toh":
                    atoms[jx][0] = "ohs"
    for ix in range(len(atoms)):
        if atoms[ix][0] == "oh":
            atoms[ix][0] = "ohy"
        elif atoms[ix][0] == "ho":
            atoms[ix][0] = "hoy"
    for ix in alsub:
        if ix in aledge:
            judge = 0
            for mx in range(len(atoms)):
                if atoms[mx][0] == "eow":
                    judge = 1
                    break
            if judge:
                atoms[ix][0] = "emgo"
            print("substi inedge", ix)


def writeinp(direction, addmol):
    # direction must be z or y, for z we insert the water into the interlayer, for y we insert the water
    # hanging up y parallel to the xz plane
    global box
    global atoms
    arr_xx = np.array(np.array(atoms)[:, 1], dtype=float)
    xx_sort = np.sort(arr_xx)
    arr_yy = np.array(np.array(atoms)[:, 2], dtype=float)
    yy_sort = np.sort(arr_yy)
    arr_zz = np.array(np.array(atoms)[:, 3], dtype=float)
    zz_sort = np.sort(arr_zz)
    xlo = (xx_sort[0] + xx_sort[len(xx_sort) - 1]) / 2 - box[0] / 2
    xhi = (xx_sort[0] + xx_sort[len(xx_sort) - 1]) / 2 + box[0] / 2
    yhi = (yy_sort[0] + yy_sort[len(yy_sort) - 1]) / 2 + box[1] / 2
    if direction == "z":
        deltaz = 0
        for ix in range(len(addmol)):
            deltaz += addmol[ix][3]*addmol[ix][1]*10/(6.02*(box[0]-2)*(box[1]-2)*addmol[ix][2])
        fin_xlo = xlo + 1
        fin_xhi = xhi - 1
        fin_ylo = yy_sort[0] + 1
        fin_yhi = yy_sort[len(yy_sort) - 1] - 1
        fin_zlo = zz_sort[len(zz_sort) - 1] + 2
        fin_zhi = zz_sort[len(zz_sort) - 1] + 2 + deltaz
        box[2] = fin_zhi - zz_sort[0]
    elif direction == "y":
        deltay = 0
        for ix in range(len(addmol)):
            deltay += addmol[ix][3]*addmol[ix][1]*10/(6.02*(box[0]-2)*box[2]*addmol[ix][2])
        fin_xlo = xlo + 1
        fin_xhi = xhi - 1
        fin_ylo = yhi + 1
        fin_yhi = yhi + 1 + deltay
        fin_zlo = zz_sort[0]
        fin_zhi = zz_sort[0] + box[2]
        box[1] += deltay
    file2 = open(direction+".inp", "w")
    print("tolerance 2.0", file=file2)
    print("output "+direction+".xyz", file=file2)
    for ix in range(len(addmol)):
        print("filetype xyz", file=file2)
        print("structure "+addmol[ix][0]+".xyz", file=file2)
        print("number ", addmol[ix][3], file=file2)
        print("inside box", fin_xlo, fin_ylo, fin_zlo, fin_xhi, fin_yhi, fin_zhi, file=file2)
        print("end structure", file=file2)
    file2.close()


def readfile2(filenm):
    global atoms
    file1 = open(filenm)
    file1.readline()
    file1.readline()
    content = file1.readlines()
    file1.close()
    for ix in content:
        tmp_1 = ix.split()[0]
        tmp_2 = float(ix.split()[1])
        tmp_3 = float(ix.split()[2])
        tmp_4 = float(ix.split()[3])
        atoms.append([tmp_1, tmp_2, tmp_3, tmp_4])


def single(filenm, leng):
    global mol_bond
    global mol_ang
    global bond_coeff
    global ang_coeff
    global dih_coeff
    atoms_t = []
    type = []
    pair = []
    ang_p = []
    dih_p = []
    imp_p = []
    bond_num = {}
    file2 = open(filenm)
    file2.readline()
    file2.readline()
    content = file2.readlines()
    file2.close()
    for ix in content:
        tmp_1 = ix.split()[0]
        tmp_2 = float(ix.split()[1])
        tmp_3 = float(ix.split()[2])
        tmp_4 = float(ix.split()[3])
        type.append(tmp_1)
        atoms_t.append([tmp_2, tmp_3, tmp_4])
    for ix in range(len(atoms_t)):
        for jx in range(ix+1, len(atoms_t)):
            if caldis(atoms_t[ix], atoms_t[jx], leng) == 0:
                if type[ix]+'-'+type[jx] in bond_coeff:
                    pair.append([type[ix]+'-'+type[jx], ix, jx])
                elif type[jx]+'-'+type[ix] in bond_coeff:
                    pair.append([type[jx]+'-'+type[ix], jx, ix])
                else:
                    print("no such bond coeff")
    for ix in range(len(pair)):
        if pair[ix][1] in bond_num:
            tmp = bond_num[pair[ix][1]]
            tmp.append(pair[ix][2])
            bond_num[pair[ix][1]] = tmp
        else:
            bond_num[pair[ix][1]] = [pair[ix][2]]
        if pair[ix][2] in bond_num:
            tmp = bond_num[pair[ix][2]]
            tmp.append(pair[ix][1])
            bond_num[pair[ix][2]] = tmp
        else:
            bond_num[pair[ix][2]] = [pair[ix][1]]
    for ix in bond_num:
        if len(bond_num[ix]) > 1:
            for jx in range(len(bond_num[ix])):
                for kx in range(jx + 1, len(bond_num[ix])):
                    if type[bond_num[ix][jx]]+'-'+type[ix]+'-'+type[bond_num[ix][kx]] in ang_coeff:
                        ang_p.append([type[bond_num[ix][jx]]+'-'+type[ix]+'-'+type[bond_num[ix][kx]],
                                      bond_num[ix][jx], ix, bond_num[ix][kx]])
                    elif type[bond_num[ix][kx]]+'-'+type[ix]+'-'+type[bond_num[ix][jx]] in ang_coeff:
                        ang_p.append([type[bond_num[ix][kx]]+'-'+type[ix]+'-'+type[bond_num[ix][jx]],
                                      bond_num[ix][kx], ix, bond_num[ix][jx]])
                    else:
                        print("no such ang coeff")
    for ix in pair:
        if len(bond_num[ix[1]]) > 1 and len(bond_num[ix[2]]) > 1:
            numn = 0
            for jx in bond_num[ix[1]]:
                for kx in bond_num[ix[2]]:
                    if jx != ix[2] and kx != ix[1]:
                        if type[jx]+'-'+type[ix[1]]+'-'+type[ix[2]]+'-'+type[kx] in dih_coeff:
                            dih_p.append([type[jx]+'-'+type[ix[1]]+'-'+type[ix[2]]+'-'+type[kx],
                                         jx, ix[1], ix[2], kx])
                            numn += 1
                        elif type[kx]+'-'+type[ix[2]]+'-'+type[ix[1]]+'-'+type[jx] in dih_coeff:
                                dih_p.append([type[jx]+'-'+type[ix[1]]+'-'+type[ix[2]]+'-'+type[kx],
                                              jx, ix[1], ix[2], kx])
                                numn += 1
                        else:
                            print("no such dih coeff")
    for ix in pair:
        ix[1] += 1
        ix[2] += 1
    for ix in ang_p:
        ix[1] += 1
        ix[2] += 1
        ix[3] += 1
    for ix in dih_p:
        ix[1] += 1
        ix[2] += 1
        ix[3] += 1
        ix[4] += 1
    mol_bond.append(pair)
    mol_ang.append(ang_p)
    mol_dih.append(dih_p)
    mol_imp.append([])


def calid():
    global atoms
    global molnum
    global moltimes
    global moltot
    global ang
    for ix in range(len(molnum)):
        tmp = 0
        for jx in range(ix+1):
            tmp += molnum[jx]*moltimes[jx]
        moltot.append(tmp)
    for ix in range(len(atoms)):
        if ix < moltot[0]:
            atoms[ix].append(1)
        elif ix >= moltot[0]:
            for jx in range(1, len(molnum)):
                if moltot[jx] > ix >= moltot[jx - 1]:
                    molexist = 0
                    for kx in range(jx):
                        molexist += moltimes[kx]
                    molid = (ix - moltot[jx - 1]) // molnum[jx] + 1 + molexist
                    atoms[ix].append(molid)


def caltype():
    global atom_turn
    global bond_turn
    global ang_turn
    global atoms
    global bond
    global ang
    for ix in range(len(atoms)):
        if atoms[ix][0] in atom_turn:
            pass
        else:
            atom_turn.append(atoms[ix][0])
    for ix in range(len(bond)):
        if bond[ix][0] in bond_turn:
            pass
        else:
            bond_turn.append(bond[ix][0])
    for ix in range(len(ang)):
        if ang[ix][0] in ang_turn:
            pass
        else:
            ang_turn.append(ang[ix][0])
    for ix in range(len(dih)):
        if dih[ix][0] in dih_turn:
            pass
        else:
            dih_turn.append(dih[ix][0])
    for ix in range(len(imp)):
        if imp[ix][0] in imp_turn:
            pass
        else:
            imp_turn.append(imp[ix][0])


def trandata(filenm):
    global mass
    global char_dic
    global pair_coeff
    global bond_coeff
    global ang_coeff
    global dih_coeff
    global imp_coeff
    global mol_bond
    global mol_ang
    global mol_dih
    global mol_imp
    mass_tmp = {}
    pair_tmp = {}
    bond_tmp = {}
    ang_tmp = {}
    dih_tmp = {}
    imp_tmp = {}
    bond_p = []
    ang_p = []
    dih_p = []
    imp_p = []
    molxyz = []
    file2 = open("./"+filenm+".data")
    file2.readline()
    file2.readline()
    file2.readline()
    whe_bond = file2.readline().split()[0]
    whe_ang = file2.readline().split()[0]
    whe_dih = file2.readline().split()[0]
    whe_imp = file2.readline().split()[0]
    while 1:
        tt = file2.readline()
        if tt == "Masses\n":
            break
    file2.readline()
    while 1:
        tt = file2.readline().split()
        if tt == []:
            break
        else:
            mass[tt[3]] = tt[1]
            mass_tmp[tt[0]] = tt[3]
    file2.readline()
    file2.readline()
    while 1:
        tt = file2.readline().split()
        if tt == []:
            break
        else:
            pair_coeff[tt[4]] = [tt[1], tt[2]]
            pair_tmp[tt[0]] = tt[4]
    if whe_bond:
        file2.readline()
        file2.readline()
        while 1:
            tt = file2.readline().split()
            if tt == []:
                break
            else:
                bond_coeff[tt[4]] = [tt[1], tt[2]]
                bond_tmp[tt[0]] = tt[4]
    if whe_ang:
        file2.readline()
        file2.readline()
        while 1:
            tt = file2.readline().split()
            if tt == []:
                break
            else:
                ang_coeff[tt[4]] = [tt[1], tt[2]]
                ang_tmp[tt[0]] = tt[4]
    if whe_dih:
        file2.readline()
        file2.readline()
        while 1:
            tt = file2.readline().split()
            if tt == []:
                break
            else:
                dih_coeff[tt[5]] = [tt[1], tt[2], tt[3]]
                dih_tmp[tt[0]] = tt[5]
    if whe_imp:
        file2.readline()
        file2.readline()
        while 1:
            tt = file2.readline().split()
            if tt == []:
                break
            else:
                imp_coeff[tt[5]] = [tt[1], tt[2], tt[3]]
                imp_tmp[tt[0]] = tt[5]
    file2.readline()
    file2.readline()
    while 1:
        tt = file2.readline().split()
        if tt == []:
            break
        else:
            molxyz.append([tt[11], tt[4], tt[5], tt[6]])
            char_dic[tt[11]] = tt[3]
    if whe_bond:
        file2.readline()
        file2.readline()
        while 1:
            tt = file2.readline().split()
            if tt == []:
                break
            else:
                bond_p.append([bond_tmp[tt[1]], int(tt[2]), int(tt[3])])
    if whe_ang:
        file2.readline()
        file2.readline()
        while 1:
            tt = file2.readline().split()
            if tt == []:
                break
            else:
                ang_p.append([ang_tmp[tt[1]], int(tt[2]), int(tt[3]), int(tt[4])])
    if whe_dih:
        file2.readline()
        file2.readline()
        while 1:
            tt = file2.readline().split()
            if tt == []:
                break
            else:
                dih_p.append([dih_tmp[tt[1]], int(tt[2]), int(tt[3]), int(tt[4]), int(tt[5])])
    if whe_imp:
        file2.readline()
        file2.readline()
        while 1:
            tt = file2.readline().split()
            if tt == []:
                break
            else:
                imp_p.append([imp_tmp[tt[1]], int(tt[2]), int(tt[3]), int(tt[4]), int(tt[5])])
    file2.close()
    mol_bond.append(bond_p)
    mol_ang.append(ang_p)
    mol_dih.append(dih_p)
    mol_imp.append(imp_p)
    file2 = open(filenm+".xyz", "w")
    print(len(molxyz), file=file2)
    print(len(molxyz), file=file2)
    for ix in range(len(molxyz)):
        print("%s %.6f %.6f %.6f" % (molxyz[ix][0], float(molxyz[ix][1]),
                                     float(molxyz[ix][2]), float(molxyz[ix][3])), file=file2)
    file2.close()


readmon()
tranohho()
subti(alfin=alsubt, sifin=sisubt)


if len(add_z1):
    for ix in range(len(add_z1)):
        molname.append(add_z1[ix][0])
        whemol.append(1)
    writeinp(direction="z", addmol=add_z1)
    for ix in range(len(add_z1)):
        file2 = open(add_z1[ix][0]+'.xyz')
        molnum.append(len(file2.readlines())-2)
        moltimes.append(add_z1[ix][3])
    os.system('~/apps/packmol/packmol < z.inp > /dev/null')
    readfile2("z.xyz")
    for ix in range(len(add_z1)):
        if os.path.exists(add_z1[ix][0]+'.data'):
            trandata(add_z1[ix][0])
        elif os.path.exists(add_z1[ix][0]+'.xyz'):
            single(add_z1[ix][0]+'.xyz', 1.6)
        else:
            print("no such file:", molname[ix])
    os.system('rm z.inp')
    os.system('rm z.xyz')

if len(add_z2):
    readmon2(addzlen)
    tranohho()
    subti(alfin=alsubt, sifin=sisubt)
    for ix in range(len(add_z2)):
        molname.append(add_z2[ix][0])
        whemol.append(1)
    writeinp(direction="z", addmol=add_z2)
    for ix in range(len(add_z2)):
        file2 = open(add_z2[ix][0]+'.xyz')
        molnum.append(len(file2.readlines())-2)
        moltimes.append(add_z2[ix][3])
    os.system('~/apps/packmol/packmol < z.inp > /dev/null')
    readfile2("z.xyz")
    for ix in range(len(add_z2)):
        if os.path.exists(add_z2[ix][0]+'.data'):
            trandata(add_z2[ix][0])
        elif os.path.exists(add_z2[ix][0]+'.xyz'):
            single(add_z2[ix][0]+'.xyz', 1.6)
        else:
            print("no such file:", molname[ix])
    os.system('rm z.inp')
    os.system('rm z.xyz')


if 'add_z3' in dir():
    readmon2(addzlen*2)
    tranohho()
    subti(alfin=alsubt, sifin=sisubt)
    for ix in range(len(add_z3)):
        molname.append(add_z3[ix][0])
        whemol.append(1)
    writeinp(direction="z", addmol=add_z3)
    for ix in range(len(add_z3)):
        file2 = open(add_z3[ix][0]+'.xyz')
        molnum.append(len(file2.readlines())-2)
        moltimes.append(add_z3[ix][3])
    os.system('~/apps/packmol/packmol < z.inp > /dev/null')
    readfile2("z.xyz")
    for ix in range(len(add_z3)):
        if os.path.exists(add_z3[ix][0]+'.data'):
            trandata(add_z3[ix][0])
        elif os.path.exists(add_z3[ix][0]+'.xyz'):
            single(add_z3[ix][0]+'.xyz', 1.6)
        else:
            print("no such file:", molname[ix])
    os.system('rm z.inp')
    os.system('rm z.xyz')


if 'add_z4' in dir():
    readmon2(addzlen*3)
    tranohho()
    subti(alfin=alsubt, sifin=sisubt)
    for ix in range(len(add_z4)):
        molname.append(add_z4[ix][0])
        whemol.append(1)
    writeinp(direction="z", addmol=add_z3)
    for ix in range(len(add_z4)):
        file2 = open(add_z4[ix][0]+'.xyz')
        molnum.append(len(file2.readlines())-2)
        moltimes.append(add_z4[ix][3])
    os.system('~/apps/packmol/packmol < z.inp > /dev/null')
    readfile2("z.xyz")
    for ix in range(len(add_z4)):
        if os.path.exists(add_z4[ix][0]+'.data'):
            trandata(add_z4[ix][0])
        elif os.path.exists(add_z4[ix][0]+'.xyz'):
            single(add_z4[ix][0]+'.xyz', 1.6)
        else:
            print("no such file:", molname[ix])
    os.system('rm z.inp')
    os.system('rm z.xyz')


if len(add_y):
    writeinp(direction="y", addmol=add_y)
    os.system('~/apps/packmol/packmol < y.inp > /dev/null')
    readfile2("y.xyz")
    for ix in range(len(add_y)):
        molname.append(add_y[ix][0])
        whemol.append(1)
    for ix in range(len(add_y)):
        file2 = open(add_y[ix][0]+'.xyz')
        molnum.append(len(file2.readlines())-2)
        moltimes.append(add_y[ix][3])
    for ix in range(len(add_y)):
        if os.path.exists(add_y[ix][0]+'.data'):
            trandata(add_y[ix][0])
        elif os.path.exists(add_y[ix][0]+'.xyz'):
            single(add_y[ix][0]+'.xyz', 1.6)
        else:
            print("no such file:", molname[ix])
    os.system('rm y.inp')
    os.system('rm y.xyz')


findang()
calid()
findwatmon()

for mm in range(0, len(molnum)):
    if whemol[mm]:
        for ix in range(moltimes[mm]):
            for jx in range(len(mol_bond[mm])):
                xx = molnum[mm]
                exit = moltot[mm-1]
                t_1 = mol_bond[mm][jx][0]
                t_2 = mol_bond[mm][jx][1] + exit + xx * ix
                t_3 = mol_bond[mm][jx][2] + exit + xx * ix
                bond.append([t_1, t_2, t_3])
            for jx in range(len(mol_ang[mm])):
                xx = molnum[mm]
                exit = moltot[mm-1]
                t_1 = mol_ang[mm][jx][0]
                t_2 = mol_ang[mm][jx][1] + exit + xx * ix
                t_3 = mol_ang[mm][jx][2] + exit + xx * ix
                t_4 = mol_ang[mm][jx][3] + exit + xx * ix
                ang.append([t_1, t_2, t_3, t_4])
            for jx in range(len(mol_dih[mm])):
                xx = molnum[mm]
                exit = moltot[mm-1]
                t_1 = mol_dih[mm][jx][0]
                t_2 = mol_dih[mm][jx][1] + exit + xx * ix
                t_3 = mol_dih[mm][jx][2] + exit + xx * ix
                t_4 = mol_dih[mm][jx][3] + exit + xx * ix
                t_5 = mol_dih[mm][jx][4] + exit + xx * ix
                dih.append([t_1, t_2, t_3, t_4, t_5])
            for jx in range(len(mol_imp[mm])):
                xx = molnum[mm]
                exit = moltot[mm-1]
                t_1 = mol_imp[mm][jx][0]
                t_2 = mol_imp[mm][jx][1] + exit + xx * ix
                t_3 = mol_imp[mm][jx][2] + exit + xx * ix
                t_4 = mol_imp[mm][jx][3] + exit + xx * ix
                t_5 = mol_imp[mm][jx][4] + exit + xx * ix
                imp.append([t_1, t_2, t_3, t_4, t_5])


caltype()


file2 = open("finallmp", "w")
print("the monte", file=file2)
print(str(len(atoms))+" atoms", file=file2)
print(str(len(bond))+" bonds", file=file2)
print(str(len(ang))+" angles", file=file2)
print(str(len(dih))+" dihedrals", file=file2)
print(str(len(imp))+" impropers", file=file2)
print(str(len(atom_turn))+" atom types", file=file2)
print(str(len(bond_turn))+" bond types", file=file2)
print(str(len(ang_turn))+" angle types", file=file2)
print(str(len(dih_turn))+" dihedral types", file=file2)
print(str(len(imp_turn))+" improper types\n", file=file2)
print("%.5f %.5f %s" % (0, box[0], " xlo xhi"), file=file2)
if add_y == []:
    print("%.5f %.5f %s" % (0, box[1], " ylo yhi"), file=file2)
else:
    print("%.5f %.5f %s" % (0, box[1]+4, " ylo yhi"), file=file2)
print("%.5f %.5f %s" % (0, box[2]+3, " zlo zhi"), file=file2)
print("\nMasses\n", file=file2)
for i in range(len(atom_turn)):
    print("%d %.6f %s" % (i+1, float(mass[atom_turn[i]]), "#"+atom_turn[i]),  file=file2)
print("\nPair Coeffs\n", file=file2)
for i in range(len(atom_turn)):
    print("%d %.10f %.10f %s" % (i+1, float(pair_coeff[atom_turn[i]][0]), float(pair_coeff[atom_turn[i]][1]),
          "#"+atom_turn[i]), file=file2)
print("\nBond Coeffs\n", file=file2)
for i in range(len(bond_turn)):
    print("%-3d %-10s " % (i+1, bond_coeff[bond_turn[i]][0]), end='', file=file2)
    for jx in range(1, len(bond_coeff[bond_turn[i]])):
        print("%-.10f " % float(bond_coeff[bond_turn[i]][jx]), end=' ', file=file2)
    print("#"+bond_turn[i], file=file2)
    # print("%d %.10f %.10f %s" % (i+1, float(bond_coeff[bond_turn[i]][0]), float(bond_coeff[bond_turn[i]][1]),
    #       "#"+bond_turn[i]), file=file2)
print("\nAngle Coeffs\n", file=file2)
for i in range(len(ang_turn)):
    print("%d %.10f %.10f %s" % (i+1, float(ang_coeff[ang_turn[i]][0]), float(ang_coeff[ang_turn[i]][1]),
          "#"+ang_turn[i]), file=file2)
if len(dih) > 0:
    print("\nDihedral Coeffs\n", file=file2)
    for i in range(len(dih_turn)):
        print("%d %.10f %.10f %.10f %s" % (i+1, float(dih_coeff[dih_turn[i]][0]),
                                           float(dih_coeff[dih_turn[i]][1]), float(dih_coeff[dih_turn[i]][2]),
                                           "#"+dih_turn[i]), file=file2)
if len(imp) > 0:
    print("\nImproper Coeffs\n", file=file2)
    for i in range(len(imp_turn)):
        print("%d %.10f %.10f %.10f %s" % (i + 1, float(imp_coeff[imp_turn[i]][0]),
                                           float(imp_coeff[imp_turn[i]][1]),
                                           float(imp_coeff[imp_turn[i]][2]), "#" + imp_turn[i]), file=file2)
print("\nAtoms\n", file=file2)
for i in range(len(atoms)):
    print("%d %d  %d %.6f %.6f %.6f %.6f " % (i + 1, atoms[i][4], atom_turn.index(atoms[i][0])+1,
                                             float(char_dic[atoms[i][0]]), float(atoms[i][1]),
                                             float(atoms[i][2]),
                                              float(atoms[i][3])), "#", atoms[i][0], file=file2)
print("\nBonds\n", file=file2)
for i in range(len(bond)):
    print("%d %d  %d %d" % (i + 1, bond_turn.index(bond[i][0])+1, bond[i][1], bond[i][2]), file=file2)
print("\nAngles\n", file=file2)
for i in range(len(ang)):
    print("%d %d  %d %d %d" % (i + 1, ang_turn.index(ang[i][0])+1,
                               ang[i][1], ang[i][2], ang[i][3]), file=file2)
if len(dih) > 0:
    print("\nDihedrals\n", file=file2)
    for i in range(len(dih)):
        print("%d %d  %d %d %d %d" % (i + 1, dih_turn.index(dih[i][0]) + 1,
                                      dih[i][1], dih[i][2], dih[i][3], dih[i][4]), file=file2)
if len(imp) > 0:
    print("\nImpropers\n", file=file2)
    for i in range(len(imp)):
        print("%d %d  %d %d %d %d" % (i + 1, imp_turn.index(imp[i][0]) + 1,
                                      imp[i][1], imp[i][2], imp[i][3], imp[i][4]), file=file2)
file2.close()
