import math
import numpy as np


spath = "./tail10dump"   # the file path to read the trajectory
file = open(spath)
arr_ele1 = np.array([])  # array to store the central atom
# arr_box = np.array([])   # array to store the simulation box
ele1 = ["ow", "hw"]  # the type name of the water atom
arr_box = np.array([[27.088, 31.28, 29.440]])


def read_box(file):
    for ix in range(4):
        file.readline()
    box = [float(jx.split()[1]) - float(jx.split()[0]) for jx in [file.readline() for ix in range(3)]]
    file.readline()
    return box


def read_xyz(file, eleme1):
    global ele1
    num_h2o = 0
    water_o = [0] * 3
    while 1:
        line = file.readline()
        if line:
            linesep = line.split()
            if linesep[0] != "ITEM:":
                if linesep[0] in ele1:
                    water_o[num_h2o % 3] = [float(linesep[1]), float(linesep[2]), float(linesep[3])]
                    num_h2o += 1
                    if num_h2o % 3 == 0:
                        eleme1.append([water_o[0], water_o[1], water_o[2]])
            else:
                return line
        else:
            return line


def readcoor():  # read the coordinate from the dump file,then save water in the array arr_multi,
    # save cations in list cations,save clay_atoms in list clays,
    # the information of box size in array arr_box
    global arr_box
    global arr_ele1
    global file
    boxs = []
    ele1s = []
    file.readline()
    eleme1 = []
    while 1:
        box1 = read_box(file)
        line = read_xyz(file, eleme1)
        boxs.append(box1)
        ele1s.append(eleme1)
        if line:
            eleme1 = []
        else:
            arr_ele1 = np.array(ele1s)
            arr_box = np.array(boxs)
            break
    file.close()
    print("all coordinates has been read...")


def get_vector(coor1, coor2, box):
    tmp1 = coor2 - coor1 - box * np.floor((coor2 - coor1) / box)
    return tmp1 - round(tmp1 / box) * box


def caldis(atom1, atom2, box, maxleng=3.5):
    tmp_arr = np.array([], dtype=np.float32)
    for ix in range(len(atom1)):
        tmpleng = np.abs(get_vector(atom1[ix], atom2[ix], box[ix]))
        if tmpleng > maxleng:
            return 0
        else:
            tmp_arr = np.append(tmp_arr, [np.float32(tmpleng)])
    finres = np.sqrt(tmp_arr.dot(tmp_arr))
    if finres <= maxleng:
        return finres
    else:
        return 0


def calang(atom1, atom2, atom3, box, maxang=30):
    edgeBA = np.array([], dtype=np.float32)
    edgeBC = np.array([], dtype=np.float32)
    for ix in range(len(atom2)):
        edgeBA = np.append(edgeBA, [np.float32(get_vector(atom2[ix], atom1[ix], box[ix]))])
        edgeBC = np.append(edgeBC, [np.float32(get_vector(atom2[ix], atom3[ix], box[ix]))])
    print(edgeBA, edgeBC)
    rescos = edgeBA.dot(edgeBC) / np.sqrt(edgeBA.dot(edgeBA)) / np.sqrt(edgeBC.dot(edgeBC)) - 0.00001
    resang = 180 * math.acos(rescos) / math.pi
    if resang <= maxang:
        return resang
    else:
        return 0


def get_midpoint(atom1, atom2, box):
    edgeBA = np.array([], dtype=np.float32)
    for ix in range(len(atom1)):
        edgeBA = np.append(edgeBA, [np.float32(get_vector(atom1[ix], atom2[ix], box[ix]) / 2 + atom1[ix])])
    return edgeBA


def get_normal_vector(atom1, atom2, atom3, box):
    vector12 = []
    vector13 = []
    for ix in range(len(atom1)):
        vector12.append(get_vector(atom1[ix], atom2[ix], box[ix]))
        vector13.append(get_vector(atom1[ix], atom3[ix], box[ix]))
    matrix1 = np.mat([vector12[:2], vector13[:2]])
    if vector12[0]*vector13[1] == vector12[1]*vector13[0]:
        resvec = [-vector13[0] / vector12[0], 1, 0]
    else:
        resvec = matrix1.I * np.mat([vector12[2], vector13[2]]).T
        resvec = [float(resvec[0][0]), float(resvec[1][0]), -1]
    return np.array(resvec)


def get_cos(vector1, vector2):
    rescos = vector1.dot(vector2) / np.sqrt(vector1.dot(vector1) * vector2.dot(vector2))
    return rescos


def get_rang_side(vector12, vector13):
    """get the other right-angle side named vector23 besides vector12"""
    vector23 = []
    delta = vector12.dot(vector13)/vector12.dot(vector12)
    for ix in range(len(vector12)):
        vector23.append(vector13[ix]-delta*vector12[ix])
    return np.array(vector23)


def get_dih_cos(atom1, atom2, atom3, atom4, box):
    vector12 = []
    vector23 = []
    vector24 = []
    for ix in range(len(atom1)):
        vector23.append(get_vector(atom2[ix], atom3[ix], box[ix]))
        vector12.append(get_vector(atom2[ix], atom1[ix], box[ix]))
        vector24.append(get_vector(atom2[ix], atom4[ix], box[ix]))
    right_side1 = get_rang_side(np.array(vector23), np.array(vector12))
    right_side2 = get_rang_side(np.array(vector23), np.array(vector24))
    return get_cos(right_side1, right_side2)


def triang_area(atom1, atom2, atom3, box):
    edgeBA = np.array([], dtype=np.float32)
    edgeBC = np.array([], dtype=np.float32)
    for ix in range(len(atom2)):
        edgeBA = np.append(edgeBA, [np.float32(get_vector(atom2[ix], atom1[ix], box[ix]))])
        edgeBC = np.append(edgeBC, [np.float32(get_vector(atom2[ix], atom3[ix], box[ix]))])
    value_bottom = edgeBC.dot(edgeBC)
    normol_edge = edgeBA.dot(edgeBC)/value_bottom*edgeBC - edgeBA
    print(normol_edge)
    value_high = np.sqrt(normol_edge.dot(normol_edge))
    print(np.sqrt(value_bottom)*value_high/2)
    return value_high * np.sqrt(value_bottom) / 2


def polygon_area(edges, time):
    fin_area = 0
    vertice = edges[0][0]
    for ix in range(len(edges)):
        if vertice not in edges[ix]:
            fin_area += triang_area(arr_ele1[time][vertice][0], arr_ele1[time][edges[ix][0]][0],
                                    arr_ele1[time][edges[ix][1]][0], arr_box[time])
    return fin_area


triang_area([0, 0, 0], [1, 0, 0], [5, 2, 0], [10, 10, 10])


def read_ice():
    global arr_ele1
    file1 = open('./ice_ih')
    ele_dic = {"O1": "ow", "O2": "ow", "H1": "hw", "H2": "hw", "H3": "hw", "H4": "hw"}
    tmpwater = []
    while 1:
        eveatom = file1.readline()
        if eveatom:
            tmpwater.append([ele_dic[eveatom.split()[0]], float(eveatom.split()[1]),
                             float(eveatom.split()[2]), float(eveatom.split()[3])])
        else:
            break
    evetime = []
    for ix in range(len(tmpwater)):
        evewater = []
        if tmpwater[ix][0] == "ow":
            evewater.append(tmpwater[ix][1:])
            for jx in range(len(tmpwater)):
                if tmpwater[jx][0] == "hw" and \
                        caldis(np.array(tmpwater[ix][1:]), np.array(tmpwater[jx][1:]), arr_box[0], maxleng=1.0):
                    evewater.append(tmpwater[jx][1:])
            evetime.append(evewater)
    arr_ele1 = np.array([evetime])
    print(arr_ele1.shape)




# readcoor()

# read_ice()

# all_dihpair = []
# for mx in range(len(arr_ele1)):
#     oxygen_pair = []
#     dih_pair = []
#     for ix in range(len(arr_ele1[mx])):
#         for jx in range(ix+1, len(arr_ele1[mx])):
#             if caldis(arr_ele1[mx][ix][0], arr_ele1[mx][jx][0], arr_box[mx]):
#                 oxygen_pair.append([ix, jx])
#     print("mx", mx)
#     for ix in range(len(oxygen_pair)):
#     # for ix in range(10):
#         tmp_dih = []
#         midpoint = get_midpoint(arr_ele1[mx][oxygen_pair[ix][0]][0], arr_ele1[mx][oxygen_pair[ix][1]][0], arr_box[0])
#         dis1 = caldis(midpoint, arr_ele1[mx][oxygen_pair[ix][0]][1], arr_box[mx])
#         dis2 = caldis(midpoint, arr_ele1[mx][oxygen_pair[ix][0]][2], arr_box[mx])
#         dis3 = caldis(midpoint, arr_ele1[mx][oxygen_pair[ix][1]][1], arr_box[mx])
#         dis4 = caldis(midpoint, arr_ele1[mx][oxygen_pair[ix][1]][2], arr_box[mx])
#         if dis1 > dis2:
#             tmp_dih = [1, oxygen_pair[ix][0], oxygen_pair[ix][1]]
#         else:
#             tmp_dih = [2, oxygen_pair[ix][0], oxygen_pair[ix][1]]
#         if dis3 > dis4:
#             tmp_dih.append(1)
#         else:
#             tmp_dih.append(2)
#         dih_pair.append(tmp_dih)
#     all_dihpair.append(dih_pair)
#
# sum_dih = 0
# allnum = 0
# for mx in range(len(all_dihpair)):
#     for ix in range(len(all_dihpair[mx])):
#         # ooh1 = get_normal_vector(arr_ele1[mx][all_dihpair[mx][ix][1]][0], arr_ele1[mx][all_dihpair[mx][ix][2]][0],
#         #                          arr_ele1[mx][all_dihpair[mx][ix][1]][all_dihpair[mx][ix][0]], arr_box[mx])
#         # ooh2 = get_normal_vector(arr_ele1[mx][all_dihpair[mx][ix][1]][0], arr_ele1[mx][all_dihpair[mx][ix][2]][0],
#         #                          arr_ele1[mx][all_dihpair[mx][ix][2]][all_dihpair[mx][ix][3]], arr_box[mx])
#         # print(get_cos(ooh1, ooh2))
#         tmp = get_dih_cos(arr_ele1[mx][all_dihpair[mx][ix][1]][all_dihpair[mx][ix][0]],
#                           arr_ele1[mx][all_dihpair[mx][ix][1]][0], arr_ele1[mx][all_dihpair[mx][ix][2]][0],
#                           arr_ele1[mx][all_dihpair[mx][ix][2]][all_dihpair[mx][ix][3]], arr_box[mx])
#         sum_dih += tmp
#         allnum += 1
#
# print("sum_dih:", sum_dih/allnum)


