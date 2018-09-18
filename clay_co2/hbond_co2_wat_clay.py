import math
import matplotlib.pyplot as plt
import numpy as np


class Basic_Calcul(object):
    @staticmethod
    def get_vector(coor1, coor2, box=100.0):
        tmp1 = coor2 - coor1 - box * np.floor((coor2 - coor1) / box)
        return tmp1 - round(tmp1 / box) * box

    def caldis(self, atom1, atom2, box, maxleng=3.5):
        tmp_arr = np.array([], dtype=np.float32)
        for ix in range(len(atom1)):
            tmpleng = np.abs(self.get_vector(atom1[ix], atom2[ix], box[ix]))
            if tmpleng > maxleng:
                return 0
            else:
                tmp_arr = np.append(tmp_arr, [np.float32(tmpleng)])
        finres = np.sqrt(tmp_arr.dot(tmp_arr))
        if finres <= maxleng:
            return finres
        else:
            return 0

    def calang(self, atom1, atom2, atom3, box, maxang=30):
        edgeBA = np.array([], dtype=np.float32)
        edgeBC = np.array([], dtype=np.float32)
        for ix in range(len(atom2)):
            edgeBA = np.append(edgeBA, [np.float32(self.get_vector(atom2[ix], atom1[ix], box[ix]))])
            edgeBC = np.append(edgeBC, [np.float32(self.get_vector(atom2[ix], atom3[ix], box[ix]))])
        rescos = edgeBA.dot(edgeBC) / np.sqrt(edgeBA.dot(edgeBA)) / np.sqrt(edgeBC.dot(edgeBC)) - 0.00001
        resang = 180 * math.acos(rescos) / math.pi
        if resang <= maxang:
            return resang
        else:
            return 0

    def triang_area(self, atom1, atom2, atom3, box):
        edgeBA = np.array([], dtype=np.float32)
        edgeBC = np.array([], dtype=np.float32)
        for ix in range(len(atom2)):
            edgeBA = np.append(edgeBA, [np.float32(self.get_vector(atom2[ix], atom1[ix], box[ix]))])
            edgeBC = np.append(edgeBC, [np.float32(self.get_vector(atom2[ix], atom3[ix], box[ix]))])
        shiftBA = np.array([edgeBA[1], edgeBA[2], edgeBA[0]])
        shiftBC = np.array([edgeBC[1], edgeBC[2], edgeBC[0]])
        final = edgeBA * shiftBC - edgeBC * shiftBA
        return np.sqrt(final.dot(final)) / 2

    def polygon_area(self, edges, time, arr_ele1, arrclayo, arr_box):
        fin_area = 0
        vertice = edges[0][0]
        for ix in range(len(edges)):
            if vertice not in edges[ix]:
                if vertice < len(arr_ele1[0]):
                    atom1 = arr_ele1[time][vertice][0]
                else:
                    atom1 = arrclayo[time][vertice-len(arr_ele1[0])]
                if edges[ix][0] < len(arr_ele1[0]):
                    atom2 = arr_ele1[time][edges[ix][0]][0]
                else:
                    atom2 = arrclayo[time][edges[ix][0]-len(arr_ele1[0])]
                if edges[ix][1] < len(arr_ele1[0]):
                    atom3 = arr_ele1[time][edges[ix][1]][0]
                else:
                    atom3 = arrclayo[time][edges[ix][1]-len(arr_ele1[0])]
                # fin_area += self.triang_area(arr_ele1[time][vertice][0], arr_ele1[time][edges[ix][0]][0],
                #                         arr_ele1[time][edges[ix][1]][0], arr_box[time])
                fin_area += self.triang_area(atom1, atom2, atom3, arr_box[time])
        return fin_area


class GetLmpData(Basic_Calcul):
    def __init__(self, filename):
        self.file = open(filename)
        self.ele1 = ["ow", "hw"]
        self.ele2 = ["oc"]
        self.clayo = ["ob", "obts"]
        self.ao = ["ao", "mgo"]
        self.arr_water = np.array([])
        self.arr_methane = np.array([])
        self.arr_clayo = np.array([])
        self.arr_box = np.array([])
        self.readcoor()

    def read_box(self):
        """read the box from lammps simulation trajectory"""
        for _ in range(4):
            self.file.readline()
        box = [float(jx.split()[1]) - float(jx.split()[0]) for jx in [self.file.readline() for _ in range(3)]]
        self.file.readline()
        return box

    def read_xyz(self, eleme1, eleme2, clay, clayao):
        """read the atom coordination from lammps trajectory"""
        num_h2o = 0
        water_o = [0] * 3
        while 1:
            line = self.file.readline()
            if line:
                linesep = line.split()
                if linesep[0] != "ITEM:":
                    if linesep[0] in self.ele1:
                        water_o[num_h2o % 3] = [float(linesep[1]), float(linesep[2]), float(linesep[3])]
                        num_h2o += 1
                        if num_h2o % 3 == 0:
                            eleme1.append([water_o[0], water_o[1], water_o[2]])
                    elif linesep[0] in self.ele2:
                        eleme2.append([float(linesep[1]), float(linesep[2]), float(linesep[3])])
                    elif linesep[0] in self.clayo:
                        clay.append([float(linesep[1]), float(linesep[2]), float(linesep[3])])
                    elif linesep[0] in self.ao:
                        clayao.append([float(linesep[1]), float(linesep[2]), float(linesep[3])])
                else:
                    return line
            else:
                return line

    def readcoor(self):
        """read the coordinate from the dump file,then save water in the array arr_water,
        save cations in list cations,save clay_atoms in arr_clay,
        the information of box size in array arr_box"""
        boxs = []
        ele1s = []
        ele2s = []
        clays = []
        eleme1 = []
        eleme2 = []
        clay = []
        clayao = []
        self.file.readline()
        while 1:
            box1 = self.read_box()
            line = self.read_xyz(eleme1, eleme2, clay, clayao)
            aozz = np.array(clayao)[:, 2]
            aozz.sort()
            allz = [aozz[0]]
            for ind1 in range(1, len(aozz)):
                gap1 = aozz[ind1] - aozz[ind1-1]
                if gap1 > 2:
                    allz.append(aozz[ind1])
            remove_mid = []
            for ind1 in range(len(clay)):
                whe_mid = 0
                for ind2 in range(len(allz)):
                    if -2.0 < self.get_vector(allz[ind2], clay[ind1][2], box1[2]) < 2.0:
                        whe_mid = 1
                        break
                if whe_mid:
                    remove_mid.append(clay[ind1])
            for ind1 in remove_mid:
                clay.remove(ind1)
            boxs.append(box1)
            ele1s.append(eleme1)
            ele2s.append(eleme2)
            clays.append(clay)
            if line:
                eleme1 = []
                eleme2 = []
                clay = []
                clayao = []
            else:
                self.arr_water = np.array(ele1s)
                self.arr_methane = np.array(ele2s)
                self.arr_clayo = np.array(clays)
                self.arr_box = np.array(boxs)
                break
        self.file.close()
        print("all coordinates has been read...")


class GetHbonds(Basic_Calcul):
    def __init__(self, getlmpdata, time):
        self.dumpdata = getlmpdata
        self.hboc = self.hb_wat_cardio(time)
        self.hbclayo = self.hb_wat_clay(time)
        self.bondleng = 3.5
        self.angleng = 45

    def hb_wat_cardio(self, time):
        hb_wat_ocs = []
        for ix in range(len(self.dumpdata.arr_methane[time])):
            for jx in range(len(self.dumpdata.arr_water[time])):
                whehb = self.whe_ow_oc(self.dumpdata.arr_water[time][jx], self.dumpdata.arr_methane[time][ix],
                                       self.dumpdata.arr_box[time])
                if whehb[0]:
                    hb_wat_ocs.append([whehb[0], whehb[1]])
        return np.array(hb_wat_ocs)

    def hb_wat_clay(self, time):
        hb_wat_clayos = []
        for ix in range(len(self.dumpdata.arr_clayo[time])):
            for jx in range(len(self.dumpdata.arr_water[time])):
                whehb = self.whe_ow_os(self.dumpdata.arr_water[time][jx], self.dumpdata.arr_clayo[time][ix],
                                       self.dumpdata.arr_box[time])
                if whehb[0]:
                    hb_wat_clayos.append([whehb[0], whehb[1]])
        return np.array(hb_wat_clayos)

    def whe_ow_oc(self, water1, carbon, box):
        if self.caldis(water1[0], carbon, box):
            # print("in distance")
            distan1 = self.caldis(water1[1], carbon, box, maxleng=2.5)
            distan2 = self.caldis(water1[2], carbon, box, maxleng=2.5)
            if distan1:
                angle1 = self.calang(water1[0], carbon, water1[1], box)
                if angle1:
                    # print("find hb1")
                    return distan1, angle1
            elif distan2:
                angle2 = self.calang(water1[0], carbon, water1[2], box)
                if angle2:
                    # print("find hb2")
                    return distan2, angle2
        return 0, 0

    def whe_ow_os(self, water1, clayo, box):
        if self.caldis(water1[0], clayo, box):
            # print("in distance2")
            distan1 = self.caldis(water1[1], clayo, box, maxleng=2.5)
            distan2 = self.caldis(water1[2], clayo, box, maxleng=2.5)
            if distan1:
                angle1 = self.calang(water1[0], clayo, water1[1], box)
                if angle1:
                    return distan1, angle1
            elif distan2:
                angle2 = self.calang(water1[0], clayo, water1[2], box)
                if angle2:
                    return distan2, angle2
        return 0, 0


def plt_hist(data, binnum, fignum):
    plt.figure(0)
    xx, yy = plt.hist(data, binnum)[:2]
    plt.figure(fignum)
    finx = [(yy[ix] + yy[ix-1])/2 for ix in range(1, len(yy))]
    plt.plot(finx, xx)

traj1 = GetLmpData(r"C:\Users\lq\Desktop\10dump")
allbondoc = []
allangoc = []
allbondcl = []
allangcl = []
for ix in range(len(traj1.arr_box)):
    hbonds = GetHbonds(traj1, ix)
    allbondoc.extend(hbonds.hboc[:, 0])
    allangoc.extend(hbonds.hboc[:, 1])
    allbondcl.extend(hbonds.hbclayo[:, 0])
    allangcl.extend(hbonds.hbclayo[:, 1])
arr_bondoc = np.array(allbondoc)
arr_angoc = np.array(allangoc)
arr_bondcl = np.array(allbondcl)
arr_angcl = np.array(allangcl)

plt_hist(arr_bondoc, 50, 1)
plt_hist(arr_angoc, 20, 2)
plt_hist(arr_bondcl, 50, 1)
plt_hist(arr_angcl, 20, 2)

plt.show()
