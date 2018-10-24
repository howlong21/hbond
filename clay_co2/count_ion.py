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
        self.iontype = ["mg", "k", "ca", "na", "ba"]
        self.arr_water = np.array([])
        self.arr_methane = np.array([])
        self.arr_clayo = np.array([])
        self.arr_box = np.array([])
        self.arr_ion = []
        self.readcoor()

    def read_box(self):
        """read the box from lammps simulation trajectory"""
        for _ in range(4):
            self.file.readline()
        box = [float(jx.split()[1]) - float(jx.split()[0]) for jx in [self.file.readline() for _ in range(3)]]
        self.file.readline()
        return box

    def read_xyz(self, eleme1, eleme2, clay, clayao, perion):
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
                    elif linesep[0] in self.iontype:
                        perion.append([linesep[0], float(linesep[1]), float(linesep[2]), float(linesep[3])])
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
        ions = []
        perion = []
        self.file.readline()
        while 1:
            box1 = self.read_box()
            line = self.read_xyz(eleme1, eleme2, clay, clayao, perion)
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
            ions.append(perion)
            if line:
                eleme1 = []
                eleme2 = []
                clay = []
                clayao = []
                perion = []
            else:
                self.arr_water = np.array(ele1s)
                self.arr_methane = np.array(ele2s)
                self.arr_clayo = np.array(clays)
                self.arr_box = np.array(boxs)
                self.arr_ion = ions
                break
        self.file.close()
        print("all coordinates has been read...")


class Count(Basic_Calcul):
    def __init__(self, lmpdata):
        self.dumpdata = lmpdata
        self.left = 10
        self.right = 50
        self.iontype = self.find_ion()
        self.count_ion(0)

    def find_ion(self):
        iontype = []
        for ix in range(len(self.dumpdata.arr_ion[0])):
            iontype.append(self.dumpdata.arr_ion[0][ix][0])
        return list(set(iontype))

    def find_bound(self, time):
        sort_ob = np.sort(self.dumpdata.arr_clayo[time, :, 2])
        indexcal = []
        layers = []
        while len(indexcal) < len(sort_ob):
            perlayer = []
            tmp = 0
            for ix in range(len(sort_ob)):
                if ix not in indexcal:
                    tmp = sort_ob[ix]
            for ix in range(len(sort_ob)):
                if np.abs(self.get_vector(sort_ob[ix], tmp, self.dumpdata.arr_box[time][2])) < 3:
                    perlayer.append(sort_ob[ix])
                    indexcal.append(ix)
            layers.append(perlayer)
        layer1min = min(layers[0])
        layer1max = max(layers[0])
        if layer1max - layer1min > self.dumpdata.arr_box[time][2]/2:
            for ix in range(len(layers[0])):
                if layers[0][ix] < layer1max - self.dumpdata.arr_box[time][2]/2:
                    layers[0][ix] += self.dumpdata.arr_box[time][2]
        layer_ave = []
        for ix in range(len(layers)):
            layer_ave.append(np.average(np.array([layers[ix]])))
        res_bound = []
        if math.fabs(layer_ave[1] - layer_ave[0]) < 7.5:
            layer_ave.append(layer_ave[0])
            layer_ave.remove(layer_ave[0])
        while len(res_bound) < len(layer_ave)/2:
            numbound = len(res_bound)
            res_bound.append([layer_ave[numbound*2], layer_ave[numbound*2+1]])
        return res_bound

    def count_ion(self, time):
        claybound = self.find_bound(time)
        fin_count = np.zeros([len(self.iontype), len(claybound)], int)
        for ix in range(len(self.dumpdata.arr_ion[time])):
            if self.left < self.dumpdata.arr_ion[time][ix][2] < self.right:
                for jx in range(len(self.iontype)):
                    if self.dumpdata.arr_ion[time][ix][0] == self.iontype[jx]:
                        for kx in range(len(claybound)):
                            layer_dis = math.fabs(self.get_vector(claybound[kx][0], claybound[kx][1],
                                                                  self.dumpdata.arr_box[time][2]))
                            updis = math.fabs(self.get_vector(self.dumpdata.arr_ion[time][ix][3], claybound[kx][1],
                                                              self.dumpdata.arr_box[time][2]))
                            downdis = math.fabs(self.get_vector(self.dumpdata.arr_ion[time][ix][3], claybound[kx][0],
                                                                self.dumpdata.arr_box[time][2]))
                            if updis < layer_dis and downdis < layer_dis:
                                fin_count[jx][kx] += 1
        for ix in range(len(fin_count)):
            print(self.iontype[ix], fin_count[ix])


traj = GetLmpData('mixdump')
t1 = Count(traj)

