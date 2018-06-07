import math
import numpy as np
import copy
import matplotlib.pyplot as plt


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
    arr_water = np.array([])
    arr_methane = np.array([])
    arr_clayo = np.array([])
    arr_box = np.array([])

    def __init__(self, filename):
        self.file = open(filename)
        self.ele1 = ["ow", "hw"]
        self.ele2 = ["co"]
        self.clayo = ["ob", "obts"]
        self.ao = ["ao", "mgo"]
        self.readcoor()

    def read_box(self):
        """read the box from lammps simulation trajectory"""
        for ix in range(4):
            self.file.readline()
        box = [float(jx.split()[1]) - float(jx.split()[0]) for jx in [self.file.readline() for ix in range(3)]]
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
            tmp1 = aozz[0]
            allz = [aozz[0]]
            for ind1 in range(1, len(aozz)):
                gap1 = aozz[ind1] - tmp1
                tmp1 = aozz[ind1]
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


class Node(object):
    # 初始化一个节点
    def __init__(self, value=None):
        self.value = value  # 节点值
        self.child_list = []    # 子节点列表

    # 添加一个孩子节点
    def add_child(self, node):
        self.child_list.append(node)


class GetHbonds(Basic_Calcul):
    def __init__(self, getlmpdata, guestnum, time):
        self.Hblist = self.per_guest(guestnum, time, getlmpdata.arr_water, getlmpdata.arr_methane,
                                              getlmpdata.arr_clayo, getlmpdata.arr_box)

    def per_guest(self, guestnum, time, arrwater, arrmeth, arrclay, arrbox):
        waterlist, claylist = self.get_watlist(guestnum, time, arrwater, arrmeth, arrclay, arrbox)
        per_bond = []
        extra_clay = []
        for jx in range(len(waterlist)):
            for kx in range(jx + 1, len(waterlist)):
                if self.whe_2_waters(arrwater[time][waterlist[jx]], arrwater[time][waterlist[kx]], arrbox[time]):
                    if waterlist[jx] < waterlist[kx]:
                        per_bond.append([waterlist[jx], waterlist[kx]])
                    else:
                        per_bond.append([waterlist[kx], waterlist[jx]])
            for kx in range(len(claylist)):
                if self.whe_ow_os(arrwater[time][waterlist[jx]], arrclay[time][claylist[kx]], arrbox[time]):
                    per_bond.append([waterlist[jx], len(arrwater[0]) + claylist[kx]])
                    extra_clay.append(claylist[kx])
                    # print(arrwater[time][waterlist[jx]], arrclay[time][claylist[kx]])
        if len(extra_clay) > 1:
            z_dis = np.array([np.abs(self.get_vector(arrclay[time][extra_clay[per_o]][2],
                                                     arrclay[time][extra_clay[0]][2]))
                              for per_o in range(len(extra_clay))])
            # print(z_dis > 2.0)
            # print(extra_clay[z_dis > 2.0])

            extra_clay2 = []
            for jx in extra_clay:
                if jx not in extra_clay2:
                    extra_clay2.append(jx)
            for jx in range(1, len(extra_clay2)):
                arrclay[time][extra_clay2[jx]][0] = arrclay[time][extra_clay2[0]][0] \
                                              + self.get_vector(arrclay[time][extra_clay2[0]][0],
                                                                arrclay[time][extra_clay2[jx]][0], arrbox[time][0])
                arrclay[time][extra_clay2[jx]][1] = arrclay[time][extra_clay2[0]][1] \
                                              + self.get_vector(arrclay[time][extra_clay2[0]][1],
                                                                arrclay[time][extra_clay2[jx]][1], arrbox[time][1])
            extra_edge = self.get_scatter_edges(arrclay[time], extra_clay2)
            for jx in range(len(extra_edge)):
                if extra_edge[jx][0] < extra_edge[jx][1]:
                    per_bond.append([extra_edge[jx][0]+len(arrwater[0]), extra_edge[jx][1]+len(arrwater[0])])
                else:
                    per_bond.append([extra_edge[jx][1] + len(arrwater[0]), extra_edge[jx][0] + len(arrwater[0])])
        print(len(per_bond), per_bond)
        return per_bond

    def get_scatter_edges(self, list1, numlist):
        numlist.sort()
        whe_in_numlist = np.array([ix in numlist for ix in range(len(list1))])
        # print(whe_in_numlist, list1)
        numlist_x = list1[whe_in_numlist][:, 0]
        min_ind = np.argmin(numlist_x)
        anglelist = []
        finedge = []
        res = []
        for indx in range(len(numlist)):
            if indx != min_ind:
                vector1 = np.array([list1[numlist[indx]][0] - list1[numlist[min_ind]][0],
                                    list1[numlist[indx]][1] - list1[numlist[min_ind]][1]])
                vector2 = np.array([0, 1])
                anglelist.append(20000*180/3.1415926*math.acos(vector1.dot(vector2) / np.sqrt(vector1.dot(vector1)))
                                 + np.sqrt(vector1.dot(vector1)))
            else:
                anglelist.append(0)
        angle2 = sorted(anglelist)
        for inde in range(len(angle2)):
            finedge.append(anglelist.index(angle2[inde]))
        for inde in range(len(finedge) - 1):
            res.append([numlist[finedge[inde]], numlist[finedge[inde + 1]]])
        if len(numlist) > 2:
            res.append([numlist[finedge[0]], numlist[finedge[-1]]])
        return res

    def get_watlist(self, guest, time, arr_ele1, arr_ele2, arr_clay, arr_box):
        per_list = []
        cl_list = []
        for jx in range(len(arr_ele1[time])):
            if self.caldis(arr_ele2[time][guest], arr_ele1[time][jx][0], box=arr_box[time], maxleng=6.5):
                per_list.append(jx)
        for jx in range(len(arr_clay[time])):
            if self.caldis(arr_ele2[time][guest], arr_clay[time][jx], box=arr_box[time], maxleng=6):
                cl_list.append(jx)
        return per_list, cl_list

    def whe_2_waters(self, water1, water2, box):
        if self.caldis(water1[0], water2[0], box):
            if self.caldis(water1[1], water2[0], box, maxleng=2.5) and self.calang(water1[0], water2[0], water1[1], box):
                return 1
            elif self.caldis(water1[2], water2[0], box, maxleng=2.5) and self.calang(water1[0], water2[0], water1[2], box):
                return 1
            elif self.caldis(water1[0], water2[1], box, maxleng=2.5) and self.calang(water2[0], water1[0], water2[1], box):
                return 1
            elif self.caldis(water1[0], water2[2], box, maxleng=2.5) and self.calang(water2[0], water1[0], water2[2], box):
                return 1
        return 0

    def whe_ow_os(self, water1, clayo, box):
        if self.caldis(water1[0], clayo, box):
            if self.caldis(water1[1], clayo, box, maxleng=2.5) and self.calang(water1[0], clayo, water1[1], box):
                return 1
            elif self.caldis(water1[2], clayo, box, maxleng=2.5) and self.calang(water1[0], clayo, water1[2], box):
                return 1
        return 0


class Cycle(Basic_Calcul):
    def __init__(self, getlmpdata, hbonds, time):
        self.arrele1 = getlmpdata.arr_water
        self.arrbox = getlmpdata.arr_box
        self.arrclayo = getlmpdata.arr_clayo
        tmptmp = self.get_atom_cycles(time, hbonds)
        self.allcycle = tmptmp[0]
        self.node = tmptmp[2]
        self.whefin = 1
        # self.whefin = tmptmp[3]
        fincycle = tmptmp[1]

        # remove the extra small cycles
        per_remove = []
        for per_cycle in self.allcycle:
            per_remove.append(self.get_left(per_cycle, fincycle))
        for per_cycle in per_remove:
            if len(per_cycle):
                self.allcycle.remove(per_cycle)
        self.cycle_edges = []
        for mx in range(len(self.allcycle)):
            for nx in range(len(self.allcycle[mx])):
                if self.allcycle[mx][nx] not in self.cycle_edges:
                    self.cycle_edges.append(self.allcycle[mx][nx])

        if len(self.allcycle):
            self.maxring = max(max([len(per_cycl) for per_cycl in self.allcycle]), len(fincycle))
            self.fincycle, self.whefin = self.get_min_fin(self.allcycle, time)
            self.clathlen = len(self.allcycle) + 1
        else:
            self.maxring = 0
            self.fincycle = []
            self.clathlen = 0
            self.whefin = 0
        # print([len(per_cycl) for per_cycl in self.allcycle], len(fincycle))
        print(self.cycle_edges)


    @staticmethod
    def get_left(list1, fincycle):
        fin_vec = list(set([jx for ix in fincycle for jx in ix]))
        tmp_vec = list(set([jx for ix in list1 for jx in ix]))
        for eve_node in tmp_vec:
            if eve_node not in fin_vec:
                return []
        share_edge = []
        for per_edge in list1:
            if per_edge in fincycle:
                share_edge.append(per_edge)
            else:
                fincycle.append(per_edge)
        for per_edge in share_edge:
            fincycle.remove(per_edge)
        return list1

    def findmax_link_edge(self, cycle_edge):
        edges = copy.deepcopy(cycle_edge)
        all_sep = []
        eve_cycle = []
        eve_nodes = []
        if len(edges):
            while len(edges):
                kk = self.judge_cycle(edges)
                all_sep.append(len(kk[2]))
                eve_cycle.append(kk[1])
                eve_nodes.append(kk[0])
            max_index = all_sep.index(max(all_sep))
            kk = [eve_nodes[max_index], eve_cycle[max_index]]
            eve_nodes[0] = eve_nodes[max_index]
            allcycles = []
            for ix in range(len(kk[1])):
                per_cycle = self.get_shortest_path(kk[0], kk[1][ix][0], kk[1][ix][1])
                allcycles.append(per_cycle)
        else:
            eve_nodes[0] = []
            allcycles = []
        print("allcycles", allcycles)

    def deep_first_search(self, cur, val, path=[]):
        path.append(cur.value)  # 当前节点值添加路径列表
        if cur.value == val:  # 如果找到目标 返回路径列表
            return path

        if cur.child_list == []:  # 如果没有孩子列表 就 返回 no 回溯标记
            return 'no'

        for node in cur.child_list:  # 对孩子列表里的每个孩子 进行递归
            t_path = copy.deepcopy(path)  # 深拷贝当前路径列表
            res = self.deep_first_search(node, val, t_path)
            if res == 'no':  # 如果返回no，说明找到头 没找到  利用临时路径继续找下一个孩子节点
                continue
            else:
                return res  # 如果返回的不是no 说明 找到了路径

        return 'no'  # 如果所有孩子都没找到 则 回溯

    def get_shortest_path(self, sear, start, end):
        # 分别获取 从根节点 到start 和end 的路径列表，如果没有目标节点 就返回no
        path1 = self.deep_first_search(sear, start, [])
        path2 = self.deep_first_search(sear, end, [])
        if path1 == 'no' or path2 == 'no':
            return '无穷大', '无节点'
        # 对两个路径 从尾巴开始向头 找到最近的公共根节点，合并根节点
        len1, len2 = len(path1), len(path2)
        for i in range(len1 - 1, -1, -1):
            if path1[i] in path2:
                index = path2.index(path1[i])
                path2 = path2[index:]
                path1 = path1[-1:i:-1]
                break
        res = path1 + path2
        # length = len(res)
        # path = '->'.join([str(ix) for ix in res])
        # return '%s:%s'%(length, path)
        if res[0] < res[-1]:
            finres = [[res[0], res[-1]]]
        else:
            finres = [[res[-1], res[0]]]
        for ix in range(len(res) - 1):
            if res[ix] < res[ix + 1]:
                finres.append([res[ix], res[ix + 1]])
            else:
                finres.append([res[ix + 1], res[ix]])
        return finres

    def xor_list(self, list1, list2):
        if len(list1) != len(list2):
            return
        res_list = []
        for ix in range(len(list1)):
            res_list.append(list1[ix] ^ list2[ix])
        return res_list

    def and_list(self, list1, list2):
        if len(list1) != len(list2):
            return
        res_list = []
        for ix in range(len(list1)):
            res_list.append(list1[ix] & list2[ix])
        return self.get_num_edge(res_list)

    def and_list2(self, list1, list2):
        if len(list1) != len(list2):
            return
        res_list = []
        for ix in range(len(list1)):
            res_list.append(list1[ix] & list2[ix])
        return res_list

    def get_bi_list(self, per_cycle, edges):
        res_bi_list = []
        for ix in range(len(edges)):
            if edges[ix] in per_cycle or [edges[ix][1], edges[ix][0]] in per_cycle:
                res_bi_list.append(1)
            else:
                res_bi_list.append(0)
        return res_bi_list

    def get_num_edge(self, per_cycle):
        num = 0
        for ix in range(len(per_cycle)):
            num += per_cycle[ix]
        return num

    def rever_bi2list(self, bi_list, edges):
        res = []
        for ix in range(len(bi_list)):
            if bi_list[ix]:
                res.append(edges[ix])
        return res

    def combine_cycle(self, allcycles, edges, time):
        for ix in range(len(allcycles)):
            for jx in range(ix + 1, len(allcycles)):
                bi_list1 = self.get_bi_list(allcycles[ix], edges)
                bi_list2 = self.get_bi_list(allcycles[jx], edges)
                if self.and_list(bi_list1, bi_list2):
                    tmp_list = self.xor_list(bi_list1, bi_list2)
                    sub_cycle = self.rever_bi2list(tmp_list, edges)
                    # print("sub_cycle", sub_cycle)
                    if self.get_num_edge(tmp_list) <= max(
                            [self.get_num_edge(bi_list1), self.get_num_edge(bi_list2)]):
                        poly_area1 = self.polygon_area(allcycles[ix], time, self.arrele1, self.arrclayo, self.arrbox)
                        poly_area2 = self.polygon_area(allcycles[jx], time, self.arrele1, self.arrclayo, self.arrbox)
                        poly_area3 = self.polygon_area(sub_cycle, time, self.arrele1, self.arrclayo, self.arrbox)
                        if self.get_num_edge(tmp_list) <= self.get_num_edge(bi_list1) and poly_area1 > poly_area3:
                            allcycles[ix] = sub_cycle
                            return 1
                        elif self.get_num_edge(tmp_list) <= self.get_num_edge(bi_list2) and poly_area2 > poly_area3:
                            allcycles[jx] = sub_cycle
                            return 1
        return 0

    def judge_cycle(self, edges):
        """find the cycles from the input edges"""
        vertices = []
        ver_nodes = []
        for kx in [jx for ix in edges for jx in ix]:
            if kx not in vertices:
                vertices.append(kx)
        allnodes = [vertices[0]]
        span_tree = [vertices[0]]
        cycle = []
        not_exam_ver = copy.deepcopy(vertices)
        for ix in range(len(vertices)):
            ver_nodes.append(Node(vertices[ix]))
        root = ver_nodes[0]
        while len(not_exam_ver) + len(span_tree) > len(set(not_exam_ver + span_tree)):
            cur_ver = not_exam_ver[0]
            del_edge_copy = []
            whe_del = 1
            for ix in edges:
                if cur_ver == ix[0]:
                    if ix[0] in span_tree:
                        if ix[1] not in span_tree:
                            ver_nodes[vertices.index(cur_ver)].add_child(ver_nodes[vertices.index(ix[1])])
                            span_tree.append(ix[1])
                            allnodes.append(ix[1])
                        else:
                            cycle.append(ix)
                        del_edge_copy.append(ix)
                    else:
                        not_exam_ver = not_exam_ver[1:] + [not_exam_ver[0]]
                        whe_del = 0
                        break
                elif cur_ver == ix[1]:
                    if ix[1] in span_tree:
                        if ix[0] not in span_tree:
                            ver_nodes[vertices.index(cur_ver)].add_child(ver_nodes[vertices.index(ix[0])])
                            span_tree.append(ix[0])
                            allnodes.append(ix[0])
                        else:
                            cycle.append(ix)
                        del_edge_copy.append(ix)
                    else:
                        not_exam_ver = not_exam_ver[1:] + [not_exam_ver[0]]
                        whe_del = 0
                        break
            if whe_del:
                for ix in del_edge_copy:
                    edges.remove(ix)
                not_exam_ver.remove(not_exam_ver[0])
        return root, cycle, allnodes

    def whe_fincycle(self, cycles, fincycle, edges, time):
        for xx in range(len(cycles)):
            bi_list1 = self.get_bi_list(cycles[xx], edges)
            bi_list2 = self.get_bi_list(fincycle, edges)
            if self.and_list(bi_list1, bi_list2):
                tmp_list = self.xor_list(bi_list1, bi_list2)
                sub_cycle = self.rever_bi2list(tmp_list, edges)
                poly_area1 = self.polygon_area(fincycle, time, self.arrele1, self.arrclayo, self.arrbox)
                poly_area3 = self.polygon_area(sub_cycle, time, self.arrele1, self.arrclayo, self.arrbox)
                if self.get_num_edge(tmp_list) <= self.get_num_edge(bi_list2) and poly_area1 > poly_area3:
                    # sharelist = self.and_list2(bi_list1, bi_list2)
                    # removelist = []
                    # for mm in range(len(sharelist)):
                    #     if sharelist[mm]:
                    #         removelist.append(edges[mm])
                    # for mm in range(len(removelist)):
                    #     edges.remove(removelist[mm])
                    return 1
        return 0

    def get_min_fin(self, allcycles, time):
        alledges = []
        for percycle in allcycles:
            for peredge in percycle:
                if peredge not in alledges and [peredge[1], peredge[0]] not in alledges:
                    alledges.append(peredge)
        while 1:
            while 1:
                allcycle_len = [len(allcycles[ix]) for ix in range(len(allcycles))]
                tmp_bi_cycle = [0 for ix in range(len(alledges))]
                for ix in range(len(allcycles)):
                    tmp_bi_cycle = self.xor_list(tmp_bi_cycle,
                                                 self.get_bi_list(allcycles[ix], alledges))
                fincycle = self.rever_bi2list(tmp_bi_cycle, alledges)
                if self.get_num_edge(tmp_bi_cycle) < max(allcycle_len):
                    allcycles[allcycle_len.index(max(allcycle_len))] = fincycle
                else:
                    break
            res_tmp = self.combine_cycle(allcycles, alledges, time)
            if not res_tmp:
                break
        if len(allcycles) > 1:
            whe_cc = self.whe_fincycle(allcycles, fincycle, alledges, time)
        else:
            whe_cc = 0
        return fincycle, whe_cc

    def get_atom_cycles(self, time, inputedges):
        # nnn = 0
        edges = copy.deepcopy(inputedges)
        all_sep = []
        eve_cycle = []
        eve_nodes = []
        if len(edges):
            while len(edges):
                kk = self.judge_cycle(edges)
                all_sep.append(len(kk[2]))
                eve_cycle.append(kk[1])
                eve_nodes.append(kk[0])
            max_index = all_sep.index(max(all_sep))
            kk = [eve_nodes[max_index], eve_cycle[max_index]]
            eve_nodes[0] = eve_nodes[max_index]
            allcycles = []
            for ix in range(len(kk[1])):
                per_cycle = self.get_shortest_path(kk[0], kk[1][ix][0], kk[1][ix][1])
                allcycles.append(per_cycle)
            if len(allcycles):
                fincycle, whe_cc = self.get_min_fin(allcycles, time)
            else:
                fincycle = []
                # whe_cc = 0
        else:
            allcycles, fincycle = [], []
            eve_nodes.append([])
        return allcycles, fincycle, eve_nodes[0]

    def output2(self, allcycle, time, filename):
        outele = []
        for t1 in range(len(allcycle)):
            for t2 in range(len(allcycle[t1])):
                if allcycle[t1][t2][0] not in outele:
                    outele.append(allcycle[t1][t2][0])
                if allcycle[t1][t2][1] not in outele:
                    outele.append(allcycle[t1][t2][1])
        res_dic = {}
        clath_atoms = []
        numnum = 0
        print("numnum", numnum, len(clath_atoms))
        lenwater = len(self.arrele1[time])
        for ix in range(len(outele)):
            index1 = outele[ix]
            res_dic[index1] = numnum
            if index1 < lenwater:
                clath_atoms.append(
                    ["ow", self.arrele1[time][index1][0][0], self.arrele1[time][index1][0][1], self.arrele1[time][index1][0][2]])
                clath_atoms.append(
                    ["hw", self.arrele1[time][index1][1][0], self.arrele1[time][index1][1][1], self.arrele1[time][index1][1][2]])
                clath_atoms.append(
                    ["hw", self.arrele1[time][index1][2][0], self.arrele1[time][index1][2][1], self.arrele1[time][index1][2][2]])
                numnum += 3
            else:
                index2 = index1 - lenwater
                clath_atoms.append(["ob", self.arrclayo[time][index2][0], self.arrclayo[time][index2][1], self.arrclayo[time][index2][2]])
                numnum += 1
        fileout = open(filename, "w")
        print(len(clath_atoms), file=fileout)
        print(len(clath_atoms), file=fileout)
        for ix in range(len(clath_atoms)):
            for jx in range(3):
                clath_atoms[ix][jx + 1] = clath_atoms[0][jx + 1] + self.get_vector(clath_atoms[0][jx + 1], clath_atoms[ix][jx + 1],
                                                                            self.arrbox[time][jx])
            print(clath_atoms[ix][0], clath_atoms[ix][1], clath_atoms[ix][2], clath_atoms[ix][3], file=fileout)
        fileout.close()
        return res_dic


class Output(Basic_Calcul):
    def __init__(self, getlmpdata, time, guestnum, node, cycle_edges, molid, fileout1, fileout2):
        self.wheclay = 0
        self.arrwater = getlmpdata.arr_water[time]
        self.arrclay = getlmpdata.arr_clayo[time]
        self.arrbox = getlmpdata.arr_box[time]
        self.center = getlmpdata.arr_methane[time][guestnum]
        self.resdic = self.output(node, fileout1)
        self.outputlabel(self.resdic, cycle_edges, molid, fileout2)

    def get_allnode(self, Node, allnodes=[]):
        allnodes.append(Node.value)
        if Node.child_list == []:
            return 0
        else:
            for ix in Node.child_list:
                self.get_allnode(ix, allnodes)
        return allnodes

    def output(self, node, fileout):
        allnodes = self.get_allnode(node, allnodes=[])
        res_dic = {}
        clath_atoms = [["CT", self.center[0], self.center[1], self.center[2]]]
        numnum = 1
        lenwater = len(self.arrwater)
        for ix in range(len(allnodes)):
            index1 = allnodes[ix]
            res_dic[index1] = numnum
            if index1 < lenwater:
                clath_atoms.append(["ow", self.arrwater[index1][0][0], self.arrwater[index1][0][1], self.arrwater[index1][0][2]])
                clath_atoms.append(["hw", self.arrwater[index1][1][0], self.arrwater[index1][1][1], self.arrwater[index1][1][2]])
                clath_atoms.append(["hw", self.arrwater[index1][2][0], self.arrwater[index1][2][1], self.arrwater[index1][2][2]])
                numnum += 3
            else:
                self.wheclay = 1
                index2 = index1 - lenwater
                clath_atoms.append(["ob", self.arrclay[index2][0], self.arrclay[index2][1], self.arrclay[index2][2]])
                numnum += 1
        print("whe_clay:", self.wheclay, file=fileout)
        print(len(clath_atoms), file=fileout)
        print(len(clath_atoms), file=fileout)
        for ix in range(len(clath_atoms)):
            for jx in range(3):
                clath_atoms[ix][jx+1] = self.center[jx] + self.get_vector(self.center[jx], clath_atoms[ix][jx+1], self.arrbox[jx])
            print(clath_atoms[ix][0], clath_atoms[ix][1], clath_atoms[ix][2], clath_atoms[ix][3], file=fileout)
        # fileout.close()
        # print(res_dic)
        return res_dic

    def outputlabel(self, resdic, cycle_edges, molid, fileout):
        print("wheclay", self.wheclay, file=fileout)
        for eachedge in cycle_edges:
            print("label add Bonds ", str(molid) + '/' + str(resdic[eachedge[0]]),
                  str(molid) + '/' + str(resdic[eachedge[1]]), file=fileout)


traj1 = GetLmpData("C:\\Users\\lq\\Desktop\\real1dump")
totaltime, totalnum, = traj1.arr_methane.shape[:2]
# totaltime, totalnum = 1, 1
fileout1 = open("xyxyz", "w")
fileout2 = open("label", "w")
for now_time in range(totaltime):
    for now_mol in range(383, 384):
# for now_time in range(0, 1):
#     for now_mol in range(0, 5):
        hbond1 = GetHbonds(traj1, guestnum=now_mol, time=now_time)
        hbedges = hbond1.Hblist
        # print("now_time, now_nol", now_time, now_mol)
        cycle1 = Cycle(traj1, hbedges, now_time)
        # print("maxring:", cycle1.maxring, "clathlen:", cycle1.clathlen)

        # print(cycle1.allcycle)
        if cycle1.maxring < 12 and cycle1.clathlen > 6:
            print("new molecule: time:", now_time, "mol:", now_mol, file=fileout1)
            print("maxring:", cycle1.maxring, "clathlen:", cycle1.clathlen,
                  "whe_fin:", cycle1.whefin, "molid: time:", now_time, "mol:", now_mol, file=fileout2, end=' ')
            Output(traj1, now_time, now_mol, cycle1.node, cycle1.cycle_edges, 0, fileout1, fileout2)
fileout1.close()
fileout2.close()
