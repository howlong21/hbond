import sys, os
import math
import matplotlib.pyplot as plt


class Hierarchical:
    def __init__(self, center, left=None, right = None, flag = None, distance = 0.0):
        self.center = center
        self.left = left
        self.right = right
        self.flag = flag
        self.distance = distance


def traverse(node):
    if node.left == None and node.right == None:
        return [node.center]
    else:
        return traverse(node.left) + traverse(node.right)


def distance(v1, v2, box):
    if len(v1) != len(v2):
        print sys.stderr, "invalid v1 and v2 !"
        sys.exit(1)
    distance = 0
    for i in range(len(v1)):
        distance += (box[i]/2 - abs(box[i]/2 - abs(v1[i] - v2[i]))) ** 2
    distance = math.sqrt(distance)
    return distance


def hcluster(data, n, length, box):
    if len(data) <= 0:
        print sys.stderr, "invalid data"
        sys.exit(1)


    clusters = [Hierarchical(data[i], flag = i) for i in range(len(data))]
    print(clusters[1].center)
    distances = {}
    min_id1 = None
    min_id2 = None
    currentCluster = -100


    while(len(clusters) > n):
        minDist = 1000000000000


        for i in range(len(clusters) - 1):
            for j in range(i + 1, len(clusters)):
                # save distance, pick up speed
                if distances.get((clusters[i].flag, clusters[j].flag)) == None:
                    alldis = []
                    for tmp1 in traverse(clusters[i]):
                        for tmp2 in traverse(clusters[j]):
                            alldis.append(distance(tmp1, tmp2, box))
                    distances[(clusters[i].flag, clusters[j].flag)] = min(alldis)


                if distances[(clusters[i].flag, clusters[j].flag)] <= minDist:
                    min_id1 = i
                    min_id2 = j
                    minDist = distances[(clusters[i].flag, clusters[j].flag)]


        if min_id1 != None and min_id2 != None and minDist != 1000000000000:
            newCenter = [(clusters[min_id1].center[i] + clusters[min_id2].center[i])/2 for i in range(len(clusters[min_id2].center))]
            newFlag = currentCluster
            currentCluster -= 1
            newCluster = Hierarchical(newCenter, clusters[min_id1], clusters[min_id2], newFlag, minDist)
            del clusters[min_id2]
            del clusters[min_id1]
            clusters.append(newCluster)
        length.append(minDist)
    finalCluster = [traverse(clusters[i]) for i in range(len(clusters))]
    return finalCluster


def loadData(filename):
    infile = open(filename, 'r')
    line = infile.readline()
    dataList = []
    tempList = []
    while line:
        lineArr = line.strip().split()
        if len(lineArr) < 2:
            line = infile.readline()
            continue
        for i in range(len(lineArr)):
            tempList.append(float(lineArr[i]))
        dataList.append(tempList)
        tempList = []
        line = infile.readline()
    return dataList


if __name__ == '__main__':
    # data = [[0, 0], [1, 0], [2.1, 0], [3.3, 0], [10, 0], [50, 0], [90, 0]]
    data = loadData('ccxyz')
    length = []
    finalCluster = hcluster(data, 140, length, box=[51.943, 130.895, 45.9719])
    print finalCluster
    print(length)
    finx = list(range(len(length)))
    print(len(finalCluster))
    num = 0
    index = 0
    for ix in range(len(finalCluster)):
        if len(finalCluster[ix]) > num:
            num = len(finalCluster[ix])
            index = ix
    print(num)
    filecc = open('caxyz', 'w')
    print >> filecc, "%d" % (len(finalCluster[index]))
    print >> filecc, "%d" % (len(finalCluster[index]))
    for ix in finalCluster[index]:
        print >> filecc, "%s %.4f %.4f %.4f" % ("co", ix[0], ix[1], ix[2])
    filecc.close()
    plt.plot(finx, length)
    plt.show()
