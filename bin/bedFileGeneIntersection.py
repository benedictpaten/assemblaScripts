import sys
import bisect
import xml.etree.ElementTree as ET
from bedFileIntersection import parseBedFile, getContainers

def parseGeneBedFile(file):
    fn = lambda (seqName, start, end, gene) : (seqName, int(start), int(end), gene)
    return [ fn(line.split()[:4]) for line in open(file, 'r').readlines() ]

def getContainment(pathIntervals, featureIntervals):
    seqName, start, end = featureIntervals[0]
    for seqName2, start2, end2 in getContainers(seqName, start, end, pathIntervals):
        print "Found interval %s %i %i %s %i %i" % (seqName, start, end, seqName2, start2, end2)
        for seqName, start, end in featureIntervals:
            print "alright %s %i %i %s %i %i" % (seqName, start, end, seqName2, start2, end2)
            if seqName != seqName2 or int(start) < int(start2) or int(end) > int(end2):
                break
        else:
            return (seqName2, start2, end2)
    return None

pathIntervals = parseBedFile(sys.argv[1])
pathIntervals.sort()
stats = ET.Element("stats", attrib={ "msaFile":sys.argv[1]})
for bedFile in sys.argv[3:]:
    samples = 0
    complete = 0
    genesToIntervals = {}
    totalLength = 0
    totalContainment = 0
    text = []
    for seqName, start, end, gene in parseGeneBedFile(bedFile):
        if not genesToIntervals.has_key(gene):
             genesToIntervals[gene] = []
        genesToIntervals[gene].append((seqName, start, end))
    for gene in genesToIntervals.keys():
        samples += 1
        interval = getContainment(pathIntervals, genesToIntervals[gene])
        length = sum([ abs(end - start + 1) for seqName, start, end in genesToIntervals[gene] ])
        totalLength += length
        if interval != None:
            complete += 1
            totalContainment += length
            text.append("_".join([ gene + "_" + "/".join([ str(j) for j in i ]) for i in genesToIntervals[gene] ]) + "_" + "/".join([ str(k) for k in interval ]))
    tag = ET.SubElement(stats, "intervals", attrib={ "featureFile":bedFile, "complete":str(complete), "samples":str(samples), "baseLength":str(totalLength), "totalComplete":str(totalContainment) })
    tag.text = " ".join(text)

fileHandle = open(sys.argv[2], "w")
ET.ElementTree(stats).write(fileHandle)
fileHandle.close()