import sys
import bisect
import xml.etree.ElementTree as ET
from bedFileIntersection import parseBedFile, getContainer

def parseGeneBedFile(file):
    fn = lambda (seqName, start, end, gene) : (seqName, int(start), int(end), gene)
    return [ fn(line.split()[:4]) for line in open(file, 'r').readlines() ]

def getContainment(pathIntervals, featureIntervals):
    seqName, start, end = featureIntervals[0]
    i = getContainer(seqName, start, end)
    if i != None:
        seqName2, start2, end2 = i
        for seqName, start, end, in featureIntervals:
            if seqName != seqName2 or int(start) < int(start2) or int(end) > int(end2):
                break
        else:
            return 1
    return 0

pathIntervals = parseBedFile(sys.argv[1])
stats = ET.Element("stats", attrib={ "msaFile":sys.argv[1]})
for bedFile in sys.argv[3:]:
    samples = 0
    complete = 0
    genesToIntervals = {}
    for seqName, start, end, gene in parseGeneBedFile(bedFile):
        if not genesToIntervals.has_key(gene):
             genesToIntervals[gene] = []
        genesToIntervals[gene].append((seqName, start, end))
    for gene in genesToIntervals.keys():
        samples += 1
        if getContainment(pathIntervals, genesToIntervals[gene]):
            complete += 1
    ET.SubElement(stats, "intervals", attrib={ "featureFile":bedFile, "complete":str(complete), "samples":str(samples) })

fileHandle = open(sys.argv[2], "w")
ET.ElementTree(stats).write(fileHandle)
fileHandle.close()