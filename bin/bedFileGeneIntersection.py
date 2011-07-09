import sys
import xml.etree.ElementTree as ET

def parseBedFile(file, fields=3):
    return [ line.split()[:fields] for line in open(file, 'r').readlines() ]

def getContainment(pathIntervals, featureIntervals):
    for seqName2, start2, end2 in pathIntervals:
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
    for seqName, start, end, gene in parseBedFile(bedFile, fields=4):
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