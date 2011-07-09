import sys
import xml.etree.ElementTree as ET

def parseBedFile(file, fields=3):
    return [ line.split()[:3] for line in open(file, 'r').readlines() ]

def getContainment(pathIntervals, featureIntervals):
    samples = 0
    complete = 0
    for seqName, start, end in featureIntervals:
        samples += 1
        for seqName2, start2, end2 in pathIntervals:
            if seqName == seqName2:
                if int(start) >= int(start2) and int(end) <= int(end2):
                    complete += 1
                    break
    return samples, complete

pathIntervals = parseBedFile(sys.argv[1])
stats = ET.Element("stats", attrib={ "msaFile":sys.argv[1]})
for bedFile in sys.argv[3:]:
    samples, complete = getContainment(pathIntervals, parseBedFile(bedFile))
    ET.SubElement(stats, "intervals", attrib={ "complete":str(complete), "samples":str(samples) })

fileHandle = open(sys.argv[2], "w")
ET.ElementTree(stats).write(fileHandle)
fileHandle.close()