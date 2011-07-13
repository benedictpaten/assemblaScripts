import sys
import bisect
import xml.etree.ElementTree as ET

def parseBedFile(file):
    fn = lambda (seqName, start, end) : (seqName, int(start), int(end))
    return [ fn(line.split()[:3]) for line in open(file, 'r').readlines() ]

def getContainer(seqName, start, end, pathIntervals):
    i = bisect.bisect_left(pathIntervals, (seqName, start, -1))
    for i in xrange(i-1, len(pathIntervals)):
        seqName2, start2, end2 = pathIntervals[i]
        if seqName < seqName2 or start < start2:
            return None
        if seqName != seqName2:
            continue
        if end <= end2:
            return (seqName2, start2, end2)
    return None
        
def getContainment(pathIntervals, featureIntervals):
    samples = 0
    complete = 0
    for seqName, start, end in featureIntervals:
        samples += 1
        if getContainer(seqName, start, end, pathIntervals) != None:
            complete += 1
    return samples, complete

pathIntervals = parseBedFile(sys.argv[1])
pathIntervals.sort()
stats = ET.Element("stats", attrib={ "msaFile":sys.argv[1]})
for bedFile in sys.argv[3:]:
    samples, complete = getContainment(pathIntervals, parseBedFile(bedFile))
    ET.SubElement(stats, "intervals", attrib={ "featureFile":bedFile, "complete":str(complete), "samples":str(samples) })

fileHandle = open(sys.argv[2], "w")
ET.ElementTree(stats).write(fileHandle)
fileHandle.close()