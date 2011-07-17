import sys
import bisect
import xml.etree.ElementTree as ET

def parseBedFile(file):
    fn = lambda (seqName, start, end) : (seqName, int(start), int(end))
    return [ fn(line.split()[:3]) for line in open(file, 'r').readlines() ]

def getContainers(seqName, start, end, pathIntervals):
    i = bisect.bisect_left(pathIntervals, (seqName, start, -1))
    if i > 0:
        i -= 1
    containers = []
    for i in xrange(i, len(pathIntervals)):
        seqName2, start2, end2 = pathIntervals[i]
        if seqName < seqName2 or (seqName == seqName2 and start < start2):
            return containers
        if seqName != seqName2:
            continue
        if end <= end2:
            assert start >= start2
            assert seqName == seqName2
            containers.append((seqName2, start2, end2))
    return containers
        
def getContainment(pathIntervals, featureIntervals):
    samples = 0
    complete = []
    totalLength = 0
    totalContainment = 0
    for seqName, start, end in featureIntervals:
        samples += 1
        containers = getContainers(seqName, start, end, pathIntervals)
        length = abs(end - start + 1)
        totalLength += length
        if len(containers) > 0:
            totalContainment += length
            complete.append("_".join([ str(j) for j in containers + [ seqName,start,end ] ]))
    return samples, complete, totalLength, totalContainment

pathIntervals = parseBedFile(sys.argv[1])
pathIntervals.sort()
pathIntervalsLength = sum( [ (end - start + 1) for seqName, start, end in pathIntervals ])
stats = ET.Element("stats", attrib={ "msaFile":sys.argv[1], "intervalsLength":str(pathIntervalsLength)})
for bedFile in sys.argv[3:]:
    samples, complete, totalLength, totalContainment = getContainment(pathIntervals, parseBedFile(bedFile))
    tag = ET.SubElement(stats, "intervals", attrib={ "featureFile":bedFile, "complete":str(len(complete)), "samples":str(samples), "baseLength":str(totalLength), "totalComplete":str(totalContainment) })
    tag.text = " ".join(complete)

fileHandle = open(sys.argv[2], "w")
ET.ElementTree(stats).write(fileHandle)
fileHandle.close()