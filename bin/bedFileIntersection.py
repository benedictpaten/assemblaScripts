import sys
import bisect
import xml.etree.ElementTree as ET
import random

def parseBedFile(file):
    fn = lambda (seqName, start, end) : (seqName, int(start), int(end))
    intervals = [ fn(line.split()[:3]) for line in open(file, 'r').readlines() ]
    intervals.sort()
    seqName, start, end = intervals[0]
    nonOverlappingIntervals = intervals[:1]
    for seqName2, start2, end2 in intervals[1:]:
        if seqName == seqName2 and start2 >= start and end2 <= end:
            continue
        seqName, start, end = seqName2, start2, end2 
        nonOverlappingIntervals.append((seqName, start, end))
    return nonOverlappingIntervals

def getContainers2(seqName, start, end, pathIntervals):
    return [ (seqName2, start2, end2) for seqName2, start2, end2 in pathIntervals if seqName2 == seqName and start >= start2 and end <= end2 ]

def getContainers(seqName, start, end, pathIntervals):
    i = bisect.bisect_left(pathIntervals, (seqName, start, -1))
    while i > 0:
        if i < len(pathIntervals):
            seqName2, start2, end2 = pathIntervals[i]
            if seqName2 < seqName:
                break
            if seqName == seqName2 and end2 < start:
                break
        i -= 1
    j = i
    containers = []
    for i in xrange(i, len(pathIntervals)):
        seqName2, start2, end2 = pathIntervals[i]
        if seqName < seqName2 or (seqName == seqName2 and start < start2):
            if random.random() > 0.999 and containers != getContainers2(seqName, start, end, pathIntervals):
                print "path intervals", pathIntervals[j:i+1]
                print "coordinates", start, end, seqName
                print "containers1", containers
                print "containers2", getContainers2(seqName, start, end, pathIntervals)
                assert 0
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
pathIntervalsLength = sum( [ (end - start + 1) for seqName, start, end in pathIntervals ])
stats = ET.Element("stats", attrib={ "msaFile":sys.argv[1], "intervalsLength":str(pathIntervalsLength)})
for bedFile in sys.argv[3:]:
    samples, complete, totalLength, totalContainment = getContainment(pathIntervals, parseBedFile(bedFile))
    tag = ET.SubElement(stats, "intervals", attrib={ "featureFile":bedFile, "complete":str(len(complete)), "samples":str(samples), "baseLength":str(totalLength), "totalComplete":str(totalContainment) })
    tag.text = " ".join(complete)

fileHandle = open(sys.argv[2], "w")
ET.ElementTree(stats).write(fileHandle)
fileHandle.close()