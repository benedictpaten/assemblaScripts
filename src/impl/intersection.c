/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactus.h"
#include "cactusMafs.h"
#include "contigPaths.h"
#include "adjacencyTraversal.h"
#include "adjacencyClassification.h"
#include "scaffoldPaths.h"
#include "assemblaCommon.h"

typedef struct _sequenceInterval {
        int32_t start;
        int32_t length;
        const char *sequenceName;
} SequenceInterval;

SequenceInterval *sequenceInterval_construct(int32_t start, int32_t length,
        const char *sequenceName) {
    SequenceInterval *sequenceInterval = st_malloc(sizeof(SequenceInterval));
    sequenceInterval->start = start;
    sequenceInterval->length = length;
    sequenceInterval->sequenceName = sequenceName;
    return sequenceInterval;
}

void sequenceInterval_destruct(SequenceInterval *sequenceInterval) {
    free(sequenceInterval);
}

stList *getContigPathIntervals(Flower *flower) {
    stList *contigPaths = getContigPaths(flower, sampleEventString, referenceEventStrings);
    stList *intervals = stList_construct3(0, (void (*)(void *))sequenceInterval_destruct);
    for (int32_t i = 0; i < stList_length(contigPaths); i++) {
        stList *contigPath = stList_get(contigPaths, i);
        Segment *_5Segment = stList_get(contigPath, 0);
        Segment *_3Segment = stList_get(contigPath, stList_length(contigPath)-1);
        Sequence *sequence = segment_getSequence(_5Segment);

        assert(sequence != NULL);
        assert(segment_getSequence(_5Segment) == segment_getSequence(_3Segment));
        assert(segment_getStrand(_5Segment));
        assert(segment_getStrand(_3Segment));
        assert(segment_getStart(_5Segment) < segment_getStart(_3Segment));

        stList_append(intervals, sequenceInterval_construct(segment_getStart(_5Segment),
                segment_getStart(_3Segment) + segment_getLength(_5Segment),
                sequence_getHeader(sequence)));
    }
    stList_destruct(contigPaths);
    return intervals;
}

stList *getBedIntervals(FILE *fileHandle) {
    stList *intervals = stList_construct3(0, (void (*)(void *))sequenceInterval_destruct);
    int32_t bufferSize = 100;
    char *line = st_malloc(bufferSize);
    int32_t lineWidth;
    while((lineWidth = benLine(&line, &bufferSize, fileHandle)) != NULL) {
        int32_t start, length;
        char *sequenceName[100];
        sscanf(line, "%s %i %i", sequenceName, &start, &length);
        stList_append(sequenceInterval_construct(start, length, sequenceName));
    }
    free(line);
    return intervals;
}


stList *getIntervalIntersection(stList *interval)

