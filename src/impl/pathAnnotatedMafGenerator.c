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

stHash *segmentsToMaximalHaplotypePaths;
stHash *maximalHaplotypePathLengths;
stHash *maximalScaffoldPathsLengths;
stList *haplotypeEventStrings;
stList *contaminationEventStrings;

static int64_t getNumberOnPositiveStrand(Block *block) {
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    Segment *segment;
    int64_t i = 0;
    while ((segment = block_getNext(it)) != NULL) {
        if (segment_getChildNumber(segment) == 0) {
            if (segment_getStrand(segment)) {
                i++;
            }
        }
    }
    block_destructInstanceIterator(it);
    return i;
}

int getSimpleCode(enum CapCode code) {
    switch (code) {
        case HAP_SWITCH:
        case HAP_NOTHING:
            return 1;
        case CONTIG_END:
        case CONTIG_END_WITH_SCAFFOLD_GAP:
        case CONTIG_END_WITH_AMBIGUITY_GAP:
            return 0;
        case SCAFFOLD_GAP:
        case AMBIGUITY_GAP:
            return 3;
        case ERROR_HAP_TO_HAP_SAME_CHROMOSOME:
        case ERROR_HAP_TO_HAP_DIFFERENT_CHROMOSOMES:
        case ERROR_HAP_TO_CONTAMINATION:
        case ERROR_HAP_TO_INSERT_TO_CONTAMINATION:
        case ERROR_HAP_TO_INSERT:
        case ERROR_HAP_TO_INSERT_AND_DELETION:
        case ERROR_HAP_TO_DELETION:
        case ERROR_CONTIG_END_WITH_INSERT:
            return 2;
        default:
            assert(0);
    }
    assert(0);
    return 2;
}

void getMAFBlock2(Block *block, FILE *fileHandle) {
    /*
     * Prints out the comment lines, then the maf blocks.
     */
    if (getNumberOnPositiveStrand(block) == 0) {
        block = block_getReverse(block);
    }
    Segment *segment;
    Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
    while ((segment = block_getNext(instanceIt)) != NULL) {
        if (strcmp(event_getHeader(segment_getEvent(segment)), assemblyEventString) == 0) { //Establish if we need a line..
            stList *maximalHaplotypePath = stHash_search(
                    segmentsToMaximalHaplotypePaths, segment);
            if (maximalHaplotypePath == NULL) {
                maximalHaplotypePath = stHash_search(
                        segmentsToMaximalHaplotypePaths, segment_getReverse(
                                segment));
            }
            if (maximalHaplotypePath != NULL) {
                assert(stHash_search(maximalHaplotypePathLengths,
                        maximalHaplotypePath) != NULL);
                int64_t length = stIntTuple_getPosition(stHash_search(
                        maximalHaplotypePathLengths, maximalHaplotypePath), 0);
                assert(stHash_search(maximalScaffoldPathsLengths,
                        maximalHaplotypePath) != NULL);
                int64_t scaffoldPathLength = stIntTuple_getPosition(
                        stHash_search(maximalScaffoldPathsLengths,
                                maximalHaplotypePath), 0);
                assert(scaffoldPathLength >= length);
                /*
                 * Codes:
                 * 0 = contig ends.
                 * 1 = correct adjacency.
                 * 2 = error adjacency.
                 * 3 = scaffold gap
                 */

                int64_t insertLength;
                int64_t deleteLength;
                Cap *otherCap;
                enum CapCode _5EndStatusNerd = getCapCode(segment_get5Cap(
                        segment), &otherCap, haplotypeEventStrings, contaminationEventStrings, &insertLength, &deleteLength, capCodeParameters);
                enum CapCode _3EndStatusNerd = getCapCode(segment_get3Cap(
                        segment), &otherCap, haplotypeEventStrings, contaminationEventStrings, &insertLength, &deleteLength, capCodeParameters);

                int64_t _5EndStatus = getSimpleCode(_5EndStatusNerd);
                int64_t _3EndStatus = getSimpleCode(_3EndStatusNerd);

                fprintf(
                        fileHandle,
                        "#HPL=%" PRIi64 " 5=%" PRIi64 " 3=%" PRIi64 " SPL=%" PRIi64 " 5NERD=%" PRIi64 " 3NERD=%" PRIi64 "\n",
                        length, _5EndStatus, _3EndStatus,
                        scaffoldPathLength,
                        (int64_t)_5EndStatusNerd, (int64_t)_3EndStatusNerd);
            }
        }
    }
    block_destructInstanceIterator(instanceIt);

    getMAFBlock(block, fileHandle);
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "pathAnnotatedMafGenerator");

    //////////////////////////////////////////////
    //Get the maximal haplotype path info.
    //////////////////////////////////////////////

    haplotypeEventStrings = getEventStrings(hap1EventString, hap2EventString);
    contaminationEventStrings = getEventStrings(contaminationEventString, NULL);
    stList *maximalHaplotypePaths = getContigPaths(flower, assemblyEventString, haplotypeEventStrings);
    segmentsToMaximalHaplotypePaths = buildSegmentToContigPathHash(
            maximalHaplotypePaths);
    maximalHaplotypePathLengths
            = buildContigPathToContigPathLengthHash(
                    maximalHaplotypePaths);
    maximalScaffoldPathsLengths = getContigPathToScaffoldPathLengthsHash(
            maximalHaplotypePaths, haplotypeEventStrings, contaminationEventStrings, capCodeParameters);

    ///////////////////////////////////////////////////////////////////////////
    // Now print the MAFs
    ///////////////////////////////////////////////////////////////////////////

    int64_t startTime = time(NULL);
    FILE *fileHandle = fopen(outputFile, "w");
    makeMAFHeader(flower, fileHandle);
    getMAFsReferenceOrdered(flower, fileHandle, getMAFBlock2);

    fclose(fileHandle);
    st_logInfo("Got the mafs in %" PRIi64 " seconds/\n", time(NULL) - startTime);

    return 0;
}
