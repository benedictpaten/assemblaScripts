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

static stHash *setOfPairs;

static void getMAFBlock2(Block *block, FILE *fileHandle) {
    /*
     * Prints out the comment lines, then the maf blocks.
     */
    if (block_getLength(block) >= minimumBlockLength) {
        Segment *segment;
        Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
        int32_t assemblyNumber = 0, hapA1Number = 0, hapA2Number = 0;
        while ((segment = block_getNext(instanceIt)) != NULL) {
            const char *segmentEvent = event_getHeader(segment_getEvent(segment));
            if (strcmp(segmentEvent, assemblyEventString) == 0) { //Establish if we need a line..
                assemblyNumber++;
            } else if (strcmp(segmentEvent, hap1EventString) == 0) {
                hapA1Number++;
            } else if (strcmp(segmentEvent, hap2EventString) == 0) {
                hapA2Number++;
            } else {
                assert(strcmp(segmentEvent, contaminationEventString) == 0);
            }
        }
        block_destructInstanceIterator(instanceIt);
        if (assemblyNumber > 0 || hapA1Number > 0 || hapA2Number > 0) {
            stIntTuple *key = stIntTuple_construct(3, hapA1Number > hapA2Number ? hapA1Number : hapA2Number,
                    hapA1Number < hapA2Number ? hapA1Number : hapA2Number, assemblyNumber);
            int32_t *value;
            if ((value = stHash_search(setOfPairs, key)) == NULL) {
                value = st_malloc(sizeof(int32_t));
                value[0] = 0;
                stHash_insert(setOfPairs, key, value);
            } else {
                stIntTuple_destruct(key);
            }
            value[0] += block_getLength(block);
        }
    }
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "copyNumberStats");

    ///////////////////////////////////////////////////////////////////////////
    // Now use the MAF printing code to generate the results..
    ///////////////////////////////////////////////////////////////////////////

    int64_t startTime = time(NULL);
    FILE *fileHandle = fopen(outputFile, "w");

    //The pairs to represent the mafs.
    setOfPairs = stHash_construct3((uint32_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn, (void(*)(void *)) stIntTuple_destruct, free);
    //Pass over the blocks.
    getMAFsReferenceOrdered(flower, fileHandle, getMAFBlock2);
    //Now calculate the linkage stats
    stList *copyNumbers = stHash_getKeys(setOfPairs);
    stList_sort(copyNumbers, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
    int32_t totalCopyNumberDeficientColumns = 0;
    int32_t totalCopyNumberDeficientBases = 0;
    int32_t totalCopyNumberDeficientColumnsGreaterThanZero = 0;
    int32_t totalCopyNumberDeficientBasesGreaterThanZero = 0;
    int32_t totalCopyNumberExcessColumns = 0;
    int32_t totalCopyNumberExcessBases = 0;
    int32_t totalColumnCount = 0;
    int32_t totalBaseCount = 0;
    for (int32_t i = 0; i < stList_length(copyNumbers); i++) {
        stIntTuple *copyNumber = stList_get(copyNumbers, i);
        int32_t *columnCount = stHash_search(setOfPairs, copyNumber);
        totalColumnCount += columnCount[0];
        totalBaseCount += stIntTuple_getPosition(copyNumber, 2) * columnCount[0];
    }
    fprintf(fileHandle, "<copy_number_stats minimumBlockLength=\"%i\" totalColumnCount=\"%i\" totalBaseCount=\"%i\">\n", minimumBlockLength, totalColumnCount, totalBaseCount);
    for (int32_t i = 0; i < stList_length(copyNumbers); i++) {
        stIntTuple *copyNumber = stList_get(copyNumbers, i);
        int32_t *columnCount = stHash_search(setOfPairs, copyNumber);
        int32_t maxHapNumber = stIntTuple_getPosition(copyNumber, 0);
        int32_t minHapNumber = stIntTuple_getPosition(copyNumber, 1);
        assert(minHapNumber >= 0);
        assert(maxHapNumber >= minHapNumber);
        int32_t assemblyNumber = stIntTuple_getPosition(copyNumber, 2);
        assert(assemblyNumber >= 0);
        assert(columnCount != NULL);
        assert(columnCount[0] >= 1);
        fprintf(
                fileHandle,
                "<copy_number_category maximumHaplotypeCopyNumber=\"%i\" minimumHaplotypeCopyNumber=\"%i\" assemblyCopyNumber=\"%i\" columnCount=\"%i\"/>\n",
                maxHapNumber, minHapNumber, assemblyNumber, columnCount[0]);
        if (assemblyNumber < minHapNumber) {
            totalCopyNumberDeficientColumns += columnCount[0];
            totalCopyNumberDeficientBases += columnCount[0] * (minHapNumber - assemblyNumber);
            if(assemblyNumber > 0) {
                totalCopyNumberDeficientColumnsGreaterThanZero += columnCount[0];
                totalCopyNumberDeficientBasesGreaterThanZero += columnCount[0] * (minHapNumber - assemblyNumber);
            }
         } else if (assemblyNumber > maxHapNumber) {
            totalCopyNumberExcessColumns += columnCount[0];
            totalCopyNumberExcessBases += columnCount[0] * (assemblyNumber - maxHapNumber);
        }
    }
    fprintf(fileHandle, "<deficientCopyNumberCounts totalColumns=\"%i\" totalBases=\"%i\" totalProportionOfColumns=\"%f\" totalProportionOfBases=\"%f\"/>",
            totalCopyNumberDeficientColumns, totalCopyNumberDeficientBases,
            ((float)totalCopyNumberDeficientColumns)/totalColumnCount, ((float)totalCopyNumberDeficientBases)/totalBaseCount);
    fprintf(fileHandle, "<deficientCopyNumberCountsGreaterThanZero totalColumns=\"%i\" totalBases=\"%i\" totalProportionOfColumns=\"%f\" totalProportionOfBases=\"%f\"/>",
                totalCopyNumberDeficientColumnsGreaterThanZero, totalCopyNumberDeficientBasesGreaterThanZero,
                ((float)totalCopyNumberDeficientColumnsGreaterThanZero)/totalColumnCount, ((float)totalCopyNumberDeficientBasesGreaterThanZero)/totalBaseCount);
    fprintf(fileHandle, "<excessCopyNumberCounts totalColumns=\"%i\" totalBases=\"%i\" totalProportionOfColumns=\"%f\" totalProportionOfBases=\"%f\"/>",
            totalCopyNumberExcessColumns, totalCopyNumberExcessBases,
            ((float)totalCopyNumberExcessColumns)/totalColumnCount, ((float)totalCopyNumberExcessBases)/totalBaseCount);
    fprintf(fileHandle, "</copy_number_stats>\n");
    fclose(fileHandle);

    st_logInfo("Got the copy number counts in %i seconds/\n", time(NULL) - startTime);

    return 0;
}
