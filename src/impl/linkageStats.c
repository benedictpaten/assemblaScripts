/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <ctype.h>
#include <math.h>

#include "sonLib.h"
#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"
#include "treeStats.h"
#include "cactusMafs.h"
#include "adjacencyTraversal.h"
#include "assemblaCommon.h"
#include "linkage.h"

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "linkageStats");

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    double bucketSize = bucketNumber / log10(upperLinkageBound);
    int32_t *correct = st_malloc(sizeof(int32_t) * bucketNumber);
    int32_t *samples = st_malloc(sizeof(int32_t) * bucketNumber);
    for (int32_t i = 0; i < bucketNumber; i++) {
        correct[i] = 0.0;
        samples[i] = 0.0;
    }

    stList *eventStrings = getEventStrings(hap1EventString, hap2EventString);
    stSortedSet *sequences = getHaplotypeSequences(flower, eventStrings);
    stSortedSetIterator *it = stSortedSet_getIterator(sequences);
    MetaSequence *metaSequence;
    while ((metaSequence = stSortedSet_getNext(it)) != NULL) {
        samplePoints(flower, metaSequence, assemblyEventString,
                sampleNumber, correct, samples,
                bucketNumber, bucketSize);
    }
    stSortedSet_destructIterator(it);
    stSortedSet_destruct(sequences);

    ///////////////////////////////////////////////////////////////////////////
    // Write it out.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(fileHandle, "<linkage_stats>\n");
    int32_t pMaxSize = 1;
    int32_t cumulativeCorrect = 0;
    int32_t cumulativeSamples = 0;
    for (int32_t i = 0; i < bucketNumber; i++) {
        if(samples[i] > 0) {
            assert(correct[i] <= samples[i]);
            assert(correct[i] >= 0);
            int32_t maxSize = pow(10, i
                    /bucketSize);
            cumulativeCorrect += correct[i];
            cumulativeSamples += samples[i];
            fprintf(fileHandle, "\t<bucket from=\"%i\" to=\"%i\" correct=\"%i\" samples=\"%i\" correct_to_samples_ratio=\"%f\" cumulative_correct=\"%i\" cumulative_samples=\"%i\" cummulate_correct_to_samples_ratio=\"%f\"/>\n",
                    pMaxSize, maxSize, correct[i], samples[i], ((float) correct[i])/samples[i], cumulativeCorrect, cumulativeSamples, ((float) cumulativeCorrect)/cumulativeSamples);
            pMaxSize = maxSize + 1;
        }
    }
    fprintf(fileHandle, "\n</linkage_stats>\n");
    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    return 0;
}

