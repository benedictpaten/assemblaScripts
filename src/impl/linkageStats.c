/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "cactus.h"
#include "cactusMafs.h"
#include "adjacencyTraversal.h"
#include "assemblaCommon.h"
#include "linkage.h"

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "contiguityStats");

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    double bucketSize = bucketNumber / log10(upperLinkageBound);
    int64_t *correct = st_malloc(sizeof(int64_t) * bucketNumber);
    int64_t *aligned = st_malloc(sizeof(int64_t) * bucketNumber);
    int64_t *samples = st_malloc(sizeof(int64_t) * bucketNumber);
    for (int64_t i = 0; i < bucketNumber; i++) {
        correct[i] = 0;
        aligned[i] = 0;
        samples[i] = 0;
    }

    stList *eventStrings = getEventStrings(hap1EventString, hap2EventString);
    stSortedSet *sequences = getMetaSequencesForEvents(flower, eventStrings);
    stSortedSetIterator *it = stSortedSet_getIterator(sequences);
    stSortedSet *sortedSegments = getOrderedSegments(flower);
    MetaSequence *metaSequence;
    while ((metaSequence = stSortedSet_getNext(it)) != NULL) {
        samplePoints(flower, metaSequence, assemblyEventString,
                sampleNumber, correct, aligned, samples,
                bucketNumber, bucketSize, sortedSegments, 1, 1.0);
    }
    stSortedSet_destructIterator(it);
    stSortedSet_destruct(sequences);
    stSortedSet_destruct(sortedSegments);

    ///////////////////////////////////////////////////////////////////////////
    // Write it out.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    fprintf(fileHandle, "<linkage_stats>\n");
    int64_t pMaxSize = 1;
    int64_t cumulativeCorrect = 0;
    int64_t cumulativeSamples = 0;
    for (int64_t i = 0; i < bucketNumber; i++) {
        if(samples[i] > 0) {
            assert(correct[i] <= samples[i]);
            assert(correct[i] >= 0);
            int64_t maxSize = pow(10, i
                    /bucketSize);
            cumulativeCorrect += correct[i];
            cumulativeSamples += samples[i];
            fprintf(fileHandle, "\t<bucket from=\"%" PRIi64 "\" to=\"%" PRIi64 "\" correct=\"%" PRIi64 "\" samples=\"%" PRIi64 "\" correct_to_samples_ratio=\"%f\" cumulative_correct=\"%" PRIi64 "\" cumulative_samples=\"%" PRIi64 "\" cummulate_correct_to_samples_ratio=\"%f\"/>\n",
                    pMaxSize, maxSize, correct[i], samples[i], ((float) correct[i])/samples[i], cumulativeCorrect, cumulativeSamples, ((float) cumulativeCorrect)/cumulativeSamples);
            pMaxSize = maxSize + 1;
        }
    }
    fprintf(fileHandle, "\n</linkage_stats>\n");
    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    return 0;
}

