/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "cactus.h"
#include "assemblaCommon.h"
#include "pathsToBeds.h"

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "linkageStats");

    ///////////////////////////////////////////////////////////////////////////
    // Get the intervals
    ///////////////////////////////////////////////////////////////////////////

    stList *haplotypeEventStrings = getEventStrings(
            treatHaplotype1AsContamination ? NULL : hap1EventString,
            treatHaplotype2AsContamination ? NULL : hap2EventString);
    stList *assemblyEventStringInList = stList_construct();
    stList_append(assemblyEventStringInList, assemblyEventString);

    stList *intervals = stList_construct3(0, (void (*)(void *))sequenceInterval_destruct);
    for(int64_t i=0; i<stList_length(haplotypeEventStrings); i++) {
        const char *hapEventString = stList_get(haplotypeEventStrings, i);
        st_logInfo("Getting contig paths for haplotype: %s", hapEventString);
        stList *contigPaths = getContigPaths(flower, hapEventString, assemblyEventStringInList);
        stList *hapIntervals = getSplitContigPathIntervals(flower, contigPaths, hapEventString,
                assemblyEventStringInList);
        stList_destruct(contigPaths);
        st_logInfo("Getting contig paths\n");
        stList_appendAll(intervals, hapIntervals);
        stList_setDestructor(hapIntervals, NULL);
        stList_destruct(hapIntervals);
    }

    st_logDebug("Got a total of %" PRIi64 " intervals\n", stList_length(intervals));

    ///////////////////////////////////////////////////////////////////////////
    // Write it out.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    for (int64_t i = 0; i < stList_length(intervals); i++) {
        SequenceInterval *sequenceInterval = stList_get(intervals, i);
        st_logDebug("We have a path interval %s %" PRIi64 " %" PRIi64 "\n", sequenceInterval->sequenceName,
                sequenceInterval->start, sequenceInterval->end);
        fprintf(fileHandle, "%s %" PRIi64 " %" PRIi64 "\n", sequenceInterval->sequenceName,
                sequenceInterval->start, sequenceInterval->end);
    }

    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    return 0;
}
