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
    stList *contaminationEventStrings =
            getEventStrings(
                    contaminationEventString,
                    treatHaplotype1AsContamination ? hap1EventString
                            : (treatHaplotype2AsContamination ? hap2EventString
                                    : NULL));
    stList *assemblyEventStringInList = stList_construct();
    stList_append(assemblyEventStringInList, assemblyEventString);

    stList *intervals = stList_construct3(0, (void (*)(void *))sequenceInterval_destruct);
    for(int32_t i=0; i<stList_length(haplotypeEventStrings); i++) {
        const char *hapEventString = stList_get(haplotypeEventStrings, i);
        st_logInfo("Getting contig paths for haplotype: %s", hapEventString);
        stList *hapIntervals;
        if (reportContigPathIntervals) {
            hapIntervals = getContigPathIntervals(flower, hapEventString,
                    assemblyEventStringInList);
            st_logInfo("Getting contig paths\n");
        } else {
            hapIntervals = getScaffoldPathIntervals(flower, hapEventString,
                    assemblyEventStringInList, contaminationEventStrings,
                    capCodeParameters);
            st_logInfo("Getting scaffold paths\n");
        }
        stList_appendAll(intervals, hapIntervals);
        stList_setDestructor(hapIntervals, NULL);
        stList_destruct(hapIntervals);
    }

    st_logDebug("Got a total of %i intervals\n", stList_length(intervals));

    ///////////////////////////////////////////////////////////////////////////
    // Write it out.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    for (int32_t i = 0; i < stList_length(intervals); i++) {
        SequenceInterval *sequenceInterval = stList_get(intervals, i);
        st_logDebug("We have a path interval %s %i %i\n", sequenceInterval->sequenceName,
                sequenceInterval->start, sequenceInterval->end);
        fprintf(fileHandle, "%s %i %i\n", sequenceInterval->sequenceName,
                sequenceInterval->start, sequenceInterval->end);
    }

    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    return 0;
}
