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
    stList *intervals;
    if (reportContigPathIntervals) {
        intervals = getContigPathIntervals(flower, assemblyEventString,
                haplotypeEventStrings);
        st_logInfo("Getting contig paths\n");
    } else {
        intervals = getScaffoldPathIntervals(flower, assemblyEventString,
                haplotypeEventStrings, contaminationEventStrings,
                capCodeParameters);
        st_logInfo("Getting scaffold paths\n");
    }

    st_logDebug("Got a total of %i intervals\n", stList_length(intervals));

    ///////////////////////////////////////////////////////////////////////////
    // Write it out.
    ///////////////////////////////////////////////////////////////////////////

    FILE *fileHandle = fopen(outputFile, "w");
    for (int32_t i = 0; i < stList_length(intervals); i++) {
        SequenceInterval *sequenceInterval = stList_get(intervals, i);
        st_logDebug("We have a path interval %s %i %i\n", sequenceInterval->sequenceName,
                sequenceInterval->start, sequenceInterval->length);
        fprintf(fileHandle, "%s %i %i\n", sequenceInterval->sequenceName,
                sequenceInterval->start, sequenceInterval->length);
    }

    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    return 0;
}
