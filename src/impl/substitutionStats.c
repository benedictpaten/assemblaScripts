/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <ctype.h>

#include "sonLib.h"
#include "cactus.h"
#include "substitutions.h"
#include "adjacencyClassification.h"
#include "scaffoldPaths.h"
#include "assemblaCommon.h"
#include "cactusMafs.h"

int64_t totalSites = 0;
double totalCorrect = 0;
int64_t totalErrors = 0;
int64_t totalCalls = 0;

int64_t totalHeterozygous = 0;
double totalCorrectInHeterozygous = 0;
int64_t totalErrorsInHeterozygous = 0;
int64_t totalCallsInHeterozygous = 0;

int64_t totalCorrectHap1InHeterozygous = 0;
int64_t totalCorrectHap2InHeterozygous = 0;

int64_t totalInOneHaplotypeOnly = 0;
double totalCorrectInOneHaplotype = 0;
int64_t totalErrorsInOneHaplotype = 0;
int64_t totalCallsInOneHaplotype = 0;

stList *indelPositions = NULL;
stList *hetPositions = NULL;

typedef struct _segmentHolder {
    Segment *segment;
    int64_t offset;
    char base1;
    char base2;
    char base3;
} SegmentHolder;

SegmentHolder *segmentHolder_construct(Segment *segment, int64_t offset, char base1, char base2, char base3) {
    SegmentHolder *segmentHolder = st_malloc(sizeof(SegmentHolder));
    segmentHolder->segment = segment;
    segmentHolder->offset = offset;
    segmentHolder->base1 = base1;
    segmentHolder->base2 = base2;
    segmentHolder->base3 = base3;
    return segmentHolder;
}

void printPositions(stList *positions, const char *substitutionType, FILE *fileHandle) {
    for (int64_t i = 0; i < stList_length(positions); i++) {
        SegmentHolder *segmentHolder = stList_get(positions, i);
        int64_t j = segment_getStart(segmentHolder->segment);
        if (segment_getStrand(segmentHolder->segment)) {
            j += segmentHolder->offset;
            assert(
                    cap_getCoordinate(segment_get5Cap(segmentHolder->segment)) == segment_getStart(
                            segmentHolder->segment));
            assert(
                    segment_getStart(segmentHolder->segment) + segment_getLength(segmentHolder->segment) - 1
                            == cap_getCoordinate(segment_get3Cap(segmentHolder->segment)));
        } else {
            j -= segmentHolder->offset;
            assert(
                    segment_getStart(segmentHolder->segment) - segment_getLength(segmentHolder->segment) + 1
                            == cap_getCoordinate(segment_get3Cap(segmentHolder->segment)));
        }

        fprintf(fileHandle, "%s: %s_%" PRIi64 " %" PRIi64 " %c %c %c\n", substitutionType,
                event_getHeader(segment_getEvent(segmentHolder->segment)),
                sequence_getLength(segment_getSequence(segmentHolder->segment)), j,
                segmentHolder->base1, segmentHolder->base2, segmentHolder->base3);
        getMAFBlock(segment_getBlock(segmentHolder->segment), fileHandle);
    }
}

static void getSnpStats(Block *block, FILE *fileHandle) {
    if (block_getLength(block) >= minimumBlockLength) {
        //Now get the column
        Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
        Segment *segment;
        char *hap1Seq = NULL;
        char *hap2Seq = NULL;
        char *assemblySeq = NULL;
        Segment *hap1Segment = NULL;
        Segment *hap2Segment = NULL;
        while ((segment = block_getNext(instanceIterator)) != NULL) {
            if (strcmp(event_getHeader(segment_getEvent(segment)), hap1EventString) == 0) {
                if (hap1Seq != NULL) {
                    goto end;
                }
                hap1Seq = segment_getString(segment);
                hap1Segment = segment;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)), hap2EventString) == 0) {
                if (hap2Seq != NULL) {
                    goto end;
                }
                hap2Seq = segment_getString(segment);
                hap2Segment = segment;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)), assemblyEventString) == 0) {
                if (assemblySeq != NULL) {
                    goto end;
                }
                assemblySeq = segment_getString(segment);
            }
        }

        assert(minimumIndentity >= 0);
        assert(minimumIndentity <= 100);
        if (hap1Seq != NULL || hap2Seq != NULL) {
            if (hap1Seq != NULL) {
                assert(strlen(hap1Seq) == block_getLength(block));
            }
            if (hap2Seq != NULL) {
                assert(strlen(hap2Seq) == block_getLength(block));
            }
            if (assemblySeq != NULL) {
                assert(strlen(assemblySeq) == block_getLength(block));
            }
            double homoMatches = 0;
            double matches = 0;
            for (int64_t i = ignoreFirstNBasesOfBlock; i < block_getLength(block) - ignoreFirstNBasesOfBlock; i++) {
                if (hap1Seq != NULL && hap2Seq != NULL) {
                    if (toupper(hap1Seq[i]) == toupper(hap2Seq[i])) {
                        homoMatches++;
                    }
                } else {
                    homoMatches = INT64_MAX;
                }
                if (assemblySeq != NULL) {
                    if (hap1Seq != NULL) {
                        if (hap2Seq != NULL) {
                            if (toupper(hap1Seq[i]) == toupper(hap2Seq[i]) && toupper(hap1Seq[i]) == toupper(
                                    assemblySeq[i])) {
                                matches++;
                            }
                        } else {
                            if (toupper(hap1Seq[i]) == toupper(assemblySeq[i])) {
                                matches++;
                            }
                        }
                    } else {
                        assert(hap2Seq != NULL);
                        if (toupper(hap2Seq[i]) == toupper(assemblySeq[i])) {
                            matches++;
                        }
                    }
                } else {
                    matches = INT64_MAX;
                }
            }
            double homoIdentity = 100.0 * homoMatches / (block_getLength(block) - 2.0 * ignoreFirstNBasesOfBlock);
            double identity = 100.0 * matches / (block_getLength(block) - 2.0 * ignoreFirstNBasesOfBlock);

            if (homoIdentity >= minimumIndentity && identity >= minimumIndentity) {
                //We're in gravy.
                for (int64_t i = ignoreFirstNBasesOfBlock; i < block_getLength(block) - ignoreFirstNBasesOfBlock; i++) {

                    if (hap1Seq != NULL) {
                        if (hap2Seq != NULL) {
                            if (toupper(hap1Seq[i]) == toupper(hap2Seq[i])) {
                                totalSites++;
                                if (assemblySeq != NULL) {
                                    totalCorrect += bitsScoreFn(assemblySeq[i], hap1Seq[i]);
                                    totalErrors += correctFn(assemblySeq[i], hap1Seq[i]) ? 0 : 1;
                                    totalCalls++;
                                }
                            } else {
                                totalHeterozygous++;
                                if (assemblySeq != NULL) {
                                    assert(toupper(hap1Seq[i]) != toupper(hap2Seq[i]));
                                    totalCorrectInHeterozygous += bitsScoreFn(assemblySeq[i], hap1Seq[i]);
                                    totalCorrectHap1InHeterozygous += bitsScoreFn(assemblySeq[i], hap1Seq[i]);
                                    totalCorrectInHeterozygous += bitsScoreFn(assemblySeq[i], hap2Seq[i]);
                                    totalCorrectHap2InHeterozygous += bitsScoreFn(assemblySeq[i], hap2Seq[i]);
                                    totalErrorsInHeterozygous += (correctFn(assemblySeq[i], hap1Seq[i]) || correctFn(
                                            assemblySeq[i], hap2Seq[i])) ? 0 : 1;
                                    totalCallsInHeterozygous++;
                                    if (!(correctFn(assemblySeq[i], hap1Seq[i])
                                            || correctFn(assemblySeq[i], hap2Seq[i]))) {
                                        stList_append(hetPositions, segmentHolder_construct(hap1Segment, i, assemblySeq[i], hap1Seq[i], hap2Seq[i]));
                                    }
                                }
                            }
                        } else {
                            totalInOneHaplotypeOnly++;
                            if (assemblySeq != NULL) {
                                totalCorrectInOneHaplotype += bitsScoreFn(assemblySeq[i], hap1Seq[i]);
                                totalErrorsInOneHaplotype += correctFn(assemblySeq[i], hap1Seq[i]) ? 0 : 1;
                                totalCallsInOneHaplotype++;
                                if (!correctFn(assemblySeq[i], hap1Seq[i])) {
                                    stList_append(indelPositions, segmentHolder_construct(hap1Segment, i, assemblySeq[i], hap1Seq[i], 'N'));
                                }
                            }
                        }
                    } else {
                        if (hap2Seq != NULL) {
                            totalInOneHaplotypeOnly++;
                            if (assemblySeq != NULL) {
                                totalCorrectInOneHaplotype += bitsScoreFn(assemblySeq[i], hap2Seq[i]);
                                totalErrorsInOneHaplotype += correctFn(assemblySeq[i], hap2Seq[i]) ? 0 : 1;
                                totalCallsInOneHaplotype++;
                                if (!correctFn(assemblySeq[i], hap2Seq[i])) {
                                    stList_append(indelPositions, segmentHolder_construct(hap2Segment, i, assemblySeq[i], 'N', hap2Seq[i]));
                                }
                            }
                        }
                    }
                }
            }
        }

        end:
        //cleanup
        if (hap1Seq != NULL) {
            free(hap1Seq);
        }
        if (hap2Seq != NULL) {
            free(hap2Seq);
        }
        if (assemblySeq != NULL) {
            free(assemblySeq);
        }
        block_destructInstanceIterator(instanceIterator);
    }
}

int main(int argc, char *argv[]) {

    assert(correctFn('Y', 'C'));
    assert(correctFn('Y', 'T'));
    assert(!correctFn('Y', 'G'));
    assert(!correctFn('Y', 'A'));

    assert(correctFn('y', 'C'));
    assert(correctFn('y', 'T'));
    assert(!correctFn('y', 'G'));
    assert(!correctFn('y', 'A'));

    assert(correctFn('Y', 'c'));
    assert(correctFn('Y', 't'));
    assert(!correctFn('Y', 'g'));
    assert(!correctFn('Y', 'a'));

    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "snpStats");

    ///////////////////////////////////////////////////////////////////////////
    // Calculate and print to file a crap load of numbers.
    ///////////////////////////////////////////////////////////////////////////

    indelPositions = stList_construct3(0, NULL);
    hetPositions = stList_construct3(0, NULL);
    FILE *fileHandle = fopen(outputFile, "w");
    getMAFs(flower, fileHandle, getSnpStats);

    ///////////////////////////////////////////////////////////////////////////
    // Print outputs
    ///////////////////////////////////////////////////////////////////////////

    fprintf(fileHandle, "<substitutionStats ");
    fprintf(fileHandle, "totalHomozygous=\"%" PRIi64 "\" "
        "totalCorrectInHomozygous=\"%f\" "
        "totalErrorsInHomozygous=\"%" PRIi64 "\" "
        "totalCallsInHomozygous=\"%" PRIi64 "\" "
        "totalHeterozygous=\"%" PRIi64 "\" "
        "totalCorrectInHeterozygous=\"%f\" "
        "totalErrorsInHeterozygous=\"%" PRIi64 "\" "
        "totalCallsInHeterozygous=\"%" PRIi64 "\" "
        "totalCorrectHap1InHeterozygous=\"%" PRIi64 "\" "
        "totalCorrectHap2InHeterozygous=\"%" PRIi64 "\" "
        "totalInOneHaplotypeOnly=\"%" PRIi64 "\" "
        "totalCorrectInOneHaplotypeOnly=\"%f\" "
        "totalErrorsInOneHaplotypeOnly=\"%" PRIi64 "\" "
        "totalCallsInOneHaplotypeOnly=\"%" PRIi64 "\" />", totalSites, totalCorrect, totalErrors, totalCalls, totalHeterozygous,
            totalCorrectInHeterozygous, totalErrorsInHeterozygous, totalCallsInHeterozygous,
            totalCorrectHap1InHeterozygous, totalCorrectHap2InHeterozygous, totalInOneHaplotypeOnly,
            totalCorrectInOneHaplotype, totalErrorsInOneHaplotype, totalCallsInOneHaplotype);

    if (printIndelPositions) {
        printPositions(indelPositions, "INDEL_SUBSTITUTION", fileHandle);
    }

    if (printHetPositions) {
        printPositions(hetPositions, "HET_SUBSTITUTION", fileHandle);
    }

    st_logInfo("Finished writing out the stats.\n");
    fclose(fileHandle);

    return 0;
}
