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

int32_t totalSites = 0;
double totalCorrect = 0;
int32_t totalErrors = 0;
int32_t totalCalls = 0;

int32_t totalHeterozygous = 0;
double totalCorrectInHeterozygous = 0;
int32_t totalErrorsInHeterozygous = 0;
int32_t totalCallsInHeterozygous = 0;

int32_t totalCorrectHap1InHeterozygous = 0;
int32_t totalCorrectHap2InHeterozygous = 0;

int32_t totalInOneHaplotypeOnly = 0;
double totalCorrectInOneHaplotype = 0;
int32_t totalErrorsInOneHaplotype = 0;
int32_t totalCallsInOneHaplotype = 0;

stList *indelPositions = NULL;
stList *hetPositions = NULL;

typedef struct _segmentHolder {
    Segment *segment;
    int32_t offset;
} SegmentHolder;

SegmentHolder *segmentHolder_construct(Segment *segment, int32_t offset) {
    SegmentHolder *segmentHolder = st_malloc(sizeof(SegmentHolder));
    segmentHolder->segment = segment;
    segmentHolder->offset = offset;
    return segmentHolder;
}

void printPositions(stList *positions, const char *substitutionType, FILE *fileHandle) {
    for (int32_t i = 0; i < stList_length(positions); i++) {
        SegmentHolder *segmentHolder = stList_get(positions, i);
        int32_t j = segment_getStart(segmentHolder->segment);
        if (segment_getStrand(segmentHolder->segment)) {
            j += segmentHolder->offset;
            assert(
                    cap_getCoordinate(segment_get5Cap(segmentHolder->segment))
                            == segment_getStart(segmentHolder->segment));
            assert(
                    segment_getStart(segmentHolder->segment)
                            + segment_getLength(segmentHolder->segment) - 1
                            == cap_getCoordinate(
                                    segment_get3Cap(segmentHolder->segment)));
        } else {
            j -= segmentHolder->offset;
            assert(
                    segment_getStart(segmentHolder->segment)
                            - segment_getLength(segmentHolder->segment) + 1
                            == cap_getCoordinate(
                                    segment_get3Cap(segmentHolder->segment)));
        }

        fprintf(
                fileHandle,
                "%s: %s_%i %i\n",
                substitutionType,
                event_getHeader(segment_getEvent(segmentHolder->segment)),
                sequence_getLength(segment_getSequence(segmentHolder->segment)),
                j);
        getMAFBlock(segment_getBlock(segmentHolder->segment), fileHandle);
    }
}

static void getSnpStats(Block *block, FILE *fileHandle) {
    if (block_getLength(block) >= minimumBlockLength) {
        //Now get the column
        Block_InstanceIterator *instanceIterator = block_getInstanceIterator(
                block);
        Segment *segment;
        char *hap1Seq = NULL;
        char *hap2Seq = NULL;
        char *assemblySeq = NULL;
        Segment *hap1Segment = NULL;
        Segment *hap2Segment = NULL;
        while ((segment = block_getNext(instanceIterator)) != NULL) {
            if (strcmp(event_getHeader(segment_getEvent(segment)),
                    hap1EventString) == 0) {
                if (hap1Seq != NULL) {
                    goto end;
                }
                hap1Seq = segment_getString(segment);
                hap1Segment = segment;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)),
                    hap2EventString) == 0) {
                if (hap2Seq != NULL) {
                    goto end;
                }
                hap2Seq = segment_getString(segment);
                hap2Segment = segment;
            }
            if (strcmp(event_getHeader(segment_getEvent(segment)),
                    assemblyEventString) == 0) {
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
            for (int32_t i = ignoreFirstNBasesOfBlock; i < block_getLength(
                    block) - ignoreFirstNBasesOfBlock; i++) {
                if (hap1Seq != NULL && hap2Seq != NULL) {
                    if (toupper(hap1Seq[i]) == toupper(hap2Seq[i])) {
                        homoMatches++;
                    }
                } else {
                    homoMatches = INT32_MAX;
                }
                if (assemblySeq != NULL) {
                    if (hap1Seq != NULL) {
                        if (hap2Seq != NULL) {
                            if (toupper(hap1Seq[i]) == toupper(hap2Seq[i])
                                    && toupper(hap1Seq[i]) == toupper(
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
                    matches = INT32_MAX;
                }
            }
            double homoIdentity = 100.0 * homoMatches / (block_getLength(block)
                    - 2.0 * ignoreFirstNBasesOfBlock);
            double identity = 100.0 * matches / (block_getLength(block) - 2.0
                    * ignoreFirstNBasesOfBlock);

            if (homoIdentity >= minimumIndentity && identity
                    >= minimumIndentity) {
                //We're in gravy.
                for (int32_t i = ignoreFirstNBasesOfBlock; i < block_getLength(
                        block) - ignoreFirstNBasesOfBlock; i++) {

                    if (hap1Seq != NULL) {
                        if (hap2Seq != NULL) {
                            if (toupper(hap1Seq[i]) == toupper(hap2Seq[i])) {
                                totalSites++;
                                if (assemblySeq != NULL) {
                                    totalCorrect += bitsScoreFn(assemblySeq[i],
                                            hap1Seq[i]);
                                    totalErrors += correctFn(assemblySeq[i],
                                            hap1Seq[i]) ? 0 : 1;
                                    totalCalls++;
                                }
                            } else {
                                totalHeterozygous++;
                                if (assemblySeq != NULL) {
                                    assert(
                                            toupper(hap1Seq[i]) != toupper(
                                                    hap2Seq[i]));
                                    totalCorrectInHeterozygous += bitsScoreFn(
                                            assemblySeq[i], hap1Seq[i]);
                                    totalCorrectHap1InHeterozygous += bitsScoreFn(
                                            assemblySeq[i], hap1Seq[i]);
                                    totalCorrectInHeterozygous += bitsScoreFn(
                                            assemblySeq[i], hap2Seq[i]);
                                    totalCorrectHap2InHeterozygous += bitsScoreFn(
                                            assemblySeq[i], hap2Seq[i]);
                                    totalErrorsInHeterozygous += (correctFn(
                                            assemblySeq[i], hap1Seq[i])
                                            || correctFn(assemblySeq[i],
                                                    hap2Seq[i])) ? 0 : 1;
                                    totalCallsInHeterozygous++;
                                    if (!(correctFn(assemblySeq[i], hap1Seq[i])
                                            || correctFn(assemblySeq[i],
                                                    hap2Seq[i]))) {
                                        stList_append(
                                                hetPositions,
                                                segmentHolder_construct(
                                                        hap1Segment, i));
                                    }
                                }
                            }
                        } else {
                            totalInOneHaplotypeOnly++;
                            if (assemblySeq != NULL) {
                                totalCorrectInOneHaplotype += bitsScoreFn(
                                        assemblySeq[i], hap1Seq[i]);
                                totalErrorsInOneHaplotype += correctFn(
                                        assemblySeq[i], hap1Seq[i]) ? 0 : 1;
                                totalCallsInOneHaplotype++;
                                if (!correctFn(assemblySeq[i], hap1Seq[i])) {
                                    stList_append(
                                            indelPositions,
                                            segmentHolder_construct(
                                                    hap1Segment, i));
                                }
                            }
                        }
                    } else {
                        if (hap2Seq != NULL) {
                            totalInOneHaplotypeOnly++;
                            if (assemblySeq != NULL) {
                                totalCorrectInOneHaplotype += bitsScoreFn(
                                        assemblySeq[i], hap2Seq[i]);
                                totalErrorsInOneHaplotype += correctFn(
                                        assemblySeq[i], hap2Seq[i]) ? 0 : 1;
                                totalCallsInOneHaplotype++;
                                if (!correctFn(assemblySeq[i], hap2Seq[i])) {
                                    stList_append(
                                            indelPositions,
                                            segmentHolder_construct(
                                                    hap2Segment, i));
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
    getMAFsReferenceOrdered(flower, fileHandle, getSnpStats);

    ///////////////////////////////////////////////////////////////////////////
    // Print outputs
    ///////////////////////////////////////////////////////////////////////////

    fprintf(fileHandle, "<substitutionStats ");
    fprintf(fileHandle, "totalHomozygous=\"%i\" "
        "totalCorrectInHomozygous=\"%f\" "
        "totalErrorsInHomozygous=\"%i\" "
        "totalCallsInHomozygous=\"%i\" "
        "totalHeterozygous=\"%i\" "
        "totalCorrectInHeterozygous=\"%f\" "
        "totalErrorsInHeterozygous=\"%i\" "
        "totalCallsInHeterozygous=\"%i\" "
        "totalCorrectHap1InHeterozygous=\"%i\" "
        "totalCorrectHap2InHeterozygous=\"%i\" "
        "totalInOneHaplotypeOnly=\"%i\" "
        "totalCorrectInOneHaplotypeOnly=\"%f\" "
        "totalErrorsInOneHaplotypeOnly=\"%i\" "
        "totalCallsInOneHaplotypeOnly=\"%i\" />", totalSites, totalCorrect,
            totalErrors, totalCalls, totalHeterozygous,
            totalCorrectInHeterozygous, totalErrorsInHeterozygous,
            totalCallsInHeterozygous,
            totalCorrectHap1InHeterozygous,
            totalCorrectHap2InHeterozygous,
            totalInOneHaplotypeOnly,
            totalCorrectInOneHaplotype, totalErrorsInOneHaplotype,
            totalCallsInOneHaplotype);

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
