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

int64_t totalHapSwitches = 0;

int64_t totalCleanEnds = 0;
int64_t totalHangingEndWithsNs = 0;

int64_t totalScaffoldGaps = 0;
int64_t totalAmbiguityGaps = 0;

int64_t totalHaplotypeLength = 0;

int64_t totalErrorsHapToHapSameChromosome = 0;
int64_t totalErrorsInterJoin = 0;
int64_t totalErrorsHapToContamination = 0;
int64_t totalErrorsHapToInsertToContamination = 0;
int64_t totalErrorsInsertion = 0;
int64_t totalErrorsDeleteion = 0;
int64_t totalErrorsInsertionAndDeletion = 0;
int64_t totalErrorsHangingInsertion = 0;
stList *insertionDistribution = NULL;
stList *deletionDistribution = NULL;

void reportHaplotypePathStatsP(Cap *cap, stList *haplotypeEventStrings, stList *contaminationEventStrings, CapCodeParameters *capCodeParameters) {
    int64_t insertLength, deleteLength;
    Cap *otherCap;
    switch (getCapCode(cap, &otherCap, haplotypeEventStrings, contaminationEventStrings, &insertLength, &deleteLength, capCodeParameters)) {
        case HAP_SWITCH:
            totalHapSwitches++;
            return;
        case HAP_NOTHING:
            return;
        case CONTIG_END:
            totalCleanEnds++;
            return;
        case CONTIG_END_WITH_SCAFFOLD_GAP:
        case CONTIG_END_WITH_AMBIGUITY_GAP:
            totalHangingEndWithsNs++;
            return;
        case SCAFFOLD_GAP:
            totalScaffoldGaps++;
            return;
        case AMBIGUITY_GAP:
            totalAmbiguityGaps++;
            return;
        case ERROR_HAP_TO_HAP_SAME_CHROMOSOME:
            totalErrorsHapToHapSameChromosome++;
            return;
        case ERROR_HAP_TO_HAP_DIFFERENT_CHROMOSOMES:
            totalErrorsInterJoin++;
            return;
        case ERROR_HAP_TO_CONTAMINATION:
            totalErrorsHapToContamination++;
            return;
        case ERROR_HAP_TO_INSERT_TO_CONTAMINATION:
            totalErrorsHapToInsertToContamination++;
            return;
        case ERROR_HAP_TO_INSERT:
            assert(insertLength > 0);
            stList_append(insertionDistribution, stIntTuple_construct1( insertLength));
            totalErrorsInsertion++;
            return;
        case ERROR_HAP_TO_DELETION:
            assert(deleteLength > 0);
            stList_append(deletionDistribution, stIntTuple_construct1( deleteLength));
            totalErrorsDeleteion++;
            return;
        case ERROR_HAP_TO_INSERT_AND_DELETION:
            assert(insertLength > 0);
            assert(deleteLength > 0);
            stList_append(insertionDistribution, stIntTuple_construct1( insertLength));
            stList_append(deletionDistribution, stIntTuple_construct1( deleteLength));
            totalErrorsInsertionAndDeletion++;
            return;
        case ERROR_CONTIG_END_WITH_INSERT:
            totalErrorsHangingInsertion++;
            return;
    }
}

static int64_t totalPathLength = 0; //sum of all blocks containing haplotype and assembly.
static stSortedSet *contigsSet;
static stSortedSet *haplotypesSet;
static stList *blockList;

int compareSequences(const void *a, const void *b) {
    return cactusMisc_nameCompare(sequence_getName((Sequence *) a), sequence_getName((Sequence *) b));
}

static int compareSequencesByLength(const void *a, const void *b) {
    return sequence_getLength((Sequence *) b) - sequence_getLength((Sequence *) a);
}

static int compareBlocksByLength(const void *a, const void *b) {
    return block_getLength((Block *) b) - block_getLength((Block *) a);
}

static stHash *maximalHaplotypePathToLength;

static int64_t getHaplotypePathLength(const void *a) {
    return stIntTuple_getPosition(stHash_search(maximalHaplotypePathToLength, (void *) a), 0);
}

static int compareMaximalHaplotypePaths(const void *a, const void *b) {
    return getHaplotypePathLength(b) - getHaplotypePathLength(a);
}

static stHash *maximalScaffoldPathToLength;

static int64_t getScaffoldPathLength(const void *a) {
    return stIntTuple_getPosition(stHash_search(maximalScaffoldPathToLength, (void *) a), 0);
}

static int compareScaffoldPaths(const void *a, const void *b) {
    return getScaffoldPathLength(b) - getScaffoldPathLength(a);
}

void accumulateBlock(Block *block, const char *assemblyEventString, stList *eventStrings) {
    if (hasCapInEvents(block_get5End(block), eventStrings)) {
        //genotypeLength += block_getLength(block);
        if (hasCapInEvent(block_get5End(block), assemblyEventString)) {
            stList_append(blockList, block);
            totalPathLength += block_getLength(block);
        }
    }
    Segment *segment;
    Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
    while ((segment = block_getNext(instanceIt)) != NULL) {
        Sequence *sequence = segment_getSequence(segment);
        assert(sequence != NULL);
        if (strcmp(event_getHeader(segment_getEvent(segment)), assemblyEventString) == 0) {
            stSortedSet_insert(contigsSet, sequence);
        }
        for(int64_t i=0; i<stList_length(eventStrings); i++) {
            if(strcmp(event_getHeader(segment_getEvent(segment)), stList_get(eventStrings, i)) == 0) {
                stSortedSet_insert(haplotypesSet, sequence);
                break;
            }
        }
        //assert(stList_length(eventStrings) == 2);
        //if (strcmp(event_getHeader(segment_getEvent(segment)), stList_get(eventStrings, 0)) == 0 ||
        //        strcmp(event_getHeader(segment_getEvent(segment)), stList_get(eventStrings, 1)) == 0) {
        //    stSortedSet_insert(haplotypesSet, sequence);
        //}
    }
    block_destructInstanceIterator(instanceIt);
}

void traverseBlocks(Flower *flower, const char *assemblyEventString, stList *eventStrings) {
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (group_getNestedFlower(group) != NULL) {
            traverseBlocks(group_getNestedFlower(group), assemblyEventString, eventStrings);
        }
    }
    flower_destructGroupIterator(groupIt);

    Flower_BlockIterator *blockIt = flower_getBlockIterator(flower);
    Block *block;
    while ((block = flower_getNextBlock(blockIt)) != NULL) {
        accumulateBlock(block, assemblyEventString, eventStrings);
    }
    flower_destructBlockIterator(blockIt);
}

int64_t getN50(int64_t genomeLength, stList *objects, int64_t(*lengthFn)(const void *)) {
    int64_t totalLength = 0;
    int64_t pJ = INT32_MAX;
    for (int64_t i = 0; i < stList_length(objects); i++) {
        int64_t j = lengthFn(stList_get(objects, i));
        assert(j <= pJ);
        pJ = j;
        totalLength += j;
        if (totalLength >= genomeLength / 2) {
            return j;
        }
    }
    return -1; //this should not happen!
}

stList *getScaffoldPathsList(stList *maximalHaplotypePaths, stList *haplotypeEventStrings, stList *contaminationEventStrings,CapCodeParameters *capCodeParameters) {
    stHash *scaffoldPaths = getScaffoldPaths(maximalHaplotypePaths, haplotypeEventStrings, contaminationEventStrings,capCodeParameters);
    stSortedSet *bucketSet = stSortedSet_construct();
    stList *scaffoldPaths2 = stList_construct();
    stHashIterator *hashIt = stHash_getIterator(scaffoldPaths);
    stList *haplotypePath;
    while ((haplotypePath = stHash_getNext(hashIt)) != NULL) {
        assert(stList_length(haplotypePath) > 0);
        stSortedSet *bucket = stHash_search(scaffoldPaths, haplotypePath);
        assert(bucket != NULL);
        if (stSortedSet_search(bucketSet, bucket) == NULL) {
            stSortedSet_insert(bucketSet, bucket);
            stList_append(scaffoldPaths2, haplotypePath);
        }
    }
    stHash_destructIterator(hashIt);
    stSortedSet_destruct(bucketSet);
    return scaffoldPaths2;
}

char *concatenateList(stList *list) {
    char **cAA = st_malloc(sizeof(char *) * (stList_length(list)/2));
    stList_sort(list, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    assert(stList_length(list) % 2 == 0);
    for (int64_t i = 0; i < stList_length(list); i+=2) {
        assert(stIntTuple_getPosition(stList_get(list, i), 0) == stIntTuple_getPosition(stList_get(list, i+1), 0));
        cAA[i/2] = stString_print("%" PRIi64 "", stIntTuple_getPosition(stList_get(list, i), 0));
    }
    char *cA = stString_join(" ", (const char **)cAA, stList_length(list)/2);
    for (int64_t i = 0; i < stList_length(list)/2; i++) {
        free(cAA[i]);
    }
    free(cAA);
    return cA;
}

void reportSamplePathStats(Flower *flower, FILE *fileHandle,
        const char *assemblyEventString,
        stList *haplotypeEventStrings, stList *contaminationEventStrings, CapCodeParameters *capCodeParameters) {
    /*
     * Gets stats on the maximal haplotype paths.
     */

    stList *maximalHaplotypePaths = getContigPaths(flower, assemblyEventString, haplotypeEventStrings);
    maximalHaplotypePathToLength = buildContigPathToContigPathLengthHash(maximalHaplotypePaths);
    maximalScaffoldPathToLength = getContigPathToScaffoldPathLengthsHash(maximalHaplotypePaths, haplotypeEventStrings, contaminationEventStrings, capCodeParameters);

    //initialise the global arrays
    insertionDistribution = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    deletionDistribution = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);

    for (int64_t i = 0; i < stList_length(maximalHaplotypePaths); i++) {
        stList *maximalHaplotypePath = stList_get(maximalHaplotypePaths, i);
        totalHaplotypeLength += contigPathLength(maximalHaplotypePath);
        for (int64_t j = 0; j < stList_length(maximalHaplotypePath); j++) {
            Segment *segment = stList_get(maximalHaplotypePath, j);
            reportHaplotypePathStatsP(segment_get5Cap(segment), haplotypeEventStrings, contaminationEventStrings, capCodeParameters);
            reportHaplotypePathStatsP(segment_get3Cap(segment), haplotypeEventStrings, contaminationEventStrings, capCodeParameters);
        }
    }

    blockList = stList_construct();
    contigsSet = stSortedSet_construct3(compareSequences, NULL);
    haplotypesSet = stSortedSet_construct3(compareSequences, NULL);

    traverseBlocks(flower, assemblyEventString, haplotypeEventStrings);

    stList *sequences = stSortedSet_getList(contigsSet);
    stList *haplotypes = stSortedSet_getList(haplotypesSet);
    stList_sort(blockList, compareBlocksByLength);
    stList_sort(sequences, compareSequencesByLength);
    stList_sort(maximalHaplotypePaths, compareMaximalHaplotypePaths);
    stList *scaffoldPaths = getScaffoldPathsList(maximalHaplotypePaths, haplotypeEventStrings, contaminationEventStrings, capCodeParameters);
    stList_sort(scaffoldPaths, compareScaffoldPaths);

    int64_t totalSequencesLength = 0;
    for (int64_t i = 0; i < stList_length(sequences); i++) {
        totalSequencesLength += sequence_getLength(stList_get(sequences, i));
    }

    int64_t averageHaplotypeLength = 0; //average length of the two haplotypes.
    for (int64_t i = 0; i < stList_length(haplotypes); i++) {
        averageHaplotypeLength += sequence_getLength(stList_get(haplotypes, i));
    }
    averageHaplotypeLength /= stList_length(haplotypeEventStrings);

    int64_t totalBlockNumber = stList_length(blockList);
    int64_t blockNG50 = getN50(averageHaplotypeLength, blockList, (int64_t(*)(const void *)) block_getLength);
    int64_t contigN50 = getN50(totalSequencesLength, sequences, (int64_t(*)(const void *)) sequence_getLength); //length of contig which appears in order (from longest to shortest) at 50% coverage.
    int64_t contigNG50 = getN50(averageHaplotypeLength, sequences, (int64_t(*)(const void *)) sequence_getLength); //length of contig which appears in order (from longest to shortest) at 50% coverage.
    int64_t haplotypePathNG50 = getN50(averageHaplotypeLength, maximalHaplotypePaths, getHaplotypePathLength); //length of maximal haplotype path which appears in order (from longest to shortest) at 50% coverage.
    int64_t scaffoldPathNG50 = getN50(averageHaplotypeLength, scaffoldPaths, getScaffoldPathLength); //length of maximal scaffold path which appears in order (from longest to shortest) at 50% coverage.

    int64_t totalContigNumber = stSortedSet_size(contigsSet); //number of contigs
    int64_t totalHaplotypePaths = stList_length(maximalHaplotypePaths); //number of haplotype paths
    int64_t totalScaffoldPaths = stList_length(scaffoldPaths); //number of scaffold paths

    int64_t totalErrors = totalErrorsHapToHapSameChromosome / 2 + totalErrorsInterJoin / 2
            + totalErrorsHapToContamination + totalErrorsHapToInsertToContamination + totalErrorsInsertion / 2
            + totalErrorsDeleteion / 2 + totalErrorsInsertionAndDeletion / 2 + totalErrorsHangingInsertion;

    double errorsPerContig = ((double) totalErrors) / totalContigNumber; //number of errors / number of contigs
    double errorsPerMappedBase = ((double) totalErrors) / totalPathLength; //number of errors / total path length
    double coverage = ((double) totalPathLength) / averageHaplotypeLength; //totalPathLength/genotypeLength


    assert(totalErrorsHapToHapSameChromosome % 2 == 0);
    assert(totalErrorsInterJoin % 2 == 0);
    assert(totalErrorsInsertion % 2 == 0);
    assert(totalErrorsDeleteion % 2 == 0);
    assert(totalErrorsInsertionAndDeletion % 2 == 0);
    assert(totalScaffoldGaps % 2 == 0);
    assert(totalAmbiguityGaps % 2 == 0);

    char *insertionDistributionString = concatenateList(insertionDistribution);
    char *deletionDistributionString = concatenateList(deletionDistribution);

    fprintf(fileHandle, "<stats totalHaplotypeSwitches=\"%" PRIi64 "\" "
        "totalScaffoldGaps=\"%" PRIi64 "\" "
        "totalAmbiguityGaps=\"%" PRIi64 "\" "
        "totalContigEnds=\"%" PRIi64 "\" "
        "totalContigEndsWithNs=\"%" PRIi64 "\" "
        "totalErrorsHaplotypeToHaplotypeSameChromosome=\"%" PRIi64 "\" "
        "totalErrorsHaplotypeToHaplotypeDifferentChromosome=\"%" PRIi64 "\" "
        "totalErrorsHaplotypeToContamination=\"%" PRIi64 "\" "
        "totalErrorsHaplotypeToInsertionToContamination=\"%" PRIi64 "\" "
        "totalErrorsHaplotypeToInsertion=\"%" PRIi64 "\" "
        "totalErrorsHaplotypeToDeletion=\"%" PRIi64 "\" "
        "totalErrorsHaplotypeToInsertionAndDeletion=\"%" PRIi64 "\" "
        "totalErrorsContigEndsWithInsert=\"%" PRIi64 "\" "
        "totalErrors=\"%" PRIi64 "\" totalPathLength=\"%" PRIi64 "\" genotypeLength=\"%" PRIi64 "\" "
        "totalContigsLength=\"%" PRIi64 "\" coverage=\"%f\" blockNG50=\"%" PRIi64 "\" contigN50=\"%" PRIi64 "\" "
        "contigNG50=\"%" PRIi64 "\" contigPathNG50=\"%" PRIi64 "\" scaffoldPathNG50=\"%" PRIi64 "\" totalBlockNumber=\"%" PRIi64 "\" "
        "totalContigNumber=\"%" PRIi64 "\" totalHaplotypePaths=\"%" PRIi64 "\" totalScaffoldPaths=\"%" PRIi64 "\" "
        "errorsPerContig=\"%f\" errorsPerMappedBase=\"%f\" "
        "insertionErrorSizeDistribution=\"%s\" "
        "deletionErrorSizeDistribution=\"%s\"/>", totalHapSwitches / 2,
            totalScaffoldGaps / 2, totalAmbiguityGaps / 2, totalCleanEnds,
            totalHangingEndWithsNs, totalErrorsHapToHapSameChromosome / 2, totalErrorsInterJoin / 2,
            totalErrorsHapToContamination, totalErrorsHapToInsertToContamination, totalErrorsInsertion / 2,
            totalErrorsDeleteion / 2, totalErrorsInsertionAndDeletion / 2, totalErrorsHangingInsertion,
            totalErrors, totalPathLength, averageHaplotypeLength, totalSequencesLength, coverage, blockNG50, contigN50,
            contigNG50, haplotypePathNG50, scaffoldPathNG50, totalBlockNumber, totalContigNumber, totalHaplotypePaths,
            totalScaffoldPaths, errorsPerContig, errorsPerMappedBase,
            insertionDistributionString, deletionDistributionString);

    stList_destruct(maximalHaplotypePaths);
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "pathStats");

    ///////////////////////////////////////////////////////////////////////////
    // Now print the haplotype path stats.
    ///////////////////////////////////////////////////////////////////////////

    int64_t startTime = time(NULL);
    FILE *fileHandle = fopen(outputFile, "w");

    assert(!(treatHaplotype1AsContamination && treatHaplotype2AsContamination));

    stList *haplotypeEventStrings = getEventStrings(treatHaplotype1AsContamination ? NULL : hap1EventString, treatHaplotype2AsContamination ? NULL : hap2EventString);
    stList *contaminationEventStrings = getEventStrings(contaminationEventString, treatHaplotype1AsContamination ? hap1EventString : (treatHaplotype2AsContamination ? hap2EventString : NULL));

    reportSamplePathStats(flower, fileHandle, assemblyEventString, haplotypeEventStrings, contaminationEventStrings, capCodeParameters);
    fclose(fileHandle);
    st_logInfo("Got the stats in %" PRIi64 " seconds/\n", time(NULL) - startTime);

    return 0;
}
