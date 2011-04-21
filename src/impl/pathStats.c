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

int32_t totalHapSwitches = 0;

int32_t totalContigEnds = 0;
int32_t totalContigEndsWithsNs = 0;

int32_t totalScaffoldGaps = 0;
int32_t totalAmbiguityGaps = 0;

int32_t totalHaplotypeLength = 0;

int32_t totalErrorsHapToHapSameChromosome = 0;
int32_t totalErrorsHapToHapDifferentChromosome = 0;
int32_t totalErrorsHapToContamination = 0;
int32_t totalErrorsHapToInsertToContamination = 0;
int32_t totalErrorsHapToInsert = 0;
int32_t totalErrorsHapToDeletion = 0;
int32_t totalErrorsHapToInsertionAndDeletion = 0;
int32_t totalErrorsContigEndsWithInsert = 0;
stList *insertionDistribution = NULL;
stList *deletionDistribution = NULL;

void reportHaplotypePathStatsP(Cap *cap, stList *eventStrings, CapCodeParameters *capCodeParameters) {
    int32_t insertLength, deleteLength;
    switch (getCapCode(cap, eventStrings, &insertLength, &deleteLength, capCodeParameters)) {
        case HAP_SWITCH:
            totalHapSwitches++;
            return;
        case HAP_NOTHING:
            return;
        case CONTIG_END:
            totalContigEnds++;
            return;
        case CONTIG_END_WITH_SCAFFOLD_GAP:
        case CONTIG_END_WITH_AMBIGUITY_GAP:
            totalContigEndsWithsNs++;
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
            totalErrorsHapToHapDifferentChromosome++;
            return;
        case ERROR_HAP_TO_CONTAMINATION:
            totalErrorsHapToContamination++;
            return;
        case ERROR_HAP_TO_INSERT_TO_CONTAMINATION:
            totalErrorsHapToInsertToContamination++;
            return;
        case ERROR_HAP_TO_INSERT:
            assert(insertLength > 0);
            stList_append(insertionDistribution, stIntTuple_construct(1, insertLength));
            totalErrorsHapToInsert++;
            return;
        case ERROR_HAP_TO_DELETION:
            assert(deleteLength > 0);
            stList_append(deletionDistribution, stIntTuple_construct(1, deleteLength));
            totalErrorsHapToDeletion++;
            return;
        case ERROR_HAP_TO_INSERT_AND_DELETION:
            assert(insertLength > 0);
            assert(deleteLength > 0);
            stList_append(insertionDistribution, stIntTuple_construct(1, insertLength));
            stList_append(deletionDistribution, stIntTuple_construct(1, deleteLength));
            totalErrorsHapToInsertionAndDeletion++;
            return;
        case ERROR_CONTIG_END_WITH_INSERT:
            totalErrorsContigEndsWithInsert++;
            return;
    }
}

static int32_t totalPathLength = 0; //sum of all blocks containing haplotype and assembly.
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

static int32_t getHaplotypePathLength(const void *a) {
    return stIntTuple_getPosition(stHash_search(maximalHaplotypePathToLength, (void *) a), 0);
}

static int compareMaximalHaplotypePaths(const void *a, const void *b) {
    return getHaplotypePathLength(b) - getHaplotypePathLength(a);
}

static stHash *maximalScaffoldPathToLength;

static int32_t getScaffoldPathLength(const void *a) {
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
        assert(stList_length(eventStrings) == 2);
        if (strcmp(event_getHeader(segment_getEvent(segment)), stList_get(eventStrings, 0)) == 0 ||
                strcmp(event_getHeader(segment_getEvent(segment)), stList_get(eventStrings, 1)) == 0) {
            stSortedSet_insert(haplotypesSet, sequence);
        }
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

int32_t getN50(int32_t genomeLength, stList *objects, int32_t(*lengthFn)(const void *)) {
    int32_t totalLength = 0;
    int32_t pJ = INT32_MAX;
    for (int32_t i = 0; i < stList_length(objects); i++) {
        int32_t j = lengthFn(stList_get(objects, i));
        assert(j <= pJ);
        pJ = j;
        totalLength += j;
        if (totalLength >= genomeLength / 2) {
            return j;
        }
    }
    return -1; //this should not happen!
}

stList *getScaffoldPathsList(stList *maximalHaplotypePaths, stList *eventStrings, CapCodeParameters *capCodeParameters) {
    stHash *scaffoldPaths = getScaffoldPaths(maximalHaplotypePaths, eventStrings, capCodeParameters);
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
    for (int32_t i = 0; i < stList_length(list); i+=2) {
        assert(stIntTuple_getPosition(stList_get(list, i), 0) == stIntTuple_getPosition(stList_get(list, i+1), 0));
        cAA[i/2] = stString_print("%i", stIntTuple_getPosition(stList_get(list, i), 0));
    }
    char *cA = stString_join(" ", (const char **)cAA, stList_length(list)/2);
    for (int32_t i = 0; i < stList_length(list)/2; i++) {
        free(cAA[i]);
    }
    free(cAA);
    return cA;
}

void reportHaplotypePathStats(Flower *flower, FILE *fileHandle,
        const char *assemblyEventString,
        stList *eventStrings, CapCodeParameters *capCodeParameters) {
    /*
     * Gets stats on the maximal haplotype paths.
     */

    stList *maximalHaplotypePaths = getContigPaths(flower, assemblyEventString, eventStrings);
    maximalHaplotypePathToLength = buildContigPathToContigPathLengthHash(maximalHaplotypePaths);
    maximalScaffoldPathToLength = getMaximalScaffoldPathLengths(maximalHaplotypePaths, eventStrings, capCodeParameters);

    //initialise the global arrays
    insertionDistribution = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    deletionDistribution = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);

    for (int32_t i = 0; i < stList_length(maximalHaplotypePaths); i++) {
        stList *maximalHaplotypePath = stList_get(maximalHaplotypePaths, i);
        totalHaplotypeLength += contigPathLength(maximalHaplotypePath);
        for (int32_t j = 0; j < stList_length(maximalHaplotypePath); j++) {
            Segment *segment = stList_get(maximalHaplotypePath, j);
            reportHaplotypePathStatsP(segment_get5Cap(segment), eventStrings, capCodeParameters);
            reportHaplotypePathStatsP(segment_get3Cap(segment), eventStrings, capCodeParameters);
        }
    }

    blockList = stList_construct();
    contigsSet = stSortedSet_construct3(compareSequences, NULL);
    haplotypesSet = stSortedSet_construct3(compareSequences, NULL);

    traverseBlocks(flower, assemblyEventString, eventStrings);

    stList *sequences = stSortedSet_getList(contigsSet);
    stList *haplotypes = stSortedSet_getList(haplotypesSet);
    stList_sort(blockList, compareBlocksByLength);
    stList_sort(sequences, compareSequencesByLength);
    stList_sort(maximalHaplotypePaths, compareMaximalHaplotypePaths);
    stList *scaffoldPaths = getScaffoldPathsList(maximalHaplotypePaths, eventStrings, capCodeParameters);
    stList_sort(scaffoldPaths, compareScaffoldPaths);

    int32_t totalSequencesLength = 0;
    for (int32_t i = 0; i < stList_length(sequences); i++) {
        totalSequencesLength += sequence_getLength(stList_get(sequences, i));
    }

    int32_t averageHaplotypeLength = 0; //average length of the two haplotypes.
    for (int32_t i = 0; i < stList_length(haplotypes); i++) {
        averageHaplotypeLength += sequence_getLength(stList_get(haplotypes, i));
    }
    averageHaplotypeLength /= 2;

    int32_t totalBlockNumber = stList_length(blockList);
    int32_t blockNG50 = getN50(averageHaplotypeLength, blockList, (int32_t(*)(const void *)) block_getLength);
    int32_t contigN50 = getN50(totalSequencesLength, sequences, (int32_t(*)(const void *)) sequence_getLength); //length of contig which appears in order (from longest to shortest) at 50% coverage.
    int32_t contigNG50 = getN50(averageHaplotypeLength, sequences, (int32_t(*)(const void *)) sequence_getLength); //length of contig which appears in order (from longest to shortest) at 50% coverage.
    int32_t haplotypePathNG50 = getN50(averageHaplotypeLength, maximalHaplotypePaths, getHaplotypePathLength); //length of maximal haplotype path which appears in order (from longest to shortest) at 50% coverage.
    int32_t scaffoldPathNG50 = getN50(averageHaplotypeLength, scaffoldPaths, getScaffoldPathLength); //length of maximal scaffold path which appears in order (from longest to shortest) at 50% coverage.

    int32_t totalContigNumber = stSortedSet_size(contigsSet); //number of contigs
    int32_t totalHaplotypePaths = stList_length(maximalHaplotypePaths); //number of haplotype paths
    int32_t totalScaffoldPaths = stList_length(scaffoldPaths); //number of scaffold paths

    int32_t totalErrors = totalErrorsHapToHapSameChromosome / 2 + totalErrorsHapToHapDifferentChromosome / 2
            + totalErrorsHapToContamination + totalErrorsHapToInsertToContamination + totalErrorsHapToInsert / 2
            + totalErrorsHapToDeletion / 2 + totalErrorsHapToInsertionAndDeletion / 2 + totalErrorsContigEndsWithInsert;

    double errorsPerContig = ((double) totalErrors) / totalContigNumber; //number of errors / number of contigs
    double errorsPerMappedBase = ((double) totalErrors) / totalPathLength; //number of errors / total path length
    double coverage = ((double) totalPathLength) / averageHaplotypeLength; //totalPathLength/genotypeLength


    assert(totalErrorsHapToHapSameChromosome % 2 == 0);
    assert(totalErrorsHapToHapDifferentChromosome % 2 == 0);
    assert(totalErrorsHapToInsert % 2 == 0);
    assert(totalErrorsHapToDeletion % 2 == 0);
    assert(totalErrorsHapToInsertionAndDeletion % 2 == 0);
    assert(totalScaffoldGaps % 2 == 0);
    assert(totalAmbiguityGaps % 2 == 0);

    char *insertionDistributionString = concatenateList(insertionDistribution);
    char *deletionDistributionString = concatenateList(deletionDistribution);

    fprintf(fileHandle, "<stats totalHaplotypeSwitches=\"%i\" "
        "totalScaffoldGaps=\"%i\" "
        "totalAmbiguityGaps=\"%i\" "
        "totalContigEnds=\"%i\" "
        "totalContigEndsWithNs=\"%i\" "
        "totalErrorsHaplotypeToHaplotypeSameChromosome=\"%i\" "
        "totalErrorsHaplotypeToHaplotypeDifferentChromosome=\"%i\" "
        "totalErrorsHaplotypeToContamination=\"%i\" "
        "totalErrorsHaplotypeToInsertionToContamination=\"%i\" "
        "totalErrorsHaplotypeToInsertion=\"%i\" "
        "totalErrorsHaplotypeToDeletion=\"%i\" "
        "totalErrorsHaplotypeToInsertionAndDeletion=\"%i\" "
        "totalErrorsContigEndsWithInsert=\"%i\" "
        "totalErrors=\"%i\" totalPathLength=\"%i\" genotypeLength=\"%i\" "
        "totalContigsLength=\"%i\" coverage=\"%f\" blockNG50=\"%i\" contigN50=\"%i\" "
        "contigNG50=\"%i\" contigPathNG50=\"%i\" scaffoldPathNG50=\"%i\" totalBlockNumber=\"%i\" "
        "totalContigNumber=\"%i\" totalHaplotypePaths=\"%i\" totalScaffoldPaths=\"%i\" "
        "errorsPerContig=\"%f\" errorsPerMappedBase=\"%f\" "
        "insertionErrorSizeDistribution=\"%s\" "
        "deletionErrorSizeDistribution=\"%s\"/>", totalHapSwitches / 2,
            totalScaffoldGaps / 2, totalAmbiguityGaps / 2, totalContigEnds,
            totalContigEndsWithsNs, totalErrorsHapToHapSameChromosome / 2, totalErrorsHapToHapDifferentChromosome / 2,
            totalErrorsHapToContamination, totalErrorsHapToInsertToContamination, totalErrorsHapToInsert / 2,
            totalErrorsHapToDeletion / 2, totalErrorsHapToInsertionAndDeletion / 2, totalErrorsContigEndsWithInsert,
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
    stList *eventStrings = getEventStrings(hap1EventString, hap2EventString);
    reportHaplotypePathStats(flower, fileHandle, assemblyEventString, eventStrings, capCodeParameters);
    fclose(fileHandle);
    st_logInfo("Got the stats in %i seconds/\n", time(NULL) - startTime);

    return 0;
}
