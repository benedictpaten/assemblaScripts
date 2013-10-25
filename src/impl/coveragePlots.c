/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten (at) gmail.com) and Dent Earl (dearl (at) soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"
#include "treeStats.h"

#include "contigPaths.h"
#include "adjacencyTraversal.h"
#include "adjacencyClassification.h"
#include "scaffoldPaths.h"
#include "assemblaCommon.h"

/*
 * For a range of block, contig and contig-path length values reports
 * presence and absence from different classes of alignment column.
 */

typedef struct _blockHolder {
    Block *block;
    int64_t blockLength;
    int64_t haplotypePathLength;
    int64_t scaffoldPathLength;
    int64_t contigLength;

    /*
     * Categories as follows:
     *
     * hap1/hap2/assembly = 0
     * hap1/hap2/!assembly = 1
     * hap1/!hap2/assembly = 2
     * hap1/!hap2/!assembly = 3
     * !hap1/hap2/assembly = 4
     * !hap1/hap2/!assembly = 5
     * !hap1/!hap2/assembly = 6
     */
    int64_t haplotypeCategory;
    /*
     * Categories as follows:
     *
     * contamination/assembly = 0
     * contamination/!assembly = 1
     * !contamination/assembly = 2
     */
    int64_t contaminationCategory;
    /*
     * Categories as follows:
     *
     * contamination/hap = 0
     * contamination/!hap = 1
     * !contamination/hap = 2
     */
    int64_t haplotypeToContaminationCategory;
} BlockHolder;

stHash *segmentsToMaximalHaplotypePaths;
stHash *maximalHaplotypePathLengths;
stHash *maximalScaffoldPathLengths;

int64_t getMaximalHaplotypePathLengthP(Block *block,
        stHash *segmentToHaplotypePath, stHash *haplotypePathLengths) {
    Segment *segment;
    Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
    int64_t maxLength = 0;
    while ((segment = block_getNext(instanceIt)) != NULL) {
        if (strcmp(event_getHeader(segment_getEvent(segment)), assemblyEventString) == 0) { //Establish if we need a line..
            stList *maximalHaplotypePath = stHash_search(
                    segmentToHaplotypePath, segment);
            if (maximalHaplotypePath == NULL) {
                maximalHaplotypePath = stHash_search(segmentToHaplotypePath,
                        segment_getReverse(segment));
            }
            if (maximalHaplotypePath != NULL) {
                assert(
                        stHash_search(haplotypePathLengths,
                                maximalHaplotypePath) != NULL);
                int64_t i = stIntTuple_get(stHash_search(
                        haplotypePathLengths, maximalHaplotypePath), 0);
                if (i > maxLength) {
                    maxLength = i;
                }
            }
        }
    }
    block_destructInstanceIterator(instanceIt);

    return maxLength;
}

int64_t getMaximalHaplotypePathLength(Block *block) {
    return getMaximalHaplotypePathLengthP(block,
            segmentsToMaximalHaplotypePaths, maximalHaplotypePathLengths);
}

int64_t getMaximalScaffoldPathLength(Block *block) {
    return getMaximalHaplotypePathLengthP(block,
            segmentsToMaximalHaplotypePaths, maximalScaffoldPathLengths);
}

int64_t getMaximalContigLength(Block *block) {
    Segment *segment;
    Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
    int64_t maxLength = 0;
    while ((segment = block_getNext(instanceIt)) != NULL) {
        if (strcmp(event_getHeader(segment_getEvent(segment)), assemblyEventString) == 0) { //Establish if we need a line..
            Sequence *sequence = segment_getSequence(segment);
            assert(sequence != NULL);
            if (sequence_getLength(sequence) > maxLength) {
                maxLength = sequence_getLength(sequence);
            }
        }
    }
    block_destructInstanceIterator(instanceIt);

    return maxLength;
}

stSortedSet *getSpecies(Block *block) {
    stSortedSet *species = stSortedSet_construct3((int(*)(const void *,
            const void *)) strcmp, free);
    Segment *segment;
    Block_InstanceIterator *segmentIterator = block_getInstanceIterator(block);
    while ((segment = block_getNext(segmentIterator)) != NULL) {
        //Event name
        Event *event = segment_getEvent(segment);
        assert(event != NULL);
        const char *eventHeader = event_getHeader(event);
        if (stSortedSet_search(species, (void *) eventHeader) == NULL) {
            stSortedSet_insert(species, stString_copy(eventHeader));
        }
    }
    block_destructInstanceIterator(segmentIterator);
    return species;
}

int64_t getHaplotypeCategory(stSortedSet *species) {
    int64_t i = stSortedSet_search(species, assemblyEventString) == NULL ? 1 : 0;
    int64_t j = stSortedSet_search(species, hap1EventString) == NULL ? 2 : 0;
    int64_t k = stSortedSet_search(species, hap2EventString) == NULL ? 4 : 0;
    return i + j + k;
}

int64_t getContaminationCategory(stSortedSet *species) {
    int64_t i = stSortedSet_search(species, assemblyEventString) == NULL ? 1 : 0;
    int64_t j = stSortedSet_search(species, contaminationEventString) == NULL ? 2 : 0;
    return i + j;
}

int64_t getHaplotypeContaminationCategory(stSortedSet *species) {
    int64_t i = stSortedSet_search(species, contaminationEventString) == NULL ? 1 : 0;
    int64_t j = (stSortedSet_search(species, hap1EventString) == NULL
            && stSortedSet_search(species, hap2EventString) == NULL) ? 2 : 0;
    return i + j;
}

static BlockHolder *blockHolder_construct(Block *block) {
    BlockHolder *blockHolder = st_malloc(sizeof(BlockHolder));

    blockHolder->block = block;

    //get block length
    blockHolder->blockLength = block_getLength(block);
    assert(blockHolder->blockLength >= 0);

    //get contig length
    blockHolder->contigLength = getMaximalContigLength(block);
    assert(blockHolder->contigLength >= 0);

    //get haplotype path length
    blockHolder->haplotypePathLength = getMaximalHaplotypePathLength(block);
    assert(blockHolder->haplotypePathLength >= 0);

    //get haplotype path length
    blockHolder->scaffoldPathLength = getMaximalScaffoldPathLength(block);
    assert(blockHolder->scaffoldPathLength >= 0);
    assert(blockHolder->scaffoldPathLength >= blockHolder->haplotypePathLength);

    //get haplotype category
    stSortedSet *species = getSpecies(block);
    blockHolder->haplotypeCategory = getHaplotypeCategory(species);
    assert(blockHolder->haplotypeCategory >= 0);
    assert(blockHolder->haplotypeCategory < 8);

    //get contamination category
    blockHolder->contaminationCategory = getContaminationCategory(species);
    assert(blockHolder->contaminationCategory >= 0);
    assert(blockHolder->contaminationCategory < 4);

    //get haplotype contamination category
    blockHolder->haplotypeToContaminationCategory = getHaplotypeContaminationCategory(species);
    assert(blockHolder->haplotypeToContaminationCategory >= 0);
    assert(blockHolder->haplotypeToContaminationCategory < 4);

    stSortedSet_destruct(species);

    return blockHolder;
}

static void blockHolder_destruct(BlockHolder *blockHolder) {
    free(blockHolder);
}

static int blockHolder_compareByBlockLength(const BlockHolder *blockHolder1,
        const BlockHolder *blockHolder2) {
    return blockHolder1->blockLength - blockHolder2->blockLength;
}

static int blockHolder_compareByHaplotypePathLength(
        const BlockHolder *blockHolder1, const BlockHolder *blockHolder2) {
    return blockHolder1->haplotypePathLength
            - blockHolder2->haplotypePathLength;
}

static int blockHolder_compareByScaffoldPathLength(
        const BlockHolder *blockHolder1, const BlockHolder *blockHolder2) {
    return blockHolder1->scaffoldPathLength - blockHolder2->scaffoldPathLength;
}

static int blockHolder_compareByContigLength(const BlockHolder *blockHolder1,
        const BlockHolder *blockHolder2) {
    return blockHolder1->contigLength - blockHolder2->contigLength;
}

static int64_t blockHolder_getBlockLength(const BlockHolder *blockHolder) {
    return blockHolder->blockLength;
}

static int64_t blockHolder_getHaplotypePathLength(
        const BlockHolder *blockHolder) {
    return blockHolder->haplotypePathLength;
}

static int64_t blockHolder_getScaffoldPathLength(const BlockHolder *blockHolder) {
    return blockHolder->scaffoldPathLength;
}

static int64_t blockHolder_getContigLength(const BlockHolder *blockHolder) {
    return blockHolder->contigLength;
}

static int64_t blockHolder_getHaplotypeCategory(const BlockHolder *blockHolder) {
    return blockHolder->haplotypeCategory;
}

static int64_t blockHolder_getContaminationCategory(const BlockHolder *blockHolder) {
    return blockHolder->contaminationCategory;
}

static int64_t blockHolder_getHaplotypeContaminationCategory(
        const BlockHolder *blockHolder) {
    return blockHolder->haplotypeToContaminationCategory;
}

static void printCumulativeLengthPlots(stList *blockHolders, int(*cmpFn)(
        const BlockHolder *, const BlockHolder *), int64_t(*getCategory)(
        const BlockHolder *), int64_t(*getLength)(const BlockHolder *),
        int64_t categoryNumber, const char **categoryNames,
        const char *outputFile) {
    FILE *fileHandle = fopen(outputFile, "w");

    stList_sort(blockHolders, (int(*)(const void *, const void *)) cmpFn);

    int64_t binNumber = 2000;
    double binSize = 8.0 / binNumber; //We go up to 100,000,000
    int64_t *cumulativeLengths = st_malloc(sizeof(int64_t) * categoryNumber
            * (binNumber + 1));
    for (int64_t i = 0; i < categoryNumber; i++) {
        cumulativeLengths[i] = 0;
    }
    for (int64_t i = 0; i < stList_length(blockHolders); i++) { //Get the start values
        BlockHolder *blockHolder = stList_get(blockHolders, i);
        cumulativeLengths[getCategory(blockHolder)]
                += blockHolder_getBlockLength(blockHolder);
    }

    int64_t j = 0; //index of bin
    //Get cumulative values for the remaining bins
    for (int64_t i = 1; i <= binNumber; i++) {
        //Copy across the cumulative values from the previous bin
        for (int64_t k = 0; k < categoryNumber; k++) {
            cumulativeLengths[i * categoryNumber + k] = cumulativeLengths[(i
                    - 1) * categoryNumber + k];
        }
        //Get minimum length for block in bin
        double minimumBlockLength = pow(10, i * binSize);
        //Discard blocks smaller than threshold
        while (j < stList_length(blockHolders)) {
            BlockHolder *blockHolder = stList_get(blockHolders, j);
            if (getLength(blockHolder) >= minimumBlockLength) {
                break;
            } else {
                j++;
            }
            if (getCategory(blockHolder) % 2 == 0) {
                cumulativeLengths[i * categoryNumber
                        + getCategory(blockHolder)]
                        -= blockHolder_getBlockLength(blockHolder);
                if (getCategory(blockHolder) < categoryNumber) {
                    //This shifts the numbers into the other column
                    cumulativeLengths[i * categoryNumber + getCategory(
                            blockHolder) + 1] += blockHolder_getBlockLength(
                            blockHolder);
                }
            }
        }
    }

    //Now print the final values..
    fprintf(fileHandle, "category\t");
    for (int64_t i = 0; i <= binNumber; i++) {
        fprintf(fileHandle, "%" PRIi64 "\t", (int64_t) pow(10, i * binSize));
    }
    fprintf(fileHandle, "\n");
    for (j = 0; j < categoryNumber; j++) {
        fprintf(fileHandle, "%s\t", categoryNames[j]);
        for (int64_t i = 0; i <= binNumber; i++) {
            fprintf(fileHandle, "%lli\t", (long long int) cumulativeLengths[i
                    * categoryNumber + j]);
        }
        fprintf(fileHandle, "\n");
    }
    free(cumulativeLengths);

    fclose(fileHandle);
}

void getBlocks(Flower *flower, stList *blockHolders) {
    Flower_BlockIterator *blockIt = flower_getBlockIterator(flower);
    Block *block;
    while ((block = flower_getNextBlock(blockIt)) != NULL) {
        if (block_getInstanceNumber(block) > 0) {
            stList_append(blockHolders, blockHolder_construct(block));
        }
    }
    flower_destructBlockIterator(blockIt);
    //Now recurse on the groups
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (group_getNestedFlower(group) != NULL) {
            getBlocks(group_getNestedFlower(group), blockHolders);
        }
    }
    flower_destructGroupIterator(groupIt);
}

int main(int argc, char *argv[]) {
    //////////////////////////////////////////////
    //Parse the inputs
    //////////////////////////////////////////////

    parseBasicArguments(argc, argv, "coveragePlots");

    ///////////////////////////////////////////////////////////////////////////
    // Calculate haplotype paths
    ///////////////////////////////////////////////////////////////////////////

    stList *haplotypeEventStrings = getEventStrings(hap1EventString, hap2EventString);
    stList *contaminationEventStrings = getEventStrings(contaminationEventString, NULL);

    stList *maximalHaplotypePaths = getContigPaths(flower, assemblyEventString, haplotypeEventStrings);
    segmentsToMaximalHaplotypePaths = buildSegmentToContigPathHash(
            maximalHaplotypePaths);
    maximalHaplotypePathLengths
            = buildContigPathToContigPathLengthHash(
                    maximalHaplotypePaths);
    maximalScaffoldPathLengths = getContigPathToScaffoldPathLengthsHash(
            maximalHaplotypePaths, haplotypeEventStrings, contaminationEventStrings, capCodeParameters);

    ///////////////////////////////////////////////////////////////////////////
    // Calculate blocks
    ///////////////////////////////////////////////////////////////////////////

    stList *blockHolders = stList_construct3(0,
            (void(*)(void *)) blockHolder_destruct);
    getBlocks(flower, blockHolders);

    const char *haplotypeCategoryNames[8] = { "hap1/hap2/assembly",
            "hap1/hap2/!assembly", "hap1/!hap2/assembly",
            "hap1/!hap2/!assembly", "!hap1/hap2/assembly",
            "!hap1/hap2/!assembly", "!hap1/!hap2/assembly", "all" };

    st_system("mkdir %s", outputFile);

    printCumulativeLengthPlots(blockHolders, blockHolder_compareByBlockLength,
            blockHolder_getHaplotypeCategory, blockHolder_getBlockLength, 8,
            haplotypeCategoryNames, stString_print(
                    "%s/blockLengthsVsCoverageOfAssemblyAndHaplotypes.txt", outputFile));
    printCumulativeLengthPlots(blockHolders,
            blockHolder_compareByHaplotypePathLength,
            blockHolder_getHaplotypeCategory,
            blockHolder_getHaplotypePathLength, 8, haplotypeCategoryNames,
            stString_print("%s/contigPathLengthsVsCoverageOfAssemblyAndHaplotypes.txt", outputFile));
    printCumulativeLengthPlots(blockHolders,
            blockHolder_compareByScaffoldPathLength,
            blockHolder_getHaplotypeCategory,
            blockHolder_getScaffoldPathLength, 8, haplotypeCategoryNames,
            stString_print("%s/scaffoldPathLengthsVsCoverageOfAssemblyAndHaplotypes.txt", outputFile));
    printCumulativeLengthPlots(blockHolders, blockHolder_compareByContigLength,
            blockHolder_getHaplotypeCategory, blockHolder_getContigLength, 8,
            haplotypeCategoryNames, stString_print(
                    "%s/contigLengthsVsCoverageOfAssemblyAndHaplotypes.txt", outputFile));

    const char *contaminationCategoryNames[4] = { "contamination/assembly", "contamination/!assembly",
            "!contamination/assembly", "all" };

    printCumulativeLengthPlots(blockHolders, blockHolder_compareByBlockLength,
            blockHolder_getContaminationCategory, blockHolder_getBlockLength, 4,
            contaminationCategoryNames, stString_print(
                    "%s/blockLengthsVsCoverageOfAssemblyAndContamination.txt", outputFile));
    printCumulativeLengthPlots(blockHolders, blockHolder_compareByContigLength,
            blockHolder_getContaminationCategory, blockHolder_getContigLength, 4,
            contaminationCategoryNames, stString_print(
                    "%s/contigLengthsVsCoverageOfAssemblyAndContamination.txt", outputFile));

    const char *contaminationHaplotypeCategoryNames[4] = { "hap/contamination", "hap/!contamination",
            "!hap/contamination", "all" };

    printCumulativeLengthPlots(blockHolders, blockHolder_compareByBlockLength,
            blockHolder_getHaplotypeContaminationCategory, blockHolder_getBlockLength,
            4, contaminationHaplotypeCategoryNames, stString_print(
                    "%s/blockLengthsVsCoverageOfHaplotypesAndContamination.txt", outputFile));
    printCumulativeLengthPlots(blockHolders, blockHolder_compareByContigLength,
            blockHolder_getHaplotypeContaminationCategory, blockHolder_getContigLength,
            4, contaminationHaplotypeCategoryNames, stString_print(
                    "%s/contigLengthsVsCoverageOfHaplotypesAndContamination.txt", outputFile));

    return 0;
}

