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

#include "sonLib.h"
#include "cactus.h"
#include "adjacencyClassification.h"
#include "assemblaCommon.h"

/*
 * Global parameters.
 */

CapCodeParameters *capCodeParameters = NULL;
char *outputFile = NULL;
Flower *flower = NULL;
CactusDisk *cactusDisk = NULL;

char *assemblyEventString = NULL;
char *hap1EventString = NULL;
char *hap2EventString = NULL;
char *contaminationEventString = NULL;

/*
 * Optional parameter used by copy number and substitution scripts.
 */
int32_t minimumBlockLength = 0;

/*
 * Parameters for the substitution script.
 */
int32_t ignoreFirstNBasesOfBlock = 0;
int32_t minimumIndentity = 0;
bool printIndelPositions = 0;

/*
 * For the linkage script.
 */
int32_t bucketNumber = 2000;
int32_t upperLinkageBound = 200000000;
int32_t sampleNumber = 1000000;

stList *getEventStrings(const char *hapA1EventString,
        const char *hapA2EventString) {
    stList *eventStrings = stList_construct3(0, NULL);
    stList_append(eventStrings, stString_copy(hapA1EventString));
    stList_append(eventStrings, stString_copy(hapA2EventString));
    return eventStrings;
}

void basicUsage(const char *programName) {
    fprintf(stderr, "%s\n", programName);
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-c --cactusDisk : The location of the cactus disk\n");
    fprintf(stderr, "-e --outputFile : The file to write the output in.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
    fprintf(
            stderr,
            "-m --minimumNsForScaffoldGap : Minimum number of Ns in assembly sequence to denote a scaffold gap\n");
    fprintf(stderr, "-n --maximumDeletionLength : Maximum length of deletion\n");
    fprintf(stderr,
            "-o --maximumInsertionLength : Maximum length of insertion\n");

    fprintf(stderr, "-p --assemblyEventString : The assembly event string\n");
    fprintf(stderr,
            "-q --haplotype1EventString : The haplotype 1 event string\n");
    fprintf(stderr,
            "-r --haplotype2EventString : The haplotype 2 event string\n");
    fprintf(stderr,
            "-s --contaminationEventString : The contamination event string\n");

    //Parameters not used by everyone..
    fprintf(stderr, "The following are not used by all scripts\n");
    fprintf(stderr, "-t --minimumBlockLength : Minimum block length\n");
    fprintf(stderr, "-u --ignoreFirstNBasesOfBlock : Minimum block length\n");
    fprintf(stderr, "-v --minimumIdentity : Minimum identity of the block\n");
    fprintf(stderr,
            "-w --printIndelPositions : Print out valid columns containing only one haplotype\n");
    fprintf(stderr, "-x --bucketNumber : Number of buckets\n");
    fprintf(stderr, "-y --upperLinkageBound : Upper linkage bound\n");
    fprintf(stderr, "-z --sampleNumber : Number of samples\n");
}

int parseBasicArguments(int argc, char *argv[], const char *programName) {
    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    int32_t k;
    capCodeParameters = capCodeParameters_construct(25, INT32_MAX, 100000);
    assemblyEventString = NULL;
    hap1EventString = NULL;
    hap2EventString = NULL;
    contaminationEventString = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument,
                0, 'c' }, { "outputFile", required_argument, 0, 'e' }, {
                "help", no_argument, 0, 'h' }, { "minimumNsForScaffoldGap",
                required_argument, 0, 'm' }, { "maximumDeletionLength",
                required_argument, 0, 'n' }, { "maximumInsertionLength",
                required_argument, 0, 'o' }, { "assemblyEventString",
                required_argument, 0, 'p' }, { "haplotype1EventString",
                required_argument, 0, 'q' }, { "haplotype2EventString",
                required_argument, 0, 'r' }, { "contaminationEventString",
                required_argument, 0, 's' }, { "minimumBlockLength",
                required_argument, 0, 't' }, { "ignoreFirstNBasesOfBlock",
                required_argument, 0, 'u' }, { "minimumIdentity",
                required_argument, 0, 'v' }, { "printIndelPositions",
                no_argument, 0, 'w' }, { "bucketNumber", required_argument, 0,
                'x' }, { "upperLinkageBound", required_argument, 0, 'y' }, {
                "sampleNumber", required_argument, 0, 'z' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:e:hm:n:o:p:q:r:s:t:u:v:wx:y:z:",
                long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'e':
                outputFile = stString_copy(optarg);
                break;
            case 'h':
                basicUsage(programName);
                return 0;
            case 'm':
                k = sscanf(optarg, "%i", &capCodeParameters->minimumNCount);
                assert(k == 1);
                break;
            case 'n':
                k = sscanf(optarg, "%i", &capCodeParameters->maxDeletionLength);
                assert(k == 1);
                break;
            case 'o':
                k
                        = sscanf(optarg, "%i",
                                &capCodeParameters->maxInsertionLength);
                assert(k == 1);
                break;
            case 'p':
                assemblyEventString = stString_copy(optarg);
                break;
            case 'q':
                hap1EventString = stString_copy(optarg);
                break;
            case 'r':
                hap2EventString = stString_copy(optarg);
                break;
            case 's':
                contaminationEventString = stString_copy(optarg);
                break;
            case 't':
                k = sscanf(optarg, "%i", &minimumBlockLength);
                assert(k == 1);
                break;
            case 'u':
                k = sscanf(optarg, "%i", &ignoreFirstNBasesOfBlock);
                assert(k == 1);
                break;
            case 'v':
                k = sscanf(optarg, "%i", &minimumIndentity);
                assert(k == 1);
                if (minimumIndentity > 100 || minimumIndentity < 0) {
                    st_errAbort(
                            "The minimum identity was not in the range [0, 100]: %i",
                            minimumIndentity);
                }
                break;
            case 'w':
                printIndelPositions = 1;
                break;
            case 'x':
                k = sscanf(optarg, "%i", &bucketNumber);
                assert(k == 1);
                if (bucketNumber < 1) {
                    st_errAbort(
                            "The number of buckets can not be less than 1: %i",
                            bucketNumber);
                }
                break;
            case 'y':
                k = sscanf(optarg, "%i", &upperLinkageBound);
                assert(k == 1);
                if (upperLinkageBound < 1) {
                    st_errAbort(
                            "The total interval of linkage is too small: %i",
                            upperLinkageBound);
                }
                break;
            case 'z':
                k = sscanf(optarg, "%i", &sampleNumber);
                assert(k == 1);
                if (sampleNumber < 0) {
                    st_errAbort(
                            "The number of samples can not be less than 0: %i",
                            sampleNumber);
                }
                break;
            default:
                st_errAbort("Unrecognised option %s", optarg);
                break;
        }
    }

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if (logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if (logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    if (outputFile == NULL) {
        st_errAbort("The output file was not specified");
    }
    if (cactusDiskDatabaseString == NULL) {
        st_errAbort("The cactus disk string was not specified");
    }
    if (assemblyEventString == NULL) {
        st_errAbort("The assembly event string was not specified");
    }
    if (hap1EventString == NULL) {
        st_errAbort("The haplotype 1 event string was not specified");
    }
    if (hap2EventString == NULL) {
        st_errAbort("The haplotype 2 event string was not specified");
    }
    if (contaminationEventString == NULL) {
        st_errAbort("The contamination event string was not specified");
    }

    assert(outputFile != NULL);
    assert(cactusDiskDatabaseString != NULL);

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Output graph file : %s\n", outputFile);
    st_logInfo("The cactus disk string : %s\n", cactusDiskDatabaseString);
    st_logInfo("The assembly event string : %s\n", assemblyEventString);
    st_logInfo("The haplotype 1 event string : %s\n", hap1EventString);
    st_logInfo("The haplotype 2 event string : %s\n", hap2EventString);
    st_logInfo("The contamination event string : %s\n",
            contaminationEventString);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(
            cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the cactus disk\n");

    //////////////////////////////////////////////
    //Cleanup
    //////////////////////////////////////////////

    free(cactusDiskDatabaseString);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    flower = cactusDisk_getFlower(cactusDisk, 0);
    assert(flower != NULL);
    st_logInfo("Parsed the top level flower of the cactus tree\n");

    return 0;
}

