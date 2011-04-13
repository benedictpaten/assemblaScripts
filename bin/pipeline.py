#!/usr/bin/env python


"""Script runs cactus to compare a bunch of assemblies against a set of two haplotypes and 
a set of contamination sequences.
"""

import os
import xml.etree.ElementTree as ET
import xml
import sys
from optparse import OptionParser

from jobTree.scriptTree.target import Target 
from jobTree.scriptTree.stack import Stack

from jobTree.src.bioio import logger
from jobTree.src.bioio import setLoggingFromOptions

from cactus.shared.config import CactusWorkflowExperiment
from cactus.shared.common import runCactusWorkflow

from sonLib.bioio import getTempFile, getTempDirectory
from sonLib.bioio import fastaRead, fastaWrite
from sonLib.bioio import system

from jobTree.test.jobTree.jobTreeTest import runJobTreeStatusAndFailIfNotComplete

def getRootPathString():
    """
    function for finding external location
    """
    import os
    import assemblathon.bin.pipeline
    i = os.path.abspath(assemblathon.bin.pipeline.__file__)
    return os.path.split(os.path.split(i)[0])[0] #os.path.split(os.path.split(os.path.split(i)[0])[0])[0]

def getCactusDiskString(alignmentFile):
    return "<st_kv_database_conf type=\"tokyo_cabinet\"><tokyo_cabinet database_dir=\"%s\"/></st_kv_database_conf>" % alignmentFile

class MakeAlignment(Target):
    """Target runs the alignment.
    """
    def __init__(self, newickTree, haplotypeSequences, 
                 assemblyFile, outputDir, configFile, options):
        Target.__init__(self, cpu=1, memory=8000000000)
        self.newickTree = newickTree
        self.haplotypeSequences = haplotypeSequences
        self.assemblyFile = assemblyFile
        self.outputDir = outputDir
        self.configFile = configFile
        self.options = options
    
    def run(self):
        cactusAlignmentName = "cactusAlignent"
        cactusAlignment = os.path.join(self.outputDir, cactusAlignmentName)
        if not os.path.exists(cactusAlignment):
            #Prepare the assembly
            #First copy it.
            if self.assemblyFile[-3:] == '.gz':
               tempAssemblyFile = getTempFile(rootDir=self.getLocalTempDir(), suffix=".gz")
               system("cp %s %s" % (self.assemblyFile, tempAssemblyFile))
               system("gunzip %s" % tempAssemblyFile)
               tempAssemblyFile = tempAssemblyFile[:-3]
               assert os.path.exists(tempAssemblyFile)
            else:
                tempAssemblyFile = getTempFile(rootDir=self.getLocalTempDir(), suffix="")
                system("cp %s %s" % (self.assemblyFile, tempAssemblyFile))
            #Make the supporting temporary files
            tempExperimentFile = getTempFile(rootDir=self.getLocalTempDir())
            tempJobTreeDir = os.path.join(self.getLocalTempDir(), "jobTree")
            #Make the experiment file
            cactusWorkflowExperiment = CactusWorkflowExperiment(
                                                 sequences=self.haplotypeSequences + [ tempAssemblyFile ], 
                                                 newickTreeString=self.newickTree, 
                                                 databaseName=cactusAlignmentName,
                                                 outputDir=self.getLocalTempDir(),
                                                 configFile=self.configFile)
            cactusWorkflowExperiment.writeExperimentFile(tempExperimentFile)
            #Now run cactus workflow
            runCactusWorkflow(experimentFile=tempExperimentFile, jobTreeDir=tempJobTreeDir, 
                              setupAndBuildAlignments=True,
                              buildTrees=False, buildFaces=False, buildReference=True,
                              batchSystem="single_machine", maxThreads=1, jobTreeStats=True)
            logger.info("Ran the workflow")
            #Check if the jobtree completed sucessively.
            runJobTreeStatusAndFailIfNotComplete(tempJobTreeDir)
            logger.info("Checked the job tree dir")
            #Compute the stats
            cactusAlignmentDir = os.path.join(self.getLocalTempDir(), cactusAlignmentName)
            tempJobTreeStatsFile = os.path.join(self.getLocalTempDir(),"jobTreeStats.xml")
            system("jobTreeStats --jobTree %s --outputFile %s" % (tempJobTreeDir, tempJobTreeStatsFile))
            #Now copy the true assembly back to the output
            system("mv %s %s/config.xml" % (tempExperimentFile, self.outputDir))
            system("mv %s %s/" % (tempJobTreeStatsFile, self.outputDir))
            system("mv %s %s/" % (cactusAlignmentDir, self.outputDir))
            assert os.path.exists(cactusAlignment)
            #We're done!
        self.addChildTarget(MakeStats1(self.outputDir, cactusAlignment, self.options))
                
class MakeAlignments(Target):
    """Makes alignments using pipeline.
    """
    def __init__(self, newickTree, haplotypeSequences, 
                 assembliesDir, outputDir, configFile, options):
        Target.__init__(self)
        self.newickTree = newickTree
        self.haplotypeSequences = haplotypeSequences
        self.assembliesDir = assembliesDir
        self.outputDir = outputDir
        self.configFile = configFile
        self.options = options
    
    def run(self):
        for assembly in os.listdir(self.assembliesDir):
            if assembly[-3:] == '.gz' or assembly[-3:] == '.fa':
                assemblyFile = os.path.join(self.assembliesDir, assembly)
                #The output directory
                outputDirForAssembly = os.path.join(self.outputDir, assembly)
                if(outputDirForAssembly[-3:] == ".gz"):
                    outputDirForAssembly = outputDirForAssembly[:-3]
                #Make the output dir if it doesn't exist
                if not os.path.exists(outputDirForAssembly):
                    os.mkdir(outputDirForAssembly)
                #Make the output file
                self.addChildTarget(MakeAlignment(newickTree=self.newickTree, haplotypeSequences=self.haplotypeSequences,
                                                  assemblyFile=assemblyFile, outputDir=outputDirForAssembly, configFile=self.configFile, options=self.options))

class MakeStats1(Target):
    """Builds basic stats and the maf alignment(s).
    """
    def __init__(self, outputDir, alignment, options, cpu=4, memory=8000000000):
        Target.__init__(self, cpu=cpu, memory=memory)
        self.alignment = alignment
        self.options = options
        self.outputDir = outputDir
        
    def runScript(self, binaryName, outputFile, specialOptions):
        if not os.path.exists(outputFile):
            tempOutputFile = getTempFile(rootDir=self.getLocalTempDir())
            os.remove(tempOutputFile)
            system("%s --cactusDisk '%s' --outputFile %s --assemblyEventString %s \
--haplotype1EventString %s --haplotype2EventString %s \
--contaminationEventString %s --minimumNsForScaffoldGap %s %s" % 
            (os.path.join(getRootPathString(), "bin", binaryName),
             getCactusDiskString(self.alignment),
             tempOutputFile, 
             self.options.assemblyEventString,
             self.options.haplotype1EventString,
             self.options.haplotype2EventString,
             self.options.contaminationEventString,
             self.options.minimumNsForScaffoldGap, specialOptions))
            system("mv %s %s" % (tempOutputFile, outputFile))
        
    def run(self):
        outputFile = os.path.join(self.outputDir, "cactusTreeStats.xml")
        if not os.path.exists(outputFile):
            system("cactus_treeStats --cactusDisk '%s' --flowerName 0 --outputFile %s --noPerColumnStats" % (getCactusDiskString(self.alignment), outputFile))
        #outputFile = "%s.maf" % self.alignment
        #if not os.path.exists(outputFile):
        #    system("cactus_MAFGenerator --cactusDisk '%s' --flowerName 0 --outputFile %s --orderByReference" % (getCactusDiskString(self.alignment), outputFile))
        outputFile = os.path.join(self.outputDir, "annotatedPaths.maf")
        self.runScript("pathAnnotatedMafGenerator", outputFile, "")
        self.addChildTarget(MakeContigPathStats(self.outputDir, self.alignment, self.options))

class MakeContigPathStats(MakeStats1):
    """Makes contig path stats.
    """
    def run(self):
        outputFile = os.path.join(self.outputDir, "pathStats.xml")
        self.runScript("pathStats", outputFile, "")
        self.addChildTarget(MakeCoveragePlots(self.outputDir, self.alignment, self.options))
        
class MakeCoveragePlots(MakeStats1):
    """Makes coverage plots.
    """
    def run(self):
        self.runScript("coveragePlots", os.path.join(self.outputDir, "coveragePlots"), "")
        self.addChildTarget(MakeSubstitutionStats(self.outputDir, self.alignment, self.options))
        
class MakeSubstitutionStats(MakeStats1):
    """Makes substitution stats.
    """
    def run(self):
        outputFile = os.path.join(self.outputDir, "substitutionStats_1000_98_5.txt")
        self.runScript("substitutionStats", outputFile, "--ignoreFirstNBases 5 --minimumBlockLength 1000 --minimumIdentity 98")
        
        outputFile = os.path.join(self.outputDir, "substitutionStats_1000_98_5_indel_positions.txt")
        self.runScript("substitutionStats", outputFile, "--ignoreFirstNBases 5 --minimumBlockLength 1000 --minimumIdentity 98 --printIndelPositions")
        
        outputFile = os.path.join(self.outputDir, "substitutionStats_0_0_0.txt")
        self.runScript("substitutionStats", outputFile, "--ignoreFirstNBases 0 --minimumBlockLength 0")
        
        self.addChildTarget(MakeCopyNumberStats(self.outputDir, self.alignment, self.options))

class MakeCopyNumberStats(MakeStats1):
    """Make copy number stats.
    """
    def run(self):
        outputFile = os.path.join(self.outputDir, "copyNumberStats_0.xml")
        self.runScript("copyNumberStats", outputFile, "")
        
        outputFile = os.path.join(self.outputDir, "copyNumberStats_1000.xml")
        self.runScript("copyNumberStats", outputFile, "--minimumBlockLength 1000")
        
        self.addChildTarget(MakeLinkageStats(self.outputDir, self.alignment, self.options))
   
class MakeLinkageStats(MakeStats1):
    """Make linkage stats.
    """
    def run(self):
        outputFile = os.path.join(self.outputDir, "linkageStats.xml")
        self.runScript("linkageStats", outputFile, "--bucketNumber 2000 --sampleNumber 100000000")
    
def main():
    ##########################################
    #Construct the arguments.
    ##########################################
    
    parser = OptionParser()
 
    parser.add_option("--haplotypeSequences", dest="haplotypeSequences")
    parser.add_option("--newickTree", dest="newickTree")
    parser.add_option("--assembliesDir", dest="assembliesDir")
    parser.add_option("--outputDir", dest="outputDir")
    parser.add_option("--configFile", dest="configFile")
    parser.add_option("--minimumNsForScaffoldGap", dest="minimumNsForScaffoldGap")
    parser.add_option("--assemblyEventString", dest="assemblyEventString")
    parser.add_option("--haplotype1EventString", dest="haplotype1EventString")
    parser.add_option("--haplotype2EventString", dest="haplotype2EventString")
    parser.add_option("--contaminationEventString", dest="contaminationEventString")
    
    Stack.addJobTreeOptions(parser)
    
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    if len(args) != 0:
        raise RuntimeError("Unrecognised input arguments: %s" % " ".join(args))

    Stack(MakeAlignments(newickTree=options.newickTree, 
                         haplotypeSequences=options.haplotypeSequences.split(), 
                         assembliesDir=options.assembliesDir, 
                         outputDir=options.outputDir, 
                         configFile=options.configFile, 
                         options=options)).startJobTree(options)
    logger.info("Done with job tree")

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    from assemblathon.bin.pipeline import *
    _test()
    main()