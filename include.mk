binPath = ${rootPath}bin
libPath = ${rootPath}lib
dataPath = ${rootPath}data
outputPath = ${rootPath}output

#Shared Parameters
#jobTreeFlags = --batchSystem parasol --logDebug --retryCount 0 --maxThreads 4
jobTreeFlags = --batchSystem singleMachine --maxThreads 20 --logDebug --retryCount 0
configFile=${libPath}/config_fast.xml
minimumNsForScaffoldGap=15

#########
#Build basic cactus alignment
#########
		
pipeline :
	rm -rf ./jobTree
	#Running pipeline to build comparisons
	python ${binPath}/pipeline.py --assemblyEventString ${assemblyEventString} --haplotype1EventString ${hap1EventString} --haplotype2EventString ${hap2EventString} --contaminationEventString ${contaminationEventString} --haplotypeSequences '${haplotypeSequences}' --newickTree '${newickTree}' --assembliesDir ${assembliesDir} --outputDir ${outputDir} --configFile ${configFile} --minimumNsForScaffoldGap ${minimumNsForScaffoldGap} --jobTree ./jobTree ${jobTreeFlags}
	jobTreeStatus --jobTree ./jobTree --failIfNotComplete
	rm -rf ./jobTree
	
basic : pipeline
