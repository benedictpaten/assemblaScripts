include ${rootPath}/include.mk
haplotypeSequences=${dataPath}/assemblathon1/simulatedHaplotypes/testHap1.fa.complete ${dataPath}/assemblathon1/simulatedHaplotypes/testHap2.fa.complete ${dataPath}/assemblathon1/simulatedHaplotypes/testHap3.fa.complete
newickTree=(((hapA1:0.002,hapA2:0.002):0.0000001,ecoli:0.002):0.000001,assembly:0.002);
hap1EventString=hapA1
hap2EventString=hapA2
assemblyEventString=assembly
contaminationEventString=ecoli
featureBedFiles=${dataPath}/tests/little/beds/test1.bed ${dataPath}/tests/little/beds/test2.bed
geneBedFiles=${dataPath}/tests/little/beds/test1.bed