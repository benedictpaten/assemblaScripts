include ${rootPath}/include.mk
haplotypeSequences=${dataPath}/assemblathon1/simulatedHaplotypes/hap1.trf.repmask.fa.complete  ${dataPath}/assemblathon1/simulatedHaplotypes/hap2.trf.repmask.fa.complete ${dataPath}/assemblathon1/simulatedHaplotypes/Escherichia_coliO83_H1_1.trf.repmask.fa.complete
newickTree=(((hapA1:0.002,hapA2:0.002):0.0000001,ecoli:0.002):0.000001,assembly:0.002);
hap1EventString=hapA1
hap2EventString=hapA2
assemblyEventString=assembly
contaminationEventString=ecoli
featureBedFiles=${dataPath}/assemblathon1/simulatedHaplotypes/bedsA1/* ${dataPath}/assemblathon1/simulatedHaplotypes/bedsA2/*
geneBedFiles=${dataPath}/assemblathon1/simulatedHaplotypes/*/*TRANSCRIPTS* ${dataPath}/assemblathon1/simulatedHaplotypes/*/*CDS* 
