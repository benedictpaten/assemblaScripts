rootPath=../../git/cactus/
include ${rootPath}include.mk
binPath=../bin

libSources = impl/*.c
libHeaders = inc/*.h
commonSources = impl/assemblaCommon.c

all : ${binPath}/coveragePlots ${binPath}/substitutionStats ${binPath}/pathAnnotatedMafGenerator ${binPath}/pathStats ${binPath}/copyNumberStats ${binPath}/linkageStats  

${binPath}/coveragePlots :  ${libSources} ${libHeaders} ${libPath}/cactusTreeStats.a  ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -I inc -o ${binPath}/coveragePlots impl/coveragePlots.c ${commonSources} ${libPath}/cactusTreeStats.a ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibs}

${binPath}/substitutionStats : ${libSources} ${libHeaders} ${libPath}/cactusMafs.a ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -I inc -o ${binPath}/substitutionStats impl/substitutionStats.c ${commonSources} ${libPath}/cactusMafs.a ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibs}

${binPath}/copyNumberStats : ${libSources} ${libHeaders} ${libPath}/cactusMafs.a ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -I inc -o ${binPath}/copyNumberStats impl/copyNumberStats.c ${commonSources} ${libPath}/cactusMafs.a ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibs}

${binPath}/linkageStats : ${libSources} ${libHeaders} ${libPath}/cactusMafs.a ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -I inc -o ${binPath}/linkageStats impl/linkageStats.c ${commonSources} ${libPath}/cactusMafs.a ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibs}

${binPath}/pathStats : ${libSources} ${libHeaders} ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -I inc -o ${binPath}/pathStats impl/pathStats.c ${commonSources} ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibs}

${binPath}/pathAnnotatedMafGenerator : ${libSources} ${libHeaders} ${libPath}/cactusMafs.a ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibsDependencies}
	${cxx} ${cflags} -I${libPath} -I inc -o ${binPath}/pathAnnotatedMafGenerator impl/pathAnnotatedMAFGenerator.c ${commonSources} ${libPath}/cactusMafs.a ${libPath}/cactusLib.a ${libPath}/assemblaLib.a ${basicLibs}

clean : 
	rm -f *.o
	rm -f ${binPath}/coveragePlots ${binPath}/pathStats ${binPath}/substitutionStats ${binPath}/pathAnnotatedMafGenerator ${binPath}/copyNumberStats ${binPath}/linkageStats