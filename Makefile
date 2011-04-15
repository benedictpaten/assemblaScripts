.PHONY: all clean cleanTest test testBig testLittle testStatic run scaffolds contigs static

all : 
	cd src && make all

clean : 
	cd src && make clean
	cd tests/big && make clean
	cd tests/little && make clean
	cd tests/static && make clean
	cd assemblathon1/scaffolds && make clean
	cd assemblathon1/contigs && make clean
	cd assemblathon1/static && make clean

cleanTest :
	cd tests/big && make clean
	cd tests/little && make clean
	cd tests/static && make clean

test : testBig testLittle testLittleStatic 

testBig :
	cd tests/big && make all

testLittle :
	cd tests/little && make all

testStatic :
	cd tests/static && make all

run : scaffolds contigs static

scaffolds :
	cd assemblathon1/scaffolds && make all

contigs : 
	cd assemblathon1/contigs && make all

static :
	cd assemblathon1/static && make all
