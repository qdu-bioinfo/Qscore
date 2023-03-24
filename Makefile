
CC:=g++
ifneq (,$(findstring Darwin,$(shell uname)))
	exist = $(shell if [ -e '/usr/local/bin/g++-11' ]; then echo "exist"; else echo "notexist"; fi;)
	ifeq ($(exist),exist)
		CC:=g++-11
	else
        	exist = $(shell if [ -e '/usr/local/bin/g++-10' ]; then echo "exist"; else echo "notexist"; fi;)
        	ifeq ($(exist),exist)
                	CC:=g++-10
		else
			CC:=g++-9
		endif
	endif
endif

HASHFLG=-Wno-deprecated
BUILDFLG=-w -ffunction-sections -fdata-sections -fmodulo-sched -msse

EXE_CMP1=bin/Qscore
EXE_CMP2=bin/Extract_16S_rRNA
EXE_CMP3=bin/Extract_WGS

tax:$(OBJ_TAX) src/Qscore.cpp
	$(CC) -o $(EXE_CMP1) src/Qscore.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_CMP2) src/Extract_16S_rRNA.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)
	$(CC) -o $(EXE_CMP3) src/Extract_WGS.cpp $(HASHFLG) $(BUILDFLG) $(OMPFLG)


	tar -zxvf Database/queryByOTU/gg_13.txt.tar.gz 			-C ./Database/queryByOTU/
	tar -zxvf Database/queryByOTU/gg_13_99.txt.tar.gz 		-C ./Database/queryByOTU/
	tar -zxvf Database/queryByOTU/silva.txt.tar.gz 			-C ./Database/queryByOTU/
	tar -zxvf Database/queryByOTU/ncbi_refseq.txt.tar.gz 	-C ./Database/queryByOTU/

	tar -zxvf Database/queryByTaxonomy/gg_13.txt.tar.gz 		-C ./Database/queryByTaxonomy/
	tar -zxvf Database/queryByTaxonomy/gg_13_99.txt.tar.gz 		-C ./Database/queryByTaxonomy/
	tar -zxvf Database/queryByTaxonomy/silva.txt.tar.gz 		-C ./Database/queryByTaxonomy/
	tar -zxvf Database/queryByTaxonomy/ncbi_refseq.txt.tar.gz 	-C ./Database/queryByTaxonomy/

clean:
	rm -rf bin/PM-* src/*.o

