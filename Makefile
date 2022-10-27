# Makefile

bin/main : src/main.cpp include/sparsematrix.h include/sparsealgs.h Makefile
	g++ -I include/ -o bin/main src/main.cpp

.PHONY : clean
clean : 
	rm -f bin/main
