# Makefile

bin/main : src/main.cpp include/sparsematrix.h include/sparsealgs.h include/randommatrix.h Makefile
	g++ -I include/ -o bin/main src/main.cpp

.PHONY : clean
clean : 
	rm -f bin/main
