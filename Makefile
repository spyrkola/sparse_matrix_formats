# Makefile

bin/main : src/main.cpp
	g++ -I include/ -o bin/main src/main.cpp

.PHONY : clean
clean : 
	rm -f bin/main
