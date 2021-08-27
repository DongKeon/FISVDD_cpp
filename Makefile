.PHONY: all clean

all: main 
	./test

test: simple_demo.cpp
	g++ -o test -Wall -Wpedantic -Wextra -std=c++11 simple_demo.cpp

clean:
	rm main
