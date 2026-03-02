LIBOMPPATH = /opt/homebrew/opt/libomp
CC = gcc
CFLAGS = -O3 -Wall -Werror

CXX = g++
CXXFLAGS = -O3 -Wall -Werror

EXES = domino-dp domino-poly domino-ff

all: $(EXES)

domino-ff: domino-ff.c
	$(CC) $< -o $@ $(CFLAGS) -march=native -Xpreprocessor -fopenmp -I$(LIBOMPPATH)/include -L$(LIBOMPPATH)/lib -lgmp -lomp

domino-poly: domino-poly.c
	$(CC) $< -o $@ $(CFLAGS) -lgmp

domino-dp: domino-dp.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) -lgmp -lgmpxx

.PHONY: clean

clean:
	@rm $(EXES)
