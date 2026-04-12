LIBOMPPATH = /opt/homebrew/opt/libomp
CXX = g++
CXXFLAGS = -O3 -Wall -Werror -Xpreprocessor -fopenmp -I$(LIBOMPPATH)/include -L$(LIBOMPPATH)/lib
LIBS = -lgmp -lgmpxx -lomp

EXES = maze maze-low-memory

all: $(EXES)

%: %.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LIBS)

clean:
	@rm $(EXES)
