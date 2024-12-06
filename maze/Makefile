CXX = g++
CXXFLAGS = -O3 -Wall -Werror
LIBS = -lgmp -lgmpxx

EXES = maze recurrence

all: $(EXES)

%: %.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LIBS)

clean:
	rm $(EXES)
