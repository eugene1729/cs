CC = gcc
CFLAGS = -Wall -Werror
LIBS = -lgmp
EXES = domino

all: $(EXES)

%: %.c
	$(CC) $< -o $@ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm $(EXES)
