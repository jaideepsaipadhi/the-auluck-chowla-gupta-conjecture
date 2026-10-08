CC ?= gcc

all: ws_fc_exact ws_fc_check

ws_fc_exact: ws_fc_exact.c
	$(CC) -O3 -march=native -Wall -o $@ $<

ws_fc_check: ws_fc_check.c
	$(CC) -O2 -fopenmp -frounding-math -o $@ $< -lm

clean:
	rm -f ws_fc_exact ws_fc_check

.PHONY: all clean
