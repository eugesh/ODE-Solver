all: plot

EXEC_NAME=odesolver
CMP=g++-12
CMP_KEYS=-Wall -Wextra -Wfloat-equal -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual -pedantic
CMP_OPT=-O3

plot: $(wildcard plot.plt)
	gnuplot plot.plt

run: build $(wildcard data.dat)
	./$(EXEC_NAME)

build: $(wildcard *.cpp *.h Makefile)
	$(CMP) -std=c++20 $(CMP_KEYS) $(CMP_OPT) *.cpp -o $(EXEC_NAME)

clear:
	rm -rf fft result.dat abs.dat
