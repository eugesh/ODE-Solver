all: plot

EXEC_NAME=odesolver
CMP=g++-12
CMP_KEYS=-Wall -Wextra -Wfloat-equal -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual -pedantic
CMP_OPT=-O3
SRC=$(wildcard src/*.cpp src/*.h)

plot: $(wildcard plot.plt) run
	gnuplot plot.plt

run: build $(wildcard data.dat)
	./$(EXEC_NAME)

build: $(SRC)
	$(CMP) -O3 -std=c++20 $(CMP_KEYS) $(CMP_OPT) $(SRC) -o $(EXEC_NAME)

clear:
	rm -rf fft result.dat abs.dat
