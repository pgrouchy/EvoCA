all: evoCA_adder evoCA_adder_test

evoCA_adder: evoCA_adder.o params.h

	g++ -O3 -Wunused-variable -o evoCA_adder evoCA_adder.o

evoCA_adder_test: evoCA_adder_test.o params.h

	g++ -O3 -Wunused-variable -o evoCA_adder_test evoCA_adder_test.o

evoCA_adder.o: evoCA_adder.cpp params.h

	g++ -O3 -Wno-unused-result -Wunused-variable -c -o evoCA_adder.o evoCA_adder.cpp

evoCA_adder_test.o: evoCA_adder_test.cpp params.h

	g++ -O3 -Wno-unused-result -Wunused-variable -c -o evoCA_adder_test.o evoCA_adder_test.cpp

clean:

	rm evoCA_adder evoCA_adder_test evoCA_adder.o evoCA_adder_test.o
