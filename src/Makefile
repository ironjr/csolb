INTEL=/opt/intel/compilers_and_libraries_2018.3.222/linux
INTEL_MKL=$(INTEL)/mkl
INTEL_COMP=$(INTEL)/compiler
CPP=g++
INC=$(INTEL_MKL)/include/
LIB=$(INTEL_MKL)/lib/intel64/
LIB2=$(INTEL_COMP)/lib/intel64/
BUILD=../build

clean:
	rm -f $(BUILD)/*
	find . -type f -name '*.o' -exec rm {} +

stress-test: stress-test-main.o solb.o
	$(CPP) -o $(BUILD)/stress-test \
		stress-test/stress-test-main.o \
		core/solb.o \
		-L$(LIB) -L$(LIB2) \
		-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm

solb.o: core/solb.cpp
	(cd core; \
		$(CPP) -Wall -I$(INC) -c solb.cpp)

stress-test-main.o: stress-test/stress-test-main.cpp
	(cd stress-test; \
		$(CPP) -Wall -fopenmp -I$(INC) -c stress-test-main.cpp)

