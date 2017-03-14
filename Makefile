include make.inc

PROJECT=hmmaw
CC=gcc
CXX=g++
LDFLAGS=
CFLAGS=

UNAME=$(shell uname -s)
ifeq ($(UNAME), Linux)
	LINUX := 1
else ifeq ($(UNAME), Darwin)
	OSX := 1
endif

ifeq ($(LINUX),1)
MEXCC=/usr/global/matlab/R2016a/bin/mex
MATLABINC=/usr/global/matlab/R2015a/extern/include/ -fPIC
endif
ifeq ($(OSX),1)
MEXCC=/Applications/MATLAB_R2016a.app/bin/mex
MATLABINC=/Applications/MATLAB_R2015a.app/extern/include/ -fPIC
endif


INCLUDES=-Iinclude/ -I$(MOSEK)/h $(CBLAS_INC) \
 -Igoogletest/googletest/include \
 -I/usr/local/Cellar/gflags/2.2.0/include \
 -I/usr/local/Cellar/glog/0.3.4_1/include
LIBRARIES=-L$(MOSEK)/bin -lmosek64.7.1 -Wl,-rpath,. -Wl,-rpath,$(MOSEK)/bin $(BLAS_LIB) \
 -L/usr/local/Cellar/glog/0.3.4_1/lib -lglog \
 -L/usr/local/Cellar/gflags/2.2.0/lib -lgflags
BUILD_DIR=build

HEADERS := $(shell find include -name '*.h')
SOURCES_CPP := $(wildcard src/*.cpp src/*/*.cpp)
OBJ := $(patsubst src/%.cpp, build/%.o, $(SOURCES_CPP))
SOURCE_CPP_WITH_MAIN=main.cpp

# OBJ_IN_FOLDER=$(subst src, obj, $(OBJ))
# OBJ = $(patsubst src/%.c,obj/%.o,$(SOURCES_C))

# Test sources:
TEST_MAIN_SRC = ./test/gtest_main.cpp
TEST_SRCS := $(shell find ./test/ -name "*_test.cpp")
TEST_SRCS := $(filter-out $(TEST_MAIN_SRC), $(TEST_SRCS))
TEST_CXX_OBJS := $(addprefix $(BUILD_DIR)/,${TEST_SRCS:.cpp=.o})
TEST_CXX_BINS := $(addsuffix .testbin,$(foreach obj,$(TEST_CXX_OBJS),$(basename $(notdir $(obj)))))
TEST_ALL_BIN := $(TEST_BIN_DIR)/test_all.testbin

GTEST_OBJ := $(addprefix $(BUILD_DIR)/, ${GTEST_SRC:.cpp=.o})
GTEST_SRC := src/gtest/gtest-all.cpp

.PHONY: clean all test

all: $(PROJECT)_train

build/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) $(INCLUDES) -MM -MT build/$*.o $< >build/$*.d
	$(CXX) -c $(CFLAGS) $(INCLUDES) $< -o $@

hmm_fit.mex: hmm_fit_jia.c $(OBJ)
	$(MEXCC) CFLAGS="\$(CFLAGS)" $< $(OBJ)

hmm_likelihood.mex: hmm_likelihood_jia.c $(OBJ)
	$(MEXCC) CFLAGS="\$(CFLAGS)" $< $(OBJ)

ifeq ($(OSX),1)
$(PROJECT)_train: $(SOURCE_CPP_WITH_MAIN) lib/libhmm.a lib/libmosek64_wrapper.dylib
	$(CXX) -o $@ $(CFLAGS) $(INCLUDES) $(LIBRARIES) $^ $(LDFLAGS)

lib/libmosek64_wrapper.dylib: build/mosek_solver.o
	$(CXX) -dynamiclib $(INCLUDES) $(LIBRARIES) -Wl,-install_name,$@ -compatibility_version 7.0 -current_version $(MOSEK_VERSION) -o $@ $< $(MOSEKLIB)
	# install_name_tool -change  @loader_path/libmosek64.$(MOSEK_VERSION).dylib  $(MOSEK)/bin/libmosek64.$(MOSEK_VERSION).dylib libmosek64_wrapper.$(MOSEK_VERSION).dylib
	# ln -sf libmosek64_wrapper.$(MOSEK_VERSION).dylib $@ 

lib/libhmm.a: $(filter-out lib/solver_mosek.o, $(OBJ))
	@echo $^
	@mkdir -p $(@D)
	ar crv $@ $^

endif

ifeq ($(LINUX),1)
$(PROJECT)_train: lib/libhmm.a lib/libmosek64_wrapper.so
	$(CXX) -o $@ $(CFLAGS) $(INCLUDES) $(LIBRARIES) $^ $(LDFLAGS)

lib/libmosek64_wrapper.so: build/mosek_solver.o
	$(CXX) -shared $(LDFLAGS) $(INCLUDES) $(LIBRARIES) -Wl,-soname,$@ -o libmosek64_wrapper.$(MOSEK_VERSION).so $< $(MOSEKLIB)
	ln -sf libmosek64_wrapper.$(MOSEK_VERSION).so $@ 

lib/libhmm.a: $(OBJ) libmosek64_wrapper.dylib
	@echo $^
	@mkdir -p $(@D)
	ar crv $@ $^

endif

clean: 
	rm -f ./build/*.o ./build/*.d

run:
	./train -i ../../data/gmm_hmm_samples_dim2_T1000.dat -m md.dat -d 2 -n 1 -s 2 -l 1000

test: $(TEST_ALL_BIN) $(TEST_CXX_BINS)

testvars:
	@echo $(HEADERS)
	@echo $(SOURCES_CPP)
	@echo $(OBJ)
	@echo $(SOURCE_CPP_WITH_MAIN)
	@echo $(MEXCC)
	@echo $(UNAME)
