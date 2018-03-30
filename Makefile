PROJECT=hmmaw
CC=gcc
CXX=g++
LDFLAGS=
CFLAGS=-std=c++11 -fPIC $(AW_DEFINES)

UNAME=$(shell uname -s)
HNAME=$(shell hostname)
IDENTIFIER=$(shell echo $(HNAME) | cut -d'.' -f 3)
IDENTIFIER1=$(shell echo $(HNAME) | cut -d'.' -f 1)
#lionxv.rcc.psu.edu
#cyberstar.psu.edu
#br005.pvt.bridges.psc.edu
ifeq ($(UNAME), Linux)
	LINUX := 1
	ifeq ($(IDENTIFIER1), lionxv)
		include make_inc/make.inc.Linux.Cyberstar
	else ifeq ($(IDENTIFIER1), cyberstar)
		include make_inc/make.inc.Linux.Cyberstar
	else ifeq ($(IDENTIFIER), bridges)
		include make_inc/make.inc.Linux
	else ifeq ($(IDENTIFIER), ESC8000)
	  include make_inc/make.inc.esc8000
	else ifneq (,$(findstring yukun-dell,$(IDENTIFIER)))
	  include make_inc/make.inc.home-dell
	endif
else ifeq ($(UNAME), Darwin)
	OSX := 1
	include make_inc/make.inc.macOS
endif


INCLUDES=-Iinclude/ -I$(MOSEK)/h $(CBLAS_INC) \
 -Itest \
 -Igoogletest/googletest/include -Igoogletest/googletest\
 -I$(GFLAGS)/include \
 -I$(GLOG)/include
LIBRARIES=-Wl,-rpath,. $(BLAS_LIB) \
 -L$(GFLAGS)/lib -lgflags -Wl,-rpath,$(GFLAGS)/lib \
 -L$(GLOG)/lib -lglog -Wl,-rpath,$(GLOG)/lib \
 -L./lib
BUILD_DIR=build

HEADERS := $(shell find include -name '*.h')
SOURCES_CPP := $(wildcard src/*.cpp src/*/*.cpp)
OBJ := $(patsubst src/%.cpp, build/%.o, $(SOURCES_CPP))
SOURCE_CPP_WITH_MAIN=main.cpp

# OBJ_IN_FOLDER=$(subst src, obj, $(OBJ))
# OBJ = $(patsubst src/%.c,obj/%.o,$(SOURCES_C))

# Test sources:
GTEST_ALL = ./googletest/googletest/src/gtest-all.cc
TEST_MAIN_SRC = ./test/gtest_main.cpp
#TEST_SRCS := $(shell find ./test -name "*_test.cpp")
TEST_SRCS := $(shell find ./test -name "*.cpp")
TEST_SRCS := $(filter-out $(TEST_MAIN_SRC), $(TEST_SRCS))
TEST_CXX_OBJS := $(patsubst ./test/%.cpp, ./build/%.o, $(TEST_SRCS))
TEST_ALL_BIN := lib/hmm_awTest

GTEST_OBJ := $(addprefix $(BUILD_DIR)/, ${GTEST_SRC:.cpp=.o})
GTEST_SRC := src/gtest/gtest-all.cpp

.PHONY: clean all test

all: $(PROJECT)_train

build/%.o: src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CFLAGS) $(INCLUDES) -MM -MT build/$*.o $< >build/$*.d
	$(CXX) -c $(CFLAGS) $(INCLUDES) $< -o $@

ifeq ($(OSX),1)
$(PROJECT)_train: $(SOURCE_CPP_WITH_MAIN) lib/libhmm.a
	$(CXX) -o $@ $(CFLAGS) $(INCLUDES) -L./lib -lhmm $(LIBRARIES) $(SOURCE_CPP_WITH_MAIN) $(LDFLAGS)

lib/libhmm.a: $(filter-out build/mosek_solver.o, $(OBJ))
	@echo $^
	@mkdir -p $(@D)
	ar crv $@ $^

endif

ifeq ($(LINUX),1)
$(PROJECT)_train: $(SOURCE_CPP_WITH_MAIN) lib/libhmm.so
	$(CXX) -o $@ $(CFLAGS) $(INCLUDES) -L./lib -lhmm -Wl,-rpath=./lib $(LIBRARIES) $(SOURCE_CPP_WITH_MAIN) $(LDFLAGS)

lib/libhmm.so: $(filter-out build/mosek_solver.o, $(OBJ))
	@echo $^
	@mkdir -p $(@D)
	$(CXX) -shared -o $@ $^ $(LIBRARIES)
	#ar crv $@ $^

endif

clean:
	rm -f ./build/*.o ./build/*.d
	rm $(PROJECT)_train

run:
	./$(PROJECT)_train --infilename ./data/gmm_hmm_samples_dim2_T1000.dat --mdfilename /tmp/md.dat --dim 2 --num 1 --statenum 2 --len 1000 --forcediag
	# ./$(PROJECT)_train -i ./data/gmm_hmm_samples_dim2_T1000.dat -m ./data/md.dat -d 2 -n 1 -s 2 -l 1000

test: $(TEST_ALL_BIN)

runtest:
	$(TEST_ALL_BIN)

$(TEST_ALL_BIN): $(TEST_MAIN_SRC) $(TEST_SRCS) lib/libhmm.so lib/libgtest.so
	@echo $@
	$(CXX) $(TEST_MAIN_SRC) $(TEST_SRCS) -o $@ $(CFLAGS) $(INCLUDES) $(LIBRARIES) -lhmm -lgtest -lm -lpthread -Wl,-rpath,./lib

testvars:
	@echo $(IDENTIFIER)
	@echo $(INCLUDES)
	@echo $(HEADERS)
	@echo $(SOURCES_CPP)
	@echo $(OBJ)
	@echo $(SOURCE_CPP_WITH_MAIN)
	@echo $(MEXCC)
	@echo $(UNAME)
	@echo $(filter-out build/mosek_solver.o, $(OBJ))
	@echo $(TEST_SRCS)
	@echo $(TEST_CXX_OBJS)
	@echo $(IN_HOME_DELL)
	@echo $(findstring yukun-dell,$(IDENTIFIER))
