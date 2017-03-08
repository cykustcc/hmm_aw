PROJECT=train
CC=gcc
LDFLAGS=-lm
CFLAGS=-g -Iinclude/

HEADERS := $(shell find include -name '*.h')
SOURCES_C := $(shell find src -name '*.c')
SOURCE_C_WITH_MAIN=$(PROJECT).c
OBJ := $(SOURCES_CPP:%.cpp=%.o) $(SOURCES_C:%.c=%.o)
# OBJ_IN_FOLDER=$(subst src, obj, $(OBJ))
# OBJ = $(patsubst src/%.c,obj/%.o,$(SOURCES_C))

# Test sources:
TEST_MAIN_SRC = ../hmmjiaTest/gtest_main.cpp
TEST_SRCS := $(shell find ../hmmjiaTest/ -name "*_test.cpp")
TEST_SRCS := $(filter-out $(TEST_MAIN_SRC), $(TEST_SRCS))
TEST_CXX_OBJS := ${TEST_SRCS:.cpp=.o}
TEST_CXX_BINS := $(addsuffix .testbin,$(foreach obj,$(TEST_CXX_OBJS),$(basename $(notdir $(obj)))))
TEST_ALL_BIN := $(TEST_BIN_DIR)/test_all.testbin

GTEST_OBJ := $(addprefix $(BUILD_DIR)/, ${GTEST_SRC:.cpp=.o})
GTEST_SRC := src/gtest/gtest-all.cpp

all: $(PROJECT)

obj:
	@mkdir -p $@

%.o:  %.c %.h
	$(CC) -o $@ $(CFLAGS) -c $+

$(PROJECT): $(OBJ)
	$(CC) -o $@ $(CFLAGS) $(SOURCE_C_WITH_MAIN) $^ $(LDFLAGS)

clean: 
	rm -f ./src/*.o

run:
	./train -i ../../data/gmm_hmm_samples_dim2_T1000.dat -m md.dat -d 2 -n 1 -s 2 -l 1000

test: $(TEST_ALL_BIN) $(TEST_CXX_BINS)

testvars:
	@echo $(SOURCES_C)
	@echo $(OBJ)
	@echo $(OBJ_IN_FOLDER)
	@echo $(SOURCE_C_WITH_MAIN)
