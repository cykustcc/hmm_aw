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

test:
	@echo $(SOURCES_C)
	@echo $(OBJ)
	@echo $(OBJ_IN_FOLDER)
	@echo $(SOURCE_C_WITH_MAIN)
