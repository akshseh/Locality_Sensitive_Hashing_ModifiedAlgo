# source https://github.com/AkselMeola/simple-c-project-template

#!/usr/bin/make -f

VERSION = $(shell cat ../VERSION)



#############################
# Local paths 				#
#############################

INCLUDES =../include
LOCAL_LIBS =../lib
OBJECTS=obj
EXE=exe
TARGET := $(EXE)/main
#############################
# Shared recources 			#
#############################

# Uncomment or add what's necessary

SHARED_CFLAGS = \
	-I$(INCLUDES) \
#	$(shell pkg-config --cflags glib-2.0)	\
#	$(shell pkg-config --cflags json-glib-1.0) \
#	$(shell pkg-config --cflags libcurl) \
#	$(shell pkg-config --cflags purple)

SHARED_LIBS = \
#	$(shell pkg-config --libs glib-2.0)	\
#	$(shell pkg-config --libs json-glib-1.0) \
#	$(shell pkg-config --libs libcurl) \
#	$(shell pkg-config --libs purple)


#############################
# Compiler flags 			#
#############################

# Uncomment or add what's necessary

CC=g++

CFLAGS= \
	$(SHARED_CFLAGS) \
	-g \
	-g3 \
	-Wall \
	-fPIC \
	-Wextra \
	-Wformat \
	-D_GLIBCXX_DEBUG\
	-std=c++ \
	-Wno-unused-parameter \
	-Werror=format-security

LDFLAGS ?= \
#	-shared \
#	-fPIC


LIBS= \
	$(SHARED_LIBS)


# Find all .h files in includes directory
_DEPS = $(wildcard *.h)
DEPS = $(patsubst %, $(INCLUDES)/%, $(_DEPS))

# Find all .c files and rename them as .o files
_OBJ = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
OBJ = $(patsubst %,$(OBJECTS)/%,$(_OBJ))

.PHONY: clean all install

$(OBJECTS)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): $(OBJ)
	$(CC) $(LDFLAGS) -o $@ $^ $(CFLAGS) $(LIBS)

# Clean up
clean:
	rm -f $(TARGET) $(OBJECTS)/*.o *~ core $(INCDIR)/*~
