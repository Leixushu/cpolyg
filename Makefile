VPATH= . src eqns

# debug flags
CFLAGS = -g -O0
# optimized flags
#CFLAGS = -O3
# turn on all warnings, use C++11
CFLAGS += -Wall -pedantic 

CPPFLAGS = -std=c++11

CC=clang
CXX=clang++

ifeq ($(CXX) , clang++)
  CPPFLAGS += -stdlib=libc++
endif

SYSTEM_INCLUDE_DIR := /usr/local/include
SYSTEM_LIB_DIR := /usr/local/lib

INCLUDES := -I$(SYSTEM_INCLUDE_DIR) -Ivoro++_2d/src -Isrc
LIBS := -L$(SYSTEM_LIB_DIR) -Llib -larmadillo -lgsl -lvoro++_2d 

CFLAGS += $(INCLUDES)

CSRC := triangle.c
SRC := main.cpp PolyMesh.cpp MeshFn.cpp Meshes.cpp Triangulation.cpp \
	   Quadrature.cpp Legendre.cpp \
	   Advection.cpp

OBJC := $(addprefix build/, $(notdir $(patsubst %.c,%.o, $(CSRC))))
OBJS := $(addprefix build/, $(notdir $(patsubst %.cpp,%.o, $(SRC))))
DEPS := $(addprefix build/, $(notdir $(patsubst %.cpp,%.d, $(SRC))))
EXEC := cpolyg

build/%.o: %.cpp Makefile
	$(CXX) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
	$(CXX) -MM -MT '$@' $(CFLAGS) $(CPPFLAGS) $< > build/$*.d

build/%.o: %.c Makefile
	$(CC) -c $(CFLAGS) $< -o $@
	$(CC) -MM -MT '$@' $(CFLAGS) $< > build/$*.d

cployg: Makefile $(OBJS) $(OBJC)
	$(CXX) $(CFLAGS) $(OBJC) $(OBJS) $(LIBS) -o $(EXEC)

clean:
	rm -f $(OBJS) $(DEPS) $(EXEC)

-include $(DEPS)
