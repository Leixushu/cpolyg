VPATH= . src eqns examples

# debug flags
CFLAGS = -g -O0
# optimized flags
#CFLAGS = -O3
# turn on all warnings, use C++11

WARNINGFLAGS = -Wall -Wshadow -pedantic

ALLWARNINGS = -Weverything -Wno-extra-semi -Wno-sign-conversion -Wno-padded \
			  -Wno-sign-compare -Wno-shorten-64-to-32 -Wno-old-style-cast -Wno-weak-vtables \
			  -Wno-c++98-compat -Wno-global-constructors -Wno-missing-prototypes

CFLAGS += $(WARNINGFLAGS) 

CPPFLAGS = -std=c++11

CC=clang
CXX=clang++

ifeq ($(CXX) , clang++)
  CPPFLAGS += -stdlib=libc++
endif

SYSTEM_INCLUDE_DIR := /usr/local/include
SYSTEM_LIB_DIR := /usr/local/lib

INCLUDES := -I$(SYSTEM_INCLUDE_DIR) -isystemvoro++_2d/src -Isrc -Ieqns
LIBS := -L$(SYSTEM_LIB_DIR) -Llib -larmadillo -lgsl -lvoro++_2d -lsuperlu

CFLAGS += $(INCLUDES)

LIBCSRC := triangle.c
LIBCOBJ := $(addprefix build/, $(notdir $(patsubst %.c,%.o, $(LIBCSRC))))
LIBSRC := PolyMesh.cpp MeshFn.cpp Meshes.cpp Triangulation.cpp Functors.cpp \
	   Quadrature.cpp Legendre.cpp MassMatrix.cpp TimeIntegration.cpp \
	   Advection.cpp Euler.cpp EulerVortex.cpp
LIBOBJS := $(addprefix build/, $(notdir $(patsubst %.cpp,%.o, $(LIBSRC))))
OUTPUT := lib/libcpolyg.a examples/adv examples/eul

ADVOBJS := build/adv.o
EULOBJS := build/eul.o

OBJS := $(LIBOBJS) $(LIBCOBJ) $(ADVOBJS) $(EULOBJS)
DEPS := $(addprefix build/, $(notdir $(patsubst %.o,%.d, $(OBJS))))

all: cpolyg adv eul Makefile

cpolyg: $(LIBOBJS) $(LIBCOBJ) Makefile
	rm -f lib/libcpolyg.a
	ar rs lib/libcpolyg.a $(LIBOBJS) $(LIBCOBJ)

adv: $(ADVOBJS) Makefile
	$(CXX) $(CFLAGS) $(ADVOBJS) $(LIBS) -lcpolyg -o examples/adv

eul: $(EULOBJS) Makefile
	$(CXX) $(CFLAGS) $(EULOBJS) $(LIBS) -lcpolyg -o examples/eul

build/%.o: %.cpp Makefile
	$(CXX) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
	$(CXX) -MM -MT '$@' $(CFLAGS) $(CPPFLAGS) $< > build/$*.d

build/%.o: %.c Makefile
	$(CC) -c $(CFLAGS) $< -o $@
	$(CC) -MM -MT '$@' $(CFLAGS) $< > build/$*.d

clean:
	rm -f $(OBJS) $(DEPS) $(OUTPUT)

-include $(DEPS)
