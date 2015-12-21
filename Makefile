VPATH= . src eqns examples src/triangle src/Timer

# debug flags
CFLAGS = -g -O0
# optimized flags
#CFLAGS = -O3
# turn on all warnings, use C++11

WARNINGFLAGS = -Wall -Wshadow -pedantic

CFLAGS += $(WARNINGFLAGS) 

CPPFLAGS = -std=c++11

CC=clang
CXX=clang++

ifeq ($(CXX) , clang++)
  CPPFLAGS += -stdlib=libc++
endif

SYSTEM_INCLUDE_DIR := /usr/local/include
SYSTEM_LIB_DIR := /usr/local/lib

INCLUDES := -I$(SYSTEM_INCLUDE_DIR) -Ivoro++_2d/src -Isrc -Ieqns
LIBS := -L$(SYSTEM_LIB_DIR) -Llib -larmadillo -lgsl -lvoro++_2d -lsuperlu

CFLAGS += $(INCLUDES)

LIBCSRC := triangle/triangle.c
LIBSRC := PolyMesh.cpp MeshFn.cpp Meshes.cpp Triangulation.cpp Functors.cpp \
	   Quadrature.cpp Legendre.cpp MassMatrix.cpp TimeIntegration.cpp \
	   Timer/CH_Timer.cpp \
	   Advection.cpp Euler.cpp EulerVortex.cpp
LIBOBJS := $(addprefix build/, $(notdir $(patsubst %.cpp,%.o, $(LIBSRC)))) \
		   $(addprefix build/, $(notdir $(patsubst %.c,%.o, $(LIBCSRC))))
OUTPUT := lib/libcpolyg.a examples/adv examples/eul

ADVOBJS := build/adv.o
EULOBJS := build/eul.o

OBJS := $(LIBOBJS) $(ADVOBJS) $(EULOBJS)
DEPS := $(addprefix build/, $(notdir $(patsubst %.o,%.d, $(OBJS))))

eul: $(LIBOBJS) $(EULOBJS) Makefile
	$(CXX) $(CFLAGS) $(LIBOBJS) $(EULOBJS) $(LIBS) -o examples/eul

adv: $(LIBOBJS) $(ADVOBJS) Makefile
	$(CXX) $(CFLAGS) $(LIBOBJS) $(ADVOBJS) $(LIBS) -o examples/adv

# all: cpolyg adv eul Makefile
# 
# cpolyg: $(LIBOBJS) Makefile
# 	rm -f lib/libcpolyg.a
# 	ar rs lib/libcpolyg.a $(LIBOBJS)
# 
# adv: $(ADVOBJS) Makefile
# 	$(CXX) $(CFLAGS) $(ADVOBJS) $(LIBS) -lcpolyg -o examples/adv
# 
# eul: $(EULOBJS) Makefile
# 	$(CXX) $(CFLAGS) $(EULOBJS) $(LIBS) -lcpolyg -o examples/eul

build/%.o: %.cpp Makefile
	$(CXX) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
	$(CXX) -MM -MT '$@' $(CFLAGS) $(CPPFLAGS) $< > build/$*.d

build/%.o: %.c Makefile
	$(CC) -c $(CFLAGS) $< -o $@
	$(CC) -MM -MT '$@' $(CFLAGS) $< > build/$*.d

clean:
	rm -f $(OBJS) $(DEPS) $(OUTPUT)

-include $(DEPS)
