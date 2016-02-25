VPATH= . src eqns examples src/triangle src/Timer src/blas

CC=clang  
CXX=clang++
#CC=gcc-5
#CXX=g++-5

ifeq ($(DEBUG), 1)
	# debug flags
	CFLAGS = -g -O0
else
	# optimized flags
	CFLAGS = -O3
endif

# turn on warnings
WARNINGFLAGS = -Wall -Wshadow -pedantic
CFLAGS += $(WARNINGFLAGS) 
# use C++11
CPPFLAGS = -std=c++11

# is using clang, need to specify the standard library
ifeq ($(CXX) , clang++)
  CPPFLAGS += -stdlib=libc++
endif

SYSTEM_INCLUDE_DIR := /usr/local/include
SYSTEM_LIB_DIR := /usr/local/lib

INCLUDES := -I$(SYSTEM_INCLUDE_DIR) -Ivoro++_2d/src -Isrc -Ieqns
LIBS := -L$(SYSTEM_LIB_DIR) -Llib -larmadillo -lgsl -lvoro++_2d -lblas -llapack
CFLAGS += $(INCLUDES)

LIBCSRC := triangle/triangle.c
LIBSRC := PolyMesh.cpp MeshFn.cpp Meshes.cpp Triangulation.cpp Functors.cpp \
	   Quadrature.cpp Legendre.cpp MassMatrix.cpp TimeIntegration.cpp \
	   Timer/CH_Timer.cpp BlockMatrix.cpp blas/blas.cpp Preconditioners.cpp \
	   Jacobian.cpp Equation.cpp BoundaryConditions.cpp\
	   Advection.cpp \
	   Euler.cpp EulerSetups.cpp
LIBOBJS := $(addprefix build/, $(notdir $(patsubst %.cpp,%.o, $(LIBSRC)))) \
		   $(addprefix build/, $(notdir $(patsubst %.c,%.o, $(LIBCSRC))))

EXAMPLES := ExpAdv ImpAdv PeriodicAdv test Convergence PeriodicAdv \
			ExpEul ImpEul Convergence KH
BINS := $(addprefix bin/, $(EXAMPLES))

OBJS := $(LIBOBJS) $(addprefix build/, $(notdir $(addsuffix .o, $(EXAMPLES))))
DEPS := $(addprefix build/, $(notdir $(patsubst %.o,%.d, $(OBJS))))

.PHONY: all clean debug

all: $(BINS)

debug:
	make -j DEBUG=1 --always-make

bin/%: build/%.o $(LIBOBJS) Makefile
	$(CXX) $(CFLAGS) $(LIBOBJS) $< $(LIBS) -o $@

build/%.o: %.cpp Makefile
	$(CXX) -c $(CFLAGS) $(CPPFLAGS) $< -o $@
	$(CXX) -MM -MT '$@' $(CFLAGS) $(CPPFLAGS) $< > build/$*.d

build/%.o: %.c Makefile
	$(CC) -c $(CFLAGS) -w $< -o $@
	$(CC) -MM -MT '$@' $(CFLAGS) $< > build/$*.d

clean:
	rm -f $(OBJS) $(DEPS) $(BINS)

-include $(DEPS)
