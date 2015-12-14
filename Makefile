VPATH= . 

# debug flags
CFLAGS = -g -O0
# optimized flags
#CFLAGS = -O3
# turn on all warnings, use C++11
CFLAGS += -Wall -std=c++11

CXX=clang++

ifeq ($(CXX) , clang++)
  CFLAGS += -stdlib=libc++
endif

SYSTEM_INCLUDE_DIR := /usr/local/include
SYSTEM_LIB_DIR := /usr/local/lib

INCLUDES := -I$(SYSTEM_INCLUDE_DIR)
LIBS := -L$(SYSTEM_LIB_DIR) -larmadillo -lgsl

CFLAGS += $(INCLUDES)

SRC := main.cpp PolyMesh.cpp MeshFn.cpp leg.cpp Meshes.cpp

OBJS := $(addprefix o/, $(notdir $(patsubst %.cpp,%.o, $(SRC))))
DEPS := $(addprefix d/, $(notdir $(patsubst %.cpp,%.d, $(SRC))))
EXEC := cpolyg

o/%.o: %.cpp Makefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > d/$*.d

cployg: Makefile $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

clean:
	rm -f $(OBJS) $(DEPS) $(EXEC)

-include $(DEPS)
