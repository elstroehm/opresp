CPP = clang++-18
CPP_FLAGS = -std=c++17 -pedantic -Wall -O2 -flto

BUILDDIR = build/$(CPP)/

SOURCEFILES = $(wildcard src/*.cpp)
HEADERFILES = $(wildcard src/*.hpp)
OBJECTFILES = $(patsubst src/%.cpp,lib/$(CPP)/%.o,$(wildcard $(SOURCEFILES)))

EIGENDIR = PATH_TO_PROJECT/external/eigen-3.4.0
CUBADIR = PATH_TO_PROJECT/external/Cuba-4.2.2
JSONDIR = PATH_TO_PROJECT/external/json
HDF5DIR = PATH_TO_PROJECT/external/hdf5-1.14.1-2/hdf5

INCLUDEDIRS = -I$(EIGENDIR) -I$(CUBADIR) -I$(JSONDIR) -I$(HDF5DIR)/include
LIBDIRS = -L$(CUBADIR) -L$(HDF5DIR)/lib

LIBS = -fopenmp -lm -Wl,--no-as-needed -ldl -lcuba -lhdf5_cpp -lhdf5

main: main.cpp $(OBJECTFILES)
	$(CPP) $(CPP_FLAGS) $(INCLUDEDIRS) $(LIBDIRS) -o $(BUILDDIR)main.out main.cpp $(OBJECTFILES) $(LIBS)

setup: setup.cpp $(OBJECTFILES)
	$(CPP) $(CPP_FLAGS) $(INCLUDEDIRS) $(LIBDIRS) -o $(BUILDDIR)setup.out setup.cpp $(OBJECTFILES) $(LIBS)

tls_response: tls_response.cpp $(OBJECTFILES)
	$(CPP) $(CPP_FLAGS) $(INCLUDEDIRS) $(LIBDIRS) -o $(BUILDDIR)tls_response.out tls_response.cpp $(OBJECTFILES) $(LIBS)

qo_enthusiast: qo_enthusiast.cpp $(OBJECTFILES)
	$(CPP) $(CPP_FLAGS) $(INCLUDEDIRS) $(LIBDIRS) -o $(BUILDDIR)qo_enthusiast.out qo_enthusiast.cpp $(OBJECTFILES) $(LIBS)

lib/$(CPP)/%.o: src/%.cpp src/%.hpp
	$(CPP) -c $(CPP_FLAGS) $(INCLUDEDIRS) $< -o $(patsubst src/%.cpp,lib/$(CPP)/%.o,$<)

# $< means "the first entry in the list of prerequisites"
# see also: "automatic variables in make"

