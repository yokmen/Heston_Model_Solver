CC = cc

# Production build: use -O3 for maximum performance with OpenMP
# Note: Changed from -O0 (no optimization) to -O3 (aggressive optimization)
# This is critical for OpenMP performance - compilation optimizations
# can provide 2-3x speedup for this numerical code
CFLAGS = -O3 -Wall -fopenmp $(shell pkg-config --cflags hdf5-serial 2>/dev/null || pkg-config --cflags hdf5)

# Debug build: slower compilation but easier debugging
CFLAGS_DEBUG = -O0 -g -Wall -fopenmp $(shell pkg-config --cflags hdf5-serial 2>/dev/null || pkg-config --cflags hdf5)

# Added -fopenmp here for the linker
LIB = \
SuiteSparse/UMFPACK/Lib/libumfpack.a \
SuiteSparse/CHOLMOD/Lib/libcholmod.a \
SuiteSparse/AMD/Lib/libamd.a \
SuiteSparse/CAMD/Lib/libcamd.a \
SuiteSparse/COLAMD/Lib/libcolamd.a \
SuiteSparse/CCOLAMD/Lib/libccolamd.a \
SuiteSparse/metis-4.0/libmetis.a \
SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a \
-lm -lblas -llapack -fopenmp \
$(shell pkg-config --libs hdf5-serial 2>/dev/null || pkg-config --libs hdf5)

IFLAGS = \
-ISuiteSparse/UMFPACK/Include \
-ISuiteSparse/SuiteSparse_config \
-ISuiteSparse/AMD/Include

OBJS = main.o prob.o grid_index.o grid.o norme.o greeks.o \
Crank_Nicolson.o dataset.o time.o umfpack.o hdf5_writer.o

default: main

main: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIB)

# Debug build: slower but easier to debug
debug: clean
	$(MAKE) CFLAGS="$(CFLAGS_DEBUG)" main

# Release build (default): optimized for speed
release: clean
	$(MAKE) CFLAGS="$(CFLAGS)" main

%.o: %.c
	$(CC) $(CFLAGS) $(IFLAGS) -c $< -o $@

clean:
	rm -f *.o main