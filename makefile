CC=clang++#g++#

ifeq ($(CC),g++)
	CFLAGS= -Wall -Wextra -g -O3 -std=c++14 -fmax-errors=5   
else
	CFLAGS= -Wall -Wextra -g -O3 -std=c++14 -ferror-limit=5
endif

LIBS= -lm \
	-lbbhutil ### to save to sdf 
#	-lhdf5 -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ ### to save to hdf5
	
LDFLAGS= #-I /usr/include/hdf5/serial/ ### to save to hdf5 

OBJECTS= main.o \
	edgb.o \
	indep_res.o \
	radial_pts.o \
	field.o \
	io_sdf.o \
	sim_params.o \
	initial_data.o \
	fd_stencils.o 

DEPS_MAIN= edgb.hpp \
	indep_res.hpp \
	radial_pts.hpp \
	field.hpp \
	io_sdf.hpp \
	sim_params.hpp \
	initial_data.hpp 

DEPS_SP= sim_params.hpp

DEPS_EDGB= edgb.hpp \
	radial_pts.hpp \
	field.hpp \
	fd_stencils.hpp 

DEPS_RP= radial_pts.hpp  

DEPS_FIELDS= field.hpp

DEPS_FD= fd_stencils.hpp

DEPS_IO_SDF= io_sdf.hpp

DEPS_ID= initial_data.hpp \
	sim_params.hpp 

DEPS_INDEP_RES= indep_res.hpp \
	field.hpp

run: $(OBJECTS)
	$(CC) -o run $(OBJECTS) $(LDFLAGS) $(LIBS) $(CFLAGS)

test: $(OBJECTS)
	$(CC) -o test $(OBJECTS) $(LDFLAGS) $(LIBS) $(CFLAGS)

main.o: main.cpp $(DEPS_MAIN)
	$(CC) -c main.cpp $(LDFLAGS) $(CFLAGS)

sim_params.o: sim_params.cpp $(DEPS_SP)
	$(CC) -c sim_params.cpp $(LDFLAGS) $(CFLAGS)

radial_pts.o: radial_pts.cpp $(DEPS_RP)
	$(CC) -c radial_pts.cpp $(LDFLAGS) $(CFLAGS)

edgb.o: edgb.cpp $(DEPS_EDGB)
	$(CC) -c edgb.cpp $(LDFLAGS) $(CFLAGS)

field.o: field.cpp $(DEPS_FIELDS)
	$(CC) -c field.cpp $(LDFLAGS) $(CFLAGS)

fd_stencils.o: fd_stencils.cpp $(DEPS_FD)
	$(CC) -c fd_stencils.cpp $(LDFLAGS) $(CFLAGS)

io_sdf.o: io_sdf.cpp $(DEPS_IO_SDF)
	$(CC) -c io_sdf.cpp $(LDFLAGS) $(CFLAGS)

indep_res.o: indep_res.cpp $(DEPS_INDEP_RES)
	$(CC) -c indep_res.cpp $(LDFLAGS) $(CFLAGS)

initial_data.o: initial_data.cpp $(DEPS_ID)
	$(CC) -c initial_data.cpp $(LDFLAGS) $(CFLAGS)
#############################################################################
.PHONY: clean_o
clean_o:
	@rm *.o 
.PHONY: clean_run
clean_run:
	@rm run 
.PHONY:clean_test
clean_TEST:
	@rm -r test 
