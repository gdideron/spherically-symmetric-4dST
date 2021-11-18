#=============================================================================
TOP= ./
BIN= $(TOP)bin/
SRC= $(TOP)src/
OBJDIR= $(TOP)obj/
OUT= $(TOP)output/

SYSLIB= -lm

VPATH= $(SRC)
#=============================================================================
CC=g++#icpc#
CFLAGS= -g -O2
#==========================================================================
ifeq ($(CC),g++)
	CFLAGS+= -std=c++14 -Wall -Wextra #-fopenmp
		#-fcheck=all  
endif

ifeq ($(CC),icpc)
	CFLAGS+= -std=c++14 -Wall
		#-check-bounds 
endif
#=============================================================================
OBJ= $(addprefix $(OBJDIR), \
	main.o \
	edgb.o \
	indep_res.o \
	radial_pts.o \
	field.o \
	io_csv.o \
	io_sdf.o \
	sim_params.o \
	initial_data.o \
	fd_stencils.o \
	)
DEPS= 	edgb.hpp \
	indep_res.hpp \
	radial_pts.hpp \
	field.hpp \
	io_csv.hpp \
	io_sdf.hpp \
	sim_params.hpp \
	initial_data.hpp \
	fd_stencils.hpp
#=============================================================================
all: default.run
test: $(TEST)
#=============================================================================
%.run: $(OBJ)
	$(CC) $^ -o $(BIN)$@ $(SYSLIB) $(CFLAGS)
#=============================================================================
$(OBJDIR)%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<
#=============================================================================
.PHONY: clean clean_out
clean:
	$(RM) $(OBJDIR)*.o 
	$(RM) $(BIN)*.run 
clean_out:
	$(RM) -r $(OUT)*
