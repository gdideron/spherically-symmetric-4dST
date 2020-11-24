#=============================================================================
TOP= ./
INC= $(TOP)include/
BIN= $(TOP)bin/
SRC= $(TOP)src/
OBJDIR= $(TOP)obj/
OUT= $(TOP)output/

SYSLIB= -lm -lbbhutil

## libraries for bbhutil (to save output as .sdf)

INCBBHUTIL= $(HOME)/rnpletal/include 
LIBBBHUTIL= $(HOME)/rnpletal/lib 

VPATH= $(BIN) $(SRC) $(INC) $(OBJDIR)
#=============================================================================
CC=g++
CFLAGS= -Wall -Wextra -g -O3 -std=c++14 -fmax-errors=5
#=============================================================================
OBJ= $(addprefix $(OBJDIR), \
	main.o \
	edgb.o \
	indep_res.o \
	radial_pts.o \
	field.o \
	io_sdf.o \
	io_csv.o \
	sim_params.o \
	initial_data.o \
	fd_stencils.o \
	)
DEPS= 	edgb.hpp \
	indep_res.hpp \
	radial_pts.hpp \
	field.hpp \
	io_sdf.hpp \
	io_csv.hpp \
	sim_params.hpp \
	initial_data.hpp \
	fd_stencils.hpp
#=============================================================================
all: default.run
test: $(TEST)
#=============================================================================
%.run: $(OBJ)
	$(CC) $^ -o $(BIN)$@ -I$(INC) -I$(INCBBHUTIL) $(SYSLIB) -L$(LIBBBHUTIL) $(CFLAGS)
#=============================================================================
$(OBJDIR)%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -I$(INC) -I$(INCBBHUTIL) -c -o $@ $<
#=============================================================================
.PHONY: clean clean_out
clean:
	$(RM) $(OBJDIR)*.o 
	$(RM) $(BIN)*.run 
clean_out:
	$(RM) -r $(OUT)*
