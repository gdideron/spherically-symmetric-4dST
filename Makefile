#=============================================================================
TOP= ./
INC= $(TOP)include/
BIN= $(TOP)bin/
SRC= $(TOP)src/
OBJDIR= $(TOP)obj/
OUT= $(TOP)output/

SYSLIB= -lm

VPATH= $(BIN) $(SRC) $(INC) $(OBJDIR)
#=============================================================================
<<<<<<< HEAD
CC=icpc#g++#
=======
CC=g++#icc#
>>>>>>> 3f67c22eac2b5916a4cd4169403aa263e75e186d
CFLAGS= -g -O2
#==========================================================================
ifeq ($(CC),g++)
	CFLAGS+= -std=c++14 -Wall -Wextra
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
	sim_params.o \
	initial_data.o \
	fd_stencils.o \
	)
DEPS= 	edgb.hpp \
	indep_res.hpp \
	radial_pts.hpp \
	field.hpp \
	io_csv.hpp \
	sim_params.hpp \
	initial_data.hpp \
	fd_stencils.hpp
#=============================================================================
all: default.run
test: $(TEST)
#=============================================================================
%.run: $(OBJ)
	$(CC) $^ -o $(BIN)$@ -I$(INC) $(SYSLIB) $(CFLAGS)
#=============================================================================
$(OBJDIR)%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -I$(INC) -c -o $@ $<
#=============================================================================
.PHONY: clean clean_out
clean:
	$(RM) $(OBJDIR)*.o 
	$(RM) $(BIN)*.run 
clean_out:
	$(RM) -r $(OUT)*
