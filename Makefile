#############################################################################
TOP= ./
INC= $(TOP)include/
BIN= $(TOP)bin/
SRC= $(TOP)src/
OBJDIR= $(TOP)obj/

SYSLIB= -lm -lbbhutil

INCBBHUTIL= /home/jripley/rnpletal/include 
LIBBBHUTIL= /home/jripley/rnpletal/lib 

VPATH= $(BIN):$(SRC):$(INC):$(OBJDIR)
#############################################################################
CC=g++
CFLAGS= -Wall -Wextra -g -O2 -std=c++11 -fmax-errors=5   
#############################################################################
OBJ= $(addprefix $(OBJDIR), \
	main.o \
	edgb.o \
	indep_res.o \
	radial_pts.o \
	field.o \
	io_sdf.o \
	sim_params.o \
	initial_data.o \
	fd_stencils.o \
	)
DEPS= 	edgb.hpp \
	indep_res.hpp \
	radial_pts.hpp \
	field.hpp \
	io_sdf.hpp \
	sim_params.hpp \
	initial_data.hpp \
	fd_stencils.hpp 
#############################################################################
RUN= $(BIN)run
TEST= $(BIN)test
all: $(RUN)
test: $(TEST)
#############################################################################
$(RUN): $(OBJ)
	$(CC) $^ -o $@ -I$(INC) -I$(INCBBHUTIL) $(SYSLIB) -L$(LIBBBHUTIL) $(CFLAGS)
$(TEST): $(OBJ)
	$(CC) $^ -o $@ -I$(INC) -I$(INCBBHUTIL) $(SYSLIB) -L$(LIBBBHUTIL) $(CFLAGS)
#############################################################################
$(OBJDIR)%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -I$(INC) -I$(INCBBHUTIL) -c -o $@ $<
#############################################################################
.PHONY: clean_o clean_run clean_test
clean_o:
	@rm $(OBJDIR)*.o 
clean_run:
	@rm $(BIN)run 
clean_test:
	@rm -r $(BIN)test 
