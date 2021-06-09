
include ./Makefile.config

# level-1 modules: self-contained
# --level-2 modules: use level-1 modules and are used by level-3 modules
#   --level-3 modules: use both level-1 and level-2 modules
#
OBJS = $(DIR_OBJ)/global.o    \
       $(DIR_OBJ)/parser.o    \
       $(DIR_OBJ)/datatype.o  \
       $(DIR_OBJ)/mmi_mpi.o   \
       $(DIR_OBJ)/string.o    \
       $(DIR_OBJ)/su.o        \
       $(DIR_OBJ)/io.o        \
       $(DIR_OBJ)/math.o      \
       $(DIR_OBJ)/pml.o       \
       $(DIR_OBJ)/source.o    \
       $(DIR_OBJ)/fdcore.o       \
       $(DIR_OBJ)/solver.o

all: dir fd2dmpi

dir:
	@-test -d $(DIR_MOD) || mkdir $(DIR_MOD)
	@-test -d $(DIR_OBJ) || mkdir $(DIR_OBJ)
#	@-test -d $(DIR_EXE) || mkdir $(DIR_EXE)

fd2dmpi: $(OBJS) $(DIR_OBJ)/fd2dmpi.o
	$(FCC) -o $@ $^ $(LDFLAGS)
	mv $@ ../$(DIR_EXE)

bin2su.x: $(OBJS) $(DIR_OBJ)/bin2su.o
	$(FCC) -o $@ $^ $(LDFLAGS)
	mv $@ ../$(DIR_EXE)

$(DIR_OBJ)/%.o: %.f90
	$(FCC) -o $@ $< $(F90FLAGS)

clean:
	rm -f $(DIR_EXE)/* $(DIR_OBJ)/*.o $(DIR_MOD)/*.mod;

