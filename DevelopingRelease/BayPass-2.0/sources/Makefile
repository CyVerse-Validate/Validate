#FC := gfortran
#FC := ifort
ifeq ($(FC),gfortran)
 FFLAGS := -O3 -fno-range-check
## -fno-range-check is to deal with intentional overflow when initializing the mt seed (mt_stream.F90: L110-L117): ifor is complient, gfortran not
 FFLAGS2 := -ffast-math -funroll-loops -O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fopenmp 
 LDFLAGS := -static
 LDFLAGS2 := 
 EXEC := g_baypass
 ifeq ($(DEBUG),Yes)
  DEBUGFLAGS := -g -Wconversion
 else
  DEBUGFLAGS :=
 endif
else
 FFLAGS :=  -O3 -r8 -traceback
 FFLAGS2 := -O3 -r8 -qopenmp -nowarn -traceback
 LDFLAGS := -static
 LDFLAGS2 := 
 EXEC := i_baypass
 ifeq ($(DEBUG),Yes)
  DEBUGFLAGS := -g -traceback
 else
  DEBUGFLAGS :=
 endif
endif

all: $(EXEC)

$(EXEC): baypass.f90 mt_stream.o f_get_coeff.o gf2xe.o
	$(FC) $(FFLAGS2) $(LDFLAGS2) $(DEBUGFLAGS) $^ -o $@

f_get_coeff.o: f_get_coeff.F90 gf2xe.o f_jump_coeff.o
	$(FC) $(FFLAGS) $(LDFLAGS) $(DEBUGFLAGS) -c $< -o $@

gf2xe.o: gf2xe.F90
	$(FC) $(FFLAGS) $(LDFLAGS) $(DEBUGFLAGS) -c $< -o $@

f_jump_coeff.o: f_jump_coeff.F90
	$(FC) $(FFLAGS) $(LDFLAGS) $(DEBUGFLAGS) -c $< -o $@

mt_stream.o: mt_stream.F90 mt_kind_defs.o
	$(FC) $(FFLAGS) $(LDFLAGS) $(DEBUGFLAGS) -c $< -o $@

mt_kind_defs.o: mt_kind_defs.F90
	$(FC) $(FFLAGS) $(LDFLAGS) $(DEBUGFLAGS) -c $< -o $@

clean:
	rm -f *.o *.mod
