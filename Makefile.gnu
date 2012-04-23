FFLAGC  =  -g -O2 -fpp -fopenmp -DHORNOS
opt     = 
profile = 
obj     = lev00
ildir   =  
input   = 

CFLAGC  = -g -O2 -fpp -fopenmp -DHORNOS

FCOMPL  = gfortran 
CCOMPL  = gcc 
# list of other directories for source files
.PREFIXES: .

.SUFFIXES:
.SUFFIXES: .f90 .c .s .o .fil 

.f90.o:
	$(FCOMPL) -c $(ildir) $(FFLAGC) $(opt) $(profile) $<

.f.fil:
	$(FCOMPL) -il $(FFLAGC) $<

.s.o:
	as $<

.c.o:
	$(CCOMPL) -c $(CFLAGC) $(profile) $<

#OBJ_mod =  atoms.o menu.o  box.o kpoints.o mendeleev.o \
#              code.o param.o dos_inc.o explore.o siesta_eig.o  

OBJECTS = atoms.o menu.o  box.o kpoints.o mendeleev.o \
code.o param.o dos_inc.o explore.o  siesta_eig.o \
device.o bastr.o do_param.o get_param_siesta.o \
hat.o read_siesta_input.o read_vasp_geom.o lev00.o \
tools.o tools_strings.o prep_dos.o dos_add.o \
mainmenu.o read_vasp_psi2.o plotting.o prep_disp.o \
dipole.o density.o write_dens.o lev_coulmb.o hould.o \
invers.o manip_dens.o plot_add.o read_density.o \
simulate.o read_siesta_pdos.o ttag.o choose_tasks.o stm_TH.o

INLINE  =

APPLIC: $(INLINE) $(OBJECTS)
	$(FCOMPL) $(FFLAGC) $(profile) -o $(obj) $(OBJECTS) $(LIBS)

test:
	@echo START TEST ON $(input) , opt = $(opt)
	@echo start test on $(input) , opt = $(opt) >> TIME.LOG
	@date >> TIME.LOG
	@( time $(obj) < $(input) > $(input).out ) 2>> TIME.LOG
	@echo - - - - - - - - - - - >> TIME.LOG

clean:
	@rm -f $(INLINE) $(OBJECTS) $(OBJ_mod) *.mod

# include file dependencies

lev00.o mainmenu.o :: param.f90 atoms.f90 code.f90 kpoints.f90 dos_inc.f90
device.o dipole.o simulate.o :: param.f90 atoms.f90 menu.f90
do_param.o get_param_siesta.o write_dens.o read_density.o read_siesta_pdos.o siesta_eig.o :: param.f90 \
                                        atoms.f90 code.f90 
read_siesta_input.o read_vasp_geom.o :: param.f90 atoms.f90 code.f90 \
                                        kpoints.f90 mendeleev.f90
prep_dos.o :: param.f90 atoms.f90 kpoints.f90 dos_inc.f90
dos_add.o :: param.f90 dos_inc.f90
read_vasp_psi2.o prep_disp.o :: param.f90 atoms.f90 dos_inc.f90
density.o :: param.f90 atoms.f90 menu.f90 box.f90
lev_coulmb.o :: param.f90 atoms.f90 menu.f90 code.f90
manip_dens.o :: param.f90 code.f90
plot_add.o :: param.f90 menu.f90 explore.f90
choose_tasks.o :: param.f90 dos_inc.f90 atoms.f90
stm_TH.o :: param.f90 atoms.f90 menu.f90
