# Makefile script for Tree Smoothed-particle-hydrodynamics codes.
#
# Author: Eraldo Pereira Marinho, PhD. at University of Sao Paulo, Brazil
#
# May 17th, 1999
#
# Original Platform: Digital Unix (OSF/1)

# Fortran 77 flags:
# FC=kf77
#FFLAGS= -ansi -extend_source -cpp -check nobounds $(DIRECTIVES) -V \
#	-cross_reference
FFLAGS=  $(DIRECTIVES) -ffixed-line-length-132 -ansi

# C-preprocessor:
CPP=/lib/cpp

# C-compiler:
# CC=gcc

# C-flags:
#CFLAGS= -ansi $(DIRECTIVES) -arch host -protect_headers all -v -w0 -fast
CFLAGS= -ansi $(DIRECTIVES) 

# Directories definitions:
MYBIN=$(HOME)/bin
MYLIB=$(HOME)/lib
MYSHARE=$(HOME)/share
MYSHLIB=$(HOME)/shlib
MYINCLUDE=$(HOME)/include
DISTDIR=.

# Main entry point:
SPH_MAIN=sph_main.o

# Base Compiled Objects:
BASE_OBJS= \
sph_advance.o \
sph_gravity.o \
sph_integral.o \
sph_intschm.o \
sph_io.o \
sph_control.o \
sph_sseeds.o \
sph_timing.o \
sph_tree.o

# MHD Subroutines:
MHD_OBJS= sph_mhd.o

#SPH Subroutines:
SPH_OBJS= \
sph_usr.o \
sph_usr_accel.o \
sph_nns.o \
sph_kernel.o \
sph_thermo.o

# Others:
OTHERS= amax.o amin.o sph_io.o

# SPH Common Data Structures
SPH_COMMON= sph_common

# Cooling Subroutines:
COOL_OBJS = \
Cool_interp.o \
sph_cool.o

# Thermodynamics objects:
THERMO = \
e_h2.o \
w_mol.o \
thermo.o

# Prefix name:
PREFIX=SPH
#
# Extra objects:
EXTRA_OBJS=

# Definite objects:
EXEC_DEPENDENCES=\
common \
$(SPH_MAIN) \
$(OTHERS) \
$(SPH_OBJS) \
$(BASE_OBJS) \
$(THERMO) \
$(EXTRA_OBJS)

# Targets:
#
nothing:
	@echo Makefile: no rule.
	@echo '	try make <instance>'
exec: $(EXEC_DEPENDENCES)
	$(FC) 	-O4 \
		-ldran \
		-o $(MYBIN)/$(PREFIX)_$(SUFFIX) \
		$(SPH_MAIN) \
		$(EXTRA_OBJS) \
		$(OTHERS) $(SPH_OBJS) \
		$(BASE_OBJS) \
		$(THERMO)
adiabatic: 
	make exec \
		'EXTRA_OBJS=$(COOL_OBJS)' \
		'SUFFIX=adiabatic' \
		'DIRECTIVES=\
		-D_VERBOSE \
		-D_MOVIE_ \
		-D_40960 \
		-D_time_resol_3 \
		-D_mid_neighbor_resol \
		-D_H_1 \
		-D_SPH_ADIABAT \
		-D_LINUX_ \
		'
adiabatic2: 
	make exec \
		'EXTRA_OBJS=$(COOL_OBJS)' \
		'SUFFIX=adiabatic2' \
		'DIRECTIVES= \
		-D_MOVIE_ \
		-D_262144 \
		-D_time_resol_3 \
		-D_mid_neighbor_resol \
		-D_H_1 \
		-D_LINUX_ \
		-D_SPH_ADIABAT'
dissp: 
	make exec \
		'EXTRA_OBJS=$(COOL_OBJS)' \
		'SUFFIX=dissp' \
		'DIRECTIVES= \
		-D_MOVIE_ \
		-D_LINUX_ \
		-D_40960 \
		-D_time_resol_3 \
		-D_mid_neighbor_resol \
		'
dissp2: 
	make exec \
		'EXTRA_OBJS=$(COOL_OBJS)' \
		'SUFFIX=dissp' \
		'DIRECTIVES=\
		-D_MOVIE_ \
		-D_262144 \
		-D_time_resol_3 \
		-D_LINUX_ \
		-D_mid_neighbor_resol \
		'
mhd-dissp: 
	make exec \
		'EXTRA_OBJS=$(MHD_OBJS) $(COOL_OBJS)' \
		'SUFFIX=mhd-dissp' \
		'DIRECTIVES= \
		-D_SPMHD \
		-D_MOVIE_ \
		-D_40960 \
		-D_time_resol_3 \
		-D_mid_neighbor_resol \
		-D_LINUX_ \
		'
mhd-dissp2: 
	make exec \
		'EXTRA_OBJS=$(MHD_OBJS) $(COOL_OBJS)' \
		'SUFFIX=mhd-dissp2' \
		'DIRECTIVES= \
		-D_SPMHD \
		-D_BACKGROUND_ \
		-D_MOVIE_ \
		-D_262144 \
		-D_time_resol_5 \
		-D_mid_neighbor_resol \
		-D_LINUX_ \
		'
mhd-ad: 
	make exec \
		'EXTRA_OBJS=$(MHD_OBJS)' \
		'SUFFIX=mhd-ad' \
		'DIRECTIVES= \
		-D_40960 \
		-D_MOVIE_ \
		-D_SPMHD \
		-D_H_2 \
		-D_LINUX_ \
		-D_SPH_ADIABAT \
		-D_BACKGROUND_ \
		-D_time_resol_3 \
		-D_mid_neighbor_resol \
		'
galaxy-2d:
	make exec \
		'EXTRA_OBJS=$(COOL_OBJS)' \
		'SUFFIX=galaxy-2d' \
		'DIRECTIVES= \
		-D_40960 \
		-D_2D \
		-D_H_1 \
		-D_MOVIE_ \
		-D_time_resol_7 \
		-D_GALAXY_2D \
		-D_SPIRAL \
		-D_SLOW_SWITCH_ON \
		-D_LINUX_ \
		-D_ROTATING_FRAME \
		-D_NON_SELFGRAVITATING \
		-D_COMPUTE_PHYSICAL_TIME \
		-D_BACKGROUND_ \
		-D_NEW_METHOD \
		-D_mid_neighbor_resol \
		'
nbody: 
	make exec \
		'PREFIX=tree' \
		'SUFFIX=nbody' \
		'THERMO=' \
		'DIRECTIVES= \
		-D_CHECK_TREE_SINGULARITY \
		-D_VERBOSE \
		-D_MOVIE_ \
		-D_N_BODY \
		-D_LINUX_ \
		-D_262144 \
		-D_time_resol_4 \
		'
common: $(SPH_COMMON).h
$(SPH_COMMON).h:
	$(CPP) $(DIRECTIVES) $(SPH_COMMON).F > $(SPH_COMMON).h
clean:
	$(RM) $(RMFLAGS) *.o *.l *.out *.cmp.f $(SPH_COMMON).h *.aux *.log *.dvi
Cooling_Table: Cooling_Table.o CO.o cooling.o numstring.o H-H2-HM.o 
	$(FC) -fast -o cooling_table \
	cooling_Table.o co.o cooling.o numstring.o h-h2-hm.o
listobj:
	@echo $(OTHERS) $(SPH_OBJS) $(BASE_OBJS) $(THERMO) $(MHD_OBJS)
vortex:=DIRECTIVES= \
		-D_SPMHD -D_H_1 -D_SPH_ADIABAT \
		-D_40960 \
		-D_mid_neighbor_resol \
		-D_LINUX_ \
		-D_time_resol_3 
vortex: common sph_vortex.o sph_tree.o sph_kernel.o sph_io.o \
	sph_mhd.o sph_usr.o sph_nns.o $(OTHERS) sph_gravity.o 
	$(FC) -fast -o $(MYBIN)/sph_vortex sph_vortex.o sph_tree.o \
					   sph_kernel.o sph_io.o \
					   sph_mhd.o sph_usr.o sph_nns.o \
					   $(OTHERS) sph_gravity.o
analysis:=DIRECTIVES= \
		-D_40960 \
		-D_time_resol_3 \
		-D_mid_neighbor_resol \
		-D_LINUX_ \
		-D_H_1 \
		-D_SPH_ADIABAT

analysis:
	make	exec \
		'SPH_MAIN=sph_analysis.o' \
		'EXTRA_OBJS=$(COOL_OBJS)' \
		'SUFFIX=analysis' \
		'DIRECTIVES=$(DIRECTIVES)'
incompress:=DIRECTIVES= \
		-D_VRM_INCOMPRESSIBLE_SPH \
		-D_LINUX_ \
		-D_MOVIE_ \
		-D_40960 \
		-D_time_resol_3 \
		-D_mid_neighbor_resol \
		-D_H_1 \
		-D_SPH_ADIABAT
incompress: 
	make exec \
		'EXTRA_OBJS=$(COOL_OBJS)' \
		'SUFFIX=incompress' \
		'DIRECTIVES=$(DIRECTIVES)'
base=iso_g

source=$(base).f

object=$(base).o

exec=$(base)

iso_g: $(exec)

$(exec): $(object)
	$(FC)  -O4 -o $(exec) $(object)

$(object):
	$(FC) -c -o $(object) -O4 $(source)
fractal:=DIRECTIVES=\
                -D_SPMHD -D_H_1 -D_SPH_ADIABAT \
		-D_LINUX_ \
                -D_262144 \
                -D_time_resol_3 \
                -D_mid_neighbor_resol
fractal: sph_common.h sph_tree.o sph_fractal.o sph_nns.o amax.o $(MHD_OBJS)
	$(FC) $(DIRECTIVES) -fast -o $(MYBIN)/sph_fractal \
			sph_tree.o \
			sph_fractal.o \
			sph_nns.o \
			amax.o \
			$(MHD_OBJS) 
Units: Units.tex Units.ps
	latex units.tex
	latex units.tex
	latex units.tex
Units.ps: Units.dvi
	dvips -o units.ps units
core:=DIRECTIVES= \
		-D_SPH_ADIABAT \
		-D_SIMPLER_DISK_MODEL \
		-D_BACKGROUND_ \
		-D_NON_SELFGRAVITATING \
		-D_MOVIE_ \
		-D_40960 \
		-D_2D \
		-D_H_1 \
		-D_LINUX_ \
		-D_time_resol_7 \
		-D_NEW_METHOD \
		-D_mid_neighbor_resol
core:
	make exec \
		'EXTRA_OBJS=$(COOL_OBJS)' \
		'SUFFIX=core' \
		'DIRECTIVES=$(DIRECTIVES)'
galaxy-2d-relax:=DIRECTIVES= \
		-D_40960 \
		-D_2D \
		-D_H_1 \
		-D_MOVIE_ \
		-D_LINUX_ \
		-D_time_resol_7 \
		-D_NON_SELFGRAVITATING \
		-D_GALAXY_2D \
		-D_NON_ROTATING_FRAME \
		-D_NON_SPIRAL \
		-D_NEW_METHOD \
		-D_COMPUTE_PHYSICAL_TIME \
		-D_mid_neighbor_resol
galaxy-2d-relax:
	make exec \
		'EXTRA_OBJS=$(COOL_OBJS)' \
		'SUFFIX=galaxy-2d-relax' \
		'DIRECTIVES=$(DIRECTIVES)'
galaxy-2d-adiabatic:=DIRECTIVES= \
		-D_40960 \
		-D_2D \
		-D_H_1 \
		-D_MOVIE_ \
		-D_time_resol_7 \
		-D_GALAXY_2D \
		-D_SPIRAL \
		-D_ROTATING_FRAME \
		-D_NON_SELFGRAVITATING \
		-D_LINUX_ \
		-D_COMPUTE_PHYSICAL_TIME \
		-D_SLOW_SWITCH_ON \
		-D_SPH_ADIABAT \
		-D_BACKGROUND_ \
		-D_mid_neighbor_resol
galaxy-2d-adiabatic:
	make exec \
		'EXTRA_OBJS=$(COOL_OBJS)' \
		'SUFFIX=galaxy-2d-adiabatic' \
		'DIRECTIVES=$(DIRECTIVES)'
disk3d:=DIRECTIVES= \
		-D_SIMPLER_DISK_MODEL \
		-D_MOVIE_ \
		-D_40960 \
		-D_LINUX_ \
		-D_3D \
		-D_time_resol_7 \
		-D_mid_neighbor_resol
disk3d:
	make exec \
		'EXTRA_OBJS=$(COOL_OBJS)' \
		'SUFFIX=disk3d' \
		'DIRECTIVES=$(DIRECTIVES)'
nbody-doc:
	@echo no rule
