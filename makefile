TARGET = simple3D

SRC =  modules.f90 3D.f90 simple3D.f90 allocation.f90  allocatecpp.f90  cadenas.f90  fe.f90  fkfun.f90  kai.f90  kinsol.f90  pxs.f90  savetodisk.f90 rands.f90 dielectric.f90 ct.f90

HOST=$(shell hostname)
$(info HOST is ${HOST})


ifeq ($(HOST),carapa)
LFLAGS = -lsundials_fkinsol -lsundials_fnvecserial -lsundials_kinsol -lsundials_nvecserial -lm
endif

ifeq ($(HOST),juno)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif

ifeq ($(HOST),spinetta)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif



ifeq ($(HOST),leomisso)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif


# some definitions
SHELL = /bin/bash
FFLAGS= -O3 # -fbacktrace -fbounds-check # -O3

LFLAGS =  -L/opt/local/sundials-2.5.0-openmpi-atlas/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1 -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../.. -lgfortran -lm -lgcc_s -lquadmath

ifeq ($(HOST),piluso.rosario-conicet.gov.ar)
LFLAGS = -L/home/mtagliazucchi.inquimae/software/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s
endif


ifeq ($(HOST),pear)
LFLAGS=-L/home/mario/software/KINSOL/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux
endif

ifeq ($(HOST),cnode01)
LFLAGS = -L/home/mtagliazucchi/software/Kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s
endif

ifeq ($(HOST),master) 
LFLAGS = -L/shared/software/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath
endif

ifeq ($(HOST),mate.bme.northwestern.edu) 
LFLAGS = -L/home/mario/software/kinsol/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath
endif

ifeq ($(HOST),quser13) 
LFLAGS = -L/home/mta183/KINSOL2.7/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lgfortranbegin -lgfortran -lm
endif

ifeq ($(HOST),quser12)
LFLAGS = -L/home/mta183/KINSOL2.7/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lgfortranbegin -lgfortran -lm
endif

ifeq ($(HOST),quser11)
LFLAGS = -L/home/mta183/KINSOL2.7/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lgfortranbegin -lgfortran -lm
endif

ifeq ($(HOST),quser10)
LFLAGS = -L/home/mta183/KINSOL2.7/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2 -L/usr/lib/gcc/x86_64-redhat-linux/4.1.2/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -lgfortranbegin -lgfortran -lm
endif

GIT_VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)
GFLAGS=-cpp -D_VERSION=\"$(GIT_VERSION)\"

FF = /opt/local/openmpi-1.8.8/bin/mpif77 #${F90}
VER = ~/bin/crystal

all:	$(TARGET)

$(TARGET): $(SRC:.f90=.o)
	$(FF) -o $(TARGET) $(SRC:.f90=.o) $(LFLAGS) $(GFLAGS)
	cp $(TARGET) $(VER)

$(SRC:.f90=.o): $(SRC)
	${FF} -c ${FFLAGS}  $(SRC) $(LFLAGS) $(GFLAGS)

install: all
	cp $(TARGET) $(VER)

clean:	
	@rm -f $(SRC:.f90=.o) $(SRC:.f90=.d) $(TARGET) *~

realclean: clean
	@rm -f .depend

depend dep:
	@$(FF)  $(CFLAGS) -MM $(SRC) > .depend 

ifeq (.depend, $(wildcard .depend))
include .depend
endif

















































