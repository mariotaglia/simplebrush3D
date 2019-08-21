TARGET = 3D

SRC = simple3D.f90 modules.f90 MCMT.f90 allocation.f90 3D.f90  allocatecpp.f90  cadenas.f90  fe.f90  fkfun.f90  kai.f90  kinsol.f90  pxs.f90  savetodisk.f90 rands.f90 dielectric.f90

# some definitions
SHELL = /bin/bash
FFLAGS= -O3 

ifeq ($(HOST),spinetta)
LFLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
LDLAGS = -lm /usr/lib/x86_64-linux-gnu/librt.so  -L/usr/local/lib  -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial ${LIBS} -Wl,-rpath,/usr/local/lib
endif

FF = mpif77 #${F90}
VER = ~/bin

all:	$(TARGET)

$(TARGET): $(SRC:.f90=.o)
	$(FF) -o $(TARGET) $(SRC:.f90=.o) $(LFLAGS) $(LDFLAGS)
$(SRC:.f90=.o):
	${FF} -c ${FFLAGS}  $(SRC) $(LFLAGS) $(LDFLAGS)

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

















































