subroutine allocation

use system
use fields_fkfun
use conformations
use chainsdat
use kinsol
use results
use mkinsol
implicit none

! fields_fkfun
ALLOCATE(xtotal(0:dimx+1, 0:dimy+1, -1:dimz+2)) ! xtotal para poor solvent
ALLOCATE(psi(0:dimx+1, 0:dimy+1, 0:dimz+1))
ALLOCATE(q(ncha))
ALLOCATE(xh(dimx, dimy, dimz))

! kinsol
ALLOCATE (xflag(2*dimx*dimy*dimz))

! results
ALLOCATE (avpol(dimx, dimy, dimz))
ALLOCATE (xpos(dimx, dimy, dimz)) ! pos ion
ALLOCATE (xneg(dimx, dimy, dimz)) ! neg ioni
ALLOCATE (qtot(dimx, dimy, dimz)) ! Carga total
ALLOCATE (xHplus(dimx, dimy, dimz)) ! H+
ALLOCATE (xOHmin(dimx, dimy, dimz)) ! OH-
ALLOCATE (fdis(dimx, dimy, dimz))
ALLOCATE (epsfcn(0:dimx+1, 0:dimy+1, 0:dimz+1))
ALLOCATE (Depsfcn(0:dimx+1, 0:dimy+1, 0:dimz+1))

! mkinsol
ALLOCATE (pp(2*dimx*dimy*dimz))

! chainsdat
allocate(posicion(ncha,2))
allocate(in1(long,3))

end subroutine
