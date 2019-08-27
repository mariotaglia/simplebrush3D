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
!ALLOCATE (xflag(3*dimx*dimy*dimz))
! results
ALLOCATE (avpol(dimx, dimy, dimz,2))
ALLOCATE (xpos(dimx, dimy, dimz)) ! pos ion
ALLOCATE (xneg(dimx, dimy, dimz)) ! neg ioni
ALLOCATE (qtot(dimx, dimy, dimz)) ! Carga total
ALLOCATE (xHplus(dimx, dimy, dimz)) ! H+
ALLOCATE (xOHmin(dimx, dimy, dimz)) ! OH-
ALLOCATE (fdis(dimx, dimy, dimz))
ALLOCATE (epsfcn(0:dimx+1, 0:dimy+1, 0:dimz+1))
ALLOCATE (Depsfcn(0:dimx+1, 0:dimy+1, 0:dimz+1))

!!G:!!!!!!!!

allocate(fdisANC(dimx, dimy, dimz)) !	fraction not charge pol-A
allocate(fdisBNC(dimx, dimy, dimz)) !	fraction not charge pol-B
allocate(fdisAas(dimx, dimy, dimz)) !	fraction associate pol-A
allocate(fdisBas(dimx, dimy, dimz)) 
allocate(fdisANa(dimx, dimy, dimz)) !	fraction not charge pol-A
allocate(fdisBCl(dimx, dimy, dimz)) !	fraction not charge pol-B
allocate(xna(dimx, dimy, dimz))
allocate(xnb(dimx, dimy, dimz))
allocate(eta(dimx, dimy, dimz))
allocate(m(dimx, dimy, dimz))
allocate(KK0check(dimx, dimy, dimz))
allocate(KKaAcheckplus(dimx, dimy, dimz))
allocate(KKaAna(dimx, dimy, dimz))
allocate(KKaBCl(dimx, dimy, dimz))
allocate(KKaBcheckmin(dimx, dimy, dimz))
!!G:!!!!!!!

! mkinsol
ALLOCATE (pp(2*dimx*dimy*dimz))

! chainsdat
allocate(posicion(ncha,2))
allocate(ct(ncha))
allocate(in1(long,3))

end subroutine
