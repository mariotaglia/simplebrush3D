module system 
real*8, parameter :: delta = 0.5
integer  dimx 
integer  dimy 
integer  dimz 
real*8 gama 
endmodule

module mkinsol
double precision, allocatable :: pp(:)
endmodule

module chainsdat
integer cuantas 
integer long 
integer ncha 
real*8, ALLOCATABLE :: in1(:,:)  ! segment positions 
real*8, ALLOCATABLE :: posicion(:,:) ! posicion graft de la cadena ncha
integer cpp
endmodule

module molecules
use system
real*8 vsol
real*8 vpol
real*8 vsalt
real*8 zpos,zneg, zpol
real*8 st
endmodule

module fields_fkfun
use system
use chainsdat
real*8, allocatable :: xtotal(:, :, :) ! xtotal para poor solvent
real*8, allocatable :: psi(:, :, :) 
real*8, allocatable :: q(:)
!pro(cuantas, cpp)
real*8, allocatable :: pro(:,:)
real*8, allocatable :: xh(:, :, :)
real*8 shift
endmodule

module conformations
integer, allocatable :: px(:,:,:)
integer, allocatable :: py(:,:,:)
integer, allocatable :: pz(:,:,:)
endmodule

module MPI
include 'mpif.h' ! librerias MPI
integer rank, size, ierr
integer flagsolver
endmodule

module kinsol
use system
integer iter
integer *4 ier ! Kinsol error flag
integer *8 neq ! Kinsol number of equations
real*8 norma
real*8, ALLOCATABLE :: xflag(:) 
endmodule

module const
real*8 dielW, dielP, dielS
real*8 constqE
real*8 dielPr, dielSr
real*8 pKw
real*8 pi 
real*8, parameter :: Na = 6.02d23 
real*8 constq
real*8 lb
integer seed
real*8 error  ! para comparar con la norma...
real*8 errel
integer itmax
integer infile
integer randominput
integer verbose
endmodule

module kai
real*8 Xu(-2:2, -2:2, -2:2)
endmodule

module results
use system
real*8, allocatable :: avpol(:,:,:)
real*8, allocatable :: epsfcn(:,:,:)
real*8, allocatable :: Depsfcn(:,:,:)
real*8, allocatable :: xpos(:,:,:) ! pos ion
real*8, allocatable :: xneg(:,:,:) ! neg ioni
real*8, allocatable :: qtot(:,:,:) ! Carga total
real*8, allocatable :: xHplus(:,:,:) ! H+
real*8, allocatable :: xOHmin(:,:,:) ! OH-
real*8, allocatable :: fdis(:,:,:)
endmodule

module bulk
real*8 expmupos,expmuneg,expmuHplus,expmuOHmin
real*8 K0
real*8 xsolbulk, xposbulk, xnegbulk, xHplusbulk,xOHminbulk
endmodule

module inputtemp
real*8 xsalt
real*8 pHbulk
real*8 pOHbulk
real*8 pKa, Ka
real*8 csalt
real*8 cHplus, cOHmin
end module
