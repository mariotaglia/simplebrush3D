subroutine fkfun(x,f,ier2)

use system
use chainsdat
use molecules
use const
use results
use bulk
use kai
use MPI
use fields_fkfun
use kinsol
use conformations
implicit none

integer*4 ier2
integer ntot
real*8 x(*),f(*)
real*8 protemp
integer i,j, ix, iy, iz, ii, ax, ay, az
integer jx, jy, jj
real*8 xpot(dimx, dimy, dimz,2)
! Charge
real*8 psitemp
! poor solvent 
real*8 sttemp
! MPI
integer tag
parameter(tag = 0)
integer err
real*8 avpol_tosend(dimx,dimy,dimz,2)
real*8 avpol_temp(dimx,dimy,dimz)
real*8 volprot(dimx,dimy,dimz)
real*8 q_tosend
real*8 gradpsi2

!-----------------------------------------------------
! Common variables

shift = 1.0
volprot = 0.0

! Jefe

if(rank.eq.0) then ! llama a subordinados y pasa vector x
   flagsolver = 1
   CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)
   CALL MPI_BCAST(x, 2*dimx*dimy*dimz , MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)
endif

!------------------------------------------------------
! DEBUG
!      if(iter.gt.2000) then
!      do i = 1, n
!      print*,i, x(i)
!      enddo
!      endif


! Recupera xh y psi desde x()

ntot = dimx*dimy*dimz ! numero de celdas
do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
     xh(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1))
     psi(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ntot)   
     xna(ix,iy,iz)=x(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ntot+ntot)   
  enddo
 enddo
enddo
      
! Condiciones de borde potencial electrostatico
! en x
 
do iy = 1, dimy
 do iz = 1, dimz
   psi(0, iy, iz) = psi(dimx, iy, iz)
   psi(dimx+1, iy, iz) = psi(1, iy, iz)
 enddo
enddo

! en y
do ix = 1, dimx
 do iz = 1, dimz
   psi(ix, 0, iz) = psi(ix, dimy, iz)
   psi(ix, dimy+1, iz) = psi(ix, 1, iz)
 enddo
enddo

! en z
do ix = 1, dimx
 do iy = 1, dimy
   psi(ix, iy, dimz+1) = 0.0  ! psibulk = 0.0
   psi(ix, iy, 0) = psi(ix, iy, 1) ! zero charge
!   psi(ix, iy, 0) = psi(ix, iy, 1)*epsfcn(ix,iy,1)/epsfcn(ix,iy,0) ! zero charge
 enddo
enddo

! aristas... importantes para lattices no cubicos...
do iz = 1, dimz
 psi(0, 0, iz) = psi(dimx, dimy, iz)
 psi(dimx+1, dimy+1, iz) = psi(1, 1, iz)
 psi(dimx+1, 0, iz) = psi(1, dimy, iz)
 psi(0, dimy+1, iz) = psi(dimx, 1, iz)
enddo

! No importan las esquinas....
! Fracciones de volumen inicial	y fdis

do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
!    avpol(ix,iy,iz)=0.0	
    xpos(ix, iy, iz) = expmupos*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zpos) ! ion plus volume fraction 
    xneg(ix, iy, iz) = expmuneg*(xh(ix, iy, iz)**vsalt)*dexp(-psi(ix, iy, iz)*zneg) ! ion neg volume fraction
    xHplus(ix, iy, iz) = expmuHplus*(xh(ix, iy, iz))*dexp(-psi(ix, iy, iz))           ! H+ volume fraction
    xOHmin(ix, iy,iz) = expmuOHmin*(xh(ix,iy,iz))*dexp(+psi(ix,iy,iz))           ! OH-  volume fraction
    !fdis(ix,iy,iz)=1.0 /(1.0 + xHplus(ix,iy,iz)/(K0*xh(ix,iy,iz)) )

!!G::
	xnb(ix,iy,iz)=1.0-xna(ix,iy,iz)-xh(ix,iy,iz)-xpos(ix,iy,iz)-xneg(ix,iy,iz)-xHplus(ix,iy,iz)-xOHmin(ix,iy,iz)
	fdisAas(ix,iy,iz)=0.0d0 													
	fdisBas(ix,iy,iz)=0.0d0			
!
if ((1.0d-6 .lt. xna(ix,iy,iz)).AND.(1.0d-6 .lt. xnb(ix,iy,iz))) THEN
	eta(ix,iy,iz)=xna(ix,iy,iz)/xnb(ix,iy,iz)

M(ix,iy,iz)=( 1.0+ (xOHmin(ix,iy,iz))/(K0B*xh(ix,iy,iz))+(xneg(ix,iy,iz)/(K0BCl*xh(ix,iy,iz)**vsalt)) )*( 1.0&
+ (xHplus(ix,iy,iz))/(K0A*xh(ix,iy,iz)) +(xpos(ix,iy,iz)/(K0ANa*xh(ix,iy,iz)**vsalt)))/(K0Eo*xna(ix,iy,iz)) 							   !!gabi: vpair = vpol!!
 	 	fdisAas(ix,iy,iz) = (1.0+eta(ix,iy,iz)+eta(ix,iy,iz)*M(ix,iy,iz))/(2.0*eta(ix,iy,iz))-(	-1.0/eta(ix,iy,iz)+(	(1.0+eta(ix,iy,iz)+eta(ix,iy,iz)*M(ix,iy,iz))/(2.0*eta(ix,iy,iz))	)**2	)**0.5 !!!!!!!!!!!!!!!!!!!!!!!
   	fdisBas(ix,iy,iz) = (1.0+eta(ix,iy,iz)+eta(ix,iy,iz)*M(ix,iy,iz))/(2.0)-eta(ix,iy,iz)*(	-1.0/eta(ix,iy,iz)+(	(1.0+eta(iz)+eta(iz)*M(ix,iy,iz))/(2.0*eta(ix,iy,iz))	)**2	)**0.5 !!!!!!!!!!!!!!!!!!!!!!!
	
	endif

	fdisANC(ix,iy,iz) = (1.0-fdisAas(ix,iy,iz) )/(1.0 + (K0A*xh(ix,iy,iz))/(xHplus(ix,iy,iz))&
+ (K0A*xpos(ix,iy,iz)* (xh(ix,iy,iz)**(1.0-vsalt)))/(K0ANa*xHplus(ix,iy,iz)) )						   !g
   fdisBNC(ix,iy,iz) = (1.0-fdisBas(ix,iy,iz) )/(1.0 + (K0B*xh(ix,iy,iz))/(xOHmin(ix,iy,iz))&
 + (K0B*xneg(ix,iy,iz)* (xh(ix,iy,iz)**(1.0-vsalt)))/(K0BCl*xOHmin(ix,iy,iz)) )						   !g
	fdisANa(ix,iy,iz) = (fdisANC(ix,iy,iz))*(K0A*xpos(ix,iy,iz)* (xh(ix,iy,iz)**(1.0-vsalt)) )/(xHplus(ix,iy,iz)*K0ANa)
	fdisBCl(ix,iy,iz) = (fdisBNC(ix,iy,iz))*(K0B*xneg(ix,iy,iz)* (xh(ix,iy,iz)**(1.0-vsalt)) )/(xOHmin(ix,iy,iz)*K0BCl)

			!KK0check(ix,iy,iz)=-dlog10( (Na/1.0d24)*fdisBas(ix,iy,iz)/(	(1.0-fdisAas(ix,iy,iz)-fdisANa(ix,iy,iz)&
!-fdisANC(ix,iy,iz))*(1.0-fdisBas(ix,iy,iz)-fdisBNC(ix,iy,iz)-fdisBCl(ix,iy,iz))*xna(ix,iy,iz) )	)/pKeo
		!	KKaAcheckplus(ix,iy,iz)= -dlog10( (xHplus(ix,iy,iz)/xh(ix,iy,iz))*((1-fdisANC(ix,iy,iz)-fdisANa(ix,iy,iz)&
!-fdisAas(ix,iy,iz))/fdisANC(ix,iy,iz))*(xsolbulk*1.0d24/(Na*vsol)))-pKaA !! esto era para chequear pkaA
		!	kkaBcheckmin(ix,iy,iz)=	 (xOhmin(ix,iy,iz)/xh(ix,iy,iz))*(1.0-fdisBas(ix,iy,iz)-fdisBCl(ix,iy,iz)-fdisBNC(ix,iy,iz))/fdisBNC(ix,iy,iz)-K0B
		!	KKaANa(ix,iy,iz)= dlog10( (xh(ix,iy,iz)/xpos(ix,iy,iz))*(fdisANa(ix,iy,iz)/(1-fdisANC(ix,iy,iz)-fdisANa(ix,iy,iz)&
!-fdisAas(ix,iy,iz)))*(xsolbulk*1.0d24/(Na*vsol)))/pKaAna !! esto era para chequear pkaA
		!	KKaBCl(ix,iy,iz)= dlog10( (xh(ix,iy,iz)/xneg(ix,iy,iz))/((1-fdisBNC(ix,iy,iz)-fdisBCl(ix,iy,iz)&
!-fdisBas(ix,iy,iz))/fdisBCl(ix,iy,iz))*(xsolbulk*1.0d24/(Na*vsol)))/pkaBCl !! esto era para chequear pkaA


!!G:
   enddo
 enddo  
enddo

! Calculo de xtotal para poor solvent
! en el lattice
do ix=1,dimx
 do iy=1,dimy
  do iz=1,dimz
   xtotal(ix, iy, iz) = 1.0-xpos(ix, iy,iz) - xneg(ix, iy, iz)- xh(ix, iy, iz) - xHplus(ix, iy, iz) - xOHmin(ix, iy, iz) ! xtotal es todo menos solvente e iones
  enddo
 enddo
enddo

do ix = 1, dimx
 do iy = 1, dimy
  xtotal(ix, iy, dimz+1) = 0.0 ! xtotal en bulk = 0.0
  xtotal(ix, iy, dimz+2) = 0.0 ! xtotal en bulk = 0.0
  xtotal(ix, iy, 0) = 0.0 ! xtotal en la superficie = 0.0
  xtotal(ix, iy, -1) = 0.0 ! xtotal en la superficie = 0.0
 enddo
enddo

do ix = 1, dimx
 do iz = 1, dimz
  xtotal(ix,0,iz) = xtotal(ix,dimy,iz)
  xtotal(ix,dimy+1,iz) = xtotal(ix,1,iz)
 enddo
enddo

do ix = 1, dimx
 do iy = 1, dimy
  xtotal(ix, iy, dimz+1) = 0.0 ! xtotal en bulk = 0.0
  xtotal(ix, iy, 0) = 0.0 ! xtotal en la superficie = 0.0
 enddo
enddo

! Compute dielectric permitivity
call dielectfcn(xtotal,volprot,epsfcn,Depsfcn)

!------------------------------------------------------------------------
! PDFs polimero
!------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! PARALELO: Cada procesador trabaja sobre una cadena...
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calcula xpot

sttemp = st*vpol*vsol

do ix=1,dimx
 do iy=1,dimy
   do iz=1,dimz

    ! xpot(ix, iy, iz) = xh(ix,iy,iz)**vpol/fdis(ix,iy,iz)*dexp(-psi(ix, iy, iz)*zpol)
   
  xpot(ix, iy, iz,1) = xh(ix, iy, iz)**vpol*dexp(-psi(ix, iy, iz)*zpolA)/(1.0-fdisAas(ix, iy, iz)-fdisANC(ix, iy, iz)-fdisANa(ix, iy, iz)) 
  xpot(ix, iy, iz,2) = xh(ix, iy, iz)**vpol*dexp(-psi(ix, iy, iz)*zpolB)/(1.0-fdisBas(ix, iy, iz)-fdisBNC(ix, iy, iz)-fdisBCl(ix, iy, iz))

     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz))**2 

    ! xpot(ix,iy,iz) = xpot(ix,iy,iz)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0*vpol)

     xpot(ix,iy,iz,1) = xpot(ix,iy,iz,1)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0*vpol)
     xpot(ix,iy,iz,2) = xpot(ix,iy,iz,2)*exp(Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0*vpol)



     protemp=0.0

     do ax = -2,2 
      do ay = -2,2
       do az = -2,2
            jx = ix+ax
            jy = iy+ay
            jx = mod(jx-1+5*dimx, dimx) + 1
            jy = mod(jy-1+5*dimy, dimy) + 1
            protemp=protemp + Xu(ax,ay,az)*sttemp*xtotal(jx, jy, iz+az)
       enddo
      enddo
     enddo

   !  xpot(ix, iy, iz) = xpot(ix, iy, iz)*dexp(protemp)

	    xpot(ix, iy, iz,1) = xpot(ix, iy, iz,1)*dexp(protemp)
    xpot(ix, iy, iz,2) = xpot(ix, iy, iz,2)*dexp(protemp)
   enddo
  enddo
enddo

!avpol_tosend = 0.0
!q = 0

q = 0.0
avpol_tosend = 0.0

do jj = 1, cpp ! cpp caden por prosc jj num procesad
   ii = jj+rank*cpp ! ii cant de cadenas
   q_tosend=0.0
   avpol_temp = 0.0
		

 do i=1,cuantas
   pro(i, jj)=shift

   do j=1,long
    ax = px(i, j, jj) ! cada uno para su cadena...
    ay = py(i, j, jj)
    az = pz(i, j, jj)         
    pro(i, jj) = pro(i, jj) * xpot(ax, ay, az,ct(ii))

   enddo
   do j=1,long
   avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj))= &
   avpol_temp(px(i,j, jj),py(i,j, jj),pz(i,j, jj))+pro(i, jj)*vpol*vsol/(delta**3)

   enddo

   q_tosend=q_tosend+pro(i, jj)

 enddo ! i
! norma 
 do ix=1,dimx
  do iy=1,dimy
   do iz=1,dimz
    avpol_tosend(ix,iy,iz,ct(ii))=avpol_tosend(ix, iy, iz,ct(ii)) + avpol_temp(ix,iy,iz)/q_tosend
    enddo
   enddo
 enddo
q(ii) = q_tosend ! no la envia ahora

enddo ! jj


!------------------ MPI ----------------------------------------------
!1. Todos al jefe


call MPI_Barrier(MPI_COMM_WORLD, err)

! Jefe
if (rank.eq.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, dimx*dimy*dimz*2, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)
endif
! Subordinados
if(rank.ne.0) then
! Junta avpol       
  call MPI_REDUCE(avpol_tosend, avpol, dimx*dimy*dimz*2, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err) 
!!!!!!!!!!! IMPORTANTE, LOS SUBORDINADOS TERMINAN ACA... SINO VER !MPI_allreduce!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  goto 3333
endif

!!!!!!!!!!!!!!!!!!!!!!! FIN MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------------------------------
!   Construye Ecuaciones a resolver 
!----------------------------------------------------------------------------------------------

! Qtot

do ix=1,dimx
   do iy=1,dimy
        do iz=1,dimz
  
         qtot(ix, iy, iz) =  (zpos*xpos(ix, iy, iz)+zneg*xneg(ix, iy, iz))/vsalt + avpol(ix, iy, iz,1)*zpolA/vpol*(1.0-fdisANC(ix, iy, iz)-fdisAas(ix, iy, iz)-fdisANa(ix, iy, iz)) + &
+ avpol(ix, iy, iz,2)*zpolB/vpol*(1.0-fdisBNC(ix, iy, iz)-fdisBas(ix, iy, iz)-fdisBCl(ix, iy, iz)) +    xHplus(ix, iy, iz) - xOHmin(ix, iy, iz)

        enddo
   enddo
enddo

! Volume fraction

do ix=1,dimx
   do iy=1,dimy
      do iz=1,dimz
      f(ix+dimx*(iy-1)+dimx*dimy*(iz-1))= avpol(ix,iy,iz,1)+avpol(ix,iy,iz,2) + xh(ix,iy,iz) + &
      xneg(ix, iy, iz) + xpos(ix, iy, iz) + xHplus(ix, iy, iz) + &
      xOHmin(ix, iy, iz) -1.000000d0

      enddo
   enddo
enddo

! Poisson eq.
!!G::  ESTE NO LO TOCO

do ix=1,dimx
   do iy=1,dimy
       do iz=1,dimz

       psitemp = epsfcn(ix,iy,iz)*(psi(ix+1, iy, iz)-2*psi(ix, iy, iz)+psi(ix-1, iy, iz)) 
       psitemp = psitemp + epsfcn(ix,iy,iz)*(psi(ix, iy+1, iz)-2*psi(ix, iy, iz)+psi(ix, iy-1, iz))
       psitemp = psitemp + epsfcn(ix,iy,iz)*(psi(ix, iy, iz+1) -2*psi(ix, iy, iz) + psi(ix, iy, iz-1))

       psitemp = psitemp + (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))*(epsfcn(ix+1,iy,iz)-epsfcn(ix-1,iy,iz))/4.0
       psitemp = psitemp + (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))*(epsfcn(ix,iy+1,iz)-epsfcn(ix,iy-1,iz))/4.0
       psitemp = psitemp + (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))*(epsfcn(ix,iy,iz+1)-epsfcn(ix,iy,iz-1))/4.0


      f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ntot)=(psitemp + &
      qtot(ix, iy, iz)*constq)/(-2.0)

      enddo
   enddo
enddo
 
! third block of f, ASSOCIATE RELATION
 
do ix=1,dimx
   do iy=1,dimy
       do iz=1,dimz
			 f(ix+dimx*(iy-1)+dimx*dimy*(iz-1)+ntot+ntot)=xna(ix,iy,iz)-avpol(ix,iy,iz,1) !	
      enddo
  enddo
enddo

norma = 0.0

do i = 1, 3*ntot
  norma = norma + (f(i))**2
enddo

iter = iter + 1
if(verbose.ge.3) then
if(rank.eq.0)print*,'fkfun:', iter, norma, q(1)
endif

3333 continue
ier2 = 0.0 

return
end
