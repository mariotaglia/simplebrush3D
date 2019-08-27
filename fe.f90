!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!    Free Energy Calculation...
!
!
!
subroutine Free_Energy_Calc(looped)

use system
use const
use fields_fkfun
use MPI
use molecules
use kai
use bulk
use results
implicit none

integer looped
real*8  q_tosend(ncha)
real*8  q0(ncha)
real*8 F_Mix_s, F_Mix_pos
real*8 F_Mix_neg, F_Mix_Hplus
real*8 Free_energy2, sumpi, sumrho, sumel, sumdiel, sum, mupol,sumas
real*8 F_Mix_OHmin, F_Conf, F_Eq, F_vdW, F_eps, F_electro
real*8 pro0(cuantas, cpp)
real*8 Free_Energy
 
! MPI
integer stat(MPI_STATUS_SIZE) 
integer source
integer dest
integer tag
parameter(tag = 0)
integer err

! Dummies
integer ix, iy, iz, i, ii, ax, ay, az, jj
integer jx, jy, iii

real*8 gradpsi2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!
!  Recupera pro(i) de todos los procesos para calculo de F
!

! Subordinados

q0 = 0
q_tosend = 0

if(rank.ne.0) then
       dest = 0
! Envia q

       do jj = 1, cpp
       iii = rank*cpp+jj
       q_tosend(iii) = q(iii) 
       enddo

        call MPI_REDUCE(q_tosend, q0, ncha, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

! Envia pro
        CALL MPI_SEND(pro, cuantas*cpp , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD,err)
      goto 888
endif


      Free_Energy = 0.0
      Free_Energy2 = 0.0

! 1. Mezcla solvente

      F_Mix_s = 0.0 

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz
      F_Mix_s = F_Mix_s + xh(ix, iy,iz)*(dlog(xh(ix, iy, iz))-1.0)
      F_Mix_s = F_Mix_s - xsolbulk*(dlog(xsolbulk)-1.0)
      enddo      
      enddo      
      enddo      
      F_Mix_s = F_Mix_s * delta**3/vsol
      Free_Energy = Free_Energy + F_Mix_s

! 2. Mezcla ion positivo

      F_Mix_pos = 0.0 

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz
      
      F_Mix_pos = F_Mix_pos + xpos(ix, iy,iz) &
      *(dlog(xpos(ix, iy, iz)/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))

      F_Mix_pos = F_Mix_pos - xposbulk &
      *(dlog(xposbulk/vsalt)-1.0-dlog(expmupos) + dlog(vsalt))

      enddo
      enddo
      enddo
      F_Mix_pos = F_Mix_pos * delta**3/vsol/vsalt
      Free_Energy = Free_Energy + F_Mix_pos

! 3. Mezcla ion negativo

      F_Mix_neg = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz


      F_Mix_neg = F_Mix_neg + xneg(ix, iy,iz) &
      *(dlog(xneg(ix, iy, iz)/vsalt)-1.0- dlog(expmuneg) + dlog(vsalt))

      F_Mix_neg = F_Mix_neg - xnegbulk &
      *(dlog(xnegbulk/vsalt)-1.0- dlog(expmuneg) + dlog(vsalt))

      enddo 
      enddo 
      enddo 
      F_Mix_neg = F_Mix_neg * delta**3/vsol/vsalt
      Free_Energy = Free_Energy + F_Mix_neg

! 4. Mezcla protones

      F_Mix_Hplus = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz


      F_Mix_Hplus = F_Mix_Hplus &
     +xHplus(ix, iy, iz)*(dlog(xHplus(ix,iy,iz))-1.0 -dlog(expmuHplus))

      F_Mix_Hplus = F_Mix_Hplus &
     -xHplusbulk*(dlog(xHplusbulk)-1.0 -dlog(expmuHplus))

      enddo
      enddo
      enddo
      F_Mix_Hplus = F_Mix_Hplus * delta**3/vsol
      Free_Energy = Free_Energy + F_Mix_Hplus

! 5. Mezcla hidroxilos

      F_Mix_OHmin = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz


      F_Mix_OHmin = F_Mix_OHmin + xOHmin(ix, iy,iz)*(dlog(xOHmin(ix, iy, iz))-1.0-dlog(expmuOHmin))

      F_Mix_OHmin = F_Mix_OHmin - xOHminbulk*(dlog(xOHminbulk)-1.0-dlog(expmuOHmin))

      enddo
      enddo
      enddo
      F_Mix_OHmin = F_Mix_OHmin * delta**3/vsol
      Free_Energy = Free_Energy + F_Mix_OHmin

! 6. Entropia interna polimero

      F_Conf = 0.0

! Jefe

       if (rank.eq.0) then ! Igual tiene que serlo, ver arriba

       do jj = 1, cpp
       iii = rank*cpp+jj
       q_tosend(iii) = q(iii) 
       enddo

        call MPI_REDUCE(q_tosend, q0, ncha, &
        MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

       do jj = 1, cpp 
       do i = 1, cuantas
       iii = jj
       
         F_Conf = F_Conf + (pro(i, jj)/q0(iii)) &
      *dlog((pro(i, jj))/q0(iii))

       enddo
       enddo 

         do ii = 2, size ! loop sobre los procesadores restantes

        source = ii-1

        call MPI_RECV(pro0, cuantas*ncha, &
        MPI_DOUBLE_PRECISION, source, tag, MPI_COMM_WORLD,stat, err)


       do jj = 1, cpp
       do i = 1, cuantas

       iii = (ii-1)*cpp+jj

         F_Conf = F_Conf + (pro0(i, jj)/q0(iii))*dlog((pro0(i, jj))/q0(iii))

       enddo
       enddo

       enddo ! ii

       endif ! rank

      Free_Energy = Free_Energy + F_Conf



! 7. Chemical Equilibrium
!comento original 
     ! F_Eq = 0.0 
            

     ! do ix  = 1, dimx
     ! do iy  = 1, dimy
     ! do iz  = 1, dimz
      
     ! F_Eq = F_Eq + fdis(ix, iy, iz)*dlog(fdis(ix, iy, iz)) &
     ! *avpol(ix,iy,iz)/vpol

  !    F_Eq = F_Eq + (1.0-fdis(ix,iy,iz)) &
   !   *dlog(1.0-fdis(ix, iy,iz))*avpol(ix, iy,iz)/vpol

!      F_Eq = F_Eq + (1.0-fdis(ix, iy,iz))*dlog(K0)*avpol(ix, iy,iz)/vpol

!      F_Eq = F_Eq + (1.0-fdis(ix, iy,iz)) &!
!      *(-dlog(expmuHplus))*avpol(ix, iy, iz)/vpol

 !     enddo
 !     enddo
 !     enddo

!      F_eq = F_eq *delta**3/vsol

 !     Free_Energy = Free_Energy + F_Eq


      F_Eq = 0.0 
            

      do ix  = 1, dimx
      do iy  = 1, dimy
      do iz  = 1, dimz
!!!APORTE POL-A
 		if (1.0d-10 < (1.0-fdisANC(ix,iy,iz)-fdisAas(ix,iy,iz)-fdisANa(ix,iy,iz))) then  !! Es posible que este if no haga falta 
 			F_Eq = F_Eq + (1.0-fdisANC(ix,iy,iz)-fdisAas(ix,iy,iz)-fdisANa(ix,iy,iz))&
*dlog(1.0-fdisANC(ix,iy,iz)-fdisAas(ix,iy,iz)-fdisANa(ix,iy,iz))*avpol(ix,iy,iz,1)/vpol 
		 endif
 		F_Eq = F_Eq + (fdisANC(ix,iy,iz))*dlog(fdisANC(ix,iy,iz))*avpol(ix,iy,iz,1)/vpol
 		F_Eq = F_Eq + (fdisANa(ix,iy,iz))*dlog(fdisANa(ix,iy,iz))*avpol(ix,iy,iz,1)/vpol
     F_Eq = F_Eq + (fdisANC(ix,iy,iz))*dlog(K0A)*avpol(ix,iy,iz,1)/vpol
 		F_Eq = F_Eq + (fdisANC(ix,iy,iz))*(-dlog(expmuHplus))*avpol(ix,iy,iz,1)/vpol
    F_Eq = F_Eq + (fdisANa(ix,iy,iz))*(dlog(K0ANa))*avpol(ix,iy,iz,1)/vpol
     F_Eq = F_Eq + (fdisANa(ix,iy,iz))*(-dlog(expmupos))*avpol(ix,iy,iz,1)/vpol

 		F_Eq = F_Eq + (fdisAas(ix,iy,iz))*(-dlog(K0Eo))*avpol(ix,iy,iz,1)/vpol

 	if (1.0d-10 < fdisAas(ix,iy,iz)) then 
 	 F_Eq = F_Eq + (fdisAas(ix,iy,iz))*dlog(fdisAas(ix,iy,iz))*avpol(ix,iy,iz,1)/vpol
 		if (1.0d-10 < avpol(ix,iy,iz,1))then 
  			F_Eq = F_Eq + (-fdisAas(ix,iy,iz))*(dlog(avpol(ix,iy,iz,1)&
*fdisAas(ix,iy,iz))-1.0)*avpol(ix,iy,iz,1)/vpol ! usando que Vpol =Vab
 		endif
 	endif

	   enddo
      enddo
      enddo
!! APORTE POL-B
      do ix  = 1, dimx
      do iy  = 1, dimy
      do iz  = 1, dimz

 		if (1.0d-10 < (1.0-fdisBNC(ix,iy,iz)-fdisBas(ix,iy,iz)-fdisBCl(ix,iy,iz))) then
     	F_Eq = F_Eq + (1.0-fdisBNC(ix,iy,iz)-fdisBas(ix,iy,iz)-fdisBCl(ix,iy,iz))&
*dlog(1.0-fdisBNC(ix,iy,iz)-fdisBas(ix,iy,iz)-fdisBCl(ix,iy,iz))*avpol(ix,iy,iz,2)/vpol
 		endif

 	  F_Eq = F_Eq + (fdisBNC(ix,iy,iz))*dlog(fdisBNC(ix,iy,iz))*avpol(ix,iy,iz,2)/vpol
    F_Eq = F_Eq + (fdisBCl(ix,iy,iz))*dlog(fdisBCl(ix,iy,iz))*avpol(ix,iy,iz,2)/vpol

    F_Eq = F_Eq + (fdisBNC(ix,iy,iz))*dlog(K0B)*avpol(ix,iy,iz,2)/vpol
   F_Eq = F_Eq + (fdisBNC(ix,iy,iz))*(-dlog(expmuOHmin))*avpol(ix,iy,iz,2)/vpol
    F_Eq = F_Eq + (fdisBCl(ix,iy,iz))*(dlog(K0BCl))*avpol(ix,iy,iz,2)/vpol
    F_Eq = F_Eq + (fdisBCl(ix,iy,iz))*(-dlog(expmuneg))*avpol(ix,iy,iz,2)/vpol

			if ((1.0d-10 < fdisBas(ix,iy,iz))) then 
			F_Eq = F_Eq +( (fdisBas(ix,iy,iz))*dlog(fdisBas(ix,iy,iz)) )*avpol(ix,iy,iz,2)/vpol
		 	endif
	   enddo
      enddo
      enddo
     F_eq = F_eq *delta**3/vsol

      Free_Energy = Free_Energy + F_Eq

!!G !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!G :F_Eq

! 8.vdW ! Ojo, los kai son negativos => atraccion

       F_vdW = 0.0

      do ix = 1, dimx
      do iy = 1, dimy
      do iz = 1, dimz

            do ax = -2,2
            do ay = -2,2
            do az = -2,2


            jx = ix+ax
            jy = iy+ay

            jx = mod(jx-1+5*dimx, dimx) + 1
            jy = mod(jy-1+5*dimy, dimy) + 1
  
      F_vdW = F_vdW - 0.5000*delta**3*xtotal(ix, iy, iz) &
        *xtotal(jx, jy,iz+az)*Xu(ax, ay, az)*st

            enddo
            enddo
            enddo


      enddo
      enddo
      enddo

      Free_Energy = Free_Energy + F_vdW

! 9. Electrostatic ! OJO

      F_electro = 0.0    

      do ix  = 1, dimx
      do iy  = 1, dimy
      do iz  = 1, dimz

!    gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))**2

!    gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))*(psi(ix+1,iy,iz)-psi(ix-1,iy,iz))
!    gradpsi2 = gradpsi2 + (psi(ix,iy+1,iz)-psi(ix,iy-1,iz))*(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))
!    gradpsi2 = gradpsi2 + (psi(ix,iy,iz+1)-psi(ix,iy,iz-1))*(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))
!    gradpsi2 = gradpsi2/4.0

!    gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))*(psi(ix+1,iy,iz)-psi(ix,iy,iz))
!    gradpsi2 = gradpsi2 + (psi(ix,iy+1,iz)-psi(ix,iy,iz))*(psi(ix,iy+1,iz)-psi(ix,iy,iz))
!    gradpsi2 = gradpsi2 + (psi(ix,iy,iz+1)-psi(ix,iy,iz))*(psi(ix,iy,iz+1)-psi(ix,iy,iz))

      F_electro = F_electro &
       + delta**3*psi(ix, iy, iz)*qtot(ix, iy, iz)/2.0/vsol

!       gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz))**2

!      F_electro = F_electro &
!       + (delta**3)/vsol*(psi(ix, iy, iz)*qtot(ix, iy, iz) - 0.5/constq*gradpsi2*epsfcn(ix,iy,iz))

      enddo
!      F_electro = F_electro + sigmaq*psi(ix, iy, 0)/2.0 ! OJO!!! REVISAR!!!!
!      F_electro = F_electro + sigmaq*psi(ix, iy, dimz+1)/2.0 ! OJO!!! REVISAR!!!!
      enddo
      enddo


      Free_Energy = Free_Energy + F_electro

      if (verbose.ge.1) then
      print*, 'Free_Energy_Calc: Free energy(1) = ', Free_energy
      endif


      F_eps = 0.0 

! minimal F

      Free_Energy2 = 0.0

        sumpi = 0.0
        sumrho=0.0
        sumel=0.0
        sumdiel = 0.0

!!G:  Término asociacion !!!!!!!!!!!G:sumas
	sumas=0.0

        do ix=1,dimx
        do iy=1,dimy
        do iz=1,dimz

           sumpi = sumpi+dlog(xh(ix, iy, iz))     
           sumpi = sumpi-dlog(xsolbulk)
     
           sumrho = sumrho + ( - xh(ix, iy, iz) -xHplus(ix, iy, iz) &
        - xOHmin(ix, iy, iz) - (xpos(ix, iy, iz)+xneg(ix, iy, iz))/vsalt)! sum over  rho_i i=+,-,s


           sumrho = sumrho - ( - xsolbulk -xHplusbulk &
       -xOHminbulk - (xposbulk+xnegbulk)/vsalt) ! sum over  rho_i i=+,-,s

           sumel = sumel - qtot(ix, iy, iz)*psi(ix, iy, iz)/2.0 
      
          gradpsi2 = (psi(ix+1,iy,iz)-psi(ix,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz))**2
          sumdiel = sumdiel + 0.5/constq*xtotal(ix,iy,iz)*gradpsi2*Depsfcn(ix,iy,iz)

		!!!G:TÉRMINO  DE PAIR AB !!!!!!!G:SUMAS.
			sumas=sumas +avpol(ix,iy,iz,1)*fdisAas(ix,iy,iz)/vpol

         enddo
         enddo
         enddo
         
!     gradpsi2 = (psi(ix+1,iy,iz)-psi(ix-1,iy,iz))**2+(psi(ix,iy+1,iz)-psi(ix,iy-1,iz))**2+(psi(ix,iy,iz+1)-psi(ix,iy,iz-1))**2
!     xpot(ix, iy, iz) = xpot(ix,iy,iz)*exp(-Depsfcn(ix,iy,iz)*(gradpsi2)/constq/2.0*vpol/4.0)


         sumpi = (delta**3/vsol)*sumpi
         sumrho = (delta**3/vsol)*sumrho
         sumel = (delta**3/vsol)*sumel
         sumdiel = (delta**3/vsol)*sumdiel

!!G: 
   		sumas= (delta**3/vsol)*sumas		
         sum = sumpi + sumrho + sumel + sumdiel+sumas

         sum = sumpi + sumrho + sumel + sumdiel


         do ii = 1, ncha
         Free_Energy2 = Free_Energy2-dlog(q0(ii)/shift) 
         enddo

!         do ii = 1, ncha
!         Free_Energy2 = Free_Energy2-dlog(qA0(ii)/shiftA) 
!         enddo
!         do ii = 1, ncha
!         Free_Energy2 = Free_Energy2-dlog(qB0(ii)/shiftB) 
!         enddo

         Free_Energy2 = Free_Energy2 + sum - F_vdW

      if (verbose.ge.1) then
      print*, 'Free_Energy_Calc: Free energy(2) = ', Free_energy2
      endif

!!diffener= Free_energy- Free_energy2 !!!!!!!g:

!!print*, 'Free energy:', Free_energy, Free_energy2,diffener !!!!G

! Guarda energia libre


        mupol = 0.0
        do ii = 1, ncha
        mupol = mupol - dlog(q0(ii)/shift)
        enddo
        mupol = mupol/ncha
!

!!!!!!!!!!!!gg!!!!!!!!!!
!!!!G: mupola  y mu pol B

!        mupolA = 0.0
 !       do ii = 1, ncha
   !     mupolA = mupolA - dlog(qA0(ii)/shiftA)
!        enddo
!        mupolA = mupolA/ncha        
!		
!       mupolB = 0.0
       ! do ii = 1, ncha
!        mupolB = mupolB - dlog(qB0(ii)/shiftB)
!        enddo
!        mupolB = mupolB/ncha
!
!!G: h sin el factor dos

! EL h  está sin el factor 2 multiplicando 
	
!h capa 
!			normhcapa=0.0

!      do ix  = 1, dimx!
!      do iy  = 1, dimy
!      do iz  = 1, dimz

!
!				normhcapa=normhcapa+(avpolA(ix,iy,iz)+avpolB(ix,iy,iz))*delta**3/vpol
!				hcapa=hcapa+(avpolA(ix,iy,iz)+avpolB(ix,iy,iz))*(ix+iy+iz)*delta*delta**3/vpol

!			enddo
!			enddo
!			enddo
!			if(1.0d-10 <normhcapa)then
!				hcapa=2.0*hcapa/normhcapa
!			!endif

!!!!!!!!!!!!!gg!!!!!!!!!!
          if(rank.eq.0) then

 
         write(301,*)looped, Free_energy
          flush(301)
         write(302,*)looped, F_Mix_s 
         write(303,*)looped, F_Mix_pos
         write(304,*)looped, F_Mix_neg
         write(305,*)looped, F_Mix_Hplus
         write(306,*)looped, F_Mix_OHmin
         write(307,*)looped, F_Conf
         write(308,*)looped, F_Eq
         write(309,*)looped, F_vdW
         write(410,*)looped, F_eps
         write(311,*)looped, F_electro

         write(312,*)looped, Free_energy2

         write(313,*)looped, mupol
!!       write(313,*)looped, mupolA
!!       write(313,*)looped, mupolB
!!			write(316,*)cc, ccc, hcapa

         endif
 
 888     call MPI_BCAST(free_energy, 1, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD, err)

         return

         end




