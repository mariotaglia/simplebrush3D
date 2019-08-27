use system
use MPI
use kinsol
use const
implicit none
integer counter, counterr
real*8 temp
real*8, external :: rands
logical flag

counter = 1
counterr = 1

call readinput
call initmpi
if(rank.eq.0)print*, 'MPI OK'

call initconst
call initall
call ctinit
call allocation   !!  G:  ver 

verbose = 5

! Calculate poor-solvent coefficients
call kais
if(rank.eq.0)print*, 'Kai OK'

call  graftpoints
if(rank.eq.0)print*, 'Graftpoints OK'

call creador ! Genera cadenas
if(rank.eq.0)print*, 'Creador OK'

if(infile.eq.0) then
   call solve
   call Free_Energy_Calc(counter)
   call savedata(counter)
else
   call retrivefromdisk(counter)
   counterr = counter
   if(rank.eq.0)print*, 'Load input from file'
   infile = 2
   call solve
   call Free_Energy_Calc(counter)
endif

call endall
end


subroutine initmpi
use MPI
use chainsdat
implicit none

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

if(mod(ncha, size).ne.0) then
  print*, 'Cannot divide', size, 'processors among ',ncha, 'chains'
  call MPI_FINALIZE(ierr) ! finaliza MPI
  stop
endif

cpp = ncha/size
call allocatecpp

end subroutine

subroutine initconst
use const
use molecules
implicit none
pi = acos(-1.0)
seed = 15615
gama = 90.0/180.0 * pi
lb = 0.714 ! bjerrum lenght in nm
zpos = 1.0
zneg = -1.0
vsol = 0.030
vsalt=((4.0/3.0)*pi*(0.2)**3)/vsol  ! volume salt in units of vsol 0.2=radius salt  
vpol= 0.11/vsol ! ((4.0/3.0)*pi*(0.2)**3)/vsol  ! volume polymer segment in units of vsol 
constq=delta*delta*4.0*pi*lb/vsol   ! multiplicative factor in poisson eq  
pKw = 14
error = 1e-4 ! para comparar con la norma...
errel=1d-6
itmax=200
end subroutine



subroutine initall
use molecules
use const
use bulk
use MPI
use chainsdat
use inputtemp
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Open common files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

if(rank.eq.0) then
       open(unit=301, file='F_tot.dat')
       open(unit=302, file='F_mixs.dat')
       open(unit=303, file='F_mixpos.dat')
       open(unit=304, file='F_mixneg.dat')
       open(unit=305, file='F_mixH.dat')
       open(unit=306, file='F_mixOH.dat')
       open(unit=307, file='F_conf.dat')
       open(unit=308, file='F_eq.dat')
       open(unit=309, file='F_vdW.dat')
       open(unit=410, file='F_eps.dat')
       open(unit=311, file='F_electro.dat')
       open(unit=312, file='F_tot2.dat')
       open(unit=314, file='F_mixpos2.dat')
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input-dependent variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!

constqE = vpol/(2.0d0*constq)
dielW = 78.54
dielPr = dielP/dielW
dielSr = dielS/dielW

zpolA = -1.0      ! charge of polyelectrolyte segment A			 G::
zpolB = 1.0      ! charge of polyelectrolyte segment B		 	G ::

pKw = 14.0				!G::
kW=10**(-pKw)		!G::
KaA=10**(-pKaA)
KaB=10**(-pKaB)
KEo=10**(-pKEo) !!!!!!!!!!!!!!!!!!!!!!!
KaANa=10**(-pKaANa)
KaBCl=10**(-pKaBCl)



cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
pOHbulk= pKw -pHbulk
cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  
xsalt=(csalt*Na/(1.0d24))*(vsalt*vsol)   ! volume fraction salt,csalt in mol/l 

if(pHbulk.le.7) then  ! pH<= 7
  xposbulk=xsalt/zpos
  xnegbulk=-xsalt/zneg+(xHplusbulk -xOHminbulk) *vsalt ! NaCl+ HCl  
else                  ! pH >7 
  xposbulk=xsalt/zpos +(xOHminbulk -xHplusbulk) *vsalt ! NaCl+ NaOH   
  xnegbulk= -xsalt/zneg 
endif

xsolbulk=1.0 -xHplusbulk -xOHminbulk -xnegbulk -xposbulk 

!K0 = (Ka*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 

K0A = (KaA*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
K0ANa = (KaANa*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
K0BCl = (KaBCl*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
K0B = (Kw/KaB*vsol/xsolbulk)*(Na/1.0d24) 
K0Eo = (KEo)*(1.0d24/Na) !!!!!!!!!!!!!!!!!!!!!!!

expmupos = xposbulk /xsolbulk**vsalt
expmuneg = xnegbulk /xsolbulk**vsalt
expmuHplus = xHplusbulk /xsolbulk   ! vsol = vHplus 
expmuOHmin = xOHminbulk /xsolbulk   ! vsol = vOHmin 

end subroutine


subroutine readinput
use molecules
use const
use bulk
use MPI
use chainsdat
use inputtemp
implicit none
character basura

!!!!!!!!!!!!!!!!!!!!!!!!!
! Read input variables
!!!!!!!!!!!!!!!!!!!!!!!!!

read(8,*), basura
read(8,*), dimx,dimy,dimz

read(8,*), basura
read(8,*), ncha

read(8,*), basura
read(8,*), long

read(8,*), basura
read(8,*), cuantas

read(8,*), basura
read(8, *), dielP, dielS


read(8, *), basura
read(8, *), csalt


read(8, *), basura
read(8, *), pKaA    ! pKaA of weak polyacid segments

read(8, *), basura
read(8, *), pKaB    ! pKaB of weak polyacid segments

read(8, *), basura
read(8, *), pKaANa    ! pKaA of weak polyacid segments

read(8, *), basura
read(8, *), pKaBCl    ! pKaB of weak polyacid segments

read(8, *), basura					!!!!!!!!!!!!!!!!!!!!!!!
read(8, *), pKEo     ! Interation 	!!!!!!!!!!!!!!!!!!!!!!!

read(8, *), basura
read(8, *), pHbulk

read(8, *), basura
read(8, *), infile

read(8, *), basura
read(8, *), st

read(8, *), basura
read(8, *), randominput
end subroutine


subroutine endall
use MPI
implicit none

!!!!!!!!!!!!!!!!!!!!!!
! Close common files
!!!!!!!!!!!!!!!!!!!!!!

close(301)
close(302)
close(303)
close(304)
close(305)
close(306)
close(307)
close(308)
close(309)
close(310)
close(311)
close(312)
close(313)

call MPI_FINALIZE(ierr) ! finaliza MPI    
stop

end subroutine



subroutine savedata(cccc)
use system
use results
use const
use molecules
use chainsdat
use kai
use fields_fkfun
use MPI
use kinsol
implicit none
integer cccc
character*20 filename
character*5  title
real*8 temp(dimx,dimy,dimz,2)
real*8 temppsi(dimx,dimy,dimz) !!chequear G!!!!!!!!!!!!!!!!!!!1

!----------------------------------------------------------
!  OUTPUT
!----------------------------------------------------------

if(rank.eq.0) then ! solo el jefe escribe a disco....
  ! Guarda infile
!  write(filename,'(A4, I3.3, A4)')'out.', cccc, '.dat'
!  open(unit=45, file=filename)
!   do i = 1, 2*n
!    write(45, *)x1(i)
!   enddo
!  close(45)

!!!!!!!!!!!!!!!!!!! Guarda archivos !!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Polimero
!!!!!!!!!!! es mejor  mantener dimensiones

  temp(:,:,:,:) = avpol(:,:,:,:)

  title = 'avpol'
  call savetodisk(temp, title, cccc)

! Solvente
!  title = 'avsol'
!  call savetodisk(xh, title, cccc)
! Cationes
!  title = 'avpos'
!  call savetodisk(xpos, title, cccc)
! Aniones
!  title = 'avneg'
!  call savetodisk(xneg, title, cccc)
! H+
!  title = 'avHpl'
!  call savetodisk(xHplus, title, cccc)
! OH-
!  title = 'avOHm'
!  call savetodisk(xOHmin, title, cccc)
! fdis
!  title = 'frdis'
!  call savetodisk(fdis, title, cccc)
! Potencial electrostatico

  temppsi(1:dimx,1:dimy, 1:dimz) = psi(1:dimx,1:dimy, 1:dimz)

  title = 'poten'
  call savetodisk(temp, title, cccc)

! system

  write(filename,'(A7, I3.3, A4)')'system.', cccc, '.dat'
  open (unit=310, file=filename)
  write(310,*)'st          = ',st ! residual size of iteration vector
  write(310,*)'fnorm       = ',norma ! residual size of iteration vector
  write(310,*)'length seg  = ',0.35 ! value see subroutine cadenas
  write(310,*)'delta       = ',delta
  write(310,*)'vsol        = ',vsol
  write(310,*)'vsalt       = ',vsalt*vsol
  write(310,*)'vpol       = ',vpol*vsol
  write(310,*)'pKw         = ',pKw
  write(310,*)'zpos        = ',zpos
  write(310,*)'zneg        = ',zneg
  write(310,*)'long        = ',long
  write(310,*)'iterations  = ',iter
  write(310,*)'sigma cad/nm2 = ',ncha/(dimx*dimy*delta*delta)
  write(310,*)'gama =          ', gama, gama*180/pi
  write(310,*)'kai =          ', Xu
  close(310)

endif

end subroutine

subroutine store2disk(counter) ! saves state to disk
use kinsol
use results
use MPI
use const
implicit none
integer counter

if(rank.eq.0) then
open (unit=8, file='out.out', form='unformatted')
write(8)counter
write(8)seed
write(8)xflag
close(8)
endif
end subroutine


subroutine retrivefromdisk(counter) ! saves state to disk
use kinsol
use results
use const
implicit none
integer counter

open (unit=8, file='in.in', form='unformatted')
read(8)counter
read(8)seed
read(8)xflag
close(8)
end subroutine


