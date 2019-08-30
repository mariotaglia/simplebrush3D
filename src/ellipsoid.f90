
subroutine randomvect(V)
use const
implicit none
real*8, external :: rands
real*8 u, w
real*8 V(3)
real*8 theta, phi
u = rands(seed)
w = rands(seed)
theta = 2*pi*u
phi = acos(2.0*w-1.0)
V(1) = cos(theta)*sin(phi)
V(2) = sin(theta)*sin(phi)
V(3) = cos(phi)
end subroutine

subroutine make_ellipsoid
use system
use ellipsoid
implicit none
integer i

! orientation vector

orient = 0
orient(1) = 1.0

AellL(1) = Aell(1)+delta
AellL(2) = Aell(2)+delta
AellL(3) = Aell(3)+delta

AellS(1) = Aell(1)-delta
AellS(2) = Aell(2)-delta
AellS(3) = Aell(3)-delta

AAA = 0.0
AAAS = 0.0
AAAL = 0.0

do i = 1,3
AAA(i,i) = 1.0/(Aell(i)**2)
AAAS(i,i) = 1.0/(AellS(i)**2)
AAAL(i,i) = 1.0/(AellL(i)**2)
enddo

end subroutine

subroutine update_matrix(flag)
use system
use ellipsoid
use ematrix
use MPI
use const
implicit none
integer npoints ! points per cell for numerical integration 
integer counter
character*5 title
real*8 temp
logical flag

call make_ellipsoid

! rotate ellipsoid matrixes according to current rotation matrix
call rotv(AAA, rotmatrix)
call rotv(AAAS, rotmatrix)
call rotv(AAAL, rotmatrix)
call rotvo(orient, rotmatrix)

gama = 90.0/180.0*pi
npoints = 10

flag = .false.

!theta = pi/5.0
!call rotx(AAA, theta)
!call rotx(AAAS, theta)

!theta = pi/15.0
!call roty(AAA, theta)
!call roty(AAAS, theta)

call integrate(AAAL,AellL, Rell,npoints, voleps ,flag)
flag = .false. ! not a problem if eps lays outside boundaries

call integrate(AAA,Aell, Rell,npoints, volprot, flag)
call integrate(AAAS,AellS, Rell,npoints, volq, flag)

title = 'avpro'
counter = 1
call savetodisk(volprot, title, counter)

if (verbose.ge.2) then
if (rank.eq.0) then
print*, 'update_matrix: Total volumen real space= ', 4.0/3.0*pi*Aell(1)*Aell(2)*Aell(3)
print*, 'update_matrix: Total discretized volumen =', sum(volprot)*delta**3
endif
endif

temp = 4.0/3.0*pi*Aell(1)*Aell(2)*Aell(3)/( sum(volprot)*delta**3)
volprot = volprot*temp

voleps = voleps-volprot
volq = volprot-volq

title = 'aveps'
counter = 1
call savetodisk(voleps, title, counter)

title = 'avcha'
counter = 1
call savetodisk(volq, title, counter)

temp = sum(volq)

volprot = volprot * 0.99
volq = volq/temp*echarge/(delta**3) ! sum(volq) is echarge

!temp = 0
!do ix = 1,dimx
!do iy = 1,dimy
!do iz = 1,dimz
!temp = temp + volq(ix,iy,iz)*delta**3
!enddo
!enddo
!enddo
!if(rank.eq.1)print*, 'Total protein charge', temp

end subroutine




subroutine integrate(AAA,Aell, Rell, npoints,volprot, flag)

use system
implicit none
integer npoints
real*8 AAA(3,3)
real*8 volprot(dimx,dimy,dimz)
real*8 Rell(3), Aell(3)
real*8 dr(3)
integer ix,iy,iz,ax,ay,az
real*8 vect
logical flagin, flagout
real*8 intcell
real*8 mmmult
integer jx,jy
integer Rpos(3)
integer maxAell
logical flag

volprot = 0.0

maxAell = int(max(Aell(1)/delta,Aell(2)/delta,Aell(3)/delta))+2
Rpos(1) = int(Rell(1)/delta)
Rpos(2) = int(Rell(2)/delta)
Rpos(3) = int(Rell(3)/delta)


! Make a list of the cells that have no ellipsoid, those that have part ellipsoid and those that have full ellipsoid
! Consider boundary conditions 

do ix = Rpos(1)-maxAell, Rpos(1)+maxAell
do iy = Rpos(2)-maxAell, Rpos(2)+maxAell
do iz = Rpos(3)-maxAell, Rpos(3)+maxAell

jx=mod(ix+dimx-1,dimx)+1
jy=mod(iy+dimy-1,dimy)+1

flagin=.false.
flagout=.false.

do ax = -1,0
do ay = -1,0
do az = -1,0

dr(1) = (ix+ax)*delta - Rell(1)
dr(2) = (iy+ay)*delta - Rell(2)
dr(3) = (iz+az)*delta - Rell(3)

vect = mmmult(dr,AAA)
if(vect.le.1.0) then           ! inside the ellipsoid
  flagin=.true.
  if(flagout.eqv..true.) then 
      if((iz.ge.1).and.(iz.le.dimz)) then
         volprot(jx,jy,iz) = intcell(AAA, Rell, ix,iy,iz, npoints)
      else
         print*,'update_matrix: Flag', iz
         flag=.true.
      endif   
      goto 999 ! one in and one out, break the cycle
  endif
else 
  flagout=.true.
  if(flagin.eqv..true.) then
      if((iz.ge.1).and.(iz.le.dimz)) then
          volprot(jx,jy,iz) = intcell(AAA, Rell, ix,iy,iz, npoints)
      else
         print*,'update_matrix: Flag', iz
          flag=.true.
      endif
      goto 999 ! one in and one out, break the cycle
  endif
endif

enddo
enddo
enddo

if((flagin.eqv..true.).and.(flagout.eqv..false.)) then 
      if((iz.ge.1).and.(iz.le.dimz)) then
         volprot(jx,jy,iz)=1.0 ! all inside
      else
         print*,'update_matrix: Flag', iz
         flag=.true.
      endif
endif
999 continue

enddo
enddo
enddo
end subroutine

double precision function intcell(AAA,Rell,ix,iy,iz,n)
use system
implicit none
real*8 AAA(3,3)
real*8 Rell(3)
integer ix,iy,iz,ax,ay,az
integer cc
real*8 vect
integer n
real*8 mmmult
real*8 dr(3)

cc = 0
do ax = 1, n
do ay = 1, n
do az = 1, n

dr(1) = ix*delta-(ax)*delta/float(n) - Rell(1)
dr(2) = iy*delta-(ay)*delta/float(n) - Rell(2)
dr(3) = iz*delta-(az)*delta/float(n) - Rell(3)

!dr(1) = mod(dr(1)+(dimx*delta),(dimx*delta))
!dr(2) = mod(dr(2)+(dimy*delta),(dimy*delta))

vect = mmmult(dr,AAA)

if(vect.le.1.0)cc=cc+1

enddo
enddo
enddo

intcell = float(cc)/(float(n)**3)
end function

double precision function mmmult(V,A)
implicit none
real*8 V(3)
real*8 A(3,3)
real*8 C(3)
C(1) = A(1,1)*V(1)+A(1,2)*V(2)+A(1,3)*V(3)
C(2) = A(2,1)*V(1)+A(2,2)*V(2)+A(2,3)*V(3)
C(3) = A(3,1)*V(1)+A(3,2)*V(2)+A(3,3)*V(3)
mmmult = V(1)*C(1) + V(2)*C(2) + V(3)*C(3)
endfunction

subroutine rotv(A, B) ! applies rotation matrix B to ellipsoid matrix A
implicit none
real*8 A(3,3)
real*8 B(3,3)
real*8 BT(3,3)
BT = TRANSPOSE(B)
A = MATMUL(A, B)
A = MATMUL(BT, A)
end subroutine

subroutine rotvo(orient, B) ! applies rotation matrix B to vector orient
implicit none
real*8 B(3,3)
real*8 orient(3)
orient = MATMUL(B, orient)
end subroutine

subroutine rotvm(A, theta, V) ! rotates the rotation matrix A theta degress round V
implicit none
real*8 A(3,3)
real*8 theta
real*8 B(3,3)
real*8 V(3)
B(1,1)=cos(theta)+V(1)*V(1)*(1.0-cos(theta))
B(1,2)=V(1)*V(2)*(1.0-cos(theta))-V(3)*sin(theta)
B(1,3)=V(1)*V(3)*(1.0-cos(theta))+V(2)*sin(theta)
B(2,1)=V(2)*V(1)*(1.0-cos(theta))+V(3)*sin(theta)
B(2,2)=cos(theta)+V(2)*V(2)*(1.0-cos(theta))
B(2,3)=V(2)*V(3)*(1.0-cos(theta))-V(1)*sin(theta)
B(3,1)=V(3)*V(1)*(1.0-cos(theta))-V(2)*sin(theta)
B(3,2)=V(3)*V(2)*(1.0-cos(theta))+V(1)*sin(theta)
B(3,3)=cos(theta)+V(3)*V(3)*(1.0-cos(theta))
A = MATMUL(B, A)
end subroutine














