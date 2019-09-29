subroutine creador

use system
use const
use system
use chainsdat
implicit none
real*8 lseg
parameter (lseg=0.35)
integer i,il,ll
integer j
real*8 indax, inday, indaz
real*8 chains(3,200,100)
real*8 altx,alty,altz,x(200),y(200),xp(200),yp(200)
real*8 theta,theta1
integer iglobal
integer nchas

il=0
iglobal=1
theta1=30.0*pi/180.0
do while (il.lt.cuantas)
  call cadenas72mr(chains,nchas)
  do i=1,nchas
    do ll=1,12
      il=il+1
      if(il.gt.cuantas) goto 100
      theta=ll*theta1
      do j=1,long
         xp(j)=chains(2,j,i)
         yp(j)=chains(3,j,i)
         x(j)=xp(j)*cos(theta)+yp(j)*sin(theta) !+ rn1
         y(j)=-xp(j)*sin(theta)+yp(j)*cos(theta) !+ rn2
         altx=x(j)
         alty=y(j)
         altz=chains(1,j,i)                 
         indax=altx
         inday=alty
         indaz=altz        
         in1(j,2)=indax
         in1(j,3)=inday
         in1(j,1)=indaz
         call pxs(il)  
      enddo
    enddo
  enddo
enddo
 100  return
end

subroutine cadenas72mr(chains,nchas)
use chainsdat
use const     
implicit none
real*8 chains(3,200,100)
integer i,state,ii,j,ive,jve
real*8 rn,state1,sitheta,cotheta,dista,lseg
real*8 siphip,cophip
character*1 test
real*8 m(3,3),mm(3,3),tt(3,3),tp(3,3),tm(3,3)
parameter (lseg=0.35)
real*8 x(3),xend(3,200),xendr(3,200)
integer nchas
real*8 rands

sitheta=sin(68.0*pi/180.0)
cotheta=cos(68.0*pi/180.0)
siphip=sin(120.0*pi/180.0)
cophip=cos(120.0*pi/180.0)
     
nchas=0

do while (nchas.eq.0) 
 x(1)=lseg
 x(2)=0.0
 x(3)=0.0
     
 xend(1,1)=lseg
 xend(2,1)=0.0
 xend(3,1)=0.0
      
 tt(1,1)=cotheta
 tt(1,2)=sitheta
 tt(1,3)=0.0
 tt(2,1)=sitheta
 tt(2,2)=-cotheta
 tt(2,3)=0.0
 tt(3,1)=0.0
 tt(3,2)=0.0
 tt(3,3)=-1.0
      
 tp(1,1)=cotheta
 tp(1,2)=sitheta
 tp(1,3)=0.0
 tp(2,1)=sitheta*cophip
 tp(2,2)=-cotheta*cophip
 tp(2,3)=siphip
 tp(3,1)=sitheta*siphip
 tp(3,2)=-cotheta*siphip
 tp(3,3)=-cophip
      
 tm(1,1)=cotheta
 tm(1,2)=sitheta
 tm(1,3)=0.0
 tm(2,1)=sitheta*cophip
 tm(2,2)=-cotheta*cophip
 tm(2,3)=-siphip
 tm(3,1)=-sitheta*siphip
 tm(3,2)=cotheta*siphip
 tm(3,3)=-cophip
      
 222  rn=rands(seed)
      
 state1=0.0
     
 m(1,1)=cotheta
 m(1,2)=sitheta
 m(1,3)=0.0
      
 m(2,1)=cos(state1)*sitheta
 m(2,2)=-cos(state1)*cotheta
 m(2,3)=sin(state1)
 m(3,1)=sin(state1)*sitheta
 m(3,2)=-sin(state1)*cotheta
 m(3,3)=-cos(state1)
      
 x(1)=m(1,1)*lseg
 x(2)=m(2,1)*lseg
 x(3)=m(3,1)*lseg
      
 xend(1,2)=lseg+x(1)
 xend(2,2)=x(2)
 xend(3,2)=x(3)
      
 do i=3,long
         rn=rands(seed)
         state=int(rn*3)
         if (state.eq.3) then 
            state=2
         endif

         if (state.eq.0) then
            call mrrrr(m,tt,mm)
          do ii=1,3
          do j=1,3
             m(ii,j)=mm(ii,j)
          enddo
          enddo
         elseif (state.eq.1) then
            call mrrrr(m,tp,mm)
          do ii=1,3
          do j=1,3
             m(ii,j)=mm(ii,j)
          enddo
          enddo
         elseif (state.eq.2) then
            call mrrrr(m,tm,mm)
          do ii=1,3
          do j=1,3
            m(ii,j)=mm(ii,j)
          enddo
          enddo
         endif
         
         x(1)=m(1,1)*lseg
         x(2)=m(2,1)*lseg
         x(3)=m(3,1)*lseg
         
         xend(1,i)=xend(1,i-1)+x(1)
         xend(2,i)=xend(2,i-1)+x(2)
         xend(3,i)=xend(3,i-1)+x(3)
 enddo       
      
 dista=0.0
 do ive=4,long
         do jve=1,ive-3
            dista=(xend(1,jve)-xend(1,ive))**(2.0)
            dista=dista+(xend(2,jve)-xend(2,ive))**(2.0)
            dista=dista+(xend(3,jve)-xend(3,ive))**(2.0)
            dista=sqrt(dista)
            if (dista.lt.lseg) then
               goto 222
            endif
         enddo
 enddo


 do i=1,300
         test='S'
         call rota36(xend,xendr,long,test)
         if (test.eq.'N')cycle
         nchas=nchas+1
         do j=1,long
            chains(1,j,nchas)=xendr(1,j)
            chains(2,j,nchas)=xendr(2,j)
            chains(3,j,nchas)=xendr(3,j)
         enddo
         if (nchas.eq.25)exit
 enddo   
enddo
return
end

subroutine rota36(xend,xendr,n,test)
      
use system
use const
implicit none      
real*8 xend(3,200),rands,xendr(3,200)
character*1 test
integer n
real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga
real*8 a,b,c
real*8 alfa, cga
integer i
real*8 gama2

fac=rands(seed)
fac1=rands(seed)
fac2=rands(seed)
alfa=fac*2*pi
cbe=2.0d0*fac1-1.0d0
gama2=fac2*2*pi

sbe=(1-cbe**2)**0.5
cal=cos(alfa)
sal=sin(alfa)
cga=cos(gama2)
sga=sin(gama2)

do i=1,n
   a=xend(1,i)
   b=xend(2,i)
   c=xend(3,i)

   xendr(1,i)=a*(-cbe*sal*sga+cal*cga)-b*(cbe*sal*cga+cal*sga)+c*sbe*sal
   xendr(2,i)=a*(cbe*cal*sga+sal*cga)+b*(cbe*cal*cga-sal*sga)-c*sbe*cal
   xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe
enddo 

do i=1,n
   if (xendr(1,i).lt.0.0) test='N'  
   if (xendr(1,i).gt.dimz*delta) then
    print*, 'rota36: Increase dimz'
    stop
   endif  
enddo

return
end
      
subroutine mrrrr(a,b,c)
implicit none
real*8 a(3,3),b(3,3),c(3,3)
integer i,j,k

do i=1,3
  do j=1,3
     c(i,j)=0
  enddo
enddo

do i=1,3
 do j=1,3
  do k=1,3
  c(i,j)=c(i,j)+a(i,k)*b(k,j)
  enddo
 enddo
enddo

return
end 
        
subroutine graftpoints
use system
use chainsdat
use const
implicit none
  
real*8 temp     
real*8 lado, postempx, postempy
real*8 lado_x, lado_y
integer ix, iy, j, ilado_x, ilado_y

seed = 86456


lado = sqrt(dfloat(dimx*dimy/(ncha/2)))
lado_x = dimx/lado  ! array
lado_y = dimy/lado 
ilado_x = int(lado_x)
ilado_y = int(lado_y)
j = 1

do ix = 1, ilado_x
do iy = 1, ilado_y

 postempx = dimx/2.0/lado_x + (ix-1)*dimx/lado_x
 postempy = dimy/2.0/lado_y + (iy-1)*dimy/lado_y

 posicion(j,1) = postempx*delta 
 posicion(j,2) = postempy*delta 

 temp = 0.0!delta!!gabi

 if(randominput.eq.1) then
	call srand(seed)
	!print *, rand()
!	posicion(j,1)=	posicion(j,1)+(rand()-0.5)*delta*2
!	posicion(j,2)=	posicion(j,2)+(rand()-0.5)*delta*2
!  temp=(rand()-0.5)*delta*2
 temp=delta*0.5
select case (j)
 case (1,5,33,37)
 posicion(j,1) = postempx*delta + temp
 posicion(j,2) = postempy*delta + temp
 case (2, 3,6,7,34,35,38,39)
 posicion(j,1) = postempx*delta + temp
 posicion(j,2) = postempy*delta 
 case (4,8,36,40)
 posicion(j,1) = postempx*delta + temp
 posicion(j,2) = postempy*delta - temp
 case (9,13,17,21,41,45,49,53)
 posicion(j,1) = postempx*delta 
 posicion(j,2) = postempy*delta + temp
 case (10, 11,14,15,18,19,22,23,42,43,46,47,50,51,54,55)
 posicion(j,1) = postempx*delta  
 posicion(j,2) = postempy*delta
 case (12,16,20,24,44,48,52,56)
 posicion(j,1) = postempx*delta 
 posicion(j,2) = postempy*delta - temp
 case (25,29,57,61)
 posicion(j,1) = postempx*delta - temp
 posicion(j,2) = postempy*delta + temp
 case (26,27,30,31,58,59,62,63)
 posicion(j,1) = postempx*delta - temp
 posicion(j,2) = postempy*delta 
 case (28,32,60,64)
 posicion(j,1) = postempx*delta - temp
 posicion(j,2) = postempy*delta - temp
 endselect

 endif

 j = j + 1

enddo
enddo

do j = ncha/2+1, ncha
posicion(j,:)=posicion(j-ncha/2,:)
enddo

end


