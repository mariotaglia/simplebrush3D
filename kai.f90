
!#####################################################################
!
! Este programa calcula los kai para poor-solvent en 3D (geometria 
! esferica) usando un metodo de MC
!
!#####################################################################


subroutine kais
use system
use const
use kai
use molecules
use MPI
implicit none

real*8 lseg ! largo del segmento
real*8 l ! medio lseg, radio del segmento
 
integer MCsteps ! numero de steps de MC
integer ix, iy , iz
real*8 x,y,z, radio
integer limit
parameter (limit = 3)
real*8 matriz(-limit:limit, -limit:limit, -limit:limit) ! matriz de kai
integer i
real*8 rands
real*8 suma

if(rank.eq.0)print*,'kais: Kai calculation'
suma = 0.0
do ix = -limit, limit
 do iy = -limit, limit
  do iz = -limit, limit
  matriz(ix, iy, iz) = 0.0
  enddo
 enddo
enddo

MCsteps = 10
if(rank.eq.0)print*, 'kais: CORREGIR MCSTEPS'

lseg=0.35
l = lseg 

do i = 1, MCsteps
 x = 3.0*(rands(seed)-0.5)*delta ! numero al azar entre -1.5*delta y 1.5*delta
 y = 3.0*(rands(seed)-0.5)*delta ! numero al azar entre -1.5*delta y 1.5*delta
 z = 3.0*(rands(seed)-0.5)*delta ! numero al azar entre -1.5*delta y 1.5*delta

 radio = sqrt(x**2 + y**2 + z**2) ! espacio real

 if(radio.gt.(1.5*delta)) cycle ! No esta dentro de la esfera del cut-off   
 if(radio.lt.l) cycle ! esta dentro de la esfera del segmento

 ! celda 
 ix = int(anint(x/delta))   ! espacio de la grilla
 iy = int(anint(y/delta))
 iz = int(anint(z/delta))
 matriz(ix, iy, iz) = matriz(ix, iy, iz) + (l/radio)**6
enddo



do ix = -2, 2
do iy = -2, 2
do iz = -2, 2
 Xu(ix, iy, iz) = matriz(ix, iy, iz)/MCsteps*((3.0*delta)**3)
 suma = suma +  matriz(ix, iy, iz)/MCsteps*((3.0*delta)**3)
enddo
enddo
enddo

if(rank.eq.0)print*, 'kais: Sum 5x5', suma

suma = 0.0
do ix = -limit, limit
do iy = -limit, limit
do iz = -limit, limit
 suma = suma +  matriz(ix, iy, iz)/MCsteps*((3.0*delta)**3)
enddo
enddo
enddo

if(rank.eq.0)print*, 'kais: Total Sum', suma
 
end



