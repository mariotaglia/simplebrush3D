!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Esta subrutina se encarga de poner a todas los segmentos dentro del slab

subroutine pxs(il)

use system
use MPI
use chainsdat
use conformations
use const
implicit none
    
integer il
integer j, ii, jj
real*8 pxtemp, pytemp
real*8 x, y

do jj = 1, cpp
  ii = rank*cpp+jj
    do j=1,long
       x = in1(j ,2)
       y = in1(j, 3)
       pxtemp = (x + posicion(ii, 1))/delta     ! en el nuevo sistema
            do while (pxtemp.gt.dimx)
             pxtemp = pxtemp - dimx
            enddo
            do while (pxtemp.lt.1.0d-20)
              pxtemp = pxtemp + dimx
            enddo     
            px(il, j, jj) = int(pxtemp) + 1 !!! 

            pytemp = (y + posicion(ii, 2))/delta ! en el nuevo sistema
            do while (pytemp.gt.dimy)
              pytemp = pytemp - dimy
            enddo
            do while (pytemp.lt.1.0d-20)
              pytemp = pytemp + dimy
            enddo
            py(il, j, jj) = int(pytemp) + 1
            pz(il,j, jj)=int((in1(j,1))/delta)+1
     enddo
enddo

return
end
      



