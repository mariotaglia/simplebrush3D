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

!do jj = 1, cpp
!  ii = rank*cpp+jj
!    do j=1,long
!       x = in1(j ,2)
!       y = in1(j, 3)
!       pxtempA = (x + posicion(ii, 1))/delta     ! en el nuevo sistema
 !           do while (pxtempA.gt.dimx)
 !            pxtempA = pxtempA - dimx
 !           enddo
 !           do while (pxtempA.lt.1.0d-20)
  !            pxtempA = pxtempA + dimx
  !          enddo     
  !          pxA(il, j, jj) = int(pxtempA) + 1 !!! 

 !           pytempA = (y + posicion(ii, 2))/delta ! en el nuevo sistema
 !           do while (pytempA.gt.dimy)
 !             pytempA = pytempA - dimy
 !           enddo
 !           do while (pytempA.lt.1.0d-20)
 !             pytempA = pytempA + dimy
 !           enddo
 !           pyA(il, j, jj) = int(pytempA) + 1
 !           pzA(il,j, jj)=int((in1(j,1))/delta)+1
!     enddo
!enddo
!do jj = 1, cpp
!  ii = rank*cpp+jj
!    do j=1,long
!       x = in1(j ,2)
!       y = in1(j, 3)
!       pxtempB = (x + posicion(ii, 1))/delta     ! en el nuevo sistema
!            do while (pxtempB.gt.dimx)
!             pxtempB = pxtempB - dimx
!            enddo
!            do while (pxtempB.lt.1.0d-20)
!              pxtempB = pxtempB + dimx
!            enddo     
!            pxB(il, j, jj) = int(pxtempB) + 1 !!! 

!            pytempB = (y + posicion(ii, 2))/delta ! en el nuevo sistema
!            do while (pytempB.gt.dimy)
!              pytempB = pytempB - dimy
!            enddo
!            do while (pytempB.lt.1.0d-20)
!              pytempB = pytempB + dimy
!            enddo
!            pyB(il, j, jj) = int(pytempB) + 1
!            pzB(il,j, jj)=int((in1(j,1))/delta)+1
!     enddo
!enddo

return
end
      



