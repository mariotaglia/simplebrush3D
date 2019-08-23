subroutine allocatecpp
use fields_fkfun
use conformations
use chainsdat
implicit none

ALLOCATE(px(cuantas, long, cpp))
ALLOCATE(py(cuantas, long, cpp))
ALLOCATE(pz(cuantas, long, cpp))
ALLOCATE(pro(cuantas, cpp))

!ALLOCATE(pxA(cuantas, long, cpp))
!ALLOCATE(pyA(cuantas, long, cpp))
!ALLOCATE(pzA(cuantas, long, cpp))
!ALLOCATE(proA(cuantas, cpp))
!ALLOCATE(pxB(cuantas, long, cpp))
!ALLOCATE(pyB(cuantas, long, cpp))
!ALLOCATE(pzB(cuantas, long, cpp))
!ALLOCATE(proB(cuantas, cpp))

end subroutine
