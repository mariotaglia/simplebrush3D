subroutine allocatecpp
use fields_fkfun
use conformations
use chainsdat
implicit none

ALLOCATE(px(cuantas, long, cpp))
ALLOCATE(py(cuantas, long, cpp))
ALLOCATE(pz(cuantas, long, cpp))
ALLOCATE(pro(cuantas, cpp))
end subroutine
