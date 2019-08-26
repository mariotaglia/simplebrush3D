subroutine ct(ii)

use system
use MPI
use chainsdat
use conformations
use const
implicit none

integer ii

do ii=1,nchas
	ct(ii)= mod(ii,2)+1
enddo

return
end
