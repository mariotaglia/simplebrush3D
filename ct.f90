subroutine ctinit
use chainsdat
implicit none
real*8 r,m
integer lado,i,j,ii

ct(1:ncha/2)=1
ct(ncha/2+1:ncha)=2

return
end

