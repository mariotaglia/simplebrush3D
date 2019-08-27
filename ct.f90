subroutine ctinit
use chainsdat
implicit none
real*8 r,m
integer lado,i,j,ii

	m=ncha*1.0
	r=sqrt(m)
	lado=int(r)
	do j=1,lado
		do i=1,lado
			ii=lado*(j-1)+i
			ct(ii)= ((-1)**(i+j)+1)/2
			print*,ii,ct(ii)
		enddo
	enddo

return
end

