module xxm
integer, parameter :: xx = 10
endmodule

module yym
use xxm
integer yy(xx)
end module

use xxm
use yym
yy(1) = 1
print*, xx
end
