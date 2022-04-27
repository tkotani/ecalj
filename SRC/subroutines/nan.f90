program ieee_naninf
use, intrinsic :: ieee_arithmetic
implicit none
real :: NaN, Inf
NaN = ieee_value(0.0, ieee_quiet_nan)
Inf = ieee_value(0.0, ieee_positive_inf)
print *, Inf, NaN
end program ieee_naninf
