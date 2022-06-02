real(8) function epslon (x)
  double precision :: x

  !     estimate unit roundoff in quantities of size x.

  double precision :: a,b,c,eps

  !     this program should function properly on all systems
  !     satisfying the following two assumptions,
  !        1.  the base used in representing floating point
  !            numbers is not a power of three.
  !        2.  the quantity  a  in statement 10 is represented to
  !            the accuracy used in floating point variables
  !            that are stored in memory.
  !     the statement number 10 and the go to 10 are intended to
  !     force optimizing compilers to generate code satisfying
  !     assumption 2.
  !     under these assumptions, it should be true that,
  !            a  is not exactly equal to four-thirds,
  !            b  has a zero for its last bit or digit,
  !            c  is not exactly equal to one,
  !            eps  measures the separation of 1.0 from
  !                 the next larger floating point number.
  !     the developers of eispack would appreciate being informed
  !     about any systems where these assumptions do not hold.

  !     this version dated 4/6/83.

  a = 4.0d0/3.0d0
10 b = a - 1.0d0
  c = b + b + b
  eps = dabs(c-1.0d0)
  if (eps == 0.0d0) go to 10
  epslon = eps*dabs(x)
  return
END function epslon

