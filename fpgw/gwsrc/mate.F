c---------------------
c mathematics library
c---------------------

      subroutine acv (c,v,
     d                n,
     o                w)

c 92.01.17
c add a constant to a vector w(i)= c + v(i)

      implicit real*8 (a-h,o-z)
      dimension v(n)
      dimension w(n)

      do       i = 1,n
        w(i)       = c + v(i)
      end do

      return
      end
c-------------------------------------------------------------------
      subroutine acm (c,v,
     d                n)

c 92.02.20
c add a constant to the diagonal elements of a matrix
c v(i,i)= c + v(i,i)

      implicit real*8 (a-h,o-z)
      dimension v(n,n)

      do       i = 1,n
        v(i,i)     = c + v(i,i)
      end do

      return
      end
c-------------------------------------------------------------------
      double precision function alagr2 (x,xi,fi)

c 92.03.02
c 92.04.10 from alagr3
c two-point interpolation
c given a function fi at two points xi, the routine interpolates
c the function at x
c f(x) = [ (x-x2)/(x1-x2) ] f1
c      + [ (x-x1)/(x2-x1) ] f2

c x  = the point at which the function is to be interpolated
c xi(2) = points where the function is given
c fi(2) = the function at xi

      implicit real*8 (a-h,o-z)
      dimension xi(2),fi(2)

      xx1        = x-xi(1)
      xx2        = x-xi(2)
      x12        = xi(1)-xi(2)
      alagr2     = (xx2*fi(1) - xx1*fi(2))/x12

      return
      end
c-----------------------------------------------------------------
      double precision function alagr3 (x,xi,fi)

c 92.03.02
c three-point interpolation with arbitrary mesh
c given a function fi at three points xi, the routine interpolates
c the function at x
c f(x) = [ { (x-x2)(x-x3) } / { (x1-x2)(x1-x3) } ] f1
c      + [ { (x-x1)(x-x3) } / { (x2-x1)(x2-x3) } ] f2
c      + [ { (x-x1)(x-x2) } / { (x3-x1)(x3-x2) } ] f3

c x  = the point at which the function is to be interpolated
c xi(3) = points where the function is given
c fi(3) = the function at xi

      implicit real*8 (a-h,o-z)
      dimension xi(3),fi(3)

      xx1        = x-xi(1)
      xx2        = x-xi(2)
      xx3        = x-xi(3)
      x12        = xi(1)-xi(2)
      x13        = xi(1)-xi(3)
      x23        = xi(2)-xi(3)
      alagr3     = xx2*xx3*fi(1)/(x12*x13)
     .           - xx1*xx3*fi(2)/(x12*x23)
     .           + xx1*xx2*fi(3)/(x13*x23)

      return
      end
c-----------------------------------------------------------------
      double precision function alagr4 (xi,f,x)

c given a function at four arbitrary points x1,x2,x3,x4,
c alagr4 gives the interpolated value at x, which can be anywhere
c x1,x2,x3,x4 do not have to be ordered
c f1,f2,f3,f4 must be the values of the function at x1,x2,x3,x4

c f(x) = [ {(x-x2)(x-x3)(x-x4)} / {(x1-x2)(x1-x3)(x1-x4)} ] f1
c      + [ {(x-x1)(x-x3)(x-x4)} / {(x2-x1)(x2-x3)(x2-x4)} ] f2
c      + [ {(x-x1)(x-x2)(x-x4)} / {(x3-x1)(x3-x2)(x3-x4)} ] f3
c      + [ {(x-x1)(x-x2)(x-x3)} / {(x4-x1)(x4-x2)(x4-x3)} ] f3

      implicit real*8(a-h,o-z)
      dimension xi(4),f(4)

      xx1        = x-xi(1)
      xx2        = x-xi(2)
      xx3        = x-xi(3)
      xx4        = x-xi(4)
      xx23       = xx2*xx3
      xx14       = xx1*xx4
      x12        = xi(1)-xi(2)
      x13        = xi(1)-xi(3)
      x14        = xi(1)-xi(4)
      x23        = xi(2)-xi(3)
      x24        = xi(2)-xi(4)
      x34        = xi(3)-xi(4)
      alagr4     = xx23*xx4*f(1)/(x12*x13*x14)
     .           - xx14*xx3*f(2)/(x12*x23*x24)
     .           + xx14*xx2*f(3)/(x13*x23*x34)
     .           - xx23*xx1*f(4)/(x14*x24*x34)

      return
      end
c-------------------------------------------------------------------
      double precision function cellvol (bas)

c calculates the volume of a parallelepiped formed by bas
c vol = |b1 . (b2 x b3)|

      implicit real*8(a-h,o-z)
      dimension bas(3,3)
      dimension w1(3)

      call cross_x  (bas(1,2),bas(1,3),w1)
      cellvol    = dabs( vdv(bas(1,1),w1,3) )

      return
      end
c-------------------------------------------------------------------
      subroutine clagr4 (xi,x,
     o c )

c given a function at four arbitrary points x1,x2,x3,x4,
c clagr4 gives the coefficients of interpolation at x,
c which can be anywhere
c x1,x2,x3,x4 do not have to be ordered but when calculating
c the function at x, f(x) = S[i=1,4] c(i) f(i)
c f1,f2,f3,f4 must be the values of the function at x1,x2,x3,x4

c f(x) = [ {(x-x2)(x-x3)(x-x4)} / {(x1-x2)(x1-x3)(x1-x4)} ] f1
c      + [ {(x-x1)(x-x3)(x-x4)} / {(x2-x1)(x2-x3)(x2-x4)} ] f2
c      + [ {(x-x1)(x-x2)(x-x4)} / {(x3-x1)(x3-x2)(x3-x4)} ] f3
c      + [ {(x-x1)(x-x2)(x-x3)} / {(x4-x1)(x4-x2)(x4-x3)} ] f3

      implicit real*8(a-h,o-z)
      dimension xi(4),c(4)

      xx1        = x-xi(1)
      xx2        = x-xi(2)
      xx3        = x-xi(3)
      xx4        = x-xi(4)
      xx23       = xx2*xx3
      xx14       = xx1*xx4
      x12        = xi(1)-xi(2)
      x13        = xi(1)-xi(3)
      x14        = xi(1)-xi(4)
      x23        = xi(2)-xi(3)
      x24        = xi(2)-xi(4)
      x34        = xi(3)-xi(4)
      c(1)       = xx23*xx4 / (x12*x13*x14)
      c(2)       =-xx14*xx3 / (x12*x23*x24)
      c(3)       = xx14*xx2 / (x13*x23*x34)
      c(4)       =-xx23*xx1 / (x14*x24*x34)

      return
      end
c-----------------------------------------------------------------
      subroutine cubstar(xp,yp,zp,
     o                   xst,yst,zst,nstar)

c finds cubic stars of a given vector
c xp,yp,zp    = the given vector
c xst,yst,zst = the stars
c nstar       = number of stars

      implicit real*8(a-h,o-z)
      dimension xst(48),yst(48),zst(48),xyz(3)
c this equivalence is used to avoid problems with the compiler
      equivalence (x,xyz(1)),(y,xyz(2)),(z,xyz(3))
      data delta/1.d-6/


c     x         = dabs(xp)
c     y         = dabs(yp)
c     z         = dabs(zp)
      xyz(1)    = dabs(xp)
      xyz(2)    = dabs(yp)
      xyz(3)    = dabs(zp)
      nx        = idnint( x/(x+delta) )
      ny        = idnint( y/(y+delta) )
      nz        = idnint( z/(z+delta) )
      i         = idnint( (x-y)/(x-y+delta) )
      j         = idnint( (y-z)/(y-z+delta) )
      k         = idnint( (z-x)/(z-x+delta) )
      inv       = iabs( i+j+k-1 )
      iper      = ( (i+j+k)/2 )*2 + 1
      nstar     = 0

ctemp
c     iftemp = ifile('TEMP')
c     write(iftemp,6000)xp,yp,zp
c     write(iftemp,6100)nx,ny,nz,i,j,k,inv,iper
c6000 format(1x,5d13.6)
c6100 format(1x,8i4)

      do 1    i = 1,inv
        do 2    j = 1,iper
          do 3   ix = nx,-nx,-2
            xi = ix
          do 3   iy = ny,-ny,-2
            yi = iy
          do 3   iz = nz,-nz,-2
            zi = iz
            nstar     = nstar + 1
c     xst(nstar)= dfloat(ix) * x
c     yst(nstar)= dfloat(iy) * y
c   3 zst(nstar)= dfloat(iz) * z
c     xst(nstar)= dble(ix) * x
c     yst(nstar)= dble(iy) * y
c     zst(nstar)= dble(iz) * z
c     write (iftemp,6000)x,y,z
c     xst(nstar)= xi * x
c     yst(nstar)= yi * y
c     zst(nstar)= zi * z
            xst(nstar)= xi * xyz(1)
            yst(nstar)= yi * xyz(2)
            zst(nstar)= zi * xyz(3)
ctemp
c     write(iftemp,6200)nstar,ix,iy,iz,
c    .xst(nstar),yst(nstar),zst(nstar)
c6200 format(1x,4i4,3d13.6)
    3     continue

          k         = nx
          xt        = x
          nx        = nz
          nz        = ny
          ny        = k
          x         = z
          z         = y
          y         = xt
    2   continue
        k         = nx
        xt        = x
        nx        = ny
        ny        = k
        x         = y
        y         = xt
    1 continue
c     write(6,601)nstar
c     write(6,600)(xst(i),yst(i),zst(i),i=1,nstar)
c 600 format(2x,3d13.5)
c 601 format(2x,i4/)
      return
      end
c--------------------------------------------------------------------
      subroutine cv (c,v,n,
     o w )

c forms w(i) = c * v(i)

      implicit real*8(a-h,o-z)
      dimension v(n)
      dimension w(n)

      do       i = 1,n
        w(i)       = c*v(i)
      end do

      return
      end
c-----------------------------------------------------------------
      subroutine cva (c,v,n,
     o w )

c forms w(i) = w(i) + c * v(i)

      implicit real*8(a-h,o-z)
      dimension v(n)
      dimension w(n)

      do       i = 1,n
        w(i)       = w(i) + c*v(i)
      end do

      return
      end
c-------------------------------------------------------------------
      subroutine cvs (c,v,n,
     o w )

c forms w(i) = c * v(i)
c single precision

      implicit real*8(a-h,o-z)
      real*4    v(n),c
      real*4    w(n)

      do       i = 1,n
        w(i)       = c*v(i)
      end do

      return
      end
c-----------------------------------------------------------------
      subroutine dbles (dv,
     d                  n,
     o                  sv)

c 92.02.07
c convert a double precision array dv to single sv

      implicit real*8 (a-h,o-z)
      real*8 dv(n)
      real*4 sv(n)

      do       i = 1,n
        sv(i)      = sngl(dv(i))
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine dmv(nr,nc,a,ndim,v,av)
      implicit real*8(a-h,o-z)

c matrix * vector
c av(m) = sum(n=1,nc) a(m,n)*v(n) ,m=1,nr
c input:
c nr = dimension of output vector
c nc = dimension summed
c a  = matrix to be multiplied
c ndim = leading dimension of a
c v  = vector to be multiplied
c av = resulting vector

      dimension a(ndim,nc),v(nc),av(nr)

      if(nr .gt. ndim)stop 'dmv: row of a too large'
      do 10    m = 1,nr
        tem        = 0.d0
        do 20    n = 1,nc
          tem        = tem + a(m,n)*v(n)
   20   continue
        av(m)      = tem
   10 continue
      return
      end
c-----------------------------------------------------------------
      subroutine finter (f,nr,
     i                   clagr,istart,nrg,
     o fint)

c interpolate a function f at gaussian points

c f   = function at equal mesh
c nr  = no. equal mesh points
c clagr  = lagrange 4-point interpolation coefficients
c istart = index for starting point of interpolation, see inter.f
c nrg    = no. gaussian points

c fint   = interpolated function

      implicit double precision(a-h,o-z)
      dimension f(nr),clagr(4,nrg),istart(nrg)
      dimension fint(nrg)

      do       i = 1,nrg
        is         = istart(i) - 1
        tem1       = 0.d0
        do       m = 1,4
          tem1       = tem1 + clagr(m,i)*f(is+m)
        end do
        fint(i)    = tem1
      end do

      return
      end
c-----------------------------------------------------------------
      subroutine gauss(n,a,b,x,w)
c
c subroutine gauss for gaussian integration   76.08.12.
c the integral of f(x) from a to b is given by the sum
c over i from 1 to n of w(i)*f(x(i)).
c maximum n is 48
c
      implicit real*8(a-h,o-z)
      dimension x(n),w(n),e06(3),w06(3),e10(5),w10(5),e12(6),w12(6),
     d          e16(8),w16(8),e20(10),w20(10),e24(12),w24(12),
     d          e32(16),w32(16),e40(20),w40(20),e48(24),w48(24)
      data e06/.238619186083197d0,.661209386466265d0,.932469514203152d0/
      data w06/.467913934572691d0,.360761573048139d0,.171324492379170d0/
      data e10/.148874338981631d0,.433395394129247d0,.679409568299024d0,
     d         .865063366688985d0,.973906528517172d0/
      data w10/.295524224714753d0,.269266719309996d0,.219086362515982d0,
     d         .149451349150581d0,.066671344308688d0/
      data e12/.125233408511469d0,.367831498998180d0,.587317954286617d0,
     d         .769902674194305d0,.904117256370475d0,.981560634246719d0/
      data w12/.249147045813403d0,.233492536538355d0,.203167426723066d0,
     d         .160078328543346d0,.106939325995318d0,.047175336386512d0/
      data e16/.0950125098376374d0,.2816035507792589d0,
     d         .4580167776572274d0,.6178762444026437d0,
     d         .7554044083550030d0,.8656312023878317d0,
     d         .9445750230732326d0,.9894009349916499d0/
      data w16/.1894506104550685d0,.1826034150449236d0,
     d         .1691565193950025d0,.1495959888165767d0,
     d         .1246289712555339d0,.0951585116824928d0,
     d         .0622535239386479d0,.0271524594117541d0/
      data e20/.0765265211334973d0,.2277858511416451d0,
     d         .3737060887154196d0,.5108670019508271d0,
     d         .6360536807265150d0,.7463319064601508d0,
     d         .8391169718222188d0,.9122344282513259d0,
     d         .9639719272779138d0,.9931285991850949d0/
      data w20/.1527533871307259d0,.1491729864726037d0,
     d         .1420961093183821d0,.1316886384491766d0,
     d         .1181945319615184d0,.1019301198172404d0,
     d         .0832767415767047d0,.0626720483341091d0,
     d         .0406014298003869d0,.0176140071391521d0/
      data e24/.0640568928626056d0,.1911188674736163d0,
     d         .3150426796961634d0,.4337935076260451d0,
     d         .5454214713888395d0,.6480936519369756d0,
     d         .7401241915785544d0,.8200019859739029d0,
     d         .8864155270044010d0,.9382745520027328d0,
     d         .9747285559713095d0,.9951872199970214d0/
      data w24/.1279381953467522d0,.1258374563468283d0,
     d         .1216704729278034d0,.1155056680537256d0,
     d         .1074442701159656d0,.0976186521041139d0,
     d         .0861901615319533d0,.0733464814110803d0,
     d         .0592985849154368d0,.0442774388174198d0,
     d         .0285313886289337d0,.0123412297999872d0/
      data e32/.0483076656877383d0,.1444719615827965d0,
     d         .2392873622521371d0,.3318686022821276d0,
     d         .4213512761306353d0,.5068999089322294d0,
     d         .5877157572407623d0,.6630442669302152d0,
     d         .7321821187402897d0,.7944837959679424d0,
     d         .8493676137325700d0,.8963211557660521d0,
     d         .9349060759377397d0,.9647622555875064d0,
     d         .9856115115452683d0,.9972638618494816d0/
      data w32/.0965400885147278d0,.0956387200792749d0,
     d         .0938443990808046d0,.0911738786957639d0,
     d         .0876520930044038d0,.0833119242269468d0,
     d         .0781938957870703d0,.0723457941088485d0,
     d         .0658222227763618d0,.0586840934785355d0,
     d         .0509980592623762d0,.0428358980222267d0,
     d         .0342738629130214d0,.0253920653092621d0,
     d         .0162743947309057d0,.0070186100094701d0/
      data e40/.0387724175060508d0,.1160840706752552d0,
     d         .1926975807013711d0,.2681521850072537d0,
     d         .3419940908257585d0,.4137792043716050d0,
     d         .4830758016861787d0,.5494671250951282d0,
     d         .6125538896679802d0,.6719566846141795d0,
     d         .7273182551899271d0,.7783056514265194d0,
     d         .8246122308333117d0,.8659595032122595d0,
     d         .9020988069688743d0,.9328128082786765d0,
     d         .9579168192137917d0,.9772599499837743d0,
     d         .9907262386994570d0,.9982377097105592d0/
      data w40/.0775059479784248d0,.0770398181642480d0,
     d         .0761103619006262d0,.0747231690579683d0,
     d         .0728865823958041d0,.0706116473912868d0,
     d         .0679120458152339d0,.0648040134566010d0,
     d         .0613062424929289d0,.0574397690993916d0,
     d         .0532278469839368d0,.0486958076350722d0,
     d         .0438709081856733d0,.0387821679744720d0,
     d         .0334601952825478d0,.0279370069800234d0,
     d         .0222458491941670d0,.0164210583819079d0,
     d         .0104982845311528d0,.0045212770985332d0/
      data e48/
     d 0.03238 01709 62869 326033d0 , 0.09700 46992 09462 698930d0,
     d 0.16122 23560 68891 718056d0 , 0.22476 37903 94689 061225d0,
     d 0.28736 24873 55455 576736d0 , 0.34875 58862 92160 738160d0,
     d 0.40868 64819 90716 729916d0 , 0.46690 29047 50958 404545d0,
     d 0.52316 09747 22233 033678d0 , 0.57722 47260 83972 703818d0,
     d 0.62886 73967 76513 623995d0 , 0.67787 23796 32663 905212d0,
     d 0.72403 41309 23814 654674d0 , 0.76715 90325 15740 339254d0,
     d 0.80706 62040 29442 627083d0 , 0.84358 82616 24393 530711d0,
     d 0.87657 20202 74247 885906d0 , 0.90587 91367 15569 672822d0,
     d 0.93138 66907 06554 333114d0 , 0.95298 77031 60430 860723d0,
     d 0.97059 15925 46247 250461d0 , 0.98412 45837 22826 857745d0,
     d 0.99353 01722 66350 757548d0 , 0.99877 10072 52426 118601d0/
      data w48/
     d 0.06473 76968 12683 922503d0 , 0.06446 61644 35950 082207d0,
     d 0.06392 42385 84648 186624d0 , 0.06311 41922 86254 025567d0,
     d 0.06203 94231 59892 663904d0 , 0.06070 44391 65893 880053d0,
     d 0.05911 48396 98398 635746d0 , 0.05727 72921 00403 215705d0,
     d 0.05519 95036 99984 162868d0 , 0.05289 01894 85193 667096d0,
     d 0.05035 90355 53854 474958d0 , 0.04761 66584 92490 474826d0,
     d 0.04467 45608 56694 280419d0 , 0.04154 50829 43464 749214d0,
     d 0.03824 13510 65830 706317d0 , 0.03477 72225 64770 438893d0,
     d 0.03116 72278 32798 088902d0 , 0.02742 65097 08356 948200d0,
     d 0.02357 07608 39324 379141d0 , 0.01961 61604 57355 527814d0,
     d 0.01557 93157 22943 848728d0 , 0.01147 72345 79234 539490d0,
     d 0.00732 75539 01276 262102d0 , 0.00315 33460 52305 838633d0/
      e1 = .5d0*(b - a)
      e2 = e1 + a
      if(n.eq. 6) go to  6
      if(n.eq.10) go to 10
      if(n.eq.12) go to 12
      if(n.eq.16) go to 16
      if(n.eq.20) go to 20
      if(n.eq.24) go to 24
      if(n.eq.32) go to 32
      if(n.eq.40) go to 40
      if(n.eq.48) go to 48
      xl = dsqrt(-e06(1))
      stop 'gauss: wrong no. points'
    6 continue
      do  7 i = 1,3
        ak = e1*e06(i)
        x(i+3) = e2 + ak
        x(4-i) = e2 - ak
        w(i+3) = w06(i)*e1
        w(4-i) = w06(i)*e1
    7 continue
      return
   10 continue
      do 11 i = 1,5
        ak = e1*e10(i)
        x(i+5) = e2 + ak
        x(6-i) = e2 - ak
        w(i+5) = w10(i)*e1
        w(6-i) = w10(i)*e1
   11 continue
      return
   12 continue
      do 13 i = 1,6
        ak = e1*e12(i)
        x(i+6) = e2 + ak
        x(7-i) = e2 - ak
        w(i+6) = w12(i)*e1
        w(7-i) = w12(i)*e1
   13 continue
      return
   16 continue
      do 17 i = 1,8
        ak = e1*e16(i)
        x(i+8) = e2 + ak
        x(9-i) = e2 - ak
        w(i+8) = w16(i)*e1
        w(9-i) = w16(i)*e1
   17 continue
      return
   20 continue
      do 21 i = 1,10
        ak = e1*e20(i)
        x(i+10) = e2 + ak
        x(11-i) = e2 - ak
        w(i+10) = w20(i)*e1
        w(11-i) = w20(i)*e1
   21 continue
      return
   24 continue
      do 25 i = 1,12
        ak = e1*e24(i)
        x(i+12) = e2 + ak
        x(13-i) = e2 - ak
        w(i+12) = w24(i)*e1
        w(13-i) = w24(i)*e1
   25 continue
      return
   32 continue
      do 33 i = 1,16
        ak = e1*e32(i)
        x(i+16) = e2 + ak
        x(17-i) = e2 - ak
        w(i+16) = w32(i)*e1
        w(17-i) = w32(i)*e1
   33 continue
      return
   40 continue
      do 41 i = 1,20
        ak = e1*e40(i)
        x(i+20) = e2 + ak
        x(21-i) = e2 - ak
        w(i+20) = w40(i)*e1
        w(21-i) = w40(i)*e1
   41 continue
      return
   48 continue
      do 49 i = 1,24
        ak      = e1*e48(i)
        x(i+24) = e2 + ak
        x(25-i) = e2 - ak
        w(i+24) = w48(i)*e1
        w(25-i) = w48(i)*e1
   49 continue
      return
      end
c-------------------------------------------------------------------
c     subroutine iinit(iw,ldim)
c
c set array iw with ldim elements to zero
c
c     implicit real*8(a-h,o-z)
c     dimension iw(ldim)
c
c     do       i = 1,ldim
c     iw(i)      = 0
c     end do

c     return
c     end
c------------------------------------------------------------------
      integer function ivsum (iv,n)

c sum the elements of an integer array
      implicit real*8(a-h,o-z)
      dimension iv(n)

      ivsum      = 0
      do       i = 1,n
        ivsum      = ivsum + iv(i)
      end do

      return
      end
c-------------------------------------------------------------------
      integer function lowtri (ldim,irow,icol)

c index for lower triangle
c    1
c    2  6
c    3  7  10
c    4  8  11  13
c    5  9  12  14  15
c e.g. ldim=5 irow=4 icol=3  ==> lowtri=11
c irow must be .ge. icol

      implicit real*8(a-h,o-z)

      if(irow .lt. icol)stop 'lowtri: irow < icol'
      lowtri     =  ldim*(icol-1) - (icol*(icol-1))/2 + irow

      return
      end
c---------------------------------------------------------------------
      subroutine madd (a,lda,b,ldb,nrow,ncol,ldc,
     o c)

c add matrices a + b = c

c lda,ldb,ldc = leading dimension of a,b and c
c nrow,ncol   = dimension of row and column of c

      implicit double precision (a-h,o-z)
      dimension a(lda,1),b(ldb,1)
      dimension c(ldc,1)

      do       j = 1,ncol
        do       i = 1,nrow
          c(i,j)     = a(i,j) + b(i,j)
        end do
      end do

      return
      end
c---------------------------------------------------------------------
      integer function maxim (num,ldim)

c find the maximum value of a set of integers

      implicit real*8(a-h,o-z)
      dimension num(ldim)

      maxim      = num(1)
      do       i = 2,ldim
        if(num(i) .gt. maxim) maxim = num(i)
      end do

      return
      end
c-------------------------------------------------------------------
      subroutine mfull (tri,n,sym,
     o full )

c 91.11.29
c recovers a matrix written in lower triangular form, see lowtri.f

c tri  = lower triangular matrix
c n    = dimension of matrix
c sym  = logical: .true.   ==> symmetric
c                 .false.  ==> anti-symmetric

c full = full matrix

      implicit double precision (a-h,o-z)

      logical   sym
      dimension tri( (n*(n+1))/2 )
      dimension full(n,n)

      i          = 0
      do      ic = 1,n
        do      ir = ic,n
          i          = i + 1
          full(ir,ic)= tri(i)
        end do
      end do

      if(sym)then
        do      ic = 1,n
          do      ir = ic+1,n
            full(ic,ir)= full(ir,ic)
          end do
        end do
      else
        do      ic = 1,n
          do      ir = ic+1,n
            full(ic,ir)= -full(ir,ic)
          end do
        end do
      endif

      return
      end
c-------------------------------------------------------------------
      subroutine mfullsd (tri,n,sym,
     o full )

c 91.11.29, 92.02.09
c recovers a matrix written in lower triangular form, see lowtri.f
c convert from single to double precision

c tri  = lower triangular matrix in single precision
c n    = dimension of matrix
c sym  = logical: .true.   ==> symmetric
c                 .false.  ==> anti-symmetric

c full = full matrix in double precision

      implicit double precision (a-h,o-z)

      logical   sym
      real*4    tri( (n*(n+1))/2 )
      dimension full(n,n)

      i          = 0
      do      ic = 1,n
        do      ir = ic,n
          i          = i + 1
          full(ir,ic)= dble( tri(i) )
        end do
      end do

      if(sym)then
        do      ic = 1,n
          do      ir = ic+1,n
            full(ic,ir)= full(ir,ic)
          end do
        end do
      else
        do      ic = 1,n
          do      ir = ic+1,n
            full(ic,ir)= -full(ir,ic)
          end do
        end do
      endif

      return
      end
c-------------------------------------------------------------------
      subroutine minv(s,ldim,n,
     w                work,ipvt,
     o si)

c invert a real (symmetric ?) matrix s

c s    = matrix to be inverted
c ldim = leading dimension of s
c n    = dimension of the problem
c work,ipvt = work arrays

c si   = inverse of s

      implicit real*8(a-h,o-z)
      dimension s(ldim,n),
     w          work(n),ipvt(n)
      dimension si(n,n)
      dimension det(2)

      do         i = 1,n
        call dcopy   (n,s(1,i),1,si(1,i),1)
      end do

      call dgefa   (si,n,n,ipvt,info)
      call dgedi   (si,n,n,ipvt,det,work,01)

      return
      end
c-------------------------------------------------------------------
      subroutine minvc (r,c,
     d                  ldim,n,
     w                  work,ipvt,w1,w2,
     o ri,ci )

c 91.11.29
c invert a complex matrix

c r,c   = real and imaginary parts of the matrix
c ldim  = leading dimension of matrix
c n     = dimension of the problem
c work,ipvt,w1,w2 = work arrays

c ri,ci = real and imaginary parts of the inverse matrix

      implicit double precision (a-h,o-z)
      dimension r(ldim,n),c(ldim,n),
     w          work(n),ipvt(n),w1(n,n),w2(n,n)
      dimension ri(n,n),ci(n,n)

c invert real part
      call minv    (r,ldim,n,
     w              work,ipvt,
     o              w1)

c real part of inverse
      call mmul    (c,ldim,w1,n,n,n,n,n,
     o              w2 )
      call mmul    (w2,n,c,ldim,n,n,n,n,
     o              ri )
      call madd    (ri,n,r,ldim,n,n,n,
     o              ci )
      call minv    (ci,n,n,
     w              work,ipvt,
     o              ri )

c imaginary part of inverse
      call mmul    (ri,n,w2,n,n,n,n,n,
     o              ci )
      call cv      (-1.d0,ci,n*n,
     o              ci )

      return
      end
c-------------------------------------------------------------------
      subroutine mms (a,lda,b,ldb,nr,nm,nc,ldc,
     o c)

c single precision matrix multiplication
c  c = a b
c lda = leading dimension of a
c ldb = leading dimension of b
c ldc = leading dimension of c
c nr  = no. rows of c
c nc  = no. columns of c
c nm  = summed dimension

      implicit real*8(a-h,o-z)
      real*4    a(lda,nm),b(ldb,nc),c(ldc,nc),tem

      if(nr .gt. lda)stop 'mms: row of a too large'
      if(nm .gt. ldb)stop 'mms: row of b too large'
      if(nr .gt. ldc)stop 'mms: row of c too large'

      do 10   i = 1,nr
      do 10   j = 1,nc
        tem       = 0.0
        do 20   k = 1,nm
          tem       = tem + a(i,k)*b(k,j)
   20   continue
        c(i,j)    = tem
   10 continue

      return
      end
c------------------------------------------------------------------
      subroutine mmul(a,lda,b,ldb,nr,nm,nc,ldc,
     o c)

c matrix multiplication
c  c = a b
c lda = leading dimension of a
c ldb = leading dimension of b
c ldc = leading dimension of c
c nr  = no. rows of c
c nc  = no. columns of c
c nm  = summed dimension

      implicit real*8(a-h,o-z)
      dimension a(lda,nm),b(ldb,nc),c(ldc,nc)

      if(nr .gt. lda)stop 'mmul: row of a too large'
      if(nm .gt. ldb)stop 'mmul: row of b too large'
      if(nr .gt. ldc)stop 'mmul: row of c too large'

c     do 10   i = 1,nr
c     do 10   j = 1,nc
c     tem       = 0.d0
c     do 20   k = 1,nm
c  20 tem       = tem + a(i,k)*b(k,j)
c  10 c(i,j)    = tem

      do      j = 1,nc

        do      i = 1,nr
          c(i,j)    = 0.d0
        enddo

        do      k = 1,nm
          bkj       = b(k,j)
          do      i = 1,nr
            c(i,j)    = c(i,j) + bkj*a(i,k)
          enddo
        enddo

      enddo

      return
      end
c------------------------------------------------------------------
      subroutine mmulc (ra,ca,lda,
     i                  rb,cb,ldb,
     i                  nrow,nmul,ncol,ldc,
     o rc,cc)

c 91.11.29
c multiply two complex matrices a b = c

c ra,ca = real and imaginary parts of a
c rb,cb =                             b
c lda,ldb,ldc = leading dimensions of a,b and c
c nrow,ncol   = no. rows and coulmns of c
c nmul  = no. contractions

c rc,cc = real and imaginary parts of c

      implicit double precision (a-h,o-z)

      dimension ra(lda,1),ca(lda,1),
     i          rb(ldb,1),cb(ldb,1)
      dimension rc(ldc,1),cc(ldc,1)

      if(nrow .gt. lda) stop 'mmulc: lda too small'
      if(nmul .gt. ldb) stop 'mmulc: ldb too small'
      if(nmul .gt. ldc) stop 'mmulc: ldc too small'
c     do      ir = 1,nrow
c     do      ic = 1,ncol
c     rsum       = 0.d0
c     csum       = 0.d0
c     do       i = 1,nmul
c     rsum       = rsum + ra(ir,i)*rb(i,ic)
c    .                  - ca(ir,i)*cb(i,ic)
c     csum       = csum + ra(ir,i)*cb(i,ic)
c    .                  + ca(ir,i)*rb(i,ic)
c     end do
c     rc(ir,ic)  = rsum
c     cc(ir,ic)  = csum
c     end do
c     end do

      do      ic = 1,ncol

        do      ir = 1,nrow
          rc(ir,ic)  = 0.d0
          cc(ir,ic)  = 0.d0
        enddo

        do       i = 1,nmul
          rbic       = rb(i,ic)
          cbic       = cb(i,ic)
          do      ir = 1,nrow
            rc(ir,ic)  = rc(ir,ic) + rbic*ra(ir,i) - cbic*ca(ir,i)
            cc(ir,ic)  = cc(ir,ic) + cbic*ra(ir,i) + rbic*ca(ir,i)
          end do
        end do

      end do

      return
      end

c------------------------------------------------------------------
      subroutine mv (a,v,
     d               ldima,nc,nr,
     o               av)

c 92.02.12
c matrix * vector
c av(m) = sum(n=1,nc) a(m,n)*v(n) ,m=1,nr

c a  = matrix to be multiplied
c v  = vector to be multiplied
c ldima = leading dimension of a
c nc = dimension summed
c nr = dimension of output vector

c av = resulting vector

      implicit real*8(a-h,o-z)
      dimension a(ldima,nc),v(nc),av(nr)

      if(nr .gt. ldima)stop 'mv: row of a too large'
      do 10    m = 1,nr
        tem        = 0.d0
        do 20    n = 1,nc
          tem        = tem + a(m,n)*v(n)
   20   continue
        av(m)      = tem
   10 continue

      return
      end
c-----------------------------------------------------------------
      subroutine ord2in(i,j,io)

c order two  integers in increasing magnitude

      implicit real*8(a-h,o-z)
      dimension io(2)

      io(2)       = max0(i,j)
      io(1)       = i + j - io(2)

      return
      end
c--------------------------------------------------------------------
      subroutine ord3in(i,j,k,io)

c order three integers in increasing magnitude

      implicit real*8(a-h,o-z)
      dimension io(3)

      maxij       = max0(i,j)
      minij       = i + j - maxij
      io(3)       = max0(maxij,k)
      maxl        = maxij + k - io(3)
      io(1)       = min0(minij,maxl)
      io(2)       = minij + maxl - io(1)

      return
      end
c--------------------------------------------------------------------
      subroutine ord4in(i,j,k,l,io)

c order four integers in increasing magnitude

      implicit real*8(a-h,o-z)
      dimension io(4)

      maxij       = max0(i,j)
      minij       = i + j - maxij
      maxkl       = max0(k,l)
      minkl       = k + l - maxkl
      io(4)       = max0(maxij,maxkl)
      maxl        = maxij + maxkl - io(4)
      io(1)       = min0(minij,minkl)
      minl        = minij + minkl - io(1)
      io(3)       = max0(maxl,minl)
      io(2)       = maxl + minl - io(3)

      return
      end
c--------------------------------------------------------------------
      subroutine ordna(ldim,frac,v,ov,indx)
      implicit real*8(a-h,o-z)

c order a set of numbers in decreasing values
c input:
c ldim  = dimension of v
c v     = array where the numbers are stored
c vo    = ordered array
c indx(old index)  = new index

      dimension v(ldim),ov(ldim),indx(ldim)

      do 20    i = 1,ldim
        indx(i)    = 0
   20 continue
c     ov(1)      = -1.d-30
      ic         = 0
 1111 ymax       = -1.d30
      do 10    i = 1,ldim
        if(indx(i) .ne. 0)goto 10
        if(v(i) .gt. ymax)then
          ymax       = v(i)
          imax       = i
        endif
   10 continue
c     if(ymax .lt. frac*ov(1))return
      ic         = ic + 1
      indx(imax) = ic
      ov(ic)     = ymax
      if(ic .lt. ldim)goto 1111
      return
      end
c---------------------------------------------------------------------
      subroutine ordnai(ldim,v,ov,indx)

c order a set of numbers in increasing values

c ldim  = dimension of v
c v     = array where the numbers are stored

c ov    = ordered array
c indx(old index)  = new index

      implicit real*8(a-h,o-z)
      dimension v(ldim),ov(ldim),indx(ldim)

      do 20    i = 1,ldim
        indx(i)    = 0
   20 continue

      ic         = 0
 1111 ymin       = 1.d30
      do 10    i = 1,ldim
        if(indx(i) .ne. 0)goto 10
        if(v(i) .lt. ymin)then
          ymin       = v(i)
          imin       = i
        endif
   10 continue
      ic         = ic + 1
      indx(imin) = ic
      ov(ic)     = ymin
      if(ic .lt. ldim)goto 1111

      return
      end
c----------------------------------------------------------------------
      subroutine ordnai1 (v,
     d                    ldim,
     o ov,indx,indxn)

c 92.01.22
c order a set of numbers in increasing values

c ldim  = dimension of v
c v     = array where the numbers are stored

c ov    = ordered array
c indx(old index)  = new index
c indxn(new) = old

      implicit real*8(a-h,o-z)
      dimension v(ldim)
      dimension ov(ldim),indx(ldim),indxn(ldim)

      do 20    i = 1,ldim
        indxn(i)   = 0
        indx(i)    = 0
   20 continue

      ic         = 0
 1111 ymin       = 1.d30
      do 10    i = 1,ldim
        if(indx(i) .ne. 0)goto 10
        if(v(i) .lt. ymin)then
          ymin       = v(i)
          imin       = i
        endif
   10 continue
      ic         = ic + 1
      indx(imin) = ic
      indxn(ic)  = imin
      ov(ic)     = ymin
      if(ic .lt. ldim)goto 1111

      return
      end
c----------------------------------------------------------------------
      double precision function pilrntz(x,delta)

c calculates pi*lorentzian
      implicit real*8(a-h,o-z)
      data tolmax/1.d10/
      pilrntz    = delta/(x*x + delta*delta)
      if (dabs(pilrntz) .gt. tolmax) stop 'fn. pilrntz: .gt. tolmax'
      return
      end
c---------------------------------------------------------------------
      subroutine piltzn (x,delta,n,
     o                   piltz)

c 92.04.19
c calculates pi*lorentzian n times

      implicit real*8 (a-h,o-z)
      dimension x(n),piltz(n)

      do       i = 1,n
        piltz(i)   = pilrntz (x(i),delta)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine pl(x,p,lmax)
      implicit real*8(a-h,o-z)
      dimension p(lmax+1)
      p(1) = 1.
      if(lmax.le.0) stop 'lmax l.t. 0'
      a2 = 0.
      p2 = a2
      p3 = 1.
      do 11 l = 1,lmax
        a1 = a2
c     a2 = dfloat(l)
        a2 = dble(l)
        p1 = p2
        p2 = p3
        p3 = ((a1 + a2)*x*p2 - a1*p1)/a2
        p(l+1) = p3
   11 continue
      return
      end
c--------------------------------------------------------------------
      subroutine plgaus (lmx,ng,
     o x,w,plg )

c generate P(x,l) at gaussian points x between -1,1

c lmx = maximum l
c ng  = no. gaussian points

c x   = gaussian points between -1,1
c w   = gaussian weights
c plg = P(x,l)

      implicit real*8(a-h,o-z)
      dimension x(ng),w(ng),plg(ng,0:lmx)

c local array
      parameter (lmax=20)
      dimension plt(0:lmax)

      if(lmx .gt. lmax) stop 'plgaus: increase lmax'
      call gauss   (ng,-1.d0,1.d0,x,w)
      do      ig = 1,ng
        call pl      (x(ig),plt(0),lmx)
        do       l = 0,lmx
          plg(ig,l)  = plt(l)
        end do
      end do

      return
      end
c-------------------------------------------------------------------
      subroutine ppp(nr,a,bpr,p1,p2,p3,result)
      implicit real*8(a-h,o-z)

c Integrate dr r*r*R1*R2*R3, with nr points including end points,
c where r = b(exp(x)-1) with x on equal mesh of size a.
c dr = (b+r)dx
c bpr = b+r/r
c p = rR, R=radial wavefunction
c With Simpson's rule, nr need not be odd.

      dimension bpr(nr),p1(nr),p2(nr),p3(nr)
      data third/0.3333333333333333d0/

      nr2        = nr/2
      nrt        = nr
      if(nr2*2 .eq. nr)nrt = nr - 1
      nr2        = nrt - 1

      even       = 0.d0
      odd        = 0.d0
      do 10    i = 2,nr2,2
        even       = even + bpr(i)*p1(i)*p2(i)*p3(i)
   10 continue
      do 20    i = 3,nr2+1,2
        odd        = odd  + bpr(i)*p1(i)*p2(i)*p3(i)
   20 continue
      result     = bpr(1)*p1(1)*p2(1)*p3(1)
     +           - bpr(nrt)*p1(nrt)*p2(nrt)*p3(nrt)
     +           + 4.d0*even + 2.d0*odd

      if(nrt .lt. nr)
     +result     = result
     +    + 0.25d0*(      -bpr(nr-2)*p1(nr-2)*p2(nr-2)*p3(nr-2)
     +               +8.d0*bpr(nr-1)*p1(nr-1)*p2(nr-1)*p3(nr-1)
     +               +5.d0*bpr(nr)  *p1(nr)  *p2(nr)  *p3(nr)  )

      result     = result*a*third
      return
      end
c-------------------------------------------------------------------
      subroutine rphi4(nr,a,r2,p1,p2,p3,p4,result)
      implicit real*8(a-h,o-z)

c integrate dr r*r*R1*R2*R3*R4 with r(i) = b*( exp((i-1)*a) - 1 )
c with Simpson's rule.
c inputs:
c  p1 = r*R1 etc.
c  r2 = (b+r)/r**2

      dimension r2(nr),p1(nr),p2(nr),p3(nr),p4(nr)
      data third/0.3333333333333333d0/

      nr2        = nr/2
      nrt        = nr
      if(nr2*2 .eq. nr)nrt = nr - 1
      nr2        = nrt - 1

      even       = 0.d0
      odd        = 0.d0
      do 10    i = 2,nr2,2
        even       = even + r2(i)  *p1(i)  *p2(i)  *p3(i)  *p4(i)
        odd        = odd  + r2(i+1)*p1(i+1)*p2(i+1)*p3(i+1)*p4(i+1)
   10 continue
      result     = 
     .           - r2(nrt)*p1(nrt)*p2(nrt)*p3(nrt)*p4(nrt)
     .           + 4.d0*even + 2.d0*odd

      if(nrt .lt. nr)
     .result     = result + 0.25d0*
     . (   -r2(nr-2)*p1(nr-2)*p2(nr-2)*p3(nr-2)*p4(nr-2)
     .+8.d0*r2(nr-1)*p1(nr-1)*p2(nr-1)*p3(nr-1)*p4(nr-1)
     .+5.d0*r2(nr)  *p1(nr)  *p2(nr)  *p3(nr)  *p4(nr)  )

      result     = result*a*third
      return
      end
c-------------------------------------------------------------------
      subroutine scudir(nsd,
     o                  wsd,xsd,ysd,zsd )

c generates special cubic directions and the weights
c I = < [dr^] f(r^) >
c   = S[i=1,nsd] wsd(i) f(i)
c where f(i) is the function calculated at the i'th direction

c nsd    = number of cubic special directions (62,114)

c wsd    = weights
c xsd,ysd,zsd, = x,y,z components of special cubic directions
c xsd^2 + ysd^2 + zsd^2 = 1

      implicit real*8(a-h,o-z)
      dimension xsd(nsd),ysd(nsd),zsd(nsd),wsd(nsd)
      dimension w(10),x(10),y(10),z(10),
     .          xst(48),yst(48),zst(48)
      save xst,yst,zst
      data pi/3.1415926535897932d0/

c-special cubic directions
      fpi           = 4.d0*pi
      if(nsd .eq. 62) goto 62
      if(nsd .eq. 114)goto 114
      stop 'scudir: wrong number of directions'

   62 w(1)          = .130612244897931d0*fpi/6.d0
      w(2)          = .128571428571554d0*fpi/8.d0
      w(3)          = .740816326530515d0*fpi/48.d0
      x(1)          = 1.d0
      y(1)          = 0.d0
      z(1)          = 0.d0
      sqrt3         = 1.d0/dsqrt(3.d0)
      x(2)          = sqrt3
      y(2)          = sqrt3
      z(2)          = sqrt3
      x(3)          = .846433804070399d0
      y(3)          = .497257813599068d0
      z(3)          = .190484860662438d0
      nshell        = 3
      if(nshell .eq. 3)goto 2222

  114 w(1)          = .076190476192774d0 *fpi/ 6.d0
      w(2)          = .137357478197258d0 *fpi/ 12.d0
      w(3)          = .344086737167612d0 *fpi/ 48.d0
      w(4)          = .442365308442356d0 *fpi/ 48.d0
      x(1)          = 1.d0
      y(1)          = 0.d0
      z(1)          = 0.d0
      x(2)          = .707106781186547d0
      y(2)          = x(2)
      z(2)          = 0.d0
      x(3)          = .733519276107007d0
      y(3)          = .570839829704020d0
      z(3)          = .368905625333822d0
      x(4)          = .909395474471327d0
      y(4)          = .385850474128732d0
      z(4)          = .155303839700451d0
      nshell        = 4

 2222 j2            = 0
      do 100      i = 1,nshell
        xt            = x(i)
        yt            = y(i)
        zt            = z(i)
        j1            = j2 + 1
        call cubstar(xt,yt,zt,
     o             xst,yst,zst,nstar)

        j2            = j1 + nstar - 1
        do 110      j = j1,j2
          xsd(j)        = xst(j-j1+1)
          ysd(j)        = yst(j-j1+1)
          zsd(j)        = zst(j-j1+1)
          wsd(j)        = w(i)
  110   continue
  100 continue
      if(j2 .ne. nsd) stop 'scudir: number of directions do not match'

ctemp
c     iftemp = ifile('TEMP')
c     do i = 1,nsd
c     write(iftemp,6000)xsd(i),ysd(i),zsd(i),wsd(i)
c     enddo
c6000 format(1x,4d13.6)
c     stop 'scudir'

      return
      end
c---------------------------------------------------------------------
      double precision function sumv (v,n)

c 92.02.14
c sum the elements of a vector

c v = vector with dimension n

      implicit real*8 (a-h,o-z)
      dimension v(n)

      sumv       = 0.d0
      do       i = 1,n
        sumv       = sumv + v(i)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine simpson(nr,a,br,f,result)
      implicit real*8(a-h,o-z)

c integrate dr f(r) with r(i) = b*( exp((i-1)*a) - 1 )
c br = b + r
c with Simpson's rule.

      dimension f(nr),br(nr)
      data third/0.3333333333333333d0/

      nr2        = nr/2
      nrt        = nr
      if(nr2*2 .eq. nr)nrt = nr - 1
      nr2        = nrt - 1

      even       = 0.d0
      odd        = 0.d0
      do 10    i = 2,nr2,2
        even       = even + br(i) * f(i)
        odd        = odd  + br(i+1) * f(i+1)
   10 continue
      result     = br(1)*f(1) - br(nrt)*f(nrt)
     +           + 4.d0*even + 2.d0*odd

      if(nrt .lt. nr)
     +result     = result
     +    + 0.25d0*(      -br(nr-2)*f(nr-2)
     +               +8.d0*br(nr-1)*f(nr-1)
     +               +5.d0*br(nr)  *f(nr)  )

      result     = result*a*third
      return
      end
c--------------------------------------------------------------------
      subroutine simpson2(nr,a,br,f,g,result)
      implicit real*8(a-h,o-z)

c integrate dr f(r)*g(r) with r(i) = b*( exp((i-1)*a) - 1 )
c and br = b+r
c with Simpson's rule.

      dimension f(nr),g(nr),br(nr)
      data third/0.3333333333333333d0/

      nr2        = nr/2
      nrt        = nr
      if(nr2*2 .eq. nr)nrt = nr - 1
      nr2        = nrt - 1

      even       = 0.d0
      odd        = 0.d0
      do 10    i = 2,nr2,2
        even       = even + br(i) * f(i) * g(i)
        odd        = odd  + br(i+1) * f(i+1) * g(i+1)
   10 continue
      result     = br(1)*f(1)*g(1) 
     +           - br(nrt)*f(nrt)*g(nrt)
     +           + 4.d0*even + 2.d0*odd

      if(nrt .lt. nr)
     +result     = result
     +    + 0.25d0*(      -br(nr-2)*f(nr-2)*g(nr-2)
     +               +8.d0*br(nr-1)*f(nr-1)*g(nr-1)
     +               +5.d0*br(nr)  *f(nr)  *g(nr)  )

      result     = result*a*third
      return
      end
c-------------------------------------------------------------------
      subroutine sl1 (l,s,ns,
     o osl1 )

c calculates 1/s^(l+1)

      implicit real*8(a-h,o-z)
      dimension s(ns)
      dimension osl1(ns)

      do       i = 2,ns
        osl1(i)    = 1.d0 / s(i)**(l+1)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine sngld (sv,
     d                  n,
     o                  dv)

c 92.02.07
c convert a single precision array sv to double dv

      implicit real*8 (a-h,o-z)
      real*4 sv(n)
      real*8 dv(n)

      do       i = 1,n
        dv(i)      = dble(sv(i))
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine specub (nrg,nsd,xmt,b,
     o                   rgaus,wr,
     o                   xsd,ysd,zsd,wsd,
     o                   xsp,ysp,zsp,wsp )

c generates special points for three dimensional integration
c e.g. I = < [d3r] f(r,r^) >
c        = < [r^] [r=0,R] r^2 f(r,r^) >
c        = < [r^] [x=0,X] (b+r)r^2 f(r,r^) >
c        = S[isd=1,nsd] wsd(isd)
c        * { S[ir=1,nrg] wr(ir) r(ir)^2 f(ir,isd) }
c        = S[i=1,nrg*nsd] wsp(i) r(i)^2 f(i)
c where f(i) is the function calculated at
c r(i) = [xsp(i),ysp(i),zsp(i)]
c note that the function should be multiplied by r^2

c nrg   = no. radial gaussian points
c nsd   = no. special cubic directions (62,114)
c xmt   = end point of equal mesh axis x in r = b(exp(x)-1)
c       = (nr-1)*a, where nr is the no. radial mesh points
c                   and a is the equal mesh size
c b     = see xmt
c xsd,ysd,zsd = x,y,z components of special cubic directions
c  xsd^2 + ysd^2 + zsd^2 = 1
c wsd   = weights of special cubic directions
c rgaus = radial gaussian points
c wr    = radial weights, can be used to integrate
c         < [r=0,R] g(r) > = < [x=0,X] (b+r) g(r) >
c                          = S[i=1,nr] wr(i) g(i)
c         g(i) is calculated at radial gaussian points rgaus(i)
c xsd,ysd,zsd = x,y,z components of special cubic directions
c xsd^2 + ysd^2 + zsd^2 = 1
c wsd   = weights of special cubic directions, generated by scudir.f

c xsp,ysp,zsp = x,y,z components of special points
c  xsp  = xsd*rgaus etc.
c wsp   = weights of special points
c       = wsd*wr
c the weight is ordered according to radial-angular

      implicit real*8(a-h,o-z)

      dimension rgaus(nrg),wr(nrg),
     o          xsd(nsd),ysd(nsd),zsd(nsd),wsd(nsd),
     o          xsp(nrg*nsd),ysp(nrg*nsd),zsp(nrg*nsd),
     o          wsp(nrg*nsd)

c generate gaussian points along radial axis
      call gauss   (nrg,0.d0,xmt,
     o              rgaus,wr)
      do       i = 1,nrg
        rgaus(i)   = b * ( dexp(rgaus(i)) - 1.d0 )
        wr(i)      = (b+rgaus(i)) * wr(i)
      end do

c generate special cubic directions
      call scudir  (nsd, 
     o              xsd,ysd,zsd,wsd)

      ic         = 0
      do       i = 1,nsd
        do       j = 1,nrg
          ic         = ic + 1
          xsp(ic)    = xsd(i)*rgaus(j)
          ysp(ic)    = ysd(i)*rgaus(j)
          zsp(ic)    = zsd(i)*rgaus(j)
          wsp(ic)    = wsd(i)*wr(j)
        end do
      end do

      return
      end
c-------------------------------------------------------------------
      subroutine specub1(nrg,nsd,xmt,b,
     i                   xsd,ysd,zsd,wsd,
     o                   xgaus,rgaus,wr,
     o                   xsp,ysp,zsp,wsp )

c generates special points for three dimensional integration
c e.g. I = < [d3r] f(r,r^) >
c        = < [r^] [r=0,R] r^2 f(r,r^) >
c        = < [r^] [x=0,X] (b+r)r^2 f(r,r^) >
c        = S[isd=1,nsd] wsd(isd)
c        * { S[ir=1,nrg] wr(ir) r(ir)^2 f(ir,isd) }
c        = S[i=1,nrg*nsd] wsp(i) r(i)^2 f(i)
c where f(i) is the function calculated at
c r(i) = [xsp(i),ysp(i),zsp(i)]
c note that the function should be multiplied by r^2

c nrg   = no. radial gaussian points
c nsd   = no. special cubic directions (62,114)
c xmt   = end point of equal mesh axis x in r = b(exp(x)-1)
c       = (nr-1)*a, where nr is the no. radial mesh points
c                   and a is the equal mesh size
c b     = see xmt
c xsd,ysd,zsd = x,y,z components of special cubic directions
c  xsd^2 + ysd^2 + zsd^2 = 1
c wsd   = weights of special cubic directions, generated by scudir.f

c rgaus = radial gaussian points
c wr    = radial weights, can be used to integrate
c         < [r=0,R] g(r) > = < [x=0,X] (b+r) g(r) >
c                          = S[i=1,nr] wr(i) g(i)
c         g(i) is calculated at radial gaussian points rgaus(i)
c xsp,ysp,zsp = x,y,z components of special points
c  xsp  = xsd*rgaus etc.
c wsp   = weights of special points
c       = wsd*wr
c the weight is ordered according to radial-angular
c it already contains the (b+r) factor

      implicit real*8(a-h,o-z)

      dimension xsd(nsd),ysd(nsd),zsd(nsd),wsd(nsd)
      dimension xgaus(nrg),rgaus(nrg),wr(nrg),
     o          xsp(nrg*nsd),ysp(nrg*nsd),zsp(nrg*nsd),
     o          wsp(nrg*nsd)

c generate gaussian points along radial axis
      call gauss   (nrg,0.d0,xmt,
     o              xgaus,wr)
      do       i = 1,nrg
        rgaus(i)   = b * ( dexp(xgaus(i)) - 1.d0 )
        wr(i)      = (b+rgaus(i)) * wr(i)
      end do

c generate special points
      ic         = 0
      do       i = 1,nsd
        do       j = 1,nrg
          ic         = ic + 1
          xsp(ic)    = xsd(i)*rgaus(j)
          ysp(ic)    = ysd(i)*rgaus(j)
          zsp(ic)    = zsd(i)*rgaus(j)
          wsp(ic)    = wsd(i)*wr(j)
        end do
      end do

      return
      end
c--------------------------------------------------------------------
      subroutine specuba(xmt,b,
     i                   xsd,ysd,zsd,wsd,
     d                   nrg,nsd,nclass,
     o                   xgaus,rgaus,wr,
     o                   xsp,ysp,zsp,wsp )

c 92.01.21
c generates special points for three dimensional integration
c for all classes
c e.g. I = < [d3r] f(r,r^) >
c        = < [r^] [r=0,R] r^2 f(r,r^) >
c        = < [r^] [x=0,X] (b+r)r^2 f(r,r^) >
c        = S[isd=1,nsd] wsd(isd)
c        * { S[ir=1,nrg] wr(ir) r(ir)^2 f(ir,isd) }
c        = S[i=1,nrg*nsd] wsp(i) r(i)^2 f(i)
c where f(i) is the function calculated at
c r(i) = [xsp(i),ysp(i),zsp(i)]
c note that the function should be multiplied by r^2

c nrg   = no. radial gaussian points
c nsd   = no. special cubic directions (62,114)
c xmt   = end point of equal mesh axis x in r = b(exp(x)-1)
c       = (nr-1)*a, where nr is the no. radial mesh points
c                   and a is the equal mesh size
c b     = see xmt
c xsd,ysd,zsd = x,y,z components of special cubic directions
c  xsd^2 + ysd^2 + zsd^2 = 1
c wsd   = weights of special cubic directions, generated by scudir.f

c rgaus = radial gaussian points
c wr    = radial weights, can be used to integrate
c         < [r=0,R] g(r) > = < [x=0,X] (b+r) g(r) >
c                          = S[i=1,nr] wr(i) g(i)
c         g(i) is calculated at radial gaussian points rgaus(i)
c xsp,ysp,zsp = x,y,z components of special points
c  xsp  = xsd*rgaus etc.
c wsp   = weights of special points
c       = wsd*wr
c the weight is ordered according to radial-angular
c it already contains the (b+r) factor

      implicit real*8(a-h,o-z)

      dimension xsd(nsd),ysd(nsd),zsd(nsd),wsd(nsd),
     i          b(nclass),xmt(nclass)

      dimension xgaus(nrg,nclass),rgaus(nrg,nclass),wr(nrg,nclass),
     o          xsp(nrg*nsd,nclass),ysp(nrg*nsd,nclass),
     o          zsp(nrg*nsd,nclass),wsp(nrg*nsd,nclass)

      do      ic = 1,nclass
        call specub1 (nrg,nsd,xmt(ic),b(ic),
     i              xsd,ysd,zsd,wsd,
     o              xgaus(1,ic),rgaus(1,ic),wr(1,ic),
     o              xsp(1,ic),ysp(1,ic),zsp(1,ic),wsp(1,ic) )
      end do

      return
      end
c--------------------------------------------------------------------
      subroutine specwa (nrg,nsd,xmt,b,nclass,
     i                   wsd,
     o                   xgaus,rgaus,wr,
     o                   wsp )

c generates the weights of special points for three
c dimensional integration for all classes
c e.g. I = < [d3r] f(r,r^) >
c        = < [r^] [r=0,R] r^2 f(r,r^) >
c        = < [r^] [x=0,X] (b+r)r^2 f(r,r^) >
c        = S[isd=1,nsd] wsd(isd)
c        * { S[ir=1,nrg] wr(ir) r(ir)^2 f(ir,isd) }
c        = S[i=1,nrg*nsd] wsp(i) r(i)^2 f(i)
c where f(i) is the function calculated at
c r(i) = [xsp(i),ysp(i),zsp(i)]
c note that the function should be multiplied by r^2

c nrg   = no. radial gaussian points
c nsd   = no. special cubic directions (62,114)
c xmt   = end point of equal mesh axis x in r = b(exp(x)-1)
c       = (nr-1)*a, where nr is the no. radial mesh points
c                   and a is the equal mesh size
c b     = see xmt
c wsd   = weights of special cubic directions, generated by scudir.f

c rgaus = radial gaussian points
c wr    = radial weights, can be used to integrate
c         < [r=0,R] g(r) > = < [x=0,X] (b+r) g(r) >
c                          = S[i=1,nr] wr(i) g(i)
c         g(i) is calculated at radial gaussian points rgaus(i)
c wsp   = weights of special points
c       = wsd*wr
c the weight is ordered according to radial-angular
c it already contains the (b+r) factor

      implicit real*8(a-h,o-z)

      dimension b(nclass),xmt(nclass),wsd(nsd)
      dimension xgaus(nrg,nclass),rgaus(nrg,nclass),wr(nrg,nclass),
     o          wsp(nrg*nsd,nclass)

      do      ic = 1,nclass

c generate gaussian points along radial axis
        call gauss   (nrg,0.d0,xmt(ic),
     o              xgaus(1,ic),wr(1,ic))
        do       i = 1,nrg
          rgaus(i,ic)= b(ic) * ( dexp(xgaus(i,ic)) - 1.d0 )
          wr(i,ic)   = (b(ic)+rgaus(i,ic)) * wr(i,ic)
        end do

c generate weights for special points
        ind        = 0
        do       i = 1,nsd
          do       j = 1,nrg
            ind        = ind + 1
            wsp(ind,ic)= wsd(i)*wr(j,ic)
          end do
        end do

      end do

      return
      end
c--------------------------------------------------------------------
      subroutine sphbes(x,lmax,a)
c
c this subroutine computes the regular spherical bessel functions j(l,x)
c for l = 0,..,lmax and returns the answer in a(i), i=1,..,lmax+1
c
      implicit real*8(a-h,o-z)
      dimension a(lmax+1)
      ax = dabs(x)
      lp = lmax + 1
      do 11 l = 1,lp
        a(l) = 0.d0
   11 continue
      a(1) = 1.d0
      if(ax.eq.0.d0) go to 99
      rx = 1.d0/x
      if(lmax.gt.0) go to 12
      a(1) = dsin(x)*rx
      go to 99
c  12 fl = dfloat(lp)
   12 fl = dble(lp)
      if(ax.gt.fl) go to 17
      f2 = 1.d0
      g2 = 1.d0
      do 13 l = 1,lmax
        g2 = g2 + 2.d0
        f2 = f2*x/g2
   13 continue
      f1 = f2*g2/x
      f3 = 1.d0
      g3 = 0.d0
      s1 = f1
      s2 = f2
      x2 = -x*x*0.5d0
      do 14 j = 1,100
        g1 = g2
        g2 = g2 + 2.d0
        g3 = g3 + 1.d0
        f4 = x2/(g2*g3)
        f1 = f1*x2/(g1*g3)
        f2 = f2*f4
        f3 = f3*f4
        s1 = s1 + f1
        s2 = s2 + f2
        if(dabs(f3).lt.1.d-17) go to 15
   14 continue
      stop 921
   15 a(lp) = s2
      a(lmax) = s1
      if(lmax.lt.2) go to 99
      lp = lmax - 1
c     fl = dfloat(2*lmax+1)
      fl = dble(2*lmax+1)
      m = lmax
      do 16 l = 1,lp
        m = m - 1
        fl = fl - 2.d0
        f3 = s1
        s1 = fl*rx*s1 - s2
        s2 = f3
        a(m) = s1
   16 continue
      go to 99
   17 f2 = rx*dcos(x)
      f3 = rx*dsin(x)
      fl = 1.d0
      a(1) = f3
      do 18 l = 1,lmax
        f1 = f2
        f2 = f3
        f3 = fl*rx*f2 - f1
        fl = fl + 2.d0
        a(l+1) = f3
   18 continue
   99 return
      end
c-----------------------------------------------------------------
      subroutine sylm(r,yl,lmx,r2s)
C- Generate unnormalized spherical harmonic polynomials
C ----------------------------------------------------------------
Ci Inputs
Ci   R,LMX, R should be normalised to unity
Co Outputs
Co   YL:  Cubic harmonics
Co   R2S: Square of vector R.
Cr Remarks
Cr   POLAR AXIS ALONG 001 AXIS. (ADAPTED FROM ASW PROGS)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer lmx
      double precision r2s
      double precision r(3),yl(*)
c     double precision r(3),yl((lmx+1)*(lmx+1))
cslut
C Local parameters
      integer i,l,lav,lavml,lavmm,lavpl,lavpm,lm1,lmm,lp1,m,mp1,n,nt
      double precision r2,st,x,y,z,z2
      integer(4),parameter:: lmxx=50
      double precision c(lmxx+1),s(lmxx+1),p(lmxx+1,lmxx+1)
      equivalence (x,c(2)),(y,s(2)),(z,p(2,1))
      data c(1),s(1),p(1,1),p(2,2) /1.d0,0.d0,1.d0,1.d0/
      n = (lmx+1)**2
      yl(1) = 1.d0
      x = r(1)
      y = r(2)
      z = r(3)
  500 st = x*x + y*y
      z2 = z*z
      r2 = st+z2
      r2s = r2
      if (n .lt. 2) return
      if (r2 .gt. 1.d-10) goto 1
      do  6  i = 2, n
        yl(i) = 0.d0
    6 continue
      return
    1 yl(2) = y
      yl(3) = z
      yl(4) = x
      nt = 1
      do  5  l = 2, lmx
        lp1 = l+1
        lm1 = l-1
        lav = l*lp1 + 1
        p(lp1,1) = ((l+lm1)*z*p(l,1) - lm1*r2*p(lm1,1)) / l
        yl(lav) = p(lp1,1)
        nt = nt+2
        p(lp1,lp1) = p(l,l)*nt
        c(lp1) = x*c(l) - y*s(l)
        s(lp1) = x*s(l) + y*c(l)
        lavpl = lav+l
        yl(lavpl) = p(lp1,lp1)*c(lp1)
        lavml = lav-l
        yl(lavml) = p(lp1,lp1)*s(lp1)
        if (st .gt. z2) goto 3
        do  2  lmm = 1, lm1
          m = l-lmm
          lavpm = lav+m
          lavmm = lav-m
          mp1 = m+1
          p(lp1,mp1) = (r2*(l+m)*p(l,mp1) - st*p(lp1,mp1+1)) / (z*(l-m))
          yl(lavpm) = p(lp1,mp1)*c(mp1)
          yl(lavmm) = p(lp1,mp1)*s(mp1)
    2   continue
        goto 5
    3   do  4  m = 1, lm1
          mp1 = m+1
          lavpm = lav+m
          lavmm = lav-m
          p(lp1,mp1) = ((lm1+m)*r2*p(l,m) - (lp1-m)*z*p(lp1,m))/st
          yl(lavpm) = p(lp1,mp1)*c(mp1)
          yl(lavmm) = p(lp1,mp1)*s(mp1)
    4   continue
    5 continue
      end
c-----------------------------------------------------------------
      subroutine sylmnc(c,lmx)
C- Normalization constants for the spherical harmonics
C ----------------------------------------------------------------
Ci Inputs
Ci   lmx
Co Outputs
Co   C
Cr Remarks
Cr   use together with sylm (from ASW package)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer lmx
      double precision c(1)
C Local parameters
      integer i,i1,i2,l,lav,lp1,m,n1,n2,n3
      double precision fn2,fpi,tlp1,tpi,y0

      tpi = 8.d0*datan(1.d0)
      fpi = 2.d0*tpi
      y0 = 1.d0/dsqrt(fpi)
      c(1) = y0
      do  2  l = 1, lmx
        lp1 = l+1
        tlp1 = l+lp1
        lav = l*lp1 + 1
        c(lav) = dsqrt(tlp1/fpi)
      do  2  m = 1, l
        n2 = lp1-m
        n1 = n2+1
        n3 = l+m
        fn2 = n2
        do  1  i = n1, n3
          fn2 = fn2*i
    1   continue
        i1 = lav+m
        i2 = lav-m
        c(i1) = dsqrt(tlp1/(fn2*tpi))
        c(i2) = c(i1)
    2 continue
      end
c------------------------------------------------------------------
      subroutine vaddv 
     i (v1,v2,n,
     o vadd )

c add two vectors vadd = v1 + v2

      implicit real*8(a-h,o-z)
      dimension v1(n),v2(n)
      dimension vadd(n)

      do       i = 1,n
        vadd(i)    = v1(i) + v2(i)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine vaddvs (v1,v2,
     d                   n,
     o vadd )

c 92.02.07
c add two vectors vadd = v1 + v2, single precision

      implicit real*8 (a-h,o-z)
      real*4    v1(n),v2(n)
      real*4    vadd(n)

      do       i = 1,n
        vadd(i)    = v1(i) + v2(i)
      end do

      return
      end
c---------------------------------------------------------------------
      double precision function vdot3(f,g,wsp,nsp)

c calculates
c S[i=1,nsp] wsp(i) f(i) g(i)

      implicit real*8(a-h,o-z)
      dimension f(nsp),g(nsp),wsp(nsp)

      vdot3      = 0.d0
      do       i = 1,nsp
        vdot3      = vdot3+ wsp(i)*f(i)*g(i)
      end do

      return
      end
c---------------------------------------------------------------------
c---------------------------------------------------------------------
ctakao
      double precision function vdv_ (v1,v2,n)
c dot product v1.v2
c n = no. elements
      implicit real*8(a-h,o-z)
      dimension v1(n),v2(n)
      vdv_       = 0.d0
      do       i = 1,n
        vdv_       = vdv_ + v1(i)*v2(i)
      end do
      return
      end
c
      double precision function vdv (v1,v2,n)
c dot product v1.v2
c n = no. elements
      implicit real*8(a-h,o-z)
      dimension v1(n),v2(n)
      vdv=sum(v1(1:n)*v2(1:n))
      return
      end

c---------------------------------------------------------------------
      subroutine vdvc (rv1,cv1,rv2,cv2,
     d                 n,
     o                 r,c)

c 92.01.20
c dot product of complex vectors v1.v2
c rv1,cv1 = real and imaginary part of v1
c rv2,cv2 = real and imaginary part of v2
c n = size of vectors

c r,c = real and imaginary part of dot product

      implicit real*8 (a-h,o-z)
      dimension rv1(n),cv1(n),rv2(n),cv2(n)

      r          = 0.d0
      c          = 0.d0
      do       i = 1,n
        r          = r + rv1(i)*rv2(i) - cv1(i)*cv2(i)
        c          = c + rv1(i)*cv2(i) + cv1(i)*rv2(i)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine vdvc1(rv1,cv1,rv2,cv2,
     d                 n,
     o                 r,c)

c 92.01.20
c dot product of complex vectors cc[v1].v2
c rv1,cv1 = real and imaginary part of v1
c rv2,cv2 = real and imaginary part of v2
c n = size of vectors

c r,c = real and imaginary part of dot product

      implicit real*8 (a-h,o-z)
      dimension rv1(n),cv1(n),rv2(n),cv2(n)

      r          = 0.d0
      c          = 0.d0
      do       i = 1,n
        r          = r + rv1(i)*rv2(i) + cv1(i)*cv2(i)
        c          = c + rv1(i)*cv2(i) - cv1(i)*rv2(i)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine veqvs (v,n,w)

c 92.04.03
c equate v to w

      implicit real*8 (a-h,o-z)
      real*4 v(n),w(n)

      do       i = 1,n
        w(i)       = v(i)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine vmc (rv,cv,ra,ca,lda,n,lva,
     o rva,cva )

c 91.12.01
c multiply a complex vector and a complex matrix
c va = v a
c va(j) =  S[i =1,n] v(i) a(i,j)

c rv,cv = real and imaginary parts of v
c ra,ca = real and imaginary parts of a
c lda   = leading dimension of a
c n     = dimension of multiplication
c lva   = dimension of resulting vector

c rva,cva = real and imaginary parts of resulting vector

      implicit double precision (a-h,o-z)
      dimension rv(n),cv(n),ra(lda,lva),ca(lda,lva)
      dimension rva(lva),cva(lva)

      do       j = 1,lva
        rsum       = 0.d0
        csum       = 0.d0
        do       i = 1,n
          rsum       = rsum + rv(i)*ra(i,j) - cv(i)*ca(i,j)
          csum       = csum + rv(i)*ca(i,j) + cv(i)*ra(i,j)
c     rsum       = rsum + rv(i)*ra(i,j) + cv(i)*ca(i,j)
c     csum       = csum + rv(i)*ca(i,j) - cv(i)*ra(i,j)
        end do
        rva(j)     = rsum
        cva(j)     = csum
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine vminv 
     i (v1,v2,n,
     o vmin )

c substract two vectors vmin = v1 - v2

      implicit real*8(a-h,o-z)
      dimension v1(n),v2(n)
      dimension vmin(n)

      do       i = 1,n
        vmin(i)    = v1(i) - v2(i)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine vminvs
     i (v1,v2,n,
     o vmin )

c single precision subtract a vector from another vmin = v1 - v2

      implicit real*8(a-h,o-z)
      real*4    v1(n),v2(n)
      real*4    vmin(n)

      do       i = 1,n
        vmin(i)    = v1(i) - v2(i)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine vmvc  (rv1,cv1,rm,cm,rv2,cv2,
     d                  n1,ldim,n2,
     o                  r,c )

c 92.03.02
c calculates S[i=1,n1] S[j=1,n2] v1(i) m(i,j) cc[v2(j)]
c v1,v2,a are complex

c rv1,cv1 = real and imaginary part of v1
c rm,cm   = real and imaginary part of m
c rv2,cv2 = real and imaginary part of v2
c n1,n2   = dimensions of v1 and v2
c ldim    = leading dimension of m

      implicit real*8 (a-h,o-z)
      dimension rv1(n1),cv1(n1),rm(ldim,n2),cm(ldim,n2),rv2(n2),cv2(n2)

      if (n1 .gt. ldim) stop 'vmvc: n1 .gt. ldim'
      r          = 0.d0
      c          = 0.d0
      do       j = 1,n2
        call vdvc    (rv1,cv1,rm(1,j),cm(1,j),
     d              n1,
     o              r1,c1)
        r          = r + r1*rv2(j) + c1*cv2(j)
        c          = c + c1*rv2(j) - r1*cv2(j)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine vmvc1 (rv1,cv1,rm,cm,rv2,cv2,
     d                  n1,ldim,n2,
     o                  r,c )

c 92.03.02
c calculates S[i=1,n1] S[j=1,n2] cc[v1(i)] m(i,j) v2(j)
c v1,v2,a are complex

c rv1,cv1 = real and imaginary part of v1
c rm,cm   = real and imaginary part of m
c rv2,cv2 = real and imaginary part of v2
c n1,n2   = dimensions of v1 and v2
c ldim    = leading dimension of m

      implicit real*8 (a-h,o-z)
      dimension rv1(n1),cv1(n1),rm(ldim,n2),cm(ldim,n2),rv2(n2),cv2(n2)

      if (n1 .gt. ldim) stop 'vmvc1: n1 .gt. ldim'
      r          = 0.d0
      c          = 0.d0
      do       j = 1,n2
        call vdvc1   (rv1,cv1,rm(1,j),cm(1,j),
     d              n1,
     o              r1,c1)
        r          = r + r1*rv2(j) - c1*cv2(j)
        c          = c + r1*cv2(j) + c1*rv2(j)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine vv (v1,v2,n,
     o v )

c forms v(i) = v1(i) v2(i)

c v1,v2 = vectors of size n
c v     = v1 v2

      implicit real*8(a-h,o-z)
      dimension v1(n),v2(n)
      dimension v(n)

      do       i = 1,n
        v(i)       = v1(i)*v2(i)
      end do

      return
      end
c-------------------------------------------------------------------
      subroutine vva (v1,v2,a,n,
     o v )

c forms v(i) = v(i) + a v1(i) v2(i)

c v1,v2 = vectors of size n

      implicit real*8(a-h,o-z)
      dimension v1(n),v2(n)
      dimension v(n)

      do       i = 1,n
        v(i)       = v(i) + a*v1(i)*v2(i)
      end do

      return
      end
c---------------------------------------------------------------------
      subroutine vvc (rv1,cv1,rv2,cv2,n,
     o rv,cv )

c 92.jan.07
c forms v(i) = cc[v1(i)] v2(i), for complex vectors

c v1,v2 = vectors of size n
c v     = v1 v2

      implicit real*8(a-h,o-z)
      dimension rv1(n),cv1(n),rv2(n),cv2(n)
      dimension rv(n),cv(n)

      do       i = 1,n
        rv(i)      = rv1(i)*rv2(i) + cv1(i)*cv2(i)
        cv(i)      = rv1(i)*cv2(i) - cv1(i)*rv2(i)
      end do

      return
      end
c-------------------------------------------------------------------
      subroutine vvc1 (rv1,cv1,rv2,cv2,n,
     o rv,cv )

c 92.jan.07
c forms v(i) = v1(i) v2(i), for complex vectors

c v1,v2 = vectors of size n
c v     = v1 v2

      implicit real*8(a-h,o-z)
      dimension rv1(n),cv1(n),rv2(n),cv2(n)
      dimension rv(n),cv(n)

      do       i = 1,n
        rv(i)      = rv1(i)*rv2(i) - cv1(i)*cv2(i)
        cv(i)      = rv1(i)*cv2(i) + cv1(i)*rv2(i)
      end do

      return
      end
c-------------------------------------------------------------------
      subroutine ylm3(y1,y2,y3,wsd,nsd, 
     o result)

c angular integration of y1*y2*y3 = result
c y1 etc = ylm, or any function of theta and phi
c nsd    = number of directions
c wsd    = weight

      implicit real*8(a-h,o-z)
      dimension y1(nsd),y2(nsd),y3(nsd),wsd(nsd)

      result     = 0.d0
      do 10    i = 1,nsd
        result     = result + wsd(i)*y1(i)*y2(i)*y3(i)
   10 continue

      return
      end
c-------------------------------------------------------------------
      subroutine ylm4(y1,y2,y3,y4,wtp,ndir,
     o result)

c angular integration of y1*y2*y3*y4 = result
c y1 etc = ylm, or any function of theta and phi
c ndir   = number of directions
c wtp    = weight

      implicit real*8(a-h,o-z)
      dimension y1(ndir),y2(ndir),y3(ndir),y4(ndir),wtp(ndir)

      result     = 0.d0
      do 10    i = 1,ndir
        result     = result + wtp(i)*y1(i)*y2(i)*y3(i)*y4(i)
   10 continue

      return
      end
c--------------------------------------------------------------------
      subroutine ylmsc(xgaus,ygaus,zgaus,nsd,
     w                 c,yl,lmax,
     o                 ylm)

c generates real ylm at special cubic directions

c nsd    = number of special directions
c xgaus,ygaus,zgaus = cartesian coordinates of special directions
c                     on unit sphere
c lmax   = max. l

c ylm    = real ylm(direction,l*l + l+1+m)

      implicit real*8(a-h,o-z)
      dimension xgaus(nsd),ygaus(nsd),zgaus(nsd),
     w          c((lmax+1)*(lmax+1)),yl((lmax+1)*(lmax+1))
      dimension ylm(nsd,(lmax+1)*(lmax+1))
      dimension rt(3)

      lmx2       = (lmax+1)*(lmax+1)
      call sylmnc(c,lmax)
ctemp
c     iftemp = ifile('TEMP')
c     do   i = 1,nsd
c     write (iftemp,6000)xgaus(i),ygaus(i),zgaus(i)
c     enddo
c     write(iftemp,*)' '
c     write(iftemp,6000)(c(i),i=1,(lmax+1)**2)
c6000 format(1x,5d13.6)
c     stop 'ylmsc'

      do       i = 1,nsd
        rt(1)      = xgaus(i)
        rt(2)      = ygaus(i)
        rt(3)      = zgaus(i)
        call sylm(rt,yl,lmax,r2s)
        do       j = 1,lmx2
          ylm(i,j)   = c(j)*yl(j)
        end do
      end do

      return
      end
c--------------------------------------------------------------------
      subroutine yyy (ylm,wsd,nsd,nl,
     o y3,ny3)

c angular integrals
c y3(i1,i2,i3) = integral ylm(i1)*ylm(i2)*ylm(i3)
c i1 = (l1,m1), l1=<2*nl-1
c i2 = (l2,n2), l2=<nl-1
c i3 = (l3,n3), l3=<nl-1

c ylm      = real ylm at special directions, see ylmsc.f
c wsd      = weight of special directions, see scudir.f
c nsd      = number of special directions

c ny3      = no. integrals

      implicit real*8(a-h,o-z)
      dimension ylm(nsd,nl*nl),wsd(nsd)
      dimension y3((2*nl-1)**2,nl*nl,nl*nl)

c angular integrals
      ll         = nl*nl
      ll2        = (2*nl-1)**2
      i          = 0
      do 30   i1 = 1,ll2
      do 30   i2 = 1,ll
      do 30   i3 = 1,ll
        call ylm3(ylm(1,i1),ylm(1,i2),ylm(1,i3),
     .          wsd,nsd,resa)
        i          = i + 1
        y3(i1,i2,i3) = resa
   30 continue
      ny3        = i
      if(ny3 .ne. nl*nl*nl*nl*(2*nl-1)**2 )
     . stop 'yyy: wrong no. integrals'

      return
      end
c--------------------------------------------------------------------
      subroutine yyyy (ylm,wsd,nsd,nl,
     i                 nyyyy,
     o yy,iyy)

c angular integrals
c yy(i1,i2,i3,i4) = integral ylm(i1)*ylm(i2)*ylm(i3)*ylm(i4)
c    i1 .le. i2 .le. i3 .le. i4
c iyy(i1,i2,i3,i4) = index, from 1 to nyyyy
c nyyy     = number of integrals, from nofyyyy.f
c ylm      = real ylm at special directions, see ylmsc.f
c wsd      = weight of special directions, see scudir.f
c nsd      = number of special directions

      implicit real*8(a-h,o-z)
      dimension ylm(nsd,nl*nl),wsd(nsd)
      dimension yy(nyyyy),iyy(nl*nl,nl*nl,nl*nl,nl*nl)

ctemp
c     iftemp = ifile('TEMP')
c     write(iftemp,6000)(wsd(i),i=1,nsd)
c     write(iftemp,*)' '
c     do j = 1,nl*nl
c     write(iftemp,6000)(ylm(i,j),i=1,nsd)
c     write(iftemp,*)' '
c     enddo
c6000 format(1x,5d13.6)
c     stop 'yyyy'

c angular integrals
      ll         = nl*nl
      i          = 0
      do 30   i1 = 1,ll
      do 30   i2 = i1,ll
      do 30   i3 = i2,ll
      do 30   i4 = i3,ll
        call ylm4(ylm(1,i1),ylm(1,i2),ylm(1,i3),ylm(1,i4),
     .          wsd,nsd,resa)
        i          = i + 1
        yy(i)      = resa
        iyy(i1,i2,i3,i4) = i
   30 continue
      if(i .ne. nyyyy)stop 'yyyy: wrong no. of integrals'

      return
      end
c-------------------------------------------------------------------
      subroutine y1234 (y1,y2,y3,y4,wsd,
     d                  nsd,nl,
     o yyyy)

c 92.02.19
c angular integrals
c yyyy(i1,i2,i3,i4) = integral ylm1(i1)*ylm2(i2)*ylm3(i3)*ylm4(i4)

c y1,y2,y3,y4 = ylm1,ylm2,ylm3,ylm4 at special directions
c wsd      = weight of special directions, see scudir.f
c nsd      = number of special directions

c yyyy     = s.o.
      implicit real*8(a-h,o-z)
      dimension y1(nsd,nl*nl),y2(nsd,nl*nl),
     i          y3(nsd,nl*nl),y4(nsd,nl*nl),wsd(nsd)
      dimension yyyy(nl*nl,nl*nl,nl*nl,nl*nl)

c angular integrals
      ll         = nl*nl
      do 30   i1 = 1,ll
      do 30   i2 = 1,ll
      do 30   i3 = 1,ll
      do 30   i4 = 1,ll
        call ylm4(y1(1,i1),y2(1,i2),y3(1,i3),y4(1,i4),
     .          wsd,nsd,resa)
        yyyy(i1,i2,i3,i4) = resa
   30 continue

      return
      end
c-------------------------------------------------------------------
