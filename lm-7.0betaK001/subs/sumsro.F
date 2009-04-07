      subroutine sumsro(rp,np,ips,a,b,rho,nttab,iax,rpos,rosum)
C- Add a r.s. superposition of spherical densities at a set of points
C ----------------------------------------------------------------------
Ci Inputs
Ci   rp    :set of np points
Ci   np    :number points
Ci   ips   :species table: density for site ib site ib found in ips(ib)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci         :a(i) is coefficient for species i
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci         :b(i) is coefficient for species i
Ci   rho   :spherical density for each species 
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   iax   :neighbor table containing pair information connecting
Ci         :site i to its neighbors
Ci   rpos  :positions of neighbors.  
Ci         :rp and rpos must be defined with respect to same origin.
Co Outputs
Co   rosum :sum of densities at each of np points
Cl Local variables
Cu Updates
Cu   21 Apr 02 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer np,ips(1),nrmx,niax
      parameter (nrmx=1501, niax=10)
      integer nttab,iax(niax,nttab)
      double precision rp(3,np),rpos(3,nttab),rho(nrmx,*),rosum(np)
      double precision a(*),b(*)
C ... Local parameters
      integer k,jb,js,ip,jx
      double precision dx(3),d,xx,y,e
      double precision rp3p,rp5p,rpp3p,rpp32,rpp5p,rppp
C     double precision rofii(nrmx),dy

C     x axis for polynomial interpolation.  Use log mesh
C     Needed only for call to polint
C      do  k = 1, nrmx
C        rofii(k) = k
C      enddo

C ... Sum the potential from nttab neighbors
      do  ip = 1, np
        rosum(ip)= 0
        do  k = 1, nttab
          dx(1) = rp(1,ip)-rpos(1,k)
          dx(2) = rp(2,ip)-rpos(2,k)
          dx(3) = rp(3,ip)-rpos(3,k)
          jb = iax(2,k)
          js = ips(jb)
          d = dsqrt(dx(1)**2+dx(2)**2+dx(3)**2)

          xx = 1 + dlog(d/b(js)+1)/a(js)
          if (int(xx) .le. nrmx-2) then
C           General polynomial interpolation.  Works, but slow.
C           jx = int(xx)-1
C           call polint(rofii,rho(1,js),nrmx,5,xx,0d0,0,jx,y,dy)

C           Interpolate in-line by 5-point formula
            jx = nint(xx)-2
C           first derivative, 3- and 5-point estimate
            rp3p = (rho(jx+3,js)-rho(jx+1,js))/2
            rp5p = (4*rp3p-(rho(jx+4,js)-rho(jx,js))/4)/(4-1)

C           Second derivative, nearest points and 2nd points
            rpp3p =  rho(jx+1,js) + rho(jx+3,js) - 2*rho(jx+2,js)
            rpp32 = (rho(jx+0,js) + rho(jx+4,js) - 2*rho(jx+2,js))/4
C           five-point estimate for second derivative
            rpp5p = (rpp3p*4 - rpp32*1) / (4 - 1)

C           Estimate for third derivative
            rppp = -6*(rp5p - rp3p)

C           Interpolate function value by polynomial
            e = xx-jx-2
            y = rho(jx+2,js)+e*rp5p+e*e/2*rpp5p+e**3/6*rppp

            rosum(ip) = rosum(ip) + y
          endif

C       Verbose debugging printout
C       print 333, k, iax(2,k), js, ip, dx, d, y, rosum(ip)
C 333   format(4i4,6f8.3)
      enddo
C     Debugging printout
C     print 334, ip, rosum(ip)
C  334 format(i4,f12.6)
      enddo

      end

