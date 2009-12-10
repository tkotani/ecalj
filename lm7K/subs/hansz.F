      subroutine hansz(rsml,nrx,lmax,nxi,lxi,exi,nr,x,y,z,shan,sbes)
C- solid hankels for a vector of points, any energy.
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsml  : vector of l-dependent smoothing radii of Hankels
Ci         : EITHER must be specified for lmin..lmax
Ci         : OR     rsmh(0) = const < 0. Implies rsmh(l) = -const for all l
Ci   nrx   :leading dimension of yl, xi, and h
Ci   nxi   :number of trial energies
Ci exi,lxi :trial energies and corresponding angular momenta 
Ci   lmax  :maximum angular momentum, must be >= {lxi(ie),ie=1,nxi} 
Ci   nr    :actual number of points, n must be =< nrx
Ci  x,y,z  :Cartesian coordinates of points at which solid functions are 
Ci         :to be evaluated 
Co Outputs
Co   shan  :solid smoothed Hankel (e<0) or Neumann (e>0) functions, H=xi*Ylm
Ci   sbes  :if e>0, solid unsmoothed Bessel functions * e^{l+1/2) 
Ci         :if e<0, not referenced 
Cl Local variables
Ci   r2    :x^2 + y^2 + z^2 
Ci   yl    :Ylm(i,ilm), the (real) spherical harmonic polynomials at (x,y,z)
Co   xi    :e>0 radial part of smoothed Neumann function / r^l
Co         :e<0 radial part of smoothed Hankel function / r^l
Co   fi    :e>0 radial part of unsmoothed Bessel function * e^{l+1/2}/r^l
Cu Updates
Cu   04 May 07 (S. Lozovoi) First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmax,nrx,nxi,nr
      integer rsml(0:lmax),lxi(nxi)
      double precision x(nrx),y(nrx),z(nrx)
      double precision exi(nxi)
      double precision shan(nrx,(lmax+1)**2,nxi)
      double precision sbes(nrx,(lmax+1)**2,nxi)
C ... Local parameters
      integer il,im,ilm,ir,job,ie,lmx
      logical lpos
      double precision e
c ... work auto-arrays
      double precision r2(nrx)
      double precision xi(nrx,0:lmax),fi(nrx,0:lmax)
      double precision yl(nrx,(lmax+1)**2)
      

      if (nrx .lt. nr) 
     . call rxi('hansz: increase nrx, needed at least ',nr)

c Spherical harmonic polynomials
      call ropyln(nr,x,y,z,lmax,nrx,yl,r2)

      job = 00
c Begin the cycle over energies
      do  ie = 1, nxi  
        lmx = lxi(ie) 
        if (lmx .gt. lmax) 
     . call rxi('hansz: increase lmax, needed at least ',lmx)
        e = exi(ie) 
        lpos = (e .gt. 0d0) 

c Make radial parts of smooth Hankels
        call hansrz(rsml,0,lmx,e,r2,nrx,nr,job,xi,fi)

c Make smooth solid Hankels
        ilm = 0
        do il = 0, lmx 
           do im = -il, il
             ilm = ilm+1
             do ir = 1, nr
               shan(ir,ilm,ie) = xi(ir,il)*yl(ir,ilm)
               if(lpos) 
     .         sbes(ir,ilm,ie) = fi(ir,il)*yl(ir,ilm)
             enddo
           enddo
        enddo

c End the cycle over energies
      enddo

      end

