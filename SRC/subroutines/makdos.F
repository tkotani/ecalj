      subroutine makdos(nqp,nband,nbmx,nsp,wgts,evl,n,w,tol,emin,emax,
     .ndos,dos)
C- Make density of states from bands
C-----------------------------------------------------------------------
Ci  Input
Ci   nqp   :number of q-points
Ci   nband :number of bands
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   wgts  :band weights
Ci   evl   :band eigenvalues
Ci   n,w   :Methfessel-Paxton order and broadening parameters
Ci   tol   :(tol>0) allowed error in DOS due to truncating the gaussian,
Ci         :        to a finite energy range (number of bins)
Ci         :(tol<0) dimensionless energy window specifying truncation
Ci         :        of gaussian.  Energy window for which gaussian is
Ci         :        taken to be nonzero is set to -tol*w
Ci   emin, emax, ndos: energy range and number of energy mesh points
Ci   nbmx  :leading dimension of evl
Co  Ouput
Co    dos: density of states
C-----------------------------------------------------------------------
C     implicit none
      integer nqp,nband,nbmx,nsp,n,ndos
      double precision wgts(nqp),evl(nbmx,nsp,nqp),dos(0:ndos-1,nsp),
     .w,emin,emax,tol,wt,emesh
      integer i,isp,iband,iq,meshpt,mesh1,mesh2,mrange,iprint
      double precision e,x,range,test,step,d,s,xx
      external delstp

      call dpzero(dos,nsp*ndos)
      step = (emax - emin) / (ndos - 1)
      if ( tol .gt. 0d0 ) then
        do  2  i = 0, ndos-1
          x = i * step / w
          call delstp(0,x,test,s,xx)
          if ( test .lt. tol ) then
            mrange = i + 1
            goto 3
          endif
    2   continue
        if (iprint() .gt. 30) print *,'makdos (warning) : tol too small'
    3   continue
        range = 2 * mrange * step
        test = tol
      else
        range = -tol * w
        mrange = range / ( 2 * step )
        call delstp(0,-tol/2,test,s,xx)
      endif
      if (iprint() .gt. 30) write (*,100) range/w,2*mrange,test
      do  7  iq = 1, nqp
        wt = abs(wgts(iq)) / nsp
        do  6  iband = 1, nband
        do  6  isp = 1, nsp
          e = evl(iband,isp,iq)
          meshpt = (e - emin) / step
          mesh1 = meshpt - mrange
          mesh2 = meshpt + mrange
          if (mesh2 .ge. ndos) mesh2 = ndos-1
          if (mesh1 .lt. 0) mesh1 = 0
          do  4  meshpt = mesh1, mesh2
            emesh = emin + meshpt * step
            x = (emesh - e) / w
            call delstp(n,x,d,s,xx)
            dos(meshpt,isp) = dos(meshpt,isp) + wt * d / w
    4     continue
    6   continue
    7 continue
  100 format(/1x,'MAKDOS :  range of gaussians is ',f5.2,
     .'W (',i4,' bins).'
     ./11x,'Error estimate in DOS : ',1pe8.2,' per state.')
      end

