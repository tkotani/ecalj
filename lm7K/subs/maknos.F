      subroutine maknos(nqp,nband,nbmx,nsp,wgts,evl,n,w,tol,emin,emax,
     .ndos,dos)
C- Make density of states from bands
C-----------------------------------------------------------------------
Ci  Input
Ci    nqp : number of q-points; nband : number of bands;
Ci    nsp : 2 for spin polarised bands, 1 otherwise;
Ci    wgts, evl : weights and bands (eigenvalues);
Ci    nbmx : first dimension of evl ;
Ci    n   : n>0 Methfessel-Paxton polynomial order
Ci        : n<0 sampling done with Fermi-Dirac statistics
Ci    w   : n>0 gaussian width in Methfessel-Paxton integration (Ry)
Ci        : n<0 Temperature for Fermi distribution (Ry)
Ci    tol : allowed error in DOS due to truncating the gaussian,
Ci          if negative on entry, range is set to -tol*W
Ci    emin, emax, ndos; energy range and number of energy mesh points
Co  Ouput
Co    dos: integrated DOS
Cu Updates
Cu   2 Nov 1995 (JEK) returns spin-polarized integrated dos
C-----------------------------------------------------------------------
C     implicit none
      integer nqp,nband,nbmx,nsp,n,ndos
      double precision wgts(nqp),evl(nbmx,nsp,nqp),dos(0:ndos-1,nsp),
     .w,emin,emax,tol,wt,emesh
      integer i,isp,iband,iq,meshpt,mesh1,mesh2,mrange,iprint,i1mach
      double precision e,x,range,test,step,d,s,xx
      external delstp

      call dpzero(dos,nsp*ndos)
      step = (emax - emin) / (ndos - 1)
      if ( tol .gt. 0d0 ) then
        do  2  i = 0, ndos-1
          x = i * step / w
          call delstp(n,x,test,s,xx)
          if ( test .lt. tol ) then
            mrange = i + 1
            goto 3
          endif
    2   continue
        if (iprint() .gt. 30) print *,'maknos (warning) : tol too small'
    3   continue
        range = 2 * mrange * step
        test = tol
      else
        range = -tol * w
        mrange = range / ( 2 * step )
        call delstp(n,-tol/2,test,s,xx)
      endif
C     if (iprint() .ge. 40) write (*,100) range/w,2*mrange,test
      if (iprint() .ge. 40) call awrit3(' MAKNOS: range = %;2dW'//
     .' (%i bins)  DOS error estimate = %;2g per state',
     .' ',80,i1mach(2),range/w,2*mrange,test)
      do  7  iq = 1, nqp
        wt = abs(wgts(iq)) / nsp
        do  6  iband = 1, nband
          do  6  isp = 1, nsp
            e = evl(iband,isp,iq)
            meshpt = (e - emin) / step
            mesh1 = meshpt - mrange
            mesh2 = meshpt + mrange
            if (mesh2 .ge. ndos) mesh2 = ndos-1
            call rxx(mesh1 .lt. 0,'MAKNOS: emin too large')
            do  4  meshpt = mesh1, mesh2
              emesh = emin + meshpt * step
              x = (emesh - e) / w
              call delstp(n,x,d,s,xx)
              dos(meshpt,isp) = dos(meshpt,isp) + wt * (1d0 - s)
    4       continue
            do  5  meshpt = mesh2+1, ndos-1
              dos(meshpt,isp) = dos(meshpt,isp) + wt
    5       continue
    6   continue
    7 continue
C     call prmx('nos',dos,ndos,ndos,nsp)
C  100 format(1x,'MAKNOS:  range of gaussians is ',f5.2,'W (',i4,
C     .  ' bins).'/10x,'Error estimate in DOS : ',1pe8.2,' per state.')
      end

