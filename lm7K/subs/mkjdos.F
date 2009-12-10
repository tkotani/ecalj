      subroutine mkjdos(nqp,nbmx,nfilm,nempm,nsp,nfilo,nfiup,nemlo,
     .  nemup,wgts,evl,n,w,tol,emin,emax,esciss,jdos,optmt,efermi,
     .  ndos,dos)
C- Make joint density of states from bands, sampling
C-----------------------------------------------------------------------
Ci  Input
Ci    nqp  :number of q-points
Ci    nbmx :first dimension of evl
Ci    nfilm:second dimension of optmt 
Ci    nempm:third dimension of optmt 
Ci    nsp  :2 for spin polarised bands, 1 otherwise
Ci    nfilo:Loop over occupied bands nfilo, nfiup
Ci    nfiup:Loop over occupied bands nfilo, nfiup
Ci    nemlo:Loop over unoccupied bands nemlo, nemup
Ci    nemup:Loop over unoccupied bands nemlo, nemup
Ci    wgts :eigenvalue weights
Ci         :NB: only abs(wtkp) is used; see bzmesh.f
Ci    evl  :eigenvalues
Ci    N, W : Methfessel-Paxton order and broadening parameters;
Ci    tol  :allowed error in DOS due to truncating the gaussian,
Ci          if negative on entry, range is set to -tol*W
Ci    emin, emax, ndos; energy range and number of energy mesh points
Ci   esciss:Shift energy of unoccupied states by esciss 
Ci          (scissors operator)
Ci   jdos  :compute joint DOS, omitting matrix elements optmt
Ci   optmt :matrix elements of gradient operator
Ci   efermi:Fermi energy
Ci   ndos  :number of energy mesh points
Co  Output
Co    dos :joint density of states, or joint DOS weighted by optmt,
Co        :added to dos
Cr  Remarks
Cr    dos accumulates all transitions between unoccupied and 
Cr    occupied states.
Cr    NB: array dos is NOT initialized in this routine.  This enables
Cr    partial contribution to total dos in multiple calls.
Cr    Caller should initialize dos in first call
Cu Updates
Cu   27 Apr 04 Take abs(wgts) in conformity with use of sign as flag
Cu             (see bzmesh 15 Sep 02 modification)
Cu   22 Feb 04 Removed initialization of dos
C-----------------------------------------------------------------------
C     implicit none
      integer nqp,nbmx,nsp,n,ndos,nfilm,nempm,nfilo,nfiup,nemlo,nemup
      double precision wgts(nqp),evl(nbmx,nsp,nqp),dos(0:ndos-1,3,nsp),
     .  optmt(3,nfilm,nempm,nsp,*),w,emin,emax,tol,wt,emesh,efermi,
     .  esciss
      logical jdos
      integer i,isp,iq,meshpt,mesh1,mesh2,mrange,iprint,ib1,ib2,
     .  lgunit,stdo,npol,k,ibf,ibm,ib2lo
      double precision e,x,range,test,step,d,s,xx
      external delstp


      stdo = lgunit(1)
      npol = 3
      if (jdos) npol = 1
C     call dpzero(dos,3*nsp*ndos)
      step = (emax - emin) / (ndos - 1)
      if (emin .lt. 0) call info0(10,1,0,
     .  ' MKJDOS: (warning) emin<0 for joint DOS')

C --- Set tolerance and/or range of gaussians ---
      if ( tol .gt. 0d0 ) then
        do  1  i = 0, ndos-1
          x = i * step / w
          call delstp(0,x,test,s,xx)
          if ( test .lt. tol ) then
            mrange = i + 1
            goto 2
          endif
    1   continue
        call info0(20,1,0,' *** Warning : tol too small in MKJDOS')
    2   continue
        range = 2 * mrange * step
        test = tol
      else
        range = -tol * w
        mrange = range / ( 2 * step )
        call delstp(0,-tol/2,test,s,xx)
      endif

      if (iprint() .gt. 30) call awrit8(
     .  ' MKJDOS:  ef=%1;6d  N.W=%1;6d;  emin=%d  emax=%d  %i bins'//
     .  '%N%10frange of gaussians=%1;4d W (%i bins)  Est DOS err=%1;6d',
     .  ' ',200,stdo,efermi,n+w,emin,emax,ndos-1,range/w,2*mrange,
     .  test)

C --- Loop over spin and k-points ---
      do  8  isp = 1, nsp
        do  7  iq = 1, nqp
C ...     Double loop over occupied bands and unoccupied bands
          ibf = 0
          do  6  ib1 = nfilo, nfiup
            ibf = ibf+1
            if (evl(ib1,isp,iq) .gt. efermi) goto 6
            ibm = 0
            ib2lo = max0(ib1+1,nemlo)
            do  5  ib2 = ib2lo, nemup
              ibm = ibm+1
              if (evl(ib2,isp,iq) .lt. efermi) goto 5
              e = evl(ib2,isp,iq) - evl(ib1,isp,iq)
              if (e .gt. emax+step+range/2) goto 5
              meshpt = (e - emin) / step
              mesh1 = max(meshpt-mrange,0)
              mesh2 = min(meshpt+mrange,ndos-1)

C ...         Loop over polarizations (if jdos=T, npol=1)
              do  4  k = 1, npol
                wt = abs(wgts(iq))/nsp
                if (.not. jdos) wt = wt*optmt(k,ibf,ibm,isp,iq)
                do  3  meshpt = mesh1, mesh2
                  emesh = emin + meshpt * step
                  x = (emesh - e) / w
                  call delstp(n,x,d,s,xx)
                  dos(meshpt,k,isp) = dos(meshpt,k,isp) + wt * d / w
    3           continue
    4         continue

    5       continue
    6     continue
    7   continue
    8 continue

      end

