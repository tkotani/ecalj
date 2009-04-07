      subroutine makbzm(nqp,nband,nbmx,nsp,evl,metal,n,w,npln,nwmx,
     .  nqmx,nw,ew,nq,wtbzm,bzm)
C- make BZ map from bands
C-----------------------------------------------------------------------
Ci Inputs
Ci   nqp : number of q-points; nband : number of bands;
Ci   nbmx : first dimension of evl ;
Ci   nsp : 2 for spin polarised bands, 1 otherwise;
Ci   evl : bands (eigenvalues); metal;
Ci   N, W : Methfessel-Paxton order and broadening parameters;
Ci   npln : total number of BZ planes;
Ci   nwmx : maximum number of energy windows;
Ci   nqmx : maximum number of k-points in one BZ plane;
Ci   nw : number of energy windows for each BZ plane;
Ci   ew(1), ew(2) : lower and upper limits of energy windows;
Ci   nq(1), nq(2) : no. of q-point divisions along x and y, see getbzp;
Ci   wtbzm : weight per k-point;
Co Ouputs
Co   bzm : integrated DOS as a function of k;
C-----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nqp,nband,nbmx,nsp,n,npln,nwmx,nqmx
      integer nw(npln),nq(2,npln)
      double precision w,wtbzm
      double precision evl(nbmx,nsp,nqp),ew(2,nwmx,npln),
     .  bzm(nqmx,nwmx,npln)
      logical metal

C Local variables
      integer i,isp,iband,iq,ip,iw,nq0,ntot
      double precision e,x,d,s,emin,emax,w1,w2,wt,xx
      external delstp

      call tcn('makbzm')

      iq = 0
      ntot = 0
      call dpzero(bzm,nqmx*nwmx*npln)

C --- Loop over BZ planes ---
      do  50  ip = 1, npln
        nq0 = nq(1,ip)*nq(2,ip)
        ntot = ntot + nq0
        call rxx(nq0 .gt. nqmx,'MAKBZM: nq0 gt nqmx')
        call rxx(ntot .gt. nqp,'MAKBZM: ntot gt nqp')
        call rxx(nw(ip) .gt. nwmx,'MAKBZM: nw gt nwmx')

C --- Loop over q-points and energy windows ---
        do  40  i = 1, nq0
          iq = iq + 1
          do  30  iw = 1, nw(ip)
            emin = ew(1,iw,ip)
            emax = ew(2,iw,ip)

C --- Accumulate integrated DOS (sampling) for each band ---
            do  20  iband = 1, nband
              do  10  isp = 1, nsp
                e = evl(iband,isp,iq)
                if (metal) then
                  x = (emin - e) / w
                  call delstp(n,x,d,s,xx)
                  w1 = 1d0 - s
                  x = (emax - e) / w
                  call delstp(n,x,d,s,xx)
                  w2 = 1d0 - s
                  wt = w2 - w1
                else
                  wt = 0d0
                  if (e .ge. emin .and. e .le. emax) wt = 1d0
                endif
                bzm(i,iw,ip) = bzm(i,iw,ip) + wt*wtbzm
   10         continue
   20       continue
   30     continue
   40   continue

   50 continue

      call rxx(ntot .ne. nqp,'MAKBZM: ntot ne nqp')
      call tcx('makbzm')

      end

