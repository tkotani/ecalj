      subroutine nosbzm(nfilem,nbmx,nsp,nspc,nchan,metal,n,w,nkp,wtbzm,
     .  eband,doswt,npln,nwmx,nqmx,nw,ew,nq,bzm)
C- make BZ maps from DOS weights
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nfilem, file handle for decomposition of the norm (MOMS or BAND)
Ci   nbmx, dimension of eband;
Ci   nsp=1 spin degenerate, =2 non-deg;
Ci   nspc: 2 for spin up and down coupled, otherwise 1
Ci   nchan, no. of channels;
Ci   metal; n, w : order and broadening; nkp no. of k-points;
Ci   wtbzm : weight per k-point
Ci   eband, work array for bands; doswt, work array for weights;
Ci   npln : total number of BZ planes
Ci   nwmx : maximum number of energy windows
Ci   nqmx : maximum number of k-points in one BZ plane
Ci   nw : number of energy windows for each BZ plane
Ci   ew(1), ew(2) : lower and upper limits of energy windows
Ci   nq(1), nq(2) : number of divisions along x and y of BZ plane
Co Outputs:
Co   bzm : integrated DOS for each k-point, spin, DOS channel,
Co         energy window, and BZ plane
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nfilem,nbmx,nsp,nspc,nchan,n,nkp,npln,nwmx,nqmx
      integer nw(npln),nq(2,npln)
      double precision w,wtbzm
      double precision eband(nbmx),doswt(nchan,nbmx,nspc),
     .  ew(2,nwmx,npln),bzm(nqmx,nsp,nchan,nwmx,npln)
      logical metal

C Local variables
      integer nspx,ntot,ip,nq0,iq,isp,nfstg,i,nevmx,iw,ib,jsp,ksp,ich
      integer iomoms
      double precision xx,emin,emax,e,x,d,s,w1,w2,wt,bwt

C --- Read past header in moments file ---
      rewind nfilem
      read(nfilem)

      nspx = nsp / nspc
      ntot = 0
      call dpzero(bzm,nqmx*nsp*nchan*nwmx*npln)

C --- Loop over BZ planes ---
      do  70  ip = 1, npln
        nq0 = nq(1,ip)*nq(2,ip)
        ntot = ntot + nq0
        call rxx(nq0 .gt. nqmx,'NOSBZM: nq0 gt nqmx')
        call rxx(ntot .gt. nkp,'NOSBZM: ntot gt nkp')
        call rxx(nw(ip) .gt. nwmx,'NOSBZM: nw gt nwmx')

C --- Get bands and weights from disc (loop over k-points and spin) ---
        do  60  iq = 1, nq0
          do  50  isp = 1, nspx
            nfstg = 11
            if (iomoms(nfilem,i,nsp,nspc,nkp,i,nfstg,nspc,1,1,nbmx,nbmx,
     .        nchan,nchan,nevmx,eband,doswt,doswt,doswt,xx,xx) .lt. 0)
     .        call rx('NOSBZM: failed to read moments')

C --- Loop over energy windows ---
            do  40  iw = 1, nw(ip)
              emin = ew(1,iw,ip)
              emax = ew(2,iw,ip)

C --- Get (sampling) weight for each band in energy window ---
              do  30  ib = 1, nevmx
                e = eband(ib)
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

                do  20  jsp = 1, nspc
C ... ksp is isp for uncoupled spins, and jsp for coupled spins ...
                  ksp = max(jsp,isp)

C --- Accumulate integrated DOS for each q-point and DOS channel ---
                  do  10  ich = 1, nchan
                    bwt = wt*wtbzm*doswt(ich,ib,jsp)
                    bzm(iq,ksp,ich,iw,ip) = bzm(iq,ksp,ich,iw,ip) + bwt
   10             continue
   20           continue

   30         continue

   40       continue

   50     continue
   60   continue

   70 continue

      call rxx(ntot .ne. nkp,'NOSBZM: ntot ne nkp')

      end

