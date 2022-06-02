subroutine nosbzm(nfilem,nbmx,nsp,nspc,nchan,metal,n,w,nkp,wtbzm, &
     eband,doswt,npln,nwmx,nqmx,nw,ew,nq,bzm)
  !- make BZ maps from DOS weights
  ! ----------------------------------------------------------------------
  !i Inputs:
  !i   nfilem, file handle for decomposition of the norm (MOMS or BAND)
  !i   nbmx, dimension of eband;
  !i   nsp=1 spin degenerate, =2 non-deg;
  !i   nspc: 2 for spin up and down coupled, otherwise 1
  !i   nchan, no. of channels;
  !i   metal; n, w : order and broadening; nkp no. of k-points;
  !i   wtbzm : weight per k-point
  !i   eband, work array for bands; doswt, work array for weights;
  !i   npln : total number of BZ planes
  !i   nwmx : maximum number of energy windows
  !i   nqmx : maximum number of k-points in one BZ plane
  !i   nw : number of energy windows for each BZ plane
  !i   ew(1), ew(2) : lower and upper limits of energy windows
  !i   nq(1), nq(2) : number of divisions along x and y of BZ plane
  !o Outputs:
  !o   bzm : integrated DOS for each k-point, spin, DOS channel,
  !o         energy window, and BZ plane
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: nfilem,nbmx,nsp,nspc,nchan,n,nkp,npln,nwmx,nqmx
  integer :: nw(npln),nq(2,npln)
  double precision :: w,wtbzm
  double precision :: eband(nbmx),doswt(nchan,nbmx,nspc), &
       ew(2,nwmx,npln),bzm(nqmx,nsp,nchan,nwmx,npln)
  logical :: metal

  ! Local variables
  integer :: nspx,ntot,ip,nq0,iq,isp,nfstg,i,nevmx,iw,ib,jsp,ksp,ich
  integer :: iomoms
  double precision :: xx,emin,emax,e,x,d,s,w1,w2,wt,bwt

  ! --- Read past header in moments file ---
  rewind nfilem
  read(nfilem)

  nspx = nsp / nspc
  ntot = 0
  call dpzero(bzm,nqmx*nsp*nchan*nwmx*npln)

  ! --- Loop over BZ planes ---
  do  70  ip = 1, npln
     nq0 = nq(1,ip)*nq(2,ip)
     ntot = ntot + nq0
     call rxx(nq0 .gt. nqmx,'NOSBZM: nq0 gt nqmx')
     call rxx(ntot .gt. nkp,'NOSBZM: ntot gt nkp')
     call rxx(nw(ip) .gt. nwmx,'NOSBZM: nw gt nwmx')

     ! --- Get bands and weights from disc (loop over k-points and spin) ---
     do  60  iq = 1, nq0
        do  50  isp = 1, nspx
           nfstg = 11
           if (iomoms(nfilem,i,nsp,nspc,nkp,i,nfstg,nspc,1,1,nbmx,nbmx, &
                nchan,nchan,nevmx,eband,doswt,doswt,doswt,xx,xx) < 0) &
                call rx('NOSBZM: failed to read moments')

           ! --- Loop over energy windows ---
           do  40  iw = 1, nw(ip)
              emin = ew(1,iw,ip)
              emax = ew(2,iw,ip)

              ! --- Get (sampling) weight for each band in energy window ---
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
                    if (e >= emin .AND. e <= emax) wt = 1d0
                 endif

                 do  20  jsp = 1, nspc
                    ! ... ksp is isp for uncoupled spins, and jsp for coupled spins ...
                    ksp = max(jsp,isp)

                    ! --- Accumulate integrated DOS for each q-point and DOS channel ---
                    do  10  ich = 1, nchan
                       bwt = wt*wtbzm*doswt(ich,ib,jsp)
                       bzm(iq,ksp,ich,iw,ip) = bzm(iq,ksp,ich,iw,ip) + bwt
10                  enddo
20               enddo

30            enddo

40         enddo

50      enddo
60   enddo

70 enddo

  call rxx(ntot .ne. nkp,'NOSBZM: ntot ne nkp')

end subroutine nosbzm

