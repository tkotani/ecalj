      subroutine v0intr(nl,nsp,lmax,z,rhozbk,rmax,a,nr,rofi,
     .  pnu,qnu,v,rhoin,rho,rhoc,g,gp,nmix,niter,qc,lfrz,avw,
     .  ekap,mode,vintra)
C- Make derivatives of total energy or C w.r.t. q_j
C ----------------------------------------------------------------
Ci Inputs
Ci   Same as atomsc, with the following additions:
Ci   mode  1  Make dC_i/dq_j
Ci         2  Make d^2 E /dq_i dq_j
Co Outputs
Co   vintra   (mode 1) dC_i/dq_j
Co            (mode 2) d^2 E /dq_i dq_j
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nr,nsp,nl,nmix,niter,lmax,mode
      double precision rofi(nr,2),v(nr,nsp),rho(nr,nsp),rhoc(nr,nsp),
     .  rhoin(nr,nsp),pnu(nl,nsp),qnu(3,nl,nsp),g(nr),gp(nr,4),rhozbk,
     .  avw,ekap,vintra(nl,nl,nsp,nsp)
      logical lfrz
C local variables
      integer isp,jsp,il,jl,k,ki,kj,n0,off,kl,ksp,ov,iprint,jmin,ncmx,
     .  nvmx
      parameter (n0=10, ncmx=50, nvmx=20)
      double precision sumec,sumtc,ekin,utot,vrmax(2),rhrmx,df,sumev,
     .  rhoeps,etot,amgm,qtot,ves,xx(2),thrpv,z,rmax,a,qc,exrmax(2),
     .  thrpvl(n0),etot0,pp(6,n0*2),cnl(2,n0),enl(2,2),fi,fj,qcut,
     .  pl(n0*2),ql(3,n0*2),ec(ncmx),ev(nvmx)
C     double precision pz(n0,2),qz(3,n0,3)
      integer idmod(n0),idmoz(n0),ovk(2),idu(5)
      parameter (qcut=.05d0)
    
C     logical aiova,ltmp
C     character*120 outs
C Heap allocation
      integer w(1)
      common /w/ w
C ... Offset of index kl,ksp for array of dimension (nl,nsp)
      off(kl,ksp) = kl + nl*(ksp-1)

      call tcn('v0intr')
      if (nl .gt. n0) call rx('v0intr: increase n0')
      call dpzero(vintra,nl*nl*nsp*nsp)
      call iinit(idu,5)
      ves = 0
      xx(1) = 0
      xx(2) = 0
      call dcopy(3*nl*nsp,qnu,1,ql,1)
      do  2  il = 1, nl
    2 idmod(il) = 0
      call iinit(idmoz,n0)
C     npan = 1
      df = .001d0
      call defdr(ov,     nr*nsp)
      call defdr(ovk(1), nr*nsp)
      call defdr(ovk(2), nr*nsp)

      call pshpr(iprint()-40)
      call getqvc(nsp,nl,lmax,z,pnu,qnu,0,0,0,0,w,qc,qtot,amgm,0d0,0d0)
      ec(1) = 0

C --- Derivatives of C_i wrt q_j ---
      if (mode .eq. 1) then
      do  10  isp = 1, nsp
      do  10  jl = 1, nl
        do  14  k = 1, 2
          if (k .eq. 1) fj= df
          if (k .eq. 2) fj=-df
          call dcopy(nl*nsp,pnu,1,pl,1)
          call dcopy(nr*nsp,v,1,w(ov),1)
          ql(1,off(jl,isp)) = ql(1,off(jl,isp)) + fj

          call atomsc(.true.,nl,nsp,lmax,z,rhozbk,0,0,0d0,rmax,a,nr,
     .      rofi,ec,ev,pl,ql,idmod,w(ov),0d0,rhoin,rho,rhoc,nmix,qc,
     .      sumec,sumtc,sumev,ekin,utot,rhoeps,etot,amgm,rhrmx,vrmax,
     .      qtot,exrmax,'pot',niter,lfrz)

          call potpar(nl,nsp,lmax,z,rmax,avw,ekap,.false.,.false.,
     .      .false.,a,nr,rofi,w(ov),pl,idmod,ves,xx,ql,idu,xx,xx,thrpv,
     .      thrpvl,g,gp,pp,xx,xx,xx,xx)

          ql(1,off(jl,isp)) = ql(1,off(jl,isp)) - fj
          do  16  il = 1, nl
   16     cnl(k,il) = pp(2,off(il,isp))
   14   continue
        do  18  il = 1, nl
   18   vintra(jl,il,isp,isp) = (cnl(1,il)-cnl(2,il))/2d0/df
   10 continue

C --- d^2 E /dq_i dq_j ---
      elseif (mode .eq. 2) then

C ... Make etot(dq=0)
      call dcopy(nl*nsp,pnu,1,pl,1)
      call dcopy(nr*nsp,v,1,w(ov),1)
      call atomsc(.true.,nl,nsp,lmax,z,rhozbk,0,0,0d0,rmax,a,nr,rofi,ec,
     .  ev,pl,ql,idmod,w(ov),0d0,rhoin,rho,rhoc,nmix,qc,sumec,sumtc,
     .  sumev,ekin,utot,rhoeps,etot0,amgm,rhrmx,vrmax,qtot,exrmax,'pot',
     .  niter,lfrz)

      do  20  isp = 1, nsp
      do  20  jsp = isp, nsp
C ... Forget the off-diagonal for now
C     if (jsp .ne. isp) goto 20

      do  22  il = 1, nl
      jmin = il
      if (jsp .ne. isp) jmin = 1
      do  22  jl = jmin, nl

        vintra(jl,il,jsp,isp) = 0
        vintra(il,jl,isp,jsp) = 0
        if (ql(1,off(il,isp))+ql(1,off(jl,jsp)) .lt. qcut/nsp) goto 22

C   ... Accumulate E++, E+-, E-+, E--
        do  24  ki = 1, 2
          if (ki .eq. 1) fi= df
          if (ki .eq. 2) fi=-df
          ql(1,off(il,isp)) = ql(1,off(il,isp)) + fi
          do  26  kj = 1, 2
            if (kj .eq. 1) fj= df
            if (kj .eq. 2) fj=-df
            call dcopy(nl*nsp,pnu,1,pl,1)
            call dcopy(nr*nsp,v,1,w(ov),1)
C       ... Make a better guess for initial v on second pass
            if (ki .eq. 2) then
              call dscal(nr*nsp,2d0,w(ov),1)
              call daxpy(nr*nsp,-1d0,w(ovk(3-kj)),1,w(ov),1)
            endif
            ql(1,off(jl,jsp)) = ql(1,off(jl,jsp)) + fj
C       ... Self-consistent total energy for these moments
            call atomsc(.true.,nl,nsp,lmax,z,rhozbk,0,0,0d0,rmax,a,nr,
     .        rofi,ec,ev,pl,ql,idmod,w(ov),0d0,rhoin,rho,rhoc,nmix,qc,
     .        sumec,sumtc,sumev,ekin,utot,rhoeps,etot,amgm,rhrmx,vrmax,
     .        qtot,exrmax,'pot',niter,.true.)
C            call awrit4('%x  ki,kj %i %i il jl %i %i',
C     .        outs,120,0,ki,kj,il,jl)
C            call awrit3('%a ql %,4;4d %,4;4d %,4;4d',
C     .        outs,120,0,ql(1,1),ql(1,2),ql(1,3))
C            if (nsp .eq. 2)
C     .        call awrit3('%a %,4;4d %,4;4d %,4;4d',
C     .        outs,120,0,ql(1,4),ql(1,5),ql(1,6))
C            call awrit1('%a e %,10;10d',outs,120,-6,etot)
C       ... Setup for a better guess for initial v for second pass
            if (ki .eq. 1) call dcopy(nr*nsp,w(ov),1,w(ovk(kj)),1)
            ql(1,off(jl,jsp)) = ql(1,off(jl,jsp)) - fj
            enl(ki,kj) = etot
   26     continue
          ql(1,off(il,isp)) = ql(1,off(il,isp)) - fi
   24     continue

C   ... E'' by finite difference E++ - E+- - E-+ + E--
        vintra(jl,il,jsp,isp) = (enl(1,1) - enl(1,2) -
     .                           enl(2,1) + enl(2,2))/(2*df)**2
        vintra(il,jl,isp,jsp) = vintra(jl,il,jsp,isp)

   22 continue
   20 continue

      else
        call fexit(-1,111,' Exit -1 v0intr: bad mode',mode)
      endif

      call rlse(ov)
      call poppr

C      ltmp = aiova(alabl,vintra,nl,nl-1,nsp,-6)
C      stop

      call tcx('v0intr')
      end

