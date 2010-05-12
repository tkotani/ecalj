      integer function iofa(spid,nxi0,nxi,exi,hfc,hfct,rsm,z,rmt,a,nr,
     .  qc,ccof,ceh,stc,rho,rhoc,v,ifi)
C- I/O for free-atom data, one species
C ----------------------------------------------------------------------
Ci Inputs
Ci   spid  :species label
Ci   nxi0  :leading dimension of hfc,hfct
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio File I/O
Cio  nxi   :number of energies used to fit tail of valence density
Cio  exi   :energies that fit tail of valence density
Cio  hfc   :coefficients for fit of tail valence density
Cio  hfct  :coefficients for fit of tail valence+core density
Cio  rsm   :smoothing radius for fit of tail valence density
Cio  z     :nuclear charge
Cio  rmt   :muffin-tin radius, in a.u.
Cio  a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Cio  nr    :number of radial mesh points
Cio  qc    :Sphere core charge
Cio  ccof  :coefficient to core density fit by unsmoothed Hankel
Cio  ceh   :hankel function energy to core density fit
Cio  stc   :core kinetic energy
Cio  rho   :valence density
Cio  rhoc  :core density
Cio  v     :spherical potential
Cr Remarks
Cu Updates
Cu   10 Jun 00 spin polarized
Cu   20 May 00 adapted from nfp rw_fa.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ifi,nr,nxi,nxi0
      character spid*8
      double precision a,ccof,ceh,qc,rmt,rsm,stc,z,
     .  hfc(nxi0,2),rho(2*nr),v(2*nr),exi(nxi0),rhoc(2*nr),hfct(nxi0,2)
C ... Local parameters
      integer i,jfi,lrel,lrel0,nglob,nsp,nsp0,ipr,iprint,stdo,lgunit
      logical isanrg
      character msg*23

      ipr    = iprint()
      stdo   = lgunit(1)
      nsp  = nglob('nsp')
      lrel = nglob('lrel')
      msg  = '         File mismatch:'
      iofa   = -1

C --- Input ---
      if (ifi .gt. 0) then

        jfi = ifi
        read(jfi,201,end=998,err=998) spid
  201   format(9x,a8/)
        read(jfi,102) z,a,nsp0,lrel0,nr,rmt,rsm
C ...   Sanity checks
        if (isanrg(lrel0, lrel,lrel,msg,'lrel',.true.)) stop
        if (isanrg(nsp0,  nsp,nsp,  msg,'nsp', .true.)) stop
        read(jfi,103) nxi
        read(jfi,205) (exi(i),i=1,nxi)
        read(jfi,205) (hfc(i,1),i=1,nxi)
        read(jfi,205) (hfct(i,1),i=1,nxi)
  205   format(5d16.8)
        if (nsp .eq. 2) read(jfi,*) (hfc(i,2),i=1,nxi)
        if (nsp .eq. 2) read(jfi,*) (hfct(i,2),i=1,nxi)
        read (jfi,210) qc,ccof,ceh,stc
  210   format(5x,4f16.7)
        call dfdump(rho,nr,ifi)
        if (nsp .eq. 2) call dfdump(rho(1+nr),nr,ifi)
        call dfdump(rhoc,nr,ifi)
        if (nsp .eq. 2) call dfdump(rhoc(1+nr),nr,ifi)
        call dfdump(v,nr,ifi)
        if (nsp .eq. 2) call dfdump(v(1+nr),nr,ifi)
      endif

C --- Output ---
      if (ifi .lt. 0)  then
        jfi = -ifi
        write(jfi,101) spid
  101   format('-------- ',a8,' ----------'/
     .     '     z       a      nsp lrel  nr     rmt          rsm')
        write(jfi,102) z,a,nsp,lrel,nr,rmt,rsm
  102   format(f8.2,f9.4,3i5,2f13.7)
        write(jfi,103) nxi
  103   format(i4)
        write(jfi,105) (exi(i),i=1,nxi)
        write(jfi,105) (hfc(i,1),i=1,nxi)
        write(jfi,105) (hfct(i,1),i=1,nxi)
        if (nsp .eq. 2) write(jfi,105) (hfc(i,2),i=1,nxi)
        if (nsp .eq. 2) write(jfi,105) (hfct(i,2),i=1,nxi)
  105   format(1p,5d16.8)
        write (jfi,110) qc,ccof,ceh,stc
  110   format(' core',4f16.7)
        call dfdump(rho,nr,-jfi)
        if (nsp .eq. 2) call dfdump(rho(1+nr),nr,-jfi)
        call dfdump(rhoc,nr,-jfi)
        if (nsp .eq. 2) call dfdump(rhoc(1+nr),nr,-jfi)
        call dfdump(v,nr,-jfi)
        if (nsp .eq. 2) call dfdump(v(1+nr),nr,-jfi)
        return
      endif

      iofa = 0
      return

C ... Error handling
  998 if (ipr .gt. 0)
     .  write(stdo,'('' iofa  : missing species id ... nothing read'')')

      end

