      subroutine atwf(mode,a,lmxa,nr,nsp,pnu,pnz,rsml,ehl,rmt,z,v0,
     .nphimx,ncore,konfig,ecore,gcore,gval)
C- Make properties related to core for one sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 return ncore, and konfig, and nphimx only;
Ci         :  see description below for contents of nphimx
Ci         :1s digit
Ci         :1 return valence wave functions
Ci         :2 return core wave functions
Ci         :3 combination of 1+2
Ci         :10s digit concerns orthogonalization
Ci         :0 do not orthogonalize
Ci         :1 return orthogonalized to valence orbitals
Ci         :2 return orthogonalized to valence orbitals
Ci         :  using large component only
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   lmxa  :augmentation l-cutoff
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   pnz   :pnu por local orbitals
Ci   rmt   :MT boundary
Ci   z     :nuclear charge      (not used if mode=0)
Ci   v0    :spherical potential (not used if mode=0)
Ci   ehl   :energy of smoothed Hankel tail for extended local orbital
Ci   rsml  :corresponding smoothing radius for sm. Hankel tail, loc. orb
Cio Inputs/Outputs
Cio  nphimx:dimensions gval.  Must be at least as large as the
Cio        :number of valence wave functions
Cio        :For mode=0, nphimx is output and is assigned to
Ci         :maximum number radial wave functions for any l channel.
Co Outputs
Co   ncore :number of core levels
Co   konfig:1s digit contains core configuration
Co         :10s digit:
Co         : 0 -> no local orbitals
Co         : 1 -> local orbital with p.q.n. < pnu
Co         : 2 -> local orbital with p.q.n. > pnu
Co   ... The following are not used if mode=0
Co   ecore :core eigenvalues
Co   gcore :core wave functions
Co   gval  :valence wave functions
Co          gval(ir,l,i,isp) radial w.f. for (ir,l,isp) and:
Co            i=0 : phi
Co            i=1 : phidot
Co            i=2 : local orbital
Cr Remarks
Cu Updates
Cu    4 Sep 04 Adapted to extended local orbitals
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
Cu   22 Apr 01 Created by MvS
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nr,nsp,lmxa,ncore,konfig(1+lmxa),n0,nrmx,nphimx
      parameter (n0=10,nrmx=1501)
      double precision rmt,z,a,v0(nr,nsp),pnu(n0,nsp),pnz(n0,nsp),
     .gval(nr*2,0:lmxa,nphimx,nsp),ecore(*),gcore(nr,2,*),
     .rsml(n0),ehl(n0)
C ... Local parameters
      logical lpz
      integer l,isp,konf,konfz,k,mode0,mode1
      double precision sumtc,sumec,e,ez,xx
C     double precision hcrl,val(5),slo(5),pi,tol
C     parameter (tol=1d-12)
      double precision rofi(nrmx),rwgt(nrmx),rhoc(nrmx,2),gp(2*nrmx*4)
      double precision phi,dphi,phip,dphip,p,phz,dphz,phzp,dphzp

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)

C --- Count number of core states ---
      lpz = .false.
      ncore = 0
      do  l = 0, lmxa
        k = l+1
        konfig(k) = pnu(k,1)
        konfz = mod(pnz(k,1),10d0)
        if (konfz .eq. 0) konfz = konfig(k)
C       Sanity check
        call isanrg(konfz,konfig(k)-1,konfig(k)+1,'atwf:','pnuz',.true.)
C       lpz = konfz .ne. konfig(k)
        do  konf = l+1, min(konfz,konfig(k))-1
          ncore = ncore+nsp
        enddo
        if (konfz .lt. konfig(k)) then
          konfig(k) = konfz + 10
          lpz = .true.
        elseif (konfz .gt. konfig(k)) then
          konfig(k) = konfig(k) + 20
          lpz = .true.
        endif
      enddo

      if (mode0 .eq. 0) then
        nphimx = 2
        if (lpz) nphimx = 3
        return
      endif

      if (nr .gt. nrmx) call rx('increase nrmx in atwf')
      call radmsh(rmt,a,nr,rofi)
      call radwgt(rmt,a,nr,rwgt)

C --- Valence wave functions ---
      if (mod(mode0,2) .eq. 1) then
        do  l = 0, lmxa
          k = l+1
          do  isp = 1, nsp
            konf = pnu(k,1)

C    ...  Make phi and phidot
C         NB: Write gdot to gp, with extra space for higher derivatives
C         nn  = konf-l-1
C         pi = 4d0*datan(1d0)
C         hcrl = 0
C         val(1) = rofi(nr)
C         slo(1) = 1 + dtan(pi*(0.5d0 - pnu(k,isp)))
C         call phidx(0,z,l,v0(1,isp),hcrl,0d0,rofi,nr,4,tol,e,val,slo,
C    .      nn,gval(1,l,1,isp),gp,xx,xx,xx,xx,pgam,xx,xx,xx,xx)
            call makrwf(0,z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnu(1,isp),4,
     .      gval(1,l,1,isp),gp,e,phi,dphi,phip,dphip,p)
C         Copy 1st derivative to passed array
            call dcopy(2*nr,gp,1,gval(1,l,2,isp),1)
C         phi,phidot already orthogonal if mode1=1
            if (mode1 .eq. 2)
     .      call ortrwf(10*(mode1-1)+2,z,l,v0(1,isp),nr,nr,nr,rofi,rwgt,
     .      e,e,ez,gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp),xx)

C     ... Make local orbital
            if (konf .ne. konfig(k)) then
              call isanrg(nphimx,3,3,'atwf:','nphimx',.true.)
              call makrwf(0,z,rofi(nr),l,v0(1,isp),a,nr,rofi,pnz(1,isp),2,
     .        gval(1,l,3,isp),gp,ez,phz,dphz,phzp,dphzp,p)

              call isanrg(mode1,0,2,'atwf:','10s digit mode',.true.)
              if (mode1 .eq. 0) then

C             Extra scaling
C              call ortrwf(0,z,l,v0(1,isp),nr,nr,nr,rofi,rwgt,e,e,ez,
C     .          gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp),xx)
C             call dscal(nr*2,1/xx,gval(1,l,3,isp),1)
C             phz = phz/xx
C             dphz = dphz/xx

                call wf2lo(l,a,nr,rofi,rwgt,phi,dphi,phip,dphip,phz,dphz,
     .          phzp,dphzp,pnz(1,isp),rsml,ehl,
     .          gval(1,l,1,isp),gval(1,l,2,isp),gval(1,l,3,isp))
              elseif (pnz(l+1,isp) .lt. 10) then
                call ortrwf(10*(mode1-1)+1,z,l,v0(1,isp),nr,nr,nr,rofi,
     .          rwgt,e,e,ez,gval(1,l,1,isp),gval(1,l,2,isp),
     .          gval(1,l,3,isp),xx)
              endif
C           call prrmsh('gz',rofi,gval(1,l,3,isp),nr,nr,2)
            endif

          enddo
        enddo

C       call prrmsh('gval',rofi,gval,nr,nr,2*(1+lmxa))

      endif

C --- Core eigenfunctions and eigenvalues ---
      if (mode0 .ge. 2) then
        call getcor(1,z,a,pnu,pnz,nr,lmxa,rofi,v0,0,0,0d0,sumec,sumtc,
     .  rhoc,ncore,ecore,gcore)
      endif

      end

      subroutine atwf2l(ifi,jfi,iclass,a,lmxa,nr,nsp,pnz,rmt,
     .nphimx,konfig,ecore,gcore,gval)
C- Translate radial wave functions from shifted to standard log mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   :file logical unit for valence w.f.
Ci   jfi   :file logical unit for core w.f.
Ci   iclass:for printout
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci         :See Remarks
Ci   lmxa  :augmentation l-cutoff
Ci   nr    :number of radial mesh points
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   pnz   :p.q.n for local orbitals
Ci   rmt   :MT boundary
Ci   nphimx:dimensions gval.  Should be 2, or 3 if any local orbitals
Ci   konfig:1s digit contains core configuration
Ci         :10s digit:
Ci         : 0 -> no local orbitals
Ci         : 1 -> local orbital with p.q.n. < pnu
Ci         : 2 -> local orbital with p.q.n. > pnu
Ci   ecore :core eigenvalues
Ci   gcore :core wave functions
Ci   gval  :valence wave functions
Ci          gval(ir,l,i,isp) radial w.f. for (ir,l,isp) and:
Ci            i=0 : phi
Ci            i=1 : phidot
Ci            i=2 : local orbital
Cr Remarks
Cr   This routine translates w.f. to a standard log mesh and writes
Cr   them in a spex-readable format.
Cr
Cr   Debugging:
Cr   set sqr = '-p -p -xe -tog -coll 1 -tog -coll 2:nc -ccat'
Cr   mc -qr out.gas $sqr -int 0 2.361911
Cu Updates
Cu   15 Jul 09 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ifi,jfi,iclass,nr,nsp,lmxa,konfig(0:lmxa),n0,nphimx
      parameter (n0=10)
      double precision rmt,a,pnz(n0),gval(nr,2,0:lmxa,nphimx,nsp),
     .ecore(*),gcore(nr,2,*)
C ... Local parameters
      integer ic,icore,iphi,iprint,ir,isp,jcore,jx,konf,l,nglob,nphi,
     .npoly,nrmx,stdo,PRTV
      double precision xxo,xx,fac,tolg,dot3
      double precision x,y,dy,dytop,dytopa  ! For interpolation
      parameter (npoly=6,nrmx=1501,PRTV=60)
      parameter (tolg=1d-8)
      double precision
     .gwk(nrmx),gwk2(nrmx,2),  ! work arrays
     .rofi(nrmx),rwgt(nrmx),   ! Points and weights for shifted log mesh
     .rwgtl(nrmx),rofil(nrmx)  ! Points and weights for stnd log mesh
      double precision onorm(2,0:lmxa,nphimx,nsp), ! norm on shifted mesh
     .norm(2,0:lmxa,nphimx,nsp)  ! norm on log mesh

C ... Setup
      stdo = nglob('stdo')
      if (nr .gt. nrmx) call rx('increase nrmx in atwf')

C ... Shifted and standard log mesh
      call radmsh(rmt,a,nr,rofi)
      call radwgt(rmt,a,nr,rwgt)
      xx = rmt
      fac = exp(-a)
      do  ir = nr, 1, -1
        rofil(ir) = xx
        rwgtl(ir) = 2d0*a*xx/3d0   ! Simpson's rule for standard log mesh
C       rwgtl(ir) = 1d0*a*xx
        xx = xx*fac
      enddo
C     Complete Simpson's rule
      do  ir = 2, nr-1, 2
        rwgtl(ir) = 2d0*rwgtl(ir)
      enddo
      rwgtl(1) = rwgtl(1)/2
      rwgtl(nr) = rwgtl(nr)/2

c     call prrmsh('standard log mesh',rofi,rofil,nr,nr,1)

C      print *, '!! testing ... int x dx',
C     .  nr,sngl(a),sngl(rmt)
C      sum1 = 0
C      sum2 = 0
C      do  ir  = 1, nr
C        sum1 = sum1 + rofi(ir)*rwgt(ir)
C        sum2 = sum2 + rofil(ir)*rwgtl(ir)
C      enddo
C      print 333, rmt**2/2,sum1-rmt**2/2,sum2-rmt**2/2
C  333 format('exact',f15.10,'  shifted mesh error',1pe10.2,
C     .  '  standard mesh error',1pe10.2)
C      stop

C --- Valence wave functions ---
C     ic = 2*(lmxa+1)*nphimx*nsp
C     call prrmsh('gval on shifted log mesh',rofi,gval,nr,nr,ic)

      write(ifi) iclass,lmxa,nsp ! 1st record in class: dimensioning
      call info2(PRTV,1,0,' Valence wave function normalizations for '
     .//'class %i%N  phi  l   spin    sqrt(norm)',iclass,isp)
      dytopa = 0
      do  isp = 1, nsp
        do  l = 0, lmxa
          nphi = 2
          if (pnz(l+1) .gt. 0) nphi = 3
          do  iphi = 1, nphi
            write(ifi) l,iphi,isp   ! indices to wave function being written
            do  ic = 1, 2  ! large, small components
              onorm(ic,l,iphi,isp) =
     .        dot3(nr,gval(1,ic,l,iphi,isp),gval(1,ic,l,iphi,isp),rwgt)

C         g(r) = phi(r)*r
C         For small r, g(r)~r^(l+1) ... so fit g(r)/r^(l+1)
              do  ir = 2, nr
                gwk(ir) = gval(ir,ic,l,iphi,isp)/rofi(ir)**(l+1)
              enddo
              jx = 0
              dytop = 0
              do  ir = 1, nr
                x = rofil(ir)
                call polint(rofi(2),gwk(2),nr-1,npoly,x,tolg,0,jx,y,dy)
                dytop = max(dytop,abs(dy))
                dytopa = max(dytopa,abs(dy))
                gwk2(ir,ic) = y*rofil(ir)**(l+1)
C           gval(ir,ic,l,iphi,isp) = gwk2(ir,ic) ! Overwrite, for debuggging
              enddo
              norm(ic,l,iphi,isp) = dot3(nr,gwk2(1,ic),gwk2(1,ic),rwgtl)
C          Debugging printout
C          call info8(10,0,0,'interp ic=%i l=%i iphi=%i isp=%i:'//
C     .      ' err ~ %;g  onorm = %;6d  nnorm = %;6d',
C     .      ic,l,iphi,isp,dytop,onorm(ic,l,iphi,isp),
C     .      norm(ic,l,iphi,isp),0)

            enddo
            write(ifi) gwk2(1:nr,1),gwk2(1:nr,2) ! Large component, followed by small component
            if (iprint() .ge. PRTV) write(stdo,345) iphi, l, isp,
     .      dsqrt(norm(1,l,iphi,isp)+norm(2,l,iphi,isp))
  345       format(2i4,i6,f15.8,f15.6)
          enddo
        enddo
      enddo
      write(ifi) -1,-1,-1 ! Flags last record for this class
C     Debugging; need to comment overwrite above
C      ic = 2*(lmxa+1)*nphimx*nsp
C      call prrmsh('gval on standard log mesh',rofil,gval,nr,nr,ic)

C --- Core eigenfunctions and eigenvalues ---
      write(jfi) iclass,lmxa,nsp ! 1st record in class: dimensioning
      dytopa = 0
      call info2(PRTV,1,0,' Core wave function normalizations for '
     .//'class %i%N   n   l   spin    sqrt(norm)        ecore',
     .iclass,isp)
      icore = 0
      do  l = 0, lmxa
        do  isp = 1, nsp
          jcore = 0
          do  konf = l+1, mod(konfig(l),10)-1
            icore = icore+1
            jcore = jcore+1

            xxo =
     .      dot3(nr,gcore(1,1,jcore),gcore(1,1,jcore),rwgt) +
     .      dot3(nr,gcore(1,2,jcore),gcore(1,2,jcore),rwgt)

            dytop = 0
            do  ic = 1, 2 ! large, small components

C           g(r) = phi(r)*r
C           For small r, g(r)~r^(l+1) ... so fit g(r)/r^(l+1)
              do  ir = 2, nr
                gwk(ir) = gcore(ir,ic,jcore)/rofi(ir)**(l+1)
              enddo
              jx = 0
              do  ir = 1, nr
                x = rofil(ir)
                call polint(rofi(2),gwk(2),nr-1,npoly,x,tolg,0,jx,y,dy)
                if (ic .eq. 1) then
                  dy = dy*rofil(ir)**(l+1)
                  dytop  = max(dytop,dy)
                  dytopa = max(dytopa,dy)
                endif
                gwk2(ir,ic) = y*rofil(ir)**(l+1)
C             gcore(ir,ic,jcore) = gwk2(ir,ic) ! Overwrite, for debuggging
              enddo
            enddo

            xx = dot3(nr,gwk2(1,1),gwk2(1,1),rwgtl)
     .      + dot3(nr,gwk2(1,2),gwk2(1,2),rwgtl)
C           Debugging printout
            call info8(99,0,0,' atwf2l interp l=%i konf=%i isp=%i:'//
     .      ' err ~ %;g  onorm = %,6;6d  nnorm-onorm = %;6d',l,
     .      konf,isp,dytop,xxo,xx-xxo,0,0)

            write(jfi) jcore, l, isp, konf, ecore(icore)
            write(jfi) gwk2(1:nr,1), gwk2(1:nr,2)

            if (iprint() .ge. PRTV) write(stdo,345) konf, l, isp,
     .      dsqrt(xx), ecore(icore)

          enddo
        enddo
      enddo
      write(jfi) -1,-1,-1,-1,-1d0 ! Flags last record for this class
      end

