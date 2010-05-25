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
      implicit none
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

