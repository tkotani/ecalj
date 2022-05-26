      subroutine sugcut(mode) !,sspec)
      use m_lmfinit,only: nspec,alat=>lat_alat,tol=>lat_tolft,n0,nkap0,nkaphh,lhh,nkapii,sspec=>v_sspec
      use m_supot,only: gv=>rv_a_ogv,ng=>lat_ng
      use m_uspecb,only:uspecb
      use m_lgunit,only:stdo
C- Find max recip for each spec and orbital block, store in struct.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 make cutoffs for standard envelope functions
Ci         :2 make cutoffs for extended local orbitals
Ci         :3 combination 1+2
Co Outputs
Co    sspec%ngcut
Cu Updates
Cu   16 Aug 04 New mode for getting cutoffs, local orbs.
Cu             Changed argument list
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu    9 May 00 Adapted from nfp su_gvcut.f
C ----------------------------------------------------------------------
      implicit none
      integer:: mode, ngcut(n0,nkap0),lh(nkap0)
      integer:: ipr,iprint,is,irep,icut,i,ik,l,lcut,nkapi,nkap1,nkap2
      real(8):: rsmh(n0,nkap0),eh(n0,nkap0),tpi,tpiba2,gg0,gg,e,rsm,gam,gmax,top
      character*8 spid
      character*1 ccc,ccl
      if (ng .eq. 0) return
c      stdo = lgunit(1)
      ipr = iprint()
      tpi = 8d0*datan(1d0)
      tpiba2 = (tpi/alat)**2
      if (ipr .ge. 20) then
        if (mode .eq. 1) write(stdo,887) tol
        if (mode .eq. 2) write(stdo,888) tol
        write(stdo,774)
  887   format(/' sugcut:  make orbital-dependent reciprocal vector',' cutoffs for tol=',1p,e9.2)
  888   format(/' sugcut:  orbital-dependent cutoffs for local',' orbitals, tol=',1p,e9.2)
      endif
      gg = -1
      do  is = 1, nspec
        spid = sspec(is)%name
        nkap1 = 1
        call uspecb(is,rsmh,eh)
        nkap2 = nkaphh(is)
        nkapi = nkapii(is)
        if (mode>1) then
          if(mode==2) nkap1 = nkapi+1
        else
          nkap2 = nkapi
        endif
        ngcut=0 
        if(mode==2 ) ngcut=sspec(is)%ngcut
        gg0 = gg
        do  ik = nkap1, nkap2
          lcut = -1
          do  l  = 0, lhh(ik,is)
            e = eh(l+1,ik)
            rsm = rsmh(l+1,ik)
            if (rsm .ne. 0) then
              if (l .lt. lhh(ik,is) .and. l .gt. lcut) then
                lcut = l-1
   12           lcut = lcut+1
                if (lcut .lt. lhh(ik,is)) then
                  if (rsmh(lcut+2,ik).eq.rsm .and. eh(lcut+2,ik).eq.e) goto 12
                endif
              endif
C     ... Get cutoff radius where exp(-gam*gmax)*gmax**l equals tol
              gam = rsm*rsm/4d0
              gmax = 1d0
              do  irep = 1, 10
                gmax = dsqrt(-dlog(tol/gmax**l)/gam)
              enddo
C     ... Find first longer vector, icut is one less
              icut = 1
              do  i = 1, ng
                gg = tpiba2*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
                if (gg .gt. gmax*gmax) goto 90
                icut = i
                gg0 = gg
              enddo
   90         continue
              top = dexp(-gam*gg0)*dsqrt(gg0)**l
              ccc = ' '
              if (icut .eq. ng) ccc = '*'
              ccl = ' '
              if (l .lt. lcut) ccl = '*'
              if (ipr .ge. 20) write(stdo,773) spid,l,ccl,rsm,e,gmax,top,icut,ccc
  773         format(2x,a,i2,a1,f7.2,f7.2,f8.3,1p,e12.2,0p,i8,a)
  774         format(' spec      l    rsm',4x,'eh',5x,'gmax',4x,'last term',4x,'cutoff')
              ngcut(l+1,ik) = icut
            endif
          enddo
        enddo
        sspec(is)%ngcut=ngcut
      enddo
      end subroutine sugcut
      
