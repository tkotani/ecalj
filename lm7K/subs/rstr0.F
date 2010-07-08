Cdwdef1 have_common_w_in rstr0 , #_of_w_access= 0
      subroutine rstr0(nxi,lxi,exi,nlmx,np,x,y,z,lmxa,iop,hl,hd)
C- Reduced strux for a vector of energies, standard definition IV-43.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nxi   :number of energies for each connecting vector
Ci   lxi   :lmax for each energy and point
Ci   exi   :list of nxi energies
Ci   nlmx  :leading dimension of hl,hd
Ci   np    :number of connecting vectors
Ci   x,y,z :connecting vectors
Ci   lmxa  :Make hl for lxi + lmxa
Ci   iop   :1: calculate h only
Ci          any other number: calculate both h and hdot
Co Outputs
Co   hl    :reduced strux for each energy, point, to lxi + lmxa
Co   hd    :reduced energy derivative
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer np,nxi,lxi(nxi,np),lmxa,nlmx,iop
      double precision hl(nlmx,nxi,np),hd(nlmx,nxi,np),
     .exi(nxi),x(np),y(np),z(np)
C ... Local parameters
      integer:: ie , ip , lmax , ilm , l , m , lmxx , lmxy
      real(8) ,allocatable :: ylm_rv(:)
      real(8) ,allocatable :: r2_rv(:)

      parameter (lmxx=20)
      double precision psi(-1:lmxx),phi(-1:lmxx),yl((lmxx+1)**2),r2,
     .rfac,xxx,fac2l(0:lmxx),psid,psil
C ... Heap
Cdwdef1       integer w(1)
Cdwdef1       common /w/ w

C --- Make ylm for all points, up to lmxy, which find ---
      lmxy = -1
      do  1  ip = 1, np
        do  1  ie = 1, nxi
          lmxy = max(lmxy,lxi(ie,ip))
    1 continue
      lmxy = lmxy+lmxa
      fac2l(0) = 1
      do  2  l = 1, lmxy
        fac2l(l) = fac2l(l-1) * (2*l-1)
    2 continue
      if (lmxy .gt. lmxx) call rxi('rstr0: lmax exceeds lmxx:',lmxy)
      allocate(ylm_rv((lmxy+1)**2*np))
Ctakao_ZeroClear_NotRequiered        if ((lmxy+1)**2*np<0) ylm_rv(:)=0.0d0

      allocate(r2_rv(np))
Ctakao_ZeroClear_NotRequiered        if (np<0) r2_rv(:)=0.0d0

      call ropyln ( np , x , y , z , lmxy , np , ylm_rv , r2_rv 
     .)


C --- For each point, do ---
      do  5  ip = 1, np
        call pvstr0 ( lmxy , ip , np , ylm_rv , r2_rv , yl , r2 
     .  )

        if (r2 .lt. 1d-10) then
          call dpzero(hl(1,1,ip),nlmx*nxi)
          if (iop .ne. 1) call dpzero(hd(1,1,ip),nlmx*nxi)
          goto 5
        endif

C --- Reduced strx hl, or hl and hd for all energies, this point ---
        if (iop .eq. 1) then
          do  10  ie = 1, nxi
            lmax = lmxa+lxi(ie,ip)
C     ... Not worth vectorizing (I think)
c         call bessl2(exi(ie)*r2,0,lmax,phi(0),psi(0))
            call besslr(exi(ie)*r2,0,-1,lmax,phi,psi)
            ilm = 0
            rfac = dsqrt(r2)
            xxx = 1d0/r2
            do  20  l = 0, lmax
              rfac = rfac*xxx
C     ... undo fac2l scaling, to recover standard MSM IV-43.
C         psil = rfac * psi(l) * fac2l(l)
              psil = rfac * psi(l)
              do  20  m = -l, l
                ilm = ilm+1
                hl(ilm,ie,ip) = psil*yl(ilm)
   20       continue
   10     continue
        else
          do  110  ie = 1, nxi
            lmax = lmxa+lxi(ie,ip)
C         call bessl2(exi(ie)*r2,-1,lmax,phi(-1),psi(-1))
            call besslr(exi(ie)*r2,0,-1,lmax,phi(-1),psi(-1))
            ilm = 0
            rfac = dsqrt(r2)
            xxx = 1d0/r2
            do  120  l = 0, lmax
              rfac = rfac*xxx
C         psil = rfac * psi(l) * fac2l(l)
C         psid = rfac * psi(l-1)*r2/(4*l-2) * fac2l(l)
              psil = rfac * psi(l)
              psid = rfac * psi(l-1)*r2/2
              do  120  m = -l, l
                ilm = ilm+1
                hl(ilm,ie,ip) = psil*yl(ilm)
                hd(ilm,ie,ip) = psid*yl(ilm)
  120       continue
  110     continue
        endif
    5 continue
      if (allocated(r2_rv)) deallocate(r2_rv)
      if (allocated(ylm_rv)) deallocate(ylm_rv)

      end
      subroutine pvstr0(lmax,ip,nd,ylm,r2,yl,rsq)
C     implicit none
      integer lmax,ip,nd,nlm,ilm
      double precision ylm(nd,1), yl(1), r2(ip), rsq

      nlm = (lmax+1)**2
      rsq = r2(ip)
      do  10  ilm = 1, nlm
        yl(ilm) = ylm(ip,ilm)
   10 continue
      end

