Cdelw1 var ['oylm', '0', 'oylm', '(lmxy+1)**2*np', 'real(8)', 'rv_w_', 'rstr0']
Cdelw1 var ['or2', '0', 'or2', 'np', 'real(8)', 'rv_w_', 'rstr0']
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
     .  exi(nxi),x(np),y(np),z(np)
C ... Local parameters
Cdelw1       integer ie,ip,lmax,ilm,l,m,lmxx,oylm,or2,lmxy
       integer:: ie , ip , lmax , ilm , l , m , lmxx , lmxy 
       real(8) ,allocatable :: rv_w_oylm(:)
       real(8) ,allocatable :: rv_w_or2(:)

      parameter (lmxx=20)
      double precision psi(-1:lmxx),phi(-1:lmxx),yl((lmxx+1)**2),r2,
     .  rfac,xxx,fac2l(0:lmxx),psid,psil
C ... Heap
      integer w(1)
      common /w/ w

C --- Make ylm for all points, up to lmxy, which find ---
      lmxy = -1
      do  1  ip = 1, np
      do  1  ie = 1, nxi
    1 lmxy = max(lmxy,lxi(ie,ip))
      lmxy = lmxy+lmxa
      fac2l(0) = 1
      do  2  l = 1, lmxy
    2 fac2l(l) = fac2l(l-1) * (2*l-1)
      if (lmxy .gt. lmxx) call rxi('rstr0: lmax exceeds lmxx:',lmxy)
Cdelw1       call defrr(oylm,(lmxy+1)**2*np)
       allocate(rv_w_oylm((lmxy+1)**2*np))
       if ((lmxy+1)**2*np<0) rv_w_oylm(:)=0.0d0

Cdelw1       call defrr(or2,np)
       allocate(rv_w_or2(np))
       if (np<0) rv_w_or2(:)=0.0d0

Cdelw1       call ropyln(np,x,y,z,lmxy,np,w(oylm),w(or2))
       call ropyln ( np , x , y , z , lmxy , np , rv_w_oylm , rv_w_or2 
     . ) 


C --- For each point, do ---
      do  5  ip = 1, np
Cdelw1       call pvstr0(lmxy,ip,np,w(oylm),w(or2),yl,r2)
       call pvstr0 ( lmxy , ip , np , rv_w_oylm , rv_w_or2 , yl , r2 
     . ) 

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
   20     hl(ilm,ie,ip) = psil*yl(ilm)
   10   continue
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
  120     continue
  110   continue
      endif
    5 continue
Cdelw1 rlse name= oylm old_list= oylm or2 
Cdelw1 rlse name= oylm new_list= (None)
Cdelw1       call rlse(oylm)
       if (allocated(rv_w_or2)) deallocate(rv_w_or2)
       if (allocated(rv_w_oylm)) deallocate(rv_w_oylm)

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
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
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      end

