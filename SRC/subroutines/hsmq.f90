subroutine hsmq(nxi,lmxa,lxi,exi,rsm,job,q,p,nrx,nlmx,wk,yl,awald,alat,qlv,nG,dlv,nD,vol,hsm,hsmp)
  use m_ropyln,only: ropyln
  use m_shortn3_plat,only: shortn3_plat,nout,nlatout,plat=>platx,qlat=>qlatx
  intent(in)::  nxi,lmxa,lxi,exi,rsm,job,q,p,nrx,nlmx,      awald,alat,qlv,nG,dlv,nD,vol
  !- Bloch-sum of smooth Hankel functions and energy derivatives at p
  !  by Ewald summation, for nxi groups of parameters (lxi,exi,rsm).
  ! ---------------------------------------------------------------
  !i Inputs:
  !i  nxi,lxi,exi,rsm:  number of (l,energy,rsm) groups, and list
  !i  dlv,nD:  direct lattice vectors, and number
  !i  qlv,nG:  reciprocal lattice vectors, and number
  !i  awald,vol: Ewald parameter and cell volume (used for Q-space sum)
  !i  lmxa:    make strux for lxi(ie)+lmxa,
  !i           skipping over any energy for which lxi(ie) < 0
  !i  job:     ones digit intended for H_smoothed - H_unsmoothed
  !i           To do, set nG=0, and the one's digit as follows:
  !i            0, generate wk(1,2,3,4)
  !i           >0  assume wk(1,3,4) and yl have been generated already
  !i           tens digit handles whether or not to add H(p=0).
  !i            0: return whole Bloch function (not allowed for rsm=0)
  !i            1: return Bloch function less H(p=0), if p=0
  !i           >1: return Bloch function less H(p=0), for any p
  !i               DOES NOT WORK!
  !i           100s digit
  !i            0: use rsm = rsm(ie)
  !i            1: use rsm = rsm(1) is the same for all energies
  !i           1000s digit
  !i            0: use standard phase convention, do not shorten p
  !i            1: use standard phase convention, but shorten p
  !i            2: scale standard phase convention by exp(-i q . p)
  !i            3: like 2, but also shorten p
  !i  q,p:     k-point and connecting vector (units of 2*pi/alat and alat)
  !i  nrx:     leading dimension of wk,yl; must be at least max(nD,nG)
  !i  nlmx:    leading dimension of hsm,hsmp
  !i  yl:      work array, dimensioned nrx*(lmax+1)**2
  !i  wk:      work array, dimensioned at least nrx*(2*lmax+10)
  !o Outputs:
  !o   yl:     ylm(1..nD,1..(lmax+1)**2) for points alat*(p-dlv(1..nD))
  !o   wk:     (*,1) holds r**2
  !o           (*,2) holds Y0 exp(-(r/rsm(nxi))**2) (rsm>0 only)
  !o           (*,3) holds cos(q.dlv)
  !o           (*,4) holds sin(q.dlv)
  !o   hsm,hsmp:smoothed Hankels and energy derivatives for
  !o           nxi energies, to lmax lxi(1..nxi).
  !b Bugs and Limitations:
  !b   Convergence is poor for p small but nonzero.  To avoid this,
  !b   hsmq SHOULD re-evaluate points for which p is small using
  !b   a power series expansion.
  !u Updates
  !u  15 Aug 00 extended to e>0; added 1000 digit job
  !u   7 May 98 handles case rsm>1/a by strict q-space summation
  !u  16 Sep 98 fixed dimensioning error arising when nD > nG
  ! ---------------------------------------------------------------
  !     implicit none
  integer :: nxi,lmxa,lmax,nG,nD,nlmx,nrx,lxi(nxi),job
  double precision :: alat,awald,vol,exi(nxi),rsm(nxi),p(3),q(3), &
       wk(nrx,*),yl(nrx,*),qlv(3,*),dlv(3,*), hsm(2,nlmx,*),hsmp(2,nlmx,*)
  ! Local variables
  integer :: ie,ir,ilm,l,m,ir1,lc,ls,job0,job1,job2,job3,lm,nlmxx, &
       lx(20),ndx
  parameter (nlmxx=(200+1)**2) !(nlmxx=(16+1)**2)
  double precision :: qdotr,y0,a2,pi,sp,gam,tpiba,tpi,rsmi, &
       x1,x2,xx,xx0,xx1(nlmxx),xx2(nlmxx),xx3(nlmxx),xx4(nlmxx), &
       r,akap=1d99,kappa,h0,h0d,a,faca,p1(3),cosp,sinp
  parameter (faca=1d0)
  double complex xxc
  logical :: dcmpre,ltmp
  real(8):: pp(3)
  dcmpre(x1,x2) = dabs(x1-x2) .lt. 1d-8

  ! --- Setup ---
  job0 = mod(job,10)
  if (nG /= 0) job0 = 0
  job1 = mod(job/10,10)
  job2 = mod(job/100,10)
  job3 = mod(job/1000,10)
  ! ... Shorten connecting vector; adjust phase later
  if (job3 == 1 .OR. job3 == 3) then
     pp=matmul(transpose(qlat),p)
     call shortn3_plat(pp)
     p1 = matmul(plat,pp+nlatout(:,1))
     !call shortn(p,p1,dlv,nD)
  else
     call dcopy(3,p,1,p1,1)
  endif

  pi = 4*datan(1d0)
  tpi = 2*pi
  ! ... lx(ie) = lxi(ie) + lmxa
  lmax = -1
  if (nxi > 20) call rx('hsmq: increase dim of lx')
  do  5  ie = 1, nxi
     lx(ie) = lxi(ie)+lmxa
     lmax = max(lmax,lx(ie))
5 enddo
  ! ... If all rsm are equal to 1/a, no real-space part
  a = awald
  nDx = nD
  if (nG > 0 .AND. a /= 0) then
     if (rsm(1) > faca/a) then
        ltmp = .true.
        do  6  ie = 1, nxi
           ltmp = ltmp .and. abs(rsm(ie)-rsm(1)) .lt. 1d-12
6       enddo
        ltmp = ltmp .or. job2 .eq. 1
        if (ltmp) then
           a = 1/rsm(1)
           nDx = 1
        endif
     endif
  endif

  if ((lmax+1)**2 > nlmx) call rx('hsmq: lxi gt ll(nlmx)')
  if ((lmax+1)**2 > nlmxx) call rx('hsmq: lxi gt ll(nlmxx)')
  if (nrx < max(nD,nG)) call rx('hsmq: nrx < nD or nG')
  y0 = 1/dsqrt(4*pi)
  a2 = a*a
  call dpzero(hsm, 2*nlmx*nxi)
  call dpzero(hsmp,2*nlmx*nxi)

  ! --- Energy-independent setup for Q-space part ---
  if (nG > 0) then
     tpiba = 2*pi/alat
     gam = 0.25d0/a2
     do  10  ir = 1, nG
        wk(ir,2) = tpiba*(q(1) + qlv(1,ir))
        wk(ir,3) = tpiba*(q(2) + qlv(2,ir))
        wk(ir,4) = tpiba*(q(3) + qlv(3,ir))
10   enddo
     call ropyln(nG,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
     do  12  ir = 1, nG
        sp = alat*(wk(ir,2)*p1(1) + wk(ir,3)*p1(2) + wk(ir,4)*p1(3))
        xxc = cdexp(dcmplx(-gam*wk(ir,1), sp))
        wk(ir,2) = -dble(xxc)
        wk(ir,3) = -dimag(xxc)
12   enddo

     !   ... Q-space part of reduced strx for all energies
     call pvhsmq(0,nxi,lmxa,lx,exi,a,vol,nG,wk,nrx,yl,nlmx, &
          wk(1,2),wk(1,3),wk(1,4),wk(1,5),wk(1,6),wk(1,7),hsm,hsmp)
  endif

  ! --- Energy-independent setup for direct-space part ---
  if (job0 == 0) then
     !   ... Put ylm in yl and alat**2*(p-dlv)**2 in wk(1)
     do  20  ir = 1, nD
        wk(ir,2) = alat*(p1(1)-dlv(1,ir))
        wk(ir,3) = alat*(p1(2)-dlv(2,ir))
        wk(ir,4) = alat*(p1(3)-dlv(3,ir))
20   enddo
     call ropyln(nD,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
     !   ... cos(q.dlv), sin(q.dlv) -> wk(3,4), Y0 exp(-(a dlv)**2) -> wk(2)
     if (dabs(q(1))+dabs(q(2))+dabs(q(3)) /= 0) then
        do  22  ir = 1, nD
           qdotr = tpi*(q(1)*dlv(1,ir)+ q(2)*dlv(2,ir)+ q(3)*dlv(3,ir))
           wk(ir,3) = dcos(qdotr)
           wk(ir,4) = dsin(qdotr)
22      enddo
     else
        do  23  ir = 1, nD
           wk(ir,3) = 1
23      enddo
        call dpzero(wk(1,4),nD)
     endif
  endif
  ! ... If we are doing Ewald sums, we need this to make H(1/a) (hansr4)
  if (nG /= 0) then
     do  24  ir = 1, nDx
        wk(ir,6) = y0*dexp(-wk(ir,1)*a2)
24   enddo
  endif

  ! --- D-space part of reduced strx for each energy: chi(l)=wk(l+lc)---
  xx0 = 0
  do  30  ie = 1, nxi
     if (lx(ie) < lmxa) goto 30
     ir1 = 1
     if (job2 /= 0) then
        rsmi = rsm(1)
     else
        rsmi = rsm(ie)
     endif
     !   ... Case connecting vector p=0
     if (dcmpre(wk(1,1),0d0)) then
        if (job1 == 1) ir1 = 2
        if (job1 == 0 .AND. dcmpre(rsmi,0d0)) &
             call rx('hsmq: 10s digit job=0 and rsm=0 not allowed')
     endif
     if (job1 > 1) ir1 = 2
     if (lx(ie) > lmax) call rx('hsmq: lxi > lmax')
     lc = 7
     ls = lx(ie)+9

     !  ...  Check whether to omit r.s. part
     nDx = nD
     if (nG > 0) then
        if (abs(rsm(ie)-1/a) < 1d-12) nDx = 1
     endif
     !        if (ndx .eq. 1) then
     !          print *, 'hi',rsm,1/awald
     !        endif
     !   ... chi = H(rsm) - H(1/a) = (H(0) - H(1/a)) - (H(0) - H(rsm))
     if (nG > 0) then
        call hansr4(wk,lx(ie),nrx,nDx,exi(ie),1/a,wk(1,6), &
             wk(1,5),wk(1,lc))
     else
        do    l = 0, lx(ie)+1
           do    ir = ir1, nDx
              wk(ir,l+lc) = 0
           enddo
        enddo
     endif
     if ( .NOT. dcmpre(rsmi,0d0)) then
        xx = 1/rsmi**2
        !     ... Remake Y0 exp(-(r/rsm)**2) if rsm has changed
        if (xx /= xx0) then
           xx0 = xx
           do  33  ir = 1, nDx
              wk(ir,2) = y0*dexp(-wk(ir,1)*xx)
33         enddo
        endif
        call hansr4(wk,lx(ie),nrx,nDx,exi(ie),rsmi,wk(1,2), &
             wk(1,5),wk(1,ls))
        !         call prm('H(0) - H(rsm)',wk(1,ls),nrx,nDx,lx(ie)+2)
        do    l = 0, lx(ie)+1
           do    ir = ir1, nDx
              wk(ir,l+lc) = wk(ir,l+lc) - wk(ir,l+ls)
           enddo
        enddo
     endif

     !  --- Special treatment of on-site term ---
     !      Subtract Ewald contribution from G-vectors, hsm(p,a).
     !      Already generated by hansr4:
     !      Case p eq 0, wk(1,lc+1) = -hsm(p->0,a)
     !                   wk(1,lc)   = -2*hsmp(p->0,a)
     !      Case p ne 0, wk(1,lc+1) = hsm(p,rsm->0) - hsm(p,a)
     !                   wk(1,lc)   = 2*hsmp(p,rsm->0) - 2*hsmp(p,a)
     !   ... For p>0,job1=2, convert w(1,lc..lc+1) to -hsm and -2*hsmp
     if (ir1 == 2 .AND. .NOT. dcmpre(wk(1,1),0d0)) then
        !     ... Subtract 2*h0d from wk(1,lc), making it into -2*hsmp(p,a)
        !         And h0 from wk(1,lc+1), making it into -hsm(p,a)
        r = dsqrt(wk(1,1))
        if (exi(ie) > 0) then
           call rx(' hsmq: unfinished job1=2 implementation for e>0')
           kappa = dsqrt(exi(ie))
           h0 = cos(kappa*r)/r
           h0d = sin(kappa*r)/r
           wk(1,lc+1) = wk(1,lc+1) - h0
           wk(2,lc+1) = wk(2,lc+1) - h0d
        else
           akap = dsqrt(-exi(ie))
           h0 = exp(-akap*r)/r
           h0d = h0*r/akap
           wk(1,lc+1) = wk(1,lc+1) - h0
           wk(1,lc)   = wk(1,lc) - h0d
        endif
        !     ... Generate -hsm(p,a,l>0) by upward recursion
        xx = 4*a*y0*dexp(-(akap/a/2)**2 - (r*a)**2)
        do  31  l = 1, lx(ie)
           wk(1,l+lc+1) = ((2*l-1)*wk(1,l+lc) -exi(ie)*wk(1,l+lc-1) + xx) /wk(1,1)
           !          print *, wk(1,l+lc+1)
           xx = 2*a**2*xx
31      enddo
        !     ... debugging check
        !         call hansmr(r,exi(ie),a,xx1,lx(ie))
        !         call dpcopy(p1,xx3,1,3,alat)
        !         call hsm_mol(xx3,1/a,exi(ie),lmax,xx1,xx2)
        !         stop
     endif
     !   ... debugging check
     !       call hansmr(0d0,exi(ie),a,xx1,lx(ie))

     !   ... Make sin(qR)*(H(rsm,r)-H(1/a,r)), cos(qR)*(H(rsm,r)-H(1/a,r))
     do  34  l = 0, lx(ie)+1
        do  35  ir = 1, nDx
           wk(ir,l+ls) = wk(ir,4)*wk(ir,l+lc)
35      enddo
        do  36  ir = 1, nDx
           wk(ir,l+lc) = wk(ir,3)*wk(ir,l+lc)
36      enddo
34   enddo
     do  32  l = 0, lx(ie)
        lm = l*l
        !     ... xx1..4 artifically m-dependent to allow unrolling of m loop
        do  38  m = 1, 2*l+1
           ilm = lm+m
           xx1(m) = 0
           xx2(m) = 0
           xx3(m) = 0
           xx4(m) = 0
           do  37  ir = 1, nDx
              xx1(m) = xx1(m) + yl(ir,ilm)*wk(ir,l+lc+1)
              xx2(m) = xx2(m) + yl(ir,ilm)*wk(ir,l+ls+1)
              xx3(m) = xx3(m) + yl(ir,ilm)*wk(ir,l+lc)
              xx4(m) = xx4(m) + yl(ir,ilm)*wk(ir,l+ls)
37         enddo
           hsm(1,ilm,ie)  = (hsm(1,ilm,ie) + xx1(m))
           hsm(2,ilm,ie)  = (hsm(2,ilm,ie) + xx2(m))
           hsmp(1,ilm,ie) = (hsmp(1,ilm,ie) + xx3(m)/2)
           hsmp(2,ilm,ie) = (hsmp(2,ilm,ie) + xx4(m)/2)
38      enddo
32   enddo

     !   ... Put in phase to undo shortening, or different phase convention
     sp = tpi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
     if (job3 >= 2) sp = sp-tpi*(q(1)*p1(1)+q(2)*p1(2)+q(3)*p1(3))
     if (sp /= 0d0) then
        cosp = dcos(sp)
        sinp = dsin(sp)
        do  40  ilm = 1, (lx(ie)+1)**2
           x1 = hsm(1,ilm,ie)
           x2 = hsm(2,ilm,ie)
           hsm(1,ilm,ie) = x1*cosp - x2*sinp
           hsm(2,ilm,ie) = x2*cosp + x1*sinp
           x1 = hsmp(1,ilm,ie)
           x2 = hsmp(2,ilm,ie)
           hsmp(1,ilm,ie) = x1*cosp - x2*sinp
           hsmp(2,ilm,ie) = x2*cosp + x1*sinp
40      enddo
     endif
30 enddo
end subroutine hsmq

subroutine hsmqe0(lmax,rsm,job,q,p,nrx,nlmx,wk,yl, &
     awald,alat,qlv,nG,dlv,nD,vol,hsm)
  use m_ropyln,only: ropyln
  use m_shortn3_plat,only: shortn3_plat,nout,nlatout,plat=>platx,qlat=>qlatx
!  use m_lattic,only: plat=>lat_plat,qlat=>lat_qlat
  !- Bloch-sum of smooth Hankel functions for energy 0
  ! ---------------------------------------------------------------
  !i Inputs:
  !i  rsm,lmax:smoothing radius, l cutoff for hsm
  !i  dlv,nD:  direct lattice vectors, and number
  !i  qlv,nG:  reciprocal lattice vectors, and number
  !i  a,vol:   Ewald parameter and cell volume (used for Q-space sum)
  !i  job:     ones digit intended for H_smoothed - H_unsmoothed
  !i           To do, set nG=0, and the one's digit as follows:
  !i            0, generate wk(1,2,3,4)
  !i           >0  assume wk(1,3,4) and yl have been generated already
  !i           tens digit handles whether or not to add H(p=0).
  !i            0: return whole Bloch function (not allowed for rsm=0)
  !i            1: return Bloch function less H(p), if p=0
  !i           >1: return Bloch function less H(p), for any p
  !i           100s digit
  !i           >0: ignore passed q, use q=0
  !i           1000s digit
  !i            0 use standard phase convention, do not shorten p
  !i            1 use standard phase convention, but shorten p
  !i            2 scale standard phase convention by exp(-i q . p)
  !i            3 like 2, but also shorten p
  !i  q,p:     k-point and connecting vector (units of 2*pi/alat and alat)
  !i  nrx:     leading dimension of wk,yl
  !i  nlmx:    leading dimension of hsm
  !i  yl:      work array, dimensioned nrx*(lmax+1)**2
  !i  wk:      work array, dimensioned at least nrx*(2*lmax+10)
  !o Outputs:
  !o   yl:     ylm(1..nD,1..(lmax+1)**2) for points alat*(p-dlv(1..nD))
  !o   wk:     (*,1) holds r**2
  !o           (*,2) holds Y0 exp(-(r/rsm)**2) (rsm>0 only)
  !o           (*,3) holds cos(q.dlv)
  !o           (*,4) holds sin(q.dlv)
  !o   hsm:    smoothed Hankels for e=0,q=0 for l=0..lmax
  !r Remarks
  !r   hsmqe0 is an adaptation of hsmq for for e=0, and possibly q=0.
  !r   For l=0, hsm(q=0,e->0) diverges.  To achieve a finite hsm(l=0,e=0),
  !r   the average lim(e->0) hsm0(e) = -sqrt(4*pi)*exp(gamma*e) / (vol*e)
  !r   is subtracted from hsm(l=0).
  !b Bugs and Limitations:
  !b   see hsmq.
  !e External routines required: hansr5, ropyln, dpzero
  !u Updates
  !u  25 Jun 03 (Kino) call pvhsmq with arrays for lmax, e=0
  !u   7 May 98 handles case rsm>1/a by strict q-space summation
  ! ---------------------------------------------------------------
  !     implicit none
  integer :: lmax,nG,nD,nlmx,nrx,job
  double precision :: alat,awald,vol,rsm,p(3),q(3), &
       wk(nrx,*),yl(nrx,*),qlv(3,*),dlv(3,*), hsm(2,nlmx)
  ! Local variables
  integer :: ir,ilm,l,m,ir1,lc,ls,job0,job1,job2,job3,lm,nlmxx,nDx, &
       lx1(1)
  parameter (nlmxx=(200+1)**2) !(nlmxx=(16+1)**2)
  double precision :: qdotr,y0,a2,pi,sp,gam,tpiba,tpi,p1(3),ex1(1), &
       x1,x2,xx,xx1(nlmxx),xx2(nlmxx),r,h0,q0(3),a,faca,cosp,sinp
  parameter (faca=1d0)
  double complex xxc
  logical :: dcmpre,lqzero
  real(8):: pp(3)
  dcmpre(x1,x2) = dabs(x1-x2) .lt. 1d-8

  ! --- Setup ---
  job0 = mod(job,10)
  if (nG /= 0) job0 = 0
  job1 = mod(job/10,10)
  job2 = mod(job/100,10)
  job3 = mod(job/1000,10)
  if ((lmax+1)**2 > nlmx) call rx('hsmqe0: lmax gt ll(nlmx)')
  if ((lmax+1)**2 > nlmxx) call rx('hsmqe0: lmax gt ll(nlmxx)')
  if (nrx < max(nD,nG)) call rx('hsmqe0: nrx < nD or nG')
  ! ... Shorten connecting vector; adjust phase later
  if (job3 == 1 .OR. job3 == 3) then
     pp=matmul(transpose(qlat),p)
     call shortn3_plat(pp)
     p1 = matmul(plat,pp+nlatout(:,1))
     !call shortn(p,p1,dlv,nD)
  else
     call dcopy(3,p,1,p1,1)
  endif
  pi = 4*datan(1d0)
  y0 = 1/dsqrt(4*pi)
  tpi = 2*pi
  call dpzero(hsm, 2*nlmx)
  if (job2 > 0) then
     call dpzero(q0,3)
  else
     call dpcopy(q,q0,1,3,1d0)
  endif
  lqzero = dcmpre(dabs(q0(1))+dabs(q0(2))+dabs(q0(3)),0d0)
  ! ... If rsm ge 1/a, set a = 1/rsm and skip r.s. part
  a = awald
  nDx = nD
  if (nG > 0) then
     if (rsm > faca/a) nDx = 1
     if (nDx == 1) then
        a = 1/rsm
     endif
  endif
  a2 = a*a

  ! --- Setup for Q-space part ---
  if (nG > 0) then
     tpiba = 2*pi/alat
     gam = 0.25d0/a2
     do  10  ir = 1, nG
        wk(ir,2) = tpiba*(q0(1) + qlv(1,ir))
        wk(ir,3) = tpiba*(q0(2) + qlv(2,ir))
        wk(ir,4) = tpiba*(q0(3) + qlv(3,ir))
10   enddo
     call ropyln(nG,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
     do  12  ir = 1, nG
        sp = alat*(wk(ir,2)*p1(1) + wk(ir,3)*p1(2) + wk(ir,4)*p1(3))
        xxc = cdexp(dcmplx(-gam*wk(ir,1), sp))
        wk(ir,2) = -dble(xxc)
        wk(ir,3) = -dimag(xxc)
12   enddo

     !   --- Q-space part hsm(1/a,e,r) of reduced strx for all energies ---
     !   ... bug fix: copy scalars to arrays for call (Kino)
     lx1(1) = lmax
     ex1(1) = 0d0
     call pvhsmq(1,1,lmax,lx1,ex1,a,vol,nG,wk,nrx,yl,nlmx, &
          wk(1,2),wk(1,3),wk(1,4),wk(1,5),wk(1,6),wk(1,7),hsm,hsm)
  endif

  ! --- Energy-independent setup for direct-space part ---
  if (job0 == 0) then
     !   ... Put ylm in yl and alat**2*(p-dlv)**2 in wk(1)
     do  20  ir = 1, nDx
        wk(ir,2) = alat*(p1(1)-dlv(1,ir))
        wk(ir,3) = alat*(p1(2)-dlv(2,ir))
        wk(ir,4) = alat*(p1(3)-dlv(3,ir))
20   enddo
     call ropyln(nDx,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
     !   ... cos(q.dlv), sin(q.dlv) -> wk(3,4), Y0 exp(-(a dlv)**2) -> wk(2)
     if ( .NOT. lqzero) then
        do  22  ir = 1, nDx
           qdotr =tpi*(q0(1)*dlv(1,ir)+q0(2)*dlv(2,ir)+q0(3)*dlv(3,ir))
           wk(ir,3) = dcos(qdotr)
           wk(ir,4) = dsin(qdotr)
22      enddo
     else
        do  23  ir = 1, nDx
           wk(ir,3) = 1
23      enddo
        call dpzero(wk(1,4),nDx)
     endif
  endif
  ! ... If we are doing Ewald sums, we need this to make H(1/a) (hansr5)
  if (nG /= 0) then
     do  24  ir = 1, nDx
        wk(ir,6) = y0*dexp(-wk(ir,1)*a2)
24   enddo
  endif

  ! --- D-space part of reduced strx for each energy: chi(l)=wk(l+lc)---
  ir1 = 1
  ! ... Case connecting vector p=0
  if (dcmpre(wk(1,1),0d0)) then
     if (job1 == 1) ir1 = 2
     if (job1 == 0 .AND. dcmpre(rsm,0d0)) &
          call rx('hsmqe0: job=0 and rsm=0 not allowed')
  endif
  if (job1 > 1) ir1 = 2
  lc = 7
  ls = lmax+8

  ! ... chi = H(rsm) - H(1/a) = (H(0) - H(1/a)) - (H(0) - H(rsm))
  if (nG > 0) then
     call hansr5(wk,lmax,nrx,nDx,1/a,wk(1,6),wk(1,5),wk(1,lc))
  else
     do   l = 0, lmax+1
        do   ir = ir1, nDx
           wk(ir,l+lc) = 0
        enddo
     enddo
  endif
  if ( .NOT. dcmpre(rsm,0d0)) then
     xx = 1/rsm**2
     !   ... Make Y0 exp(-(r/rsm)**2)
     do  33  ir = 1, nDx
        wk(ir,2) = y0*dexp(-wk(ir,1)*xx)
33   enddo
     call hansr5(wk,lmax,nrx,nDx,rsm,wk(1,2),wk(1,5),wk(1,ls))
     !       call prm('H(0) - H(rsm)',wk(1,ls),nrx,nDx,lmax+1)
     do    l = 0, lmax
        do    ir = ir1, nDx
           wk(ir,l+lc) = wk(ir,l+lc) - wk(ir,l+ls)
        enddo
     enddo
  endif

  ! --- Special treatment of on-site term ---
  !     Subtract Ewald contribution from G-vectors, hsm(p,a).
  !     Already generated by hansr5:
  !     Case p eq 0, wk(1,lc) = -hsm(p->0,a)
  !     Case p ne 0, wk(1,lc) =  hsm(p,rsm->0) - hsm(p,a)
  ! ... For p>0,job=2, convert w(1,lc) to -hsm
  if (ir1 == 2 .AND. .NOT. dcmpre(wk(1,1),0d0)) then
     !   ... Subtract h0(0..1) from wk(1,lc..lc+1), making it into -hsm(p,a)
     r = dsqrt(wk(1,1))
     h0 = 1/r
     wk(1,lc) = wk(1,lc) - h0
     if (lmax > 0) wk(1,lc+1) = wk(1,lc+1) - h0/wk(1,1)
     !   ... Generate -hsm(p,a,l>0) by upward recursion
     xx = (2*a**2)*4*a*y0*dexp(-(r*a)**2)
     !  ... debugging check
     !      call hansmr(r,0d0,a,xx1,lmax)
     do  31  l = 2, lmax
        wk(1,l+lc) = ((2*l-1)*wk(1,l+lc-1) + xx)/wk(1,1)
        !       print *, wk(1,l+lc)
        xx = 2*a**2*xx
31   enddo
  endif

  ! ... Make sin(qR)*(H(rsm,r)-H(1/a,r)), cos(qR)*(H(rsm,r)-H(1/a,r))
  if ( .NOT. lqzero) then
     do  32  l = 0, lmax
        do  35  ir = 1, nDx
           wk(ir,l+ls) = wk(ir,4)*wk(ir,l+lc)
35      enddo
        do  34  ir = 1, nDx
           wk(ir,l+lc) = wk(ir,3)*wk(ir,l+lc)
34      enddo
32   enddo
  endif
  sp = tpi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
  if (job3 >= 2) sp = sp-tpi*(q(1)*p1(1)+q(2)*p1(2)+q(3)*p1(3))
  cosp = dcos(sp)
  sinp = dsin(sp)
  do  36  l = 0, lmax
     lm = l*l
     !   ... xx1..4 artifically m-dependent to allow unrolling of m loop
     if (lqzero) then
        do  138  m = 1, 2*l+1
           ilm = lm+m
           xx1(m) = 0
           xx2(m) = 0
           do  137  ir = 1, nDx
              xx1(m) = xx1(m) + yl(ir,ilm)*wk(ir,l+lc)
137        enddo
           hsm(1,ilm)  = (hsm(1,ilm) + xx1(m))
           hsm(2,ilm)  = 0d0
138     enddo
     else
        do  38  m = 1, 2*l+1
           ilm = lm+m
           xx1(m) = 0
           xx2(m) = 0
           do  37  ir = 1, nDx
              xx1(m) = xx1(m) + yl(ir,ilm)*wk(ir,l+lc)
              xx2(m) = xx2(m) + yl(ir,ilm)*wk(ir,l+ls)
37         enddo
           hsm(1,ilm)  = (hsm(1,ilm) + xx1(m))
           hsm(2,ilm)  = (hsm(2,ilm) + xx2(m))
38      enddo
     endif

     !   ... Put in phase to undo shortening
     if (sp /= 0d0) then
        do  40  m = 1, 2*l+1
           ilm = lm+m
           x1 = hsm(1,ilm)
           x2 = hsm(2,ilm)
           hsm(1,ilm) = x1*cosp - x2*sinp
           hsm(2,ilm) = x2*cosp + x1*sinp
40      enddo
     endif

36 enddo

  ! --- Add extra term for l=0 when e=0 and q=0 ---
  if (lqzero) hsm(1,1) = hsm(1,1) + rsm**2/(4*vol*y0)

end subroutine hsmqe0
subroutine pvhsmq(job,nx,lmxa,lxi,exi,a,vol,n,qsq,nyl,yl,nlmx, &
     cs,sn,cs2,sn2,cs3,sn3,hl,hlp)
  !- Makes H_L(1/a,e,r) by summation over reciprocal lattice vectors.
  !  job=1 : make hl but not hlp
  !  Skips any ie for which lxi(ie)<lmxa.
  !  This is a kernel called by hsmq.
  !     implicit none
  integer :: n,nx,lmxa,lxi(nx),nlmx,job,nyl
  double precision :: cof,e,gam,pi,vol
  double precision :: cs(n),sn(n),cs2(n),sn2(n),cs3(n),sn3(n), &
       a,exi(1),qsq(n),yl(nyl,*),hl(2,nlmx,*),hlp(2,nlmx,*) !kino
  ! Local variables
  integer :: i,i1,ie,ilm,l,lx,m,je,lm,nlmxx
  parameter (nlmxx=(200+1)**2)! (nlmxx=(16+1)**2)
  double precision :: c,s,xx,x1,x2
  double precision :: c2(nlmxx),s2(nlmxx),c3(nlmxx),s3(nlmxx)
  logical :: dcmpre
  dcmpre(x1,x2) = dabs(x1-x2) .lt. 1d-8
  gam = 0.25d0/a**2
  pi = 4*datan(1d0)
  ! --- For each exi ne 0 and ilm do ---
  do  20  ie = 1, nx
     e = exi(ie)
     lx = lxi(ie)
     if (lxi(ie) < lmxa) goto 20

     !   ... Copy existing hl,hlp if already calculated for this exi
     do  30  je = 1, ie-1
        if (dabs(exi(je)-exi(ie)) > 1d-8) goto 30
        if (lxi(ie) > lxi(je)) goto 30
        hl (1:2,1:(lxi(ie)+1)**2,ie)= hl (1:2,1:(lxi(ie)+1)**2,je)
        hlp(1:2,1:(lxi(ie)+1)**2,ie)= hlp(1:2,1:(lxi(ie)+1)**2,je)
        goto 20
30   enddo

     i1 = 1
     if (dcmpre(e-qsq(1),0d0)) i1 = 2
     do  22  i = i1, n
        xx = 1/(e-qsq(i))
        c = cs(i)*xx
        cs2(i) = c
        cs3(i) = c*xx
        s = sn(i)*xx
        sn2(i) = s
        sn3(i) = s*xx
22   enddo
     cof = -4d0*pi*dexp(gam*e)/vol
     ilm = 0
     do  21  l = 0, lx, 2
        cof = -cof
        lm = l*l
        !     ... Do the even and odd l's separately
        !     ... NB c2,s2,c3,s3 m-dependent to allow m loop to be unrolled
        if (job == 0) then
           do  24  m = 1, 2*l+1
              !           ilm = ilm+1
              ilm = lm+m
              c2(ilm) = 0
              s2(ilm) = 0
              c3(ilm) = 0
              s3(ilm) = 0
              do  23  i = i1, n
                 c2(ilm) = c2(ilm) + yl(i,ilm)*cs2(i)
                 c3(ilm) = c3(ilm) + yl(i,ilm)*cs3(i)
                 s2(ilm) = s2(ilm) + yl(i,ilm)*sn2(i)
                 s3(ilm) = s3(ilm) + yl(i,ilm)*sn3(i)
23            enddo
              hl(1,ilm,ie)  = hl(1,ilm,ie)  + cof*c2(ilm)
              hlp(1,ilm,ie) = hlp(1,ilm,ie) + cof*(gam*c2(ilm)-c3(ilm))
              hl(2,ilm,ie)  = hl(2,ilm,ie)  + cof*s2(ilm)
              hlp(2,ilm,ie) = hlp(2,ilm,ie) + cof*(gam*s2(ilm)-s3(ilm))
24         enddo
        else
           do  124  m = 1, 2*l+1
              !           ilm = ilm+1
              ilm = lm+m
              c2(ilm) = 0
              s2(ilm) = 0
              do  123  i = i1, n
                 c2(ilm) = c2(ilm) + yl(i,ilm)*cs2(i)
                 s2(ilm) = s2(ilm) + yl(i,ilm)*sn2(i)
123           enddo
              hl(1,ilm,ie)  = hl(1,ilm,ie)  + cof*c2(ilm)
              hl(2,ilm,ie)  = hl(2,ilm,ie)  + cof*s2(ilm)
124        enddo
        endif
        if (l+1 <= lx) then
           lm = (l+1)*(l+1)
           if (job == 0) then
              do  26  m = 1, 2*l+3
                 !             ilm = ilm+1
                 ilm = lm+m
                 c2(ilm) = 0
                 s2(ilm) = 0
                 c3(ilm) = 0
                 s3(ilm) = 0
                 do  27  i = i1, n
                    c2(ilm) = c2(ilm) + yl(i,ilm)*cs2(i)
                    c3(ilm) = c3(ilm) + yl(i,ilm)*cs3(i)
                    s2(ilm) = s2(ilm) + yl(i,ilm)*sn2(i)
                    s3(ilm) = s3(ilm) + yl(i,ilm)*sn3(i)
27               enddo
                 hl(1,ilm,ie)  = hl(1,ilm,ie)  + cof*s2(ilm)
                 hlp(1,ilm,ie) = hlp(1,ilm,ie) + cof*(gam*s2(ilm)-s3(ilm))
                 hl(2,ilm,ie)  = hl(2,ilm,ie)  - cof*c2(ilm)
                 hlp(2,ilm,ie) = hlp(2,ilm,ie) - cof*(gam*c2(ilm)-c3(ilm))
26            enddo
           else
              do  126  m = 1, 2*l+3
                 !             ilm = ilm+1
                 ilm = lm+m
                 c2(ilm) = 0
                 s2(ilm) = 0
                 do  127  i = i1, n
                    c2(ilm) = c2(ilm) + yl(i,ilm)*cs2(i)
                    s2(ilm) = s2(ilm) + yl(i,ilm)*sn2(i)
127              enddo
                 hl(1,ilm,ie)  = hl(1,ilm,ie)  + cof*s2(ilm)
                 hl(2,ilm,ie)  = hl(2,ilm,ie)  - cof*c2(ilm)
126           enddo
           endif
        endif
21   enddo

     ! --- Add extra term for l=0 when e=0 and q=0 ---
     if (dcmpre(e-qsq(1),0d0)) hl(1,1,ie) = hl(1,1,ie) - &
          sqrt(4d0*pi)/vol*gam
20 enddo
end subroutine pvhsmq

subroutine hansr4(rsq,lmax,nrx,nr,e,rsm,wk,wk2,chi)
  !- Difference between true,smoothed hankels for l=-1...lmax
  !  for e nonzero.  Vectorizes for negative e.
  ! ---------------------------------------------------------------
  !i Inputs
  !i   rsq,nr:vector of points r**2, and number of points.
  !i          Only the first value of rsq may be zero (See Remarks)
  !i   lmax:  highest l for which to evaluate xi.
  !i   e,rsm: smoothing radius and energy
  !i   nrx:   leading dimension of chi
  !i   wk:    array containing y0*dexp(-(r/rsm)**2)
  !i   wk2:   a work array of length nr.
  !o Outputs:
  !o   chi:  Difference between smoothed and true Hankel for l=-1..lmax.
  !r Remarks
  !r   Smooth Hankels for e<=0 defined in J. Math. Phys. 39, 3393 (1998).
  !r   Notes on conventions:
  !r  *  akap^2 = -e with akap>0
  !r   JMP defines akap in contradistinction to usual convention for kappa:
  !r     kappa = sqrt(e), Im(kappa) >= 0.
  !r   It is related to kappa (defined according to usual conventions) as
  !r     akap = -i kappa (real and positive for e<0)
  !r  *JMP defines
  !r      u(+/-) = exp(-/+ akap*r) erfc(akap*rsm/2 -/+ r/rsm)
  !r             = exp(+/- i kappa r) erfc(-i*kappa*rsm/2 -/+ r/rsm)
  !r   Note: old codes (e.g. rcnsl,hsmbld) define
  !r      uplus = u-  and umins = u+
  !r   The smoothed Hankels for l=0,-1 are:
  !r      h^s_0 (r) = (u+ - u-) / 2r    = (umins - uplus) / 2r
  !r      h^s_-1(r) = (u+ + u-) / 2akap = (umins + uplus) / 2akap
  !r
  !r   We can a new set of functions
  !r      U+/- = 1/2 exp(+/- i kappa r) erfc(r/rsm +/- i*kappa*rsm/2)
  !r   It is convenient to use
  !r      erfc(-x^*) = 2-erfc^*(x).
  !r
  !r  *For e<0, i*kappa is real and U+/- are also real:
  !r      U+   = exp( i kappa r ) - u+/2
  !r      U-   = u-/2
  !r   The difference in smoothed, unsmoothed Hankels is
  !r      h_0  - h^s_0  = (exp(-akap r) - umins/2 + uplus/2) /r
  !r                    = (exp(i kappa r ) - u+/2 + u-/2) /r
  !r                    = (U+ + U-) /r
  !r      h_-1 - h^s_-1 = (exp(-akap r) - umins/2 - uplus/2) /akap
  !r                    = (exp(i kappa r ) - u+/2 - u-/2) /akap
  !r                    = (U+ - U-) /(-i kappa)
  !r
  !r  *For e>0, kappa is real; thus U+ = U-*  since
  !r     erfc(x*) = erfc*(x)
  !r   Keeping the same conventions in h-h^s for e>0 we get
  !r      h_0  - h^s_0  = (U+ + U-) /r = 2 Re (U+) /r
  !r      h_-1 - h^s_-1 = (U+ - U-) / (-i kappa) = -2 Im (U+) / kappa
  !r   The differences are thus real.
  !r
  !r  *The first point may have r->0.  In this case, chi(1,l)=0 except
  !r   chi(1,0), chi(1,-1) which are returned as the (-) value of the
  !r   smoothed hankels, i.e. the infinite unsmoothed part is discarded.
  !r
  !r   In the limit r->0
  !r     h_0  -> 1/r + i*kappa

  !r     U+/- -> 1/2 erfc(+/- i*kappa*rsm/2) = 1 - erfc(i*kappa*rsm/2)
  !r     U+   -> 1/2 erfc(i*kappa*rsm/2)
  !r     U-   -> 1/2 erfc(-i*kappa*rsm/2)
  !r          =  1 - 1/2 erfc(i*kappa*rsm/2) = 1 - U+
  !r     U+ + U- -> 1 + terms of order r . Calculate:
  !r     d/dr (erfc(r/rsm+/-i*kappa*rsm/2)) -> 2/srpi/rsm *
  !r                                           exp(-(i*kappa*rsm/2)^2)
  !r     dU+/- /dr |r=0 -> +/- i*kappa U+/-(r->0)
  !r                        - 1/srpi/rsm * exp(-(i*kappa*rsm/2)^2)
  !r     U+ + U- -> 1 + r*i*kappa( U+(0) - U-(0))
  !r                  - r/srpi/rsm * exp(-(i*kappa*rsm/2)^2))
  !r              = 1 + r*i*kappa( 2*U+(0) - 1)
  !r                  - r/srpi/rsm * exp(-(i*kappa*rsm/2)^2))
  !r     h^s_0  = h0 - U+ - U-
  !r           -> 1/r + i*kappa - (1/r + i*kappa( 2*U+(0) - 1))
  !r                  + 1/srpi/rsm * exp(-(i*kappa*rsm/2)^2))
  !r            = i*kappa*erfc(i*kappa*rsm/2) +
  !r              1/srpi/rsm * exp(-(i*kappa*rsm/2)^2))
  !r     h_-1  = exp(i kappa r) /(-i kappa) -> 1/(-i kappa)
  !r     h^s_-1 = 1/(-i kappa) - (U+(0) - U-(0))/(-i kappa)
  !r            = 1/(-i kappa) - (2*U+(0) - 1)/(-i kappa)
  !r            = 1/(-i kappa) - (1 - 2*U-(0))/(-i kappa)
  !r            = 2*U-(0)/(-i kappa)
  !r   For e<0   i*kappa = -akap (real)
  !r     h^s_0  = akap*erfc(akap*rsm/2) + 1/srpi/rsm*exp(-(akap*rsm/2)^2))
  !r     h^s_-1 = 2*U-(0) = erfc(akap*rsm/2)
  !r   For e>0   Use the expressions above
  !r             Warning! hs is NOT REAL ... returns only real part
  !r             h^s_-1 is discontinuous there:
  !r             As e->0, for e>0, 2*U-(0)/(-i kappa) -> i/kappa
  !r             As e->0, for e<0, 2*U-(0)/(-i kappa) -> 1/akap
  !r
  !r  Limiting behavior e<0:
  !r    If both (akap*rsm/2 -/+ r/rsm) >> 1, we have
  !r   -log u(+/-) -> (akap*rsm/2 -/+ r/rsm)^2 -/+ akap*r
  !r               =  (akap*rsm/2)^2 + (r/rsm)^2 >> 1
  !r    u(+/-)     -> exp[-(akap*rsm/2)^2 - (r/rsm)^2]
  !r                  is negligible compared to 1
  !r  Also, if akap*r >> 1,
  !r    h_0 -> exp(-akap*r)/r > h_s^0  is negligible compared to 1
  !u Updates
  !u  15 Aug 00 extended to e>0 (but erfc not vectorized)
  ! ---------------------------------------------------------------
  !     implicit none
  integer :: nrx,nr,lmax
  double precision :: rsq(1),e,chi(nrx,-1:*),rsm,wk(nr),wk2(nr)
  ! Local variables
  logical :: lpos
  integer :: l,ir,ir1
  real(8):: sre,r2,xx,ra,h0,arsm,earsm,kappa=0d0,ta=1d99, akap,a,r,um,up,x,facgl
  complex(8):: zikap=-99999d0,zerfc,eikr,zuplus

  ! ... erfc(x) is evaluated inline as a ratio of polynomials,
  !     to a relative precision of <10^-15 for x<5.
  !     Different polynomials are used for x<1.3 and x>1.3.
  !     Numerators and denominators are t,b respectively.
  double precision :: w,f1,f2, &
       t10,t11,t12,t13,t14,t15,t16,t17,b11,b12,b13,b14,b15,b16,b17,b18, &
       t20,t21,t22,t23,t24,t25,t26,t27,b21,b22,b23,b24,b25,b26,b27,b28
  parameter ( &
       t10=2.1825654430601881683921d0, t20=0.9053540999623491587309d0, &
       t11=3.2797163457851352620353d0, t21=1.3102485359407940304963d0, &
       t12=2.3678974393517268408614d0, t22=0.8466279145104747208234d0, &
       t13=1.0222913982946317204515d0, t23=0.3152433877065164584097d0, &
       t14=0.2817492708611548747612d0, t24=0.0729025653904144545406d0, &
       t15=0.0492163291970253213966d0, t25=0.0104619982582951874111d0, &
       t16=0.0050315073901668658074d0, t26=0.0008626481680894703936d0, &
       t17=0.0002319885125597910477d0, t27=0.0000315486913658202140d0, &
       b11=2.3353943034936909280688d0, b21=1.8653829878957091311190d0, &
       b12=2.4459635806045533260353d0, b22=1.5514862329833089585936d0, &
       b13=1.5026992116669133262175d0, b23=0.7521828681511442158359d0, &
       b14=0.5932558960613456039575d0, b24=0.2327321308351101798032d0, &
       b15=0.1544018948749476305338d0, b25=0.0471131656874722813102d0, &
       b16=0.0259246506506122312604d0, b26=0.0061015346650271900230d0, &
       b17=0.0025737049320207806669d0, b27=0.0004628727666611496482d0, &
       b18=0.0001159960791581844571d0, b28=0.0000157743458828120915d0)
  ! ... f1(w=x-1/2) is erfc(x) for x<1.3, if xx is y0*dexp(-x*x)
  f1(w) = xx*(((((((t17*w+t16)*w+t15)*w+t14)*w+t13)*w+t12)* &
       w+t11)*w+t10)/((((((((b18*w+b17)*w+b16)*w+b15)*w+b14)* &
       w+b13)*w+b12)*w+b11)*w+1)
  ! ... f2(w=x-2) is erfc(x) for x>1.3, if xx is y0*dexp(-x*x)
  f2(w) = xx*(((((((t27*w+t26)*w+t25)*w+t24)*w+t23)*w+t22)* &
       w+t21)*w+t20)/((((((((b28*w+b27)*w+b26)*w+b25)*w+b24)* &
       w+b23)*w+b22)*w+b21)*w+1)

  ! --- Setup ---
  if (lmax < 0 .OR. nr == 0) return
  a = 1/rsm
  lpos = e .gt. 0d0
  if (lpos) then
     kappa  =  dsqrt(e)
     zikap  = dcmplx(0d0,1d0)*kappa
     ta = 2*a
     !       facgl = 4*a*exp(e/(2*a)**2)
     facgl = 4*a*exp(e/ta**2)
  else
     akap = dsqrt(-e)
     arsm = akap*rsm/2
     earsm = dexp(-arsm**2)/2
     !       facgl = 4*a*exp(e/(2*a)**2)
     facgl = 8*a*earsm
  endif
  ! --- chi(*,-1), chi(*,0) = true - sm hankel for each r ---
  ! ... chi(r,-1) = (h0 - um - up)*r/akap   and
  !     chi(r,0)  = (h0 - um + up)          with
  !     h0 = exp(-akap*r)/r
  !     up = erfc(akap*rsm/2+r/rsm)*exp(akap*r)/(2r)
  !     um = erfc(akap*rsm/2-r/rsm)/exp(akap*r)/(2r)
  !     um,up correspond to umins/(2r),uplus/(2r) in hansmr
  ir1 = 1
  ! ... See Remarks for derivation of expressions in this case
  if (rsq(1) < 1d-12) then
     ir1 = 2
     if (lpos) then
        chi(1,0)  = zerfc(zikap/ta)*zikap-facgl/dsqrt(16d0*datan(1d0))
        chi(1,-1) = -dimag(zerfc(zikap/ta))/kappa
     else
        !     ... make h0 = erfc(arsm) = erfc(akap*rsm/2)
        xx = earsm/dsqrt(4d0*datan(1d0))
        if (arsm > 1.3d0) then
           h0 = f2(arsm-2d0)
        else
           h0 = f1(arsm-.5d0)
        endif
        !         chi(-1) -> erfc(akap/2a)/akap for r=0
        chi(1,0)  = akap*h0 - 4*a*xx
        chi(1,-1) = -h0/akap
     endif
  endif
  ! ... Make chi(r,rsm->0,l) - chi(r,rsm,l) for l=-1, 0
  if (lpos) then
     do  22  ir = ir1, nr
        r2 = rsq(ir)
        r = dsqrt(r2)
        ra = r*a
        eikr = exp(zikap*r)
        !         h_0  - h^s_0  =  Re (U+) /r;  h_-1 - h^s_-1 = -Im (U+) / kappa
        zuplus  = zerfc(r*a + zikap/ta)*eikr

        chi(ir,0)  =  dble(zuplus)/r
        chi(ir,-1) = -dimag(zuplus)/kappa
        wk2(ir) = facgl*wk(ir)
22   enddo
  else
     do  20  ir = ir1, nr
        r2 = rsq(ir)
        r = dsqrt(r2)
        ra = r*a
        sre = akap*r
        h0 = dexp(-sre)/r
        xx = earsm*wk(ir)/r
        !     ... Evaluate um; use (akap/2a-ra)^2 = (akap/2a)^2+(ra)^2-akap*r
        x = ra - arsm
        if (x > 1.3d0) then
           um = h0 - f2(x-2d0)
        elseif (x > 0) then
           um = h0 - f1(x-.5d0)
        elseif (x > -1.3d0) then
           um = f1(-x-.5d0)
        else
           um = f2(-x-2d0)
        endif
        !     ... Evaluation of up assumes x gt 0
        x = ra + arsm
        if (x > 1.3d0) then
           up = f2(x-2d0)
        else
           up = f1(x-.5d0)
        endif
        !     ... Make chi(r,rsm->0,l) - chi(r,rsm,l) for l=-1, 0
        chi(ir,-1) = (h0 - um - up)*r/akap
        chi(ir,0)  =  h0 - um + up
        wk2(ir) = facgl*wk(ir)
20   enddo
  endif
  ! --- chi(ir,l) for l>1 by upward recursion ---
  facgl = 2*a**2
  do  31  l = 1, lmax
     chi(1,l) = 0
31 enddo
  do  30  l = 1, lmax
     xx = 2*l-1
     do  7  ir = ir1, nr
        chi(ir,l) = (xx*chi(ir,l-1) - e*chi(ir,l-2) + wk2(ir))/rsq(ir)
        wk2(ir) = facgl*wk2(ir)
7    enddo
30 enddo
end subroutine hansr4

subroutine hansr5(rsq,lmax,nrx,nr,rsm,wk,wk2,chi)
  !- Difference between true,smoothed hankels for l=0...lmax, e=0
  ! ---------------------------------------------------------------
  !i Inputs
  !i   rsq,nr:vector of points r**2, and number of points.
  !i          Only the first value of rsq may be zero (See Remarks)
  !i   lmax:  highest l for which to evaluate xi.
  !i   rsm:   smoothing radius and energy
  !i   nrx:   leading dimension of chi
  !i   wk:    array containing y0*dexp(-(r/rsm)**2)
  !i   wk2:   a work array of length nr.
  !o Outputs:
  !o   chi:  Difference between smoothed and true Hankel for l=0..lmax.
  !r Remarks
  !r   The first point may have r->0.  In this case, chi(1,l)=0 except
  !r   chi(1,0), which is returned as the -hsm(e=0,r=0) i.e.
  !r   the infinite unsmoothed part is discarded.
  ! ---------------------------------------------------------------
  !     implicit none
  integer :: nrx,nr,lmax
  double precision :: rsq(1),chi(nrx,0:lmax),rsm,wk(nr),wk2(nr)
  double precision :: xx,ra,a,r,facgl
  integer :: l,ir,ir1
  ! ... erfc(x) is evaluated as a ratio of polynomials,
  !     to a relative precision of <10^-15 for x<5.
  !     Different polynomials are used for x<1.3 and x>1.3.
  !     Numerators and denominators are t,b respectively.
  double precision :: w,f1,f2, &
       t10,t11,t12,t13,t14,t15,t16,t17,b11,b12,b13,b14,b15,b16,b17,b18, &
       t20,t21,t22,t23,t24,t25,t26,t27,b21,b22,b23,b24,b25,b26,b27,b28
  parameter ( &
       t10=2.1825654430601881683921d0, t20=0.9053540999623491587309d0, &
       t11=3.2797163457851352620353d0, t21=1.3102485359407940304963d0, &
       t12=2.3678974393517268408614d0, t22=0.8466279145104747208234d0, &
       t13=1.0222913982946317204515d0, t23=0.3152433877065164584097d0, &
       t14=0.2817492708611548747612d0, t24=0.0729025653904144545406d0, &
       t15=0.0492163291970253213966d0, t25=0.0104619982582951874111d0, &
       t16=0.0050315073901668658074d0, t26=0.0008626481680894703936d0, &
       t17=0.0002319885125597910477d0, t27=0.0000315486913658202140d0, &
       b11=2.3353943034936909280688d0, b21=1.8653829878957091311190d0, &
       b12=2.4459635806045533260353d0, b22=1.5514862329833089585936d0, &
       b13=1.5026992116669133262175d0, b23=0.7521828681511442158359d0, &
       b14=0.5932558960613456039575d0, b24=0.2327321308351101798032d0, &
       b15=0.1544018948749476305338d0, b25=0.0471131656874722813102d0, &
       b16=0.0259246506506122312604d0, b26=0.0061015346650271900230d0, &
       b17=0.0025737049320207806669d0, b27=0.0004628727666611496482d0, &
       b18=0.0001159960791581844571d0, b28=0.0000157743458828120915d0)
  ! ... f1(w=x-1/2) is erfc(x) for x<1.3, if xx is y0*dexp(-x*x)
  f1(w) = xx*(((((((t17*w+t16)*w+t15)*w+t14)*w+t13)*w+t12)* &
       w+t11)*w+t10)/((((((((b18*w+b17)*w+b16)*w+b15)*w+b14)* &
       w+b13)*w+b12)*w+b11)*w+1)
  ! ... f2(w=x-2) is erfc(x) for x>1.3, if xx is y0*dexp(-x*x)
  f2(w) = xx*(((((((t27*w+t26)*w+t25)*w+t24)*w+t23)*w+t22)* &
       w+t21)*w+t20)/((((((((b28*w+b27)*w+b26)*w+b25)*w+b24)* &
       w+b23)*w+b22)*w+b21)*w+1)
  ! --- Setup ---
  if (lmax < 0 .OR. nr == 0) return
  a = 1/rsm
  facgl = 4*a
  ! --- chi(*0), chi(*,1) ---
  ir1 = 1
  if (rsq(1) < 1d-6) then
     ir1 = 2
     !   ... hsm(e=0,r=0) = 4*a*y0
     chi(1,0) = -2d0*a/dsqrt(4d0*datan(1d0))
  endif
  do  20  ir = ir1, nr
     r = dsqrt(rsq(ir))
     ra = r*a
     xx = wk(ir)/r
     if (ra > 1.3d0) then
        chi(ir,0) = f2(ra-2d0)
     else
        chi(ir,0) = f1(ra-.5d0)
     endif

     wk2(ir) = facgl*wk(ir)
20 enddo
  ! --- chi(ir,l) for l>0 by upward recursion ---
  facgl = 2*a**2
  do  31  l = 1, lmax
     chi(1,l) = 0
31 enddo
  do  30  l = 1, lmax
     xx = 2*l-1
     do  7  ir = ir1, nr
        chi(ir,l) = (xx*chi(ir,l-1) + wk2(ir))/rsq(ir)
        wk2(ir) = facgl*wk2(ir)
7    enddo
30 enddo
end subroutine hansr5

