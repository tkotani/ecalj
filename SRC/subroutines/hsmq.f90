!>Bloch sum of smooth Hankel functions. JMP39: https://doi.org/doi:10.1063/1.532437.
module m_hsmq !Bloch sum of smooth Hankel. lm expansion
  ! JMP39: Bott, E., M. Methfessel, W. Krabs, and P. C. Schmidt.
  ! “Nonsingular Hankel Functions as a New Basis for Electronic Structure Calculations.”
  ! Journal of Mathematical Physics 39, no. 6 (June 1, 1998): 3393–3425.
  ! https://doi.org/doi:10.1063/1.532437.
  public:: hsmq,hsmqe0
contains
  subroutine hsmq(nxi,lmxa,lxi,exi,rsm,job1,q,p,nrx,nlmx,yl,awald,alat,qlv,nG,dlv,nD,vol,hsm,hsmp) !Bloch-sum of smooth Hankel functions and energy derivatives at p by Ewald summation, for nxi groups of parameters (lxi,exi,rsm).
    use m_ropyln,only: ropyln
    use m_shortn3_plat,only: shortn3_plat,nout,nlatout,plat=>platx,qlat=>qlatx
    implicit none
    intent(in)::  nxi,lmxa,lxi,exi,rsm,job1,q,p,nrx,nlmx,   awald,alat,qlv,nG,dlv,nD,vol
    intent(out)::                                        yl,                             hsm,hsmp
    !i Inputs:
    !i  nxi,lxi,exi,rsm:  number of (l,energy,rsm) groups, and list
    !i  dlv,nD:  direct lattice vectors, and number
    !i  qlv,nG:  reciprocal lattice vectors, and number
    !i  awald,vol: Ewald parameter and cell volume (used for Q-space sum)
    !i  lmxa:    make strux for lxi(ie)+lmxa,
    !i           skipping over any energy for which lxi(ie) < 0
    !i  job1 only:  
    !i  q,p:     k-point and connecting vector (units of 2*pi/alat and alat)
    !i  nrx:     leading dimension of wk,yl; must be at least max(nD,nG)
    !i  nlmx:    leading dimension of hsm,hsmp
    !i  yl:      work array, dimensioned nrx*(lmax+1)**2
    !o Outputs:
    !o   hsm,hsmp:smoothed Hankels and energy derivatives for
    !o           nxi energies, to lmax lxi(1..nxi).

    !  yl:     ylm(1..nD,1..(lmax+1)**2) for points alat*(p-dlv(1..nD))
    !
    !     rsq  : holds r**2
    !     y0fac: Y0 exp(-(r/rsm(nxi))**2) (rsm>0 only)
    !     phase 

    !b Bugs and Limitations:
    !b   Convergence is poor for p small but nonzero.  To avoid this,
    !b   hsmq SHOULD re-evaluate points for which p is small using
    !b   a power series expansion.
    !u Updates
    !u  15 Aug 00 extended to e>0; added 1000 digit job
    !u   7 May 98 handles case rsm>1/a by strict q-space summation
    !u  16 Sep 98 fixed dimensioning error arising when nD > nG
    ! ---------------------------------------------------------------
    integer :: nxi,lmxa,lmax,nG,nD,nrx,lxi(nxi),job,nlmx, ie,ir,ilm,l,m,ir1,lc,ls,job0,job1,job2,job3,lm, lx(20),ndx,le,i1,i2
    real(8),parameter:: pi = 4d0*datan(1d0), tpi = 2d0*pi,y0 = 1/dsqrt(4d0*pi),faca=1d0
    real(8):: qdotr,a2,sp,gam,tpiba,rsmi, pp(3),x1,x2,xx,xx0,r,akap=1d99,kappa,h0,h0d,a,p1(3),cosp,sinp,&
         alat,awald,vol,exi(nxi),rsm(nxi),p(3),q(3), yl(nrx,*),qlv(3,*),dlv(3,*),y0fac(nrx),rsq(nrx)
    complex(8):: hsmp(nlmx,nxi), hsm(nlmx,nxi),xxc,phase(nrx),img=(0d0,1d0) !,phasex
    if(nrx < max(nD,nG)) call rx('hsmq: nrx < nD or nG')
    if(nxi > 20)         call rx('hsmq: increase dim of lx')
    p1=p
    lx(1:nxi) = lxi(1:nxi)+lmxa ! ... lx(ie) = lxi(ie) + lmxa
    lmax = maxval(lx(1:nxi))
    a = awald ! ... If all rsm are equal to 1/a, no real-space part
    nDx = nD
    hsm=0d0  
    hsmp=0d0
    if(nG>0 .AND. a/=0d0 .and. rsm(1)>faca/a .and. all(abs(rsm(1:nxi)-rsm(1))<1d-12) ) then
       a = 1/rsm(1)
       nDx = 1
    endif
    a2 = a*a
    qqblock: block ! --- Energy-independent setup for Q-space part ---
      real(8):: qq(nG,3)
      complex(8):: xxcc(nG)
      if (nG > 0) then 
         tpiba = 2d0*pi/alat
         gam = 0.25d0/a2
         do  ir = 1, nG
            qq(ir,:) = tpiba*(q(:) + qlv(:,ir))
         enddo
         call ropyln(nG,qq(1,1),qq(1,2),qq(1,3),lmax,nrx,yl,rsq)
         xxcc = -[(exp(-gam*rsq(ir)+ img*alat*sum(qq(ir,:)*p1) ),ir=1,nG)]
         !   ... Q-space part of reduced strx for all energies
         call pvhsmq(0,nxi,lmxa,lx,exi,a,vol,nG,rsq,nrx,yl,(lmax+1)**2,xxcc, hsm,hsmp)
      endif
    endblock qqblock
    pplblock: block ! --- Energy-independent setup for Real-space part ---
      real(8):: pp(nD,3)
      do  ir = 1, nD 
         pp(ir,:) = alat*(p1(:)-dlv(:,ir))
      enddo
      call ropyln(nD,pp(1,1),pp(1,2),pp(1,3),lmax,nrx,yl,rsq) !rsq=alat**2*(p-dlv)**2
    endblock pplblock
    if (nG /= 0) y0fac(1:nDx) = y0*dexp(-rsq(1:nDx)*a2) !we need this to make H(1/a) (hansr4)
    phase(1:nDx)=[(exp(img*tpi*sum(q*dlv(:,ir))), ir=1,nDx)]
    ! ...  D-space part of reduced strx for each energy: chi(l)=wk(l+lc)---
    xx0 = 0
    do  30  ie = 1, nxi
       if (lx(ie) < lmxa) cycle
       hansblock: block
         real(8):: wkc(nDx,0:lx(ie)+1),wks(nDx,0:lx(ie)+1)  
         ir1 = 1
         rsmi = rsm(ie)   !rsmi=0 for hvcc
         if (dcmpre(rsq(1))) then!   ... Case connecting vector p=0
            if(job1==1) ir1 = 2
            if(job1==0 .AND. dcmpre(rsmi))call rx('hsmq: 10s digit job=0 and rsm=0 not allowed')
         endif
         if (lx(ie) > lmax) call rx('hsmq: lxi > lmax')
         nDx = nD
         !   ... chi = H(rsm) - H(1/a) = (H(0) - H(1/a)) - (H(0) - H(rsm))
         if(nG>0) then
            if(abs(rsm(ie)-1/a)<1d-12) nDx = 1 !  ...  Check whether to omit r.s. part
            call hansr4(rsq,lx(ie),nDx,nDx,exi(ie),1/a,y0fac,wkc) !H(1/a)
         else
            wkc(ir1:nDx,:) = 0d0
         endif
         if ( .NOT. dcmpre(rsmi)) then
            xx = 1/rsmi**2
            y0fac(1:nDx) = y0*dexp(-rsq(1:nDx)*xx)! Y0 exp(-(r/rsm)**2) if rsm has changed 
            call hansr4(rsq,lx(ie),nDx,nDx,exi(ie),rsmi,y0fac,wks)
            wkc(ir1:nDx,:) = wkc(ir1:nDx,:) - wks(ir1:nDx,:) ! H(1/a)-H(rsmi)
         endif
         !  --- Special treatment of on-site term --- Subtract Ewald contribution from G-vectors, hsm(p,a).
         !      Already generated by hansr4: lc=0
         !      Case p eq 0, wk(1,lc+1) = -hsm(p->0,a)
         !                   wk(1,lc)   = -2*hsmp(p->0,a)
         !      Case p ne 0, wk(1,lc+1) = hsm(p,rsm->0) - hsm(p,a)
         !                   wk(1,lc)   = 2*hsmp(p,rsm->0) - 2*hsmp(p,a)
         !   ... For p>0,job1=2, convert w(1,lc..lc+1) to -hsm and -2*hsmp
         if (ir1 == 2 .AND. .NOT. dcmpre(rsq(1))) then
            if (exi(ie) > 0) call rx(' hsmq: not implement e>0')
            !     ... Subtract 2*h0d from wk(1,lc), making it into -2*hsmp(p,a) And h0 from wk(1,lc+1), making it into -hsm(p,a)
            r = dsqrt(rsq(1))
            akap = dsqrt(-exi(ie))
            wkc(1,0:1) = wkc(1,0:1) - [exp(-akap*r)/akap, exp(-akap*r)/r] ![h0,h0d]
            xx = 4d0*a*y0*dexp(-(akap/a/2)**2 - (r*a)**2)
            do l = 1, lx(ie) ! Generate -hsm(p,a,l>0) by upward recursion
               wkc(1,l+1) = ((2*l-1)*wkc(1,l) -exi(ie)*wkc(1,l-1) + xx*(2*a**2)**(l-1)) /r**2
            enddo
         endif
         !   ... Make exp(iqR)*(H(rsm,r)-H(1/a,r))
         do  l = 0, lx(ie)
            i1=l*l+1
            i2=l*l+2*l+1
            do ilm=i1,i2
               hsm(ilm,ie)  = hsm(ilm,ie)+ sum([(yl(ir,ilm)*wkc(ir,l+1)*phase(ir),ir=1,nDx)])
               hsmp(ilm,ie) = hsmp(ilm,ie) + sum([(yl(ir,ilm)*wkc(ir,l)  *phase(ir),ir=1,nDx)])/2d0
            enddo   
         enddo
         le=(lx(ie)+1)**2
       endblock hansblock
30  enddo
  end subroutine hsmq
  subroutine hsmqe0(lmax,rsm,job1,q,p,nrx,nlmx,yl,awald,alat,qlv,nG,dlv,nD,vol,hsm)! Bloch-sum of smooth Hankel functions for energy 0
    use m_ropyln,only: ropyln
    use m_shortn3_plat,only: shortn3_plat,nout,nlatout,plat=>platx,qlat=>qlatx
    implicit none
    intent(in)::    lmax,rsm,job1,q,p,nrx,nlmx,   awald,alat,qlv,nG,dlv,nD,vol
    intent(out)::                              yl,                             hsm
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
    !o   hsm:    smoothed Hankels for e=0,q=0 for l=0..lmax
    ! wk:     (*,1) holds r**2
    !         (*,2) holds Y0 exp(-(r/rsm)**2) (rsm>0 only)
    !         (*,3) holds cos(q.dlv)
    !         (*,4) holds sin(q.dlv)
    !r Remarks
    !r   hsmqe0 is an adaptation of hsmq for for e=0, and possibly q=0.
    !r   For l=0, hsm(q=0,e->0) diverges.  To achieve a finite hsm(l=0,e=0), the average lim(e->0) hsm0(e) = -sqrt(4*pi)*exp(gamma*e) / (vol*e)
    !r   is subtracted from hsm(l=0).
    !b Bugs and Limitations:  see hsmq.
    ! ---------------------------------------------------------------
    integer :: lmax,nG,nD,nlmx,nrx,job, ir,ilm,l,m,ir1,lc,ls,job0,job1,job2,job3,lm,nDx, lx1(1),i1,i2,ix
    real(8) :: qdotr,a2,sp,gam,tpiba,p1(3),ex1(1), x1,x2,xx,xx1((lmax+1)**2),xx2((lmax+1)**2),r,h0,q0(3),a,faca
    real(8):: alat,awald,vol,rsm,p(3),q(3), wk(nrx,(2*lmax+10)),yl(nrx,*),qlv(3,*),dlv(3,*),y0fac(nrx),rsq(nrx)
    complex(8):: hsm((lmax+1)**2)  
    parameter (faca=1d0)
    complex(8):: xxc,img=(0d0,1d0),phasex,phase(nrx)
    real(8):: pp(3)
    real(8),parameter:: pi = 4d0*datan(1d0), y0 = 1/dsqrt(4d0*pi),  tpi = 2d0*pi
    if (nrx < max(nD,nG)) call rx('hsmqe0: nrx < nD or nG')
    ! ... Shorten pp=matmul(transpose(qlat),p); call shortn3_plat(pp);p1 = matmul(plat,pp+nlatout(:,1))
    p1=p 
    hsm=0d0 
    q0=q
    a = awald ! ... If rsm ge 1/a, set a = 1/rsm and skip r.s. part
    nDx = nD
    if (nG > 0.and. rsm > faca/a) then
       nDx = 1
       a = 1/rsm
    endif
    a2 = a*a
    qqblock: block ! --- Setup for Q-space part ---
      real(8):: qq(nG,3)
      complex(8):: xxcc(nG)
      if (nG > 0) then
         tpiba = 2d0*pi/alat
         gam = 0.25d0/a2
         forall(ir = 1:nG) qq(ir,:) = tpiba*(q0(:) + qlv(:,ir))
         call ropyln(nG,qq(1,1),qq(1,2),qq(1,3),lmax,nrx,yl,rsq)
         xxcc = -[(exp(-gam*rsq(ir)+ img*alat*sum(qq(ir,:)*p1) ),ir=1,nG)]
         lx1(1) = lmax
         ex1(1) = 0d0
         call pvhsmq(1,1,lmax,lx1,ex1,a,vol,nG,rsq,nrx,yl,nlmx,xxcc, hsm,hsm) !Q-space part hsm(1/a,e,r) of reduced strx for all energies ---
      endif
    endblock qqblock
    ppblock: block ! --- Energy-independent setup for direct-space part ---
      real(8):: pp(nDx,3)
      forall(ir = 1:nDx) pp(ir,:) = alat*(p1(:)-dlv(:,ir))
      call ropyln(nDx,pp(1,1),pp(1,2),pp(1,3),lmax,nrx,yl,rsq)
    endblock ppblock
    phase(1:nD)=[(exp(img*tpi*sum(q0*dlv(:,ir))), ir=1,nD)]
    ! ... If we are doing Ewald sums, we need this to make H(1/a) (hansr5)
    if (nG /= 0) y0fac(1:nDx) = [(y0*exp(-rsq(ir)*a2),ir=1,nDx)]
    ! --- D-space part of reduced strx for each energy: chi(l)=wk(l+lc)---
    ir1 = 1
    if (dcmpre(rsq(1))) then ! ... Case connecting vector p=0
       if (job1 == 1) ir1 = 2
       if (job1 == 0 .AND. dcmpre(rsm)) call rx('hsmqe0: job=0 and rsm=0 not allowed')
    endif
    lc = 1
    ls = lmax+8
    if (nG > 0) then
       call hansr5(rsq,lmax,nrx,nDx,1/a,y0fac,wk(1,lc)) ! ... chi = H(rsm) - H(1/a) = (H(0) - H(1/a)) - (H(0) - H(rsm))
    else
       wk(ir1:nDx,lc:lc+lmax+1) = 0d0
    endif
    if(.NOT. dcmpre(rsm)) then  !wk(1:nDx,2) = y0*dexp(-wk(1:nDx,1)/rsm**2)!   ... Make Y0 exp(-(r/rsm)**2) !  dcmpre(rsm,0d0)=T for hvccfp0
       y0fac(1:nDx)=[(y0*exp(-rsq(ir)/rsm**2),ir=1,nDx)]
       call hansr5(rsq,lmax,nrx,nDx,rsm,y0fac,wk(1,ls)) !wk(1,ls)='H(0) - H(rsm)'
       do    l  = 0, lmax
          do ir = ir1, nDx
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
    if (ir1 == 2 .AND. .NOT. dcmpre(rsq(1))) then !   ... Subtract h0(0..1) from wk(1,lc..lc+1), making it into -hsm(p,a)
       r = dsqrt(rsq(1))
       h0 = 1/r
       wk(1,lc) = wk(1,lc) - h0
       if (lmax > 0) wk(1,lc+1) = wk(1,lc+1) - h0/rsq(1)
       xx = 4d0*a*y0*dexp(-(r*a)**2)
       do l = 2, lmax
          wk(1,l+lc) = ((2*l-1)*wk(1,l+lc-1) + xx*(2*a**2)**(l-1))/rsq(1) !Generate -hsm(p,a,l>0) by upward recursion
       enddo
    endif
    ! ... Make sin(qR)*(H(rsm,r)-H(1/a,r)), cos(qR)*(H(rsm,r)-H(1/a,r))
    llooop: do l = 0, lmax
       i1=l*l+1
       i2=l*l+2*l+1
       do ix=i1,i2
          hsm(ix)  = hsm(ix) + sum([(yl(ir,ix)*phase(ir)*wk(ir,l+lc),ir=1,nDx)])
       enddo   
    enddo llooop
    if (dcmpre(sum(dabs(q0)))) hsm(1) = hsm(1) + rsm**2/(4d0*vol*y0) !! --- Add extra term for l=0 when e=0 and q=0 ---
  end subroutine hsmqe0
  subroutine hansr4(rsq,lmax,nrx,nr,e,rsm,wk,chi)! Difference between true,smoothed hankels for l=-1...lmax  for e nonzero.  Vectorizes for negative e.
    use m_ftox
    ! we need to simplify description here.....
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
    implicit none
    real(8),parameter::   pi = 4d0*datan(1d0),  y0 = 1/dsqrt(4d0*pi)
    integer :: nrx,nr,lmax
    real(8):: rsq(nr),e,chi(nrx,-1:*),rsm,wk(nr),wk2(nr),erfcee
    logical :: lpos
    integer :: l,ir,ir1
    real(8):: sre,r2,xx,ra,h0,arsm,earsm,kappa=0d0,ta=1d99, akap,a,r,um,up,x,facgl
    if (e>0d0) call rx('EH >0 not supported')
    if (lmax < 0 .OR. nr == 0) return
    ! --- chi(*,-1), chi(*,0) = true - sm hankel for each r ---
    ! ... chi(r,-1) = (h0 - um - up)*r/akap   and
    !     chi(r,0)  = (h0 - um + up)          with
    !     h0 = exp(-akap*r)/r
    !     up = erfc(akap*rsm/2+r/rsm)*exp(akap*r)/(2r)
    !     um = erfc(akap*rsm/2-r/rsm)/exp(akap*r)/(2r)
    !     um,up correspond to umins/(2r),uplus/(2r) in hansmr
    ! ... See Remarks for derivation of expressions in this case
    a = 1/rsm
    akap = dsqrt(-e)
    arsm = akap*rsm/2
    earsm = dexp(-arsm**2)/2
    facgl = 8*a*earsm
    ir1 = 1
    if(rsq(1) < 1d-12) then
       ir1 = 2
       chi(1,-1:0) = [-erfc(arsm)/akap,akap*erfc(arsm) - 4d0*a*earsm/dsqrt(4d0*datan(1d0))] !  chi(-1) -> erfc(akap/2a)/akap for r=0
    endif
    do ir = ir1, nr !Make chi(r,rsm->0,l) - chi(r,rsm,l) for l=-1, 0
       r2 = rsq(ir)
       r = dsqrt(r2)
       ra = r*a
       sre = akap*r
       h0 = dexp(-sre)/r
       xx = earsm*wk(ir)/r
       x = ra - arsm
       um = merge(h0-xx*erfcee(x),xx*erfcee(x),x>0d0)
       up = xx*erfcee(ra + arsm) ! = xx* erfc(x )/y0/exp(-x**2) !     ... Evaluation of up assumes x gt 0
       chi(ir,-1:0) = [(h0 - um - up)*r/akap, h0 - um + up] ! Make chi(r,rsm->0,l) - chi(r,rsm,l) for l=-1, 0
       wk2(ir) = facgl*wk(ir)
    enddo
    ! --- chi(ir,l) for l>1 by upward recursion ---
    facgl = 2d0*a**2
    chi(1,1:lmax) = 0
    do l = 1, lmax
       chi(ir1:nr,l) = [(((2*l-1)*chi(ir,l-1) - e*chi(ir,l-2) + wk2(ir))/rsq(ir),ir=ir1,nr)]
       wk2(ir1:nr) = facgl*wk2(ir1:nr)
    enddo
  end subroutine hansr4
  subroutine hansr5(rsq,lmax,nrx,nr,rsm,wk,chi) !Difference between true,smoothed hankels for l=0...lmax, e=0
    use m_ftox
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
    !r   The first point may have r->0.  In this case, chi(1,l)=0 except chi(1,0), which is returned as the -hsm(e=0,r=0) i.e.
    !r   the infinite unsmoothed part is discarded.
    implicit none
    integer :: nrx,nr,lmax,l,ir,ir1
    real(8) :: rsq(nr),chi(nrx,0:lmax),rsm,wk(nr),wk2(nr), xx,ra,a,r,facgl,erfcee
    real(8),parameter::   pi = 4d0*datan(1d0),  y0 = 1/dsqrt(4d0*pi)
    if (lmax < 0 .OR. nr == 0) return
    a = 1/rsm
    ir1 = 1
    if (rsq(1) < 1d-6) then
       ir1 = 2
       chi(1,0) = -2d0*a/dsqrt(4d0*datan(1d0))!   ... hsm(e=0,r=0) = 4*a*y0
    endif
    do  ir = ir1, nr
       r = rsq(ir)**.5d0
       chi(ir,0)= wk(ir)/r * erfcee(r*a) !=xx*erfc(ra)/y0/exp(-ra**2)
       wk2(ir) = 4d0*a*wk(ir)
    enddo
    ! --- chi(ir,l) for l>0 by upward recursion ---
    chi(1,1:lmax)=0d0
    do   l = 1, lmax
       chi(ir1:nr,l) = [(((2*l-1)*chi(ir,l-1) + wk2(ir))/rsq(ir),ir=ir1,nr)]
       wk2(ir1:nr)   = 2*a**2*wk2(ir1:nr)
    enddo
  end subroutine hansr5
  subroutine pvhsmq(job,nx,lmxa,lxi,exi,a,vol,n,qsq,nyl,yl,nlmx,cssn, hl,hlp) ! Makes H_L(1/a,e,r) by summation over reciprocal lattice vectors.
    !o  hl,hlp
    !  job=1 : make hl only (not hlp)
    !  Skips any ie for which lxi(ie)<lmxa. This is a kernel called by hsmq.
    implicit none
    integer :: n,nx,lmxa,lxi(nx),nlmx,job,nyl,nlm,i,i1,ie,ilm,l,lx,m,je,lm,li,le
    double precision ::e,gam,vol, a,exi(1),qsq(n),yl(nyl,*), c,s,xx,x1,x2
    complex(8):: hl(nlmx,*),hlp(nlmx,*)
    complex(8):: cssn(n),cof0,cof
    real(8),parameter::pi = 4d0*datan(1d0)
    gam = 0.25d0/a**2
    do  20  ie = 1, nx !For each exi ne 0 and ilm do ---
       e = exi(ie)
       lx = lxi(ie)
       if (lxi(ie) < lmxa) cycle
       do  je = 1, ie-1 ! Copy existing hl,hlp if already calculated for this exi
          if (dabs(exi(je)-exi(ie))> 1d-8 .or. lxi(ie)>lxi(je)) cycle
          hl (1:(lxi(ie)+1)**2,ie)= hl (1:(lxi(ie)+1)**2,je)
          hlp(1:(lxi(ie)+1)**2,ie)= hlp(1:(lxi(ie)+1)**2,je)
       enddo
       nlm=(lx+1)**2
       i1 = 1
       if (dcmpre(e-qsq(1))) i1 = 2
       nlmb: block
         complex(8):: img=(0d0,1d0), c2(nlm),c3(nlm)
         c2 = [(sum(yl(i1:n,ilm)*cssn(i1:n)/(e-qsq(i1:n))   ),ilm=1,nlm)]
         if(job/=1) c3 = [(sum(yl(i1:n,ilm)*cssn(i1:n)/(e-qsq(i1:n))**2),ilm=1,nlm)]
         cof0 = 4d0*pi*dexp(gam*e)/vol
         ilm = 0
         do l = 0, lx
            cof = (-img)**l*cof0
            li= l**2+1
            le= l**2+2*l+1
            hl(li:le,ie)  = hl(li:le,ie)  + cof*c2(li:le)
            if(job/=1) hlp(li:le,ie) = hlp(li:le,ie) + cof*(gam*c2(li:le)-c3(li:le))
         enddo
       endblock nlmb
       if(dcmpre(e-qsq(1))) hl(1,ie)= hl(1,ie) - sqrt(4d0*pi)/vol*gam !Add extra term for l=0 when e=0 and q=0 
20  enddo
  end subroutine pvhsmq
  pure logical function dcmpre(x1)
    real(8),intent(in):: x1
    dcmpre = dabs(x1) < 1d-8
  end function dcmpre
end module m_hsmq

