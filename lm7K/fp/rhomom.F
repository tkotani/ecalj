Cdelw1 var ['orofi', '0', 'orofi', 'nrmx', 'real(8)', 'rv_w_', 'rhomom']
Cdelw1 var ['orwgt', '0', 'orwgt', 'nrmx', 'real(8)', 'rv_w_', 'rhomom']
Cdelw1 var ['oh', '0', 'oh', 'nrmx', 'real(8)', 'rv_w_', 'rhomom']
Cdelw1 var ['ov', '0', 'ov', 'nrmx', 'real(8)', 'rv_w_', 'rhomom']
Cdelw1 undel [['call'], 'orhoat', 'rhomom']
Cgetarg...info...           structure ['rhomom', 'site', 'ssite']
Cgetarg...info...           structure ['rhomom', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['rhomom', 'm_struc_def']
      subroutine rhomom(nbas,ssite,sspec,orhoat,qmom,vsum)
       
       use m_struc_def  !Cgetarg

C- Multipole moments of valence sphere densities
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct containing site-specific information
Ci     Elts read: spec
Ci   sspec :struct containing species-specific information
Ci     Elts read: lmxl z qc a nr rmt rg
Ci     Passed to: corprm
Ci   orhoat:vector of offsets containing site density
Co Outputs
Co   qmom  :Multipole moments, stored as a single long vector.
Co         := integral r^l Y_L (rho1-rho2) + core spillout
Co   vsum  :sum over all sites ib of difference vs1_ib - vs2_ib
Co         :vs1_ib = integral in sphere of estat potential[true rho]
Co         :vs2_ib = integral in sphere of estat potential[sm rho]
Cr Remarks
Cu Updates
Cu   01 Jul 05 handle sites with lmxl=-1 -> no augmentation
Cu   11 Jun 00 spin polarized
Cu   22 Apr 00 Adapted from nfp rhomom.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
Cdelw1 do not change  [['call'], 'orhoat', 'rhomom']
      integer nbas,orhoat(3,1)
Cgetarg       double precision ssite(1),sspec(1),qmom(1),vsum
       real(8):: qmom(1) , vsum 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)

C ... Local parameters
Cdelw1       integer stdo,ipr,iprint,lgunit,orofi,nrmx,orwgt,oh,ov,nsp,nglob,
Cdelw1      .  j1,ib,is,igetss,lmxl,nr,nlml,ilm,j,lfoc
       integer:: stdo , ipr , iprint , lgunit , nrmx , nsp , nglob , 
     . j1 , ib , is , igetss , lmxl , nr , nlml , ilm , j , lfoc 
       real(8) ,allocatable :: rv_w_orofi(:)
       real(8) ,allocatable :: rv_w_orwgt(:)
       real(8) ,allocatable :: rv_w_oh(:)
       real(8) ,allocatable :: rv_w_ov(:)

      parameter( nrmx=1501)
      double precision fpi,y0,z,qc,a,rmt,qcorg,qcorh,qsc,cofg,cofh,rg,
     .  ceh,rfoc,vs1,vs2
C ... Heap
      integer w(1)
      common /w/ w

C --- Setup ---
      ipr  = iprint()
      stdo = lgunit(1)
      fpi  = 16d0*datan(1d0)
      y0   = 1d0/dsqrt(fpi)
      nsp  = nglob('nsp')
Cdelw1       call defrr(orofi,nrmx)
       allocate(rv_w_orofi(nrmx))
Ctakao_ZeroClear_NotRequiered        if (nrmx<0) rv_w_orofi(:)=0.0d0

Cdelw1       call defrr(orwgt,nrmx)
       allocate(rv_w_orwgt(nrmx))
Ctakao_ZeroClear_NotRequiered        if (nrmx<0) rv_w_orwgt(:)=0.0d0

Cdelw1       call defrr(oh,   nrmx)
       allocate(rv_w_oh(nrmx))
Ctakao_ZeroClear_NotRequiered        if (nrmx<0) rv_w_oh(:)=0.0d0

Cdelw1       call defrr(ov,   nrmx)
       allocate(rv_w_ov(nrmx))
Ctakao_ZeroClear_NotRequiered        if (nrmx<0) rv_w_ov(:)=0.0d0

      if (ipr .ge. 40) write(stdo,221)

C --- Loop over sites ---
      j1 = 1
      vsum = 0d0
      do  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
Cgetarg...info...           real(8) :: qc , real(8) :: sspec%qc 1
Cgetarg         call upack('spec lmxl z qc',sspec,is,lmxl,z,qc,0)
         
         lmxl=sspec(is)%lmxl
         z=sspec(is)%z
         qc=sspec(is)%qc

Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg         call upack('spec a nr rmt',sspec,is,a,nr,rmt,0)
         
         a=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt

Cgetarg...info...           real(8) :: rg , real(8) :: sspec%rg 1
Cgetarg         call upack2('spec rg',sspec,is,rg)
         
         rg=sspec(is)%rg

        if (lmxl .eq. -1) goto 10
        call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
        qc = qcorg+qcorh
        nlml = (lmxl+1)**2
Cdelw1         call radmsh(rmt,a,nr,w(orofi))
         call radmsh ( rmt , a , nr , rv_w_orofi ) 

Cdelw1         call radwgt(rmt,a,nr,w(orwgt))
         call radwgt ( rmt , a , nr , rv_w_orwgt ) 

Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1 3 is digit, dropped
Cdelw1         call pvrhom(nlml,nr,nsp,w(orofi),w(orwgt),w(orhoat(1,ib)),
Cdelw1      .    w(orhoat(2,ib)),w(orhoat(3,ib)),cofg,cofh,rg,ceh,rfoc,z,w(oh),
Cdelw1      .    qmom(j1))
         call pvrhom ( nlml , nr , nsp , rv_w_orofi , rv_w_orwgt , w ( 
     .   orhoat ( 1 , ib ) ) , w ( orhoat ( 2 , ib ) ) , w ( orhoat ( 
     .   3 , ib ) ) , cofg , cofh , rg , ceh , rfoc , z , rv_w_oh , qmom 
     .   ( j1 ) ) 

        if (ipr .ge. 40) then
          write(stdo,220) ib,1,qmom(j1),qmom(j1)/y0,qc,z
          do  ilm = 2, nlml
            j = j1+ilm-1
            if (dabs(qmom(j)).gt.1d-6) write(stdo,220) ib,ilm,qmom(j)
          enddo
        endif
  220   format(i13,i6,f12.6,f12.6,2f9.2)
  221   format(/' rhomom:   ib   ilm      qmom',8x,'Qval',7x,
     .     'Qc',8x,'Z')
Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1 3 is digit, dropped
Cdelw1         call pvrhm2(z,a,nr,nlml,nsp,w(orofi),w(orwgt),w(orhoat(1,ib)),
Cdelw1      .    w(orhoat(2,ib)),w(orhoat(3,ib)),qmom(j1),cofg,cofh,
Cdelw1      .    rg,ceh,rfoc,w(oh),w(ov),vs1,vs2)
         call pvrhm2 ( z , a , nr , nlml , nsp , rv_w_orofi , rv_w_orwgt 
     .   , w ( orhoat ( 1 , ib ) ) , w ( orhoat ( 2 , ib ) ) , w ( orhoat 
     .   ( 3 , ib ) ) , qmom ( j1 ) , cofg , cofh , rg , ceh , rfoc , 
     .   rv_w_oh , rv_w_ov , vs1 , vs2 ) 


        vsum = vsum+vs1-vs2
        j1 = j1+nlml

   10   continue
      enddo

Cdelw1 rlse name= orofi old_list= orofi orwgt oh ov 
Cdelw1 rlse name= orofi new_list= (None)
Cdelw1       call rlse(orofi)
       if (allocated(rv_w_ov)) deallocate(rv_w_ov)
       if (allocated(rv_w_oh)) deallocate(rv_w_oh)
       if (allocated(rv_w_orwgt)) deallocate(rv_w_orwgt)
       if (allocated(rv_w_orofi)) deallocate(rv_w_orofi)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [orhoat]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine rhomom 


      subroutine pvrhom(nlml,nr,nsp,rofi,rwgt,rho1,rho2,rhoc,
     .  cofg,cofh,rg,ceh,rfoc,z,h,qmom)

C- Multipole moments for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlml  :L-cutoff for charge
Ci   nr    :number of radial mesh points
Ci   nsp   :number of spins
Ci   rofi  :radial mesh points
Ci   rwgt  :radial integration weights
Ci   rho1  :local true density*r**2, tabulated on a radial mesh
Ci   rho2  :local smoothed density*r**2, tabulated on a radial mesh
Ci   rhoc  :core density
Ci   cofg  :coefficient to Gaussian part of pseudocore density (corprm)
Ci   cofh  :coefficient to Hankel part of pseudocore density (corprm)
Ci   rg    :smoothing radius for compensating gaussians
Ci   ceh   :energy of hankel function to fit core tail
Ci   rfoc  :smoothing radius for hankel head fitted to core tail
Ci   z     :nuclear charge
Ci   h     :work array of dimension nr
Co Outputs
Co   qmom  :multipole moments for one site
Cr Remarks
Cr   Q_L = integral r^l (rho1-rho2) + l=0 contr. from core spillout
Cr   The core spillout term is:
Cr      qcore(rhoc)-z  - sm_qcore-sm_qnuc
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlml,nr,nsp
      double precision ceh,cofg,cofh,rfoc,rg,z
      double precision rofi(1),rwgt(1),h(1),qmom(nlml),rhoc(nr,nsp),
     .  rho1(nr,nlml,nsp),rho2(nr,nlml,nsp)
C ... Local parameters
      integer n0,i,ilm,l,m,lmxl,ll,isp
      parameter (n0=10)
      double precision ag,delq,fac,gnu,pi,qcor1,qcor2,qnuc2,r,
     .  rhochs,rhocsm,rhonsm,srfpi,sum,sumg,xi(0:n0),qcor1s
C     double precision sumh,samh,y0

      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      lmxl = ll(nlml)
      do  i = 1, nr
        h(i) = rwgt(i)
      enddo
      ilm = 0
      do  l = 0, lmxl
        do m = -l, l
          ilm = ilm+1
          sum = 0d0
          do  isp = 1, nsp
          do  i = 1, nr
            sum = sum + h(i)*(rho1(i,ilm,isp)-rho2(i,ilm,isp))
          enddo
          enddo
          qmom(ilm) = sum
        enddo
        do  i = 1, nr
          h(i) = h(i)*rofi(i)
        enddo
      enddo

C ... l=0 includes core; some might be spilled out
      ag = 1d0/rg
      fac = 4*pi*(ag*ag/pi)**1.5d0

C ... Renormalize gaussian
      sumg = 0d0
      do  i = 2, nr
        r = rofi(i)
        gnu = fac* r*r * dexp(-ag*ag*r*r)
        sumg = sumg + rwgt(i)*gnu
      enddo
      fac = fac/sumg

C     Make: 
C     qcor1 = true core charge inside rmax
C     qcor2 = smooth core charge inside rmax * fac
C     qnuc2 = smooth nuclear charge inside rmax * fac
C     delq = (true core q-z) - (sm core q - sm z) * fac
C     sumh = 0d0
      qcor2 = 0d0
      qnuc2 = 0d0
      qcor1 = 0d0
      qcor1s= 0d0
      do  i = 2, nr
        r = rofi(i)
        gnu = fac * r*r * dexp(-ag*ag*r*r)
        call hansmr(r,ceh,1d0/rfoc,xi,1)
        rhonsm = -z*gnu
        rhochs = srfpi*cofh*xi(0)*r*r
        rhocsm = srfpi*cofg*gnu + rhochs
C       sumh = sumh + rwgt(i)*rhochs
        qcor2 = qcor2 + rwgt(i)*rhocsm
        qnuc2 = qnuc2 + rwgt(i)*rhonsm
        qcor1 = qcor1 + rwgt(i)*rhoc(i,1)
        qcor1s= qcor1s + rwgt(i)*rhoc(i,nsp)
      enddo
      if (nsp .eq. 2) qcor1 = qcor1+qcor1s
C     samh = -y0*cofh*4d0*pi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh
      delq = qcor1-z - qcor2-qnuc2
      qmom(1) = qmom(1) + delq/srfpi

C      y0 = 1d0/srfpi
C      samh = -y0*cofh*4d0*pi*dexp(ceh*rfoc*rfoc*0.25d0)/ceh
C      write(stdo,942) samh,sumh,samh-sumh
C  942 format(' integral in smH core:',2f12.6/
C     .   ' Core spill-out charge is',f12.6)
C      write(stdo,821) delq
C  821 format('delq=',f12.6)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pvrhom 


      subroutine pvrhm2(z,a,nr,nlml,nsp,rofi,rwgt,rho1,rho2,rhoc,qmom,
     .  cofg,cofh,rg,ceh,rfoc,h,v,vsum1,vsum2)

C- Integral over electrostatic potential, to fix energy origin
C ----------------------------------------------------------------------
Ci Inputs
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge
Ci   nsp   :number of spins
Ci   rofi  :radial mesh points
Ci   rwgt  :radial integration weights
Ci   rho1  :local true density, tabulated on a radial mesh
Ci   rho2  :local smoothed density, tabulated on a radial mesh
Ci   rhoc  :core density
Ci   qmom  :multipole moments for one site (only l=0 term used)
Ci   cofg  :coefficient to Gaussian part of pseudocore density (corprm)
Ci         := y0 * pseudocore charge
Ci   cofh  :coefficient to Hankel part of pseudocore density (corprm)
Ci   rg    :smoothing radius for compensating gaussians
Ci   ceh   :energy of hankel function to fit core tail
Ci   rfoc  :smoothing radius for hankel head fitted to core tail
Ci   h     :work array of dimension nr
Ci   v     :work array (holds spherical estat potential - poiss0.f)
Co Outputs
Co   vsum1 :integral in sphere of electrostatic potential for true rho
Co   vsum2 :integral in sphere of electrostatic potential for smooth rho
Co         :including compensating gaussians
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,nlml,nsp
      double precision a,ceh,cofg,cofh,rfoc,rg,vsum1,vsum2,z,rofi(1),
     .  rwgt(1),rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),h(1),v(nr),qmom(1),
     .  rhoc(nr,nsp)
C ... Local parameters
      integer n0,i
      parameter (n0=10)
      double precision af,ag,b,cg,fac,facs,fpi,gnu,pi,q1,q2,r,srfpi,
     .  vhrho,vsum,y0,xi(0:n0)

      pi = 4d0*datan(1d0)
      fpi = 4d0*pi
      srfpi = dsqrt(fpi)
      y0 = 1d0/srfpi
      b = rofi(nr)/(dexp(a*nr-a)-1d0)

C ... True density
      q1 = 0d0
      facs = 1d0/(3-nsp)
      do  i = 1, nr
        h(i) = facs*(srfpi*(rho1(i,1,1)+rho1(i,1,nsp))
     .                    + rhoc(i,1) + rhoc(i,nsp))
        q1 = q1 + rwgt(i)*h(i)
      enddo
      call poiss0(z,a,b,rofi,h,nr,0d0,v,vhrho,vsum,1)
      vsum1 = 0d0
      do  i = 2, nr
        r = rofi(i)
        vsum1 = vsum1 + 4*pi*rwgt(i)*r*r*(v(i)-2*z/r)
      enddo

C ... Smooth density, including compensating gaussians
      cg = qmom(1) + cofg - y0*z
      ag = 1d0/rg
      af = 1d0/rfoc
      fac = fpi*(ag*ag/pi)**1.5d0
      q2 = 0d0
      facs = 1d0/(3-nsp)
      do  i = 1, nr
        r = rofi(i)
        gnu = fac* r*r * dexp(-ag*ag*r*r)
        call hansmr(r,ceh,af,xi,1)
        h(i) = srfpi*(facs*(rho2(i,1,1)+rho2(i,1,nsp))
     .                + cg*gnu + cofh*xi(0)*r*r)
        q2 = q2 + rwgt(i)*h(i)
      enddo
      call poiss0(0d0,a,b,rofi,h,nr,0d0,v,vhrho,vsum,1)
      vsum2 = 0d0
      do  i = 2, nr
        r = rofi(i)
        vsum2 = vsum2 + 4*pi*rwgt(i)*r*r*v(i)
      enddo

C      write (stdo,400) q1,q2,vsum1,vsum2
C  400 format(' qsph=',2f12.6,'   vsum=',2f12.6)
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pvrhm2 


