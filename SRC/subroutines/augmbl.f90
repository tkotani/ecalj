module m_augmbl !Add augmentation part to H and S. aughsoc add SO part to H.
                !Inputs are Site integrals, sig,tau,pi,hso 
  public augmbl,aughsoc
  private
contains
  subroutine augmbl(isp,q,osig,otau,oppi,ndimh, h,s)  !Add augmentation part to H and S. 
    use m_lmfinit,only: nsp,nlmto, sspec=>v_sspec
    use m_lmfinit,only: nbas,nkaph,alat=>lat_alat,ispec
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol,rv_a_opos
    use m_bstrux,only: Bstrux_set, bstr
    use m_orbl,only: Orblib, norb,ltab,ktab,offl
    use m_struc_def,only: s_cv1,s_rv1,s_rv4
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   isp   :current spin channel
    !i   q     :Bloch wave number
    !i   osig  :overlap matrix of P_kL
    !i         :NB: also head-head, head-tail contributions; see augmat.f
    !i   otau  :kinetic energy matrix of P_kL
    !i         :NB: also head-head, head-tail contributions; see augmat.f
    !i         (otau is not needed because folded into ppi already)
    !i   oppi  :kinetic energy + potential matrix of P_kL
    !i         :NB: also head-head, head-tail contributions; see augmat.f
    !i   ndimh :dimension of h and s
    !i   napw  :number of PWs in APW part of basis
    !i   igapw :PWs in units of reciprocal lattice vectors
    !o Outputs
    !o   h     :augmentation part of hamiltonian matrix added to h
    !o   s     :augmentation part of overlap matrix added to s
    !l Local variables
    !l   nkaph :number of orbital types for a given L quantum no. in basis
    !l         :at augmentation site ia, including local orbitals
    !l   nlmto :number of lmto basis functions
    !r Remarks
    !r   Some expressions labelled JMP refer to J.Math.Phys39, 3393 (1998)
    ! ----------------------------------------------------------------------
    implicit none
    type(s_cv1),target :: oppi(3,nbas)
    type(s_rv1),target :: otau(3,nbas)
    type(s_rv4),target :: osig(3,nbas)
    integer:: isp , ndimh , napw 
    real(8):: q(3)
    complex(8):: h(ndimh,ndimh),s(ndimh,ndimh)
    integer,parameter :: ktop0=20, nlmbx=49, nlmax=49
    complex(8),allocatable:: b(:,:,:),bx(:,:,:),bb(:,:,:)
    integer:: ibas , isa , kmax , lmxa , lmxb ,  nglob , nlma, nlmb
    real(8):: rsma,pa(3),xx 
    integer:: initbas, endbas,lm,iq,nh,np,i1,i2,ilm1,ilm2,k1,k2
    logical:: debug=.false.
    complex(8),pointer:: ppi1(:),ppi2(:),ppi3(:),Lm1(:),Lm2(:),Lm3(:),Lz1(:),Lz2(:),Lz3(:)
    real(8),pointer:: sig1(:),sig2(:),sig3(:)
    !--------------------------
    call tcn ('augmbl')
    do ibas = 1,nbas
       isa = ispec(ibas) 
       pa  = rv_a_opos(:,ibas) 
       lmxa=sspec(isa)%lmxa !max l of augmentation
       lmxb=sspec(isa)%lmxb !max l of basis
       kmax=sspec(isa)%kmxt !max of radial k
       nlmb = (lmxb+1)**2
       nlma = (lmxa+1)**2
       if (lmxa == -1) cycle
       call bstrux_set(ibas,q) !--- Set bloch sum of structure constant bstr
       ! to expand all orbitals at site ia. (Bloch sum of C^i_akL (C.1) in Ref.[1])
       allocate( b(0:kmax,nlma,ndimh) )
       do lm=1,nlma
          b(:,lm,:) = transpose(bstr(:,lm,:))
       enddo
       addaug: block
         complex(8) :: sig(0:kmax,0:kmax,0:lmxa,nsp),ppi(0:kmax,0:kmax,nlma,nlma,nsp),g(0:kmax,nlma)
         complex(8):: ppihh(nkaph,nkaph,nlmb,nlmb,nsp), ppihp(nkaph,0:kmax,nlmb,nlma,nsp)
         real(8)::    sighh(nkaph,nkaph,0:lmxb,nsp),    sighp(nkaph,0:kmax,0:lmxb,nsp)
         integer:: i1,i2,ilm,ll,l1,ik1,l2,ik2,ilm1,ilm2,j,iorb,jorb
         sighh=reshape(osig(3,ibas)%v,shape(sighh)) !augmentation head-head overlap matrix 
         sighp=reshape(osig(2,ibas)%v,shape(sighp)) !augmentation head-Pkl overlap matrix  
         ppihh=reshape(oppi(3,ibas)%cv,shape(ppihh)) !augmentation head-head potential matrix
         ppihp=reshape(oppi(2,ibas)%cv,shape(ppihp)) !augmentation head-Pkl potential matrix
         call orblib(ibas) !See use section. Return norb,ltab,ktab,offl
         do  iorb = 1, norb
            l1  = ltab(iorb)
            ik1 = ktab(iorb)
            do  ilm1 = l1**2+1, (l1+1)**2
               i1 = offl(iorb)+ilm1-l1**2  !Two-center terms
               s(i1,:) = s(i1,:) + [(       sum(sighp(ik1,:,l1,isp)*b(:,ilm1,j)),  j=1,ndimh)]
               s(:,i1) = s(:,i1) + [(dconjg(sum(b(:,ilm1,j)*sighp(ik1,:,l1,isp))), j=1,ndimh)]
               h(i1,:) = h(i1,:) + [(       sum(ppihp(ik1,:,ilm1,:,isp)*b(:,:,j)), j=1,ndimh)]
               h(:,i1) = h(:,i1) + [(dconjg(sum(b(:,:,j)*ppihp(ik1,:,ilm1,:,isp))),j=1,ndimh)]
            enddo
            do  ilm1 = l1**2+1, (l1+1)**2
               i1 = offl(iorb)+ilm1-l1**2  !Two-center terms
               do  jorb = 1, norb !one center terms
                  l2  = ltab(jorb)
                  ik2 = ktab(jorb)
                  do  ilm2 = l2**2+1, (l2+1)**2
                     i2 = offl(jorb)+ilm2-l2**2
                     h(i1,i2) = h(i1,i2) + ppihh(ik1,ik2,ilm1,ilm2,isp)
                     if (ilm1 == ilm2) s(i1,i2) = s(i1,i2) + sighh(ik1,ik2,l1,isp)
                  enddo
               enddo
            enddo
         enddo
         ppi=reshape(oppi(1,ibas)%cv,shape(ppi))
         do  i2 = 1, ndimh
            do  ilm = 1, nlma
               g(0:kmax,ilm) = [(sum(ppi(k1,:,ilm,:,isp)*b(:,:,i2)),k1=0,kmax)]
            enddo
            h(:,i2) = h(:,i2) + [(sum(dconjg(b(:,:,i1))*g(:,:)),i1=1,ndimh)]
         enddo
         sig= reshape(osig(1,ibas)%v,shape(sig))
         do i2 = 1, ndimh
            do ilm = 1, nlma
               g(:,ilm) = matmul(sig(:,:,ll(ilm),isp),b(:,ilm,i2))
            enddo
            s(1:i2,i2) = s(1:i2,i2) + [(sum( dconjg(b(:,:,i1))*g(:,:) ),i1=1,i2)]
         enddo
       endblock addaug
       deallocate(b)
    enddo
    call tcx ('augmbl')
  end subroutine augmbl
  subroutine aughsoc(qp,ohsozz,ohsopm,ndimh, hso)
    use m_orbl,only: Orblib, norb,ltab,ktab,offl
    use m_struc_def,only: s_cv1,s_rv1,s_sblock
    use m_lmfinit,only: nsp, lsox=>lso, nbas, nkaph, ispec, sspec=>v_sspec,socaxis
    use m_bstrux,only: Bstrux_set, bstr
    use m_lattic,only: plat=>lat_plat,qlat=>lat_qlat
    !!- Spin-orbit couping matrix hso
    !i   qp    :Bloch wave number
    !i   hsozz,hsopm : atomic parts of SOC (Lz and L-)
    !i   ndimh :dimension of halimtonian.
    !o   hso   :spin diagonal and off-diagonal block of spin-orbit hamiltonian
    ! note  'shorbz need to be improved in future (the method in shortn3)'.
    ! note  We obtain Lz,L+,and L- (Lzz Lmm Lpp) in this routine. From their linear combinatios, we have hso.
    implicit none
    type(s_sblock),target :: ohsozz(3,nbas),ohsopm(3,nbas)
    integer:: isp, ndimh, ibas, isa , kmax , lmxa , lmxb ,  nglob , nlma , nlmb,lso
    integer:: initbas, endbas,lm,iq,nh,np,isp1,isp2,nspx
    real(8):: q(3),qp(3),fac
    complex(8):: hso(ndimh,ndimh,3)
    complex(8),allocatable:: b(:,:,:)
    complex(8),pointer:: ppi1(:),ppi2(:),ppi3(:)
    type(s_sblock),pointer:: Lzz(:),Lmp(:)
    complex(8),allocatable:: SSbPP(:),SSbHP(:),SSbPH(:),SSbHH(:) !tailtail, headtail, headhead
    complex(8):: img=(0d0,1d0), facso(3,3), f1,f2,f3
    real(8)::d2
    logical,save:: init=.true.
    logical:: cmdopt0
    call tcn ('aughsoc')
    lso=lsox
    if(cmdopt0('--socmatrix')) lso=1
    if(lso==1) then
       if( sum(abs(socaxis-[0d0,0d0,1d0]))  < 1d-6) then
          !     Mixing matrix for Spin-block facso based on (Lz,L-,L+)
          !     (001)                              Lz   L-   L+
          facso(:,1) = [complex(8)::  1d0, 0d0, 0d0] ! diagonal part isp=1,isp=1
          facso(:,2) = [complex(8):: -1d0, 0d0, 0d0] ! diagonal part isp=2,isp=2
          facso(:,3) = [complex(8)::  0d0, 1d0, 0d0] ! off-diag part isp=1,isp=2
          facso=0.5d0* facso  ! prefactor 1/2
       elseif( sum(abs(socaxis-[1d0,1d0,0d0])) <1d-6) then
          !     (110)                               Lz         L-          L+
          !     Lx= 1/2 (L- + L+), Ly=1/2i (-L- + L+)
          !     Lx+Ly = (1/2-1/2i)L-  + (1/2+1/2i)L+
          !     Lx-Ly = (1/2+1/2i)L-  + (1/2-1/2i)L+
          d2= dsqrt(2d0)/4d0
          facso(:,1) = [complex(8)::  0d0,  d2-d2/img,  d2+d2/img ]
          facso(:,2) = [complex(8)::  0d0, -d2+d2/img, -d2-d2/img ]
          facso(:,3) = [complex(8):: -1d0,  d2*img+d2,  d2*img-d2 ]
          facso=0.5d0* facso  ! prefactor 1/2
       else
          call rx('Given HAM_SOCAXIS is not yet implemented')
       endif
    endif

    ! so=0  sumev=      -27.066983
    ! sumev=      -27.066983  val*vef=    -388.462096   sumtv=     361.395114

    ! so=2 ==> sumev=      -27.074354  val*vef=    -388.462096   sumtv=     361.387743
    ! sumev=      -27.074354  val*vef=    -388.462096   sumtv=     361.387743

    ! so=1 SO (001)          Lz   L-   L+
    !      facso(:,1) = [ 1d0, 0d0, 0d0]  ! diagonal part isp=1,isp=1
    !      facso(:,2) = [-1d0, 0d0, 0d0]  ! diagonal part isp=2,isp=2
    !      facso(:,3) = [ 0d0, 1d0, 0d0]  ! off-diag part isp=1,isp=2
    !      facso=0.5d0* facso        ! prefactor 1/2
    ! sumev=      -27.089234  val*vef=    -388.462096   sumtv=     361.372862

    ! so=1 BZ_SOCAXIS=1,1,0
    ! (110)             Lz      L-          L+  ! Taken from (A8) in Liqin2019,PhysRevB.99.054418
    !      d2= dsqrt(2d0)/4d0
    !      facso(:,1) = [complex(8)::  0d0,  d2-d2/img, d2+d2/img]
    !      facso(:,2) = [complex(8)::  0d0, -d2+d2/img,-d2-d2/img]
    !      facso(:,3) = [complex(8):: -1d0,  d2*img+d2, d2*img-d2]
    !      facso=0.5d0* facso        ! prefactor 1/2
    ! sumev=      -27.088951  val*vef=    -388.462096   sumtv=     361.373145

    hso=0d0
    q=qp 
    do ibas = 1,nbas
       isa =ispec(ibas) 
       lmxa=sspec(isa)%lmxa !max l of augmentation
       lmxb=sspec(isa)%lmxb !max l of basis
       kmax=sspec(isa)%kmxt !max of radial k
       nlmb = (lmxb+1)**2
       nlma = (lmxa+1)**2
       if (lmxa == -1) cycle
       call bstrux_set(ibas,q) !Make strux b to expand all orbitals at site ia
       if(allocated(b)) deallocate(b)
       allocate( b(0:kmax,nlma,ndimh) )
       do lm=1,nlma
          b(:,lm,:) = transpose(bstr(:,lm,:))
       enddo
       nh= nkaph*nlmb     ! size of head nh
       np= (kmax+1)*nlma  ! size of tail np
       !! Get Lzz,Lmp,Lmp(spinfliped)= (Lz,L-,L+)  See mkpot-locpot-augmat-gaugm-pvagm1,pvaglc to generate hsozz,hsopm
       Lzz => ohsozz(:,ibas)   ! Lz block  1:P*P, 2:H*P, 3:H*H for up and dn
       if(lso==1) Lmp => ohsopm(:,ibas) ! <up|L-|dn> for isp=1 , <dn|L+|up> for isp=2. See gaugm.F, pvagm1,pvaglc
       allocate(SSbPP(np*np),SSbHP(np*nh),SSbPH(np*nh),SSbHH(nh*nh))
       nspx=nsp
       if(lso==1) nspx=3
       do isp=1,nspx
          ! hso(:,:,isp=1) is (1,1) block in (A8) in in Liqin2019,PhysRevB.99.054418
          ! hso(:,:,isp=2) is (1,1) block
          ! hso(:,:,isp=3) is (1,2) block  =L- in the case of 001 spin axis.
          if(isp/=3)  isp1=isp
          if(isp==3)  isp1=1
          if(lso==1) then      !P*P, H*P, H*H,
             ! WARN!  Except 001 case, we need to assume --phispinsym (radial functions are the same in both spins).
             !     We currently calculate only spin-diagonal Sz, and <up|L-|dn>, <dn|L+|up> (See text arount the
             !     end of augmat).
             !  Folloing construction of SSbPP is generally under the assumption of --phispinsym.

             ! SSb* is the atomic site contribution from ibas (augmentation parts. see m_bandcal_init->hambl->augmbl)
             f1=facso(1,isp); f2=facso(2,isp); f3=facso(3,isp)
             !                Lz                      L-                          L+
             SSbPP(:)= f1*Lzz(1)%sdiag(:,isp1) +f2*Lmp(1)%soffd(:,1)  + f3*Lmp(1)%soffd(:,2)
             SSbHP(:)= f1*Lzz(2)%sdiag(:,isp1) +f2*Lmp(2)%soffd(:,1)  + f3*Lmp(2)%soffd(:,2)
             SSbPH(:)= f1*dconjg(Lzz(2)%sdiag(:,isp1)) +f2*dconjg(Lmp(2)%soffd(:,2))+f3*dconjg(Lmp(2)%soffd(:,1))
             SSbHH(:)= f1*Lzz(3)%sdiag(:,isp1) +f2*Lmp(3)%soffd(:,1)  + f3*Lmp(3)%soffd(:,2)
          else
             fac = 1.5d0-isp
             SSbPP(:) = fac*Lzz(1)%sdiag(:,isp)
             SSbHP(:) = fac*Lzz(2)%sdiag(:,isp)
             SSbPH(:) = fac*dconjg(Lzz(2)%sdiag(:,isp))
             SSbHH(:) = fac*Lzz(3)%sdiag(:,isp)
          endif
          augq2zhso: block 
            integer::l1,ik1,i1,j,l2,ik2,i2,k,ilm1,jlm1,ilm2,iorb,jorb,k1
            complex(8):: &
                 hsohh(nkaph,nkaph,  nlmb,nlmb),& !HH
                 hsohp(nkaph,0:kmax, nlmb,nlma),&! HP
                 hsoph(nkaph,0:kmax, nlmb,nlma),&! PH (index ordering is transposed. the same as HP)
                 hsopp(0:kmax,0:kmax,nlma,nlma),& ! PP
                 g(0:kmax,nlma)
            !call augq2zhso(ibas,nkaph,lmxb,nlmb,kmax,nlma,b,ndimh, SSbHH,SSbHP,SSbPH,SSbPP, hso(:,:,isp))
            hsohh=reshape(SSbHH,shape(hsohh))
            hsohp=reshape(SSbHP,shape(hsohp))
            hsoph=reshape(SSbPH,shape(hsoph))
            hsopp=reshape(SSbPP,shape(hsopp))
            !write(6,*)'hhhhhhhsss',sum(hsohh),sum(hsohp),sum(hsoph),sum(hsopp)
            call orblib(ibas) !return norb,ltab,ktab,offl...
            do iorb = 1, norb
               l1  = ltab(iorb)
               ik1 = ktab(iorb)
               i1 = offl(iorb)
               do  ilm1 = l1**2+1, (l1+1)**2
                  i1 = i1+1
                  do  j = 1, ndimh
                     hso(i1,j,isp) = hso(i1,j,isp) + sum([(sum(hsohp(ik1,k,ilm1,:)*b(k,:,j))        ,k=0,kmax)])
                     hso(j,i1,isp) = hso(j,i1,isp) + sum([(sum(dconjg(b(k,:,j))*hsoph(ik1,k,ilm1,:)),k=0,kmax)])
                  enddo
                  do  jorb = 1, norb !one center
                     l2  = ltab(jorb)
                     ik2 = ktab(jorb)
                     i2 = offl(jorb)
                     do  ilm2 = l2**2+1, (l2+1)**2
                        i2 = i2+1
                        hso(i1,i2,isp) = hso(i1,i2,isp) + hsohh(ik1,ik2,ilm1,ilm2)
                     enddo
                  enddo
               enddo
            enddo
            do  i2 = 1, ndimh
               do  jlm1 = 1, nlma
                  g(0:kmax,jlm1) = [(sum(hsopp(k1,:,jlm1,:)*b(:,:,i2)),k1=0,kmax)]
               enddo
               do i1 = 1, ndimh
                  hso(i1,i2,isp) = hso(i1,i2,isp) + sum(dconjg(b(:,:,i1))*g(:,:))
               enddo
            enddo
          endblock augq2zhso
       enddo
       deallocate(b,SSbPP,SSbHP,SSbPH,SSbHH)
    enddo
    call tcx ('aughsoc')
  end subroutine aughsoc
end module m_augmbl



! subroutine augq2zhso(ia,nkaph,lmxb,nlmb,kmax, nlma,b,ndimh, hsohh,hsohp,hsoph,hsopp,  hso)
  !   use m_orbl,only: Orblib, norb,ltab,ktab,offl
  !   !- Add one and two-center terms to h,s for complex potential
  !   !o   hso   :1- and 2- center spin up-down spin orbit block
  !   !r Remarks
  !   !r  The hsopmhh(i,i,i,i,1), hsopmhh(i,i,i,i,2) are the head-head matrix
  !   !r  elements of LxSx+LySy. The ppihp(i,i,i,i,3), ppihp(i,i,i,i,4) are
  !   !r  the corresponding head-tail elements.
  !   ! takao  NOTE: LzSz is alreay added h by locpot-augmat-gaugm
  !   !r  The 2c term has the form h_{i,j} = Sum_kL(conjg(b_{i;k,L})*p_{j;k,L})+
  !   !r   Sum_kL(p_{i;k,L}*p_{j;k,L}); To get the second term for spin orbit
  !   !r   we rely on the hermicity of the ppi_{LxSx+LySy} block.
  !   !r   Symbolically:
  !   !r   hso_{i,j,u,d} =  Sum_kL[p_{i,j,u,d}*b_{j} + conjg(p_{j,i,d,u}*b_{i})]
  !   !r   where u = spin-up and d = spin-down.
  !   !r   If the structure constants become noncollinear, additional terms have
  !   !r   to be added in the matrix element above.
  !   ! ----------------------------------------------------------------------
  !   implicit none
  !   integer :: mode,ia,isp,kmax,nkaph,nlma,lmxb,nlmb,ndimh
  !   complex(8):: hsohh(nkaph,nkaph,nlmb,nlmb), & !HH
  !        hsohp(nkaph,0:kmax,nlmb,nlma),&! HP
  !        hsoph(nkaph,0:kmax,nlmb,nlma),&! PH (index ordering is transposed. the same as HP)
  !        hsopp(0:kmax,0:kmax,nlma,nlma) ! PP
  !   double complex b(0:kmax,nlma,ndimh), hso(ndimh,ndimh), g(0:kmax,nlma)
  !   integer :: iorb,ik1,j,k,ilma,i1,i2,ilm1,ilm2,l1,jorb,ik2,l2,jsp,ksp,isp1,isp2
  !   double precision :: xx
  !   integer :: jlm1,jlm2,k1,k2
  !   call tcn ('augq2zhs0')
  !   call orblib(ia) !return norb,ltab,ktab,offl...
  !   do iorb = 1, norb
  !      l1  = ltab(iorb)
  !      ik1 = ktab(iorb)
  !      i1 = offl(iorb)
  !      do  ilm1 = l1**2+1, (l1+1)**2
  !         i1 = i1+1
  !         do  j = 1, ndimh
  !            hso(i1,j) = hso(i1,j) + sum([(sum(hsohp(ik1,k,ilm1,:)*b(k,:,j))        ,k=0,kmax)])
  !            hso(j,i1) = hso(j,i1) + sum([(sum(dconjg(b(k,:,j))*hsoph(ik1,k,ilm1,:)),k=0,kmax)])
  !         enddo
  !         do  jorb = 1, norb !one center
  !            l2  = ltab(jorb)
  !            ik2 = ktab(jorb)
  !            i2 = offl(jorb)
  !            do  ilm2 = l2**2+1, (l2+1)**2
  !               i2 = i2+1
  !               hso(i1,i2) = hso(i1,i2) + hsohh(ik1,ik2,ilm1,ilm2)
  !            enddo
  !         enddo
  !      enddo
  !   enddo
  !   do  i2 = 1, ndimh
  !      do  jlm1 = 1, nlma
  !         g(0:kmax,jlm1) = [(sum(hsopp(k1,:,jlm1,:)*b(:,:,i2)),k1=0,kmax)]
  !      enddo
  !      do i1 = 1, ndimh
  !         hso(i1,i2) = hso(i1,i2) + sum(dconjg(b(:,:,i1))*g(:,:))
  !      enddo
  !   enddo
  !   call tcx ('augq2zhso')
  ! end subroutine augq2zhso

!       ppi1=>oppi(1,ibas)%cv !pi integral
!       ppi2=>oppi(2,ibas)%cv
!       ppi3=>oppi(3,ibas)%cv 
!       sig1=>osig(1,ibas)%v  !sigma integral
!       sig2=>osig(2,ibas)%v
!       sig3=>osig(3,ibas)%v
       !!  --- Add 1-center and 2-center terms ---
       !call augq2z(ibas,isp,nkaph,lmxb,nlmb,kmax,nlma, b,b,ndimh, sig3,sig2,ppi3,ppi2, s,h)
       !call augq3z(kmax,nlma,ndimh,isp,  b,                ppi1,           h) !B+ ppi B to h
       !call augqs3(kmax,lmxa,nlma,ndimh,isp, b,b, sig1,  s) !B+ sig B to s
       ! augq3z: block ! this works but very slow in gfortran 2023jan -----------
       !   complex(8):: ppi(0:kmax,0:kmax,nlma,nlma,nsp)
       !   ppi=reshape(ppi1,shape(ppi))
       !   do  i2 = 1, ndimh
       !      do i1 = 1, ndimh
       !         h(i1,i2) = h(i1,i2) +&
       !              sum([(        sum([(        sum([(          sum([( &
       !              dconjg(b(k1,ilm1,i1))*ppi(k1,k2,ilm1,ilm2,isp)*b(k2,ilm2,i2) &
       !              ,k1=0,kmax)]) ,k2=0,kmax)]) ,ilm2=1,nlma)]) ,ilm1=1,nlma)])
       !      enddo
       !   enddo
       ! endblock augq3z

  ! subroutine augq2z(ia,isp,nkaph,lmxb,nlmb,kmax,nlma,bL,bR,ndimh,sighh,sighp,ppihh,ppihp,s,h)
  !   use m_orbl,only: Orblib, norb,ltab,ktab,offl
  !   !- Add one and two-center terms to h,s for complex potential
  !   ! ----------------------------------------------------------------------
  !   !i Inputs
  !   !i   ia    : augmentation site about which strux are expanded
  !   !i   isp   : current spin channel
  !   !i   nkaph : dimensions augmentation matrices
  !   !i   nlmb : dimensions augmentation potential matrix at site a
  !   !i   lmxb : dimensions sighh at site a
  !   !i   kmax  :polynomial cutoff
  !   !i   nlma  :augmentation L-cutoff
  !   !i   sighh :augmentation head-head overlap matrix
  !   !i   ppihh :augmentation head-head potential matrix
  !   !i   sighp :augmentation head-Pkl overlap matrix
  !   !i   ppihp :augmentation head-Pkl potential matrix
  !   !i   b     :Bloch strux connecting site ia to all sites
  !   !i   ndimh :hamiltonian dimension
  !   !o Outputs
  !   !o   h     :1- and 2- center augmentation part of ham. added to h
  !   !o   s     :1- and 2- center augmentation part of ovlp added to s
  !   !r Remarks
  !   !r  In this implementation, the augmentation matrices and the row
  !   !r  dimension of the structure constants b follow normal L order.
  !   !r  The ppihh(i,i,i,i,3), ppihh(i,i,i,i,4) are the head-head matrix
  !   !r  elements of LxSx+LySy. The ppihp(i,i,i,i,3), ppihp(i,i,i,i,4) are
  !   !r  the corresponding head-tail elements.
  !   !r  The 2c term has the form h_{i,j} = Sum_kL(conjg(b_{i;k,L})*p_{j;k,L})+
  !   !r   Sum_kL(p_{i;k,L}*p_{j;k,L})
  !   ! ----------------------------------------------------------------------
  !   implicit none
  !   integer :: mode,ia,isp,kmax,nkaph,nlma,lmxb,nlmb,ndimh
  !   real(8)::  sighh(nkaph,nkaph,0:lmxb,1),sighp(nkaph,0:kmax,0:lmxb,1),xx
  !   complex(8):: ppihh(nkaph,nkaph,nlmb,nlmb,*), ppihp(nkaph,0:kmax,nlmb,nlma,*),&
  !        bL(0:kmax,nlma,ndimh),bR(0:kmax,nlma,ndimh), &
  !        s(ndimh,ndimh),h(ndimh,ndimh)
  !   integer :: iorb,ik1,j,k,ilma,i1,i2,ilm1,ilm2,l1,jorb,ik2,l2,jsp,ksp
  !   complex(8),allocatable:: tso(:,:,:,:)
  !   call tcn ('augq2z')
  !   call orblib(ia) !See use section. Return norb,ltab,ktab,offl
  !   do  iorb = 1, norb
  !      l1  = ltab(iorb)
  !      ik1 = ktab(iorb)
  !      do  ilm1 = l1**2+1, (l1+1)**2
  !         i1 = offl(iorb)+ilm1-l1**2  !Two-center terms
  !         s(i1,:) = s(i1,:) + [(       sum(sighp(ik1,:,l1,isp)*bR(:,ilm1,j)),  j=1,ndimh)]
  !         s(:,i1) = s(:,i1) + [(dconjg(sum(bL(:,ilm1,j)*sighp(ik1,:,l1,isp))), j=1,ndimh)]
  !         h(i1,:) = h(i1,:) + [(       sum(ppihp(ik1,:,ilm1,:,isp)*bR(:,:,j)), j=1,ndimh)]
  !         h(:,i1) = h(:,i1) + [(dconjg(sum(bL(:,:,j)*ppihp(ik1,:,ilm1,:,isp))),j=1,ndimh)]
  !      enddo   
  !      do  ilm1 = l1**2+1, (l1+1)**2
  !         i1 = offl(iorb)+ilm1-l1**2  !Two-center terms
  !         do  jorb = 1, norb !one center terms
  !            l2  = ltab(jorb)
  !            ik2 = ktab(jorb)
  !            do  ilm2 = l2**2+1, (l2+1)**2
  !               i2 = offl(jorb)+ilm2-l2**2
  !               h(i1,i2) = h(i1,i2) + ppihh(ik1,ik2,ilm1,ilm2,isp)
  !               if (ilm1 == ilm2) s(i1,i2) = s(i1,i2) + sighh(ik1,ik2,l1,isp)
  !            enddo
  !         enddo
  !      enddo
  !   enddo
  !   call tcx ('augq2z')
  ! end subroutine augq2z
  ! subroutine augq3z(kmax,nlma,ndimh,isp,b,ppi, h)
  !   implicit none 
  !   integer :: kmax,nlma,ndimh,isp!,mode1
  !   complex(8):: ppi(0:kmax,0:kmax,nlma,nlma,*), &
  !        b(0:kmax,nlma,ndimh),h(ndimh,ndimh), &
  !        g(0:kmax,nlma),csum,csum1 !,hso(ndimh,ndimh),gso(0:kmax,nlma)
  !   integer :: i1,i2,jlm1,jlm2,k1,k2,kjlm !,kjtop
  !   call tcn ('augq3z')
  !   do  i2 = 1, ndimh
  !      do  jlm1 = 1, nlma
  !         g(0:kmax,jlm1) = [(sum(ppi(k1,:,jlm1,:,isp)*b(:,:,i2)),k1=0,kmax)]
  !      enddo
  !      do i1 = 1, ndimh
  !         h(i1,i2) = h(i1,i2) + sum(dconjg(b(:,:,i1))*g(:,:))
  !      enddo
  !   enddo
  !   call tcx ('augq3z')
  ! end subroutine augq3z
  ! subroutine augqs3(kmax,lmxa,nlma,ndimh,isp,bL,bR,sig,s)
  !   !- Add B+ sig B to s for L-diagonal sig
  !   ! ----------------------------------------------------------------------
  !   !i Inputs
  !   !i   kmax  :polynomial cutoff
  !   !i   lmxa  :dimensions sig at site a
  !   !i   nlma  :augmentation L-cutoff
  !   !i   ndimh :hamiltonian dimension
  !   !i   isp   :current spin channel
  !   !i   sig   :augmentation Pkl-Pkl overlap matrix
  !   !i   b     :Bloch structure constants (hxpbl)
  !   !o Outputs
  !   !o   s     :overlap matrix
  !   ! ----------------------------------------------------------------------
  !   implicit none
  !   integer :: kmax,lmxa,nlma,ndimh,isp
  !   double precision :: sig(0:kmax,0:kmax,0:lmxa,isp)
  !   complex(8):: bL(0:kmax,nlma,ndimh),bR(0:kmax,nlma,ndimh),&
  !        s(ndimh,ndimh),g(0:kmax,nlma),csum
  !   integer :: i1,i2,ilm,k1,k2,l,kjlm, ll
  !   call tcn ('augqs3')
  !   do i2 = 1, ndimh
  !      do ilm = 1, nlma
  !         g(:,ilm) = matmul(sig(:,:,ll(ilm),isp),bR(:,ilm,i2))
  !      enddo
  !      s(1:i2,i2) = s(1:i2,i2) + [(sum( dconjg(bL(:,:,i1))*g(:,:) ),i1=1,i2)]
  !   enddo
  !   call tcx ('augqs3')
  ! end subroutine augqs3
