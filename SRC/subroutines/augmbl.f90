subroutine aughsoc(qp,ohsozz,ohsopm, ndimh, hso)
  use m_struc_def,only: s_cv1,s_rv1,s_sblock
  use m_lmfinit,only: nsp, lsox=>lso, nbas, nkaph, ssite=>v_ssite, sspec=>v_sspec,socaxis
  use m_bstrux,only: Bstrux_set, bstr
  use m_lattic,only: plat=>lat_plat,qlat=>lat_qlat
  !!- Spin-orbit couping matrix hso
  !i   qp    :Bloch wave number
  !i   hsozz,hsopm : atomic parts of SOC (Lz and L-)
  !i   ndimh :dimension of halimtonian.
  !o   hso   :spin diagonal and off-diagonal block of spin-orbit hamiltonian
  ! note  qp is shortened to be q as in the same manner of hambl.
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
     !     if(init) then
     !     write(6,*)'fff lso=',lso
     !     write(6,"('fff facso 1=',6(2f9.4,2x))")facso(:,1)
     !     write(6,"('fff facso 2=',6(2f9.4,2x))")facso(:,2)
     !     write(6,"('fff facso 3=',6(2f9.4,2x))")facso(:,3)
     !     init=.false.
     !     endif
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
  !sss call shorbz(qp,q,qlat,plat) !is this fine?
  q=qp !sss 
  do ibas = 1,nbas
     isa =ssite(ibas)%spec
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
     !        print *,' np nh nsp=',np,nh,nsp,np*np*nsp
     !        allocate(SSbPP(np*np*nsp),SSbHP(np*nh*nsp),SSbPH(np*nh*nsp),SSbHH(nh*nh*nsp))
     allocate(SSbPP(np*np),SSbHP(np*nh),SSbPH(np*nh),SSbHH(nh*nh))
     nspx=nsp
     if(lso==1) nspx=3
     do isp=1,nspx
        ! hso(:,:,isp=1) is (1,1) block in (A8) in in Liqin2019,PhysRevB.99.054418
        ! hso(:,:,isp=2) is (1,1) block
        ! hso(:,:,isp=3) is (1,2) block  =L- in the case of 001 spin axis.
        if(isp/=3)  isp1=isp
        if(isp==3)  isp1=1
        !           if(isp1==1) isp2=2
        !           if(isp1==2) isp2=1
        if(lso==1) then      !P*P, H*P, H*H,
           ! WARN!  Except 001 case, we need to assume --phispinsym (radial functions are the same in both spins).
           !     We currently calculate only spin-diagonal Sz, and <up|L-|dn>, <dn|L+|up> (See text arount the
           !     end of augmat).
           !  Folloing construction of SSbPP is generally under the assumption of --phispinsym.

           ! SSb* is the atomic site contribution from ibas (augmentation parts. see m_bandcal_init->hambl->augmbl)
           f1=facso(1,isp); f2=facso(2,isp); f3=facso(3,isp)
           !                                Lz                      L-                          L+
           SSbPP(:)= f1*Lzz(1)%sdiag(:,isp1) +f2*Lmp(1)%soffd(:,1)  + f3*Lmp(1)%soffd(:,2)
           SSbHP(:)= f1*Lzz(2)%sdiag(:,isp1) +f2*Lmp(2)%soffd(:,1)  + f3*Lmp(2)%soffd(:,2)
           SSbPH(:)= f1*dconjg(Lzz(2)%sdiag(:,isp1)) &
                +     f2*dconjg(Lmp(2)%soffd(:,2)) + f3*dconjg(Lmp(2)%soffd(:,1))
           SSbHH(:)= f1*Lzz(3)%sdiag(:,isp1) +f2*Lmp(3)%soffd(:,1)  + f3*Lmp(3)%soffd(:,2)
        else
           fac = 1.5d0-isp
           SSbPP(:) = fac*Lzz(1)%sdiag(:,isp)
           SSbHP(:) = fac*Lzz(2)%sdiag(:,isp)
           SSbPH(:) = fac*dconjg(Lzz(2)%sdiag(:,isp))
           SSbHH(:) = fac*Lzz(3)%sdiag(:,isp)
        endif
        call augq2zhso(ibas,nkaph,lmxb,nlmb,kmax,nlma,b,ndimh, SSbHH,SSbHP,SSbPH,SSbPP, hso(:,:,isp))
     enddo
     deallocate(b,SSbPP,SSbHP,SSbPH,SSbHH)
  enddo
  call tcx ('aughsoc')
end subroutine aughsoc

subroutine lmlps(n,din,dout)
  integer:: n
  complex(8):: din(n,2),dout(n,2)
  dout(:,2) = din(:,1)
  dout(:,1) = din(:,2)
end subroutine lmlps

subroutine augmbl(ssite,sspec,isp, &
     q , sv_p_osig , sv_p_otau , sv_p_oppi, ndimh , h , s )  !,ohsozz,ohsopm
  use m_lmfinit,only: nsp,nlmto!,lso
  use m_struc_def
  use m_lmfinit,only: rv_a_ocg , iv_a_oidxcg , iv_a_ojcg , rv_a_ocy
  use m_lmfinit,only: nbas,nkaph,lat_alat
  use m_lattic,only: lat_qlat, lat_vol,rv_a_opos
  use m_bstrux,only: Bstrux_set, bstr
  ! this is used for lso=0 only now (aug2021), but lso=2 should work.
  !- Adds augmentation part of H and S
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode1 :  0 do not compute hso
  ! x         :  1 compute hso.
  ! x         :    Note: only a portion of hso is computed for a
  ! x         :    particular isp.  The total hso is assembled
  ! x         :    after isp loops from 1..2.  hso should not be
  ! x         :    initialized between isp=1 and isp=2 loops.
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
  ! x   hso   :spin off-diagonal block of spin-orbit hamiltonian
  !o   s     :augmentation part of overlap matrix added to s
  !l Local variables
  !l   nkaph :number of orbital types for a given L quantum no. in basis
  !l         :at augmentation site ia, including local orbitals
  !l   nlmto :number of lmto basis functions
  !r Remarks
  !r   Some expressions labelled JMP refer to J.Math.Phys39, 3393 (1998)
  !b Bugs
  !b   Not really a bug, but an inefficiency:
  !b   Right now, strux are kept for all orbitals in the basis, including
  !b   expansions coffs for local orbitals (which are set to zero).
  !b   Better to condense strux to reduce computational effort for 2-
  !b   and 3-center terms.
  !m MPI
  !m   See remarks in hsibl. Buffers for h and s are taken from the heap.
  !m   In addition a buffer the same size as as h and s for ALLREDUCE.
  !u Updates
  !u   05 Jul 08 (T. Kotani) output density for new PW part
  !u             Option to accumulate energy-weighted output density
  !u   08 Sep 06 (WRL) updated MPI to work with SO coupling
  !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
  !u   03 Feb 05 (A. Chantis) calculate hso
  !u    1 Sep 04 Adapted to handle complex ppi.  S.O. folded into ppi
  !u   25 Aug 04 Adapted to extended local orbitals
  !u   29 Jun 04 (A. Chantis) Include LzSz spin-orbit coupling
  !u   14 Aug 02 Added overlap-only option
  !u   15 Feb 02 (ATP) Added MPI parallelization
  !u   11 Jan 02 Adapted for f90 constructs
  !u   27 Aug 01 Extended to local orbitals.
  !u   17 Jun 00 spin polarized
  !u   18 May 00 Adapted from nfp augm_q.f
  !u   1998      (DLN) parallel version for SGI
  ! ----------------------------------------------------------------------
  implicit none
  type(s_cv1),target :: sv_p_oppi(3,nbas) !, ohsozz(3,nbas),ohsopm(3,nbas)
  type(s_rv1) :: sv_p_otau(3,nbas)
  type(s_rv1) :: sv_p_osig(3,nbas)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  integer:: isp , ndimh , napw ,i_copy_size,numprocs !lcplxp ,
  real(8):: q(3)
  double complex h(ndimh,ndimh),s(ndimh,ndimh)!,hso(ndimh,ndimh)
  integer :: nlmbx,nlmax,ktop0
  parameter (ktop0=20, nlmbx=49, nlmax=49)
  complex(8),allocatable:: b(:,:,:),bx(:,:,:),bb(:,:,:)
  integer:: ibas , isa , kmax , lmxa , lmxb ,  nglob , nlma &
       , nlmb !,   lso
  double precision :: rsma,pa(3),xx,alat,qlat(3,3),vol
  integer:: initbas, endbas,lm,iq,nh,np
  logical:: debug=.false.
  complex(8):: zxx(1)
  complex(8),pointer:: ppi1(:),ppi2(:),ppi3(:),Lm1(:),Lm2(:),Lm3(:),Lz1(:),Lz2(:),Lz3(:)
  !--------------------------
  call tcn ('augmbl')
  alat=lat_alat
  qlat=lat_qlat
  vol=lat_vol
  do ibas = 1,nbas
     isa =ssite(ibas)%spec
     pa  =rv_a_opos(:,ibas) !ssite(ibas)%pos
     lmxa=sspec(isa)%lmxa !max l of augmentation
     lmxb=sspec(isa)%lmxb !max l of basis
     kmax=sspec(isa)%kmxt !max of radial k
!     rsma=sspec(isa)%rsma
     nlmb = (lmxb+1)**2
     nlma = (lmxa+1)**2
     if (lmxa == -1) cycle
     !   --- Make strux to expand all orbitals at site ia ---
     call bstrux_set(ibas,q)
     if(allocated(b)) deallocate(b)
     allocate( b(0:kmax,nlma,ndimh) )
     do lm=1,nlma
        b(:,lm,:) = transpose(bstr(:,lm,:))
     enddo
     !! Add LzSz contribution, spin diagonal part
     !        if(lso==0) then
     ppi3=>sv_p_oppi(3,ibas)%cv !pi integral
     ppi2=>sv_p_oppi(2,ibas)%cv
     ppi1=>sv_p_oppi(1,ibas)%cv
     !        else                    !Add LzSz contribution to pi-integral ppi (moved from gaugm).
     !           call rx('augmbl is only for lso=0')
     !           nh=nkaph*nlmb
     !           np=(kmax+1)*nlma
     !           Lz1 => ohsozz(1,ibas)%cv ! Lz block  P*P
     !           Lz2 => ohsozz(2,ibas)%cv !           H*P
     !           Lz3 => ohsozz(3,ibas)%cv !           H*H
     !           allocate(ppi3(nh**2*nsp),ppi2(nh*np*nsp),ppi1(np**2*nsp))
     !           call hsozzadd(np**2, sv_p_oppi(1,ibas)%cv,Lz1,nsp,isp, ppi1)
     !           call hsozzadd(nh*np, sv_p_oppi(2,ibas)%cv,Lz2,nsp,isp, ppi2)
     !           call hsozzadd(nh**2, sv_p_oppi(3,ibas)%cv,Lz3,nsp,isp, ppi3)
     !        endif
     !!  --- Add 1-center and 2-center terms ---
     call augq2z(ibas , isp , nkaph , lmxb, nlmb , kmax , nlma ,  b, ndimh, &
          sv_p_osig(3,ibas)%v, sv_p_osig(2,ibas)%v, ppi3,ppi2,         s,h)
     call augq3z(kmax,nlma,ndimh,isp,  b,                ppi1,           h) !B+ ppi B to h
     call augqs3( kmax, lmxa, nlma, ndimh, isp, b, sv_p_osig(1,ibas)%v,  s) !B+ sig B to s
     deallocate(b)
  enddo
  !        if(lso/=0) deallocate(ppi1,ppi2,ppi3)
  call tcx ('augmbl')
end subroutine augmbl


subroutine augq2z(ia,isp,nkaph,lmxb,nlmb,kmax, &
     nlma,b,ndimh,sighh,sighp,ppihh,ppihp,s,h)
  use m_orbl,only: Orblib, norb,ltab,ktab,offl
  !- Add one and two-center terms to h,s for complex potential
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   do not compute hso
  !i   ia    :augmentation site about which strux are expanded
  !i   isp   :current spin channel
  !i   nkaph :dimensions augmentation matrices
  !i   nlmb :dimensions augmentation potential matrix at site a
  !i   lmxb :dimensions sighh at site a
  !i   kmax  :polynomial cutoff
  !i   nlma  :augmentation L-cutoff
  !i   sighh :augmentation head-head overlap matrix
  !i   ppihh :augmentation head-head potential matrix
  !i   sighp :augmentation head-Pkl overlap matrix
  !i   ppihp :augmentation head-Pkl potential matrix
  !i   b     :Bloch strux connecting site ia to all sites
  !i   ndimh :hamiltonian dimension
  !o Outputs
  !o   h     :1- and 2- center augmentation part of ham. added to h
  !o   s     :1- and 2- center augmentation part of ovlp added to s
  !r Remarks
  !r  In this implementation, the augmentation matrices and the row
  !r  dimension of the structure constants b follow normal L order.
  !r  The ppihh(i,i,i,i,3), ppihh(i,i,i,i,4) are the head-head matrix
  !r  elements of LxSx+LySy. The ppihp(i,i,i,i,3), ppihp(i,i,i,i,4) are
  !r  the corresponding head-tail elements.

  !r  The 2c term has the form h_{i,j} = Sum_kL(conjg(b_{i;k,L})*p_{j;k,L})+
  !r   Sum_kL(p_{i;k,L}*p_{j;k,L})
  ! ----------------------------------------------------------------------
  implicit none
  integer :: mode,ia,isp,kmax,nkaph,nlma,lmxb,nlmb,ndimh
  double precision :: &
       sighh(nkaph,nkaph,0:lmxb,1),sighp(nkaph,0:kmax,0:lmxb,1)
  complex(8):: ppihh(nkaph,nkaph,nlmb,nlmb,*), ppihp(nkaph,0:kmax,nlmb,nlma,*),&
       b(0:kmax,nlma,ndimh),s(ndimh,ndimh),h(ndimh,ndimh)
  ! ... Local parameters
  integer :: iorb,ik1,j,k,ilma,i1,i2,ilm1,ilm2,l1,n0,nkap0,jorb,ik2,l2,jsp,ksp
  parameter (n0=10,nkap0=3)
  !      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
  double precision :: xx
  double complex cadd,cadd1
  complex(8),allocatable:: tso(:,:,:,:)
  call tcn ('augq2z')
  ! --- Loop over basis functions at site ia (augentation index) ---
  call orblib(ia) !See use section. Return norb,ltab,ktab,offl
  do  iorb = 1, norb
     l1  = ltab(iorb)
     ik1 = ktab(iorb)
     i1 = offl(iorb)
     do  ilm1 = l1**2+1, (l1+1)**2
        i1 = i1+1
        do  j = 1, ndimh !Two-center terms
           do  k = 0, kmax
              cadd = sighp(ik1,k,l1,isp)*b(k,ilm1,j)
              s(i1,j) = s(i1,j) + cadd
              s(j,i1) = s(j,i1) + dconjg(cadd)
              do  ilma = 1, nlma
                 cadd = ppihp(ik1,k,ilm1,ilma,isp)*b(k,ilma,j)
                 h(i1,j) = h(i1,j) + cadd
                 h(j,i1) = h(j,i1) + dconjg(cadd)
              enddo
           enddo
        enddo
        do  jorb = 1, norb !one center terms
           l2  = ltab(jorb)
           ik2 = ktab(jorb)
           i2 = offl(jorb)
           do  ilm2 = l2**2+1, (l2+1)**2
              i2 = i2+1
              h(i1,i2) = h(i1,i2) + ppihh(ik1,ik2,ilm1,ilm2,isp)
              if (ilm1 == ilm2) s(i1,i2) = s(i1,i2) + sighh(ik1,ik2,l1,isp)
           enddo
        enddo
     enddo
  enddo
  call tcx ('augq2z')
end subroutine augq2z
! sssssssssssssssssssssssss
subroutine augq2zhso(ia,nkaph,lmxb,nlmb,kmax, nlma,b,ndimh, hsohh,hsohp,hsoph,hsopp,  hso)
  use m_orbl,only: Orblib, norb,ltab,ktab,offl
  !- Add one and two-center terms to h,s for complex potential
  !o   hso   :1- and 2- center spin up-down spin orbit block
  !r Remarks
  !r  The hsopmhh(i,i,i,i,1), hsopmhh(i,i,i,i,2) are the head-head matrix
  !r  elements of LxSx+LySy. The ppihp(i,i,i,i,3), ppihp(i,i,i,i,4) are
  !r  the corresponding head-tail elements.

  ! takao  NOTE: LzSz is alreay added h by locpot-augmat-gaugm

  !r  The 2c term has the form h_{i,j} = Sum_kL(conjg(b_{i;k,L})*p_{j;k,L})+
  !r   Sum_kL(p_{i;k,L}*p_{j;k,L}); To get the second term for spin orbit
  !r   we rely on the hermicity of the ppi_{LxSx+LySy} block.
  !r   Symbolically:
  !r   hso_{i,j,u,d} =  Sum_kL[p_{i,j,u,d}*b_{j} + conjg(p_{j,i,d,u}*b_{i})]
  !r   where u = spin-up and d = spin-down.
  !r   If the structure constants become noncollinear, additional terms have
  !r   to be added in the matrix element above.
  ! ----------------------------------------------------------------------
  implicit none
  integer :: mode,ia,isp,kmax,nkaph,nlma,lmxb,nlmb,ndimh
  complex(8):: &
       hsohh(nkaph,nkaph,nlmb,nlmb), & !HH
  hsohp(nkaph,0:kmax,nlmb,nlma),&! HP
  hsoph(nkaph,0:kmax,nlmb,nlma),&! PH (index ordering is transposed. the same as HP)
  hsopp(0:kmax,0:kmax,nlma,nlma) ! PP
  double complex b(0:kmax,nlma,ndimh), hso(ndimh,ndimh), g(0:kmax,nlma)
  ! ... Local parameters
  integer :: iorb,ik1,j,k,ilma,i1,i2,ilm1,ilm2,l1,n0,nkap0,jorb,ik2,l2,jsp,ksp,isp1,isp2
  parameter (n0=10,nkap0=3)
  !      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
  double precision :: xx
  double complex cadd,cadd1,fac
  integer :: jlm1,jlm2,k1,k2
  call tcn ('augq2zhs0')
  call orblib(ia) !return norb,ltab,ktab,offl...
  do iorb = 1, norb
     l1  = ltab(iorb)
     ik1 = ktab(iorb)
     i1 = offl(iorb)
     do  ilm1 = l1**2+1, (l1+1)**2
        i1 = i1+1
        do  j = 1, ndimh
           do  k = 0, kmax
              do  ilma = 1, nlma
                 hso(i1,j) = hso(i1,j) + hsohp(ik1,k,ilm1,ilma)*b(k,ilma,j)
                 hso(j,i1) = hso(j,i1) + dconjg(b(k,ilma,j))*hsoph(ik1,k,ilm1,ilma)
              enddo
           enddo
        enddo
        do  jorb = 1, norb !one center
           l2  = ltab(jorb)
           ik2 = ktab(jorb)
           i2 = offl(jorb)
           do  ilm2 = l2**2+1, (l2+1)**2
              i2 = i2+1
              hso(i1,i2) = hso(i1,i2) + hsohh(ik1,ik2,ilm1,ilm2)
           enddo
        enddo
     enddo
  enddo
  do  i2 = 1, ndimh
     g = 0d0
     do  jlm1 = 1, nlma
        do  k1 = 0, kmax
           g(k1,jlm1) = g(k1,jlm1) + sum(hsopp(k1,:,jlm1,:)*b(:,:,i2))
        enddo
     enddo
     do i1 = 1, ndimh
        hso(i1,i2) = hso(i1,i2) + sum(dconjg(b(:,:,i1))*g(:,:))
     enddo
  enddo
  call tcx ('augq2z')
end subroutine augq2zhso


! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine augq3z(kmax,nlma,ndimh,isp,b,ppi, h)
  implicit none
  integer :: kmax,nlma,ndimh,isp!,mode1
  complex(8):: ppi(0:kmax,0:kmax,nlma,nlma,*), &
       b(0:kmax,nlma,ndimh),h(ndimh,ndimh), &
       g(0:kmax,nlma),csum,csum1 !,hso(ndimh,ndimh),gso(0:kmax,nlma)
  integer :: i1,i2,jlm1,jlm2,k1,k2,kjlm !,kjtop
  call tcn ('augq3z')
  do  i2 = 1, ndimh
     g = 0d0
     do  jlm1 = 1, nlma
        do  k1 = 0, kmax
           g(k1,jlm1) = g(k1,jlm1) + sum(ppi(k1,:,jlm1,:,isp)*b(:,:,i2))
        enddo
     enddo
     do i1 = 1, ndimh
        h(i1,i2) = h(i1,i2) + sum(dconjg(b(:,:,i1))*g(:,:))
     enddo
  enddo
  call tcx ('augq3z')
end subroutine augq3z
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine augqs3(kmax,lmxa,nlma,ndimh,isp,b,sig,s)
  !- Add B+ sig B to s for L-diagonal sig
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   kmax  :polynomial cutoff
  !i   lmxa  :dimensions sig at site a
  !i   nlma  :augmentation L-cutoff
  !i   ndimh :hamiltonian dimension
  !i   isp   :current spin channel
  !i   sig   :augmentation Pkl-Pkl overlap matrix
  !i   b     :Bloch structure constants (hxpbl)
  !o Outputs
  !o   s     :overlap matrix
  ! ----------------------------------------------------------------------
  implicit none
  integer :: kmax,lmxa,nlma,ndimh,isp
  double precision :: sig(0:kmax,0:kmax,0:lmxa,isp)
  double complex b(0:kmax,nlma,ndimh),s(ndimh,ndimh),g(0:kmax,nlma),csum
  integer :: i1,i2,ilm,k1,k2,l,kjlm
  integer :: ll
  call tcn ('augqs3')
  do i2 = 1, ndimh
     do ilm = 1, nlma
        l = ll(ilm)
        g(:,ilm) = matmul(sig(:,:,l,isp),b(:,ilm,i2))
     enddo
     do  i1 = 1, i2
        s(i1,i2) = s(i1,i2) + sum( dconjg(b(:,:,i1))*g(:,:) )
     enddo
  enddo
  call tcx ('augqs3')
end subroutine augqs3


subroutine hsozzadd(nx,ppi,hsozz,nsp,isp,ppiz)
  integer:: nx,isp,nsp
  complex(8):: ppi(nx,nsp),hsozz(nx,nsp),ppiz(nx,nsp)
  ppiz(:,isp) = ppi(:,isp)+(1.5d0-isp)*hsozz(:,isp)
end subroutine hsozzadd
subroutine hsozzadd2(nx,ppi,hsozz,nsp,ppiz)
  integer:: nx,isp,nsp
  complex(8):: ppi(nx,nsp),hsozz(nx,nsp),ppiz(nx,nsp)
  do isp=1,nsp
     ppiz(:,isp) = ppi(:,isp)+(1.5d0-isp)*hsozz(:,isp)
  enddo
end subroutine hsozzadd2

!$$$!sssssssssssssssssssssssss
!$$$      subroutine augq3zhso(kmax,nlma,ndimh,ppi,b,hso)
!$$$C- Add B+ ppi B to H for non-L-diagonal, complex matrix ppi
!$$$Co   hso   :3-center from this augmentation site added to hso
!$$$C ----------------------------------------------------------------------
!$$$      implicit none
!$$$      integer kmax,nlma,ndimh,isp!,mode1
!$$$      complex(8):: ppi(0:kmax,0:kmax,nlma,nlma,3),
!$$$     & b(0:kmax,nlma,ndimh),!h(ndimh,ndimh),
!$$$     & csum,gso(0:kmax,nlma),csum1,hso(ndimh,ndimh)
!$$$      integer i1,i2,jlm1,jlm2,k1,k2,kjlm !,kjtop
!$$$!  ...  Make LxSx+LySy part of SO
!$$$      do i2 = 1, ndimh
!$$$        gso=0d0
!$$$        do  jlm1 = 1, nlma
!$$$        do  k1 = 0, kmax
!$$$            gso(k1,jlm1) = gso(k1,jlm1) + sum(ppi(k1,:,jlm1,:,3)*b(:,:,i2))
!$$$        enddo
!$$$        enddo
!$$$        do i1 = 1, ndimh
!$$$           hso(i1,i2) = hso(i1,i2) + .5d0*sum(dconjg(b(:,:,i1))*gso(:,:))
!$$$        enddo
!$$$      enddo
!$$$      end subroutine augq3zhso
