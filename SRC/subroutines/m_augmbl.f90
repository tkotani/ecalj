module m_augmbl !Add augmentation part to H and S. aughsoc add SO part to H.
  use m_lmfinit,only: lmxa_i=>lmxa,lmxb_i=>lmxb,kmxt_i=>kmxt
  !Inputs are Site integrals, sig,tau,pi,hso See JPSJ.kotani
  use m_ll,only:ll
  public augmbl,aughsoc
  private
contains
  subroutine augmbl(isp,q,osig,otau,oppi,ndimh, h,s)  !Add augmentation part to H and S. 
    use m_lmfinit,only: nsp,nlmto
    use m_lmfinit,only: nbas,alat=>lat_alat,ispec
    use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol,rv_a_opos
    use m_bstrux,only: Bstrux_set, bstr
    use m_orbl,only: Orblib, norb,ltab,ktab,offl
    use m_struc_def,only: s_cv1,s_rv1,s_rv4,s_cv5
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
    type(s_cv5),target :: oppi(3,nbas)
    type(s_rv4),target :: otau(3,nbas)
    type(s_rv4),target :: osig(3,nbas)
    integer:: isp,ndimh,napw, initbas, endbas,lm,iq,nh,np,i1,i2,ilm1,ilm2,k1,k2, ibas,isa,kmax,lmxa,lmxb, nglob,nlma, nlmb
    real(8):: q(3),rsma,pa(3),xx 
    complex(8):: h(ndimh,ndimh),s(ndimh,ndimh)
    integer,parameter :: ktop0=20, nlmbx=49, nlmax=49
    logical:: debug=.false.
    call tcn ('augmbl')
    do ibas = 1,nbas
       isa = ispec(ibas) 
       pa  = rv_a_opos(:,ibas) 
       lmxa= lmxa_i(isa) !max l of augmentation
       lmxb= lmxb_i(isa) !max l of basis
       kmax= kmxt_i(isa) !max of radial k
       nlmb = (lmxb+1)**2
       nlma = (lmxa+1)**2
       if (lmxa == -1) cycle
       call bstrux_set(ibas,q) !Get structure constant bstr for given ibas and q for expanding orbitals at site ia.(Bloch sum of C^i_akL (C.1) in Ref.[1])
       addaug: block
         integer:: i1,i2,ilm,l1,ik1,l2,ik2,ilm1,ilm2,j,iorb,jorb
         complex(8):: b(0:kmax,nlma,ndimh),g(0:kmax,nlma)
         !complex(8)::sig(0:kmax,0:kmax,0:lmxa,nsp),ppi(0:kmax,0:kmax,nlma,nlma,nsp), ppihh(nkaph,nkaph,nlmb,nlmb,nsp), ppihp(nkaph,0:kmax,nlmb,nlma,nsp)
         !real(8):: sighh(nkaph,nkaph,0:lmxb,nsp),    sighp(nkaph,0:kmax,0:lmxb,nsp)
         b = reshape(bstr,shape=shape(b),order=[3,2,1])
         associate(&
              sig  =>osig(1,ibas)%v, & !sig= reshape(osig(1,ibas)%v,shape(sig))
              sighp=>osig(2,ibas)%v, & !reshape(osig(2,ibas)%v,shape(sighp)) !augmentation head-Pkl overlap matrix  
              sighh=>osig(3,ibas)%v, & !reshape(,shape(sighh)) !augmentation head-head overlap matrix 
              ppi  =>oppi(1,ibas)%cv,& ! ppi=reshape(oppi(1,ibas)%cv,shape(ppi))
              ppihp=>oppi(2,ibas)%cv,& !reshape(oppi(2,ibas)%cv,shape(ppihp)) !augmentation head-Pkl potential matrix
              ppihh=>oppi(3,ibas)%cv)  !reshape(oppi(3,ibas)%cv,shape(ppihh)) !augmentation head-head potential matrix
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
         ! do  i2 = 1, ndimh
         !    do  ilm = 1, nlma
         !       g(0:kmax,ilm) = [(sum(ppi(k1,:,ilm,:,isp)*b(:,:,i2)),k1=0,kmax)]
         !    enddo
         !    h(:,i2) = h(:,i2) + [(sum(dconjg(b(:,:,i1))*g(:,:)),i1=1,ndimh)]
         ! enddo
         ! do i2 = 1, ndimh
         !    do ilm = 1, nlma
         !       g(:,ilm) = matmul(sig(:,:,ll(ilm),isp),b(:,ilm,i2))
         !    enddo
         !    s(1:i2,i2) = s(1:i2,i2) + [(sum( dconjg(b(:,:,i1))*g(:,:) ),i1=1,i2)]
         ! enddo
         ! MO The above two loops were placed in the following blas_mode block 2024-11-07
         blas_mode: block
           use m_blas, only: zmm => zmm_h, m_op_C
           complex(8), allocatable :: ppib(:,:,:), ppi_isp(:,:,:,:), sigb_ilm(:,:),b_ilm(:,:), csig(:,:)
           integer :: istat
           allocate(ppib(0:kmax,nlma,ndimh), ppi_isp(0:kmax,nlma,0:kmax,nlma))
           ppi_isp = reshape(source=ppi(0:kmax,0:kmax,1:nlma,1:nlma,isp), shape=[kmax+1,nlma,kmax+1,nlma], order=[1,3,2,4])
           istat = zmm(ppi_isp, b, ppib, m=(kmax+1)*nlma, n=ndimh, k=(kmax+1)*nlma)
           istat = zmm(b, ppib, h, m=ndimh, n=ndimh, k=(kmax+1)*nlma, beta=(1d0,0d0), opA=m_op_C)
           deallocate(ppib, ppi_isp)
           allocate(sigb_ilm(0:kmax,ndimh), b_ilm(0:kmax,ndimh), csig(0:kmax,0:kmax))
           do ilm=1,nlma
             b_ilm(0:kmax,1:ndimh) = b(0:kmax,ilm,1:ndimh)
             csig(0:kmax,0:kmax) = sig(0:kmax,0:kmax,ll(ilm),isp) !convert to complex from real
             istat = zmm(csig, b_ilm, sigb_ilm, m=kmax+1, n=ndimh, k=kmax+1)
             istat = zmm(b_ilm, sigb_ilm, s, m=ndimh, n=ndimh, k=kmax+1, beta=(1d0,0d0), opA=m_op_C)
           enddo
           deallocate(sigb_ilm, b_ilm, csig)
         endblock blas_mode
         endassociate
       endblock addaug
    enddo
    call tcx ('augmbl')
  end subroutine augmbl
  subroutine aughsoc(qp,ohsozz,ohsopm,ndimh, hso) ! Spin-orbit-couping matrix hso
    use m_orbl,only: Orblib, norb,ltab,ktab,offl
    use m_struc_def,only: s_cv1,s_rv1,s_sblock
    use m_lmfinit,only: nsp, lsox=>lso, nbas, nkaphh, ispec, socaxis
    use m_bstrux,only: Bstrux_set, bstr
    use m_lattic,only: plat=>lat_plat,qlat=>lat_qlat
    !i   qp    :Bloch wave number
    !i   hsozz,hsopm : atomic parts of SOC (Lz and L-)
    !i   ndimh :dimension of halimtonian.
    !o   hso   :spin diagonal and off-diagonal block of spin-orbit hamiltonian
    ! note  'shorbz need to be improved in future (the method in shortn3)'.
    ! note  We obtain Lz,L+,and L- (Lzz Lmm Lpp) in this routine. From their linear combinatios, we have hso.
    implicit none
    type(s_sblock),target :: ohsozz(3,nbas),ohsopm(3,nbas)
    integer:: isp, ndimh, ibas, isa,kmax,lmxa,lmxb, nglob,nlma,nlmb,lso,nkaph
    integer:: initbas, endbas,lm,iq,nh,np,isp1,isp2,nspx
    real(8):: q(3),qp(3),fac
    complex(8):: hso(ndimh,ndimh,3)
    complex(8),allocatable:: b(:,:,:)
    complex(8),pointer:: ppi1(:),ppi2(:),ppi3(:)
    type(s_sblock),pointer:: Lzz(:),Lmp(:)
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
          call rx('Given HAM_SOCAXIS is not yet implemented. Modify facso matrix in subrouitne aughsoc in m_augmbl.f90.')
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
       lmxa=lmxa_i(isa) !max l of augmentation
       lmxb=lmxb_i(isa) !max l of basis
       kmax=kmxt_i(isa) !max of radial k
       nlmb = (lmxb+1)**2
       nlma = (lmxa+1)**2
       if (lmxa == -1) cycle
       call bstrux_set(ibas,q) !Make strux b to expand all orbitals at site ia
       if(allocated(b)) deallocate(b)
!       write(6,*)'kkkkkkkkkkkkk',kmax,nlma,ndimh
       allocate( b(0:kmax,nlma,ndimh) )
       b = reshape(bstr,shape(b),order=[3,2,1])
       nkaph=nkaphh(isa)
       nh= nkaph*nlmb     ! size of head nh
       np= (kmax+1)*nlma  ! size of tail np
       !! Get Lzz,Lmp,Lmp(spinfliped)= (Lz,L-,L+)  See mkpot-locpot-augmat-gaugm-pvagm1,pvaglc to generate hsozz,hsopm
       Lzz => ohsozz(:,ibas)   ! Lz block  1:P*P, 2:H*P, 3:H*H for up and dn
       if(lso==1) Lmp => ohsopm(:,ibas) ! <up|L-|dn> for isp=1 , <dn|L+|up> for isp=2. See gaugm.F, pvagm1,pvaglc
       nspx=nsp
       if(lso==1) nspx=3
       do isp=1,nspx
          ! hso(:,:,isp=1) is (1,1) block in (A8) in in Liqin2019,PhysRevB.99.054418
          ! hso(:,:,isp=2) is (1,1) block
          ! hso(:,:,isp=3) is (1,2) block  =L- in the case of 001 spin axis.
          if(isp/=3)  isp1=isp
          if(isp==3)  isp1=1
          augq2zhso: block 
            integer::l1,ik1,i1,j,l2,ik2,i2,k,ilm1,jlm1,ilm2,iorb,jorb,k1
            complex(8):: &
                 hsohh(nkaph,nkaph,  nlmb,nlmb),& !HH
                 hsohp(nkaph,0:kmax, nlmb,nlma),&! HP
                 hsoph(nkaph,0:kmax, nlmb,nlma),&! PH (index ordering is transposed. the same as HP)
                 hsopp(0:kmax,0:kmax,nlma,nlma),&! PP
                 g(0:kmax,nlma)
            if(lso==1) then      !P*P, H*P, H*H,
               ! WARN!  Except 001 case, we need to assume --phispinsym (radial functions are the same in both spins).
               !     We currently calculate only spin-diagonal Sz, and <up|L-|dn>, <dn|L+|up> (See text arount the
               !     end of augmat).
               !  Folloing hsofoobar, we assume --phispinsym.
               f1=facso(1,isp); f2=facso(2,isp); f3=facso(3,isp)
               !                Lz                      L-                          L+
               hsopp= f1*Lzz(1)%sdiag(:,:,:,:,isp1)         +f2*Lmp(1)%soffd(:,:,:,:,1)        + f3*Lmp(1)%soffd(:,:,:,:,2) 
               hsohp= f1*Lzz(2)%sdiag(:,:,:,:,isp1)         +f2*Lmp(2)%soffd(:,:,:,:,1)        + f3*Lmp(2)%soffd(:,:,:,:,2)
               hsoph= f1*dconjg(Lzz(2)%sdiag(:,:,:,:,isp1)) +f2*dconjg(Lmp(2)%soffd(:,:,:,:,2))+ f3*dconjg(Lmp(2)%soffd(:,:,:,:,1))
               hsohh= f1*Lzz(3)%sdiag(:,:,:,:,isp1)         +f2*Lmp(3)%soffd(:,:,:,:,1)        + f3*Lmp(3)%soffd(:,:,:,:,2) 
            else
               fac = 1.5d0-isp
               hsopp= fac*Lzz(1)%sdiag(:,:,:,:,isp)
               hsohp= fac*Lzz(2)%sdiag(:,:,:,:,isp)
               hsoph= fac*dconjg(Lzz(2)%sdiag(:,:,:,:,isp))
               hsohh= fac*Lzz(3)%sdiag(:,:,:,:,isp) 
            endif
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
       deallocate(b)
    enddo
    call tcx ('aughsoc')
  end subroutine aughsoc
end module m_augmbl
