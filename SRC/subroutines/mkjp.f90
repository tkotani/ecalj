module  m_vcoulq
  use m_cputm,only:cputm
  use m_mpi,only: ipr,mpi__rank
  use m_lgunit,only: stdo
  use m_ftox
  public vcoulq_4,mkjb_4,mkjp_4,genjh, ajr, a1r, vcoul_termb
  private
  character(1024):: aaaw
  ! Tables of one group of atoms (same radial mesh and lx), made by mkjp_4 for the first atom of the group
  ! (hasBessel=F) and used by the other atoms of the group and by vcoul_termb, which must run before the next group
  ! replaces them.  ajr(r,ig,l) = j_l(|q+G| r) r/|q+G|^l; for eee/=0, a1r(r,ig,l) of sigkernel, which includes the
  ! Simpson weights times dr/di (fac_integral), so sum_r a1r(r,ig,l) f(r) is already the radial integral.
  ! In the GPU build only the device copies are set.
  real(8), allocatable  :: ajr(:,:,:), a1r(:,:,:)
contains
  subroutine vcoulq_4(q,nbloch,ngc,nbas,lx,lxx,nx,nxx,alat,qlat,vol,ngvecc, & !Coulmb matrix for each q
       strx,rojp,rojb,sgbb,sgpb,fouvb,ngb,bas,rmax, eee, nr,nrx,rofi, vcoulb,   vcoul)
    use m_ll,only: ll
    use m_blas, only: m_op_T,m_op_C
#ifdef __GPU
    use m_blas, only: zmm=>zmm_d
#else
    use m_blas, only: zmm=>zmm_h
#endif
    !i strx:  Structure factors
    !i nlx corresponds to (lx+1)**2 . lx corresponds to 2*lmxax.
    !i rho-type integral
    !i  ngvecc     : q+G vector
    !i  rojp rojb  : rho-type integral
    !i  sigma-type onsite integral
    !i  nx(l,ibas) : max number of radial function index for each l and ibas.
    !i               Note that the definition is a bit different from nx in basnfp.
    !i  nxx        : max number of nx among all l and ibas.
    !i  lx(nbas)   : max number of l for each ibas.
    !o vcoul:  in a.u. You have to multiply by 2 for vcoul in Ry.
    !---------------------------------------------------------------------------
    !  rojp = <j_aL(r) | P(q+G)_aL > where
    !         |P(q+G)_aL> : the projection of exp(i (q+G) r) to aL channnel.
    !         |j_aL>      : \def r^l/(2l+1)!! Y_L.  The spherical bessel functions near r=0.  Energy-dependence is omitted.
    implicit none
    integer :: nbloch, ngb, nbas, lxx,lx(nbas), nxx, nx(0:lxx,nbas)
    integer :: ibl1, ibl2,ig1,ig2,ibas,ibas1,ibas2, l,m,n, n1,l1,m1,lm1,n2,l2,m2,lm2,ipl1,ipl2
    integer :: ibasbl(nbloch), nbl(nbloch), lbl(nbloch), mbl(nbloch), lmbl(nbloch)
    integer :: ngc, ngvecc(3,ngc),nrx,nr(nbas),lm
    real(8) :: vol,q(3), qlat(3,3),alat,absqg2(ngc),qg(3), rojb(nxx, 0:lxx, nbas)
    real(8) :: sgbb(nxx,  nxx,  0:lxx,      nbas) !i sigma-type onsite integral
    real(8) :: fpivol,tpiba, bas(3,nbas),r2s,rmax(nbas)
    real(8) ::  fkk(0:lxx),fkj(0:lxx),fjk(0:lxx),fjj(0:lxx),sigx(0:lxx),radsig(0:lxx) !,radsig(0:lxx,nbas),fjj(0:lxx,nbas)
    real(8) :: eee, rofi(nrx,nbas)
    real(8),allocatable  :: cy(:),yl(:)
    real(8),parameter:: pi=4d0*datan(1d0),fpi=4d0*pi
    complex(8) :: rojp(ngc, (lxx+1)**2, nbas)   !rho-type onsite integral
    complex(8) :: strx((lxx+1)**2, nbas, (lxx+1)**2,nbas) !structure constant. The multicenter expantion of 1/|r-r'|
    complex(8) :: sgpb(ngc,  nxx,  (lxx+1)**2, nbas)
    complex(8) :: fouvb(ngc,  nxx, (lxx+1)**2, nbas) ,vcoul(ngb, ngb) !<exp(i q+G r)|xxx>
    complex(8) :: vcoulb(*)   !eee/=0: the onsite parts of <P_G1|v|P_G2> of all groups (vcoul_termb), on the device;
                              !packed lower triangle, vcoulb(ig1*(ig1-1)/2+ig2) for ig2<=ig1
    complex(8),allocatable :: pjyl_(:,:),phase(:,:)
    complex(8) :: img=(0d0,1d0)
    integer :: istat,lm2x
    integer,allocatable :: llx(:)
    write(aaaw,'(" vcoulq_4: ngb  nbloch ngc nrx procid=",5i6)') ngb,nbloch,ngc,nrx,mpi__rank
    call cputm(stdo,aaaw)
    fpivol = 4*pi*vol
    allocate( pjyl_((lxx+1)**2,ngc),phase(ngc,nbas),source=(0d0,0d0) )
    allocate( cy((lxx+1)**2), yl((lxx+1)**2),source=0d0)
    allocate( llx((lxx+1)**2),source=0)
    do lm =1,(lxx+1)**2
      llx(lm) = ll(lm)
    enddo
    call sylmnc(cy,lxx)
    tpiba = 2*pi/alat
    do ig1 = 1,ngc
      qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvecc(1:3,ig1))) !q+G in a.u.
      absqg2(ig1)  = sum(qg(1:3)**2)+1d-32
      phase(ig1,:) = exp( img*matmul(qg(1:3),bas(1:3,:))*alat  )
      call sylm(qg/sqrt(absqg2(ig1)),yl,lxx,r2s) !spherical factor Y( q+G )
      do lm =1,(lxx+1)**2
        l = ll(lm)
        pjyl_(lm,ig1) = fpi*img**l *cy(lm)*yl(lm)  * sqrt(absqg2(ig1))**l  ! <jlyl | exp i q+G r> projection of exp(i q+G r) to jl yl on MT
      enddo
    enddo
    ibl1 = 0
    do ibas= 1, nbas
      do l   = 0, lx(ibas) !-- index (mx,nx,lx,ibas) order.
        do n   = 1, nx(l,ibas)
          do m   = -l, l
            ibl1  = ibl1 + 1
            ibasbl(ibl1) = ibas
            nbl   (ibl1) = n
            lbl   (ibl1) = l
            mbl   (ibl1) = m
            lmbl  (ibl1) = l**2 + l+1 +m
          enddo
        enddo
      enddo
    enddo
    if(ibl1/=nbloch) call rx(' vcoulq: error ibl1/=nbloch', ibl1, nbloch)
    !$acc enter data copyin(rojb) create(vcoul) copyin(ibasbl, nbl, lbl, mbl, lmbl)   ! strx, rojp, sgpb, fouvb: on the device (hvccfp0)
    !$acc kernels
    vcoul(:,:) = 0d0
    !$acc end kernels
    !-- <B|v|B> block
    !$acc kernels loop independent collapse(2)
    BvB: do ibl2= 1, nbloch
      do ibl1= 1, nbloch
        ibas1= ibasbl(ibl1)
        n1   = nbl (ibl1)
        l1   = lbl (ibl1)
        m1   = mbl (ibl1)
        lm1  = lmbl(ibl1)
        ibas2= ibasbl(ibl2)
        n2   = nbl (ibl2)
        l2   = lbl (ibl2)
        m2   = mbl (ibl2)
        lm2  = lmbl(ibl2)
        vcoul(ibl1,ibl2) = rojb(n1, l1, ibas1) *strx(lm1,ibas1,lm2,ibas2) *rojb(n2, l2, ibas2)    ! offsite Coulomb
        if(ibas1==ibas2.AND.lm1==lm2) vcoul(ibl1,ibl2) = vcoul(ibl1,ibl2) + sgbb(n1,n2,l1, ibas1) ! sigma-type onsite parts
      enddo
    enddo BvB
    !$acc end kernels

    ! <P_G|v|B>
    PvB_dev_mo:block
      write(aaaw,ftox)' vcoulq_4: goto PvB procid=', mpi__rank
      call cputm(stdo,aaaw)
      PvB2: block
        complex(8) :: strxx(1:(lxx+1)**2,1:nbas,nbloch)
        complex(8) :: crojp_ibas(ngc,(lxx+1)**2,nbas)
        !$acc data create(strxx, crojp_ibas)
        !$acc kernels
        crojp_ibas(1:ngc,1:(lxx+1)**2,1:nbas) = dconjg(rojp(1:ngc,1:(lxx+1)**2, 1:nbas))
        !$acc end kernels
        !$acc kernels loop present(strx, rojb, nbl, lbl, ibasbl, lmbl)
        do ibl2=1,nbloch
          strxx(:,:,ibl2)= -strx(:,:,lmbl(ibl2),ibasbl(ibl2)) * rojb(nbl (ibl2),lbl (ibl2),ibasbl(ibl2))
        enddo
        !$acc end kernels
        !$acc host_data use_device(vcoul)
        istat = zmm(crojp_ibas, strxx,  vcoul(nbloch+1,1), m=ngc, n=nbloch, k=(lxx+1)**2*nbas,LdC=ngb) !not LdC is needed
        !$acc end host_data
        !$acc end data
      endblock PvB2
      !$acc data present(fouvb, sgpb)
      !$acc kernels loop present(vcoul, ibasbl, nbl, lbl, lmbl)
      PvB: do ibl2= 1, nbloch
        ibas2= ibasbl(ibl2)
        n2   = nbl (ibl2)
        l2   = lbl (ibl2)
        lm2  = lmbl(ibl2)
        !m2   = mbl (ibl2)
        vcoul(     nbloch+1:nbloch+ngc,ibl2) =&
             vcoul(nbloch+1:nbloch+ngc,ibl2) &
             +fouvb(1:ngc, n2, lm2, ibas2) - sgpb(1:ngc, n2, lm2, ibas2)   !<exp(i(q+G)r)|v|B_n2L2> !punch out onsite part
      enddo PvB
      !$acc end kernels
      !$acc end data
    endblock PvB_dev_mo
    !$acc exit data delete(rojb, ibasbl, nbl, lbl, mbl, lmbl)

    ! <P_G|v|P_G>
    PvP_dev_mo: block
      complex(8) :: rojpstrx((lxx+1)**2,nbas,ngc)
      integer :: ibas_order(nbas), isrt, jsrt, ktmp, itype_start, itype_end, ib_next
      complex(8) :: cPhi
      complex(8), allocatable :: vcoul_termA(:,:)
      write(aaaw,ftox) " vcoulq_4: goto PvP procid ngc lxx nrx=", mpi__rank,ngc,lxx,nrx
      call cputm(stdo,aaaw)

      lm2x= (lxx+1)**2

      ! Sort atoms by (nr, lx) to maximize Bessel/wronkj reuse
      do isrt = 1, nbas; ibas_order(isrt) = isrt; enddo
      do isrt = 1, nbas-1
        do jsrt = 1, nbas-isrt
          if(nr(ibas_order(jsrt)) > nr(ibas_order(jsrt+1)) .or. &
             (nr(ibas_order(jsrt)) == nr(ibas_order(jsrt+1)) .and. lx(ibas_order(jsrt)) > lx(ibas_order(jsrt+1)))) then
            ktmp = ibas_order(jsrt); ibas_order(jsrt) = ibas_order(jsrt+1); ibas_order(jsrt+1) = ktmp
          endif
        enddo
      enddo

      !$acc data create(rojpstrx) copyin(absqg2)

      !$acc host_data use_device(strx, rojp)
      istat = zmm(strx, rojp, rojpstrx, m=nbas*(lxx+1)**2, n=ngc, k=nbas*(lxx+1)**2, opA=m_op_T, opB=m_op_C)
      !$acc end host_data

      ! --- Term A: sum over all atoms via single BLAS call ---
      ! vcoul_A(ig1,ig2) = sum_{lm,ibas} rojpstrx(lm,ibas,ig1)*rojp(ig2,lm,ibas)
      allocate(vcoul_termA(ngc, ngc))
      !$acc data create(vcoul_termA)
      !$acc host_data use_device(rojpstrx, rojp, vcoul_termA)
      istat = zmm(rojpstrx, rojp, vcoul_termA, m=ngc, n=ngc, k=lm2x*nbas, opA=m_op_T, opB=m_op_T)
      !$acc end host_data
      !$acc kernels
      do ig1 = 1, ngc
        do ig2 = 1, ig1
          vcoul(nbloch+ig1, nbloch+ig2) = vcoul(nbloch+ig1, nbloch+ig2) + vcoul_termA(ig1, ig2)
        enddo
      enddo
      !$acc end kernels
      !$acc end data
      deallocate(vcoul_termA)

      if(eee/=0d0) then   ! Term B: the onsite parts of each group, made by vcoul_termb
        write(aaaw,ftox) " vcoulq_4: add the onsite parts (vcoul_termb)", mpi__rank
        call cputm(stdo,aaaw)
        !$acc parallel loop gang vector collapse(2) present(vcoulb(1:(ngc*(ngc+1))/2))
        do ig1 = 1, ngc
          do ig2 = 1, ngc
            if(ig2 <= ig1) vcoul(nbloch+ig1,nbloch+ig2) = vcoul(nbloch+ig1,nbloch+ig2) + vcoulb((ig1*(ig1-1))/2+ig2)
          enddo
        enddo
      else   ! eee=0: Term B on the host, once per group of atoms (same nr, lx and rofi; "type" below)
        write(aaaw,ftox) " vcoulq_4: goto igig loop (type-batched)", mpi__rank
        call cputm(stdo,aaaw)
        itype_start = 1
        do while(itype_start <= nbas)
          ibas = ibas_order(itype_start)
          ! Find end of this atom type (same nr, lx, rofi)
          itype_end = itype_start
          do while(itype_end < nbas)
            ib_next = ibas_order(itype_end + 1)
            if(nr(ib_next) /= nr(ibas) .or. lx(ib_next) /= lx(ibas)) exit
            if(.not. all(abs(rofi(1:nr(ibas),ibas) - rofi(1:nr(ib_next),ib_next)) < 1d-10)) exit
            itype_end = itype_end + 1
          enddo
          write(aaaw,ftox) " vcoulq_4: type atoms", itype_start, '-', itype_end, 'nr=', nr(ibas), 'lx=', lx(ibas), 'procid=', mpi__rank
          call cputm(stdo,aaaw)
          ! CPU path: wronkj + sigintpp once per type, multiply by Phi_type
          !$acc update self(vcoul)
          do ig1 = 1, ngc
            do ig2 = 1, ig1
              fjj(0:lxx) = 0d0   ! wronkj sets 0:lx; the sum below runs to lxx
              call wronkj( absqg2(ig1), absqg2(ig2), rmax(ibas), lx(ibas), fkk, fkj, fjk, fjj)
              call sigintpp( absqg2(ig1)**.5d0, absqg2(ig2)**.5d0, lx(ibas), rmax(ibas), sigx)
              radsig(0:lxx) = 0d0
              forall(l = 0:lx(ibas)) radsig(l) = fpi/(2*l+1) * sigx(l)
              cPhi = (0d0, 0d0)
              do jsrt = itype_start, itype_end
                cPhi = cPhi + dconjg(phase(ig1, ibas_order(jsrt))) * phase(ig2, ibas_order(jsrt))
              enddo
              vcoul(nbloch+ig1,nbloch+ig2) = vcoul(nbloch+ig1,nbloch+ig2) &
                + cPhi * sum( dconjg(pjyl_(1:lm2x,ig1)) * pjyl_(1:lm2x,ig2) &
                  * ((fpi/(absqg2(ig1)-eee)+fpi/(absqg2(ig2)-eee))*fjj(llx(1:lm2x)) + radsig(llx(1:lm2x))) )
            enddo
          enddo
          !$acc update device(vcoul)
          itype_start = itype_end + 1
        enddo ! type loop
      endif

      !$acc kernels
      do ig1 = 1, ngc
        vcoul(nbloch + ig1, nbloch + ig1) = vcoul(nbloch + ig1, nbloch + ig1) + fpivol/(absqg2(ig1) - eee) !eee is negative
      end do
      !$acc end kernels

      !$acc end data
    endblock PvP_dev_mo

    !$acc exit data copyout(vcoul)

    RightUpperPartOFvcoul: do ipl1=1, nbloch+ngc
      do ipl2=1, ipl1-1
        vcoul(ipl2,ipl1) = dconjg(vcoul(ipl1,ipl2))
      enddo
    enddo RightUpperPartOFvcoul
!    do ix = 1,nbloch+ngc
!      if((mod(ix,20)==1 .OR. ix>nbloch+ngc-10).and.ipr) write(6,"(' Diagonal Vcoul =',i5,2d18.10)") ix,vcoul(ix,ix)
!    enddo
#ifdef __GPU
    GPUmemoryRelease: block
      use openacc
      call acc_clear_freelists()
    end block GPUmemoryRelease
#endif    
  end subroutine vcoulq_4

  
  ! ptest=.False.
  ! PlaneWavetest: if(ptest) then !check Coulomb by plane wave expansion.
  !   if(ipr) write(6,*) ' --- plane wave Coulomb matrix check 1---- '
  !   write(197,*) ' --- off diagonal ---- '
  !   nblochngc = nbloch+ngc
  !   allocate(matp(nblochngc),matp2(nblochngc))
  !   do ig1 = 1,ngc
  !     matp = 0d0
  !     do ibl2= 1, nbloch
  !       ibas2= ibasbl(ibl2)
  !       n2   = nbl (ibl2)
  !       l2   = lbl (ibl2)
  !       m2   = mbl (ibl2)
  !       lm2  = lmbl(ibl2)
  !       matp(ibl2) = fouvb(ig1, n2, lm2, ibas2)*absqg2(ig1)/fpi
  !     enddo
  !     matp(nbloch+ig1) = 1d0
  !     ig2=ig1
  !     !      do ig2 = 1,ngc !off diagnal
  !     matp2 = 0d0
  !     do ibl2= 1, nbloch
  !       ibas2= ibasbl(ibl2)
  !       n2   = nbl (ibl2)
  !       l2   = lbl (ibl2)
  !       m2   = mbl (ibl2)
  !       lm2  = lmbl(ibl2)
  !       matp2(ibl2) = fouvb(ig2, n2, lm2, ibas2)*absqg2(ig2)/fpi
  !     enddo
  !     matp2(nbloch+ig2) = 1d0
  !     xxx= sum( matmul(matp(1:nblochngc),vcoul(1:nblochngc,1:nblochngc)) *dconjg(matp2(1:nblochngc))  )
  !     if(ig1/=ig2) then  !off diagnal
  !       if(abs(xxx)>1d-1 ) then
  !         write(197,'(2i5, 2d13.6)') ig1,ig2, xxx
  !         write(197,'("    matpp ", 2d13.6)') vcoul(nbloch+ig1,nbloch+ig2)
  !         write(197,*)
  !       endif
  !     else
  !       write(196,'(2i5," exact=",3d13.6,"q ngsum=",3f8.4,i5)') &
  !            ig1,ig2,fpi*vol/absqg2(ig1) , fpi*vol/absqg2(ig2),absqg2(ig1), q(1:3) , sum(ngvecc(1:3,ig1)**2)
  !       write(196,'("           cal  =", 2d13.6)') xxx
  !       write(196,'("           vcoud=", 2d13.6)') vcoul(nbloch+ig1,nbloch+ig2)
  !       write(196,*)
  !     endif
  !   enddo
  !   deallocate(matp,matp2)
  ! endif PlaneWavetest
  !end subroutine vcoulq_4

  subroutine mkjp_4(q,ngc,ngvecc,alat,qlat,lxx,lx,nxx,nx,bas,a,b,rmax,nr,nrx,rprodx,eee,rofi,rkpr,rkmr, rojp,sgpb,fouvb,hasBessel)! Integrals@MT and fouvb
    ! The integrals rojp, fouvb,fouvp are for  J_L(r)= j_l(sqrt(e) r)/sqrt(e)**l Y_L, which behaves as r^l/(2l+1)!! near r=0.
    ! oniste integral is based on 1/|r-r'| = \sum 4 pi /(2k+1) \frac{r_<^k }{ r_>^{k+1} } Y_L(r) Y_L(r')
    ! See PRB34 5512(1986) for sigma type integral
    ! hasBessel=F: first atom of a group (atoms with the same radial mesh and lx), which makes the tables ajr and a1r
    ! used by the next atoms of the group and by vcoul_termb.
    use m_ll,only: ll
    use m_bessl, only: bessl2 => bessl, wronkj2 => wronkj
    implicit none
    integer:: ngc,ngvecc(3,ngc), lxx, lx, nxx,nx(0:lxx),nr,nrx, nlx,ig1,l,n,ir,lm
    real(8):: q(3),bas(3), rprodx(nrx,nxx,0:lxx),a,b,rmax,alat, qlat(3,3)
    real(8):: pi,fpi,tpiba, qg1(3), fkk(0:lx),fkj(0:lx),fjk(0:lx),fjj(0:lx),absqg1, &
         phi(0:lx),psi(0:lx),r2s,sig
    real(8):: rofi(nrx),rkpr(nrx,0:lxx),rkmr(nrx,0:lxx),eee
    real(8),allocatable::cy(:),yl(:)
    real(8),allocatable ::a1(:,:,:), qg(:,:),absqg(:), rofi_nr(:), fac_integral(:)
    complex(8) :: rojp(ngc, (lxx+1)**2)        ! rho-type onsite integral
    complex(8) :: sgpb(ngc,  nxx,  (lxx+1)**2) !sigma-type onsite integral
    complex(8) :: fouvb(ngc,  nxx, (lxx+1)**2)
    complex(8) :: img =(0d0,1d0),phase
    complex(8),allocatable :: pjyl(:,:)
    logical, intent(in) :: hasBessel
    nlx = (lx+1)**2
    allocate(qg(3,ngc),absqg(ngc), pjyl((lx+1)**2,ngc) )

    pi    = 4d0*datan(1d0)
    fpi   = 4*pi
    tpiba = 2*pi/alat
    allocate(cy((lx+1)**2),yl((lx+1)**2))
    call sylmnc(cy,lx)
    !... q+G and <J_L | exp(i q+G r)>  J_L= j_l/sqrt(e)**l Y_L
    do ig1 = 1,ngc
      qg(1:3,ig1) = tpiba * (q(1:3)+ matmul(qlat, ngvecc(1:3,ig1)))
      qg1(1:3) = qg(1:3,ig1)
      absqg(ig1)  = sqrt(sum(qg1(1:3)**2))
      absqg1   = absqg(ig1) +1d-32
      phase = exp( img*sum(qg1(1:3)*bas(1:3))*alat  )
      call sylm(qg1/absqg1,yl,lx,r2s) !spherical factor Y( q+G )
      do lm =1,nlx
        l = ll(lm)
        pjyl(lm,ig1) = fpi*img**l *cy(lm)*yl(lm) *phase  *absqg1**l ! <jlyl | exp i q+G r> projection of exp(i q+G r) to jl yl  on MT
      enddo
    enddo 
    write(aaaw,ftox)' mkjp_4: goto rojploop'
    call cputm(stdo,aaaw)

    allocate(rofi_nr, source = rofi(1:nr))
    allocate(fac_integral(nr))
    do ir = 1, nr    ! Simpson weights times dr/di
      fac_integral(ir) = a*b*dexp(a*(ir-1))/3d0
      if(ir /= 1 .and. ir /= nr) fac_integral(ir) = fac_integral(ir)*merge(4d0,2d0,mod(ir,2)==0)
    enddo
    if(eee==0d0) then
      allocate(a1(1:nr,0:lx,ngc))
      do ig1 = 1,ngc
        call sigintAn1( absqg(ig1), lx, rofi_nr, nr,a1(1:nr, 0:lx,ig1) )
      enddo
    endif
    write(aaaw,ftox)' mkjp_4: goto dev_mo block. nx', nx(:)
    call cputm(stdo,aaaw)
    dev_mo: block
#ifdef __GPU
      use m_blas, only: dmm => dmm_d, m_op_T
#else
      use m_blas, only: dmm => dmm_h, m_op_T
#endif
      real(8), allocatable :: sigg(:,:,:), radintg(:,:,:), rprodw(:,:,:)
      integer :: llist(nlx), istat

      llist(1:nlx) = [(ll(lm), lm=1, nlx)]
      ! rojp, sgpb and fouvb stay on the device for vcoulq_4; rkpr, rkmr and rprodx are there for all q (hvccfp0)
      !$acc data copyin(absqg, rofi_nr, fac_integral, pjyl, llist, nx) present(rkpr, rkmr, rprodx, rojp, sgpb, fouvb)
      !$acc parallel loop gang vector private(fkk(0:lx), fkj(0:lx), fjk(0:lx), fjj(0:lx))
      rojploop: do ig1 = 1, ngc
        call wronkj2( absqg(ig1)**2, eee, rmax,lx, fkk,fkj,fjk,fjj)
        do lm = 1, (lxx+1)**2
          if(lm <= nlx) then
            rojp(ig1,lm) = (-fjj(llist(lm)))* pjyl(lm,ig1)
          else
            rojp(ig1,lm) = 0d0
          endif
        enddo
      enddo rojploop
      setTables: if(.not.hasBessel) then ! tables of the group of this atom (else the previous atom has the same mesh and lx)
        if(allocated(ajr)) then
          !$acc exit data delete(ajr)
          deallocate(ajr)
        endif
        allocate(ajr(1:nr,ngc,0:lx))
        !$acc enter data create(ajr)
        !$acc parallel loop collapse(2) private(phi(0:lx), psi(0:lx))
        do ig1 = 1, ngc
          do ir = 1, nr
            call bessl2(absqg(ig1)**2*rofi_nr(ir)**2,lx,phi,psi)
            do l = 0, lx
              ajr(ir,ig1,l) = phi(l)* rofi_nr(ir) **(l +1 )  ! ajr = j_l(sqrt(e) r) * r / (sqrt(e))**l
              !  Sperical Bessel j_l(r) \propto r**l/ (2l+1)!! near r=0.
            enddo
          enddo
        enddo
        !$acc end parallel
        if(eee/=0d0) then
          if(allocated(a1r)) then
            !$acc exit data delete(a1r)
            deallocate(a1r)
          endif
          allocate(a1r(1:nr,ngc,0:lx))
          !$acc enter data create(a1r)
          !$acc parallel loop gang vector collapse(2) present(ajr, a1r)
          do l = 0, lx
            do ig1 = 1, ngc
              call sigkernel(nr, a, b, rofi_nr, rkpr(1,l), rkmr(1,l), ajr(1,ig1,l), fac_integral, a1r(1,ig1,l))
            enddo
          enddo
        endif
      endif setTables

      if(eee==0d0) then
        do lm = 1, nlx
          l = llist(lm)
          do n = 1, nx(l)      
            do ig1 = 1,ngc
              call gintxx(a1(1,l,ig1),rprodx(1,n,l),A,B,NR, sig )
              sgpb(ig1,n,lm) = dconjg(pjyl(lm,ig1))* sig/(2*l+1)*fpi
            enddo
          enddo
        enddo
        !$acc update device(sgpb)
      else
        allocate(sigg(ngc,nxx,0:lx))
        !$acc data create(sigg)
        !$acc host_data use_device(a1r, rprodx, sigg)
        do l = 0, lx
          if(nx(l) == 0) cycle
          istat = dmm(a1r(1,1,l), rprodx(1,1,l), sigg(1,1,l), m=ngc, n=nx(l), k=nr, opA=m_op_T, ldB=nrx)
        enddo
        !$acc end host_data
        !$acc kernels loop independent collapse(2)
        do lm = 1, nlx
          do n = 1, nxx
            l = llist(lm)
            if(n > nx(l)) cycle
            sgpb(1:ngc,n,lm) = dconjg(pjyl(lm,1:ngc))* sigg(1:ngc,n,l)/(2*l+1)*fpi
          enddo
        enddo
        !$acc end kernels
        !$acc end data
        deallocate(sigg)
      endif
      allocate(radintg(ngc,nxx,0:lx), rprodw(nr,nxx,0:lx))
      !$acc data create(radintg, rprodw)
      !$acc parallel loop collapse(3)
      do l = 0, lx
        do n = 1, nxx
          do ir = 1, nr
            rprodw(ir,n,l) = rprodx(ir,n,l)*fac_integral(ir)
          enddo
        enddo
      enddo
      !$acc host_data use_device(ajr, rprodw, radintg)
      do l = 0, lx
        if(nx(l) == 0) cycle
        istat = dmm(ajr(1,1,l), rprodw(1,1,l), radintg(1,1,l), m=ngc, n=nx(l), k=nr, opA=m_op_T)
      enddo
      !$acc end host_data
      !$acc kernels
      fouvb(:,:,:) = 0d0
      !$acc end kernels
      !$acc kernels loop independent collapse(2)
      do lm = 1, nlx
        do n = 1, nxx
          l = llist(lm)
          if(n > nx(l)) cycle
          fouvb(1:ngc, n, lm) = fpi/(absqg(1:ngc)**2-eee) *dconjg(pjyl(lm,1:ngc))*radintg(1:ngc,n,l)
        enddo
      enddo
      !$acc end kernels
      !$acc end data
      deallocate(radintg, rprodw)

      !$acc end data
    endblock dev_mo

    deallocate(absqg, qg, pjyl, cy, yl, rofi_nr, fac_integral)
  end subroutine mkjp_4
  subroutine vcoul_termb(q,ngc,ngvecc,alat,qlat,lxg,natg,basg,rmaxg,nr,eee, vcoulb) ! Onsite part of <P_G1|v|P_G2> of one group of atoms
    ! vcoulb(igg) += cPhi(G1,G2) sum_{l<=lxg} sum_m conj(pjyl(lm,G1)) pjyl(lm,G2)
    !                  * [ (4pi/(|q+G1|^2-e)+4pi/(|q+G2|^2-e)) fjj_l + 4pi/(2l+1) sigx_l ],   igg = ig1*(ig1-1)/2+ig2, ig2<=ig1,
    ! cPhi = sum_a conj(exp(i(q+G1)R_a)) exp(i(q+G2)R_a) over the atoms basg(:,1:natg) of the group (the same radial mesh
    ! and lx=lxg).  fjj_l is fjj of wronkj from the Bessel values and slopes at rmax of each G; sigx_l = a1r^T ajr of the
    ! tables mkjp_4 made for the group.  One l at a time, so the work memory is one ngc x ngc matrix (sx).  The sum
    ! stops at lxg, as rojp (zero for l>lx) and the eee=0 path of vcoulq_4 do.
    use m_ll,only: ll
    use m_bessl, only: radkj2 => radkj
#ifdef __GPU
    use m_blas, only: dmm => dmm_d, m_op_T
#else
    use m_blas, only: dmm => dmm_h, m_op_T
#endif
    implicit none
    integer,intent(in):: ngc, ngvecc(3,ngc), lxg, natg, nr
    real(8),intent(in):: q(3), alat, qlat(3,3), basg(3,natg), rmaxg, eee
    complex(8):: vcoulb((ngc*(ngc+1))/2)
    real(8),parameter:: pi=4d0*datan(1d0), fpi=4d0*pi
    complex(8),parameter:: img=(0d0,1d0)
    integer:: ig, ig1, ig2, igg, l, m, lm, ia, istat, nlmg
    real(8):: tpiba, qg(3), r2s, e1w, e2w, rw, rjw, fjj
    real(8):: akw(lxg+2), ajw(lxg+2), dkw(lxg+2), djw(lxg+2)
    real(8),allocatable:: absqg2(:), cy(:), yl(:), ajg(:,:), djg(:,:), aje(:,:), dje(:,:), sx(:,:)
    complex(8),allocatable:: pjyl(:,:), phase(:,:)
    complex(8):: pp, cphi
    nlmg = (lxg+1)**2
    allocate(absqg2(ngc), cy(nlmg), yl(nlmg), pjyl(nlmg,ngc), phase(ngc,natg))
    allocate(ajg(0:lxg,ngc), djg(0:lxg,ngc), aje(0:lxg,ngc), dje(0:lxg,ngc), sx(ngc,ngc))
    call sylmnc(cy,lxg)
    tpiba = 2*pi/alat
    do ig = 1, ngc  ! q+G, the phase of each atom and pjyl as in vcoulq_4; Bessel values and slopes at rmax
      qg(1:3) = tpiba*(q(1:3) + matmul(qlat, ngvecc(1:3,ig)))
      absqg2(ig) = sum(qg(1:3)**2)+1d-32
      phase(ig,1:natg) = exp(img*matmul(qg(1:3),basg(1:3,1:natg))*alat)
      call sylm(qg/sqrt(absqg2(ig)),yl,lxg,r2s)
      do lm = 1, nlmg
        l = ll(lm)
        pjyl(lm,ig) = fpi*img**l*cy(lm)*yl(lm)*sqrt(absqg2(ig))**l
      enddo
      call radkj2(absqg2(ig), rmaxg, lxg, akw, ajw, dkw, djw, 0)
      ajg(0:lxg,ig) = ajw(1:lxg+1)
      djg(0:lxg,ig) = djw(1:lxg+1)
      aje(0:lxg,ig) = 0d0
      dje(0:lxg,ig) = 0d0
      if (dabs(absqg2(ig)) > 1d-6) then   ! job 1 divides by e; used only for equal nonzero energies
        call radkj2(absqg2(ig), rmaxg, lxg, akw, ajw, dkw, djw, 1)
        aje(0:lxg,ig) = ajw(1:lxg+1)
        dje(0:lxg,ig) = djw(1:lxg+1)
      endif
    enddo
    rw = rmaxg
    !$acc data copyin(absqg2, pjyl, phase, ajg, djg, aje, dje) create(sx) present(vcoulb, ajr, a1r)
    do l = 0, lxg
      !$acc host_data use_device(a1r, ajr, sx)
      istat = dmm(a1r(1,1,l), ajr(1,1,l), sx, m=ngc, n=ngc, k=nr, opA=m_op_T)   ! sigx_l
      !$acc end host_data
      !$acc parallel loop gang vector collapse(2) private(e1w, e2w, rjw, fjj, pp, cphi, igg, lm, ia, m)
      do ig1 = 1, ngc
        do ig2 = 1, ngc
          if(ig2 > ig1) cycle
          e1w = absqg2(ig1); e2w = absqg2(ig2)
          if (dabs(e1w) <= 1d-6 .and. dabs(e2w) <= 1d-6) then   ! the formulas of wronkj for fjj
            rjw = 1d0/rw
            do m = 0, l
              rjw = rjw*rw/(2*m+1)
            enddo
            fjj = -rjw*rjw*(rw*rw*rw)/(2*l+3)
          elseif (dabs(e1w-e2w) > 1d-6) then
            fjj = 1d0/(e2w-e1w)*rw*rw*(ajg(l,ig1)*djg(l,ig2)-djg(l,ig1)*ajg(l,ig2))
          else
            fjj = rw*rw*(ajg(l,ig1)*dje(l,ig1)-djg(l,ig1)*aje(l,ig1))
          endif
          pp = 0d0
          do lm = l*l+1, (l+1)**2
            pp = pp + dconjg(pjyl(lm,ig1))*pjyl(lm,ig2)
          enddo
          cphi = 0d0
          do ia = 1, natg
            cphi = cphi + dconjg(phase(ig1,ia))*phase(ig2,ia)
          enddo
          igg = (ig1*(ig1-1))/2+ig2
          vcoulb(igg) = vcoulb(igg) + cphi*pp*((fpi/(e1w-eee)+fpi/(e2w-eee))*fjj + fpi/(2*l+1)*sx(ig1,ig2))
        enddo
      enddo
    enddo
    !$acc end data
    deallocate(absqg2, cy, yl, pjyl, phase, ajg, djg, aje, dje, sx)
  end subroutine vcoul_termb
  real(8) function fac2m(i)   ! A table of (2l-1)!! data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
    integer:: i,l
    logical,save::  init=.true.
    real(8),save:: fac2mm(0:100)
    if(init) then
      fac2mm(0)=1d0
      do l=1,100
        fac2mm(l)=fac2mm(l-1)*(2*l-1)
      enddo
    endif
    fac2m=fac2mm(i)
  END function fac2m
  !=====================================================================
  subroutine genjh(eee,nr,a,b,lx,nrx,lxx, rofi,rkpr,rkmr) ! Generate radial mesh rofi, spherical bessel, and hankel functions
    ! rkpr, rkmr are real fucntions 
    !i eee=E= -kappa**2 <0
    ! rkpr = (2l+1)!! * j_l(i sqrt(abs(E)) r) * r / (i sqrt(abs(E)))**l
    ! rkmr = (2l-1)!! * h_l(i sqrt(abs(E)) r) * r * i*(i sqrt(abs(E)))**(l+1)
    ! rkpr reduced to be r**l*r      at E \to 0
    ! rkmr reduced to be r**(-l-1)*r at E \to 0
    implicit none
    integer:: nr,lx, nrx,lxx,ir,l
    real(8):: a,b,eee,psi(0:lx),phi(0:lx), rofi(nrx),rkpr(nrx,0:lxx),rkmr(nrx,0:lxx) 
    rofi(1)    = 0d0
    do ir      = 1, nr
      rofi(ir) = b*( exp(a*(ir-1)) - 1d0)
    enddo
    if(eee==0d0) then
      do l = 0,lx
        rkpr(1:nr,l) = rofi(1:nr)**(l +1)
        rkmr(2:nr,l) = rofi(2:nr)**(-l-1 +1)
        rkmr(1,l)    = rkmr(2,l)
      enddo
    else
      do ir  = 1, nr
        call bessl(eee*rofi(ir)**2,lx,phi(0:lx),psi(0:lx))
        do l = 0,lx    !fac2m(l)= (2l-1)!!
          rkpr(ir,l) = phi(l)* rofi(ir)**(l +1) *fac2m(l+1)
          if(ir/=1) rkmr(ir,l) = psi(l)* rofi(ir) **(-l ) /fac2m(l)
        enddo
      enddo
      rkmr(1,0:lx) = rkmr(2,0:lx)
    endif
  end subroutine genjh
  !=============================================================
  subroutine mkjb_4( lxx,lx,nxx,nx,a,b,nr,nrx,rprodx,rofi,rkpr,rkmr, rojb,sgbb) ! make integrals in each MT. and the Fourier matrix.
    implicit none
    integer:: lxx, lx, nxx, nx(0:lxx),nr,nrx, l,n,ir,n1,n2,l1
    real(8):: q(3), rprodx(nrx,nxx,0:lxx),a,b
    real(8):: rojb(nxx, 0:lxx)      !i rho-type onsite integral
    real(8):: sgbb(nxx, nxx, 0:lxx) !i sigma-type onsite integral
    real(8):: fac, xxx,sig, rofi(nrx),rkpr(nrx,0:lxx),rkmr(nrx,0:lxx)
    real(8),parameter:: pi = 4d0*datan(1d0), fpi = 4d0*pi
    rojb=0d0
    sgbb=0d0
    ! rojb
    fac = 1d0
    rojbloop: do l = 0,lx
      fac = fac/(2*l+1)
      do n = 1,nx(l)
        call gintxx(rkpr(1,l), rprodx(1,n,l), a,b,nr, rojb(n,l) )
      enddo
      rojb(1:nx(l),l) = fac*rojb(1:nx(l),l)
    enddo rojbloop
    sgbbloop: do l  = 0,lx
      do n1 = 1,nx(l)
        do n2 = 1,nx(l)
          call sigint_4(rkpr(1,l),rkmr(1,l),lx,a,b,nr,rprodx(1,n1,l),rprodx(1,n2,l), rofi,sig )
          sgbb(n1, n2, l)=sig/(2*l+1)*fpi
        enddo
      enddo
    enddo sgbbloop
  end subroutine mkjb_4
  subroutine sigint_4(rkp,rkm,kmx,a,b,nr,phi1,phi2,rofi, sig)
    implicit none
    integer:: nr,kmx,k,ir
    real(8):: a,b, a1(nr),a2(nr),b1(nr),rkp(nr),rkm(nr), int1x(nr),int2x(nr), phi1(nr), phi2(nr),rofi(nr),sig
    real(8),parameter:: fpi = 4d0*3.14159265358979323846d0
    a1(1) = 0d0;  a1(2:nr) = rkp(2:nr)
    a2(1) = 0d0;  a2(2:nr) = rkm(2:nr)
    b1(1:nr) = phi1(1:nr)
    call intn_smpxxx(a1,b1,int1x,a,b,rofi,nr)
    call intn_smpxxx(a2,b1,int2x,a,b,rofi,nr)
    a1(1) = 0d0; a1(2:nr) = rkm(2:nr) *( int1x(1)-int1x(2:nr) )+ rkp(2:nr) * int2x(2:nr)
    b1(1:nr) = phi2(1:nr)
    call gintxx(a1,b1,A,B,NR, sig )
  end subroutine sigint_4
  subroutine intn_smpxxx(g1,g2,intg,a,b,rofi,nr) ! Intergral of two wave function. used in ppdf
    ! int(r) = \int_(r)^(rmax) u1(r') u2(r') dr' Simpson rule ,and with higher rule for odd devision.
    IMPLICIT none
    integer :: nr,ir,lr0
    real(8) :: g1(nr),g2(nr),intg(nr),a,b,rofi(nr),w1,w2,w3,ooth,foth
    ! if(mod(nr,2) == 0) call rx( ' INTN: nr should be odd for simpson integration rule')
    intg(1)=0d0
    do ir = 3,nr,2
      intg(ir)=intg(ir-2) &
           + 1d0/3d0*G1(IR-2)*G2(IR-2)*( a*(b+rofi(ir-2)) ) &
           + 4d0/3d0*G1(IR-1)*G2(IR-1)*( a*(b+rofi(ir-1)) ) &
           + 1d0/3d0*G1(IR)  *G2(IR)  *( a*(b+rofi(ir)) )
    enddo
    do ir = 2,nr-1,2 ! We use the three-point interpolation used in the Simpson rule. !Checked by bing 2024-11-8
      intg(ir)=intg(ir-1) &
           + 5d0/12d0 *G1(IR-1)*G2(IR-1)*( a*(b+rofi(ir-1)) ) &
           + 2d0/3d0  *G1(IR)  *G2(IR)*  ( a*(b+rofi(ir)  ) ) &
           - 1d0/12d0 *G1(IR+1)*G2(IR+1)*( a*(b+rofi(ir+1)) )
    enddo
    do ir=1,nr
      intg(ir)=intg(nr)-intg(ir)
    enddo
  end subroutine intn_smpxxx
  subroutine sigkernel(nr,a,b,rofi,rkp,rkm,aj,fac,a1) ! a1(r)= fac(r)*( rkm(r) \int_0^r rkp*aj dr' + rkp(r) \int_r^rmax rkm*aj dr' )
    !$acc routine seq
    ! The integrals of intn_smpxxx (Simpson rule at odd points, three-point rule at even points; nr odd) in two
    ! passes without work arrays: the total of rkm*aj, then the running integrals.  aj(1)=0 (r=0), and a1(1)=0.
    implicit none
    integer :: nr, ir
    real(8) :: a, b, rofi(nr), rkp(nr), rkm(nr), aj(nr), fac(nr), a1(nr)
    real(8) :: s2, p1, p2, q1, q2, w, h1o, h2o, h1m, h2m, h1p, h2p
    real(8), parameter :: c3 = 1d0/3d0, c43 = 4d0/3d0, c512 = 5d0/12d0, c23 = 2d0/3d0, c112 = 1d0/12d0
    s2 = 0d0; h2o = 0d0
    do ir = 3, nr, 2
      h2m = rkm(ir-1)*aj(ir-1)*(a*(b+rofi(ir-1)))
      h2p = rkm(ir)*aj(ir)*(a*(b+rofi(ir)))
      s2 = s2 + c3*h2o + c43*h2m + c3*h2p
      h2o = h2p
    enddo
    a1(1) = 0d0
    p1 = 0d0; p2 = 0d0; h1o = 0d0; h2o = 0d0
    do ir = 3, nr, 2
      w = a*(b+rofi(ir-1)); h1m = rkp(ir-1)*aj(ir-1)*w; h2m = rkm(ir-1)*aj(ir-1)*w
      w = a*(b+rofi(ir));   h1p = rkp(ir)*aj(ir)*w;     h2p = rkm(ir)*aj(ir)*w
      q1 = p1 + c512*h1o + c23*h1m - c112*h1p
      q2 = p2 + c512*h2o + c23*h2m - c112*h2p
      a1(ir-1) = (rkm(ir-1)*q1 + rkp(ir-1)*(s2-q2))*fac(ir-1)
      p1 = p1 + c3*h1o + c43*h1m + c3*h1p
      p2 = p2 + c3*h2o + c43*h2m + c3*h2p
      a1(ir) = (rkm(ir)*p1 + rkp(ir)*(s2-p2))*fac(ir)
      h1o = h1p; h2o = h2p
    enddo
  end subroutine sigkernel
  subroutine sigintAn1( absqg, lx, rofi, nr, a1int) ! a1int(r')= r' * \int_0^a r^2 {r_{<}}^l / (r_{>})^{l+1} * j_l(absqg r)/absqg**l
    implicit none
    integer:: nr,l,ir,lx
    real(8):: a1int(nr,0:lx), rofi(nr),absqg
    real(8):: ak(0:lx) ,aj(0:lx), dk(0:lx), dj(0:lx), aknr(0:lx),ajnr(0:lx),dknr(0:lx),djnr(0:lx), phi(0:lx),psi(0:lx)
    if(absqg<1d-10) call rx( "sigintAn1: absqg=0 is not supported yet. Improve here.")
    call radkj(absqg**2, rofi(nr),lx,aknr,ajnr,dknr,djnr,0)
    a1int(1,:) = 0d0
    do ir = 2,nr
      call radkj(absqg**2, rofi(ir),lx,ak,aj,dk,dj,0)
      do l = 0,lx
        a1int(ir,l) = ((2*l+1)* aj(l) -((l+1)* ajnr(l)+ rofi(nr)*djnr(l))* (rofi(ir)/rofi(nr))**l)/absqg**2 *rofi(ir)
      enddo
    enddo
  end subroutine sigintAn1
  subroutine sigintpp( absqg1, absqg2, lx, rmax, sig)! sig(l) =\int_0^a r^2 {r_{<}}^l / (r_{>})^{l+1} *j_l(absqg1 r)/absqg1**l *j_l(absqg2 r)/absqg2**l
    ! e1\ne0 e2\ne0
    implicit none
    integer:: l,lx
    real(8)::  rmax,sig(0:lx), absqg1,absqg2, e1,e2, ak1(0:lx) ,aj1(0:lx), dk1(0:lx), dj1(0:lx), &
         ak2(0:lx) ,aj2(0:lx), dk2(0:lx), dj2(0:lx), fkk(0:lx),fkj(0:lx),fjk(0:lx),fjj(0:lx)
    e1 = absqg1**2
    e2 = absqg2**2
    call wronkj( e1,e2, rmax,lx,   fkk,fkj,fjk,fjj )
    call  radkj( e1,    rmax,lx,   ak1,aj1,dk1,dj1,0)
    call  radkj( e2,    rmax,lx,   ak2,aj2,dk2,dj2,0)
    do l = 0,lx
      sig(l)= (-l*(l+1)*rmax*aj1(l)*aj2(l) +rmax**3*dj1(l)*dj2(l) +0.5d0*rmax**2*(aj1(l)*dj2(l)+aj2(l)*dj1(l)) -fjj(l)*(2*l+1)*(e1+e2)/2d0)&
           /(e1*e2)
    enddo
  end subroutine sigintpp
endmodule m_vcoulq
