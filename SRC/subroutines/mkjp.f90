module  m_vcoulq
  use m_cputm,only:cputm
  use m_mpi,only: ipr,mpi__rank
  use m_lgunit,only: stdo
  use m_ftox
  public vcoulq_4,mkjb_4,mkjp_4,genjh, ajr
  private
  character(1024):: aaaw
  real(8), allocatable  :: ajr(:,:,:)
contains
  subroutine vcoulq_4(q,nbloch,ngc,nbas,lx,lxx,nx,nxx,alat,qlat,vol,ngvecc, & !Coulmb matrix for each q
       strx,rojp,rojb,sgbb,sgpb,fouvb,ngb,bas,rmax, eee, aa,bb,nr,nrx,rkpr,rkmr,rofi,    vcoul)
    use m_ll,only: ll
    use m_blas, only: m_op_T,m_op_C, zmv_h !, zmm_h
#ifdef __GPU
    use m_blas, only: dmm => dmm_d, zmm=>zmm_d
#else
    use m_blas, only: dmm => dmm_h, zmm=>zmm_h
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
    integer :: ngc, ngvecc(3,ngc),nblochngc,nev,nmx,ix,nrx,nr(nbas),ir,ig,lm
    real(8) :: egtpi,vol,q(3), qlat(3,3),alat,absqg2(ngc),qg(3), rojb(nxx, 0:lxx, nbas)
    real(8) :: sgbb(nxx,  nxx,  0:lxx,      nbas) !i sigma-type onsite integral
    real(8) :: fpivol,tpiba, bas(3,nbas),r2s,rmax(nbas)
    real(8) ::  fkk(0:lxx),fkj(0:lxx),fjk(0:lxx),fjj(0:lxx),sigx(0:lxx),radsig(0:lxx) !,radsig(0:lxx,nbas),fjj(0:lxx,nbas)
    real(8) :: eee , int1x(nrx),int2x(nrx),phi(0:lxx),psi(0:lxx) &
         ,aa(nbas),bb(nbas),rkpr(nrx,0:lxx,nbas),rkmr(nrx,0:lxx,nbas),rofi(nrx,nbas)
    real(8),allocatable  :: eb(:),cy(:),yl(:), a1(:,:,:)
    real(8),parameter:: pi=4d0*datan(1d0),fpi=4d0*pi
    complex(8) :: rojp(ngc, (lxx+1)**2, nbas)   !rho-type onsite integral
    complex(8) :: strx((lxx+1)**2, nbas, (lxx+1)**2,nbas) !structure constant. The multicenter expantion of 1/|r-r'|
    complex(8) :: sgpb(ngc,  nxx,  (lxx+1)**2, nbas)
    complex(8) :: fouvb(ngc,  nxx, (lxx+1)**2, nbas) ,vcoul(ngb, ngb) !<exp(i q+G r)|xxx>
    !complex(8),allocatable :: hh(:,:),oo(:,:),zz(:,:),matp(:),matp2(:),pjyl_(:,:),phase(:,:)
    complex(8),allocatable :: hh(:,:),oo(:,:),zz(:,:),matp(:),matp2(:),pjyl_(:,:),phase(:,:),pjyl_p(:,:) !,pjyl_p(:,:,:)
    complex(8) :: xxx, img=(0d0,1d0), fouvp_ig1_ig2, fouvp_ig2_ig1, sgpp_ig1_ig2
    integer :: istat,lm2x
    integer,allocatable :: llx(:)
    write(aaaw,'(" vcoulq_4: ngb  nbloch ngc nrx procid=",5i6)') ngb,nbloch,ngc,nrx,mpi__rank
    call cputm(stdo,aaaw)
    fpivol = 4*pi*vol
    allocate( pjyl_((lxx+1)**2,ngc),pjyl_p((lxx+1)**2,ngc),phase(ngc,nbas),source=(0d0,0d0) )!,pjyl_p((lxx+1)**2,ngc,nbas)
    allocate( cy((lxx+1)**2), yl((lxx+1)**2),source=0d0)
    allocate( llx((lxx+1)**2),source=0)
    !allocate( pjyl_((lxx+1)**2,ngc),phase(ngc,nbas), cy((lxx+1)**2), yl((lxx+1)**2))
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
    !$acc enter data copyin(strx, rojb, rojp) create(vcoul) copyin(ibasbl, nbl, lbl, mbl, lmbl)
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
      !$acc data copyin(fouvb, sgpb) 
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
      use m_bessl, only: bessl2 => bessl, wronkj2 => wronkj
      use m_keyvalue,only: getkeyvalue
      ! real(8), allocatable :: sigx_tmp(ngc,ngc,0:lxx), a1g(nrx,ngc), aabb_by3
      ! real(8) :: ajr_tmp(nrx,ngc), phi_rg(nrx,ngc,0:lxx), rofi_tmp(1:nrx) !  complex(8) :: crojp((lxx+1)**2,nbas,ngc)
      real(8), allocatable :: fac_integral(:), a1g(:,:), ajr_tmp(:,:), phi_rg(:,:,:), rofi_tmp(:)
      complex(8) :: rojpstrx((lxx+1)**2,nbas,ngc)
      logical :: hasBessel, keepWronkj
      real(8), allocatable ::  keep_fjj(:,:), keep_sigx(:,:), sigx_tmp(:,:)
      integer, allocatable :: iggtable(:,:)
      integer :: nggc, igg
      ! Get integral coefficients of int (a*b) G_1(ir) G_2(ir) exp(a*r))
      ! simpson rule is used. nr(ibas) was set as odd number
      !   sigx_tmp(ig1,ig2,l) is int dr (aa(ibas)*bb(ibas)) a1g(r,g1)* ajr(r,l,ibas,g2) exp(aa(ibas)*r))
      call getkeyvalue("GWinput","KeepWronkj",keepWronkj,default=.true.)
      write(aaaw,ftox) " vcoulq_4: goto PvP procid ngc lxx nrx=", mpi__rank,ngc,lxx,nrx
      call cputm(stdo,aaaw)

      !make table for ig1,ig2 from one dimensional index (igg = 1, ..., ngg)
      nggc = (ngc*(ngc+1))/2
      allocate(iggtable(2,nggc))
      igg = 0
      do ig1 = 1,ngc
        do ig2 = 1, ig1
          igg = igg + 1
          iggtable(1,igg) = ig1
          iggtable(2,igg) = ig2
        enddo
      enddo

      !$acc data create(rojpstrx,pjyl_p) copyin(rofi, rkpr, rkmr, aa, bb, nr, absqg2, pjyl_, phase, iggtable)

      !$acc host_data use_device(strx, rojp)
      istat = zmm(strx, rojp, rojpstrx, m=nbas*(lxx+1)**2, n=ngc, k=nbas*(lxx+1)**2, opA=m_op_T, opB=m_op_C)
      !$acc end host_data

      write(aaaw,ftox) " vcoulq_4: goto igig loop", mpi__rank
      call cputm(stdo,aaaw)
      lm2x= (lxx+1)**2

      igigLoopSlow: do ibas= 1, nbas
        !$acc kernels loop collapse(2)
        do ig1 = 1,ngc
          do lm2=1,(lx(ibas)+1)**2
            pjyl_p(lm2,ig1)=pjyl_(lm2,ig1)*phase(ig1,ibas)
          enddo
        enddo
        !$acc end kernels

        if(eee==0d0) then
          !copy GPU -> CPU 
          !$acc update self(pjyl_p, rojpstrx(1:(lxx+1)**2,ibas,1:ngc), vcoul)
          do ig1 = 1,ngc !this loop is slow for large system, but maybe vcoulq_4 is rather the critical step 
            do ig2 = 1,ig1
              call wronkj( absqg2(ig1), absqg2(ig2), rmax(ibas),lx(ibas), fkk,fkj,fjk,fjj)
              call sigintpp( absqg2(ig1)**.5, absqg2(ig2)**.5, lx(ibas), rmax(ibas), sigx)
              radsig(0:lxx) = 0d0 
              forall(l = 0:lx(ibas)) radsig(l) = fpi/(2*l+1) * sigx(l)
              vcoul(nbloch+ig1,nbloch+ig2) =  vcoul(nbloch+ig1,nbloch+ig2) + sum( rojpstrx(1:lm2x,ibas,ig1)*rojp(ig2, 1:lm2x, ibas) &
                   + dconjg(pjyl_p(1:lm2x,ig1))*pjyl_p(1:lm2x,ig2)* &
                   ( (fpi/(absqg2(ig1)-eee)+fpi/(absqg2(ig2)-eee)) *fjj(llx(1:lm2x)) + radsig(llx(1:lm2x)) )   )
            enddo
          enddo
          !$acc update device(vcoul)
        else !eee is nonzero

          hasBessel = .false.
          if(ibas > 1) then
            if(nr(ibas) == nr(ibas-1)) then
              if( all(abs(rofi(1:nr(ibas),ibas) - rofi(1:nr(ibas-1),ibas-1)) < 1d-10) .and. lx(ibas) == lx(ibas-1) ) hasBessel = .true.
            endif
          endif

          setBessel: if(.not.hasBessel) then
            allocate(phi_rg(nr(ibas), ngc, 0:lx(ibas)))
            allocate(rofi_tmp(1:nr(ibas)), fac_integral(1:nr(ibas)), a1g(nr(ibas),ngc), ajr_tmp(nr(ibas),ngc))
            allocate(sigx_tmp(ngc,ngc))
            !$acc data create(phi_rg, ajr_tmp, a1g, rofi_tmp, fac_integral, sigx_tmp)

            !$acc parallel loop collapse(2) private(phi(0:lxx), psi(0:lxx))
            do ig = 1, ngc
              do ir = 1, nr(ibas)
                call bessl2(absqg2(ig)*rofi(ir,ibas)**2,lx(ibas),phi, psi)
                phi_rg(ir,ig,0:lx(ibas)) = phi(0:lx(ibas))
              enddo
            enddo
            !$acc end parallel

            if(keepWronkj) then
              if(allocated(keep_fjj)) then
                !$acc exit data delete(keep_fjj)
                 deallocate(keep_fjj)
              endif
              allocate(keep_fjj(0:lx(ibas),nggc))
              !$acc enter data create(keep_fjj)
              !$acc parallel loop private(fkk(0:lxx), fkj(0:lxx), fjk(0:lxx), fjj(0:lxx))
              do igg = 1, nggc
                ig1 = iggtable(1,igg)
                ig2 = iggtable(2,igg)
                call wronkj2( absqg2(ig1), absqg2(ig2), rmax(ibas),lx(ibas), fkk,fkj,fjk,fjj)
                keep_fjj(0:lx(ibas),igg) = fjj(0:lx(ibas))
              enddo
              !$acc end parallel
            endif

            if(allocated(keep_sigx)) then
              !$acc exit data delete(keep_sigx)
              deallocate(keep_sigx)
            endif
            allocate(keep_sigx(0:lx(ibas),nggc))
            !$acc enter data create(keep_sigx)

            !$acc kernels
            do ir = 1, nr(ibas)
              fac_integral(ir) = aa(ibas)*bb(ibas)*dexp(aa(ibas)*(ir-1))/3d0
              if( ir /= 1 .and. ir /= nr(ibas)) fac_integral(ir) = fac_integral(ir)*merge(4d0,2d0,mod(ir,2)==0)
            enddo
            !$acc end kernels
            do l = 0, lx(ibas)
              !$acc kernels
              rofi_tmp(1:nr(ibas)) = rofi(1:nr(ibas),ibas)**(l+1)
              !$acc end kernels
              !$acc kernels loop independent private(int1x, int2x)
              do ig = 1, ngc
                ajr_tmp(1:nr(ibas),ig) = phi_rg(1:nr(ibas),ig,l)*rofi_tmp(1:nr(ibas))
                call intn_smpxxx( rkpr(1,l,ibas), ajr_tmp(1,ig),int1x,aa(ibas),bb(ibas),rofi(1,ibas),nr(ibas))
                call intn_smpxxx( rkmr(1,l,ibas), ajr_tmp(1,ig),int2x,aa(ibas),bb(ibas),rofi(1,ibas),nr(ibas))
                ! a1g(1:nr(ibas),ig) = [0d0,(rkmr(2:nr(ibas),l,ibas) *( int1x(1)-int1x(2:nr(ibas)) ) &
                !                          + rkpr(2:nr(ibas),l,ibas) *  int2x(2:nr(ibas)))* fac_integral(2:nr(ibas))]  ! error in GPU version 08/18/2025
                a1g(1,ig) = 0d0
                do ir = 2, nr(ibas)
                  a1g(ir,ig) = (rkmr(ir,l,ibas) * (int1x(1) - int1x(ir)) + rkpr(ir,l,ibas) * int2x(ir)) * fac_integral(ir)
                enddo
              enddo
              !$acc end kernels
              istat = dmm(a1g, ajr_tmp, sigx_tmp, m=ngc, n=ngc, k=nr(ibas), opA=m_op_T)
              !$acc kernels
              do igg = 1, nggc
                ig1 = iggtable(1,igg)
                ig2 = iggtable(2,igg)
                keep_sigx(l,igg) = sigx_tmp(ig1,ig2)
              enddo
              !$acc end kernels
            enddo

            !$acc end data
            deallocate(ajr_tmp, a1g, rofi_tmp, fac_integral, phi_rg, sigx_tmp)
          endif setBessel

          write(aaaw,ftox) " vcoulq_4:  igig loop procid ibas nr lx, hasBessel=", mpi__rank,ibas, nr(ibas), lx(ibas), hasBessel
          call cputm(stdo,aaaw)

          !$acc parallel loop private(fkk(0:lxx), fkj(0:lxx), fjk(0:lxx), fjj(0:lxx), sigx(0:lxx), radsig(0:lxx)) present(keep_sigx)
          do igg = 1, nggc
            ig1 = iggtable(1,igg)
            ig2 = iggtable(2,igg)
            if(keepWronkj) then
              fjj(0:lx(ibas)) = keep_fjj(0:lx(ibas),igg)
            else
              call wronkj2( absqg2(ig1), absqg2(ig2), rmax(ibas),lx(ibas), fkk,fkj,fjk,fjj)
            endif
            sigx(0:lx(ibas)) = keep_sigx(0:lx(ibas),igg)
            radsig(0:lxx) = 0d0 
            forall(l = 0:lx(ibas)) radsig(l) = fpi/(2*l+1) * sigx(l)
            vcoul(nbloch+ig1,nbloch+ig2) =  vcoul(nbloch+ig1,nbloch+ig2) + sum( rojpstrx(1:lm2x,ibas,ig1)*rojp(ig2, 1:lm2x, ibas) &
                 + dconjg(pjyl_p(1:lm2x,ig1))*pjyl_p(1:lm2x,ig2)* &
                 ( (fpi/(absqg2(ig1)-eee)+fpi/(absqg2(ig2)-eee)) *fjj(llx(1:lm2x)) + radsig(llx(1:lm2x)) )   )
          enddo
          !$acc end parallel

        endif
      enddo igigLoopSlow
      if(allocated(keep_fjj)) then
        !$acc exit data delete(keep_fjj)
        deallocate(keep_fjj)
      endif
      if(allocated(keep_sigx)) then
        !$acc exit data delete(keep_sigx)
        deallocate(keep_sigx)
      endif
      deallocate(iggtable)

      !$acc kernels
      do ig1 = 1, ngc
        vcoul(nbloch + ig1, nbloch + ig1) = vcoul(nbloch + ig1, nbloch + ig1) + fpivol/(absqg2(ig1) - eee) !eee is negative
      end do
      !$acc end kernels

      !$acc end data
    endblock PvP_dev_mo

    !$acc exit data copyout(vcoul) delete(strx, rojp)

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
    use m_ll,only: ll
    use m_bessl, only: bessl2 => bessl
    implicit none
    integer:: ngc,ngvecc(3,ngc), lxx, lx, nxx,nx(0:lxx),nr,nrx, nlx,ig1,ig2,l,n,ir,n1,n2,lm !, ibas
    real(8):: q(3),bas(3), rprodx(nrx,nxx,0:lxx),a,b,rmax,alat, qlat(3,3)
    real(8):: pi,fpi,tpiba, qg1(3), fkk(0:lx),fkj(0:lx),fjk(0:lx),fjj(0:lx),absqg1,absqg2, &
         fac,radint,radsigo(0:lx),radsig(0:lx),phi(0:lx),psi(0:lx),r2s,sig,sig1,sig2,sigx(0:lx),sig0(0:lx) ,qg2(3)
    real(8):: rofi(nrx),rkpr(nrx,0:lxx),rkmr(nrx,0:lxx),eee,qg1a(3)
    real(8),allocatable::cy(:),yl(:)
    real(8),allocatable ::a1(:,:,:), qg(:,:),absqg(:), rofi_nr(:)
    complex(8) :: rojp(ngc, (lxx+1)**2)        ! rho-type onsite integral
    complex(8) :: sgpb(ngc,  nxx,  (lxx+1)**2) !sigma-type onsite integral
    complex(8) :: fouvb(ngc,  nxx, (lxx+1)**2)
    complex(8) :: img =(0d0,1d0),phase
    complex(8),allocatable :: pjyl(:,:)
    logical, intent(in) :: hasBessel
    nlx = (lx+1)**2
    ! allocate(ajr(1:nr,0:lx,ngc),a1(1:nr,0:lx,ngc), qg(3,ngc),absqg(ngc), pjyl((lx+1)**2,ngc) )
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
    ! rojp
    write(aaaw,ftox)' mkjp_4: goto rojploop'
    call cputm(stdo,aaaw)
    rojp = (0d0, 0d0)
    rojploop: do ig1 = 1,ngc
      call wronkj( absqg(ig1)**2, eee, rmax,lx, fkk,fkj,fjk,fjj)
      do lm = 1,nlx
        l = ll(lm)
        rojp(ig1,lm) = (-fjj(l))* pjyl(lm,ig1)
      enddo
    enddo rojploop

    allocate(rofi_nr, source = rofi(1:nr))
    !$acc enter data copyin(absqg(1:ngc), rofi_nr(1:nr))

    setBessel: if(.not.hasBessel) then
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
    endif setBessel
    !-------------------------
    if(eee==0d0) then
      allocate(a1(1:nr,0:lx,ngc))
      do ig1 = 1,ngc
        call sigintAn1( absqg(ig1), lx, rofi_nr, nr,a1(1:nr, 0:lx,ig1) )
      enddo
      !      else       ! We need to implement a version of sigintAn1 to treat eee/=0 case...
    endif
    write(aaaw,ftox)' mkjp_4: goto dev_mo block. size of ajr:,nx', size(ajr), nx(:)
    call cputm(stdo,aaaw)
    dev_mo: block
#ifdef __GPU
      use m_blas, only: dmm => dmm_d, m_op_T
#else
      use m_blas, only: dmm => dmm_h, m_op_T
#endif
      real(8):: a1work(nr), a2work(nr), int1x(nr), int2x(nr), a1g(nr,ngc), fac_integral(1:nr)
      real(8), allocatable :: sigg(:,:,:), radintg(:,:,:)
      integer :: llist(nlx)
      integer:: istat

      llist(1:nlx) = [(ll(lm), lm=1, nlx)]
      fac_integral(1) = a*b/3d0
      do ir = 2, nr
        fac_integral(ir) = fac_integral(ir-1)*dexp(a)
      enddo
      forall(ir=2:nr-1) fac_integral(ir) = fac_integral(ir)*merge(4d0,2d0,mod(ir,2)==0)

      !$acc data copyin(fac_integral, rkpr, rkmr, pjyl, rprodx, llist, nx) create(a1g)
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
      else
        allocate(sigg(ngc,nxx,0:lx))
        !$acc data create(sigg) copyout(sgpb)
        do l = 0, lx
          if(nx(l) == 0) cycle
          !$acc kernels loop independent private(int1x, int2x, a1work, a2work)
          do ig1 = 1, ngc
            a1work(1) = 0d0;  a1work(2:nr) = rkpr(2:nr,l)
            a2work(1) = 0d0;  a2work(2:nr) = rkmr(2:nr,l)
            call intn_smpxxx(a1work,ajr(1,ig1,l),int1x,a,b,rofi_nr,nr)
            call intn_smpxxx(a2work,ajr(1,ig1,l),int2x,a,b,rofi_nr,nr)
            a1g(1,ig1) = 0d0
            a1g(2:nr,ig1) = rkmr(2:nr,l) *( int1x(1)-int1x(2:nr) )+ rkpr(2:nr,l) * int2x(2:nr)
            a1g(1:nr,ig1) = a1g(1:nr,ig1)*fac_integral(1:nr)
          enddo
          !$acc end kernels
          istat = dmm(a1g, rprodx(1,1,l), sigg(1,1,l), m=ngc, n=nx(l), k=nr, opA=m_op_T, ldB=nrx)
        enddo
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
!      if(ipr) write(stdo,ftox)' mkjp_4: fouvb dev block'
      allocate(radintg(ngc,nxx,0:lx))
      !$acc data create(radintg) copyout(fouvb)
      do l = 0, lx
        if(nx(l) == 0) cycle
        !$acc kernels loop independent present(ajr)
        do ig1 = 1, ngc
          a1g(1:nr,ig1) = ajr(1:nr,ig1,l)*fac_integral(1:nr)
        enddo
        !$acc end kernels
        istat = dmm(a1g, rprodx(1,1,l), radintg(1,1,l), m=ngc, n=nx(l), k=nr, opA=m_op_T, ldB=nrx)
      enddo
      !$acc kernels
      fouvb(:,:,:) = 0d0
      !$acc end kernels
      !$acc kernels loop independent collapse(2) present(absqg)
      do lm = 1, nlx
        do n = 1, nxx
          l = llist(lm)
          if(n > nx(l)) cycle
          fouvb(1:ngc, n, lm) = fpi/(absqg(1:ngc)**2-eee) *dconjg(pjyl(lm,1:ngc))*radintg(1:ngc,n,l)
        enddo
      enddo
      !$acc end kernels
      !$acc end data
      deallocate(radintg)

      !$acc end data
    endblock dev_mo

    !$acc exit data delete(absqg, rofi_nr)
    deallocate(absqg, qg, pjyl, cy, yl, rofi_nr)
  end subroutine mkjp_4
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
    !$acc routine seq
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
