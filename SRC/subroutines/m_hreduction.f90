module m_hreduction
use m_cmdopt_registry, only: c0_gs, c0_mlo_diagnorm, c0_mlo_feb4, c0_mlo_ortho, c0_mlo_orthonorm
contains
  subroutine Hreduction(mlomethod,iprx,ndimPMT,hamm,ovlm,ndimMTO,ix,fff1, hammout,ovlmout, qp, cmlo,nev, zMLO) !> Reduce H(ndimPMT) to H(ndimMTO)
    ! cmlo= <Psi^MPT i|F^MLO j>
   use m_zhev,only:zhev_tk4
   use m_nvfortran, only: findloc
   use m_readqplist,only: eferm
!   use m_HamPMT,only: GramSchmidt!,epsovl
   use m_lgunit,only:stdo
   use m_lmfinit,only:oveps
   use m_keyvalue,only: getkeyvalue
   use m_GWinput, only: gwinput_init, gwinput_loaded, &
                        tg_mlo_nskip => mlo_nskip, tg_mlo_eww => mlo_eww, &
                        tg_mlo_emax => mlo_emax
   implicit none
   integer::i,j,ndimPMT,ndimMTO,nx,nmx,ix(ndimMTO),nev,nxx,jj,ndimPMTx,nvpmt,mlomethod,nskip,nskipin
   real(8)::beta,emu,val,wgt(ndimPMT),evlmto(ndimMTO),evl(ndimPMT),evlx(ndimPMT),qp(3),eww,eadd
   complex(8):: evecmto(ndimMTO,ndimMTO),evecpmt(ndimPMT,ndimPMT)
   complex(8):: ovlmx(ndimPMT,ndimPMT),hammx(ndimPMT,ndimPMT),fac(ndimPMT,ndimMTO),ddd(ndimMTO,ndimMTO)
   complex(8):: hamm(ndimPMT,ndimPMT),ovlm(ndimPMT,ndimPMT)
   complex(8):: hammout(ndimMTO,ndimMTO),ovlmout(ndimMTO,ndimMTO)
   complex(8),optional,intent(out):: cmlo(ndimPMT,ndimMTO) !<Psi^PMT_i|F^MLO_k>, PMT-eigenstate-basis coefficients
   complex(8),optional,intent(out):: zMLO(ndimPMT,ndimMTO) !|F_MLO_k> = sum_m |chi^PMT_m> zMLO(m,k)
   complex(8):: cmlo_loc(ndimPMT,ndimMTO) !internal working array
!   complex(8),optional:: zcplz(ndimPMT,ndimMTO)
   complex(8),allocatable :: Amat(:,:)
   real(8):: fff1,fff !epsovl=1d-8 epsovlm=0d0 ,
   logical:: iprx
   ! choosed MTO Hamiltonian. Get evlmto,evecmto
   ovlmx= ovlm
   hammx= hamm
   nmx = ndimMTO !  write(stdo,*)'Start Hreduction: 111'
   call zhev_tk4(ndimMTO,hamm(ix(1:ndimMTO),ix(1:ndimMTO)),ovlm(ix(1:ndimMTO),ix(1:ndimMTO)), nmx,nev, evlmto, evecmto, oveps)
   if(nev/=ndimMTO) call rx('Hreduction: nev/=ndimMTO We didnot get eigenfuncitons of ndimMTO. Linear dependency problem?')
   ! PMT Hamiltonian
   ovlm= ovlmx
   hamm= hammx
   nmx = ndimPMT !  write(stdo,*)'Start Hreduction: 222'
   call zhev_tk4(ndimPMT,hamm(1:ndimPMT,1:ndimPMT),ovlm(1:ndimPMT,1:ndimPMT), nmx,nev, evl,evecpmt, oveps) !PMT
   ovlm=ovlmx
   ndimPMTx=nev !obtained. oveps may reduce ndimPMT to be ndimPMTx
   fac = (0d0,0d0)
   do j=1,ndimMTO !Amat is corrected matrix element of fac=<psi_PMT|psi_MTO>
      do i=1,nev
         fac(i,j)= sum(dconjg(evecpmt(:,i))*matmul(ovlmx(:,ix(1:ndimMTO)),evecmto(1:ndimMTO,j))) !<Psi_PMT|Psi_MTO>
      enddo
   enddo
   ModifyMatrixElements :block
      use m_nvfortran,only : findloc
      use m_ftox
      integer:: ie,nidxevlmto,nidxevl,ibx,jx,idxevlmto(ndimMTO),idxevl(ndimPMT),jbx,nval,nnn,imx,nbx,ii
      real(8):: eee,fffx,ecut,xxx,rydberg,facww,sss,fff,epscore,emax,alpha,emin,ww(ndimPMTx),dex,ddd !,ewcutf
      real(8),allocatable::mulfac(:,:),mulfacw(:,:)
      complex(8):: imag=(0d0,1d0)
      ! Assert block for normalization check
      do j=1,ndimMTO 
        if(abs(sum(abs(fac(:,j))**2)-1d0)>1d-4) call rxi('Hreduction: normalization error band index=',j)
      enddo
      if(iprx) then
        do j=1,ndimMTO !Amat is corrected matrix element of <psi_PMT|psi_MTO>
          do i=1,ndimPMTx
            if(abs(fac(i,j))**2>.1) write(stdo,ftox)'fac matrix ',j,i,ftof(abs(fac(i,j))**2)
          enddo
        enddo
      endif
      
      ! Determine nskip, eigenfunctions PMT(1:nskip), semicores, are removed.
      epscore=0.5d0
      nskipin = findloc( sum(abs(fac(:,:))**2,dim=2) > epscore, value=.true.,dim=1)-1 !semicore level skip by LO. Or skip evec outside of MTOa
      call gwinput_init()
      if (gwinput_loaded) then
         nskip = tg_mlo_nskip
         if (nskip == -huge(0)) nskip = nskipin   ! sentinel ⇒ key absent
      else
         call rx('m_GWinput: legacy GWinput reader is disabled. GWinput.toml is required.')
!         call getkeyvalue("GWinput","mlo_nskip",nskip,default=nskipin) !nskip is LO bands. This will be automatic
      endif
      write(stdo,ftox) 'nnnnn nskip',nskip !,ftof(sum(abs(fac(:,:))**2,dim=2))

      ! === Usage ===
      ! Simple version ---> Set mlo_method 0 with mlo_emax for Semiconductor (\lesssim VBM) or Al2O3_Cr (7eV or higher)
      !                     Not necessary for NiO, Ru2O3.  
      ! Generally speaking, we need emax for localized bands, while we need mlomethod2 for semiconductors (broad sp bands, smooth cutoff).
      ! 1. Only localized bands, I think no switch needed.
      ! 2. For semiconductors, set emax = Efermi (or even -9999) around ( ---> then mlomethod0 is close to mlomethod2).
      ! 3. For Al2O3_Cr (sp and d bandd), we need to set emax a little about the localized bands.
      ! --------------------
      !
      ! For sp bands smooth cutoff.
      !  Semiconductors: Si, GaAs,
      !    + We have to set mlo_emax 0 or something. For Al2O3_Cr, we need to set mlo_emax as 15 eV or so.
      !    (For Al2O3_Cr, we found mlomethod 1 with emax= 7 eV works well).
      !    
      !  NiO, Ru2O3 
      !    mlo_method 1 auto emax, or mlo_method 0 auto emax. auto emax 
      ! 
      !  Extract 3d or 4f bands
      !     + mlomethod 0 works.
      !     mlomethod 2 works for 4f extraction. No emax
      !
      !  We have to set mlo_emax, up to which we have to include i for fac=<Psi_PMT(i)|Psi_MTO(j)> for semiconductors or broad band included.
      !

      ! When P = \sum_i \sum_j |Psi^PMT_i> <Psi^PMT_i|Psi^MTO_j> <Psi^MTO_j|, we have |F^MTO_k>=  P| F^MTO_k>. That is P is identical operator.
      ! Instead of <Psi^PMT_i|Psi^MTO_j>, we use Amat which is a modified version.
      if (gwinput_loaded) then
         eww = tg_mlo_eww
      else
         call rx('m_GWinput: legacy GWinput reader is disabled. GWinput.toml is required.')
!         call getkeyvalue("GWinput","mlo_eww",eww,default=0.2d0) !smoothing cutoff
      endif
      emax = evl(ndimMTO+nskip) - eferm   ! emax is the max of evl at ndimMTO+nskip. This is mainly useful for localized bands range.
!      emax = evlmto(ndimMTO) - eferm
      if (gwinput_loaded) then
         eee = tg_mlo_emax
         if (eee == huge(0d0)) eee = emax*rydberg()  ! sentinel ⇒ key absent
      else
         call rx('m_GWinput: legacy GWinput reader is disabled. GWinput.toml is required.')
!         call getkeyvalue("GWinput","mlo_emax",eee,default=emax*rydberg())  !eV relative to Ef.
      endif
      emax=eee/rydberg()+eferm

      !Amat is a modification of fac(ndimPMTx,ndimMTO), which is <Psi_PMT_i |Psi_MTO j>. Psi are eigenfunctions.
      allocate(Amat(ndimPMTx,ndimMTO),source=(0d0,0d0))!this is to avoid bug in ifort18.0.5
      mloloop : do j=1,ndimMTO 
        if(mlomethod==0) then          ! Determine ecut for j to determine maxmum i index for PMT.
          ecut = max(emax, evlmto(j))  !  emax(relative to ef) is the rigid limit for localized MTOs
        elseif(mlomethod==1) then
          ecut = emax 
        elseif(mlomethod==2) then  
          ecut = evlmto(j)  
        endif
        pmtloop: do i=nskip+1,ndimPMTx 
          Amat(i,j)= fac(i,j) * fermidist( (evl(i) - ecut) /eww)
        enddo pmtloop
      enddo mloloop
      Amat(1:nskip,:)=0d0
! do we need GramSchmidt orthogonalizaition? We expect lower is enphasized more for mode0 and for mode2.
      if(c0_gs) call GramSchmidt(ndimPMTx,ndimMTO,Amat) !Amat= ¥bar{<Psi_PMT_i |Psi_MTO j>}
      
      ! cmlo(i,k) = \sum_j ¥bar{<Psi_PMT_i |Psi_MTO j>} <Psi_MTO j|F_MTO k> = <Psi_PMT_i|F_MLO k>
      ! |F^MLO_k> = P |F_MTO_k> = (\sum_{i,j} |Psi_PMT_i> ¥bar{<Psi_PMT_i |Psi_MTO j>} <Psi_MTO j|) |F_MTO k>=  |Psi_PMT> C * zMTO 
      ! Here P is a projector-like (probably not a projector because of GramSchmidt).
      nx = ndimPMTx
      cmlo_loc(ndimPMTx+1:ndimPMT,1:ndimMTO)=0d0
      cmlo_loc(1:ndimPMTx,1:ndimMTO) = matmul(Amat(1:ndimPMTx,1:ndimMTO),&
           matmul(transpose(dconjg(evecmto(:,:))),ovlmx(ix(1:ndimMTO),ix(1:ndimMTO)))) ! where <Psi_MTO j|MTO_k> = (evecmto*) @ ovlmx

      ! Per-orbital diagonal normalization: |F^MLO_i> -> |F^MLO_i>/sqrt(<F^MLO_i|F^MLO_i>).
      ! Diagonal-only; off-diagonal overlap is left untouched. Use --mlo_diagnormalization
      ! to enable. Matches Feb 2026 commit 464a2d510 behavior when on.
      ! (--mlo_ortho below performs full Lowdin orthogonalization.)
      MLODiagonalNormalize: if (c0_mlo_diagnorm .or. c0_mlo_feb4) then
         do i = 1, ndimMTO
            ddd = sum(dconjg(cmlo_loc(1:nx,i))*cmlo_loc(1:nx,i)) !<F^MLO|F^MLO>
            cmlo_loc(1:nx,i) = cmlo_loc(1:nx,i)/sqrt(ddd)
         enddo
      endif MLODiagonalNormalize
      MLOLowdinOrthogonalization:if(c0_mlo_orthonorm .or. c0_mlo_ortho) then
        block
          use m_lapack, only: zhev => zhev_h
          complex(8) :: ovlm_mlo(ndimMTO,ndimMTO), evl_ovl_buf(ndimMTO,ndimMTO), sinv_half(ndimMTO, ndimMTO)
          real(8) :: eval(ndimMTO), einv_half
          real(8), parameter :: eps = 1d-12, eps_ovl_chk = 1d-8
          integer :: istat
          ovlm_mlo = matmul(dconjg(transpose(cmlo_loc)), cmlo_loc)
          istat = zhev(ovlm_mlo, n=ndimMTO, evl=eval)
          do i = 1, ndimMTO
            einv_half = merge(0d0, 1d0/sqrt(eval(i)), eval(i) < eps)
            evl_ovl_buf(:,i) = ovlm_mlo(:,i)*einv_half
          enddo
          sinv_half = matmul(evl_ovl_buf, transpose(dconjg(ovlm_mlo)))
          cmlo_loc = matmul(cmlo_loc, sinv_half)
          !check
          ovlm_mlo = matmul(dconjg(transpose(cmlo_loc)), cmlo_loc)
          forall(i=1:ndimMTO) ovlm_mlo(i,i) = ovlm_mlo(i,i) - 1d0
          if (any(abs(ovlm_mlo) > eps_ovl_chk)) call rx('Hreduction: LowdinOrthogonalization FAILD')
        endblock
      endif MLOLowdinOrthogonalization


      ! |F^MLO j'>= |F^PMT_i'> z^PMT_i'i cmlo(i,j)
      do i=1,ndimMTO
        do j=1,ndimMTO
          hammout(i,j)= sum( dconjg(cmlo_loc(1:nx,i))*evl(1:nx)*cmlo_loc(1:nx,j)) !|F^MLO_i> = |Psi^PMT_j> cmlo(j,i)
          ovlmout(i,j)= sum( dconjg(cmlo_loc(1:nx,i))*cmlo_loc(1:nx,j) ) !<F^MLO|F^MLO>
        enddo
      enddo
      if(present(cmlo)) cmlo = cmlo_loc
      if(present(zMLO)) zMLO = matmul(evecpmt(:,1:nx), cmlo_loc(1:nx,:)) !PMT-basis-function coefficients

!       GetZCZ: if(present(zcplz)) then ! This is for k mesh of GWinput. zcplz(ndimPMT, ndimMLO), ndimMLO=ndimMTO
!         diagonalizeMTO: block         !  |F^MLO j> = |Psi^PMT i> zcplz(i,j)  !MLO eigenfunctions constructed from the PMT basis.
!           real(8):: evlx(ndimMTO),oveps=0d0
!           complex(8):: evecmlo(ndimMTO,ndimMTO),hh(ndimMTO,ndimMTO),oo(ndimMTO,ndimMTO)
!           hh=hammout
!           oo=ovlmout
!           call zhev_tk4(ndimMTO,hh,oo,ndimMTO,nev, evlx,evecmlo, oveps)
!           ! nmx=ndimMTO.
!           ! CAUTION for zhev_tk4: If nmx=0, only eigenvalues returned. Diangonalize (hamm- evl ovlm) z=0
!           ! evecmlo= z^MLO_jj' convert to eigenfunction Psi
! !          zcplz = matmul(matmul(evecpmt(1:ndimPMT,1:nx), cmlo(1:nx,1:ndimMTO)), evecmlo) !zPMT* C * zMTO
!           zcplz = cmlo(1:nx,1:ndimMTO), evecmlo) !  |F^MLO j> = |Psi^PMT i> cpmo(i,j)
!                                             !MLO eigenfunctions constructed from the PMT basis.
!           !   do i=1,ndimMTO
!           !     write(stdo,ftox)'eigen111',i,ftof(qp,3),'  ',ftof(evl(i+nskip)),' ',ftof(evlx(i))
!           !   enddo
!         endblock diagonalizeMTO
!       endif GetZCZ
       
     endblock ModifyMatrixElements
     return
   end subroutine Hreduction
   real(8) function fermidist(x)
     real(8),intent(in) :: x
     if(x>100d0) then
       fermidist=0d0
     elseif(x<-100d0) then
       fermidist=1d0
     else
       fermidist=1d0/(exp(x)+1)
     endif
   end function fermidist
   subroutine GramSchmidt(nv,n,zmel)
     integer:: igb=1,it,itt,n,nv
     complex(8):: ov(n),vec(nv),dnorm2(nv),zmel(nv,n)
     real(8):: dnorm
     do it = 1,n
       vec(:)= zmel(:,it)
       do itt = 1,it-1
         ov(itt) = sum( dconjg(zmel(:,itt))*vec(:))
       enddo
       vec = vec - matmul(zmel(:,1:it-1),ov(1:it-1))
       zmel(:,it) = vec/sum(dconjg(vec)*vec)**.5d0
     enddo
   end subroutine GramSchmidt
end module m_hreduction
