!>Crystal symmetry data are stored by call m_mksym_init. NOTE:nbas (atomic sites)-> nspec (species) -> nclass (class)
module m_mksym 
  public :: m_mksym_init
  integer,allocatable,protected :: oics(:)    ! ispec= ics(iclass) gives spec for iclass.
  real(8),allocatable,protected :: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:),dlmm(:,:,:,:)
  integer,allocatable,protected :: invgx(:),miat(:,:), oistab (:,:)   ! j= istab(i,ig): site i is mapped to site j by grp op ig
  integer,allocatable,protected :: iclasst(:)
  integer,protected::              nclasst !number of equivalent class
  integer,protected :: ngrp  !# of lattice symmetry 
  integer,protected :: npgrp !# of lattice symmetry with artificial inversion (as time-reversal symmetry) if lc=1 =>mkqp
  logical,protected:: AFmode
  integer,protected:: ngrpAF
  !NOTE: We have ngrp+ngrpAF symmetry for symops... ngrp for lattice, ngrpAF is extra symmetry for AF (spin-flip symmetry)
  !  ngrp+ngrpAF is the all symmetry for SYMGRP+SYMGRPAF  (SYMGRPAF needs AF index for AFpairs
  !  SpinFlip symmetry: y=matmul(symopsF(:,:,ig),x)+ ag_af(:,ig), ig=ngrp+1,ngrp+ngrpAF.
  !NOTE for Antiferro symmetric self-energy mode.
  !To obtain self-consistency, it may be useful to keep AF condition during iteration. We need to set SYMGRPAF.
  ! This is a sample for NiO
  ! SYMGRPAF i:(1,1,1) ! translation + inversion for spin-flip symmetry.
  ! SYMGRP r3d          ! this keeps spin axis. (probably larger symmetry allowed for SO=0)
  ! STRUC  ALAT={a} PLAT= 0.5 0.5 1.0  0.5 1.0 0.5  1.0 0.5 0.5
  !        NBAS= 4  NSPEC=3
  ! SITE   ATOM=Niup POS=  .0   .0   .0   AF=1    <--- AF symmetric pair
  !        ATOM=Nidn POS= 1.0  1.0  1.0   AF=-1   <--- 
  !        ATOM=Mnup POS=  .0   .0   .0   AF=2   (if another AF pairs)
  !        ATOM=Mndn POS= 1.0  1.0  1.0   AF=-2 
  !        ATOM=O POS=  .5   .5   .5
  !        ATOM=O POS= 1.5  1.5  1.5
contains
  subroutine m_mksym_init(prgnam) 
    use m_mpitk,only:   master_mpi
    use m_lgunit,only:  stdo
    use m_lmfinit,only: nspec,nbas, sstrnsymg,symgaf,ips=>iv_a_oips,slabl,iantiferro,addinv,lmxax
    use m_lattic,only:  plat=>lat_plat,rv_a_opos
    use m_mksym_util,only: mksym,mptauof,rotdlmm
    use m_ftox
    implicit none
    integer,parameter::  ngmx = 48
    character,intent(in)::  prgnam*(*)
    integer:: ibas,lc,j,iprint,nclass,ngrpTotal,k,npgrpAll,ig
    integer,parameter::recln=511
    logical ::cmdopt0,ipr10=.false.
    character strn*(recln),strn2*(recln),outs(recln)
    real(8):: osymgr (3,3,ngmx),oag(3,ngmx)
    real(8),parameter:: tol=1d-4
    integer,allocatable:: iclasstAll(:)
    call tcn('m_mksym_init')
    if(.NOT.prgnam=='LMFA') ipr10= iprint()>10 !this is only for master
    strn = 'find'
    if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
    if(cmdopt0('--nosym') .OR. cmdopt0('--pdos') ) strn = ' '
    lc=merge(1,0,addinv) ! Add inversion to get sampling k points. phi*. When we have TR with keeping spin \sigma, psi_-k\sigm(r) = (psi_k\sigma(r))^* 
    !lmxax=lmxax
    if(master_mpi) call pshpr(60)
    allocate(iclasst(nbas),oics(nbas),oistab(nbas,ngmx))
    if(ipr10) write(stdo,"(a)")   'SpaceGroupSym of Lattice: ========start========================== '
    write(stdo,"(a)") ' SYMGRP = '//trim(strn)
    call mksym(lc,slabl,strn,ips, iclasst,nclasst,npgrp,ngrp,oag,osymgr,oics,oistab)
    if(ipr10) write(stdo,"(a)") 'SpaceGroupSym of Lattice: ========end =========================== '
    allocate(symops,source=osymgr)
    allocate(ag,source=oag)
    AFmodeBlock: block
      real(8):: osymgrAll(3,3,ngmx),oagAll(3,ngmx)
      integer:: oicsAll(nbas),oistabAll(nbas,ngmx),ipsAF(nbas),iga,igall,nclassAll,ig
      AFmode=len_trim(symgaf)>0 
      if(AFmode) then
         strn2=trim(strn)//' '//trim(symgaf)
         if(ipr10) then
            write(stdo,*)
            write(stdo,"(a)") 'Add SpaceGroupSym ops by AF symmetry===start========= '
            write(stdo,"(a)") 'AF: Antiferro mode: SYMGRPAF    = '//trim(symgaf)
            write(stdo,"(a)") 'AF:                 SYMGRPAF all= '//trim(strn2)
            write(stdo,"(a,2i3)")  ('AF:  ibas,AF=',j,iantiferro(j),j=1,nbas)
         endif
         ipsAF = ips
         do j=1,nbas
            do k=j,nbas
               if( iantiferro(j)+iantiferro(k)==0) then
                  ipsAF(k) = ipsAF(j) !to drive mksym for AF mode (Assuming AF pairs with the same spec).
                  exit
               endif
            enddo
         enddo
         if(master_mpi) call pshpr(60)
         allocate(iclasstAll(nbas))
         if(ipr10) write(stdo,"(a)")   'SpaceGroupSym of Lattice+AF: ========start========================== '
         call mksym(lc,slabl,strn2,ipsAF, iclasstAll,nclassAll, npgrpAll,ngrpTotal,oagAll,osymgrAll,oicsAll,oistabAll) !Big symmetry for lattice+AF
         ngrpAF=ngrpTotal-ngrp
         if(ipr10) write(stdo,"(a)")   'SpaceGroupSym of Lattice+AF: ========end========================== '

!         do ig=1,ngrpTotal
!            write(stdo,ftox)'symall ig=',ig,ftof(reshape(osymgrAll(:,:,ig),[9]),2),' ',oagAll(:,ig)
!         enddo
      
         iga=ngrp
         do igall=1,ngrpTotal !Pick up symmetry by AF
            if(any( [( sum(abs(symops(:,:,ig)-osymgrALL(:,:,igall)))+sum(abs(ag(:,ig)-oagAll(:,igall)))<tol,ig=1,ngrp )]  )) cycle
            iga=iga+1
            symops(:,:,iga)= osymgrALL(:,:,igall)
            ag(:,iga)      = oagALL(:,igall)
            oistab(:,iga)  = oistabAll(:,igall)
         enddo
         if(iga/=ngrpAF+ngrp) call rxiii('ngrpAF+nggp/=ngrpTotal',ngrp,ngrpAF,iga)
      else
         ngrpTotal=ngrp
         ngrpAF=0
      endif
    endblock AFmodeBlock
    MiatTiatDlmm:block
      allocate(miat(nbas,ngrpTotal),tiat(3,nbas,ngrpTotal),invgx(ngrpTotal),shtvg(3,ngrpTotal),&
           dlmm(-lmxax:lmxax,-lmxax:lmxax,0:lmxax,ngrpTotal))
      call            mptauof(symops,             ngrp,  plat,nbas,rv_a_opos,iclasst,   miat,tiat,invgx,shtvg)  !for ig=1,ngrp
      if(AFmode) call mptauof(symops(:,:,ngrp+1:),ngrpAF,plat,nbas,rv_a_opos,iclasstAll, & !  ig=ngrp+1,ngrpAF
           miat(:,ngrp+1:),tiat(:,:,ngrp+1:),invgx(ngrp+1:),shtvg(:,ngrp+1:),afmode) !mapping of sites by spacegrope ops
!      write(stdo,ftox)'mmmmm iclasst=',iclasst
!      do ig=1,ngrp
!         write(stdo,ftox)'mmmm ig=',ig, 'miat=',miat(1:nbas,ig)
!      enddo
!      write(stdo,ftox)'mmmmm iclasstAll=',iclasstAll
!      do ig=ngrp+1,ngrp+ngrpAF
!         write(stdo,ftox)'mmmm ig=',ig,'symops=',ftof(reshape(symops(:,:,ig),[9])),'miat=',miat(1:nbas,ig)
!      enddo
      call rotdlmm(symops,ngrpTotal, lmxax+1, dlmm) ! Get rotation matrix Dlmm in real spherical harmonics.  !for sigm mode, dlmm needed.
      if(ipr10) write(stdo,*)
    endblock MiatTiatDlmm
    call tcx('m_mksym_init')
  end subroutine m_mksym_init
end module m_mksym
