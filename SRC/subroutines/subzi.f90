module m_subzi ! Obtain weight wtkb(ib,isp,iq) for brillowine zone integation
  real(8),allocatable,protected :: rv_a_owtkb(:,:,:) ! wtkb : tetrahedron integration weights. it might be from wkp.*
  integer,protected:: nevmx=0 !for band mode plbnd=1. --->nev=min(nevmx,ndimhx) =0
!  integer,protected:: lswtk=0
contains
  subroutine m_subzi_init()
    use m_lgunit,only:stdo
    use m_ext,only: sname
    use m_lmfinit, only: nsp,nspc, nevmxin=>bz_nevmx, lmet=>bz_lmet, qbg=>zbak,lso
    use m_mkqp,only: ntet=> bz_ntet ,bz_nkp
    use m_suham,only: ndham=>ham_ndham,ndhamx=>ham_ndhamx
    use m_mkpot,only:  qval
!    use m_MPItk,only: master_mpi
    !- Brillouin-integration setup
    !i Inputs
    !i   lmet  :See Remarks
    !i         :0 assume insulator
    !i         :3 always make weights
    !i   ltet  :T for tetrahedron weights
    !i   ndham :leading dimension of owtkb
    !i         :Hamiltonian should not exceed this dimension
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nkp   :number of irreducible k-points (bzmesh.f)
    ! o Inputs/Outputs
    ! o  nevmx :On input, maximum number of eigenvectors to find
    ! o         input nevmx<0 => do not generate eigenvectors
    ! o         input nevmx=0 => choose a default value
    !xxxo   swtk  :memory is allocated for spin weights, nspc=2
    !!
    !r   To integrate the output density in the Brillouin Zone, integration
    !r   weights are needed for each band and qp, but they are not known
    !r   until all the bands are obtained.  The problem is solved in one of
    !r   the following ways:
    !r -------------------------------------
    !r   zval  :total valence charge
    !r     lmet=0 system assumed to be an insulator; weights known a priori
    !r
    !r     (removed) lmet=1 eigenvectors are written to disk, in which case the
    !r            integration for the charge density can be deferred until
    !r            all the bands are obtained
    !r
    !r     (removed) lmet=2 integration weights are assumed from a prior band pass
    !r
    !r     lmet=3 After m_bandcal_init, we calcutlate wtkb, which is passed to m_bandcal_2nd to accumulate quantities of sum in BZ.
    !r
    !      (removed) lmet=4 is removed now ------------.
    implicit none
    logical :: ltet,lwt,tdos,cmdopt0,PROCARon,fullmesh
    integer :: nkp,mpsord
    double precision :: zval,ef0,def,esmear
    integer :: ifi,lerr,iprint,n
    character(11) :: strni(2)
    integer :: nevx0,nq0,nsp0
    data strni /'sampling','tetrahedron'/
    call tcn('m_subzi_init')
    if(cmdopt0('--cls') .OR. cmdopt0('--tdos') .OR. cmdopt0('--mkprocar'))  nevmx = ndhamx
    fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface')
    if(cmdopt0('--band') .OR. fullmesh) return
    nkp=bz_nkp
    if(allocated(rv_a_owtkb)) deallocate(rv_a_owtkb)
    zval=qval-qbg
    ltet = ntet>0
    nevmx=nevmxin !initial contdition except  cmdopt0('--band').or.fullmesh
    if(lmet>0) then
      if(lso==1) allocate(rv_a_owtkb(ndhamx,1,nkp))
      if(lso/=1) allocate(rv_a_owtkb(ndham,nsp,nkp))
    endif  
    if(nevmx == 0) then
       nevmx = (int(zval) + 1)/2
       if (lmet /= 0) nevmx = max(nevmx+nevmx/2,9)
       nevmx = min(nevmx,ndham)
       if (nspc == 2) nevmx = 2*nevmx
       nevmx=nevmx+5 !+5 is for safer setting. At least +1 is required...
    endif
    call tcx('m_subzi_init')
  end subroutine m_subzi_init
  subroutine m_subzi_bzintegration(evlall, eferm,sev,sumqv,vmag)
    use m_ftox
    use m_lgunit,only:stdo
    use m_MPItk,only:  numprocs=>nsize
    use m_lmfinit,only: lso,nsp,ham_scaledsigma,nlibu,lmaxu,bz_w,lmet=>bz_lmet,nbas,epsovl=>ham_oveps,nspc,bz_n
    use m_lmfinit,only: bz_fsmommethod,qbg=>zbak,fsmom=>bz_fsmom,ndos=>bz_ndos
    use m_mkqp,only: nkabc=> bz_nabc,ntet=> bz_ntet,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet
    use m_lmfinit,only: nchan=>pot_nlma, nvl=>pot_nlml
    use m_qplist,only: m_qplist_init, nkp,xdatt,labeli,labele,dqsyml,etolc,etolv, nqp2n_syml,nqp_syml,nqpe_syml,nqps_syml,nsyml
    use m_mkpot,only: qval
    use m_suham,only: ndham=>ham_ndham,ndhamx=>ham_ndhamx,nspx=>ham_nspx
    use m_ext,only: sname
    use m_bzwts,only: bzwtsf,bzwtsf2
    implicit none
    intent(in)::                   evlall
    intent(out)::                          eferm,sev,sumqv,vmag
    logical:: lfill=.false.,ltet
    logical:: debug,cmdopt0 
    integer:: ierr,ifimag,i,ifi,unlink !,iobzwt
    real(8):: dosrng,evlall(*),sev,sumqv(3,*),eferm,vmag,ef0,bz_ef
    real(8),parameter::    NULLR =-99999
    call tcn('m_subzi_bzintegration')
    write(stdo,ftox)'m_subzi_bzintegration:'
    sev=0d0
    ltet = ntet>0
    debug = cmdopt0('--debugbndfp')
    ! --- BZ integration for fermi level, band sum and qp weights ---
    dosrng = 8
    if(bz_n<0) dosrng = 16
    if(bz_fsmommethod == 1) then ! vmag (in Ry) contains magnetic field.  For eigenvalus, add -vmag/2 for isp=1, and +vmag/2 for isp=2.
       call bzwtsf2 ( ndham,ndham,nsp,nspc,nkabc(1),nkabc(2),nkabc(3),nkp,ntet,iv_a_oidtet,qval-qbg,& ! & note qval is output
            fsmom,lmet.ne.0,ltet,bz_n,ndos,bz_w,dosrng,rv_a_owtkp,evlall,eferm,sev,rv_a_owtkb,sumqv(1,2),lfill,vmag)!, lwtkb&! lswtk,swtk &
    else
       call bzwtsf ( ndham,ndham,nsp,nspc,nkabc(1),nkabc(2),nkabc(3),nkp,ntet,iv_a_oidtet,qval-qbg,&
            fsmom,lmet.ne.0,ltet,bz_n,ndos,bz_w,dosrng,rv_a_owtkp,evlall,eferm,sev,rv_a_owtkb,sumqv(1,2),lfill,vmag) !, lwtkb  & ! lswtk,swtk &
    endif
    sumqv(:,1) = sumqv(:,2)
    call tcx('m_subzi_bzintegration')
  end subroutine m_subzi_bzintegration
end module m_subzi
