module m_subzi ! Obtain weight wtkb(ib,isp,iq) for brillowine zone integation
  real(8),allocatable,protected :: rv_a_owtkb(:,:,:) ! wtkb : tetrahedron integration weights. it might be from wkp.*
  integer,protected:: nevmx=0 !for band mode plbnd=1. --->nev=min(nevmx,ndimhx) =0
contains
  subroutine m_subzi_init() !! Brillouin-integration setup
    use m_lgunit,only:stdo
    use m_ext,only: sname
    use m_lmfinit, only: nsp,nspc, nevmxin=>bz_nevmx, lmet=>bz_lmet, qbg=>zbak,lso,nspx
    use m_mkqp,only: ntet=> bz_ntet ,bz_nkp
    use m_suham,only: ndham=>ham_ndham,ndhamx=>ham_ndhamx
    use m_mkpot,only:  qval
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
    !r   To integrate the output density in the Brillouin Zone, integration
    !r   weights are needed for each band and qp, but they are not known
    !r   until all the bands are obtained.  The problem is solved in one of
    !r   the following ways:
    !r -------------------------------------
    !r   zval  :total valence charge
    !r     lmet=0 system assumed to be an insulator; weights known a priori
    !r     lmet=3 After m_bandcal_init, we calcutlate wtkb, which is passed to m_bandcal_2nd to accumulate quantities of sum in BZ.
    !r
    !r     (removed) lmet=1 eigenvectors are written to disk, in which case the
    !r            integration for the charge density can be deferred until
    !r            all the bands are obtained
    !r     (removed) lmet=2 integration weights are assumed from a prior band pass
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
    zval = qval-qbg
    ltet = ntet>0
    nevmx= nevmxin !initial contdition except  cmdopt0('--band').or.fullmesh
    if(lmet>0) allocate(rv_a_owtkb(ndhamx,nspx,nkp))
    if(nevmx == 0) then
       nevmx = (int(zval) + 1)/2
       if(lmet /= 0) nevmx = max(nevmx+nevmx/2,9)
       nevmx = min(nevmx,ndham)
       if(nspc == 2) nevmx = 2*nevmx
       nevmx=nevmx+5 !+5 is for safer setting. At least +1 is required...
    endif
    call tcx('m_subzi_init')
  end subroutine m_subzi_init
  subroutine m_subzi_bzintegration(evlall, eferm,sev,sumqv,vmag)
    use m_ftox
    use m_lgunit,only: stdo
    use m_lmfinit,only: qbg=>zbak
    use m_mkpot,only: qval
    use m_mkqp,only: nkabc=> bz_nabc,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet
    use m_qplist,only: m_qplist_init, nkp,xdatt,labeli,labele,dqsyml,etolc,etolv, nqp2n_syml,nqp_syml,nqpe_syml,nqps_syml,nsyml
    use m_suham,only: ndham=>ham_ndham!,ndhamx=>ham_ndhamx,nspx=>ham_nspx
    use m_bzwts,only: bzwtsf2
    implicit none
    intent(in)::                   evlall
    intent(out)::                          eferm,sev,sumqv,vmag
    logical:: lfill=.false.,ltet, debug,cmdopt0 ,wtsf2
    integer:: ierr,ifimag,i,ifi,unlink !,iobzwt
    real(8):: dosrng,evlall(*),sev,sumqv(3,*),eferm,vmag,ef0,bz_ef
    real(8),parameter::    NULLR =-99999
    call tcn('m_subzi_bzintegration')
    write(stdo,ftox)'m_subzi_bzintegration:'
    call bzwtsf2(ndham,nkabc(1),nkabc(2),nkabc(3),nkp,iv_a_oidtet,qval-qbg,& ! & note qval is output
         rv_a_owtkp, evlall,eferm,sev,rv_a_owtkb,sumqv(1:2,2),lfill,vmag)!, lwtkb&! lswtk,swtk &
    sumqv(1:2,1) = sumqv(1:2,2)
    call tcx('m_subzi_bzintegration')
  end subroutine m_subzi_bzintegration
end module m_subzi
