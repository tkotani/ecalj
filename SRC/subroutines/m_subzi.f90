module m_subzi ! Obtain weight wtkb(ib,isp,iq) for brillowine zone integation
  use m_ftox
  use m_lgunit,only: stdo
  real(8),allocatable,protected :: wtkb(:,:,:) ! wtkb : tetrahedron integration weights. it might be from wkp.*
  integer,protected:: nevmx=0 !for band mode plbnd=1. --->nev=min(nevmx,ndimhx) =0
contains
  subroutine m_subzi_init() !! Brillouin-integration setup
    use m_ext,only: sname
    use m_lmfinit, only: nsp,nspc, nevmxin=>bz_nevmx, lmet=>bz_lmet, qbg=>zbak,lso,nspx
    use m_mkqp,only: ntet=> bz_ntet ,bz_nkp
    use m_suham,only: ndham=>ham_ndham,ndhamx=>ham_ndhamx
    use m_mkpot,only:  qval
    !i Inputs
    !   lmet  :See Remarks
    !         :0 assume insulator
    !         :3 always make weights
    !     (we have removed lmet=1,2,4)
    !   ltet  :T for tetrahedron weights
    !   ndham :leading dimension of wtkb
    !   nsp   :2 for spin-polarized case, otherwise 1
    !   nkp   :number of irreducible k-points (bzmesh.f)
    !Outputs
    !o  nevmx :On input, maximum number of eigenvectors to find
    !          input nevmx<0 => do not generate eigenvectors
    !          input nevmx=0 => choose a default value
    !r   zval  :total valence charge
    !r   lmet:
    !r     lmet=0 system assumed to be an insulator; weights known a priori
    !r     lmet=3 After m_bandcal_init, we calcutlate wtkb, which is passed to m_bandcal_2nd to accumulate quantities of sum in BZ.
    implicit none
    logical :: ltet,lwt,tdos,cmdopt0,PROCARon,fullmesh
    integer :: nkp,mpsord,ifi,lerr,iprint,n,nevx0,nq0,nsp0
    real(8) :: zval,ef0,def,esmear
    call tcn('m_subzi_init')
    if(cmdopt0('--cls') .OR. cmdopt0('--tdos') .OR. cmdopt0('--mkprocar'))  nevmx = ndhamx
    fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface')
    if(lso==1.and.(cmdopt0('--band') .OR. fullmesh)) then
       nevmx=ndhamx
       return
    endif
    if(cmdopt0('--band') .OR. fullmesh) return
    nkp=bz_nkp
    if(allocated(wtkb)) deallocate(wtkb)
    zval = qval-qbg
    ltet = ntet>0
    nevmx= nevmxin !initial contdition except  cmdopt0('--band').or.fullmesh
    if(lmet>0) allocate(wtkb(ndhamx,nspx,nkp))
    if(nevmx == 0) then
       nevmx = (int(zval) + 1)/2
       if(lmet /= 0) nevmx = max(nevmx+nevmx/2,9)
       nevmx = min(nevmx,ndham)
       if(nspc == 2) nevmx = 2*nevmx
       nevmx=nevmx+5 !+5 is for safer setting. At least +1 is required...
    endif
    call tcx('m_subzi_init')
  end subroutine m_subzi_init
  subroutine m_subzi_bzintegration(evlall,efermi,sev,qvalm,vmag)
    use m_bzintegration2,only: bzintegration2
    implicit none
    intent(in)::                   evlall
    intent(out)::                         efermi,sev,qvalm,vmag
    real(8):: evlall(*),sev, qvalm(2),efermi,vmag
    call bzintegration2(evlall, efermi,sev,wtkb,qvalm,vmag)
  end subroutine m_subzi_bzintegration
end module m_subzi
