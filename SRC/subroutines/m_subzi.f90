module m_subzi ! Obtain weight wtkb(ib,isp,iq) for brillowine zone integation
  use m_ftox
  use m_lgunit,only: stdo
  real(8),allocatable,protected :: wtkb(:,:,:) ! wtkb : tetrahedron integration weights. it might be from wkp.*
  integer,protected:: nevmx
contains
  subroutine m_subzi_init() ! Set nevmx and allocate wtkb.
    use m_ext,only: sname
    use m_lmfinit, only: nspc,lmet=>bz_lmet, qbg=>zbak,nspx
    use m_mkqp,only: ntet=> bz_ntet ,bz_nkp
    use m_suham,only: ndhamx=>ham_ndhamx
    use m_mkpot,only:  qval
    !   lmet/=0 : allocate tetrahedron weight wtkb
    !   ndhamx : leading dimension of wtkb
    !   nsp    : 2 for spin-polarized case, otherwise 1
    !   nkp    : number of irreducible k-points (bzmesh.f)
    !   nevmx  : maximum number of eigenvectors to find 
    implicit none
    logical :: cmdopt0
    integer :: nkp
    real(8) :: zval
    call tcn('m_subzi_init')
    if(lmet>0) then
      nkp  = bz_nkp
      if(allocated(wtkb)) deallocate(wtkb)       
      allocate(wtkb(ndhamx,nspx,nkp))
    endif
    if(cmdopt0('--tdos').or. cmdopt0('--band').or.cmdopt0('--fermisurface')) then !nevmx=0 implies eigenvalue-only mode
      nevmx= 0
    elseif(cmdopt0('--pdos').or.cmdopt0('--mkprocar').or.cmdopt0('--zmel0').or.cmdopt0('--cls')) then
      nevmx= ndhamx  !all bands
    else  !above occipied bands. (tetrahedron method may require a little more than zval/2)
      zval = qval-qbg
      nevmx = ceiling(zval)/2
      if(lmet /= 0) nevmx = max(nevmx+nevmx/2,9) !probably safer setting nevmx for metal. Rough estimation.
      nevmx = min(nevmx*nspc, ndhamx) !nspc=2 for lso=1
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
