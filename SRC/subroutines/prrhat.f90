subroutine prrhat ( nbas , ssite , sspec , sv_p_orhoat )
  use m_struc_def
  use m_lmfinit,only: nsp
  use m_lgunit,only:stdo
  integer:: nbas
  type(s_rv1) :: sv_p_orhoat(3,nbas)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  ! ... Local parameters
  integer:: ib , lmxl , nr , nlml , is , igetss
  real(8) ,allocatable :: rofi_rv(:)
  real(8) ,allocatable :: rwgt_rv(:)
  double precision :: a,rmt
  do  ib = 1, nbas
     is = int(ssite(ib)%spec)
     lmxl=sspec(is)%lmxl
     a=sspec(is)%a
     nr=sspec(is)%nr
     rmt=sspec(is)%rmt
     if (lmxl == -1) goto 10
     allocate(rofi_rv(nr))
     allocate(rwgt_rv(nr))
     call radmsh ( rmt , a , nr , rofi_rv )
     call radwgt ( rmt , a , nr , rwgt_rv )
     nlml = (lmxl+1)**2
     write(stdo,200) ib,rmt,nr,nlml
200  format(/' Density at site',i3,'   rmt=',f8.4, &
          '   nr=',i5,'   nlml=',i3)
     call prlrho ( 'true density' , nr , nlml , nsp , rofi_rv , rwgt_rv &
          , sv_p_orhoat( 1 , ib )%v )
     call prlrho ( 'smooth density' , nr , nlml , nsp , rofi_rv , &
          rwgt_rv , sv_p_orhoat( 2 , ib )%v )
     if (allocated(rwgt_rv)) deallocate(rwgt_rv)
     if (allocated(rofi_rv)) deallocate(rofi_rv)
10   continue
  enddo
end subroutine prrhat

subroutine prlrho(str,nr,nlm,nsp,rofi,rwgt,rho)
  use m_lgunit,only:stdo
  !- Print info about site density
  !     implicit none
  ! ... Passed parameters
  integer :: nr,nlm,nsp
  double precision :: rho(nr,nlm,nsp),rofi(nr),rwgt(nr)
  character*(*) str
  ! ... Local parameters
  integer :: ilm,l,ll,itop,ibot,i
  double precision :: xx,pi,srfpi,top,bot,sum
  !      stdo = lgunit(1)
  pi = 4d0*datan(1d0)
  srfpi = dsqrt(4d0*pi)
  write(stdo,601) str
  do  ilm = 1, nlm
     l = ll(ilm)
     top = -1d10
     bot = 1d10
     sum = 0d0
     itop = 0
     ibot = 0
     do  i = 1, nr
        xx = (rho(i,ilm,1)+rho(i,ilm,nsp))/(3-nsp)
        if (xx > top) then
           itop = i
           top = xx
        endif
        if (xx < bot) then
           ibot = i
           bot = xx
        endif
        sum = sum + rwgt(i)*xx*(rofi(i)+1d-32)**l
     enddo
     xx = dmax1(dabs(top),dabs(bot))
     if (xx > 1d-6) then
        if (ilm == 1) then
           write(stdo,600) ilm,bot,ibot,top,itop,sum,srfpi*sum
        else
           write(stdo,600) ilm,bot,ibot,top,itop,sum
        endif
     endif
600  format(i4,f15.8,i5,f15.8,i5,f15.8,f12.6)
601  format(/' prlrho: ',a/'  ilm',6x,'rhomin    pnt', &
          7x,'rhomax    pnt       moment')
  enddo
end subroutine prlrho


