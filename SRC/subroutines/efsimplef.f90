subroutine efsimplef2ax ( legas, esmr, &
     valn,ef)
  use m_READ_BZDATA,only: &
       nqbz,nqibz,ginv,  qibz,wibz,qbz
  use m_genallcf_v3,only: nspin,z,natom,nclass,iclass,nl,konfig=>konf
  use m_readeigen, only: readeval
  use m_readhbe,only: nband
  use m_mpi, only: mpi__root
  use m_hamindex,only:   zbak
  !!== Calculate efermi for discrete sum. (not for tetrahedron method) ==
  !! You need to call init_reaeigen before you call this.
  !! user readeval (readeigen.f) to get eigenvalues.
  ! nspin   = 1, paramagnetic
  !           2, ferromagnetic
  ! ef      = fermi level
  ! nband   = no. states
  ! nqbz    = no. k-points
  ! valn    = number of valence electron.

  ! -------------------
  !      e(iband) < efm : occupation is one
  ! efm< e(iband) < efp : occupation is wfacef.
  ! efp< e(iband)       : occupation is zero

  implicit none
  integer(4):: is,iq,ik,isig,kpx,ifev(2) !nspin,nqibz
  integer(4):: ieaord(nband*nqibz*nspin),mbytes,mwords,iwksize, &
       iqibz  !,iindxk
  real(8)   :: ekt(nband, nqibz,nspin), ektx(nband*nqibz*nspin)
  real(8)   :: wgt(nband, nqibz,nspin), wgtx(nband*nqibz*nspin)
  real(8)   :: qx(3),qbas(3,3),wwg !qbzx(3), ,ginv(3,3),

  integer(4):: ncore,l,ia,ic ,ierr
  !      integer(4):: konfig(0:nl-1,nclass) !konf(nl,nclass)
  real(8)   :: valn,ef

  integer(4) :: nbnqnsp,ix,ikx=-9999,ikini,nne
  real(8)    :: ew1,ew2,ein,valx,enumef_gauss,esmr, efini &
       ,eee2,wwg2 ,enumef
  !      real(8) :: efp,efm,wwgo,wfacef
  logical :: legas,autoew,GaussSmear=.true. !is external

  integer(4):: ifile_handle,if8301,if8302 !nqbz,
  !      integer(4):: n_index_qbz,index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)
  !     real(8)   :: qbz(3,nqbz)
  !      integer:: ifzbak
  !      real(8)::zbak
  !--------------------------------------------------------------------
  autoew =.false.
  if(GaussSmear) then
     write(6,*)' efsimplef2(gaussian mode):start'
  else
     write(6,*)' efsimplef2(rectangular mode):start'
  endif
  if(esmr<=0d0) autoew= .TRUE. 
  ! total valence charge
  if(legas) then
     write(6,*)' efsimplef2: legas=T use given valn = ',valn
  else
     valn    = 0d0
     do ia   = 1,natom
        ic    = iclass(ia)
        valn  = valn + z(ic)
        write(6,*)' ia z(ic)=',ia, z(ic)
        do    l = 0,nl-1
           write(6,*)' l (konfig(l+1,ic)-l-1) 2*(2l+1)=',l,(konfig(l+1,ic)-l-1),( 2*l +1)*2
           valn  = valn - (konfig(l+1,ic)-l-1) *( 2*l +1)*2
        end do
     end do
  endif
  !! read ZBAK file
  !      open(newunit=ifzbak,file='ZBAK')
  !      read(ifzbak,*)zbak
  valn=valn-zbak
  !      close(ifzbak)
  !      valn = valn+zbak
  !      write(6,*)' valn=valn+zbak =',valn

  do is = 1,nspin
     do iq = 1,nqibz
        ekt(:,iq,is) = readeval(qibz(:,iq),is)
     enddo
  enddo

  if(abs(sum(wibz(1:nqibz))-2d0)>1d-10) then
     write(6,*) 'sum (wibz)=', sum(wibz(1:nqibz))
     ! top2rx 2013.08.09 kino        stop 'efsimplef2: wibzsumerr'
     call rx( 'efsimplef2: wibzsumerr')
  endif
  do is = 1,nspin
     do iq = 1,nqibz
        wgt(1:nband,iq,is) = wibz(iq)
        if(nspin==2) wgt(1:nband,iq,is) = wgt(1:nband,iq,is)/2d0
     enddo
  enddo

  ! ekt and wgt
  call dcopy ( nband*nqibz*nspin, ekt,1, ektx ,1)
  call dcopy ( nband*nqibz*nspin, wgt,1, wgtx ,1)
  call sortea( ektx,ieaord,nband*nqibz*nspin,isig)

  !      write(6,*)nband,nqibz,nspin,nband*nqibz*nspin
  !      do ik  = 1, nband*nqibz*nspin
  !        write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wgtx(ieaord(ik))
  !      enddo

  if (mpi__root) then
     !         if8301=ifile_handle()
     open(newunit=if8301,file = "DOSACC.lda")
     !         if8302=ifile_handle()
     open(newunit=if8302,file = "DOSACC2.lda")
  endif
  wwg = 0d0
  wwg2= 0d0
  eee2= -1d99
  ikini= 1
  ierr= 1
  nne = nband*nqibz*nspin
  do ik = 1, nne
     !---
     if(eee2 +1d-4< ektx(ieaord(ik)) .OR. ik==nne ) then
        ! degeneracy check
        if (mpi__root) then
           if(ik/=1) write(if8302,"(2i6,2d23.15)") ikini,ik-1,eee2,wwg2
        endif
        wwg2 = wgtx(ieaord(ik))
        eee2 = ektx(ieaord(ik))
        ikini =ik
     else
        wwg2= wwg2 + wgtx(ieaord(ik))
     endif
     !---
     wwg = wwg + wgtx(ieaord(ik))
     if(wwg<valn+2d0) write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wwg

     if (mpi__root) then
        write(if8301,"(2i6,3d23.15)") &
             ik,ieaord(ik),ektx(ieaord(ik)),wwg,wgtx(ieaord(ik))
     endif

     if( wwg>valn-1d-8 .AND. ierr==1 ) then
        write(6,*)
        efini = .5*(ektx(ieaord(ik+1))+ ektx(ieaord(ik)))
        if(autoew) then
           ! top2rx 2013.08.09 kino            if(ik<3) stop ' efsimplef2: ik<3'
           if(ik<3) call rx( ' efsimplef2: ik<3')
           esmr  = ektx(ieaord(ik)) - ektx(ieaord(ik-1))
        endif
        ierr=0
     endif
     if( wwg > valn+1d0) ikx=ik
  enddo

  ! top2rx 2013.08.09 kino      if(ierr==1) stop ' efsimplef2: ierr=1 given nval is too large'
  if(ierr==1) call rx( ' efsimplef2: ierr=1 given nval is too large')

  nbnqnsp = nband*nqibz*nspin
  ! c gaussian
  if(GaussSmear) then
     valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp)) &
          ,efini,esmr,nbnqnsp)
  else
     valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp)) &
          ,efini,esmr,nbnqnsp)
  endif
  write(6,*) 'valx at efini=',efini,valx
  if(abs(valx-valn)<1d-8) then
     ef=efini
     goto 8891
  endif

  ew1= ektx(ieaord(1))-0.01d0
  ew2= ektx(ieaord(ikx))
  nbnqnsp =nband*nqibz*nspin
  do ix = 1,100
     ein = 0.5d0*(ew1+ew2)
     if(GaussSmear) then
        valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp)) &
             ,ektx(ieaord(1:nbnqnsp)),ein,esmr,nbnqnsp)
     else
        valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp)) &
             ,ein,esmr,nbnqnsp)
     endif
     if(valx>valn) ew2=ein
     if(valx<valn) ew1=ein
     if(abs(ew1-ew2)<1d-15) exit
  enddo
  ef = 0.5d0*(ew1+ew2)

8891 continue
  if (mpi__root) then
     write(if8301,*) " ef=",ef
     close(if8301)
     write(if8302,*) " ef=",ef
     close(if8302)
  endif

  !      write(6,*)' esmr        =',esmr
  !      write(6,*)' determined ef =',ef
  !----------------------------------
  !      wwg = 0d0
  !      do ik  = 1, nband*nqibz*nspin
  !       wwgo = wwg
  !       wwg  = wwg + wgtx(ieaord(ik))
  !       if( abs(wwg-valn)<1d-6) then
  !         ef = 0.5d0*( ektx(ieaord(ik))+ektx(ieaord(ik+1)) )
  !         efp = ef + 0.25d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efp is just above the fermi
  !         efm = ef - 0.25d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efm is just below the fermi
  !        elseif(wwg>valn) then
  !         ef      = ektx(ieaord(ik))
  !         wfacef  = (valn-wwgo)/wgtx(ieaord(ik))
  !         efp = ef + 0.5d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efp is just above the fermi
  !         efm = ef - 0.5d0*(ektx(ieaord(ik))-ektx(ieaord(ik-1)))  !efm is just below the fermi
  !          write(6,*)' determined ef    =',ef
  !          write(6,*)'            efp   =',efp
  !          write(6,*)'            efm   =',efm
  !          write(6,*)'           wfacef =',wfacef
  !          return
  !        endif
  !     enddo
  if(GaussSmear) then
     write(6,*)' efsimplef2ax(gauss):end'
  else
     write(6,*)' efsimplef2ax:end'
  endif
end subroutine efsimplef2ax
!----------------------------------------------------------------------------
subroutine efsimplef2a (nspin,wibz,qibz,ginv, &
     nband,nqibz, konfig,z,nl,natom,iclass,nclass, &
     valn, legas, esmr, &! & input only for empty case
  qbz,nqbz, &! & index_qbz, n_index_qbz,
  ef)
  use m_readeigen, only: readeval
  use m_mpi, only: mpi__root
  use m_hamindex,only:   zbak
  !!== Calculate efermi for discrete sum. (not for tetrahedron method) ==
  !! You need to call init_reaeigen before you call this.
  !! user readeval (readeigen.f) to get eigenvalues.
  ! nspin   = 1, paramagnetic
  !           2, ferromagnetic
  ! ef      = fermi level
  ! nband   = no. states
  ! nqbz    = no. k-points
  ! valn    = number of valence electron.

  ! -------------------
  !      e(iband) < efm : occupation is one
  ! efm< e(iband) < efp : occupation is wfacef.
  ! efp< e(iband)       : occupation is zero

  implicit none
  integer(4):: is,iq,nspin,nqibz,ik,isig,kpx,nband,ifev(2)
  integer(4):: ieaord(nband*nqibz*nspin),mbytes,mwords,iwksize, &
       iqibz  !,iindxk
  real(8)   :: ekt(nband, nqibz,nspin), ektx(nband*nqibz*nspin)
  real(8)   :: wgt(nband, nqibz,nspin), wgtx(nband*nqibz*nspin)
  real(8)   :: qbzx(3), qx(3),qbas(3,3),ginv(3,3),wwg

  integer(4):: nclass,natom,nl,ncore,l,ia,ic !,indxk(*)
  real(8)   :: wibz(nqibz),valn,ef, z(nclass),qibz(3,nqibz)
  integer(4):: iclass(natom),konfig(0:nl-1,nclass),ierr

  integer(4) :: nbnqnsp,ix,ikx=-99999,ikini,nne
  real(8)    :: ew1,ew2,ein,valx,enumef_gauss,esmr, efini &
       ,eee2,wwg2 ,enumef
  !      real(8) :: efp,efm,wwgo,wfacef
  logical :: legas,autoew,GaussSmear=.true. !is external

  integer(4):: nqbz,if8301,if8302 !ifile_handle,
  !      integer(4):: n_index_qbz,index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)
  real(8)   :: qbz(3,nqbz)
  !      integer:: ifzbak
  !      real(8)::zbak
  !--------------------------------------------------------------------
  autoew =.false.
  if(GaussSmear) then
     write(6,*)' efsimplef2(gaussian mode):start'
  else
     write(6,*)' efsimplef2(rectangular mode):start'
  endif
  if(esmr<=0d0) autoew= .TRUE. 
  ! total valence charge
  if(legas) then
     write(6,*)' efsimplef2: legas=T use given valn = ',valn
  else
     valn    = 0d0
     do ia   = 1,natom
        ic    = iclass(ia)
        valn  = valn + z(ic)
        write(6,*)' ia z(ic)=',ia, z(ic)
        do    l = 0,nl-1
           write(6,*)' l (konfig(l,ic)-l-1) 2*(2l+1)=',l,(konfig(l,ic)-l-1),( 2*l +1)*2
           valn  = valn - (konfig(l,ic)-l-1) *( 2*l +1)*2
        end do
     end do
  endif
  ! ccccccccccccccccccccccccccccccccccccccccc
  !      valn=29d0
  ! cccccccccccccccccccccccccccccccccccccccc
  write(6,*)' valn=',valn
  !! read ZBAK file
  !      open(newunit=ifzbak,file='ZBAK')
  !      read(ifzbak,*)zbak
  valn=valn-zbak
  !      close(ifzbak)
  !      valn = valn+zbak
  !      write(6,*)' valn=vanl+zbak =',valn

  do is = 1,nspin
     do iq = 1,nqibz
        ekt(:,iq,is) = readeval(qibz(:,iq),is)
     enddo
  enddo

  if(abs(sum(wibz(1:nqibz))-2d0)>1d-10) then
     write(6,*) 'sum (wibz)=', sum(wibz(1:nqibz))
     call rx( 'efsimplef2: wibzsumerr')
  endif
  do is = 1,nspin
     do iq = 1,nqibz
        wgt(1:nband,iq,is) = wibz(iq)
        if(nspin==2) wgt(1:nband,iq,is) = wgt(1:nband,iq,is)/2d0
     enddo
  enddo
  ! ekt and wgt
  call dcopy ( nband*nqibz*nspin, ekt,1, ektx ,1)
  call dcopy ( nband*nqibz*nspin, wgt,1, wgtx ,1)
  call sortea( ektx,ieaord,nband*nqibz*nspin,isig)
  if (mpi__root) then
     !         if8301=ifile_handle()
     open(newunit=if8301,file = "DOSACC.lda")
     !         if8302=ifile_handle()
     open(newunit=if8302,file = "DOSACC2.lda")
  endif
  wwg = 0d0
  wwg2= 0d0
  eee2= -1d99
  ikini= 1
  ierr= 1
  nne = nband*nqibz*nspin
  do ik = 1, nne
     !---
     if(eee2 +1d-4< ektx(ieaord(ik)) .OR. ik==nne ) then
        ! degeneracy check
        if (mpi__root) then
           if(ik/=1) write(if8302,"(2i6,2d23.15)") ikini,ik-1,eee2,wwg2
        endif
        wwg2 = wgtx(ieaord(ik))
        eee2 = ektx(ieaord(ik))
        ikini =ik
     else
        wwg2= wwg2 + wgtx(ieaord(ik))
     endif
     !---
     wwg = wwg + wgtx(ieaord(ik))
     if(wwg<valn+2d0) write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wwg

     if (mpi__root) then
        write(if8301,"(2i6,3d23.15)") &
             ik,ieaord(ik),ektx(ieaord(ik)),wwg,wgtx(ieaord(ik))
     endif

     if( wwg>valn-1d-8 .AND. ierr==1 ) then
        write(6,*)
        efini = .5*(ektx(ieaord(ik+1))+ ektx(ieaord(ik)))
        if(autoew) then
           ! top2rx 2013.08.09 kino            if(ik<3) stop ' efsimplef2: ik<3'
           if(ik<3) call rx( ' efsimplef2: ik<3')
           esmr  = ektx(ieaord(ik)) - ektx(ieaord(ik-1))
        endif
        ierr=0
     endif
     if( wwg > valn+1d0) ikx=ik
  enddo

  ! top2rx 2013.08.09 kino      if(ierr==1) stop ' efsimplef2: ierr=1 given nval is too large'
  if(ierr==1) call rx( ' efsimplef2: ierr=1 given nval is too large')

  nbnqnsp = nband*nqibz*nspin
  ! c gaussian
  if(GaussSmear) then
     valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp)) &
          ,efini,esmr,nbnqnsp)
  else
     valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp)) &
          ,efini,esmr,nbnqnsp)
  endif
  write(6,*) 'valx at efini=',efini,valx
  if(abs(valx-valn)<1d-8) then
     ef=efini
     goto 8891
  endif

  ew1= ektx(ieaord(1))-0.01d0
  ew2= ektx(ieaord(ikx))
  nbnqnsp =nband*nqibz*nspin
  do ix = 1,100
     ein = 0.5d0*(ew1+ew2)
     if(GaussSmear) then
        valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp)) &
             ,ektx(ieaord(1:nbnqnsp)),ein,esmr,nbnqnsp)
     else
        valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp)) &
             ,ein,esmr,nbnqnsp)
     endif
     if(valx>valn) ew2=ein
     if(valx<valn) ew1=ein
     if(abs(ew1-ew2)<1d-15) exit
  enddo
  ef = 0.5d0*(ew1+ew2)

8891 continue
  if (mpi__root) then
     write(if8301,*) " ef=",ef
     close(if8301)
     write(if8302,*) " ef=",ef
     close(if8302)
  endif

  !      write(6,*)' esmr        =',esmr
  !      write(6,*)' determined ef =',ef
  !----------------------------------
  !      wwg = 0d0
  !      do ik  = 1, nband*nqibz*nspin
  !       wwgo = wwg
  !       wwg  = wwg + wgtx(ieaord(ik))
  !       if( abs(wwg-valn)<1d-6) then
  !         ef = 0.5d0*( ektx(ieaord(ik))+ektx(ieaord(ik+1)) )
  !         efp = ef + 0.25d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efp is just above the fermi
  !         efm = ef - 0.25d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efm is just below the fermi
  !        elseif(wwg>valn) then
  !         ef      = ektx(ieaord(ik))
  !         wfacef  = (valn-wwgo)/wgtx(ieaord(ik))
  !         efp = ef + 0.5d0*(ektx(ieaord(ik+1))-ektx(ieaord(ik)))  !efp is just above the fermi
  !         efm = ef - 0.5d0*(ektx(ieaord(ik))-ektx(ieaord(ik-1)))  !efm is just below the fermi
  !          write(6,*)' determined ef    =',ef
  !          write(6,*)'            efp   =',efp
  !          write(6,*)'            efm   =',efm
  !          write(6,*)'           wfacef =',wfacef
  !          return
  !        endif
  !     enddo
  if(GaussSmear) then
     write(6,*)' efsimplef2(gauss):end'
  else
     write(6,*)' efsimplef2:end'
  endif
end subroutine efsimplef2a
!----------------------------------------------------------------------------
!------------------------------------------------------
real(8) function enumef_gauss( wgtx,ektx,ein,esmr,nbnqnsp)
  implicit none
  integer(4):: nbnqnsp,ik
  real(8) :: ektx(nbnqnsp),wgtx(nbnqnsp),wwg, &
       derfcx,ein,esmr
  wwg = 0d0
  do ik = 1, nbnqnsp
     wwg= wwg + wgtx(ik) &
          *0.5d0* derfcx( -(ein-ektx(ik))/sqrt(2d0)/esmr )
  enddo
  enumef_gauss = wwg
  !     write(6,*)' ein enumef=', ein, enumef
END function enumef_gauss

real(8) function derfcx(a)
  real(8):: a, derfc,ax, amx=12d0
  ax = a
  if( abs(a)>amx) ax= a/abs(a) * amx
  !      write(6,*)' xxx ',ax
  !      write(6,*)' yyy ',derfc(ax)
  derfcx=derfc(ax)
END function derfcx

!------------
real(8) function enumef( wgtx,ektx,ein,esmr,nbnqnsp)
  implicit real*8(a-h,o-z)
  integer:: nbnqnsp,ik
  real(8) :: ektx(nbnqnsp),wgtx(nbnqnsp)
  !     write(6,*) esmr
  wwg = 0d0
  do ik = 1, nbnqnsp
     !       write(6,*)'ik=',ik,ektx(ik),wgtx(ik)
     if    (  ektx(ik) + 0.5d0*esmr < ein ) then
        wwg  = wwg + wgtx(ik)
     elseif(  ektx(ik) - 0.5d0*esmr < ein ) then
        wwg  = wwg + wgtx(ik)*(ein- (ektx(ik)-0.5d0*esmr))/esmr
     endif
  enddo
  enumef = wwg
  !     write(6,*)' ein enumef=', ein, enumef
END function enumef



!---------------------------------------------------------------
subroutine findemaxmin(nband,qbz,nqbz,nspin, &
     emax,emin)
  use m_readeigen, only: readeval
  implicit none
  integer(4) :: nband,nqbz,nspin,isp,kx,i !,ifev(2)
  real(8)::emax,emin,qbz(3,nqbz),eee
  real(8),allocatable:: ekxxx(:,:,:)
  allocate( ekxxx(nband,nqbz,nspin))
  Emax=-1d9
  do isp =1, nspin
     do kx = 1, nqbz
        ekxxx(1:nband,kx,isp) = readeval(qbz(:,kx), isp)
        do i=1,nband
           eee= ekxxx(i,kx,isp)
           !            print *,i,eee
           if(eee>Emax .AND. eee<1d9) Emax=eee !not eee<1d9 corresponds to 1d20 for padding in lmf2gw.F and sugw.Fago
        enddo
     enddo
  enddo
  !      Emax = maxval(ekxxx(1:nband)
  Emin = minval(ekxxx)
  deallocate(ekxxx)
end subroutine findemaxmin
!$$$      subroutine findemaxmin(ifev,nband,nqbz,nspin,
!$$$     o   emax,emin)
!$$$      implicit none
!$$$      integer(4) :: nband,nqbz,nspin,isp,kx,ifev(2)
!$$$      real(8)::emax,emin
!$$$      real(8),allocatable:: ekxxx(:,:,:)
!$$$      allocate( ekxxx(nband,nqbz,nspin))
!$$$      do isp =1, nspin
!$$$      do kx = 1, nqbz
!$$$        call rwdd1 (ifev(isp), kx, nband, ekxxx(1:nband,kx,isp) )
!$$$      enddo
!$$$      enddo
!$$$      Emax = maxval(ekxxx)
!$$$      Emin = minval(ekxxx)
!$$$      deallocate(ekxxx)
!$$$      end

!$$$      subroutine readomgc(ifinin,omg_c)
!$$$      real(8)::omg_c,blank
!$$$      integer::ifinin,is,iopen,iclose
!$$$      ifinin=iopen('GWIN_V2',1,0,0)
!$$$      read(ifinin,*)
!$$$      read(ifinin,*)
!$$$      read(ifinin,*) blank,omg_c !omg_c is frequency parameter
!$$$      ! freq2(iw)=dw(iw-1)+dw**2(iw-1)**2/2/omg_c
!$$$      ! quadratic term is essential for energy > omg_c
!$$$      is=iclose('GWIN_V2')
!$$$      end