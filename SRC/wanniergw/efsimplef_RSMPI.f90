subroutine efsimplef2a_RSMPI(ispin,wibz,qibz,ginv, &
     nband,nqibz, &
     konfig,z,nl,natom,iclass,nclass, &
     valn, legas, esmr, qbz,nqbz, ef)
  use m_readeigen,only:readeval
  use rsmpi !RS
  !- Calculate efermi for discrete sum. (not for tetrahedron method)
  !r user readeval (readeigen.f) to get eigenvalues.
  !r You need to call init_reaeigen before you call this.
  ! ispin   = 1, paramagnetic
  !           2, ferromagnetic
  ! ef      = fermi level
  ! nband   = no. states
  ! nqbz    = no. k-points
  ! valn    = number of valence electron.

  ! -------------------
  !      e(iband) < efm : occupation is one
  ! efm< e(iband) < efp : occupation is wfacef.
  ! efp< e(iband)       : occupation is zero

  !------------------------------------------------------------
  ! May 2007 Rei Sakuma : modified for MPI version
  !    only process with
  !        Is_IO_Root_RSMPI()==.true.(see gwsrc/RSMPI_mod.F)
  !      prints  messages and open output files.
  !------------------------------------------------------------

  implicit none
  integer(4):: is,iq,ispin,nqibz,ik,isig,kpx,nband,ifev(2)
  integer(4):: ieaord(nband*nqibz*ispin),mbytes,mwords,iwksize, &
       iqibz  !,iindxk
  real(8)   :: ekt(nband, nqibz,ispin), ektx(nband*nqibz*ispin)
  real(8)   :: wgt(nband, nqibz,ispin), wgtx(nband*nqibz*ispin)
  real(8)   :: qbzx(3), qx(3),qbas(3,3),ginv(3,3),wwg

  integer(4):: nclass,natom,nl,ncore,l,ia,ic !,indxk(*)
  real(8)   :: wibz(nqibz),valn,ef, z(nclass),qibz(3,nqibz)
  integer(4):: iclass(natom),konfig(0:nl-1,nclass),ierr

  integer(4) :: nbnqnsp,ix,ikx,ikini,nne
  real(8)    :: ew1,ew2,ein,valx,enumef_gauss,esmr, efini &
       ,eee2,wwg2 ,enumef
  !      real(8) :: efp,efm,wwgo,wfacef
  logical :: legas,autoew ,GaussSmear=.true. !is external

  integer(4):: nqbz,if8301,if8302,ifile_handle
  !      integer(4):: n_index_qbz,index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)
  real(8)   :: qbz(3,nqbz)

  !--------------------------------------------------------------------
  autoew =.false.
  if (Is_IO_Root_RSMPI()) then
     if(GaussSmear) then
        write(6,*) ' efsimplef2(gaussian mode):start'
     else
        write(6,*) ' efsimplef2(rectangular mode):start'
     endif
  endif

  if(esmr<=0d0) autoew= .TRUE. 
  ! total valence charge
  if(legas) then
     ! RS: legas mode is not implemented..
     call RSMPI_Stop(' efsimplef2_RSMPI: legas=T')
  else
     valn    = 0d0
     do ia   = 1,natom
        ic    = iclass(ia)
        valn  = valn + z(ic)
        if (Is_IO_Root_RSMPI()) write(6,*) ' ia z(ic)=',ia, z(ic)
        do    l = 0,nl-1
           if (Is_IO_Root_RSMPI()) then
              write(6,*) &
                   ' l (konfig(l,ic)-l-1) 2*(2l+1)=', &
                   l,(konfig(l,ic)-l-1),( 2*l +1)*2
           endif
           valn  = valn - (konfig(l,ic)-l-1) *( 2*l +1)*2
        end do
     end do
     if (Is_IO_Root_RSMPI()) write(6,*) ' valn=',valn
  endif

  do is = 1,ispin
     do iq = 1,nqibz
        ekt(:,iq,is) = readeval(qibz(:,iq),is)
     enddo
  enddo

  if(abs(sum(wibz(1:nqibz))-2d0)>1d-10) then
     if(Is_IO_Root_RSMPI()) then
        write(6,*) 'sum (wibz)=', sum(wibz(1:nqibz))
     endif
     call RSMPI_Stop('efsimplef2: wibzsumerr')
  endif
  do is = 1,ispin
     do iq = 1,nqibz
        wgt(1:nband,iq,is) = wibz(iq)
        if(ispin==2) wgt(1:nband,iq,is) = wgt(1:nband,iq,is)/2d0
     enddo
  enddo

  ! ekt and wgt
  call dcopy ( nband*nqibz*ispin, ekt,1, ektx ,1)
  call dcopy ( nband*nqibz*ispin, wgt,1, wgtx ,1)
  call sortea( ektx,ieaord,nband*nqibz*ispin,isig)

  !      write(6,*)nband,nqibz,ispin,nband*nqibz*ispin
  !      do ik  = 1, nband*nqibz*ispin
  !        write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wgtx(ieaord(ik))
  !      enddo


  ! RS: only ioroot opens these files
  if (Is_IO_Root_RSMPI()) then
     if8301=ifile_handle()
     if8302=ifile_handle()
     open(if8301,file = "DOSACC.lda")
     open(if8302,file = "DOSACC2.lda")
  endif

  wwg = 0d0
  wwg2= 0d0
  eee2= -1d99
  ikini= 1
  ierr= 1
  nne = nband*nqibz*ispin
  do ik = 1, nne
     !---
     if(eee2 +1d-4< ektx(ieaord(ik)) .OR. ik==nne ) then
        ! degeneracy check
        if((ik/=1) .AND. Is_IO_Root_RSMPI()) then
           write(if8302,"(2i6,2d23.15)") ikini,ik-1,eee2,wwg2
        endif
        wwg2 = wgtx(ieaord(ik))
        eee2 = ektx(ieaord(ik))
        ikini =ik
     else
        wwg2= wwg2 + wgtx(ieaord(ik))
     endif
     !---
     wwg = wwg + wgtx(ieaord(ik))
     !        if(wwg<valn+2d0) write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wwg
     if((wwg<valn+2d0) .AND. (Is_IO_Root_RSMPI())) then
        write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wwg
     endif

     if (Is_IO_Root_RSMPI()) write(if8301,"(2i6,3d23.15)") &
          ik,ieaord(ik),ektx(ieaord(ik)),wwg,wgtx(ieaord(ik))

     if( wwg>valn-1d-8 .AND. ierr==1 ) then
        if (Is_IO_Root_RSMPI()) write(6,*)
        efini = .5*(ektx(ieaord(ik+1))+ ektx(ieaord(ik)))
        if(autoew) then
           if(ik<3) call RSMPI_Stop(' efsimplef2: ik<3')
           esmr  = ektx(ieaord(ik)) - ektx(ieaord(ik-1))
        endif
        ierr=0
     endif
     if( wwg > valn+1d0) ikx=ik
  enddo

  if(ierr==1) then
     call RSMPI_Stop(' efsimplef2: ierr=1 given nval is too large')
  endif


  nbnqnsp = nband*nqibz*ispin
  ! c gaussian
  if(GaussSmear) then
     valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp)) &
          ,efini,esmr,nbnqnsp)
  else
     valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp)) &
          ,efini,esmr,nbnqnsp)
  endif
  if (Is_IO_Root_RSMPI()) write(6,*) 'valx at efini=',efini,valx
  if(abs(valx-valn)<1d-8) then
     ef=efini
     goto 8891
  endif

  ew1= ektx(ieaord(1))-0.01d0
  ew2= ektx(ieaord(ikx))
  nbnqnsp =nband*nqibz*ispin
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

8891 if (Is_IO_Root_RSMPI()) then
     write(if8301,*) " ef=",ef
     close(if8301)
     write(if8302,*) " ef=",ef
     close(if8302)
  endif
  !      write(6,*)' esmr        =',esmr
  !      write(6,*)' determined ef =',ef
  !----------------------------------
  !      wwg = 0d0
  !      do ik  = 1, nband*nqibz*ispin
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
  if (Is_IO_Root_RSMPI()) then
     if(GaussSmear) then
        write(6,*) ' efsimplef2(gauss):end'
     else
        write(6,*) ' efsimplef2:end'
     endif
  endif
end subroutine efsimplef2a_RSMPI
!----------------------------------------------------------------------------
