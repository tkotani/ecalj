subroutine getwemax(lqall,wemax)
  !!> this routine is just in order to get |e_ip-ef| on real space integration
  !! too complicated ---> need to fix in future
  use m_readhbe,only: nband
  use m_read_bzdata,only: read_bzdata, &
       nqibz,qibz,ginv,qbz,nqbz,wibz
  !     & dq_,qbz,wbz,qibz,wibz,qbzw,
  !     & idtetf,ib1bz,idteti,
  !     & nstar,irk,nstbz
  use m_genallcf_v3,only: &
       nspin, konf,z,nl,natom,iclass,nclass,esmr,deltaw!,dw
  use m_keyvalue,only:getkeyvalue
  use m_readeigen,only: readeval !init* is alreaday called.
  use m_ReadEfermi,only: ef !ef is set at main routine
  use m_anf,only:laf,anfcond

  implicit none
  logical,intent(in):: lqall
  real(8),intent(out):: wemax
  integer:: ntq,i,nq,ifqpnt,ret,iaf,ipx1,itx1,ipx2,itx2
  integer,allocatable:: itq(:)
  real(8),allocatable:: q(:,:),eqt(:) !,q0p(:,:)
  logical ::  legas = .false., tetra,lqallxxx
  real(8):: omegav,omegac,eee,efnew,emaxv,eminc,ffac,we,valn
  integer:: ifief,ib,ip,iq,iqall,it,k,is,nspinmx
  !      logical :: tetra=.false.

  lqallxxx=.true.
  call anfcond()
  if( .NOT. lqall) then
     call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
     call readx   (ifqpnt,10)
     iaf=-1
     read (ifqpnt,*) iqall,iaf
     if (iqall == 1) then
        lqallxxx = .true.
     else
        lqallxxx = .false.
        call readx   (ifqpnt,100)
        read (ifqpnt,*) ntq
        allocate( itq(ntq) )
        read (ifqpnt,*) (itq(i),i=1,ntq)
        call readx   (ifqpnt,100)
        read (ifqpnt,*) nq
        allocate(q(3,nq))
        do       k = 1,nq
           read (ifqpnt,*) i,q(1,k),q(2,k),q(3,k)
           write(6,'(i3,3f13.6)') i,q(1,k),q(2,k),q(3,k)
        enddo
     endif
     close(ifqpnt) !Walter fixed this Jul2008
  endif
  if(lqallxxx) then
     ntq = nband
     allocate (itq(ntq))
     do i = 1, ntq
        itq(i) = i
     enddo
     nq  = nqibz
     allocate(q(3,nq))
     q = qibz(1:3,1:nqibz) !call dcopy   (3*nqibz,qibz,1,q,1)
  endif

  nspinmx = nspin
  if (laf .OR. iaf==1) nspinmx =1

  !      if(tetra) then
  !        open( ifief,file='EFERMI')
  !        read( ifief,*) ef
  !        close(ifief)
  !        write(6,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef
  !      else

  !! comment out followings because esmr is not for hx0fp0, hx0fp0.sc.
  !! sep2017takao
  !! This is a bug; however, we have no changes for previous results, because
  !! a lucky another accicent that scissors_x0() in switches.F called readefermi() (this again read 'EFERMI').
  !!
  !$$$      call readefermi() !readin ef
  !$$$      if(esmr/=0d0) then
  !$$$          call efsimplef2a(nspin,wibz,qibz,ginv,
  !$$$     i        nband,nqibz
  !$$$     i       ,konf,z,nl,natom,iclass,nclass
  !$$$     i       ,valn, legas, esmr,  !!! valn is input for legas=T, output otherwise.
  !$$$     i       qbz,nqbz ! index_qbz, n_index_qbz,
  !$$$     o       ,efnew)
  !$$$          ef = efnew
  !$$$cc- check total ele number -------
  !$$$c          ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz)
  !$$$c          write(6,*)' ef    =',ef
  !$$$c          write(6,*)' esmr=',esmr
  !$$$c          write(6,*)' valn  =',valn
  !$$$c          write(6,*)' ntot  =',ntot
  !$$$      endif



  !      endif
  !$$$c ---  q near zero
  !$$$c      write(6,*) 'reading QOP'
  !$$$      open (ifq0p,file='Q0P')
  !$$$      read (ifq0p,"(i5)") nq0it
  !$$$c      write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0it
  !$$$      allocate( q0i(1:3,1:nq0it) )
  !$$$      do i=1,nq0it
  !$$$        read (101,* ) wdummy,q0i(1:3,i)
  !$$$      enddo
  !$$$c      write(6,*) ' k in Q0P =', nq0it
  !$$$c      write(6,"(i3, 3f14.6)" )(i, q0i(1:3,i),i=1,nq0it)
  !$$$      close(ifq0p)
  !$$$c --- qibze(3,nqbze) qbze(3,nqibze)
  !$$$      nqbze  = nqbz *(1 + nq0it)
  !$$$      allocate(  qbze(3,nqbze) )
  !$$$      call dcopy(3*nqbz, qbz, 1, qbze,1)
  !$$$      do i = 1,nq0it
  !$$$        ini = nqbz*(1 + i -1)
  !$$$        do ix=1,nqbz
  !$$$          qbze (:,ini+ix) = q0i(:,i) + qbze(:,ix)
  !$$$        enddo
  !$$$      enddo


  !! for 1shot GW deltaw id for d\Sigma/d_omega
  allocate(eqt(nband))
  omegav =  0d0 !1d99 fixed oct.2003 takao
  omegac =  0d0 !-1d99
  do is = 1,nspinmx
     do ip = 1,nq
        eqt = readeval(q(1,ip),is)
        do it=1,ntq
           eee = eqt(itq(it)) - 2d0*deltaw  !2.d0*(-1d0-shtw)*deltaw
           !            write(6,*)' is ip it eee=',eee,eqt(itq(it))
           if(eee>=1d20-1d10) cycle !takao jun2009
           if( eee <ef .AND. eee< omegav )  then
              omegav = eee
              ipx1 = ip
              itx1 = it
           endif
           eee = eqt(itq(it)) + 2d0*deltaw  !+ 2.d0*(1d0-shtw)*deltaw
           !            write(6,*)' is ip it eee=',eee,eqt(itq(it))
           if( eee >ef .AND. eee> omegac )  then
              omegac = eee
              ipx2 = ip
              itx2 = it
           endif
        enddo
     enddo
  enddo
  !! for Gaussian smearing
  !      if(GaussSmear()) then
  ffac=10d0 !This is OK?
  !      else
  !        ffac=0.5d0
  !      endif
  emaxv =  0d0 !-1d99 fixed oct.2003 takao
  eminc =  0d0 !1d99
  do is = 1, nspinmx
     do iq = 1, nqbz
        eqt= readeval(qbz(1,iq),is)
        !          print *
        do ib=1,nband
           eee = eqt(ib)
           !            if(ib<6) write(6,*)' iq ib eee=',iq,ib,eee
           if( eee <ef+ffac*esmr .AND. eee>emaxv) emaxv = eee
           if( eee >ef-ffac*esmr .AND. eee<eminc) eminc = eee
        enddo
     enddo
  enddo
  deallocate(eqt)
  we  = max(abs(emaxv - ef), abs(omegav-ef),abs(omegac- ef) , abs(eminc-ef) )
  wemax= we+ffac*esmr
  !      nw  = idint (we/2d0/dw) + 3
  !      write(6,*)' --------------------------------'
  !      write(6,*)' emaxv= ',emaxv
  !      write(6,*)' eminc= ',eminc
  !      write(6,*)' omegav ip it=',omegav ,ipx1,itx1
  !      write(6,*)' omegac ip it=',omegac ,ipx2,itx2
  !      write(6,*)' we max for <ef', emaxv - omegav
  !      write(6,*)' we max for >ef', omegac- eminc
  !      write(6,("(' wemax=  ',f13.4)") wemax
  !      write(6,*)"write nw to NW"
  !      open(1101,file='NW')
  !      write(1101,*) nw
  !      close(1101)
  !      write(6, *) ' --- computational conditions --- '
  !      write(6,'("    deltaw  =",f13.6)') deltaw
  !      write(6,'("    esmr    =",f13.6)') esmr
  !      write(6,'("    niw nw dw   =",2i6,f13.6)') niw,nw,dw
end subroutine getwemax
