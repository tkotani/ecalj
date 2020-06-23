      subroutine getwemax(lqall,wemax)
!!> this routine is just in order to get |e_ip-ef| on real space integration
!! too complicated ---> need to fix in future
      use m_zmel,only:
     & nband 
      use m_read_bzdata,only: read_bzdata,
     & nqibz,qibz,ginv,qbz,nqbz,wibz
c     & dq_,qbz,wbz,qibz,wibz,qbzw,
c     & idtetf,ib1bz,idteti,
c     & nstar,irk,nstbz
      use m_genallcf_v3,only:
     i  nspin, konf,z,nl,natom,iclass,nclass,esmr,deltaw!,dw
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
      integer:: ifile_handle,ifief,ib,ip,iq,iqall,it,k,is,nspinmx
c      logical :: tetra=.false.

      lqallxxx=.true.
      call anfcond()
      if(.not.lqall) then
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
      if (laf.or.iaf==1) nspinmx =1

c      if(tetra) then
c        ifief=ifile_handle()
c        open( ifief,file='EFERMI')
c        read( ifief,*) ef
c        close(ifief)
c        write(6,"(a,f12.6)")' --- READIN ef from EFERMI. ef=',ef
c      else

!! comment out followings because esmr is not for hx0fp0, hx0fp0.sc.
!! sep2017takao
!! This is a bug; however, we have no changes for previous results, because
!! a lucky another accicent that scissors_x0() in switches.F called readefermi() (this again read 'EFERMI').
!! 
c$$$      call readefermi() !readin ef
c$$$      if(esmr/=0d0) then
c$$$          call efsimplef2a(nspin,wibz,qibz,ginv, 
c$$$     i        nband,nqibz
c$$$     i       ,konf,z,nl,natom,iclass,nclass
c$$$     i       ,valn, legas, esmr,  !!! valn is input for legas=T, output otherwise.
c$$$     i       qbz,nqbz ! index_qbz, n_index_qbz,
c$$$     o       ,efnew)
c$$$          ef = efnew
c$$$cc- check total ele number -------
c$$$c          ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz)
c$$$c          write(6,*)' ef    =',ef
c$$$c          write(6,*)' esmr=',esmr
c$$$c          write(6,*)' valn  =',valn
c$$$c          write(6,*)' ntot  =',ntot
c$$$      endif



c      endif
c$$$c ---  q near zero
c$$$c      write(6,*) 'reading QOP'
c$$$      open (ifq0p,file='Q0P')
c$$$      read (ifq0p,"(i5)") nq0it
c$$$c      write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0it
c$$$      allocate( q0i(1:3,1:nq0it) )
c$$$      do i=1,nq0it
c$$$        read (101,* ) wdummy,q0i(1:3,i)
c$$$      enddo
c$$$c      write(6,*) ' k in Q0P =', nq0it
c$$$c      write(6,"(i3, 3f14.6)" )(i, q0i(1:3,i),i=1,nq0it)
c$$$      close(ifq0p)
c$$$c --- qibze(3,nqbze) qbze(3,nqibze)
c$$$      nqbze  = nqbz *(1 + nq0it)
c$$$      allocate(  qbze(3,nqbze) )
c$$$      call dcopy(3*nqbz, qbz, 1, qbze,1)
c$$$      do i = 1,nq0it
c$$$        ini = nqbz*(1 + i -1)
c$$$        do ix=1,nqbz
c$$$          qbze (:,ini+ix) = q0i(:,i) + qbze(:,ix)
c$$$        enddo
c$$$      enddo


!! for 1shot GW deltaw id for d\Sigma/d_omega
      allocate(eqt(nband))
      omegav =  0d0 !1d99 fixed oct.2003 takao
      omegac =  0d0 !-1d99
      do is = 1,nspinmx
        do ip = 1,nq
          call readeval(q(1,ip),is,eqt)
          do it=1,ntq
            eee = eqt(itq(it)) - 2d0*deltaw  !2.d0*(-1d0-shtw)*deltaw
c            write(6,*)' is ip it eee=',eee,eqt(itq(it))
            if(eee>=1d20-1d10) cycle !takao jun2009
            if( eee <ef .and. eee< omegav )  then
              omegav = eee
              ipx1 = ip
              itx1 = it
            endif
            eee = eqt(itq(it)) + 2d0*deltaw  !+ 2.d0*(1d0-shtw)*deltaw
c            write(6,*)' is ip it eee=',eee,eqt(itq(it))
            if( eee >ef .and. eee> omegac )  then
              omegac = eee
              ipx2 = ip
              itx2 = it
            endif
          enddo
        enddo
      enddo
!! for Gaussian smearing
c      if(GaussSmear()) then
        ffac=10d0 !This is OK?
c      else
c        ffac=0.5d0
c      endif
      emaxv =  0d0 !-1d99 fixed oct.2003 takao
      eminc =  0d0 !1d99
      do is = 1, nspinmx
        do iq = 1, nqbz
          call readeval(qbz(1,iq),is,eqt)
c          print *
          do ib=1,nband
            eee = eqt(ib)
c            if(ib<6) write(6,*)' iq ib eee=',iq,ib,eee
            if( eee <ef+ffac*esmr .and. eee>emaxv) emaxv = eee
            if( eee >ef-ffac*esmr .and. eee<eminc) eminc = eee
          enddo
        enddo
      enddo
      deallocate(eqt)
      we  = max(abs(emaxv - ef), abs(omegav-ef),abs(omegac- ef) , abs(eminc-ef) )
      wemax= we+ffac*esmr
c      nw  = idint (we/2d0/dw) + 3
c      write(6,*)' --------------------------------'
c      write(6,*)' emaxv= ',emaxv
c      write(6,*)' eminc= ',eminc
c      write(6,*)' omegav ip it=',omegav ,ipx1,itx1
c      write(6,*)' omegac ip it=',omegac ,ipx2,itx2
c      write(6,*)' we max for <ef', emaxv - omegav
c      write(6,*)' we max for >ef', omegac- eminc
c      write(6,("(' wemax=  ',f13.4)") wemax
c      write(6,*)"write nw to NW"
c      open(1101,file='NW')
c      write(1101,*) nw
c      close(1101)
c      write(6, *) ' --- computational conditions --- '
c      write(6,'("    deltaw  =",f13.6)') deltaw
c      write(6,'("    esmr    =",f13.6)') esmr
c      write(6,'("    niw nw dw   =",2i6,f13.6)') niw,nw,dw
      end
