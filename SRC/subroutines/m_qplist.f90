!> Get q points list.
module m_qplist 
  use m_ftox
  use m_lgunit,only:stdo
  public :: m_qplist_init,m_qplist_qspdivider,qshortn
  integer,protected,public::   napwmxqp
  integer,allocatable,public:: igv2qp(:,:,:),igv2revqp(:,:,:,:),napwkqp(:)
  
  integer,parameter,private:: nsymlmax=100
  integer,protected,public:: nkp,ngpmx ,nqi,iqibzmax ,  iqini,iqend     !for current rank
  integer,protected,public:: nsyml=0, nqp_syml(nsymlmax),nqps_syml(nsymlmax),nqpe_syml(nsymlmax),nqp2n_syml(nsymlmax)
  integer,allocatable,protected,public:: ngplist(:),iprocq(:,:),ispp(:),ngvecp(:,:,:), kpproc(:)
  integer,protected,allocatable,public:: iqproc(:),isproc(:)
  integer,protected,public:: niqisp
  real(8),protected,public:: dqsyml(nsymlmax),etolv,etolc
  real(8),allocatable,protected,public:: qplist(:,:), xdatt(:)
  character*20,protected,public ::labeli(nsymlmax),labele(nsymlmax)
  private
  real(8),allocatable:: qplistss(:,:)
  real(8)::tolq=1d-12
contains
  subroutine m_qplist_init(plbnd,llmfgw)
    use m_lmfinit,only: nspx,pwemax,alat,pwmode
    use m_mkqp,only: rv_p_oqp,rv_a_owtkp
    use m_MPItk,only: procid,master,master_mpi,mlog
    use m_mkqp,only: nkabc=> bz_nabc,bz_nkp
    use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
    use m_ext,only: sname
    intent(in)::             plbnd,llmfgw
    integer:: nqp2_syml(nsymlmax),nqp2s_syml(nsymlmax),nqp2e_syml(nsymlmax)
    logical:: masslineon(nsymlmax),llmfgw
    integer:: plbnd,ifqplist,nkk1,nkk2,nkk3,ik1,ik2,ik3,iq,ifi,ifisyml,isyml,i,ierr
    logical:: cmdopt0,fullmesh,PROCARon,fsmode
    integer:: ikp,ifbnd,nsymln,onesp=0,jobgw
    real(8):: qps_syml(3,nsymlmax), qpe_syml(3,nsymlmax),rq
    real(8)::  totxdatt
    integer:: fxsts,ii,ifiproc
    character(512)::schar
    character(3):: charnum3
    character(50)::infoq
    integer:: nqnum, nqbz,iqq,isp,ifiqg,irr,ngp ,imx,nqibz,iqibz
    real(8):: QpGcut_psi,q(3)
    real(8)::qq(3)
    call tcn('m_qplist_init')
    if(master_mpi) write(stdo,ftox)'m_qplistinit:start'
    nkk1=nkabc(1)
    nkk2=nkabc(2)
    nkk3=nkabc(3)
    fullmesh = cmdopt0('--fullmesh').or. cmdopt0('--fermisurface') !full mesh stop just after do 2010 iq loop.
    fsmode   = cmdopt0('--fermisurface') !full mesh stop just after do 2010 iq loop.
    PROCARon = cmdopt0('--mkprocar')
    if(allocated(qplist)) deallocate(qplist)
    plbndmode: if( .NOT. master_mpi) then
       continue
    elseif(llmfgw) then ! GW driver mode !Read QGpsi Get ngplist,ngvecp in addition to qplist.
       open(newunit=ifiqg,file='QGpsi',form='unformatted',status='old')
       read(ifiqg) nqnum, ngpmx ,QpGcut_psi,nqbz ,nqi ,imx,nqibz
       allocate(qplist(3,nqi), ngplist(nqi), ngvecp(3,ngpmx,nqi))
       ngvecp=0
       open(newunit=ifqplist,file='QPLIST.lmfgw.chk')
       nkp=0
       iqibz=0
       do iq=1,nqnum
          read(ifiqg) q,ngp,irr
          if(irr==0) then
             read(ifiqg)
             cycle
          endif
          nkp = nkp+1 !we count only irr=1 case, that is, nkp=nqirr is returned.
          qplist(:,nkp) = q
          ngplist(nkp)  = ngp
          iqibz=iqibz+1
          if(iqibz==nqibz) iqibzmax=nkp
          read(ifiqg) ngvecp(1:3,1:ngp,nkp)
          write(ifqplist,"(i5,3f23.15)")nkp, qplist(:,nkp) !,ngplist(nkp)
       enddo
       close(ifiqg)
       close(ifqplist)
    elseif(plbnd==0) then
       nkp=bz_nkp
       allocate(qplist(3,nkp))
       qplist= rv_p_oqp !call dcopy(3*nkp, rv_p_oqp,1,qplist,1)
       !open(newunit=ifqplist,file='QPLIST.IBZ')
       !do iq=1,nkp
       !   write(ifqplist,"(i5,3f23.15,3x,3f23.15)")iq, qplist(:,iq), rv_a_owtkp(iq)/2d0
       !enddo
       !close(ifqplist)
       !! plbnd/=1 band plot mode
    else
       if(cmdopt0('--onesp') .AND. nspx==1) onesp = 1
       if(fullmesh) then
          nkp = nkk1*nkk2*nkk3
          allocate(qplist(3,nkp))
          if(fsmode) then !!! Fermi surface version for xcrysden
             iq=0
             do ik1=1,nkk1    !ordering is differnt from procaron case...
                do ik2=1,nkk2
                   do ik3=1,nkk3
                      iq=iq+1
                      qplist(:,iq)   =  qlat(:,1)*dble(ik1-1)/(nkk1-1) &
                           +   qlat(:,2)*dble(ik2-1)/(nkk2-1) &
                           +   qlat(:,3)*dble(ik3-1)/(nkk3-1)
                   enddo
                enddo
             enddo
          endif
          !! pdos (--mkprocar and --fullmesh)
          if(procaron) then !to fit to the tetirr.F
             iq=0
             do ik3=1,nkk3
                do ik2=1,nkk2
                   do ik1=1,nkk1
                      iq=iq+1
                      qplist(:,iq) =  qlat(:,1)*dble(ik1-1)/nkk1 &
                           +   qlat(:,2)*dble(ik2-1)/nkk2 &
                           +   qlat(:,3)*dble(ik3-1)/nkk3
                   enddo
                enddo
             enddo
          endif
          !! syml direct read for plbnd mode. See "call writeband" below. feb2015
       else
          !            readeferm=.false.
          !            if(master_mpi) write(stdo,*)' --- Readin efermi.lmf --- '
          !            open(newunit=ifi,file='efermi.lmf',status='old',err=1012)
          !            read(ifi,*,err=1012) eferm
          !            readeferm=.true.
          !            close(ifi)
          !            goto 1013
          ! 1012       continue
          !            call rx('No efermi.lmf!: Copy it, or run lmf-MPIK (sc mode) to get efermi.lmf.')
          ! 1013       continue

          !! --- example of syml file ---
          ! ndiv qleft(1:3) qright(1:3) llabel rlabel  ndiv2 ninit2 nend2 etolv(Ry) etolc(Ry)
          ! 5  0 0 0   .5 .5  .5        Gamma  L       1025  1  16     0.1      0.01
          ! 5  0 0 0    1.  0  0        Gamma  X
          ! 5  0 0 0   .75 .75 0        Gamma  K       1025  8  38     0.1      0.01
          !! As this shows you can add or not add line after ndiv2. These are for highly resolved calculations.
          write(stdo,*)' --- Readin syml file --- '
          open(newunit=ifisyml,status='old',file='syml.'//trim(sname))
          nsyml=0
          nsymln=0
          nqp2_syml=0
          do
             if(nsyml+1>nsymlmax) call rx('bndfp: Enlarge nsymlmax')
             read(ifisyml,"(a)",end=1015) schar
             if(len(trim(schar))==0 .OR. schar(1:1)=='#' .OR. schar(1:1)=='!' .OR. schar(1:1)=='%') cycle !comment line
             read(schar,*,err=1014,end=1014) &
                  nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1), &
                  labeli(nsyml+1),labele(nsyml+1), nqp2_syml(nsyml+1),nqp2s_syml(nsyml+1),nqp2e_syml(nsyml+1),etolv,etolc
             masslineon(nsyml+1)=.true.
             nqp2n_syml(nsyml+1)= nqp2e_syml(nsyml+1)-nqp2s_syml(nsyml+1)+1
             nqps_syml(nsyml+1)=1
             nqpe_syml(nsyml+1)=nqp_syml(nsyml+1)
             write(stdo,"(' ',i4,3f9.4,' ',3f9.4,' ',a,' ',a,'  Massl:div,init,end=',3i5)") &
                  nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1), &
                  trim(labeli(nsyml+1)),trim(labele(nsyml+1)), &
                  nqp2_syml(nsyml+1),nqp2s_syml(nsyml+1),nqp2e_syml(nsyml+1)
             goto 1025

1014         continue
             read(schar,*,err=1015,end=1015) &
                  nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1), &
                  labeli(nsyml+1),labele(nsyml+1)
             masslineon(nsyml+1)=.false.
             nqps_syml(nsyml+1)=1
             nqpe_syml(nsyml+1)=nqp_syml(nsyml+1)
             write(stdo,"(' ',i4,3f9.4,' ',3f9.4,' ',a,' ',a)") &
                  nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1), &
                  trim(labeli(nsyml+1)),trim(labele(nsyml+1))
             nsymln = nsymln+1
             nqp2n_syml(nsyml+1)= 0
1025         continue
             if(nqp_syml(nsyml+1)==0) exit
             nsyml = nsyml + 1
          enddo
1015      continue
          close(ifisyml)
          nkp = sum(nqp_syml(1:nsyml)+nqp2n_syml(1:nsyml))
          allocate(qplist(3,nkp))
          if(allocated(xdatt)) deallocate(xdatt)
          allocate(xdatt(nkp))
          totxdatt=0d0
          ikp=0
          do isyml=1,nsyml
             dqsyml(isyml) = dsqrt(sum((qpe_syml(1:3,isyml) -qps_syml(1:3,isyml))**2))
             do i=1,nqp_syml(isyml)+nqp2n_syml(isyml)
                ikp= ikp+1
                if(i<=nqp_syml(isyml)) then
                   rq = dble(i-1)/(nqp_syml(isyml)-1)
                else
                   ii= i-nqp_syml(isyml)-1
                   rq = dble(nqp2s_syml(isyml)-1+ii)/(nqp2_syml(isyml)-1)
                endif
                qplist(:,ikp)= (1d0-rq)*qps_syml(1:3,isyml) +rq*qpe_syml(1:3,isyml)
                xdatt(ikp) = totxdatt + dqsyml(isyml)*rq
             enddo
             totxdatt = totxdatt + dqsyml(isyml)
          enddo
          write(stdo,"('nsyml nkp=',3i5)") nsyml,nkp
       endif
       if (nkp <= 0) call rx('bndfp: nkp<=0') ! quit if nkp==0
    endif plbndmode
    
    if (master_mpi .AND. nsyml>0) then !plbnd mode
       open(newunit=ifqplist,file='QPLIST')
       print *,'-------- qplist --------',nsyml
       iq=0
       do isyml=1,nsyml
          do i=1,nqp_syml(isyml) + nqp2n_syml(isyml)
             iq=iq+1
             infoq=''
             if(i==1) infoq=' <-- isyml= '//charnum3(isyml)
             if(i==nqp_syml(isyml)+1) infoq=' <-- isyml Mass= '//charnum3(isyml)
             write(stdo,"(i5,3f8.3,' ',a)")iq,qplist(:,iq),trim(infoq)
             write(ifqplist,"(i5,3f23.15,x,f12.6,' ',a)")iq,qplist(:,iq),xdatt(iq),trim(infoq)
          enddo
          write(ifqplist,*)
       enddo
       close(ifqplist)
    endif
    ! broadcase nkp and qplist
    call mpibc1_int( nkp,1,  'qplit_nkp')
    call mpibc1_int( onesp,1,'qplist_onesp')
    if( .NOT. master_mpi) allocate(qplist(3,nkp))
    call mpibc1_real(qplist, 3*nkp , 'qplist_qp'  )
    if(llmfgw) then
       if( .NOT. master_mpi) allocate(ngplist(nkp))
       call mpibc1_int(ngpmx,1,      'qplist_ngpmx')
       if( .NOT. master_mpi) allocate(ngvecp(3,ngpmx,nkp))
       call mpibc1_int(ngplist, nkp, 'qplist_ngp' )
       call mpibc1_int( ngvecp,  3*ngpmx*nkp,  'qplist_ngvecp' )
       call mpibc1_int( iqibzmax,1, 'qplist_iqibzmax')
    endif
    if(procaron .AND. nsyml==0 .AND. master_mpi) then !xdatt is dummy
       if(allocated(xdatt)) deallocate(xdatt)
       allocate(xdatt(nkp))
       xdatt=0d0
    endif
    if(PROCARon) then
       if( .NOT. master_mpi) allocate(xdatt(nkp))
       call mpibc1( xdatt, nkp , 4 , .false. , 'bndfp' , 'xdatt'  )
    endif
    allocate(qplistss(3,nkp))
    qplistshortened: block  !2023-4-27
      use m_shortn3_qlat,only: shortn3_qlat,nout,nlatout
      real(8):: qfrac(3),qs(3)
      do ikp=1,nkp !
         qfrac = matmul(transpose(plat),qplist(:,ikp)) !qfrac: qplist in fractional coordinate
         call shortn3_qlat(qfrac)
         qplistss(:,ikp)=matmul(qlat, qfrac+nlatout(1:3,1)) !shortened qplist
      enddo
    endblock qplistshortened
    !
    Gindexqplist: block
      integer:: nnn(3),ig,ikt,i
      real(8):: qqq(3),dum,pwgmax
      PMTmodeOnly: if(mod(pwmode,10)/=0 .and. pwemax>1d-8) then
         pwgmax = pwemax**.5
         allocate(napwkqp(nkp))
         do ikt=1,nkp
            qqq = merge(qplist(:,ikt), 0d0, mod(pwmode/10,10)==1)
            call getgv2(alat,plat,qlat,qqq, pwgmax,1, napwkqp(ikt),dum)
         enddo
         napwmxqp=maxval(napwkqp)
         allocate(igv2qp(3,napwmxqp,nkp))
         do ikt = 1,nkp
            qqq = merge(qplist(:,ikt),0d0,mod(pwmode/10,10) == 1)
            call getgv2(alat,plat,qlat,qqq, pwgmax,2, napwkqp(ikt),igv2qp(:,:,ikt)) 
         enddo
         imx = maxval([(maxval(abs(igv2qp(1:3,1:napwkqp(ikt),ikt))),ikt=1,nkp)])
         allocate( igv2revqp(-imx:imx,-imx:imx,-imx:imx,nkp),source=999999 ) !Reverse table of igv2 
         do ikt = 1,nkp
            do ig  = 1,napwkqp(ikt)
               nnn = igv2qp(1:3, ig, ikt)
               igv2revqp( nnn(1), nnn(2),nnn(3), ikt) = ig
            enddo
         enddo
      endif PMTmodeOnly
!      write(6,*)'napwkqp=',napwkqp
!      stop 'vvvvvvvvvvvvvvvvvvvvvvv'
    endblock Gindexqplist
    call tcx('m_qplist_init')
  end subroutine m_qplist_init
  function qshortn(q) result(qs) !shortest q vectror. module of qlat 2023-4-27
    intent(in):: q
    integer:: ikp,iq
    real(8):: q(3),qs(3)
    if(sum(abs(q))<tolq) then
       qs=0d0
    else
       iq = findloc([(sum(abs(qplist(:,ikp)-q))<tolq,ikp=1,nkp)],value=.true.,dim=1)
       if(iq==0) call rx('qshortn: can not find iq')
       qs = qplistss(:,iq) !write(aaa,ftox)'q=',ftof(q,3),'qshortn=',ftof(qs,3),'qdiff=',ftof(q-qs,3),iq
    endif   
  end function qshortn
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine m_qplist_qspdivider() ! MPIK k point divider. From iqini to iqend for each processor. ! Set iqini,iqend ispx for each rank (procid)
    use m_MPItk,only: procid,master,master_mpi,mlog, numprocs=>nsize
    use m_lmfinit,only: nsp,nspc,afsym
    use m_ext,only: sname
    use m_dstrbp,only: dstrbp
    implicit none
    integer:: iqq,isp,ispx,icount,iqs,ncount,iqsi,iqse,iprint,nspxx,idat,i,nsize
    logical:: cmdopt0
    call tcn('m_qplist_qpsdivider')
    nspxx=nsp/nspc !nspc is 2 for spin-coupled case
!    if(cmdopt0('--afsym')) nspxx=1
    if((.not.cmdopt0('--jobgw')).and.(.not.cmdopt0('--writeham')).and.afsym) nspxx=1
    allocate(kpproc(0:numprocs))
    call dstrbp(nkp*nspxx, numprocs,1,kpproc(0))
    ! i=1,nkp*nspxx is divided into [kpproc(procid),kpporc(procid+1)-1] for each procid.
    iqsi = kpproc(procid)
    iqse = kpproc(procid+1)-1
    niqisp=iqse-iqsi+1
    ! (iq,isp) is ordered as (1,1),(1,2),(2,1),(2,2),(3,1),(3,2),(4,1),(4,2).....  
    if(niqisp>0) then
       allocate(iqproc(niqisp),isproc(niqisp))
       do i=iqsi,iqse
          iqproc(i-iqsi+1) =   (i-1)/nspxx +1
          isproc(i-iqsi+1) = mod(i-1,nspxx)+1
       enddo
       iqini=iqproc(1)
       iqend=iqproc(niqisp)
    else
       iqini=0
       iqend=-1
    endif
!    do idat=1,niqisp
!       write(stdo,ftox)'qspdivider: procid=',procid,' idat iq isp=',idat,iqproc(idat),isproc(idat)
!    enddo
    call tcx('m_qplist_qpsdivider')
  end subroutine m_qplist_qspdivider
end module m_qplist
! !mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
! module m_read_lmfgw_kdivider
!   use m_lgunit,only:stdo
!   character*256,allocatable:: extp(:)
!   integer,allocatable,protected:: iprocq(:,:)
! contains
!   subroutine readlmfgw_kdivider()
!     integer:: ifiproc,nqixx, nspxx, nrank,ispiqq,iqqxx,ispxx,ixxx,procid,iqq,isp
!     character*256:: ext,extn
!     !! === readin lmfgw_kdivider, and get extensions
!     open(newunit=ifiproc,file= 'lmfgw_kdivider')
!     read(ifiproc,*) ext  !extension is read here
!     read(ifiproc,*) nqixx, nspxx, nrank
!     allocate(iprocq(nqixx,nspxx))
!     do isp=1,nspxx
!        do iqq=1,nqixx
!           read(ifiproc,*) iqqxx, ispxx, ixxx
!           if(iqqxx/=iqq) stop 'iqqxx/=iqq'
!           if(ispxx/=isp) stop 'ispxx/=isp'
!           iprocq(iqq,isp) = ixxx
!           write(stdo,"('iqq isp irank=',i8,i2,i6)") iqq,isp, iprocq(iqq,isp)
!        enddo
!     enddo
!     close(ifiproc)
!     !! for mpi files.
!     allocate(extp(0:nrank-1))
!     extp(0) = trim(ext)
!     write(stdo,"('  0 ext= ',a,a)") trim(extp(0)),' ----------'
!     do procid=1,nrank-1
!        write(extn,"(i10)") procid
!        extp(procid)=trim(adjustl(ext))//'_'//trim(adjustl(extn))
!        write(stdo,"(i3,' ext= ',a,a)") procid,trim(extp(procid)),' ----------'
!     enddo
!   end subroutine readlmfgw_kdivider
! end module m_read_lmfgw_kdivider
! !mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
module m_readqplist 
  use m_lgunit,only:stdo
  integer,protected:: ndat
  real(8),allocatable,protected:: xdat(:),qplistsy(:,:)
  real(8),protected:: eferm !<-- temporary use. just as a memo.
contains
  subroutine readqplistsy()
    implicit none
    integer:: ifqplistsy,nnn,ix
    open(newunit=ifqplistsy,file='qplist.dat')
    nnn=10000
    if(allocated(xdat)) deallocate(xdat)
    if(allocated(qplistsy)) deallocate(qplistsy)
    allocate(xdat(nnn),qplistsy(1:3,nnn))
    read(ifqplistsy,*) eferm
    ix=0
    do
       ix=ix+1
       read(ifqplistsy,*,end=1011) xdat(ix),qplistsy(1:3,ix)
       !        write(stdo,'(" qplist ix xdata q=",i5,f9.4,x,3f9.4)')ix,xdat(ix),qplistsy(:,ix)
    enddo
1011 continue
    ndat=ix-1
    close(ifqplistsy)
  end subroutine readqplistsy
end module m_readqplist
