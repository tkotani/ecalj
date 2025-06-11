!>output routines such as writeband used in bndfp
module m_writeband
  use m_MPItk,only: comm
  use m_ftox
  real(8),external:: rydberg
  public writeband,writefs,writepdos,writedossawada
  private
contains
  subroutine writeband(evlall,eferm,vesav,evtop,ecbot,spinweightsoc) !write band file. bnd* and bandplot.isp*.glt
    use m_lgunit,only:stdo
    use m_lmfinit,only:nsp,alat=>lat_alat,lso,nspx
    use m_qplist,only: nkp,nsyml,xdatt,nqp_syml,nqp2n_syml,qplist,labeli,labele,nqps_syml,nqpe_syml,dqsyml,etolv,etolc
    use m_bandcal,only:nevls
    use m_ext,only: sname,dirname
    implicit none
    real(8),intent(in):: eferm,vesav,evtop,ecbot ! evtop is max of n-th band. !evbot is bottom of bands upper than n+1
    integer:: ifbndo,ikp,isyml,jsp
    character*300::filenm(2),bchar
    character*3::charnum3
    integer::ifbndsp(2),iprint,iq,i,ifglt,ifbndsp_nearef(2),ifmass(2),ifmglt,ifmglt2,ndhamx
    character*300::aaa,addx
    real(8):: emin=-20d0,emax=20d0 !eV for default plotting.
    integer:: ikps,ne,ifi,ix,nee,ibb,ilt,imin,imax
    real(8),allocatable::diffeb(:),diff2eb(:)
    real(8)::polinta,tpiba,eee,eqq,qqq,dEdkatef
    character*100::fname,massfile,fname2,sss=''
    integer,allocatable::ikpoff(:)
    real(8),allocatable::disoff(:)
    character(100),allocatable:: fnameb(:,:),fnamem(:,:,:)
    logical:: initii
    logical:: semiconband,metalband
    integer:: idat,ifglts(2),iqplist,isp,nx(3)
    character(100)::acrossef
    character(13)::massd,mass2d,labelp
    logical:: scd,cmdopt0
    real(8)::kef
    real(8),allocatable:: kabs(:)
    real(8)::  evlall(:,:,:),basel
    logical:: eszero
    real(8):: spinweightsoc(:,:,:)
    nx=shape(evlall)
    ndhamx=nx(1)
    scd= evtop <ecbot
    tpiba = 2d0*4d0*datan(1d0)/alat
    !! ikpoffset
    allocate(ikpoff(nsyml+1),disoff(nsyml+1))
    ikpoff(1)=0
    disoff(1)=0
    do isyml = 2,nsyml+1
       ikpoff(isyml) = ikpoff(isyml-1)+ nqp_syml(isyml-1)+nqp2n_syml(isyml-1)
       disoff(isyml) = disoff(isyml-1)+dqsyml(isyml-1)
       !write(stdo,"('ikpoff=',2i5)") isyml,ikpoff(isyml)
    enddo
    !! write bandplot.glt for gnuplot
    eszero = cmdopt0('--eszero')
    basel= merge(vesav,eferm,eszero)
    allocate(fnameb(nsyml,nspx))
    do jsp = 1, nspx
       fname='bandplot.isp'//char(48+jsp)
       open(newunit=ifglt, file=trim(fname)//'.glt')
       ifglts(jsp)=ifglt
       write(ifglt,'(a)')'set terminal postscript enhanced color eps'
       write(ifglt,'(a)')'set output "'//trim(fname)//'.band.eps"'
       write(ifglt,'(a)')'set xzeroaxis'
       write(ifglt,'(a)')'set grid'
       if(eszero) write(ifglt,'(a)')'set ylabel "Energy-Eesav(eV)"'
       if(.not.eszero) write(ifglt,'(a)')'set ylabel "Energy-Efermi(eV)"'
       write(ifglt,'(a)')'# This is given written in subroutine writeband in lm7K/fp/bndfp.F'
       if(nspx==1) addx=''
       if(nspx==2) addx=' isp='//char(48+jsp)
       write(ifglt,'(a)')'set title "Band '//trim(sname)//trim(addx)//' at '//trim(dirname)//'"'
       write(ifglt,'(a,F12.5,a,F12.5,a)') 'set yrange [',emin,':',emax,']'
       write(ifglt,'(a,F12.5,a)') 'set xrange [0.0:',disoff(nsyml+1),']'
       !!
       write(ifglt,'(a,$)') 'set xtics ('
       do isyml = 1,nsyml   !symmetry line start
          labelp=labeli(isyml)
          if(isyml>1) then
             if(trim(labele(isyml-1))/=labeli(isyml)) labelp=trim(labele(isyml-1))//'|'//trim(labeli(isyml))
          endif   
          write(ifglt,'(a,F15.10,",\")')"'"//trim(labelp)//"'",disoff(isyml)
       enddo
       write(ifglt,'(a,F15.10,")")') "'"//trim(labele(nsyml))//"'",disoff(nsyml+1)
!       write(ifglt,'(a,$)') 'set x2tics ('
!       do isyml = 2,nsyml   !symmetry line end
!          write(ifglt,'(a,F15.10,",\")')"'"//trim(labele(isyml-1))//"'",disoff(isyml)
!       enddo
!       write(ifglt,'(a,F15.10,")")')"'"//trim(labele(nsyml))//"'",disoff(nsyml+1)
       write(ifglt,'(a,F12.5)') 'tpia=2*3.1415926/',alat
       write(ifglt,"('plot \')")
       do isyml = 1,nsyml
          fnameb(isyml,jsp)='bnd'//charnum3(isyml)//'.'//'spin'//char(48+jsp)
          if(isyml/=1) write(ifglt,'(",\")')
          write(ifglt,"(a,$)")  '"'//trim(fnameb(isyml,jsp))//'" u ($2):($3) lt 1 pt 1 w lp'
       enddo
    enddo

    !! Write bnd**.isp*.glt, open qplist.dat
    open(newunit=iqplist,file='qplist.dat')
    do 4111 isyml = 1,nsyml
       ne = nqp_syml(isyml) !+nqp2n_syml(isyml)
       do iq = 1,ne
          ikp = ikpoff(isyml)+iq
          do jsp = 1, nspx
             if (iprint() > 20) then
                write(stdo,'(" bndfp: kpt",i5," of",i5, " k jsp=",3f9.5,i2," nev=",i5)') &
                     ikp,nkp,qplist(:,ikp),jsp,nevls(ikp,jsp)
             endif
             if(iprint() >=35) then
                write(stdo,"(9f8.4)") (evlall(i,jsp,ikp), i=1,nevls(ikp,jsp))
             endif
          enddo
       enddo
       do 4113 jsp = 1, nspx   ! ispx index.
          open(newunit=ifbndsp(jsp),file=fnameb(isyml,jsp))
          if(lso==1) sss='                '//'spinup     spindn'
          write(ifbndsp(jsp),"('#',i5,'         ',f10.5,2x,'QPE(ev)',9x,'1st-deri',14x,'qvec',a)") nkp,eferm,trim(sss)
          if(.not.eszero) write(ifbndsp(jsp),ftox)'#base=eferm  ', ftof(eferm*rydberg(),8),ftof(vesav*rydberg(),8),' ! efermi Vesav (eV)'
          if(     eszero) write(ifbndsp(jsp),ftox)'#base=estatic', ftof(eferm*rydberg(),8),ftof(vesav*rydberg(),8),' ! efermi Vesav (eV)'
          ibb=0
          do 5113 i=1,minval(nevls(:,:)) !band index
             !! take derivatives
             allocate(diffeb(ne))
             diffeb(1)  = 0d0
             diffeb(ne) = 0d0
             do iq = 2,ne-1
                ikp = ikpoff(isyml)+iq
                diffeb(iq)= (evlall(i,jsp,ikp+1)-evlall(i,jsp,ikp-1))/(xdatt(ikp+1)-xdatt(ikp-1))/tpiba
             enddo
             do iq = 1,ne
                ikp = ikpoff(isyml)+iq
                if(lso==1) write(sss,"(2f11.6)") spinweightsoc(i,1:2,ikp)
                write(ifbndsp(jsp),"(i5,(1x,f12.8), 2(1x,f16.10),x,3f9.5,a)") &
                     i,xdatt(ikp),(evlall(i,jsp,ikp)-basel)*rydberg(),diffeb(iq), qplist(:,ikp), trim(sss)
                !! write qplist.dat (xaxis, q, efermi) used for band plot
                if(i==1 .AND. jsp==1) then
                   if(iq==1 .AND. isyml==1) write(iqplist,"(f16.10,' ! ef or estaticav')") basel
                   ikp = ikpoff(isyml)+iq
                   write(iqplist,"(f16.10,x,3f9.5,' !x q')") &
                        xdatt(ikp), qplist(:,ikp)
                endif
             enddo
             deallocate(diffeb)
             write(ifbndsp(jsp),*)
5113      enddo
          close(ifbndsp(jsp))
4113   enddo
4111 enddo
    close(iqplist)

    !! Write Band*Syml*Spin*.dat : These are between [VBM(left)-etolv(Ry), CBM(left) + etolc(Ry)]
    allocate(fnamem(minval(nevls(:,:)),nsyml,nspx))
    do jsp = 1, nspx
       do isyml = 1,nsyml
          ibb=0
          do i=1,minval(nevls(:,:))
             ibb = ibb+1
             fnamem(ibb,isyml,jsp)= &
                  'Band'//charnum3(ibb)//'Syml'//charnum3(isyml)//'Spin'//char(48+jsp)//'.dat'
          enddo
       enddo
    enddo
    do 2111 isyml = 1,nsyml
       ne = nqp2n_syml(isyml)
       if(ne==0) cycle
       allocate(diffeb(ne),diff2eb(ne))
       !! diffeb and diff2eb are 1st and 2nd derivatives,
       !! dEi(k)/dk d^2Ei(k)/d^2k along a symmetry line (isyml index).
       !! But be careful to determine effective mass, because it is not analytic at Gamma point.
       !! E.g, see Kittel's book.
       do 2112 jsp = 1, nspx   ! ispx index.
          ibb=0
          do 2113 i = 1,minval(nevls(:,:)) !band index
             !! take derivatives
             diffeb(1)  = 0d0
             diffeb(ne) = 0d0
             diff2eb(1) = 0d0
             diff2eb(ne)= 0d0
             do iq = 2,ne-1
                ikp = ikpoff(isyml)+ iq + nqp_syml(isyml)
                diffeb(iq)= (evlall(i,jsp,ikp+1)-evlall(i,jsp,ikp-1))/(xdatt(ikp+1)-xdatt(ikp-1))/tpiba
             enddo
             !! NOTE: effective mass. 2d0 is because emass in Ry unit is 1/2
             ikps= ikpoff(isyml)+1 + nqp_syml(isyml)
             eee = evlall(i,jsp,ikps)
             !            print *,'jjjjjjjj',jsp,isyml,ikps,i,eee
             semiconband = scd.and.evtop-etolv<eee .and. eee<ecbot+etolv
             ! heck i-th band is neare VCM and CBM
             metalband = ichangesign(evlall(i,jsp,ikps:ikps+ne-1)-basel,ne) >-1
             ! or metal crosspoint point across basel
             ! heck i-th band is near VCM and CBM
             if(semiconband .OR. metalband) then
                ibb = ibb+1
                if(metalband) then
                   imin=max(2,idat-1)  !diffeb is given from 2 to ne-1
                   imax=min(ne-1,idat+2)
                   allocate( kabs(imin:imax))
                   kabs(imin:imax) = (xdatt(ikps+imin-1:ikps+imax-1)-disoff(isyml))*tpiba ! tpiba is 2pi/alat. |k|. left-end is zero.
                   kef = polinta(0d0, evlall(i,jsp,ikps+imin-1:ikps+imax-1)-basel, kabs, imax-imin+1) !we have k at Ef |k|=qef
                   dEdkatef = polinta(kef, kabs, diffeb(imin:imax), imax-imin+1)
                   deallocate(kabs)
                   !                 write(stdo,"('Effective Mass for metal: ',a,
                   !     &            ' iq isyml |k|_ef/(2pi/alat) dEdk 2*k/dEdk=', 2i4,' ',f8.3,' ',f8.3,' ',f8.3)")
                   !     &            trim(fnamem(ibb,isyml,jsp)),idat,isyml, kef/tpiba, dEdkatef, 2d0*kef/dEdkatef
                endif
                !               print *,'ffff ', ibb,isyml,jsp,trim(fnamem(ibb,isyml,jsp))
                open(newunit=ifmass(jsp),file=trim(fnamem(ibb,isyml,jsp)))
                write(ifmass(jsp),"(a)") '# mass1 is for metal, mass2 is for semiconductor'
                write(ifmass(jsp),"(a)") '# isyml,ib,iq,isp, |k|/(2pi*alat), QPE-EF, QPE-QPE(start),' &
                     //' mass2=2*(2*(QPE-QPE(start))/|k|**2), mass1=2*|k|/(dE/dk)'
                do ix=1,ne
                   acrossef=''
                   if(ix==idat)   then
                      write(bchar,"(f13.8,' ',f8.3)") kef/tpiba,2d0*kef/dEdkatef
                      acrossef=' <--- Ef. kef/(2pi/alat), 2*|k|/(dE/dk)_kef= '//trim(bchar)
                   endif
                   if(ix==idat+1) acrossef=' <--- Ef.'
                   qqq= (xdatt(ikps+ix-1)-disoff(isyml))*tpiba ! tpiba is 2pi/alat. |k|. left-end is zero.
                   eqq= (evlall(i,jsp,ikps+ix-1)-evlall(i,jsp,ikps))
                   ! in rydberg !eqq is relative to the eval at ikps(left end point)
                   if(ix==1 .OR. ix==ne) then
                      massd= '   --------  '
                      mass2d='   --------  '
                   else
                      write(massd,"(f13.8)")  2d0*qqq/diffeb(ix)
                      write(mass2d,"(f13.8)") 2d0/(2d0*eqq/qqq**2)
                   endif
                   write(ifmass(jsp),"(3i4,i2,' ',f13.8,' ',f13.8,' ',f13.8,' ',a,' ',a,a)") &
                        isyml,i,ix,jsp,(xdatt(ikps+ix-1)-disoff(isyml))*tpiba, &
                        (evlall(i,jsp,ikps+ix-1)-basel)*rydberg(), &
                        (evlall(i,jsp,ikps+ix-1)-evlall(i,jsp,ikps))*rydberg(), &
                        mass2d,  massd,  trim(acrossef)
                enddo
                close(ifmass(jsp))
                write(aaa,"(a,f13.5,a)") '" u ($5/tpia+',disoff(isyml),'):($6) lt'
                write(bchar,"(a,i2,a,i2,a,i3,a)") '"'//trim(fnamem(ibb,isyml,jsp))//trim(aaa), &
                     ibb+1,' pt ',ibb+1,' w lp'
                write(ifglts(jsp),'(",\")')
                write(ifglts(jsp),"(a,$)")trim(bchar)
             endif
2113      enddo
2112   enddo
       deallocate(diffeb,diff2eb)
2111 enddo
    do jsp=1,nspx
       ifglt=ifglts(jsp)
       write(ifglt,*)
       write(ifglt,'(a)')'set terminal x11'
       write(ifglt,'(a)')'replot'
       close(ifglts(jsp))
    enddo
    !   !! bnds.* is only for backward compatibility.
    !   open(newunit=ifbndo,file='bnds.'//trim(sname))
    !   write(ifbndo,"(i5,f10.5,i6)") sum(nqp_syml(1:nsyml)),basel,0
    !   do 3111 isyml = 1,nsyml
    !      ne = nqp_syml(isyml)
    !      write(ifbndo,"(i5)") ne
    !      do iq = 1,ne
    !         ikp = ikpoff(isyml) + iq
    !         do jsp = 1, nspx
    !            write(ifbndo,"(3f10.5,i6)") qplist(:,ikp),nevls(ikp,jsp)
    !            write(ifbndo,"(10f8.4)")(evlall(i,jsp,ikp),i=1,nevls(ikp,jsp))
    !         enddo
    !      enddo
    ! 3111 enddo
    !   write(ifbndo,"(i6)") 0
    !   close(ifbndo)
  end subroutine writeband
  pure integer function ichangesign(a,n)
    implicit none
    intent(in)::                    a,n
    integer:: i,n
    real(8):: a(n)
    ichangesign=-1
    do i=1,n-1
       if(a(i)*a(i+1) <0) then
          ichangesign=i
          exit
       endif
    enddo
  end function ichangesign
  subroutine writefs(evlall,eferm,spinweightsoc)!Fermi surface mode (eigenvalues in full BZ). No output variables.
    use m_lmfinit, only: nsp,nspc,lso,nspx
    use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
    use m_mkqp,only: bz_nabc
    use m_qplist,only:nkp,qplist
    use m_shortn3_qlat,only: shortn3_qlat,nout,nlatout
    implicit none
    logical:: cmdopt0,allband
    real(8):: ppin(3),eferm
    integer:: ip,i,isp,ififm,nbxx,iq,ib,nkk1,nkk2,nkk3,ifi,nx(3),ndhamx
    real(8):: rlatp(3,3),xmx2(3),vadd,qshort(3),evlall(:,:,:)
    real(8):: spinweightsoc(:,:,:)
    character*100::sss=''
    integer:: iout 
    nkk1=bz_nabc(1)
    nkk2=bz_nabc(2)
    nkk3=bz_nabc(3)
    allband  = cmdopt0('--allband')
    nx=shape(evlall)
    ndhamx=nx(1)
    do isp=1,nsp/nspc
       if(isp==1) open(newunit=ifi, file='fermiup.bxsf')
       if(isp==2) open(newunit=ifi, file='fermidn.bxsf')
       if(isp==1) open(newunit=ififm, file='fermiup.data')
       if(isp==2) open(newunit=ififm, file='fermidn.data')
       write(ififm,*) " # Generated by bndfp.F with --fermisurface"
       write(ififm,*) " # q-point is reduced in 1st BZ"
       write(ifi,*) "BEGIN_INFO"
       write(ifi,*) " # usage xcrysden --bxsf fermiup.bxsf"
       write(ifi,*) " # http://www.xcrysden.org/doc/XSF.html#2l.16"
       write(ifi,'(a,f9.5)') "  Fermi Energy:",eferm
       write(ififm,'(a,f9.5)') "  # Fermi Energy [eV]:",eferm
       if(lso==1) sss='                           spinup     spindn'
       write(ififm,*) " # qshort(3) band_energy[eV]", trim(sss)
       write(ifi,*) "END_INFO"
       write(ifi,*)"BEGIN_BLOCK_BANDGRID_3D"
       write(ifi,*)"  this_is_for_xcrysden"
       write(ifi,*)"  BEGIN_BANDGRID_3D_simple_test"
       nbxx=0
       do ib=1,ndhamx
          if(allband .OR. (minval(evlall(ib,isp,:))<eferm+0.5 .AND.maxval(evlall(ib,isp,:))>eferm-0.5)) then
             nbxx=nbxx+1
          endif
       enddo
       write(ifi,"(4x,i8)") nbxx
       write(ifi,"(4x,3i8)") nkk1,nkk2,nkk3
       write(ifi,"(4x,3(f9.5,1x))") 0d0,0d0,0d0
       write(ifi,"(4x,3(f9.5,1x))") qlat(:,1)
       write(ifi,"(4x,3(f9.5,1x))") qlat(:,2)
       write(ifi,"(4x,3(f9.5,1x))") qlat(:,3)
       !call shortn3_initialize(qlat)
       do ib=1,ndhamx
          if(allband .OR. (minval(evlall(ib,isp,:))<eferm+0.5 .AND. maxval(evlall(ib,isp,:))>eferm-0.5)) then
             write(ifi,"(a,i8)")"  BAND: ",ib
             write(ififm,"(a,i8)")"  # BAND: ",ib
             write(ifi,"(3x,10(x,f9.5))") (evlall(ib,isp,iq), iq=1,nkk1*nkk2*nkk3)
             ! to reduce q-point to qshort
             do iq=1,nkk1*nkk2*nkk3
                ppin=matmul(transpose(plat),qplist(:,iq))
                call shortn3_qlat(ppin) !bug (qlist(:,iq) ==>ppin 2022-6-8 tkotani)
                iout=1
                qshort(1:3)=  matmul(qlat(:,:), ppin+nlatout(:,iout))
                if(lso==1) write(sss,"(2f11.6)") spinweightsoc(ib,1:2,iq)
                write(ififm,"(4f13.5,a)") qshort(:),evlall(ib,isp,iq), trim(sss) !q-point reduced to 1stBZ
             enddo
          endif
       enddo
       write(ifi,*)
       write(ifi,*)"  END_BANDGRID_3D"
       write(ifi,*)"END_BLOCK_BANDGRID_3D"
    enddo
    close(ifi)
    close(ififm)
  end subroutine writefs
  subroutine writepdos(ext)! Readin pdosinput, and print out pdos files.
    use m_dstrbp,only: dstrbp
!    use m_lmfinit,only:lso
    implicit none
    integer:: ifip,ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,ndos,nkp &
         ,ibas,jsp,ifi,init,iend,ipts,j,ndos_,ichan,isp,itet,ksp,i,ib
    integer,allocatable::idtete(:,:),ipqe(:,:,:)
    real(8),allocatable:: evlall(:,:,:),dwgtall(:,:,:,:,:),pdosp(:,:),pdosalla(:,:,:,:)
    real(8)::eminp,emaxp,ef0,eee,eminp_,emaxp_
    character*3::charnum3
    real(8):: bin,eigen(4),vvv,wt,bin2
    character strn*120
    character*(*)::ext
    character(8)::xt
    logical::cmdopt2!,mlog
    integer, dimension(:),allocatable :: kpproc
    integer::numprocs,procid,ierr,itete,iteti,ispx,lso
    logical :: cmdopt0, idwmode
    real(8), allocatable :: dwgt4(:,:,:,:)
    integer :: iq, idt, idw
    include "mpif.h"
    idwmode = cmdopt0('--writedw')
    if(idwmode) print *, 'writedw mode ON'
    print *,' pdosdata file=','pdosdata.'//trim(ext)
    open(newunit=ifip,form='unformatted',file='pdosdata.'//trim(ext))
    read(ifip) ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,nkp,lso !ndos,nkp mar2015
    allocate(idtete(0:4,6*nkp),ipqe(nkk1,nkk2,nkk3))
    allocate(evlall(ndhamx,nspx,nkp))
    if(.not.idwmode) allocate(dwgtall(nchanp,nbas,ndhamx,nsp,nkp))
    if(idwmode) allocate(dwgt4(nchanp,nbas,ndhamx,4))
    read(ifip) idtete
    read(ifip) evlall
    if(.not.idwmode) read(ifip) dwgtall
    read(ifip) ef0
    close(ifip)
    call MPI_COMM_RANK( comm, procid, ierr )
    call MPI_COMM_SIZE( comm, numprocs, ierr )
    allocate (kpproc(0:numprocs), stat=ierr)
    call dstrbp(ntete,numprocs,1,kpproc(0))
    iteti = kpproc(procid)
    itete = kpproc(procid+1)-1
    eminp =-25.0/rydberg() !default values
    emaxp = 30.0/rydberg()
    ndos  = 5500
    if (cmdopt2('-emin=',strn)) then
       read(strn,*) eminp_
       eminp = eminp_/rydberg()
    endif
    if (cmdopt2('-emax=',strn)) then
       read(strn,*) emaxp_
       emaxp = emaxp_/rydberg()
    endif
    if(cmdopt2('-ndos=',strn))  then
       read(strn,*) ndos_
       ndos = ndos_
    endif
    bin = (emaxp - eminp) / (ndos - 1)
    vvv = ( 3d0  -  nsp ) / ( nkk1 * nkk2 * nkk3 * 6d0 )/4d0
    allocate(pdosalla(ndos,nsp,nchanp,nbas))
    tetrehedronloop:do itet = iteti, itete
       do isp = 1, nsp
          if(idwmode) then
             do idt = 1, 4
                 iq = idtete(idt,itet)
                 open(newunit=idw,file='dwgt.dir/dwgtk'//trim(xt(iq))//trim(xt(isp)),form='unformatted')
                 read(idw) dwgt4(:,:,:,idt)
                 close(idw)
             enddo
          endif
          ispx= merge(1,isp,lso==1)
          do ib = 1, nevmin
             eigen(1:4) = evlall(ib,ispx,idtete(1:4,itet))
             if( minval(eigen) > emaxp+ef0 ) cycle
             do ibas = 1,nbas
                do ichan = 1, nchanp
                   if(idwmode) then
                     wt = sum(dwgt4(ichan,ibas,ib,1:4)) * idtete(0,itet) * vvv
                   else
                     wt = sum(dwgtall(ichan,ibas,ib,isp,idtete(1:4,itet))) * idtete(0,itet) * vvv
                   endif
                   call slinz(wt,eigen,eminp+ef0,emaxp+ef0,pdosalla(1,isp,ichan,ibas),ndos)
                enddo
             enddo
          enddo
       enddo
    enddo tetrehedronloop
    call mpibc2_real( pdosalla, ndos*nsp*nchanp*nbas, 'writepdos_pdosalla' )
    if(procid==0) then
       allocate(pdosp (ndos,nchanp))
       bin2 = 2d0 * bin
       do isp =1,nsp
          do ibas=1,nbas
             open(newunit=ifi,file=trim('dos.isp'//char(48+isp)//'.site'//charnum3(ibas))//'.'//trim(ext))
             write(ifi,"('#lm ordering. See the end of lmdos. relative to efermi')")
             DOSfromfinitedifferenceofNOS:do ichan=1,nchanp
                do i = 2, ndos - 1
                   pdosp(i,ichan)=(pdosalla(i+1,isp,ichan,ibas) - pdosalla(i-1,isp,ichan,ibas))/bin2
                enddo
                pdosp(1,ichan)    = pdosp(2,ichan)
                pdosp(ndos,ichan) = pdosp(ndos-1,ichan)
             enddo DOSfromfinitedifferenceofNOS
             do ipts=1,ndos
                eee = eminp+ (ipts-1d0)*(emaxp-eminp)/(ndos-1d0)
                write(ifi,"(255(f13.5,x))")eee,pdosp(ipts,1:nchanp)
             enddo
             close(ifi)
          enddo
       enddo
       deallocate(pdosp,pdosalla)
    endif
  end subroutine writepdos
  subroutine writedossawada()
    use m_dstrbp,only: dstrbp
    !! DOS generator
    ! step 0.
    !     prepare ctrl.cubic
    !     ctrl.cubic should contain plat= 1 0 0 0 1 0 0 0 1.
    !     Set nk1=10 nk2=10 nk3=10 or someghing in ctrl.cubic
    !     Anything fine, because we just like to generate tetraf.dat qlistf.dat. Eg. ctrl.cubic I have.
    ! step 1.
    !    > job_pdos  cubic -np 4 --tetraw
    !    gives tetradata.dat and qlistf.dat (-np 1 is fine)
    ! step 2.
    !     Caluluate eigenvalues for given qlistf.dat
    !       Write eigenf.dat. One line for eigenvalues for q in qplistf.dat
    !       First line is for 'dimension, minimum energy, maximum energy, division'
    !     See how to read eigenf.dat below.
    ! step 3.
    !     Now we have tetraf.dat and eigenf.dat  .
    !     Run this binary generated from program writedossawada.
    !     >lmf-MPIK --wdsawada (see lmv7.F)
    !     Then you get dosf.dat file
    !! == readin dos input, print out dos
    implicit none
    integer:: ndhamx,nsp,nspx,nevmin,nchanp,nkk1,nkk2,nkk3,ntete,ndos,nkp &
         ,ibas,jsp,ifi,init,iend,ipts,j,isp,itet,ksp,i,ib,ifip,nbas
    integer,allocatable::idtete(:,:)
    real(8),allocatable:: evlall(:,:),pdosp(:,:),pdosalla(:,:),dwgtall(:,:,:)
    real(8)::eminp,emaxp,ef0,eee,eminp_,emaxp_
    character*100::strn
    !      character*100::filenm(2)
    real(8):: eigen(4),wt,tot,bin,bin2
    logical:: mlog
    integer, dimension(:),allocatable :: kpproc
    complex(8),allocatable:: ham(:,:,:)
    integer::numprocs,procid,ierr,itete,iteti,ikp
    include "mpif.h"
    open(newunit=ifip,form='unformatted',file='tetraf.dat')
    read(ifip) ndhamx,nkp,ntete
    allocate(idtete(0:4,6*nkp))
    read(ifip) idtete
    close(ifip)
    !!
    ndos=1000
    open(newunit=ifip,form='formatted',file='eigenf.dat')
    read(ifip,*) nevmin
    nbas=nevmin
    allocate(evlall(nevmin,nkp),dwgtall(nbas,nevmin,nkp))
    do ikp=1,nkp
       read(ifip,*) evlall(1:nevmin,ikp)
       !         read(ifip,*) dwgtall(ibas,1:nevmin,ikp)
       !         enddo
       !         write(6,"(i5,1000f10.3)") ikp,evlall(1:nevmin,ikp)
    enddo
    close(ifip)
    !!
    open(newunit=ifip,form='formatted',file='hamiltonian.dat')
    read(ifip,*)
    allocate(ham(nevmin,nevmin,nkp))
    read(ifip,*)ham(1:nevmin,1:nevmin,1:nkp)
    close(ifip)

    eminp =  minval(evlall)-0.5
    emaxp  = maxval(evlall)+0.5
    write(6,*) 'read eigenvalue data from eigenf.dat nkp=',nkp
    ef0=0d0
    mlog=.false.
    !      mlog = cmdopt('--mlog',6,0,strn) !--mlog here is taken by getarg.
    call MPI_COMM_RANK( comm, procid, ierr )
    call MPI_COMM_SIZE( comm, numprocs, ierr )
    allocate (kpproc(0:numprocs), stat=ierr)
    call dstrbp(ntete,numprocs,1,kpproc(0))
    iteti = kpproc(procid)
    itete = kpproc(procid+1)-1
    !! dostet ---
    !      call dostet(ndhamx,nsp,nspx,nevmin,nchanp*nbas,nkk1,nkk2,nkk3,ntete,idtete,evlall,
    !     &  dwgtall, ndos, eminp+ef0, emaxp+ef0,.false.,wkd,pdosall)
    bin = (emaxp - eminp) / (ndos - 1)
    !      vvv = ( 3d0  -  nsp ) / ( nkk1 * nkk2 * nkk3 * 6d0 )/4d0
    allocate(pdosalla(ndos,0:nbas))
    !! --- Loop over tetrahedra ---
    do itet = iteti, itete
       do ib = 1, nevmin
          eigen(1:4) = evlall(ib,idtete(1:4,itet))
          if( minval(eigen) > emaxp+ef0 ) cycle
          if( maxval(eigen) < eminp+ef0 ) cycle
          do ibas = 0,nbas
             if(ibas==0) wt = idtete(0,itet)
             if(ibas/=0) wt = sum(dwgtall(ibas,ib,idtete(1:4,itet))) * idtete(0,itet)
             call slinz(wt,eigen,eminp+ef0,emaxp+ef0,pdosalla(:,ibas),ndos)
          enddo
       enddo
    enddo
    call mpibc2_real( pdosalla, ndos, 'writedossawada_pdosalla' )
    if(procid==0) then        !master only
       allocate(pdosp (ndos,0:nbas))
       bin2 = 2d0 * bin
       open(newunit=ifi,file='dosf.dat')
       do ibas=0,nbas
          do i = 2, ndos - 1
             pdosp(i,ibas)=(pdosalla(i+1,ibas) - pdosalla(i-1,ibas))/bin2
          enddo
       enddo
       pdosp(1,:)    = pdosp(2,:)
       pdosp(ndos,:) = pdosp(ndos-1,:)
       tot= sum(pdosp(1:ndos,:))*bin
       do ipts=1,ndos
          eee = eminp+ (ipts-1d0)*(emaxp-eminp)/(ndos-1d0)
          write(ifi,"(10000(f13.5,x))")eee,(pdosp(ipts,ibas)/tot,ibas=0,nbas)
       enddo
       close(ifi)
    endif
  end subroutine writedossawada
endmodule m_writeband
