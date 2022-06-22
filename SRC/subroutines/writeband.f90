subroutine writeband(eferm,evtop,ecbot)!, etolv,etolc)
  !      subroutine writeband(stdo,nkp,nsyml,nspx,nsp,ndhamx,
  !     & xdatt,nqp_syml,nqp2n_syml,nevls,evlall,qplist,labeli,labele,alat,
  !     & eferm,evtop,ecbot,sname,nqps_syml,nqpe_syml,nsymln,dqsyml,
  !     &     etolv,etolc)
  use m_lmfinit,only:stdo,nsp,alat=>lat_alat
  use m_qplist,only: nkp,nsyml,xdatt,nqp_syml,nqp2n_syml,qplist,labeli,labele, &
       nqps_syml,nqpe_syml,dqsyml,etolv,etolc
  use m_suham,only: ndham=>ham_ndham, ndhamx=>ham_ndhamx,nspx=>ham_nspx
  use m_bandcal,only:nevls,evlall
  use m_ext,only: sname
  !! == write band file ==
  implicit none
  !      integer,intent(in):: nkp,nsyml,nspx,ndhamx,nsp,  nevls(nkp,nspx),stdo,nsymln
  !      real(8),intent(in):: evlall(ndhamx,nspx,nkp),qplist(3,nkp),xdatt(nkp),alat  !vnow,
  !      integer,intent(in):: nqp_syml(nsyml),nqp2n_syml(nsyml),nqps_syml(nsyml),nqpe_syml(nsyml)
  real(8),intent(in):: eferm,evtop,ecbot ! evtop is max of n-th band. !evbot is bottom of bands upper than n+1
  !                                             ! if evtop <evbot ---> insulator
  !      real(8),intent(in):: dqsyml(nsyml),etolv,etolc
  !      character*20,intent(in)::labeli(nsyml),labele(nsyml)
  !      character(*),intent(in)::sname

  integer:: ifbndo,ikp,isyml,jsp
  character*300::filenm(2),bchar
  real(8):: rydberg=13.6058d0 !,vadd
  character*3::charnum3
  integer::ifbndsp(2),iprint,iq,i,ifile_handle,ifglt,ifbndsp_nearef(2),ifmass(2),ifmglt,ifmglt2
  character*300::aaa,addx
  real(8):: emin=-20d0,emax=20d0 !eV for default plotting.
  integer:: ikps,ne,ifi,ix,nee,ibb,ilt,imin,imax
  real(8),allocatable::diffeb(:),diff2eb(:)
  real(8)::polinta,tpiba,eee,eqq,qqq,dEdkatef
  character*100::fname,massfile,fname2
  integer,allocatable::ikpoff(:)
  real(8),allocatable::disoff(:)
  character(100),allocatable:: fnameb(:,:),fnamem(:,:,:)
  logical:: initii
  logical:: semiconband,metalband
  integer:: idat,ichangesign,ifglts(2),iqplist
  character(100)::acrossef
  character(13)::massd,mass2d
  logical:: scd
  real(8)::kef
  real(8),allocatable:: kabs(:)

  scd= evtop <ecbot
  tpiba = 2d0*4d0*datan(1d0)/alat

  !! ikpoffset
  allocate(ikpoff(nsyml+1),disoff(nsyml+1))
  ikpoff(1)=0
  disoff(1)=0
  do isyml = 2,nsyml+1
     ikpoff(isyml) = ikpoff(isyml-1)+ nqp_syml(isyml-1)+nqp2n_syml(isyml-1)
     disoff(isyml) = disoff(isyml-1)+dqsyml(isyml-1)
     write(stdo,"('ikpoff=',2i5)") isyml,ikpoff(isyml)
  enddo

  !! write bandplot.glt for gnuplot
  allocate(fnameb(nsyml,nspx))
  do jsp = 1, nspx
     ifglt =ifile_handle()
     ifglts(jsp)=ifglt
     fname='bandplot.isp'//char(48+jsp)
     open(unit=ifglt, file=trim(fname)//'.glt')
     write(ifglt,'(a)')'set terminal postscript enhanced color eps'
     write(ifglt,'(a)')'set output "'//trim(fname)//'.band.eps"'
     write(ifglt,'(a)')'set xzeroaxis'
     write(ifglt,'(a)')'set grid'
     write(ifglt,'(a)')'set ylabel "Energy-Efermi(eV)"'
     write(ifglt,'(a)')'# This is given written in subroutine writeband in lm7K/fp/bndfp.F'
     if(nspx==1) addx='"'
     if(nspx==2) addx=' isp='//char(48+jsp)//'"'
     write(ifglt,'(a)')'set title "Band '//trim(sname)//trim(addx)
     write(ifglt,'(a,F12.5,a,F12.5,a)') 'set yrange [',emin,':',emax,']'
     write(ifglt,'(a,F12.5,a)') 'set xrange [0.0:',disoff(nsyml+1),']'
     !!
     write(ifglt,'(a,$)') 'set xtics ('
     do isyml = 1,nsyml   !symmetry line start
        write(ifglt,'(a,F15.10,",\")')"'"//trim(labeli(isyml))//"'",disoff(isyml)
     enddo
     write(ifglt,'(a,F15.10,")")') "'"//trim(labele(nsyml))//"'",disoff(nsyml+1)
     !!
     write(ifglt,'(a,$)') 'set x2tics ('
     do isyml = 2,nsyml   !symmetry line end
        write(ifglt,'(a,F15.10,",\")')"'"//trim(labele(isyml-1))//"'",disoff(isyml)
     enddo
     write(ifglt,'(a,F15.10,")")')"'"//trim(labele(nsyml))//"'",disoff(nsyml+1)

     write(ifglt,'(a,F12.5)') 'tpia=2*3.1415926/',alat
     write(ifglt,"('plot \')")
     do isyml = 1,nsyml
        fnameb(isyml,jsp)='bnd'//charnum3(isyml)//'.'//'spin'//char(48+jsp)
        if(isyml/=1) write(ifglt,'(",\")')
        write(ifglt,"(a,$)")  '"'//trim(fnameb(isyml,jsp))//'" u ($2):($3) lt 1 pt 1 w lp'
     enddo
  enddo

  !! Write bnd**.isp*.glt, open qplist.dat
  iqplist = ifile_handle()
  open(iqplist,file='qplist.dat')
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
        ifbndsp(jsp) =ifile_handle()
        open(unit=ifbndsp(jsp),file=fnameb(isyml,jsp))
        write(ifbndsp(jsp),"('#  ',i5,f10.5,15x,'QPE(ev)',11x,'1st-deri')") nkp,eferm
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
              write(ifbndsp(jsp),"(i5, (1x,f12.8), 2(1x,f16.10),x,3f9.5)") &
                   i,xdatt(ikp),(evlall(i,jsp,ikp)-eferm)*rydberg,diffeb(iq), qplist(:,ikp)
              !! write qplist.dat (xaxis, q, efermi) used for band plot
              if(i==1 .AND. jsp==1) then
                 if(iq==1 .AND. isyml==1) write(iqplist,"(f16.10,' ! ef')") eferm
                 ikp = ikpoff(isyml)+iq
                 write(iqplist,"(f16.10,x,3f9.5,' !x q')") &
                      xdatt(ikp), qplist(:,ikp)
              endif
           enddo
           deallocate(diffeb)
           write(ifbndsp(jsp),*)
5113    enddo
        close(ifbndsp(jsp))
4113 enddo
4111 enddo
  close(iqplist)

  !$$$
  !$$$!! write massplotband.glt for gnuplot
  !$$$!! write massplot.glt for gnuplot (curvature, not the slope for metal).
  !$$$      allocate(fnamem(minval(nevls(:,:)),nsymln+1:nsyml,nspx))
  !      do jsp = 1, nspx
  !$$$c        ifmglt=ifile_handle()
  !$$$c        fname='massplot.isp'//char(48+jsp)
  !$$$c        open(unit=ifmglt,file=trim(fname)//'.glt')
  !$$$c        ifmglt2=ifile_handle()
  !$$$        fname2='massplotband.isp'//char(48+jsp)
  !$$$c        open(unit=ifmglt2,file=trim(fname2)//'.glt')
  !$$$c   set terminal postscript enhanced color eps
  !$$$c   set output 'temp.eps'
  !$$$c        write(ifmglt,'(a)')'set terminal postscript enhanced color eps'
  !$$$c        write(ifmglt,'(a)')'set output "'//trim(fname)//'.eps"'
  !$$$c        write(ifmglt,'(a,i1)')'set multiplot layout 1,',nsyml-nsymln
  !$$$c        write(ifmglt,'(a)')'set xzeroaxis'
  !$$$c        write(ifmglt,'(a)')'set grid'
  !$$$c        write(ifmglt,'(a)')'set ylabel "mass/mass(electron)"'
  !$$$c        write(ifmglt,'(a)')'set xlabel " |q|/(2pi/alat) "'
  !$$$c        write(ifmglt,'(a)')'# If necessary, see subroutine writeband given in lm7K/fp/bndfp.F for your purpose!'
  !$$$
  !$$$c$$$        write(ifglt,'(a)')'set terminal postscript enhanced color eps'
  !$$$c$$$        write(ifglt,'(a)')'set output "'//trim(fname2)//'.eps"'
  !$$$c$$$        write(ifglt,'(a,i1)')'set multiplot layout 1,',nsyml-nsymln
  !$$$c$$$        write(ifglt,'(a)')'set xzeroaxis'
  !$$$c$$$        write(ifglt,'(a)')'set grid'
  !$$$c$$$        write(ifglt,'(a)')'set ylabel " e - E(fermi) (eV)"'
  !$$$c$$$        write(ifglt,'(a)')'set xlabel " |q| " '
  !$$$c$$$        write(ifglt,'(a)')'# If necessary, see subroutine writeband given in lm7K/fp/bndfp.F for your purpose!'
  !$$$c$$$        if(nspx==1) addx=''
  !$$$c$$$        if(nspx==2) addx=' isp='//char(48+jsp)
  !$$$c$$$c        write(ifmglt,*)
  !        do isyml = nsymln+1,nsyml
  !$$$          ne = nqpe_syml(isyml)-nqps_syml(isyml)+1
  !$$$c$$$c          if(scd)         aaa='# insul u ($5):($8)'
  !$$$c$$$c          if(.not.scd)    aaa='# metal u ($5):($9)'
  !$$$c$$$c          if(.not.scd)    aaa='# metal u ($5):($8)'
  !$$$c$$$c          write(ifmglt,'(a)')'set title "Mass '//trim(sname)//trim(addx)//
  !$$$c$$$c     &    ' '//trim(labeli(isyml))//'--'//trim(labele(isyml))//' '//trim(aaa)//'"'
  !$$$c$$$c          write(ifmglt,'(a,F12.5,a)') 'set xrange [0.0:',xdatt(ikpoff(isyml+1))-disoff(isyml),']'
  !$$$c$$$c          write(ifmglt,"('plot \')")
  !$$$c$$$          write(ifglt,'(a)')'set title "Band '//trim(sname)//trim(addx)//
  !$$$c$$$     &    ' '//trim(labeli(isyml))//'--'//trim(labele(isyml))//'"'
  !$$$c$$$          write(ifglt,'(a,F12.5,a)') 'set xrange [0.0:',xdatt(ikpoff(isyml+1))-disoff(isyml),']'
  !$$$          write(ifglt,'(a,F12.5)') 'tpia=2*3.1415926/',alat
  !$$$          write(ifglt,"('plot \')")
  !$$$          ibb=0
  !$$$          ikps= ikpoff(isyml)+1
  !$$$          initii=.true.
  !          do i=1,minval(nevls(:,:))
  !$$$            eee = evlall(i,jsp,ikps)
  !$$$            semiconband = scd .and. evtop-etolv<eee .and. eee<ecbot+etolv .and. isyml>nsymln !check i-th band is neare VCM and CBM
  !$$$            idat= ichangesign(evlall(i,jsp,ikps:ikps+ne-1)-eferm,ne) !for metal crosspoint point across eferm
  !$$$            metalband  = idat>-1 .and. isyml>nsymln !check i-th band is near VCM and CBM
  !$$$            if(semiconband.or.metalband) then
  !               ibb = ibb+1
  !               fnamem(ibb,isyml,jsp)='Band'//charnum3(ibb)//'Syml'//charnum3(isyml)//'Spin'//char(48+jsp)//'.mass'
  !$$$c               if(.not.initii) write(ifmglt,'(",\")')
  !$$$c               if(scd)         aaa='" u ($5):($8) lt'
  !$$$c               if(.not.scd)    aaa='" u ($5):($9) lt'
  !$$$c               write(bchar,"(a,i2,a,i2,a,i3,a)") '"'//trim(fnamem(ibb,isyml,jsp))//trim(aaa),
  !$$$c     &         ibb,' pt ',ibb,' w lp ti "band=',i,'"'
  !$$$c               write(ifmglt,"(a,$)")trim(bchar)
  !$$$               if(.not.initii) write(ifglt,'(",\")')
  !$$$               aaa='" u (tpia*$5):($6) lt'
  !$$$               write(bchar,"(a,i2,a,i2,a,i3,a)") '"'//trim(fnamem(ibb,isyml,jsp))//trim(aaa),
  !$$$     &         ibb,' pt ',ibb,' w lp ti "band=',i,'"'
  !$$$               write(ifglt,"(a,$)")trim(bchar)
  !$$$               initii=.false.
  !$$$            endif
  !          enddo
  !$$$c          write(ifmglt,*)
  !$$$c          write(ifmglt,'(a)')'unset ylabel'
  !$$$          write(ifglt,*)
  !$$$          write(ifglt,'(a)')'unset ylabel'
  !$$$        enddo
  !$$$c        close(ifmglt)
  !$$$c        close(ifmglt2)
  !      enddo


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
           metalband = ichangesign(evlall(i,jsp,ikps:ikps+ne-1)-eferm,ne) >-1
           ! or metal crosspoint point across eferm
           ! heck i-th band is near VCM and CBM
           if(semiconband .OR. metalband) then
              ifmass(jsp) = ifile_handle()
              ibb = ibb+1
              if(metalband) then
                 imin=max(2,idat-1)  !diffeb is given from 2 to ne-1
                 imax=min(ne-1,idat+2)
                 allocate( kabs(imin:imax))
                 kabs(imin:imax) = (xdatt(ikps+imin-1:ikps+imax-1)-disoff(isyml))*tpiba ! tpiba is 2pi/alat. |k|. left-end is zero.
                 kef = polinta(0d0, evlall(i,jsp,ikps+imin-1:ikps+imax-1)-eferm, kabs, imax-imin+1) !we have k at Ef |k|=qef
                 dEdkatef = polinta(kef, kabs, diffeb(imin:imax), imax-imin+1)
                 deallocate(kabs)
                 !                 write(stdo,"('Effective Mass for metal: ',a,
                 !     &            ' iq isyml |k|_ef/(2pi/alat) dEdk 2*k/dEdk=', 2i4,' ',f8.3,' ',f8.3,' ',f8.3)")
                 !     &            trim(fnamem(ibb,isyml,jsp)),idat,isyml, kef/tpiba, dEdkatef, 2d0*kef/dEdkatef
              endif
              !               print *,'ffff ', ibb,isyml,jsp,trim(fnamem(ibb,isyml,jsp))
              open(unit=ifmass(jsp),file=trim(fnamem(ibb,isyml,jsp)))
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
                      (evlall(i,jsp,ikps+ix-1)-eferm)*rydberg, &
                      (evlall(i,jsp,ikps+ix-1)-evlall(i,jsp,ikps))*rydberg, &
                      mass2d,  massd,  trim(acrossef)
              enddo
              close(ifmass(jsp))
              write(aaa,"(a,f13.5,a)") '" u ($5/tpia+',disoff(isyml),'):($6) lt'
              write(bchar,"(a,i2,a,i2,a,i3,a)") '"'//trim(fnamem(ibb,isyml,jsp))//trim(aaa), &
                   ibb+1,' pt ',ibb+1,' w lp'
              write(ifglts(jsp),'(",\")')
              write(ifglts(jsp),"(a,$)")trim(bchar)
           endif
2113    enddo
2112 enddo
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
!   write(ifbndo,"(i5,f10.5,i6)") sum(nqp_syml(1:nsyml)),eferm,0
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

