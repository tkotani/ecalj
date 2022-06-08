subroutine eibzgen(nqibz,symgg,ngrp,qibze,iqxini,iqxend,qbz,nqbz,timereversal,ginv,iprintx, &
     nwgt,igx,igxt,eibzsym,timerout)
!  use m_shortn3,only: shortn3_initialize,shortn3
  use m_shortn3_qlat,only: shortn3_qlat,nout,nlatout
  use m_genallcf_v3,only:  plat
  use m_hamindex,only: pwmode
  implicit none
  !! === Obtain info for eibz symmetrization. See PRB81,125102-9 ===
  !! For GaAs 4x4x4 (with timereversal), we have IBZxBZ=10x640 is reduced to be IBZxEBZ=286.
  integer:: nqibz,ngrp,iqxini,iqxend,iq,ig,igxx,neibzx,ieibz,ibz,nwgtsum,itimer,ntimer
  integer:: neibz(iqxini:iqxend),nwgt(nqbz,iqxini:iqxend),ik,nqbz,i,it
  logical:: timereversal
  real(8):: symgg(3,3,ngrp),qibze(3,iqxini:iqxend),qbz(3,nqbz),qeibz(3,nqbz,iqxini:iqxend), &
       q(3),qx(3),qlat(3,3),qdiff(3),ginv(3,3),qxx(3),ddd,timer,qxi(3) , &
       sss(3,3),sumcheck1 !qout(3,nqbz,iqxini:iqxend,nqbz)
  integer:: igx(ngrp*2,nqbz,iqxini:iqxend),igxt(ngrp*2,nqbz,iqxini:iqxend),imm,im
  integer:: eibzsym(ngrp,-1:1,iqxini:iqxend),ibzxi(nqbz),ibzx,eallow(ngrp,-1:1),immxx
  logical:: iprintx,timerout

  integer,allocatable:: immx(:),mqbz(:)
  integer,allocatable:: nxx(:,:,:),ibznxx(:,:),nww(:,:),igx_(:,:),igxt_(:,:)

  real(8)::    epd=1d-5
  integer:: nkey(3),isig,kkk,kkk3(3),ik1,ik2,ik3,ix,ifi,ifiese
  integer,allocatable:: ieord(:)
  integer,allocatable:: key(:,:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
  real(8),allocatable:: qbzl(:,:)
  logical :: inotime,nexist,debug=.false.
!  integer,parameter:: noutmx=48
  integer:: iout,iapw !,nout,nlatout(3,noutmx)
  real(8)::ppin(3)
  real(8):: rlatp(3,3),xmx2(3),tolq=1d-8
  call minv33tp(plat,qlat)
  ! ccccccccccccccccccc
  !      print *,'ginv test'
  !      ginv=ginv/2d0
  !      write(6,"('ginv=',3f9.4)")ginv(:,:)
  ! ccccccccccccccccccc
  if(iprintx) write(6,"('eibzgen: TimeReversalSwitch ngrp= ',l,i3)") timereversal,ngrp
  ntimer=1
  if(timereversal) ntimer=-1
  sss(1,:)=(/1d0,0d0,0d0/)
  sss(2,:)=(/0d0,1d0,0d0/)
  sss(3,:)=(/0d0,0d0,1d0/)
  sumcheck1= sum(abs(symgg(:,:,1)-sss))
  if(sumcheck1>tolq) call rx( 'eibzgen: symgg(:,:,1) is not E')

  !!=== Get key and nkey for each ix. (similar with those in readqg,readeigen) ===
  allocate(qbzl(3,nqbz),key(3,0:nqbz),ieord(nqbz))
  key(:,0)=0 !dummy
  do ibz=1,nqbz
     call rangedq(matmul(ginv,qbz(:,ibz)), qbzl(1,ibz))
     ! cccccccccccccccc
     !            write(6,"('ibz=',i5,3f8.4)")ibz,qbzl(:,ibz)
     ! cccccccccccccccc
  enddo
  do ix =1,3
     call sortea(qbzl(ix,:),ieord,nqbz,isig)
     ik=0
     do i=1,nqbz
        kkk= ( qbzl(ix,ieord(i))+0.5*epd )/epd  !kkk is digitized by 1/epd
        if(i==1 .OR. key(ix,ik)<kkk) then
           ik=ik+1
           key(ix,ik) = kkk
           !               write(6,*)ix, ik,i, key(ix,ik), qbzl(ix,ieord(i))
        elseif (key(ix,ik)>kkk) then
           write(6,*)ix, ik,i, key(ix,ik), qbzl(ix,ieord(i))
           call rx( 'iqindx: bug not sorted well')
        endif
     enddo
     nkey(ix)=ik
  enddo
  deallocate(ieord)
  !!  key is reallocated. inverse mattping, iqkkk
  allocate( kk1(nkey(1)),kk2(nkey(2)),kk3(nkey(3)) )
  kk1(:) = key(1,1:nkey(1))
  kk2(:) = key(2,1:nkey(2))
  kk3(:) = key(3,1:nkey(3))
  deallocate(key)
  allocate( iqkkk(nkey(1),nkey(2),nkey(3)) )
  if(debug) write(6,*)' initqqq nqbz=',nqbz
  do i=1,nqbz
     kkk3= (qbzl(:,i)+0.5*epd)/epd !kkk is digitized by 1/epd
     call tabkk(kkk3(1), kk1,nkey(1), ik1)
     call tabkk(kkk3(2), kk2,nkey(2), ik2)
     call tabkk(kkk3(3), kk3,nkey(3), ik3)
     iqkkk(ik1,ik2,ik3)=i
     if(debug) write(6,"(' iqkkk ik1,ik2,ik3,=',i4,3i3,x,3i6)")i,ik1,ik2,ik3,kkk3(:)
  enddo
  if(debug) then
     write(6,*)'kk1=',kk1
     write(6,*)'kk2=',kk2
     write(6,*)'kk3=',kk3
  endif
  deallocate(qbzl)
  !! -----------------------------------------------------------
  !      open(newunit=ifiese,file='PWMODE')
  !      read(ifiese,*)pwmode
  !     close(ifiese)
!  if(pwmode>0 .AND. pwmode<10) call shortn3_initialize(qlat)

  !! ===  main iq loop ===
  eibzsym = 0
  do iq = iqxini,iqxend
     if(debug) write(6,"('iq =',i7,'  xxxxxxxxxxxxxxxxxxxxxxx')") iq
     q = qibze(:,iq) !q means k in eq.50 of PRB81,125102-9
     !! Allowed operation to keep q,   eallow(ig,it) \in EIBZ(q)
     eallow=0d0
     do it=1,ntimer,-2
        timer=dble(it)
        do ig=1,ngrp
           qdiff= timer*matmul(symgg(:,:,ig),q) - q
           call rangedq(matmul(ginv,qdiff), qxx)
           if(debug) write(6,*)'ig qdiff',ig,qdiff,q,ginv,qxx
           ! cccccccccccccccccccccccccccccccccccccccccccc
           if(pwmode>0 .AND. pwmode<10) then
              !     check wheter q is on 1st BZ boundary or not
              ppin=matmul(transpose(plat),q) !qlat-based fractional coodinate
              call shortn3_qlat(ppin) !shortn3(ppin,noutmx, nout,nlatout)
              if(ig/=1 .AND. nout>1) then
                 do iout=1,nout
                    if(debug) write(*,"(a,3i5,f10.4,3f8.4)")'rrrrn1 =',nlatout(:,iout), &
                         sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2), &
                         matmul(qlat(:,:),ppin+nlatout(:,iout))
                 enddo
                 print *
                 cycle
              endif
           endif
           ! ccccccccccccccccccccccccccccccccccccccccccc
           if(sum(abs(qxx))<tolq) then
              eallow(ig,it)=1
           endif
        enddo
     enddo
     eibzsym(:,:,iq)=eallow ! === eibzsym is eallow -----

     ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !        write(6,*)'xxx nqbz=',nqbz,ngrp,ntimer, 8*nqbz*ngrp/1000000000.
     ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     !! === Generate all nxx (integered q with epditit) generated from ibz. ===
     allocate(nxx(3,nqbz,ngrp*2),immx(nqbz)) !,nww(nbz,ngrp*2))
     allocate(igx_(ngrp*2,nqbz),igxt_(ngrp*2,nqbz))
     immxx=0
     do ibz=1,nqbz
        if(debug) write(6,"(' qbz=',i3,' ',3f9.4,' ',3f9.4)")ibz,qbz(:,ibz),matmul(ginv,qbz(:,ibz))
        imm=0
        do it=1,ntimer,-2 !ntimer is 1 or -1.
           timer = dble(it)
           do ig=1,ngrp
              if(eallow(ig,it)/=1) cycle !eallow=1 keeps q invariant. Allowed operation of EIBZ(q)
              qx = timer*matmul(symgg(:,:,ig),qbz(:,ibz))
              qxi= matmul(ginv,qx)
              call rangedq(qxi, qxx) ! qxx(i) takes from 0 to 1.
              do im=1,imm
                 if(sum(abs(nxx(:,ibz,im) - (qxx(:)+.5*epd)/epd)) <2) then
                    !                 nww(ibz,im) = nww(ibz,im)+1
                    goto 3012      !<2 for safe
                 endif
              enddo
              imm=imm+1
              nxx(:,ibz,imm)= (qxx+.5*epd)/epd  !integer nxx(i) takes from 0 to 1/epd.
              igx_ (imm,ibz) = ig  !inefficient memory usage, I think.
              igxt_(imm,ibz) = it  !
              if(debug) write(6,"('imm ibz it ig nxx qxx qxi=',4i7,x,3i6,x, 3f8.4,x,3f8.4,' qx=',3f8.4)") &
                   imm,ibz,ig,it ,nxx(:,ibz,imm), qxx, qxi, qx
              !            write(6,*)'qbz=',qbz(:,ibz)
              !            write(6,*)'qx =',qx
              !            write(6,*)'qxi=',qxi
              !! inversion table:
              !!             ibz<---from qbzrid
              !! we can skip such ibz.
3012          continue
           enddo
        enddo
        immx(ibz)=imm
        if(imm>immxx) immxx=imm
     enddo

     if(debug) write(6,*)'eibz aaaaaaa'
     !        call cputid(0)

     !! inequivalent points.
     allocate(mqbz(nqbz))
     mqbz=0

     !! may2015 there zero-cleare is moved from next ix loop. (bug fix). correct?
     nwgt(:,iq)=0
     igx(:,:,iq) =999991
     igxt(:,:,iq)=999992
     !!
     do ix=1,ntimer,-2 !try no time reversal for ix=1
        if(ix== 1) inotime= .TRUE. 
        if(ix==-1) inotime= .FALSE. 
        do ibz=1,nqbz
           if(debug) write(6,"('1111 ix,ibz',i3,i5,i5)")ix,ibz
           if(mqbz(ibz)==1) cycle
           do im=1,immx(ibz)
              if(inotime) then
                 if(igxt_(im,ibz) == -1) cycle
              elseif(inotime) then
                 if(igxt_(im,ibz) ==  1) cycle
              endif
              kkk3 = nxx(:,ibz,im)
              if(debug) write(6,*)'uuuuu ibz,im kkk3=',ibz,im,kkk3
              call tabkk(kkk3(1), kk1,nkey(1), ik1)
              call tabkk(kkk3(2), kk2,nkey(2), ik2)
              call tabkk(kkk3(3), kk3,nkey(3), ik3)
              if(ik1==-999999 .OR. ik2==-999999 .OR. ik3==-999999) cycle !! june2016
              if(debug) write(6,*)'uuuuu ik=',ik1,ik2,ik3
              ibzx = iqkkk(ik1,ik2,ik3)
              if(mqbz(ibzx)==1) then
                 cycle
              else
                 mqbz(ibzx)=1
                 nwgt(ibz,iq)= nwgt(ibz,iq) + 1
                 igx (nwgt(ibz,iq),ibz,iq) = igx_(im,ibz)
                 igxt(nwgt(ibz,iq),ibz,iq) = igxt_(im,ibz)
              endif
           enddo
        enddo
        if(debug) then
           do ibz=1,nqbz
              write(6,"('2222 ix,ibz mqbz nqbz=',i2,i3,i5,i5)")ix,ibz,mqbz(ibz),nwgt(ibz,iq)
           enddo
        endif
        if(sum(nwgt(:,iq))==nqbz) then
           if(ix==1) eibzsym(:,-1,iq)=0
           goto 1202
        endif
        write(6,*)' mptauo: we use time reversal'
     enddo
     call rx( 'eibzgen: bug. nwgt sum is not nqbz')
1202 continue

     ! cccccccccccccccccccccccccccccccccc
     !         do ibz=1,200 !nqbz
     !           if(nwgt(ibz,iq)/=0) then
     !             write(6,"('yyy0: ',i8,2x,25(i3,i2))") ibz,(igx(i,ibz,iq),igxt(i,ibz,iq),i=1,nwgt(ibz,iq))
     !           endif
     !         enddo
     ! cccccccccccccccccccccccccccccccccc

     !          write(6,*)'pppppp nwgt(ibz,iq),nqbz=',iq,sum(nwgt(:,iq)),nqbz
     !$$$!! === Search other ibz generated from ibz. ===
     !$$$!!   ibznxx(ibz,imm) means "ibz rotated by igx(imm,ibz,iq),igxt(imm,ibz,iq)"
     !$$$        allocate(mqbz(nqbz))
     !$$$        do ibzx=1,nqbz
     !$$$          do ibz=1,ibzx !search equivalent (ibzx,1) with first appeard (ibz,im) .
     !$$$          do im=1,immx(ibz)
     !$$$             if(sum(abs(nxx(:,ibzx,1) - nxx(:,ibz,im))) <2) then
     !$$$               nwgt(ibz,iq) = nwgt(ibz,iq) +  1
     !$$$               igx (nwgt(ibz,iq),ibz,iq) = igx_(im,ibz)
     !$$$               igxt(nwgt(ibz,iq),ibz,iq) = igxt_(im,ibz)
     !$$$               goto 3312
     !$$$             endif
     !$$$          enddo
     !$$$          enddo
     !$$$ 3312     continue
     !$$$        enddo
     deallocate(igx_,igxt_)
     !        write(6,*)'eibz bbbbbb'
     !        call cputid(0)

     neibzx=0
     do ibz=1,nqbz
        if(nwgt(ibz,iq) /=0) then
           neibzx=neibzx+1
           qeibz(:,neibzx,iq)=qbz(:,ibz)
           write(6,"('   ibz qeibz =',i7,3f10.5, ' nwgt=',i7)")ibz,qeibz(:,neibzx,iq),nwgt(ibz,iq)
        endif
     enddo
     neibz(iq) = neibzx !number of eibz for iq
     if( iprintx ) then
        write(6,"('iq=',i8,' # of EIBZ: Full=',i8, &
             ' Used(TimeR 1 or -1)=',i3,'=',i2,'+',i2,' neibz= ',i7)")iq,sum(abs(eallow(:,:))), &
             sum(eibzsym(:,:,iq)),sum(eibzsym(:,1,iq)),sum(eibzsym(:,-1,iq)),neibz(iq)
        !          write(6,"('eibz: iq neibz nqbz= ',i3,3f11.5,3i7)") iq,q,neibz(iq),nqbz
     endif
     deallocate(nxx,immx,mqbz)
  enddo
  !!
  nwgtsum = sum(nwgt(1:nqbz,iqxini:iqxend))
  timerout=timereversal
  if(minval(igxt)==1) timerout= .FALSE. 

  if(iprintx) then
     write(6,"(' nqbz,  sum(neibz(iq)), sum(ngwt)=sum(nqbz)= ',10i7)") &
          nqbz,sum(neibz(iqxini:iqxend)), nwgtsum, nqbz*(iqxend-iqxini+1)
  endif
  !      stop 'test stop -------- end of eibzgen ----------'
end subroutine eibzgen

!!------------------------------------------------
subroutine tabkk(kkin, kktable,n, nout)
  integer, intent(in) :: n,kkin, kktable(n)
  integer, intent(out) :: nout
  integer:: i,mm,i1,i2
  i1=1
  i2=n
  if(kkin==kktable(1)) then
     nout=1
     return
  elseif(kkin==kktable(n)) then
     nout=n
     return
  endif
  do i=1,n
     mm=(i1+i2)/2
     if(kkin==kktable(mm)) then
        nout=mm
        return
     elseif(kkin>kktable(mm)) then
        i1=mm
     else
        i2=mm
     endif
  enddo
  !$$$      do i=1,n
  !$$$         if(kkin==kktable(i)) then
  !$$$            nout=i
  !$$$            return
  !$$$         endif
  !$$$      enddo
  ! un2016
  !      call rx( 'takk: error')
  !      write(6,*) i1,i2,kkin
  !      write(6,*) kktable(i1),kktable(i2)
  nout=-999999
end subroutine tabkk
