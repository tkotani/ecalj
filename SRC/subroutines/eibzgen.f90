! === Obtain info for eibz symmetrization. See PRB81,125102-9 ===
! For GaAs 4x4x4 (with timereversal), we have IBZxBZ=10x640 is reduced to be IBZxEBZ=286.
subroutine eibzgen(nqibz,symgg,ngrp,qibze,iqxini,iqxend,qbz,nqbz,timereversal,ginv,iprintx, &
     nwgt,igx,igxt,eibzsym,timerout)
  use m_genallcf_v3,only:  plat
  use m_hamindex,only: pwmode
  use m_nvfortran,only: findloc
  implicit none
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
  real(8),parameter::    epd=1d-5
  integer:: nkey(3),isig,kkk,kkk3(3),ik1(1),ik2(1),ik3(1),ix,ifi,ifiese
  integer,allocatable:: ieord(:)
  integer,allocatable:: key(:,:),kk1(:),kk2(:),kk3(:),iqkkk(:,:,:)
  real(8),allocatable:: qbzl(:,:)
  logical :: inotime,nexist,debug=.false.
  integer:: iout,iapw 
  real(8)::ppin(3)
  real(8):: rlatp(3,3),xmx2(3),tolq=1d-8
  call minv33tp(plat,qlat)
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
  enddo
  do ix =1,3
     call sortea(qbzl(ix,:),ieord,nqbz,isig)
     ik=0
     do i=1,nqbz
        kkk= nint( qbzl(ix,ieord(i))/epd)  !kkk is digitized by 1/epd
        if(i==1 .OR. key(ix,ik)<kkk) then
           ik=ik+1
           key(ix,ik) = kkk !    write(6,*)ix, ik,i, key(ix,ik), qbzl(ix,ieord(i))
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
     kkk3= nint(qbzl(:,i)/epd) !kkk is digitized by 1/epd
     ik1= findloc(kk1,value=kkk3(1))
     ik2= findloc(kk2,value=kkk3(2))
     ik3= findloc(kk3,value=kkk3(3))
     iqkkk(ik1(1),ik2(1),ik3(1))=i
     if(debug) write(6,"(' iqkkk ik1,ik2,ik3,=',i4,3i3,x,3i6)")i,ik1,ik2,ik3,kkk3(:)
  enddo
  deallocate(qbzl)
  !! ===  main iq loop ===
  eibzsym = 0
  do iq = iqxini,iqxend
     q = qibze(:,iq) !q means k in eq.50 of PRB81,125102-9
     write(6,"('iq =',i7,3f10.5,' xxxxxxxxxxxxxxxxx')") iq,q
     !! Allowed operation to keep q,   eallow(ig,it) \in EIBZ(q)
     eallow=0d0
     do it=1,ntimer,-2
        timer=dble(it)
        do ig=1,ngrp
           qdiff= timer*matmul(symgg(:,:,ig),q) - q
           call rangedq(matmul(ginv,qdiff), qxx)
           if(sum(abs(qxx))<tolq) then
              eallow(ig,it)=1
           endif
        enddo
     enddo
     eibzsym(:,:,iq)=eallow ! === eibzsym is eallow -----
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
                 if(sum(abs(nxx(:,ibz,im) - nint(qxx/epd))) <2) then
                    !                 nww(ibz,im) = nww(ibz,im)+1
                    goto 3012      !<2 for safe
                 endif
              enddo
              imm=imm+1
              nxx(:,ibz,imm)= nint(qxx/epd)  !integer nxx(i) takes from 0 to 1/epd.
              igx_ (imm,ibz) = ig  !inefficient memory usage, I think.
              igxt_(imm,ibz) = it  !
              !if(debug) write(6,"('immx ibz it ig nxx qxx qxi=',4i7,x,3i6,xx,3f26.14)") &
              !     imm,ibz,ig,it ,nxx(:,ibz,imm), qxx!qxi,qx
              !! inversion table:      ibz<---from qbzrid    !! we can skip such ibz.
3012          continue
           enddo
        enddo
        immx(ibz)=imm
        if(imm>immxx) immxx=imm
     enddo
     if(debug) write(6,*)'eibz aaaaaaa'
     !! inequivalent points.
     allocate(mqbz(nqbz))
     mqbz=0
     nwgt(:,iq)=0
     igx(:,:,iq) =999991
     igxt(:,:,iq)=999992
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
              ik1= findloc(kk1,value=kkk3(1))
              ik2= findloc(kk2,value=kkk3(2))
              ik3= findloc(kk3,value=kkk3(3))
              if(ik1(1)+ik2(1)+ik3(1)==0) cycle 
              if(debug) write(6,*)'uuuuu ik=',ik1,ik2,ik3
              ibzx = iqkkk(ik1(1),ik2(1),ik3(1))
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
     deallocate(igx_,igxt_)
     neibzx=0
     do ibz=1,nqbz
        if(nwgt(ibz,iq) /=0) then
           neibzx=neibzx+1
           qeibz(:,neibzx,iq)=qbz(:,ibz)
           write(6,"('   ibz qeibz =',i6,3f10.5, '  nwgt=',i4)")ibz,qeibz(:,neibzx,iq),nwgt(ibz,iq)
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
  nwgtsum = sum(nwgt(1:nqbz,iqxini:iqxend))
  timerout=timereversal
  if(minval(igxt)==1) timerout= .FALSE. 
  if(iprintx) then
     write(6,"(' nqbz,  sum(neibz(iq)), sum(ngwt)=sum(nqbz)= ',10i7)") &
          nqbz,sum(neibz(iqxini:iqxend)), nwgtsum, nqbz*(iqxend-iqxini+1)
  endif   !      stop 'test stop -------- end of eibzgen ----------'
end subroutine eibzgen
