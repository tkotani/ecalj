module m_hamindex   !hamiltonian index read/write for successive GW calculaitons.
  use NaNum,only: NaN       !for initialization, but not working well
  use m_lgunit,only:stdo
  use m_ftox
  public:: Readhamindex
  integer,protected,public:: ngrp_original,pwmode,ndham !ngrpaf,
  integer,protected,public:: nqi=NaN, nqnum=NaN, ngrp=NaN, lxx=NaN, kxx=NaN,norbmto=NaN, &
       nqtt=NaN, ndimham=NaN, napwmx=NaN, lxxa=NaN, ngpmx=NaN, imx=NaN,nbas=NaN
  integer,allocatable,protected,public:: & ! offH (:), & !iclasstaf(:),
       ltab(:),ktab(:),offl(:), offlrev(:,:,:),ibastab(:), & !iclasst(:),
       iqimap(:),iqmap(:),igmap(:),invgx(:),miat(:,:),ibasindex(:), &
       igv2(:,:,:),napwk(:),igv2rev(:,:,:,:)! for rotation of evec       
  real(8),allocatable,protected,public:: & !symops_af(:,:,:), ag_af(:,:), &
       symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:), dlmm(:,:,:,:),qq(:,:), qtt(:,:),qtti(:,:)
  real(8),protected,public:: plat(3,3)=NaN,qlat(3,3)=NaN,zbak
  logical,protected,public:: readhamindex_init=.false., AFmode
  private
  logical,private:: debug=.false.
contains
  subroutine readhamindex() !Read info of PMT Hamiltoninan for wave rotation.
    logical::qpgexist
    integer:: ifi,nkt
!    logical::pmton
!    logical,save:: done=.false.
!    if(done) call rx('readhamindex is already done')
!    done=.true.
    readhamindex_init=.true.
    open(newunit=ifi,file='HAMindex',form='unformatted')
    read(ifi)ngrp,nbas,kxx,lxx,imx,ngpmx,norbmto,pwmode,zbak,ndham,AFmode
    allocate(symops(3,3,ngrp),ag(3,ngrp))
    allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp))
    write(stdo,ftox) 'ngrp=',ngrp
    read(ifi)symops,ag,invgx,miat,tiat,shtvg
    allocate( ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto) )
    allocate( offlrev(nbas,0:lxx,kxx))
    read(ifi) lxxa
    allocate( dlmm(-lxxa:lxxa, -lxxa:lxxa, 0:lxxa, ngrp))
    read(ifi) dlmm
    read(ifi) ibastab,ltab,ktab,offl,offlrev
    read(ifi) qpgexist
    if(.not.qpgexist) then
       close(ifi)
       return
    endif
    read(ifi) nqtt,nqi,nqnum
    allocate( qq(3,nqnum),qtt(3,nqtt),qtti(3,nqi)) !this was qq(3,nqnum*2) until Aug2012 when shorbz had been used.
    allocate(iqmap(nqtt),igmap(nqtt),iqimap(nqtt))
    read(ifi)qq,qtt,qtti,iqmap,igmap,iqimap !,ngvecp,ngvecprev
    read(ifi)plat,qlat,napwmx
    if(napwmx/=0)then !for APW rotation used in rotwvigg
       nkt=nqtt
       allocate( igv2(3,napwmx,nkt) )
       allocate( napwk(nkt))
       allocate( igv2rev(-imx:imx,-imx:imx,-imx:imx,nkt) )
       read(ifi) igv2,napwk,igv2rev
    endif
    close(ifi)
  end subroutine readhamindex
end module m_hamindex


module m_hamindexW   !Write hamiltonian index file 'HAMindex' for rdsigm2 and GW parts
  public m_hamindexW_init
contains
  subroutine m_hamindexW_init() !Set up m_hamiltonian. Index for Hamiltonian. --
    use NaNum,only: NaN        !for initialization, but not working well
    use m_lmfinit,only: pwmode=>ham_pwmode,pwemax,ldim=>nlmto,noutmx,nsp,stdo, &
         alat=>lat_alat,nl,nbas,ispec,sspec=>v_sspec,n0,nkap0,zbak,slabl,z
    use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat,rv_a_opos
    use m_suham,only: ndham=>ham_ndham !max dimension of hamiltonian +napwad (for so=0,2)
    use m_lmfinit,only:norbx,ltabx,ktabx,offlx
    use m_MPItk,only: master_mpi
    use m_lgunit,only:stdo
    use m_ftox
    use m_mksym,only: rv_a_osymgr,rv_a_oag,lat_nsgrp,iclasst,AFmode_mksym=>AFmode,& 
         ngrp,kxx,lxx,norbmto,AFmode,symops,ag,invgx,miat,tiat,shtvg,lxxa,dlmm,ltab,ktab,offl,offlrev,ibastab
    implicit none
    integer:: ngrp_original
    integer:: nqi, nqnum,nqtt, ndimham, napwmx, ngpmx, imx
    integer,allocatable::  offH (:), iqimap(:),iqmap(:),igmap(:),ibasindex(:), igv2(:,:,:),napwk(:),igv2rev(:,:,:,:)
    real(8),allocatable:: qq(:,:), qtt(:,:),qtti(:,:)
    integer:: ibas,k,l,ndim,ipr,nglob,off,offs,iorb,offsi,ib,is
    integer:: nkabc(3),nkp,lshft(3),napwx,ig,nini,nk1,nk2,nk3,ik1,ik2,ik3,ikt,nkt
    integer:: ifi, ifisym,i,ifiqibz,igg,iqq,iqi,irr,iqi_,jobgw
    integer:: iapw ,iprint,ngadd,igadd
    integer:: ngp, ifiqg,iq,nnn(3),ixx,ndummy,nqbz___ ,ifatomlist
    integer,allocatable:: iqtt(:), kv(:)
    real(8):: pwgmax, QpGcut_psi,qxx(3),qtarget(3),platt(3,3),q(3),qx(3),qqx(3)
    real(8):: dum,qb(3,3),ddd(3),ppin(3), rlatp(3,3),xmx2(3),qqq(3),diffs,ddf
    real(8),allocatable:: symtmp(:,:,:)
    logical:: qpgexist
    character(8)::  spid(nbas)
    integer:: ndima,lmxax,npqn,ificlass,nat,lmaxa,ipqn,ifinlaindx,isp,konf
    logical,save:: done=.false.
    real(8):: tolq
    character(1):: lorb(1:3)=['p','d','l'],dig(1:9)=['1','2','3','4','5','6','7','8','9'],&
         lsym(0:n0-1)=['s','p','d','f','g','5','6','7','8','9']
    call tcn('m_hamindex_init')
    inquire(file='QGpsi',EXIST=qpgexist) 
    QGpsimodeWriteHamindex: if(qpgexist) then
       open(newunit=ifiqg,file='QGpsi',form='unformatted') ! q on mesh and shortened q.
       read(ifiqg) nqnum, ngpmx ,QpGcut_psi, nqbz___, nqi
       if(allocated(qq)) deallocate(qq)
       nqtt=nqnum 
       nkt=nqtt
       allocate( qtti(3,nqi), qq(3,nqtt),iqtt(nqtt) )
       iqi=0
       do  iq = 1, nqnum
          read(ifiqg)  qxx,ngp,irr  ! read q and number of G vectors (irr=1 meand irreducible points)
          if(irr/=0) then
             iqi=iqi+1
             qtti(:,iqi)=qxx
             iqtt(iqi)=iq
          endif
          read(ifiqg)
          qq(:,iq)=qxx
          if(master_mpi)write(stdo,"(' qq=',i5,3f10.5)") iq,qq(:,iq)
       enddo
       close(ifiqg)
       allocate(iqmap(nqtt),igmap(nqtt),iqimap(nqtt))
       platt= transpose(plat) !inverse of qlat
       allocate(qtt(3,nqtt))
       qtt(:,1:nqtt)=qq(:,1:nqtt)
       do i=1,nqtt !Generate info for rotwv and write 
          qtarget(:)=qtt(:,i)
          do iqi=1,nqi
             q=qtti(:,iqi)
             iqq=iqtt(iqi)
             iqi_=iqi
             do ig=1,ngrp
                call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,ig),q)) ), qx)
                if(sum(abs(qx))<tolq()) then
                   igg=ig
                   ddf=sum(abs(matmul(platt,(qtarget-matmul(symops(:,:,ig),q)))))
                   if(ddf-nint(ddf)>1d-8)write(stdo,ftox)'qxqx',ftof(q),ftof(qtarget),ftof(qtarget-matmul(symops(:,:,ig),q))
                   goto 2012
                endif
             enddo
          enddo
          errorexitq: if(master_mpi) then !error exit
             write(stdo,"(a,3f7.3)")'gen_ham: qtarget cannot found. Need to add SYMGRP explicitly (for SO=1), '//&
                  'or You have to delete inconsistent QGpsi. qtarget=',qtarget
             write(stdo,*)'gen_hamindex: qtarget can not found by SYMOPS.'
             write(stdo,"('qq20 ',3d16.8,2x,3d16.8,2x,3d16.8)")q,qtarget,matmul(platt,(qtarget-matmul(symops(:,:,ig),q)))
             do ig=1,ngrp
                call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,ig),q)) ), qx)
                write(stdo,"('qqqq2 ',3d16.8,2x,3d16.8,2x,3d16.8)") qtarget-matmul(symops(:,:,ig),q),qx
             enddo
             call rx('gen_hamindex: you may need to repeat echo 1|qg4gw, when you changed SYMOPS.')
          endif errorexitq
2012      continue
          iqmap(i)=iqq
          iqimap(i)=iqi_
          igmap(i)=igg
       enddo
       ! === For rotation of APW.===
       if(master_mpi) call pshpr(0) !print index is pushed to be zero
       if(mod(pwmode,10)==0 .OR. pwemax<1d-8) then
          allocate(napwk(nkt),source=0)
          napwmx=0
       else ! for APW rotation.  ! ... Get igv2(3,iapw,ikt). pwmode>=10 only
          if(master_mpi) print *,' gen_hamindex goto APW part: pwmode pwemax=',pwmode,pwemax
          allocate(napwk(nkt))
          pwgmax = pwemax**.5
          do ikt=1,nkt
             qqq = merge(qq(:,ikt), 0d0, mod(pwmode/10,10)==1)
             call getgv2(alat,plat,qlat,qqq, pwgmax,1, napwk(ikt),dum)
          enddo
          napwmx=maxval(napwk)
          allocate( igv2(3,napwmx,nkt))
          do ikt = 1,nkt
             qqq = merge(qq(:,ikt),0d0,mod(pwmode/10,10) == 1)
             call getgv2(alat,plat,qlat,qqq, pwgmax,2, napwk(ikt),igv2(:,:,ikt)) 
          enddo
          imx = maxval([(maxval(abs(igv2(1:3,1:napwk(ikt),ikt))),ikt=1,nkt)])
          allocate( igv2rev(-imx:imx,-imx:imx,-imx:imx,nkt),source=999999 ) !Reverse table of igv2 --->igv2rev
          do ikt = 1,nkt
             do ig  = 1,napwk( ikt )
                nnn  = igv2(1:3, ig, ikt)
                igv2rev( nnn(1), nnn(2),nnn(3), ikt) = ig
             enddo
          enddo
       endif
    endif QGpsimodeWriteHamindex
    if(master_mpi) then  !          call writehamindex()
       if(done) call rx('writehamindex is already done')
       done=.true.
       open(newunit=ifi,file='HAMindex',form='unformatted')
       write(ifi)ngrp,nbas,kxx,lxx,imx,ngpmx,norbmto,pwmode,zbak,ndham,AFmode
       write(ifi)symops,ag,invgx,miat,tiat,shtvg
       write(ifi)lxxa
       write(ifi)dlmm
       write(ifi)ibastab,ltab,ktab,offl,offlrev !for rotation of MTO. recovered sep2012 for EIBZ for hsfp0
       write(ifi) qpgexist
       if(qpgexist) then
         write(ifi)nqtt,nqi,nqnum
         write(ifi)qq,qtt,qtti,iqmap,igmap,iqimap
         write(ifi)plat,qlat,napwmx
         if(napwmx/=0) write(ifi) igv2,napwk,igv2rev !for APW rotation used in rotwvigg
         write(ifi) alat,rv_a_opos
         close(ifi)
      endif
      call poppr !print index is poped.
    endif
    call tcx('m_hamindex_init')
  end subroutine m_hamindexW_init
end module m_hamindexW
