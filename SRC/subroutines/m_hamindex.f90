!>Read HAMindex 
module m_hamindex   
  public:: Readhamindex
  integer,protected,public:: pwmode,ndham,lmxax, nqi, ngrp, lxx, kxx,norbmto, nqtt, ndimham, napwmx, ngpmx, imx,nbas
  integer,allocatable,protected,public:: ltab(:),ktab(:),offl(:), offlrev(:,:,:),ibastab(:),iqimap(:),iqmap(:),igmap(:)
  integer,allocatable,protected,public:: invgx(:),miat(:,:),ibasindex(:), igv2(:,:,:),napwk(:),igv2rev(:,:,:,:)
  real(8),allocatable,protected,public:: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:), dlmm(:,:,:,:),qq(:,:), qtt(:,:),qtti(:,:)
  real(8),protected,public:: plat(3,3),qlat(3,3),zbak
  logical,protected,public:: readhamindex_init=.false., AFmode
  integer,protected,public:: ngrpAF
  private
  logical,private:: debug=.false.
  integer,private::ngall
contains
  subroutine readhamindex() !Read info of PMT Hamiltoninan for wave rotation.
    use m_lgunit,only:stdo
    use m_ftox
    logical::qpgexist
    integer:: ifi
    readhamindex_init=.true.
    open(newunit=ifi,file='HAMindex',form='unformatted')
    read(ifi)ngrp,nbas,kxx,lxx,norbmto,pwmode,zbak,ndham,AFmode,ngrpAF
    ngall=ngrp+ngrpAF
    allocate(symops(3,3,ngall),ag(3,ngall),invgx(ngall),miat(nbas,ngall),tiat(3,nbas,ngall),shtvg(3,ngall))
    read(ifi)symops,ag,invgx,miat,tiat,shtvg
    allocate( ltab(norbmto),ktab(norbmto),offl(norbmto),ibastab(norbmto),offlrev(nbas,0:lxx,kxx))
    read(ifi) lmxax
    allocate( dlmm(-lmxax:lmxax, -lmxax:lmxax, 0:lmxax, ngall))
    read(ifi) dlmm
    read(ifi) ibastab,ltab,ktab,offl,offlrev
    read(ifi) qpgexist
    if(.not.qpgexist) then
       close(ifi)
       return
    endif
    read(ifi) nqtt,nqi,ngpmx
    allocate(qtt(3,nqtt),qtti(3,nqi)) 
    allocate(iqmap(nqtt),igmap(nqtt),iqimap(nqtt))
    read(ifi)qtt,qtti,iqmap,igmap,iqimap !,ngvecp,ngvecprev
    read(ifi)plat,qlat,napwmx,imx
    if(napwmx/=0)then !for APW rotation used in rotwvigg
       nqtt=nqtt
       allocate( igv2(3,napwmx,nqtt), napwk(nqtt),igv2rev(-imx:imx,-imx:imx,-imx:imx,nqtt))
       read(ifi) igv2,napwk,igv2rev
    endif
    close(ifi)
  end subroutine readhamindex
end module m_hamindex

!>Write 'HAMindex'
module m_hamindexW   
  public m_hamindexW_init
contains
  subroutine m_hamindexW_init() !Set up m_hamiltonian. Index for Hamiltonian. --
    use m_lmfinit,only: pwmode=>ham_pwmode,pwemax,ldim=>nlmto,noutmx,nsp,alat=>lat_alat,nbas,ispec,n0,nkap0,zbak,slabl,z
    use m_lmfinit,only: kxx,lxx,norbmto,lmxax,ltab,ktab,offl,offlrev,ibastab
    use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat,rv_a_opos
    use m_suham,only: ndham=>ham_ndham !max dimension of hamiltonian +napwad (for so=0,2)
    use m_lmfinit,only:norbx,ltabx,ktabx,offlx,lmxax
    use m_MPItk,only: master_mpi,comm
    use m_lgunit,only:stdo
    use m_ftox
    use m_mksym,only: rv_a_osymgr=>symops,rv_a_oag=>ag,iclasst,AFmode_mksym=>AFmode,ngrp,AFmode,symops,ag
    use m_mksym,only: invgx,miat,tiat,shtvg,dlmm,ngrpAF
!    use m_mpi,only: MPI__barrier
    implicit none
    integer:: nqi, nqtt, ndimham, napwmx, ngpmx, imx
    integer,allocatable::  offH (:), iqimap(:),iqmap(:),igmap(:),ibasindex(:), igv2(:,:,:),napwk(:),igv2rev(:,:,:,:)
    real(8),allocatable:: qtt(:,:),qtti(:,:)
    integer:: ibas,k,l,ndim,ipr,nglob,off,offs,iorb,offsi,ib,is,nkabc(3),nkp,lshft(3),napwx,ig,nini,nk1,nk2,nk3,ik1,ik2,ik3,ikt
    integer:: ifi, ifisym,i,ifiqibz,igg,iqq,iqi,irr,iqi_,jobgw,iapw ,iprint,ngadd,igadd
    integer:: ngp, ifiqg,iq,nnn(3),ixx,ndummy,nqbz___ ,ifatomlist
    integer,allocatable:: iqtt(:), kv(:)
    real(8):: pwgmax, QpGcut_psi,qxx(3),qtarget(3),platt(3,3),q(3),qx(3),qqx(3)
    real(8):: dum,qb(3,3),ddd(3),ppin(3), rlatp(3,3),xmx2(3),qqq(3),diffs,ddf
    real(8),allocatable:: symtmp(:,:,:)
    logical:: qpgexist
    character(8)::  spid(nbas)
    integer:: ndima,npqn,ificlass,ipqn,ifinlaindx,isp,konf,ngall,info
    logical,save:: done=.false.
    character(1):: lorb(1:3)=['p','d','l'],dig(1:9)=['1','2','3','4','5','6','7','8','9']
    character(1):: lsym(0:n0-1)=['s','p','d','f','g','5','6','7','8','9']
    character(256)::aaa
    include 'mpif.h'
    call tcn('m_hamindex_init')
    if(.not.master_mpi) goto 9999
    if(done) call rx('writehamindex is already done')
    done=.true.
    open(newunit=ifi,file='HAMindex',form='unformatted')
    write(ifi)ngrp,nbas,kxx,lxx,norbmto,pwmode,zbak,ndham,AFmode,ngrpAF
    ngall=ngrp+ngrpAF
    write(ifi)symops(:,:,1:ngall),ag(:,1:ngall),invgx(1:ngall),miat(:,1:ngall),tiat(:,:,1:ngall),shtvg(:,1:ngall)
    write(ifi)lmxax
    write(ifi)dlmm
    write(ifi)ibastab,ltab,ktab,offl,offlrev !for rotation of MTO. recovered sep2012 for EIBZ for hsfp0
    inquire(file='QGpsi',EXIST=qpgexist)  !qpgexist is for GW drivermode
    write(ifi) qpgexist
    if(.not.qpgexist) goto 9998
    QGpsimodeWriteHamindex:block
      real(8):: tolq
      open(newunit=ifiqg,file='QGpsi',form='unformatted') ! q on mesh and shortened q.
      read(ifiqg) nqtt, ngpmx ,QpGcut_psi, nqbz___, nqi
      allocate(qtt(3,nqtt), qtti(3,nqi), iqtt(nqtt),iqmap(nqtt),igmap(nqtt),iqimap(nqtt))
      iqi=0
      do iq = 1, nqtt
         read(ifiqg) qtt(:,iq),ngp,irr  ! read q and number of G vectors (irr=1 meand irreducible points)
         if(irr/=0) then
            iqi=iqi+1
            qtti(:,iqi)= qtt(:,iq)
            iqtt(  iqi)= iq
         endif
         read(ifiqg) 
      enddo
      close(ifiqg)
      platt= transpose(plat) !inverse of qlat
      do i=1,nqtt !Generate info for rotwv and write 
         qtarget(:)=qtt(:,i)
         do iqi=1,nqi
            q   =qtti(:,iqi)
            iqq =iqtt(iqi)
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
         Errorexitq: if(master_mpi) then
            do ig=1,ngrp
               call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,ig),q)) ), qx)
               write(stdo,"(3d16.8,2x,3d16.8,2x,3d16.8)") qtarget-matmul(symops(:,:,ig),q),qx
            enddo
            write(aaa,ftox) 'm_hamindexW: no qtarget by symops:',ftof(q),' ',ftof(qtarget)
            call rx('m_hamindex:'//trim(aaa))
         endif Errorexitq
2012     continue
         iqmap(i)=iqq
         iqimap(i)=iqi_
         igmap(i)=igg
      enddo
      ! === For rotation of APW.===
      if(master_mpi) call pshpr(0) !print index is pushed to be zero
      if(mod(pwmode,10)==0 .OR. pwemax<1d-8) then
         allocate(napwk(nqtt),source=0)
         napwmx=0
      else ! for APW rotation.  ! ... Get igv2(3,iapw,ikt). pwmode>=10 only
         if(master_mpi) print *,' gen_hamindex goto APW part: pwmode pwemax=',pwmode,pwemax
         allocate(napwk(nqtt))
         pwgmax = pwemax**.5
         do ikt=1,nqtt
            qqq = merge(qtt(:,ikt), 0d0, mod(pwmode/10,10)==1)
            call getgv2(alat,plat,qlat,qqq, pwgmax,1, napwk(ikt),dum)
         enddo
         napwmx=maxval(napwk)
         allocate(igv2(3,napwmx,nqtt))
         do ikt = 1,nqtt
            qqq = merge(qtt(:,ikt),0d0,mod(pwmode/10,10) == 1)
            call getgv2(alat,plat,qlat,qqq, pwgmax,2, napwk(ikt),igv2(:,:,ikt)) 
         enddo
         imx = maxval([(maxval(abs(igv2(1:3,1:napwk(ikt),ikt))),ikt=1,nqtt)])
         allocate(igv2rev(-imx:imx,-imx:imx,-imx:imx,nqtt),source=999999 ) !Reverse table of igv2 --->igv2rev
         do ikt = 1,nqtt
            do ig=1,napwk(ikt)
               igv2rev(igv2(1,ig,ikt), igv2(2,ig,ikt),igv2(3,ig,ikt), ikt) = ig
            enddo
         enddo
      endif
      write(ifi)nqtt,nqi,ngpmx
      write(ifi)qtt,qtti,iqmap,igmap,iqimap
      write(ifi)plat,qlat,napwmx,imx
      if(napwmx/=0) write(ifi) igv2,napwk,igv2rev !for APW rotation used in rotwvigg
      write(ifi) alat,rv_a_opos
    endblock QGpsimodeWriteHamindex
9998 continue
    close(ifi)
    call poppr !print index is poped.
9999 continue
    call MPI_Barrier( comm, info )
    call tcx('m_hamindex_init')
  end subroutine m_hamindexW_init
end module m_hamindexW
