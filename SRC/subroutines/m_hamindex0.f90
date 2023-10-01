!> m_hAMIndex0 contains informatio of SYMOPS,LATTC,CLASS,NLAindx.
module m_hamindex0 
  use m_lmfinit,only: ham_pwmode,pwemax,ldim=>nlmto,noutmx,nsp_in=>nsp,lat_alat,ctrl_nbas=>nbas,ispec,n0,nkap0
  use m_lmfinit,only: zbak_read=>zbak,slabl,z,lmxa_i=>lmxa
  use m_lattic,only: lat_qlat,lat_plat,rv_a_opos
  use NaNum,only: NaN       !for initialization, but not working well
  integer,protected,public:: pwmode,ngrp=NaN, nbas=NaN
  integer,allocatable,protected,public:: ltab(:),ktab(:),offl(:),offlrev(:,:,:),ibastab(:), &
       iqimap(:),iqmap(:),igmap(:),invgx(:),miat(:,:),ibasindex(:), igv2(:,:,:),napwk(:),igv2rev(:,:,:,:),iclasst(:)
  real(8),allocatable,protected,public:: symops(:,:,:),ag(:,:),tiat(:,:,:),shtvg(:,:), dlmm(:,:,:,:),qq(:,:), &
       qtt(:,:),qtti(:,:),zz(:)
  real(8),protected,public:: plat(3,3)=NaN,qlat(3,3)=NaN,zbak
  real(8),protected,public::alat
  integer,protected,public::lmxax,nsp,ndima,norb,npqn,nclass,nphimx
  integer,allocatable,public:: konft(:,:,:),iqnum(:,:),lmxa(:),nlindx(:,:,:),pqn(:)
  character(9),allocatable,public:: caption(:)
  character(8),allocatable,public::  spid(:)
  integer,allocatable,public:: nindx(:),lindx(:),ibasindx(:)
  public:: M_hamindex0_init,Readhamindex0
  private
  logical,private:: debug=.false.
contains
  subroutine m_hamindex0_init() !-- Set up m_hamiltonian. Index for Hamiltonian. --
    ! m_hamindex0_init sets all module variables in m_hamindex0, then write to HAMindex0
    ! readhamindex0 read all module variables in m_hamindex0 from HAMindex0
    !
    ! No floating orbital case  !    ldim  = ham_ldham(1)
    use m_mksym,only: osymgr=>symops,oag=>ag,ngrpin=>ngrp,iclasstin=>iclasst 
    use m_mksym_util,only:mptauof
    use m_MPItk,only: master_mpi
    use m_density,only: pnzall,pnuall
    !r As you see in subroutine rotwvigg, the index for Hamiltonian reads as;
    !      do iorb=1,norbmto             !orbital-blocks are specified by ibas, l, and k.
    !        ibas  = ibastab(iorb)
    !         l    = ltab(iorb)
    !         k    = ktab(iorb)        !kappa (diffent MTO index for each l)
    !        init1 = offl(iorb)+1      !starting index for the block iorb
    !        iend1 = offl(iorb)+2*l+1  !end of the block for the iorb
    !      enddo
    implicit none
    integer:: ibas,k,l,ndim,ipr,nglob,off,offs,specw,fieldw,iorb,offsi,ib,is
    integer:: nkabc(3),nkp,lshft(3),napwx,ig,nini,nk1,nk2,nk3,ik1,ik2,ik3,ikt,nkt
    integer:: i_spacks,i_spackv,ifi,nbas_in,ifisym,i,ifiqibz,igg,iqq,iqi,irr,iqi_,jobgw
    integer:: iout,nout,nlatout(3,noutmx),iapw ,iprint,ngadd,igadd
    integer:: ngp, ifiqg,iq,nnn(3),ixx,ndummy,nqbz___ ,ifatomlist
    integer,allocatable:: ltabx(:,:),ktabx(:,:),offlx(:,:),iqtt(:), kv(:)
    real(8):: pwgmax, pwgmin, QpGcut_psi,qxx(3),qtarget(3),platt(3,3),q(3),qx(3),qqx(3)
    real(8):: dum,qb(3,3),ddd(3),ppin(3), tolq, rlatp(3,3),xmx2(3),qqq(3),diffs,ddf
    real(8),allocatable:: symtmp(:,:,:)
    logical:: siginit, qpgexist,debug=.false., llmfgw,prpushed
    integer:: ificlass,lmaxa,ipqn,ifinlaindx,isp,konf
    real(8) ::xx(5),pnu(n0,2),pnz(n0,2)
    integer:: ipb(ctrl_nbas),ipc(ctrl_nbas),ipcx(ctrl_nbas),mxint
    character lsym(0:n0-1)*1, lorb(3)*6, dig(9)*1, strn4*9
    data dig /'1','2','3','4','5','6','7','8','9'/
    data lsym /'s','p','d','f','g','5','6','7','8','9'/
    data lorb /'phi','phidot','phiz'/
    ! forget floating orbital case  !    ldim  = ham_ldham(1)
    nbas=ctrl_nbas !nbas_in
    ngrp=ngrpin   !note nsgrp given in mksym.F is without inversion.
    nsp=nsp_in
    plat=lat_plat
    qlat=lat_qlat
    allocate(symops,source=osymgr)
    allocate(ag,    source=oag)
    allocate(invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp),spid(nbas),lmxa(nbas),zz(nbas))
    do ib=1,nbas
       is =ispec(ib) 
       spid(ib) =slabl(is) 
       lmxa(ib) =lmxa_i(is) 
       zz(ib)   =z(is)
    enddo
    call mptauof(symops,ngrp,plat,nbas,rv_a_opos, iclasstin,miat,tiat,invgx,shtvg ) ! translation informations in miat tiat invgx, shtvg
    ndima = 0
    norb=0
    do  ib = 1, nbas
       pnz(:,1:nsp) = pnzall(:,1:nsp,ib)
       do l = 0, lmxa(ib)
          npqn = merge(3,2,pnz(l+1,1) /= 0)
          ndima = ndima + npqn*(2*l+1)
          norb= norb + npqn
       enddo
    enddo
    ! ndima is the number of augmented wave. total number of (principle,L,ibas)
    lmxax = maxval(lmxa)  
    allocate(konft(0:lmxax,nbas,nsp))
    do ib = 1, nbas
       pnu(:,1:nsp)=pnuall(:,1:nsp,ib) 
       pnz(:,1:nsp)=pnzall(:,1:nsp,ib) 
       do    isp = 1, nsp
          do  l  = 0, lmxa(ib)
             konft(l,ib,isp) = merge(mod(pnz(l+1,isp),10d0),pnu(l+1,isp), mod(pnz(l+1,isp),10d0)<pnu(l+1,isp).AND.pnz(l+1,isp)>0)
          enddo
       enddo
    enddo
    open(newunit=ifinlaindx,file='NLAindx.chk')
    write(ifinlaindx,'(''----NLAindx start---------------''/I6)') ndima
    npqn=3
    allocate(nlindx(npqn,0:lmxax,nbas),nindx(ndima),lindx(ndima),ibasindx(ndima),caption(ndima),pqn(ndima))
    iorb=0
    nlindx=-1
    ndima = 0
    do  ipqn = 1, 3
       do  ib = 1, nbas
          pnu=pnuall(:,1:nsp,ib) 
          pnz=pnzall(:,1:nsp,ib)
          do  l = 0, lmxa(ib)
             npqn = merge(3,2,pnz(l+1,1) /= 0)
             if(ipqn>npqn) cycle
             iorb=iorb+1
             konf= merge(mod(pnz(l+1,1),10d0), pnu(l+1,1),ipqn == 3)
             strn4 = dig(konf)//lsym(l)//'_'//lorb(ipqn)
             nlindx(ipqn,l,ib)=ndima
             nindx   (ndima+1:ndima+2*l+1)=ipqn
             lindx   (ndima+1:ndima+2*l+1)=l
             ibasindx(ndima+1:ndima+2*l+1)=ib
             caption (ndima+1:ndima+2*l+1)=strn4
             pqn(ndima+1:ndima+2*l+1)=konf !principle quantum number
             nphimx=max(nphimx,ipqn)
             write(ifinlaindx,'(i6,i3,i4,i6,4x,a)')ipqn,l,ib,     ndima,strn4
             ndima = ndima + (2*l+1)
          enddo
       enddo
    enddo
    close(ifinlaindx)
    if(norb/=iorb) call rx('m_hamindex0:norb/=iorb')
    alat=lat_alat
    npqn=3
    nclass=maxval(iclasstin)
    call wkonfchk(alat,plat,nbas,lmxax,lmxa,nsp,konft)
    open(newunit=ifi,file='HAMindex0',form='unformatted')
    write(ifi) alat,plat,qlat,nbas,lmxax,nsp,ngrp,ndima,norb,npqn,nclass,nphimx
    write(ifi) konft(0:lmxax,1:nbas,1:nsp),lmxa(1:nbas),nlindx(1:npqn,0:lmxax,1:nbas)
    write(ifi) iclasstin(1:nbas),spid(1:nbas),zz(1:nbas)
    write(ifi) nindx(1:ndima),lindx(1:ndima),ibasindx(1:ndima),caption(1:ndima),pqn(1:ndima)
    write(ifi) symops(1:3,1:3,1:ngrp),invgx(1:ngrp),shtvg(1:3,1:ngrp)
    close(ifi)
  end subroutine m_hamindex0_init
  subroutine readhamindex0()
    implicit none
    integer:: ifi,ibas,i
    open(newunit=ifi,file='HAMindex0',form='unformatted')
    read(ifi) alat,plat,qlat,nbas,lmxax,nsp,ngrp,ndima,norb,npqn,nclass,nphimx
    allocate( konft(0:lmxax,1:nbas,1:nsp),lmxa(1:nbas),nlindx(1:npqn,0:lmxax,1:nbas))
    allocate( iclasst(1:nbas),spid(1:nbas),zz(1:nbas))
    allocate( nindx(1:ndima),lindx(1:ndima),ibasindx(1:ndima),caption(1:ndima),pqn(1:ndima))
    allocate( symops(1:3,1:3,1:ngrp),invgx(1:ngrp),shtvg(1:3,1:ngrp))
    read(ifi) konft(0:lmxax,1:nbas,1:nsp),lmxa(1:nbas),nlindx(1:npqn,0:lmxax,1:nbas)
    read(ifi) iclasst(1:nbas),spid(1:nbas),zz(1:nbas)
    read(ifi) nindx(1:ndima),lindx(1:ndima),ibasindx(1:ndima),caption(1:ndima),pqn(1:ndima)
    read(ifi) symops(1:3,1:3,1:ngrp),invgx(1:ngrp),shtvg(1:3,1:ngrp)
    close(ifi)
  end subroutine readhamindex0
  subroutine wkonfchk(alat,plat,nbas,lmxax,lmxa,nsp,konf) ! Write konf.chk file
    ! We have no floating orbitals added.
    !i         :given site index including those orbitals
    !i   lmxax :global maximum of lmxa,  lmxa  :augmentation l-cutoff
    !i   konf  :principal quantum numbers
    implicit none
    integer :: nbas,nsp,lmxax,lmxa(nbas),konf(0:lmxax,nbas,nsp)
    real(8) :: plat(3,3),alat
    integer :: ifi,ib,isp
    open(newunit=ifi,file='konf.chk')
    write(ifi,*) ' ------------------------------------------- '
    write(ifi,"(2i4,' ! nbas lmxax=(max l for augmentaion)')") nbas,lmxax
    write(ifi,*) ' ------------------------------------------- '
    do isp = 1, nsp
       write(ifi,"(' -- ibas lmxa konf(s) konf(p) konf(d)... ', ' isp=',2i2)")isp
       do  ib = 1, nbas
          write(ifi,"('   ',99i4)")ib,lmxa(ib),konf(0:lmxa(ib),ib,isp)
       enddo
    enddo
    close(ifi)
  end subroutine wkonfchk
end module m_hamindex0
