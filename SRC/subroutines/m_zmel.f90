!> Get the matrix element zmel =  ZO^-1 <MPB psi|psi> , where ZO is ppovlz(inverse of overlap matrix) !  "call get_zmel_init" return zmel 
!  All dependencies (use foobar below ) are inputs (must be protected).
module m_zmel
  use m_genallcf_v3,only: nclass,natom,nspin,nl,nn,nnv,nnc, nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, niw
  use m_genallcf_v3,only: alat,delta,deltaw,esmr,iclass,nlnmv,nlnmc,icore,ncore,plat,pos,z,ecore,mnl=>nlnm,nl,nn,nlnmx,il,in,im
  use m_hamindex,only: ngrp, symgg=>symops,invg=>invgx
  use m_rdpp,only: Rdpp, nxx,lx,nx,mdimx,nbloch,cgr,ppbrd,nblocha,done_rdpp
  use m_read_bzdata,only: nqbz,nqibz,  qlat,ginv,qbz,qibz,wbz, done_read_bzdata
  use m_readhbe,only: nband
  use m_readQG,only: ngpmx,ngcmx,Readqg
  use m_readVcoud,only: zcousq,ngc,ngb !! zcousq is the eigenfuncition of the Coulomb matrix
  public:: get_zmel_init, Mptauof_zmel,Setppovlz,Setppovlz_chipm!,get_zmel_init1,get_zmel_init2 ! Call mptauof_zmel and setppovlz in advance to get_zmel_init
  complex(8),allocatable,protected,public :: zmel(:,:,:) ! OUTPUT: zmel(nbb,nmtot, nqtot) ,nbb:mixproductbasis, nmtot:middlestate, nqtot:endstate

  
  !private
  integer,protected:: nbb 
  integer,allocatable,protected :: miat(:,:)
  real(8),allocatable,protected :: tiat(:,:,:),shtvg(:,:), ppbir(:,:,:,:,:,:)
  complex(8),allocatable,protected :: ppovlz(:,:)
  

contains
  subroutine setppovlz(q,matz) ! Set ppovlz for given q
    intent(in)::       q,matz
    logical:: matz
    !    ppolvz(igb,ivcou)= (1    0 ) \times  zcousq(igb, ivcou)
    !                       (0 ppovl)
    !    If matz=F, no multiplication by ivcou.  Thus we have ppolz(igb,igb)
    real(8) :: q(3)
    complex(8),allocatable :: ppovl_(:,:),ppovl(:,:)!,ppovlzinv(:,:) !    logical:: eibz4x0
    integer:: i,ngc_r,ippovl0
    real(8):: qx(3),tolq=1d-8
    if(allocated(ppovlz)) deallocate(ppovlz)
    if(allocated(ppovl)) deallocate(ppovl)
    allocate( ppovl(ngc,ngc),ppovlz(ngb,ngb))
    open(newunit=ippovl0,file='PPOVL0',form='unformatted') !inefficient search for PPOVLO for given q
    do 
       read(ippovl0) qx,ngc_r
       if(sum(abs(qx-q))<tolq) then
          if(ngc_r/=ngc) call rx( 'readin ppovl: ngc_r/=ngc')
          read(ippovl0) ppovl
          exit
       endif
    enddo
    close(ippovl0)
    ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
    ppovlz(nbloch+1:nbloch+ngc,:)=matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
    deallocate(ppovl)
    nbb=ngb
  end subroutine setppovlz
  subroutine setppovlz_chipm(zzr,nmbas1) !Set ppovlz for chipm case
    intent(in)::             zzr,nmbas1
    integer::nmbas1
    complex(8):: zzr(ngb,nmbas1)
    if(allocated(ppovlz)) deallocate(ppovlz)
    allocate(ppovlz(ngb,nmbas1))
    ppovlz= zzr
    nbb=nmbas1
  end subroutine setppovlz_chipm
  subroutine mptauof_zmel(symops,ng)! Set miat,tiat,invgx,shtvg, and then get ppbir
    use m_mksym_util,only:mptauof
    use m_hamindex0,only: Readhamindex0,iclasst
    intent(in)::          symops,ng
    integer:: ng
    real(8):: symops(9,ng)
    integer,allocatable ::  invgx(:)
    call readhamindex0()
    allocate(invgx(ng),miat(natom,ng),tiat(3,natom,ng),shtvg(3,ng))
    call mptauof(symops,ng,plat,natom,pos,iclasst,miat,tiat,invgx,shtvg ) !Set miat,tiat,shtvg for space group.
    deallocate(invgx)
    call rdpp(ng,symops)  !return ppbrd:radial integrals and cgr:rotated cg coeffecients. 
    ppbafp_v2_zmel: block 
      integer :: is,irot,lmxax, ic, i,lb,nb,mb,lmb,i1,ibas,i2, np,lp,mp,lmp,n,l,m,lm
      lmxax=nl-1
      allocate(ppbir(nlnmx,nlnmx,mdimx,nclass,ng,nspin)) ! ppbir is rotated <Phi(SLn,r) Phi(SL'n',r) B(S,i,rot^{-1}(r))> by rotated cg coefficients cgr
      do irot = 1,ng
         do is = 1,nspin 
            do concurrent (ic=1:nclass)
               ibas = ic
               i = 0 !i = product basis index.
               do lb  = 0, lx (ibas)
                  do nb  = 1, nx (lb,ibas)
                     do mb  = -lb, lb
                        i    = i+1           !The number of product basis is  =(i at the end of loop).
                        lmb  = lb*lb + lb + mb + 1
                        do concurrent (i2 = 1:mnl(ic),i1 = 1:mnl(ic)) !phi1 phi2 index
                           np= in(i2,ic);  lp= il(i2,ic);  mp= im(i2,ic);  lmp=lp*lp + lp + mp + 1
                           n = in(i1,ic);  l = il(i1,ic);  m = im(i1,ic);  lm = l*l + l + m + 1
                           ppbir(i1,i2,i,ic,irot,is)=cgr(lm,lmp,lmb,irot)*ppbrd(l,n,lp,np,lb,nb,is+nspin*(ic-1))
                           !note ppbir is not necesssairy if i1,i2 are different spins
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
    endblock ppbafp_v2_zmel
  end subroutine mptauof_zmel
  
  subroutine get_zmel_init(q,kvec,irot,rkvec, ns1,ns2,ispm, nqini,nqmax,ispq, nctot,ncc, iprx,zmelconjg)
    use m_readeigen,only: readcphif 
    use m_readeigen,only: readgeigf
    use m_itq,only: itq,ntq
    intent(in)::           q,kvec,irot,rkvec, ns1,ns2,ispm, nqini,nqmax,ispq, nctot,ncc, iprx,zmelconjg
    !note: ncc can be =nctot or 0 
    ! ns1:ns2 is the range of middle states (count including core index (i=1,nctot). Thus kth valence is at nctot+k.
    !! kvec is in the IBZ, rk = Rot_irot(kvec)
    !! \parameter all inputs
    !! \parameter matrix <MPB psi|psi>
    !! igb: index of mixed product basis       at rkvec (or written as rk)
    !!   igb=1,ngb
    !!   ngb=nbloch+ngc  ngb: # of mixed product basis
    !!                   nbloch: # of product basis (within MTs)
    !!                   ngc: # of IPW for the Screened Coulomb interaction.
    !!                   igc is for given
    ! zmelt(itp|it,ib)= ZO^-1 <MPB(rkvec,ngb) phim(q-rkvec, ns1:ns2, ispm)|phiq(q,ncc+nqmax,ispq)> , or dconjg( ) when zmelconjg=T
    ! zmel = transpose(ppovlz,zmelt(:,it,itp))
    !
    ! matmul(symgg(:,:,irot),kvec)-rkvec can have difference of reciprocal vectors.
    !
    real(8),parameter::tolq=1d-8
    complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
    integer:: isp,ns1,ns2,nqmax,irot,ispq,ispm,nqini, nctot,ncc,ncnv,ncorec,nccc,mdim
    integer:: it,ia,kx,verbose,nstate,ispqk
    integer:: ngp1, ngp2, ngvecpB1(3,ngpmx),ngvecpB2(3,ngpmx),nadd(3)
    integer:: i,iap,ias,ib,ic,icp,nc,nc1,nv,ics,itp,iae,ims,ime
    real(8):: quu(3),q(3), kvec(3),rkvec(3),qkt(3),qt(3), qdiff(3)
    real(8) :: ppb(nlnmx,nlnmx,mdimx,nclass) ! ppb= <Phi(SLn,r-R)_q,isp1 |Phi(SL'n',r-R)_qk,isp2 B_k(S,i,rot^{-1}(r-R))>
    logical:: iprx
    logical:: zmelconjg
    integer,allocatable:: ngveccR(:,:)
    complex(8)::cphiq(nlmto,nband), cphim(nlmto,nband)
    complex(8):: geigq(ngpmx,nband),dgeigqk(ngpmx,nband)
  integer:: nmini, nmmax,invr,nt0,ntp0,nmtot,nqtot
  integer:: iasx(natom),icsx(natom),iatomp(natom),imdim(natom)
  real(8)::tr(3,natom)
  real(8)::qk(3),symope(3,3),shtv(3)
    ! nblocha     = number of optimal product basis functions for each class
    ! nlnmx     = maximum number of l,n,m
    ! nctot      = total no. allowed core states
    ! nbloch     = total no. product basis within MT
    !           if(mdimx /= maxval(mdim) ) call rx( 'psi2b_v3: wrong mdimx')
    !           if(sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v3: wrong nbloch')
    !           if(sum(nlnmv(iclass(1:natom)))/=nlmto) call rx( ' psi2b_v3:sum(nlnmv)/= nlmto')
    !           if(sum(ncore(iclass(1:natom)))/= nctot) call rx( "psicb_v3:sum(ncore) wrong")
    !           if(ncc/=0 .AND. ncc/=nctot) call rx( "psicb_v3: ncc/=0 and ncc/=ncctot")
      !! zmelp(igc(qi),it(qk),itp(q)) = <igc it(for q2+G2) |itp(for q1+G1)>
      !! NOTE: shtv = g(delta_{g^-1})
      !!-- zmelp(igc,it,itp) = <igc it(for G2)|itp(for G1)|> matrix element.
      !!   zmelp0(igc,it,itp)= <Gc' G2|G1> geig(G1,itp) geig^*(G2,it)
      !!   zmelp(igc,it,itp) = <Gc|Gc'>^-1 zmelp0(Gc',it,itp) 
      !!   (<Gc|Gc'>^-1 = ppovlinv)
      !! New ggg matrix <Gc |G1 G2> is introduced.
      !!
      !!    <Gc G2|G1> is equivalent to <-Gc+G1-G2>; described by ggg
      !! Readin input
      !!    ggg(1:nggg) = <Gc+G2-G1>
      !!    nvggg(3,1:nggg)   for Gc+G2-G1
      !!    nvgcgp2(3,ngcgp2) for Gc+G2
      !!    ppovlinv(ngc,ngc) <Gc|Gc> matrix
      !
      ! ggitp(Gc+G2)= \sum_G1 <(Symope(Gc)+G2)-G1> geigq1(G1,itp)*exp(-i*G1*shtv)*exp(-i(q-Gadd)*shtv)
      ! NOTE: nvgcgp2(:,igcgp2) means symope(Gc)+ G2
      !========================================================
      ! NOTE \bfr'= g (\bfr) +\delta_g. Then mapping is ROT[f(\bfr)]= f(g^-1(\bfr)+\delta_{g^-1})
    ! zmelp0(igc'(Gc'),it(G2),itp(G1)) = <Gc'G2|G1> geigq(G1,itp) geigqk*(G2,it) = <Gc' itp(G2)|it(G1)>
    
    if(allocated(zmel)) deallocate(zmel)
    qk =  q - rkvec ! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
    associate(cphitemp=> readcphif(q,ispq))    
      cphiq(1:nlmto,1:ntq) = cphitemp(1:nlmto,itq(1:ntq)) 
    endassociate
    cphim = readcphif(qk, ispm) 
    symope= symgg(:,:,irot)
    if(ngc/=0) then
      call readqg('QGpsi',q,     qt, ngp1, ngvecpB1) !q is mapped to qt in BZ
      call readqg('QGpsi',qk,   qkt, ngp2, ngvecpB2)
      qdiff = matmul(symope,kvec)-qt+qkt ! rkvec + qkt - qt is not zero. <M(rkvec) Phi(qk) |Phi(q)>
      nadd = nint(matmul(transpose(plat), qdiff)) !nadd: difference in the unit of reciprocal lattice vectors.
      block
        use m_read_ppovl,only: getppx2, ngvecc,ngcread
        integer:: igcgp2,nn(3),iggg,igp1,itp,igc,igp2,igcgp2i_(ngc,ngp2)
        call getppx2(qlat,kvec) ! read and allocate ppovlinv
        if(ngc/=ngcread) call rx( 'melpln2t: ngc/= ngcx by getppx:PPOVLG')
        allocate(ngveccR(1:3,1:ngc))
        call rotgvec(symope, 1, ngc, [ngc], qlat, ngvecc, ngveccR)
      endblock
      geigq   = readgeigf(q, ispq) !read IPW part at q   !G1 for ngp1
      dgeigqk = readgeigf(qk,ispm) !read IPW part at qk  !G2 for ngp2
      dgeigqk = dconjg(dgeigqk)
    endif

    
    ppb = ppbir(:,:,:,:,irot,ispq)           !MPB has no spin dependence
    invr  = invg(irot)       !invrot (irot,invg,ngrp) ! Rotate atomic positions invrot*R = R' + T
    tr    = tiat(:,:,invr)
    iatomp= miat(:,invr)
    shtv  = matmul(symope,shtvg(:,invr))
    imdim = [( sum(nblocha(iclass(1:ia-1)))+1  ,ia=1,natom)]
    iasx=[(sum(nlnmv(iclass(1:ia-1)))+1,ia=1,natom)]
    icsx=[(sum(ncore(iclass(1:ia-1))),ia=1,natom)]
    nmini = ns1          !starting index of middle state  (nctot+nvalence order)
    nmmax = ns2-nctot    !end      index of middle state  (nctot+nvalence order)
    nt0= nmmax-nmini+1
    ntp0=nqmax-nqini+1
    nmtot  = nctot + nt0     ! = phi_middle nmtot=ns2-ns1+1
    nqtot  = ncc   + ntp0    ! = phi_end
    ZmelBlock:block
      complex(8):: zmelt(1:nbloch+ngc,nmtot,nqtot)
      zmelt=0d0
      ZmelWithinMT: block !- Calculates <psi_q(itp) |psi_qk(it) B_k(rot(r-R))> 
        complex(8):: phasea(natom) 
        phasea = [(exp(-img *tpi* sum(kvec*tr(:,ia))),ia=1,natom)] 
        iatomloop: do concurrent(ia = 1:natom)
          ic    = iclass(ia)
          nc    = nlnmc(ic) !nlnmc      = number of l,n,m for core states
          nv    = nlnmv(ic) !nlnmv      = number of l,n,m for valence
          nc1   = nc + 1
          ncnv  = nc+nv  !ncore + nvalence
          iap   = iatomp(ia)  
          icp   = iclass(iap)
          ias   = iasx(ia)  !start of nlnmv for ia
          iae   = ias+nv-1  !end  
          ics   = icsx(ia)   
          ims   = imdim(iap) !start of PB for ia
          ime   = ims-1+nblocha(icp)
          mdim  = nblocha(icp)
          ncorec= merge(ncore(ic),0,nctot>0)
          nccc  = merge(ncore(ic),0,ncc>0)
          associate(      &
               dcphiqk => dconjg(cphim(ias:iae,nmini:nmmax)),&
               tdcphiqk=> transpose(dconjg(cphim(ias:iae,nmini:nmmax))),&
               cphiq   => cphiq(ias:iae,nqini:nqmax),&
               ppbc    => ppb(nc1:ncnv,icore(1:ncorec,ic),:,icp)) !all readin but we need only bands nmini: or nqini:
            do concurrent(i=1:mdim) !valence-valence  !=== this may be time-consuming block ==================
              zmelt(i-1+ims,nctot+1:nctot+nt0,ncc+1:ncc+ntp0)=phasea(ia) &
                   *matmul(tdcphiqk,matmul(transpose(ppb(nc1:ncnv,nc1:ncnv,i,icp)),cphiq))
            enddo
            do concurrent(it=1:ncorec) !core-valence
              zmelt(ims:ime, ics+it, ncc+1:ncc+ntp0)=phasea(ia)* matmul(transpose(ppbc(:,it,1:mdim)),cphiq(:,1:ntp0))
            enddo
            do concurrent(itp=1:nccc) !valence-core
              zmelt(ims:ime, nctot+1:nctot+nt0, ics+itp)=phasea(ia)*matmul(transpose(ppbc(:,itp,1:mdim)),dcphiqk(:,1:nt0))
            enddo                              !^^^^^^^^^phasea(ia) right? 2024-1-7 or phasea(ia) or dcongj(phasea(ia))?
          endassociate
        enddo iatomloop
      endblock ZmelWithinMT
      if(ngc/=0)then
        ZmelIPW:block  !> Mattrix elements <Plane psi |psi> from interstitial plane wave.
          use m_read_ppovl,only: igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze,nvgcgp2,ngcgp,ggg,ppovlinv
          integer:: igcgp2,nn(3),iggg,igp1,itp,igc,igp2,igcgp2i_(ngc,ngp2)
          complex(8):: zmelp0(ngc,nt0,ntp0),phase(ngc) , ggitp(ngcgp,ntp0),gggmat(ngcgp,ngp1)
          do concurrent(igcgp2=1:ngcgp,igp1=1:ngp1) !G synthesized 
            nn = ngvecpB1(:,igp1)- nvgcgp2(:,igcgp2) - nadd ! G1 -(Gc+G2) - Gadd.  Note that -Gadd= -rk + qt -qkt
            if(nn(1)<nxi .OR. nxe<nn(1) .OR. nn(2)<nyi .OR. nye<nn(2) .OR. nn(3)<nzi .OR. nze<nn(3)) cycle
            iggg = igggi(nn(1),nn(2),nn(3)) !inversion table
            if(iggg>=0) gggmat(igcgp2,igp1)=ggg(iggg)
          enddo
          do concurrent (igc=1:ngc,igp2=1:ngp2) !igc for B !G synthesized 
            nn = ngveccR(:,igc) + ngvecpB2(:,igp2)
            igcgp2i_(igc,igp2)=igcgp2i(nn(1),nn(2),nn(3))
          enddo
          phase(:)=[(exp( -img*tpi*sum((matmul(symope,kvec)+matmul(qlat,ngveccR(:,igc)))*shtv) ),igc=1,ngc)]
          !          associate(&
          !               geigq   => readgeigf(q, ispq),&         !read IPW part at q   !G1 for ngp1
          !               dgeigqk => dconjg(readgeigf(qk,ispqk))) !read IPW part at qk  !G2 for ngp2
            ggitp(:,:)= matmul(gggmat,geigq(1:ngp1,itq(nqini:nqmax)))
            do concurrent (itp= 1:ntp0) !=== this may be time-consuming block (or maynot)==================
              associate( ggitp_=>reshape([((ggitp(igcgp2i_(igc,igp2),itp),igc=1,ngc),igp2=1,ngp2)],shape=[ngc,ngp2]))
                zmelp0(:,:,itp)= matmul(ggitp_,dgeigqk(1:ngp2,nmini:nmmax))
              endassociate
            enddo
            forall(igc=1:ngc) zmelp0(igc,:,:)=phase(igc)*zmelp0(igc,:,:) 
            !          endassociate
          call matm(ppovlinv,zmelp0,zmelt(nbloch+1:nbloch+ngc,nctot+1:nctot+nt0,ncc+1:ncc+ntp0),ngc,ngc,ntp0*nt0)
        endblock ZmelIPW
      endif
      allocate(zmel(nbb,ns1:ns2, nqtot))
      call matm(dconjg(transpose(ppovlz)), zmelt, zmel,nbb,ngb,nmtot*nqtot) !MultiplePPOVLZ
      if(zmelconjg) zmel=dconjg(zmel)
    endblock ZmelBlock
  end subroutine get_zmel_init

    ! this is a copy of get_zmel_init to use blas
  subroutine get_zmel_init_gpu(q,kvec,irot,rkvec, ns1,ns2,ispm, nqini,nqmax,ispq, nctot,ncc, iprx,zmelconjg)
    use m_readeigen,only: readcphif 
    use m_readeigen,only: readgeigf
    use m_itq,only: itq,ntq
    use m_blas, only: m_op_c, m_op_n, m_op_t, zmm, zmm_batch !CPU or GPU versions specifed by macro __GPU
    implicit none
    intent(in)::           q,kvec,irot,rkvec, ns1,ns2,ispm, nqini,nqmax,ispq, nctot,ncc, iprx,zmelconjg
    real(8),parameter::tolq=1d-8
    complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
    integer:: isp,ns1,ns2,nqmax,irot,ispq,ispm,nqini, nctot,ncc,ncnv,ncorec,nccc,mdim
    integer:: it,ia,kx,verbose,nstate,ispqk
    integer:: ngp1, ngp2, ngvecpB1(3,ngpmx),ngvecpB2(3,ngpmx),nadd(3)
    integer:: i,iap,ias,ib,ic,icp,nc,nc1,nv,ics,itp,iae,ims,ime
    real(8):: quu(3),q(3), kvec(3),rkvec(3),qkt(3),qt(3), qdiff(3)
    real(8) :: ppb(nlnmx,nlnmx,mdimx,nclass) ! ppb= <Phi(SLn,r-R)_q,isp1 |Phi(SL'n',r-R)_qk,isp2 B_k(S,i,rot^{-1}(r-R))>
    logical:: iprx
    logical:: zmelconjg
    integer,allocatable:: ngveccR(:,:)
    complex(8)::cphiq(nlmto,nband), cphim(nlmto,nband)
    complex(8),allocatable:: geigq(:,:),dgeigqk(:,:)
    integer:: nmini, nmmax,invr,nt0,ntp0,nmtot,nqtot
    integer:: iasx(natom),icsx(natom),iatomp(natom),imdim(natom)
    real(8)::tr(3,natom)
    real(8)::qk(3),symope(3,3),shtv(3)
    integer :: ierr

    if(allocated(zmel)) then
      !$acc exit data delete(zmel)
      deallocate(zmel)
    endif

    SetByCPU :block
      qk =  q - rkvec ! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
      associate(cphitemp=> readcphif(q,ispq))    
        cphiq(1:nlmto,1:ntq) = cphitemp(1:nlmto,itq(1:ntq)) 
      endassociate
      cphim = readcphif(qk, ispm) 
      symope= symgg(:,:,irot)
      allocate(geigq(ngpmx,nband),dgeigqk(ngpmx,nband))
      if(ngc/=0) then
        call readqg('QGpsi',q,     qt, ngp1, ngvecpB1) !q is mapped to qt in BZ
        call readqg('QGpsi',qk,   qkt, ngp2, ngvecpB2)
        qdiff = matmul(symope,kvec)-qt+qkt ! rkvec + qkt - qt is not zero. <M(rkvec) Phi(qk) |Phi(q)>
        nadd = nint(matmul(transpose(plat), qdiff)) !nadd: difference in the unit of reciprocal lattice vectors.
        block
          use m_read_ppovl,only: getppx2, ngvecc,ngcread
          integer:: igcgp2,nn(3),iggg,igp1,itp,igc,igp2,igcgp2i_(ngc,ngp2)
          call getppx2(qlat,kvec) ! read and allocate ppovlinv
          if(ngc/=ngcread) call rx( 'melpln2t: ngc/= ngcx by getppx:PPOVLG')
          allocate(ngveccR(1:3,1:ngc))
          call rotgvec(symope, 1, ngc, [ngc], qlat, ngvecc, ngveccR)
        endblock
        geigq   = readgeigf(q, ispq) !read IPW part at q   !G1 for ngp1
        dgeigqk = readgeigf(qk,ispm) !read IPW part at qk  !G2 for ngp2
        dgeigqk = dconjg(dgeigqk)
      endif

      ppb = ppbir(:,:,:,:,irot,ispq)           !MPB has no spin dependence
      invr  = invg(irot)       !invrot (irot,invg,ngrp) ! Rotate atomic positions invrot*R = R' + T
      tr    = tiat(:,:,invr)
      iatomp= miat(:,invr)
      shtv  = matmul(symope,shtvg(:,invr))
      imdim = [( sum(nblocha(iclass(1:ia-1)))+1  ,ia=1,natom)]
      iasx=[(sum(nlnmv(iclass(1:ia-1)))+1,ia=1,natom)]
      icsx=[(sum(ncore(iclass(1:ia-1))),ia=1,natom)]
      nmini = ns1          !starting index of middle state  (nctot+nvalence order)
      nmmax = ns2-nctot    !end      index of middle state  (nctot+nvalence order)
      nt0= nmmax-nmini+1
      ntp0=nqmax-nqini+1
      nmtot  = nctot + nt0     ! = phi_middle nmtot=ns2-ns1+1
      nqtot  = ncc   + ntp0    ! = phi_end
    end block SetByCPU

    ZmelBlock:block
      complex(8):: zmelt(1:nbloch+ngc,nmtot,nqtot)
      complex(8), allocatable :: zmelt_d(:,:,:)
#ifdef __GPU
      attributes(device) :: zmelt, zmelt_d
#endif
      !$acc kernels
      zmelt=0d0
      !$acc end kernels
      ZmelWithinMT: block !- Calculates <psi_q(itp) |psi_qk(it) B_k(rot(r-R))> 
        complex(8):: phasea(natom) 
        complex(8), allocatable :: ppbvphiq_d(:,:,:), cphim_d(:,:), cphiq_d(:,:), ppbc_d(:,:,:), ppbv_d(:,:,:)
#ifdef __GPU
        attributes(device) :: ppbvphiq_d, cphim_d, cphiq_d, ppbc_d, ppbv_d
#endif
        phasea = [(exp(-img *tpi* sum(kvec*tr(:,ia))),ia=1,natom)] 
        iatomloop: do ia = 1, natom
          ic    = iclass(ia)
          nc    = nlnmc(ic) !nlnmc      = number of l,n,m for core states
          nv    = nlnmv(ic) !nlnmv      = number of l,n,m for valence
          nc1   = nc + 1
          ncnv  = nc+nv  !ncore + nvalence
          iap   = iatomp(ia)  
          icp   = iclass(iap)
          ias   = iasx(ia)  !start of nlnmv for ia
          iae   = ias+nv-1  !end  
          ics   = icsx(ia)   
          ims   = imdim(iap) !start of PB for ia
          ime   = ims-1+nblocha(icp)
          mdim  = nblocha(icp)
          ncorec= merge(ncore(ic),0,nctot>0)
          nccc  = merge(ncore(ic),0,ncc>0)

          !copy from CPU to GPU
          if(nt0>0) allocate(cphim_d(nv,nt0), source = cphim(ias:iae,nmini:nmmax))

          allocate(cphiq_d(nv,ntp0), source = cphiq(ias:iae,nqini:nqmax))
          ! valence-valence part
          if(nt0 > 0) then
            allocate(ppbvphiq_d(nv,ntp0,mdim))
            allocate(zmelt_d(nt0,ntp0,mdim))
            allocate(ppbv_d(nv,nv,mdim))
            ppbv_d(1:nv,1:nv,1:mdim) = dcmplx(ppb(nc1:ncnv,nc1:ncnv,1:mdim,icp)) 
            ierr = zmm_batch(ppbv_d, cphiq_d, ppbvphiq_d, nv, ntp0, nv, mdim, opA = m_op_T, sameB = .true.)
            ierr = zmm_batch(cphim_d, ppbvphiq_d, zmelt_d, nt0, ntp0, nv, mdim, alpha = phasea(ia), opA = m_op_C, sameA = .true.)
            !$acc kernels
            do i = 1, mdim
              zmelt(i-1+ims,nctot+1:nctot+nt0,ncc+1:ncc+ntp0) = zmelt_d(1:nt0,1:ntp0,i)
            enddo
            !$acc end kernels
            deallocate(zmelt_d, ppbvphiq_d, ppbv_d)
          endif

          ! core-valence part NOT TESTED
          allocate(zmelt_d(mdim,ntp0,max(1,ncorec)), ppbc_d(mdim,nv,max(1,ncorec)))
          !$acc kernels
          do i = 1, mdim
            ppbc_d(i,1:nv,1:max(1,ncorec)) = dcmplx(ppb(nc1:ncnv,1:(max(1,ncorec)),i,icp))
          enddo
          !$acc end kernels
          ierr = zmm_batch(ppbc_d, cphiq_d, zmelt_d, mdim, ntp0, nv, ncorec, alpha = phasea(ia), sameB = .true.)
          !$acc kernels
          do it = 1, ncorec
            zmelt(ims:ime,ics+it,ncc+1:ncc+ntp0) = zmelt_d(1:mdim,1:ntp0,it)
          enddo
          !$acc end kernels
          deallocate(zmelt_d)
          ! valence-core part NOT TESTED
          if(nt0 > 0) then
            allocate(zmelt_d(mdim,nt0,max(nccc,1)))
            ierr =  zmm_batch(ppbc_d, cphim_d, zmelt_d, mdim, nt0, nv, nccc, alpha = phasea(ia), sameB = .true.)
            !$acc kernels
            do itp = 1, nccc
              zmelt(ims:ime,nctot+1:nctot+nt0,ics+itp) = zmelt_d(1:mdim,1:nt0,itp)
            enddo
            !$acc end kernels
            deallocate(zmelt_d)
          endif
          if(nt0>0) deallocate(cphim_d)
          deallocate(cphiq_d, ppbc_d)
        enddo iatomloop
      endblock ZmelWithinMT

      if(ngc/=0 .and. nt0 > 0)then
        ZmelIPW:block  !> Mattrix elements <Plane psi |psi> from interstitial plane wave.
          use m_read_ppovl,only: igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze,nvgcgp2,ngcgp,ggg,ppovlinv,nnxi,nnxe, nnyi,nnye,nnzi,nnze,nggg
          integer:: igcgp2,nn(3),iggg,igp1,itp,igc,igp2,igcgp2i_(ngc,ngp2)
          complex(8):: zmelp0(ngc,nt0,ntp0),phase(ngc) , ggitp(ngcgp,ntp0),gggmat(ngcgp,ngp1)
          complex(8):: ggitp_all(ngc, ngp2, ntp0)

          phase(:)=[(exp( -img*tpi*sum((matmul(symope,kvec)+matmul(qlat,ngveccR(:,igc)))*shtv) ),igc=1,ngc)]  !prepared by CPU
          !$acc data copyin(dgeigqk, geigq, phase, ngvecpB1, ngvecpB2, ngveccR, nadd, ggg(1:nggg), nvgcgp2(1:3,1:ngcgp), &
          !$acc             igggi(nxi:nxe,nyi:nye,nzi:nze), igcgp2i(nnxi:nnxe,nnyi:nnye,nnzi:nnze), ppovlinv(1:ngb,1:nbb)) &
          !$acc      create(zmelp0, nn, gggmat, igcgp2i_, ggitp, ggitp_all)
          !$acc kernels loop independent collapse(2)
          do igcgp2 = 1, ngcgp
            do igp1 =1, ngp1
              nn(1:3) = ngvecpB1(1:3,igp1) - nvgcgp2(1:3,igcgp2) - nadd(1:3)
              if(nn(1)<nxi .OR. nxe<nn(1) .OR. nn(2)<nyi .OR. nye<nn(2) .OR. nn(3)<nzi .OR. nze<nn(3)) cycle
              iggg = igggi(nn(1),nn(2),nn(3))
              if(iggg>=0) gggmat(igcgp2,igp1)=ggg(iggg)
            enddo
          enddo
          !$acc end kernels
          !$acc kernels loop independent collapse(2)
          do igc = 1, ngc
           do igp2 = 1, ngp2
              nn(1:3) = ngveccR(1:3,igc) + ngvecpB2(1:3,igp2)
              igcgp2i_(igc,igp2) = igcgp2i(nn(1),nn(2),nn(3))
            enddo
          enddo
          !$acc end kernels
          ierr = zmm(gggmat, geigq(1,itq(nqini)), ggitp, ngcgp, ntp0, ngp1, ldB = ngpmx)

          !$acc kernels loop independent collapse(2)
          do igp2 = 1, ngp2
            do igc = 1, ngc
              ggitp_all(igc, igp2, 1:ntp0) = ggitp(igcgp2i_(igc,igp2), 1:ntp0)
            enddo
          enddo
          !$acc end kernels
          ierr = zmm_batch(ggitp_all, dgeigqk(1,nmini), zmelp0, ngc, nt0, ngp2, ntp0, ldB = ngpmx, sameB = .true.)
          !$acc kernels
          do igc = 1, ngc
            zmelp0(igc,1:nt0,1:ntp0) = phase(igc)*zmelp0(igc,1:nt0,1:ntp0)
          enddo
          !$acc end kernels
          allocate(zmelt_d(ngc,nt0,ntp0))
          ierr = zmm(ppovlinv, zmelp0, zmelt_d, ngc, ntp0*nt0, ngc) 
          !$acc kernels
          zmelt(nbloch+1:nbloch+ngc,nctot+1:nctot+nt0,ncc+1:ncc+ntp0) = zmelt_d(1:ngc,1:nt0,1:ntp0)
          !$acc end kernels
          !$acc end data
          deallocate(zmelt_d)
        endblock ZmelIPW
      endif

      allocate(zmel(nbb,ns1:ns2, nqtot))
      !$acc enter data create(zmel)
      !$acc host_data use_device(zmel)
      !$acc data copyin(ppovlz(1:ngb,1:nbb))
      ierr = zmm(ppovlz, zmelt, zmel, nbb, nmtot*nqtot, ngb, opA = m_op_C)
      !$acc end data
      if(zmelconjg) then
        !$acc kernels
        zmel = dconjg(zmel)
        !$acc end kernels
      endif
      !$acc end host_data

    endblock ZmelBlock
  end subroutine get_zmel_init_gpu

end module m_zmel
