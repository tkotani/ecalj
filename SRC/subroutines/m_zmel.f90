!> Get the matrix element zmel =  ZO^-1 <MPB psi|psi> , where ZO is ppovlz(inverse of overlap matrix)
!  "call get_zmel_init" return zmel 
!  All dependencies (use foobar below ) are inputs (must be protected).
module m_zmel
  use m_genallcf_v3,only: nclass,natom,nspin,nl,nn,nnv,nnc, nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, niw
  use m_genallcf_v3,only: alat,delta,deltaw,esmr,iclass,nlnmv,nlnmc,icore,ncore,plat,pos,z,ecore,mnl=>nlnm,nl,nn,nlnmx,il,in,im
  use m_hamindex,only: ngrp, symgg=>symops,invg=>invgx
  use m_rdpp,only: Rdpp, nxx,lx,nx,mdimx,nbloch,cgr,ppbrd,nblocha,done_rdpp
  use m_readeigen,only: Readcphif 
  use m_read_bzdata,only: nqbz,nqibz,  qlat,ginv,qbz,qibz,wbz, done_read_bzdata
  use m_readhbe,only: nband
  use m_itq,only: itq,ntq
  use m_readQG,only: ngpmx,ngcmx,Readqg
  use m_readVcoud,only: zcousq,ngc,ngb !! zcousq is the eigenfuncition of the Coulomb matrix
  ! OUTPUT: zmel(nbb,nmtot, nqtot) ,nbb:mixproductbasis, nmtot:middlestate, nqtot:endstate
  complex(8),allocatable,protected,public :: zmel(:,:,:)
  public:: Get_zmel_init, Mptauof_zmel,Setppovlz,Setppovlz_chipm ! Call mptauof_zmel and setppovlz in advance to get_zmel_init
  private
  integer,allocatable,private :: miat(:,:)
  real(8),allocatable,private :: tiat(:,:,:),shtvg(:,:)
  real(8),allocatable,private :: ppbir(:,:,:,:,:,:)
  complex(8),allocatable,private :: ppovlz(:,:)
  real(8),private:: q_bk(3)=1d10,qk_bk(3)=1d10
  logical,private:: init=.true.
  complex(8),allocatable,private :: cphiq(:,:), cphim(:,:)
  integer:: nkmin, nkqmin, isp_k, isp_kq,nmtot,nqtot,ispq_bk,ispm_bk
  logical:: debug=.false. 
  integer:: nbb,nbbx=0!nbb:1st dimension of zmel. MPB
  real(8), parameter :: pi=3.1415926535897932D0
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
    open(newunit=ippovl0,file='PPOVL0',form='unformatted')
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
  end subroutine setppovlz
  subroutine setppovlz_chipm(zzr,nmbas1) !Set ppovlz for chipm case
    intent(in)::             zzr,nmbas1
    integer::nmbas1
    complex(8):: zzr(ngb,nmbas1)
    if(allocated(ppovlz)) deallocate(ppovlz)
    allocate(ppovlz(ngb,nmbas1))
    ppovlz= zzr
    nbbx=nmbas1
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
      allocate(ppbir(nlnmx,nlnmx,mdimx,nclass,ng,nspin)) ! ppbir is rotated <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)> by rotated cg coefficients cgr
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
    use m_readeigen,only : readgeigf
    intent(in)::           q,kvec,irot,rkvec, ns1,ns2,ispm, nqini,nqmax,ispq, nctot,ncc, iprx,zmelconjg
    ! ns1:ns2 is the range of middle states (count including core index (i=1,nctot). Thus kth valence is at nctot+k.
    !! Get zmel= <phiq(q,ncc+nqmax,ispq) |phim(q-rkvec, ns1:ns2, ispm) MPB(rkvec,ngb)> ZO^-1
    !! kvec is in the IBZ, rk = Rot_irot(kvec)
    !! \parameter all inputs
    !! \parameter matrix <MPB psi|psi>
    !! igb: index of mixed product basis       at rkvec (or written as rk)
    !!   igb=1,ngb
    !!   ngb=nbloch+ngc  ngb: # of mixed product basis
    !!                   nbloch: # of product basis (within MTs)
    !!                   ngc: # of IPW for the Screened Coulomb interaction.
    !!                   igc is for given
    ! zmelt= <itp|it,ib>
    ! zmel(igb,it,itp) = transpose(dcojng(ppovlz))*zmelt(:,it,itp)
    !   matmul(symgg(:,:,irot),kvec)-rkvec can have difference of reciprocal vectors.
    !
    !note: ncc can be =nctot or 0 
    real(8),parameter::tolq=1d-8
    complex(8),parameter:: img=(0d0,1d0),tpi= 8d0*datan(1d0)
    integer:: isp,nmmax,ns1,ns2,nqmax,irot,ispq,ispm,nmini,nqini, nctot,ncc,ncnv
    integer:: it,ia,kx,verbose,nstate,imdim(natom),nt0,ntp0,invr, iatomp(natom),ispqk
    integer:: ngp1, ngp2, ngvecpB1(3,ngpmx),ngvecpB2(3,ngpmx),nadd(3)
    integer:: iasx(natom),i,iap,ias,ib,ic,icp,nc,nc1,nv,ics,itp,icsx(natom),iae,ims,ime
    real(8):: quu(3),q(3), kvec(3),rkvec(3),symope(3,3),shtv(3),tr(3,natom),qk(3),qkt(3),qt(3), qdiff(3)
    complex(8):: phase(natom) , geig1(ngpmx,nband),geig2(ngpmx,nband)
    real(8),allocatable :: ppb(:,:,:,:)
    complex(8),allocatable::  zmelt(:,:,:)
    logical:: debug=.false.,iprx
    logical,optional:: zmelconjg
    nmini = ns1          !starting index of middle state  (nctot+nvalence order)
    nmmax = ns2-nctot    !end      index of middle state  (nctot+nvalence order)
    qk =  q - rkvec ! qk = q-rk. rk is inside 1st BZ, not restricted to the irreducible BZ
    invr  = invg(irot)       !invrot (irot,invg,ngrp) ! Rotate atomic positions invrot*R = R' + T
    tr    = tiat(:,:,invr)
    iatomp= miat(:,invr)
    symope= symgg(:,:,irot)
    shtv  = matmul(symope,shtvg(:,invr))
    allocate( ppb(nlnmx,nlnmx,mdimx,nclass)) ! ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
    ppb = ppbir(:,:,:,:,irot,ispq) !rotated ppb. MPB has no spin dependence
    imdim = [( sum(nblocha(iclass(1:ia-1)))+1  ,ia=1,natom)]
    phase = [(exp(img *tpi* sum(kvec*tr(:,ia))),ia=1,natom)] 
    nt0= nmmax-nmini+1
    ntp0=nqmax-nqini+1
    nmtot  = nctot + nt0     ! = phi_middle nmtot=ns2-ns1+1
    nqtot  = ncc   + ntp0    ! = phi_end   
    if(init) then
       allocate( cphiq(nlmto,nband), cphim(nlmto,nband))
       init=.false.
    endif
    ReadcphiqATq: if(sum(abs(q-q_bk))>tolq .OR. ispq/=ispq_bk)  then  
       associate(cphitemp=> readcphif(q,ispq))    
         cphiq(1:nlmto,1:ntq) = cphitemp(1:nlmto,itq(1:ntq))
       endassociate
       q_bk=q
       ispq_bk=ispq
    endif ReadcphiqATq
    ReadcphimATqk: if(sum(abs(qk-qk_bk))>tolq .OR. ispm/=ispm_bk) then
       cphim= readcphif(qk, ispm) !all readin but we need only bands nmini:nmmax
       qk_bk= qk
       ispm_bk= ispm
    endif ReadcphimATqk
    iasx=[(sum(nlnmv(iclass(1:ia-1)))+1,ia=1,natom)]
    icsx=[(sum(ncore(iclass(1:ia-1))),ia=1,natom)]
    allocate(zmelt(1:nbloch+ngc,nmtot,nqtot),source=(0d0,0d0))
    ZmelWithinMT: block !- Calculates <psi (k+q,t) |core(k,t) B(R,ibloch)>  also <core(k+q,t) |psi (k,t) B(R,ibloch)>
      !r B(R,i)   = Mixed basis. Bloch orthonormal product basis for atom R
      !r core(k,t)= core states
      !r <psi(k+q,t') | core(k,t) B(R,i)> = S[RLn]  cphik(RLn,k+q,t')^*  * ppb
      ! nt0=valence num of middle states at k! ntp0=valence num of end states at kq
      ! cphik, cphikq:  coefficients of MT part or valence eigenfunction. ! middle states at k, end state at k+q
      ! icore     = index for core states,         !i ncore      = no. core states in each class
      ! ppb      = <phi(RLn) phi(RL'n') B(R,i)>
      ! mdim      = number of optimal product basis functions for each class
      ! nbloch    = total no. optimal product basis
      ! nlnmx     = maximum number of l,n,m
      ! nctot      = total no. allowed core states
      ! nbloch     = total no. product basis within MT
      !           complex(8):: zppb(nlnmx,nlnmx),zz(nlnmx,ntp0)
      !           if(mdimx /= maxval(mdim) ) call rx( 'psi2b_v3: wrong mdimx')
      !           if(sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v3: wrong nbloch')
      !           if(sum(nlnmv(iclass(1:natom)))/=nlmto) call rx( ' psi2b_v3:sum(nlnmv)/= nlmto')
      !           if(sum(ncore(iclass(1:natom)))/= nctot) call rx( "psicb_v3:sum(ncore) wrong")
      !           if(ncc/=0 .AND. ncc/=nctot) call rx( "psicb_v3: ncc/=0 and ncc/=ncctot")
      associate(nmmax=>ntp0, mdim=>nblocha, cphik=>cphim(:,nmini:), cphikq=>cphiq(:,nqini:), zpsi2b =>zmelt(1:nbloch,:,:))
        zpsi2b = 0d0
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
           ime   = ims-1+mdim(icp)
           prodloop: do concurrent(i=1:mdim(icp)) !<psi(k+q,t') | psi(k,t) B(i)>=sum(Ln) bkq(Ln,t') * <phi(Ln) phi(L'n') B(i)>
              zpsi2b(i-1+ims,nctot+1:nctot+nt0,ncc+1:ncc+ntp0)= phase(ia) * matmul( transpose(cphik(ias:ias+nv-1,1:nt0)), &
                   dconjg( matmul(transpose(ppb(nc1:nc+nv,nc1:nc+nv,i,icp)),cphikq(ias:ias+nv-1,1:ntp0)) ) )
           enddo prodloop
           if(nctot==0) cycle
           !           do concurrent(i=1:mdim(icp)) !itp=1:ntp0, it=1:ncore(ic))
!           bpp =ppb(nc1:ncnv,icore(1:ncore(ic),ic),:,icp),order=[3,2,1]
           do concurrent(itp=1:ntp0, it=1:ncore(ic))
              zpsi2b(ims:ime,ics+it,ncc+itp)=phase(ia)*&
                   dconjg(matmul(cphikq(ias:iae,itp),ppb(nc1:ncnv,icore(it,ic),1:mdim(icp),icp)))
!              zpsi2b(i-1+ims:ime,ics+1:ics+ncore(ic),ncc+1:ncc+ntp0)=phase(ia)*&
!                   dconjg(matmul(cphikq(ias:iae,ncc+1:ncc+ntp0),ppb(nc1:ncnv,icore(1:ncore(ic),ic),i,icp)))
           enddo
           if(ncc>0) cycle !note ncc=nctot or 0
           do concurrent(itp=1:ncore(ic), it=1:nt0) 
              zpsi2b(ims:ime,nctot+it,ics+itp)=dconjg(phase(ia))&
                   *matmul(cphik(ias:iae,it),ppb(nc1:ncnv,icore(itp,ic),1:mdim(icp),icp))
           enddo
        enddo iatomloop
      endassociate
    endblock ZmelWithinMT
    if(ngc==0)goto 9898
    ZmelIPW:block ! calculates <psi(k',t') | psi(k,t) B(R,i)> for all R ! Valence x Valence = <IPW phi_valence |phi_valence>
      ispqk=ispm
      call readqg('QGpsi',q,     qt, ngp1, ngvecpB1) !q is mapped to qt in BZ
      call readqg('QGpsi',qk,   qkt, ngp2, ngvecpB2)
      geig1= readgeigf(q, ispq)   !read IPW part at q
      geig2= readgeigf(qk,ispqk)  !read IPW part at qk
      qdiff = matmul(symope,kvec) - qt + qkt ! rkvec + qkt - qt is not zero. <M(rkvec) Phi(qk) |Phi(q)>
      nadd = nint(matmul(transpose(plat),  qdiff)) !nadd: difference in the unit of reciprocal lattice vectors.
      !! zmelp(igc(qi),it(q2),itp(q1)) = <itp(for q1+G1)| it(for q2+G2) igc>
      !! NOTE: shtv = g(delta_{g^-1})
      !!--- zmelp(igc,it,itp) = <itp(for G1)|it(for G2) igc> matrix element.
      !!   zmelp0(igc,it,itp) = <G1|G2 Gc'> geig^*(G1,itp) geig(G2,it)
      !!   zmelp(igc,it,itp) =   = zmelp0(Gc',it,itp) <Gc'|Gc>^-1
      !!   (<Gc'|Gc>^-1 is dconjg(ppovlinv)
      !! New ggg matrix <Gc |G1 G2> is introduced.
      !!
      !!    <Gc G2|G1> is equivalent to <-Gc+G1-G2>; described by ggg
      !! Readin input
      !!    ggg(1:nggg) = <Gc+G2-G1>
      !!    nvggg(3,1:nggg)   for Gc+G2-G1
      !!    nvgcgp2(3,ngcgp2) for Gc+G2
      !!    ppovlinv(ngc,ngc) <Gc|Gc> matrix
      !!
      !! ggitp(Gc+G2)= \sum_G1 <(Symope(Gc)+G2)-G1> geigq1(G1,itp)*exp(-i*G1*shtv)*exp(-i(q-Gadd)*shtv)
      !! NOTE: nvgcgp2(:,igcgp2) means symope(Gc)+ G2
      melpln2t: block !> Mattrix elements <Plane psi |psi> from interstitial plane wave.
        use m_read_ppovl,only: getppx2,igggi,igcgp2i,nxi,nxe,nyi,nye,nzi,nze,nvggg,nvgcgp2,ngvecc
        use m_read_ppovl,only: nggg,ngcgp,ngcread, ggg,ppovlinv, ngc2,ngvecc2
        integer:: ngveccR(1:3,1:ngc),igcgp2,nn(3),iggg,igp1,itp,igc,igp2
        complex(8),allocatable::ggitp(:,:) 
        complex(8):: zmelp0(ngc,nt0,ntp0),ggitp_(ngp2,ngc),phase(ngc)!,ggitp(ntp0,ngcgp)
        call getppx2(qlat,kvec) ! read and allocate PPOVL*
        if(ngc/=ngcread) call rx( 'melpln2t: ngc/= ngcx by getppx:PPOVLG')!  write(6,*)qt,ngcread,ngc
        call rotgvec(symope, 1, ngc, [ngc], qlat, ngvecc, ngveccR)
        allocate(ggitp(ntp0,ngcgp),source=(0d0,0d0))
        do concurrent (igcgp2=1:ngcgp) !for ngc+ngp2
           do igp1  = 1,ngp1   
              nn = ngvecpB1(:,igp1)- nvgcgp2(:,igcgp2) - nadd ! G1 -(Gc+G2) - Gadd ! -Gadd= -rk + qt -qkt
              if(nn(1)<nxi .OR. nxe<nn(1) .OR. nn(2)<nyi .OR. nye<nn(2) .OR. nn(3)<nzi .OR. nze<nn(3)) cycle
              iggg = igggi(nn(1),nn(2),nn(3)) !inversion table
              if(iggg<0) cycle ! ggg(iggg) = <qt+G1 -(rk+Gc) -(qk+G2) >, where
              ggitp(:,igcgp2)=ggitp(:,igcgp2)+ggg(iggg)*geig1(igp1,itq(nqini:nqmax)) !time-consuing
           enddo
        enddo
        ggitp = dconjg(ggitp)
        ! zmelp <=  \sum_G2 ggitp(Gc+G2) geigqg2(G2)) !! note \bfr'= g (\bfr) +\delta_g  (\bfr= {\bf r})
        ! mapping of function g[f(\bfr)]= f(g^-1(\bfr)+\delta_{g^-1})
        ! zmelp0(igc'(Gc'),it(G2),itp(G1)) = <G1|G2 Gc'> geig*(G1,itp)geig(G2,it) = <itp(G1)|it(G2) Gc'>
        !zmelp0=0d0
        phase(:)=[(exp( img*2d0*pi*sum((matmul(symope,kvec)+matmul(qlat,ngveccR(:,igc)))*shtv) ),igc=1,ngc)]
        do concurrent (itp= 1:ntp0)
           do concurrent (igc=1:ngc,igp2=1:ngp2)
              nn = ngveccR(:,igc) + ngvecpB2(:,igp2)
              igcgp2 = igcgp2i(nn(1),nn(2),nn(3))
              ggitp_(igp2,igc) = phase(igc)*ggitp(itp,igcgp2)
           enddo
           zmelp0(:,:,itp)= matmul(transpose(ggitp_),geig2(1:ngp2,nmini:nmmax))
        enddo
        call matm(dconjg(ppovlinv),zmelp0,zmelt(nbloch+1:nbloch+ngc,nctot+1:nctot+nt0,ncc+1:ncc+ntp0),ngc,ngc,ntp0*nt0)
        deallocate(ggitp)
      endblock melpln2t
    endblock ZmelIPW
9898 continue
    nbb=ngb
    if(nbbx/=0) nbb=nbbx
    if(allocated(zmel)) deallocate(zmel)
    MultiplePPOVLZ: block
      allocate(zmel(nbb,ns1:ns2, nqtot))
      if(nt0==0) then ! core only for middle states. fast mode 
         zmel=0d0
         do concurrent(ia=1:natom)
            ics = sum(ncore(iclass(1:ia-1)))
            iap = iatomp(ia)
            icp = iclass(iap)
            do concurrent(it=1:ncore(iclass(ia)))
               if(ics+it<ns1.or.ns2<ics+it) cycle
               associate(nb1 => imdim(iap),nb2 => imdim(iap)+nblocha(icp)-1)
                 zmel(1:nbb,ics+it,1:nqtot) = dconjg(matmul(transpose(ppovlz(nb1:nb2,1:nbb)),zmelt(nb1:nb2,ics+it,1:nqtot)))
               endassociate
            enddo
         enddo
      else
         call matm(dconjg(transpose(ppovlz)), dconjg(zmelt), zmel,nbb,ngb,nmtot*nqtot)
      endif
    endblock MultiplePPOVLZ
    deallocate(zmelt)
    if(present(zmelconjg)) zmel=dconjg(zmel)
  end subroutine get_zmel_init
end module m_zmel

!  subroutine GramSchmidt_zmel() !oct15 2021 !experimenal
!    integer:: igb=1,it,itt
!    complex(8):: ov(nmtot),vec(nqtot),dnorm2(nmtot)
!    real(8):: dnorm
!    do it = 1,nmtot        ! occ
!       vec(:)= zmel(igb,it,:)
!       do itt = 1,it-1
!          ov(itt) = sum( dconjg(zmel(igb,itt,:))*vec(:))/dnorm2(itt)
!       enddo
!       vec = vec - matmul(ov(1:it-1),zmel(igb,1:it-1,:))
!       dnorm2(it) = sum(dconjg(vec)*vec)
!       zmel(igb,it,:) = vec
!    enddo
!  end subroutine GramSchmidt_zmel
