module m_hsibl
  public hsibl,hsibl1
  private
contains
subroutine hsibl(ssite,sspec,k1,k2,k3,vsm,isp,q,ndimh,napw,igapw,h)
  use m_struc_def
  use m_lmfinit,only: lat_alat,nspec,nbas
  use m_lattic,only: lat_qlat
  use m_lattic,only: lat_vol
  use m_supot,only: lat_nabc
  use m_supot,only: lat_ng
  use m_supot,only: lat_gmax
  use m_lattic,only:  lat_plat
  use m_uspecb,only:uspecb
  use m_orbl,only: Orblib1,Orblib2,ktab1,ltab1,offl1,norb1,ktab2,ltab2,offl2,norb2
  use m_ftox
  use m_sugcut,only: ngcut
  use m_lmfinit,only: nkaph,nl,ndimx
  !- Interstitial ME of smooth Bloch Hankels, smooth potential.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ssite :struct containing site-specific information
  !i     Elts read: spec pos
  !i     Passed to:
  !i   sspec :struct containing species-specific information
  !i     Elts read: ngcut
  !i     Passed to: tbhsi uspecb
  !i   slat  :struct containing information about the lattice
  !i     Elts read: alat plat qlat gmax nabc ng vol
  !i   k1,k2,k3 dimensions of vsm
  !i   vsm   :smoothed potential, real-space mesh
  !i   isp   :current spin channel (1 or 2)
  !i   q     :Bloch vector
  !i   ndimh :dimension of hamiltonian
  !i   napw  :number of augmented PWs in basis
  !i   igapw :vector of APWs, in units of reciprocal lattice vectors
  !o Outputs
  !o   h     :interstitial matrix elements of vsm added to h
  !r Remarks
  !r  *How orbital is extracted and employed.
  !r   See Remarks in smhsbl.f
  !m MPI
  !m   Parallelise over the main loop over nbas. In the serial code, h
  !m   is added to in each pass through the loop. Furthermore h is non
  !m   zero on entry to hsibl. This leads to a problem for MPI because
  !m   each process cannot simply add to the existing array h and pool
  !m   results with ALLREDUCE: this would lead to double counting. Instead
  !m   the increment to h from each call to hsibl must be pooled after
  !m   the main loop and then added to h. This requires allocating a
  !m   workspace of the same dimension as h. A second workspace of the
  !m   same length is needed as a buffer for ALLREDUCE. This second
  !m   work array can be dispensed with once MPI-2 is implemented with
  !m   the MPI_IN_PLACE feature of ALLREDUCE. Because these two work
  !m   arrays are large, they are taken from the heap rather than
  !m   ALLOCATEd using F90. Note that allocation of one work array the
  !m   size of h from the heap does not increase memory load because the
  !m   workspace for the eigenvectors is not yet allocated.
  !u Updates
  !u   05 Jul 08 (T. Kotani) new APW part of basis
  !u   12 Aug 04 First implementation of extended local orbitals
  !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
  !u   15 Feb 02 (ATP) Added MPI parallelization
  !u   27 Aug 01 Extended to local orbitals.
  !u             At present, must compile with SGI with local orbitals!
  !u   12 Oct 00 Use q-dependent list of G vectors
  !u   22 May 00 Adapted from nfp hsi_q.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: k1,k2,k3,isp,ndimh,napw,igapw(3,napw)
  real(8):: q(3)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  double complex h(ndimh,ndimh),vsm(k1,k2,k3,isp)
  integer :: n0,npmx,nkap0,nkape,nermx,ngabc(3),nlmto
  parameter (n0=10,nkap0=3,nkape=2,nermx=100)
  integer:: ltop , n1 , n2 , n3 , net , ng , nglob , nlmtop &
       , nrt , iprint
  real(8) ,allocatable :: g_rv(:)
  real(8) ,allocatable :: g2_rv(:)
  real(8) ,allocatable :: gv_rv(:)
  real(8) ,allocatable :: he_rv(:)
  real(8) ,allocatable :: hr_rv(:)
  integer ,allocatable :: kv_iv(:)
  real(8) ,allocatable :: yl_rv(:)
  double precision :: alat,plat(3,3),qlat(3,3),vol,gmax,q0(3)
  double precision :: etab(nermx),rtab(nermx)
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
  integer:: i , ib1 , ib2 , ie , ofh1 , ofh2 , ip , ir , is1 , &
       is2 , j , mp , nlm1 , nlm2 , ik1 , ik2 , l1 , iorb1 , l2 , l2t &
       , iorb2 , jorb2 , osin1 , osin2 , oc1 , ocf1 , oc2 , ocf2 , ocos1 &
       , ocos2 , of , owk , ndim1 , ndim2 , nkap1 , nkap2
  integer ,allocatable :: iv_iv(:),ncuti(:)
  complex(8) ,allocatable :: wk2_zv(:)
  integer :: iprt(n0,nkap0,nermx),ipet(n0,nkap0,nermx), ncut(n0,nkap0)!,lh(nkap0)
  integer:: blks1(n0*nkap0),ntab1(n0*nkap0)
  integer:: blks2(n0*nkap0),ntab2(n0*nkap0)
  double precision :: eh1(n0,nkap0),rsmh1(n0,nkap0)
  double precision :: eh2(n0,nkap0),rsmh2(n0,nkap0)
  double precision :: xx(n0),p1(3),p2(3)
  integer:: xxxx(nkap0)
  integer :: ig1,i1,ig2,i2,igx(3),igx1,igx2,igx3,oiv1, iloop
  complex(8) ,allocatable :: h_zv(:)
  integer:: nnn
  complex(8),allocatable:: w_oc1(:,:),w_ocf1(:,:),w_oc2(:,:),w_ocf2(:,:),w_of(:,:)
  real(8),allocatable:: w_ocos1(:,:), w_osin1(:,:),w_ocos2(:,:), w_osin2(:,:),w_owk(:,:)
  real(8):: gmin=0d0
  integer:: ibini,ibend
  integer:: nnrl,lmri,li,nnrlx,nnrli,ik,ib,ndim
  call tcn('hsibl')
  nlmto = ndimh - napw
  if (nspec > nermx) call rx('hsibl: increase nermx')
  alat=lat_alat
  plat=lat_plat
  qlat=lat_qlat
  gmax=lat_gmax
  ngabc=lat_nabc
  ng=lat_ng
  vol=lat_vol
  ! --- <MTO | V | MTO and < MTO | V PW> parts of h ---
  if (nlmto > 0) then
     ! ... Setup for q-dependent gv ... also makes kv, gv+q and iv
     !     NB: gv generated by gvlst2 has q already added to it!
     call tcn('gvlst2')
     call pshpr(iprint()-30)
     ! OTE difference of arguments between gvlist and gvlst2. gmin and mshlst
     !        call gvlist(alat,plat,q,n1,n2,n3,gmax,500,0,ng,xx,xx,xx,xx)
     call gvlst2(alat,plat,q,n1,n2,n3,gmin,gmax,0,500,0,ng,xx,xx,xx,xx)
     allocate(gv_rv(ng*3))
     allocate(kv_iv(ng*3))
     allocate(iv_iv(ng*3))
     call gvlst2(alat, plat, q, n1, n2, n3, gmin, gmax, 0, 509, ng, ng, kv_iv, gv_rv, iv_iv, iv_iv)
     call poppr
     call tcx('gvlst2')
     ! ... Tables of energies, rsm, indices to them
     call tbhsi(sspec,nspec,nermx,net,etab,ipet,nrt,rtab,iprt,ltop)
     !     ndimx = maximum hamiltonian dimension for any site (in m_lmfinit now)
     nlmtop = (ltop+1)**2
     allocate(g_rv(ng*3))
     allocate(yl_rv(ng*nlmtop))
     allocate(g2_rv(ng))
     allocate(he_rv(ng*net))
     allocate(hr_rv(ng*nrt))
     call dpzero(q0,3)
     call hsibl1 ( net , etab , nrt , rtab , ltop , alat , q0 , ng &
          , gv_rv , g_rv , g2_rv , yl_rv , he_rv , hr_rv )
     mp = 1
     nnn=min(mp,nbas)
     allocate( w_oc1( ng*ndimx,nnn), w_ocf1(ndimx,nnn))
     allocate( w_oc2(nnn, ng*ndimx), w_ocf2(nnn,  ndimx), &
          w_ocos1( ng,nnn),  w_osin1( ng,nnn), &
          w_ocos2( ng,nnn),  w_osin2( ng,nnn), &
          w_owk  ( ng,nnn),  w_of(k1*k2*k3,nnn))
     ibini=1
     ibend=nbas
     do  iloop = ibini,ibend
        ib1=iloop
        ip = 1
        if (nbas < mp) ip = ib1
        ndim1 = 0
        is1=ssite(ib1)%spec
        p1=ssite(ib1)%pos
        call suphas ( q , p1 , ng , iv_iv , n1 , n2 , n3 , qlat , w_ocos1(1,ip) , w_osin1(1,ip) )
        call orblib1(ib1) !norb1,ltab1,ktab1,xx,offl1,xx)
        ofh1 = offl1(1)
        call uspecb(is1,rsmh1,eh1)!       Block routines into groups with common (e,rsm)
        call gtbsl1(7+16,norb1,ltab1,ktab1,rsmh1,eh1,ntab1,blks1) ![1,1,1,0,1]
        do  iorb1 = 1, norb1
           if (blks1(iorb1) == 0) cycle
           l1   = ltab1(iorb1)
           ik1  = ktab1(iorb1)
           nlm1 = l1**2+1
           nlm2 = nlm1 + blks1(iorb1)-1
           ie   = ipet(l1+1,ik1,is1)
           ir   = iprt(l1+1,ik1,is1)
           call hsibl3 ( ie , ir , etab , rtab , vol , nlm1 , nlm2 , ndim1 &
                , ng , yl_rv , he_rv , hr_rv , w_ocos1(1,ip) , w_osin1(1,ip) &
                , w_owk(1,ip) , w_oc1(1,ip) , w_ocf1(1,ip) )
           ndim1 = ndim1 + max(blks1(iorb1),0)
        enddo
        !   ... Multiply potential into wave functions for orbitals in ib1
        call hsibl4(n1,n2,n3,k1,k2,k3,vsm(1,1,1,isp),w_of(1,ip),ng,kv_iv,ndim1,w_oc1(1,ip) )
        !   ... Loop over second of (ib1,ib2) site pairs
        do 1010 ib2 = ib1, nbas
           is2 =ssite(ib2)%spec
           p2  =ssite(ib2)%pos
           ncut=ngcut(:,:,is2)
           call suphas (q,p2,ng, iv_iv , n1 , n2 , n3 , qlat , w_ocos2(1,ip) , w_osin2(1,ip) )
           call orblib2(ib2) !norb2,ltab2,ktab2,offl2
           ofh2 = offl2(1)
           call uspecb(is2,rsmh2,eh2) ! Block into groups with consecutive l and common (e,rsm)
           call gtbsl1(7+16,norb2,ltab2,ktab2,rsmh2,eh2,ntab2,blks2) ![1,1,1,0,1]
           ndim2 = 0
           do  iorb2 = 1, norb2
              if (blks2(iorb2) == 0) cycle
              nlm1 = l2**2+1
              nlm2 = nlm1 + blks2(iorb2)-1
              ndim2 = ndim2 + max(blks2(iorb2),0)
           enddo
           allocate(ncuti(ndim2))
           ndim2 = 0
           do  iorb2 = 1, norb2
              if (blks2(iorb2) == 0) cycle
              jorb2 = ntab2(iorb2)
              l2t  = ltab2(jorb2)
              l2   = ltab2(iorb2)
              ik2  = ktab2(iorb2)
              nlm1 = l2**2+1
              nlm2 = nlm1 + blks2(iorb2)-1
              ie   = ipet(l2+1,ik2,is2)
              ir   = iprt(l2+1,ik2,is2)
              call hsibl3 ( ie , ir , etab , rtab , vol , nlm1 , nlm2 , ndim2 &
                   , ng , yl_rv , he_rv , hr_rv , w_ocos2(1,ip) , w_osin2(1,ip) &
                   , w_owk(1,ip) , w_oc2(1,ip) , w_ocf2(1,ip) )
              ncuti(ndim2+1:ndim2+nlm2-nlm1+1)=ncut(l2t+1,ik2)
              ndim2 = ndim2 + max(nlm2-nlm1+1,0)
           enddo
           !     ... Scalar products phi1*vsm*phi2 for all orbitals in (ib1,ib2)
           allocate(wk2_zv(ndim1*ndim2))
           wk2_zv=0d0
           ! ncuti are only at Gamma point; thus symmetry can not be kept well for other k points.
           call ncutcorrect ( ncuti , ndim2 , gv_rv , ng )
           call hsibl2 ( ndim1 , ndim2 , ng , ncuti , w_oc1(1,ip) , w_ocf1(1,ip) &
                , w_oc2(1,ip) , w_ocf2(1,ip) , ndim1 , 0 , 0 , wk2_zv )
           deallocate(ncuti)
           call hsibl5 ( norb1 , blks1 , offl1 , ndim1 , norb2 , blks2 , &
                offl2 , ndim2 , ndimh , wk2_zv , h )
           deallocate(wk2_zv)
1010    enddo 
        !   ... Matrix elements <Hsm| Vsm |PW>
        call hsibl6 ( ndimh , nlmto , norb1 , blks1 , offl1 , ng , iv_iv &
             , napw , igapw , w_oc1(1,ip) , w_ocf1(1,ip) , h )
     enddo 
     deallocate(hr_rv, he_rv, g2_rv, yl_rv, g_rv, iv_iv, kv_iv, gv_rv, &
          w_oc1,w_ocf1, w_oc2, w_ocf2, w_ocos1,w_osin1,w_ocos2,w_osin2,w_owk,w_of)
  endif
  ! --- <e^i qpG | V |e^i qpG'>/vol = V(G'-G) ---
  if (napw > 0) then
     call fftz3(vsm(1,1,1,isp),n1,n2,n3,k1,k2,k3,1,0,-1)
     do  ig1 = 1, napw
        i1 = ig1+nlmto
        do  ig2 = ig1, napw
           i2 = ig2+nlmto
           igx = igapw(:,ig1) - igapw(:,ig2)
           igx1 = mod(igx(1)+10*n1,n1)
           igx2 = mod(igx(2)+10*n2,n2)
           igx3 = mod(igx(3)+10*n3,n3)
           if (igx1<0 .OR. igx2<0 .OR. igx3<0) call rx('igx<0')
           h(i1,i2) = h(i1,i2) + vsm(igx1+1,igx2+1,igx3+1,isp)
        enddo
     enddo
     call fftz3(vsm(1,1,1,isp),n1,n2,n3,k1,k2,k3,1,0,1)
  endif
  do  i = 1, ndimh
     do  j = i, ndimh
        h(j,i) = dconjg(h(i,j)) ! ... Occupy second half of matrix
     enddo
  enddo
  call tcx('hsibl')
end subroutine hsibl


subroutine hsibl1(net,et,nrt,rt,ltop,alat,q,ng,gv,g,g2,yl,he,hr)
  use m_ropyln,only: ropyln
  !- Make yl's, energy and rsm factors for list of G vectors
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   net   :size of table et
  !i   et    :table of all inequivalent energies
  !i   nrt   :size of table rt
  !i   rt    :table of all inequivalent smoothing radii
  !i   ltop  :largest l at any site
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   q     :Bloch wave number
  !i   ng    :number of G-vectors
  !i   gv    :list of reciprocal lattice vectors G (gvlist.f)
  !o Outputs
  !o   g     :2*pi/alat * (q+gv) for all g-vectors
  !o   g2    :g**2
  !o   yl    :Y_L
  !o   he    :1/(et-g2) for all inequivalent e's and g-vectors
  !o   hr    :dexp(-(rsm/2)**2*g2(i)) for all inequivalent rsm and g-vecs
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ltop,net,ng,nrt
  double precision :: alat,q(3),gv(ng,3),g(ng,3),yl(ng,1),he(ng,net), &
       hr(ng,nrt),g2(ng),et(net),rt(nrt)
  ! ... Local parameters
  integer :: i,ie,ir
  double precision :: pi,tpiba,gam

  ! ... Make (2*pi/alat)*(gv+q) in g
  pi = 4d0*datan(1d0)
  tpiba = 2d0*pi/alat
  do  i = 1, ng
     g(i,1) = tpiba*(gv(i,1)+q(1))
     g(i,2) = tpiba*(gv(i,2)+q(2))
     g(i,3) = tpiba*(gv(i,3)+q(3))
  enddo
  ! ... Make the yl's and g2
  call ropyln(ng,g(1,1),g(1,2),g(1,3),ltop,ng,yl,g2)
  ! ... Make the energy factors
  do  ie = 1, net
     do  i = 1, ng
        he(i,ie) = 1d0/(et(ie)-g2(i))
     enddo
  enddo
  ! ... Make the rsm factors
  do  ir = 1, nrt
     gam = 0.25d0*rt(ir)*rt(ir)
     do  i = 1, ng
        hr(i,ir) = dexp(-gam*g2(i))
     enddo
  enddo
end subroutine hsibl1

subroutine hsibl2(n1,n2,ng,ncut2,c1,cf1,c2,cf2,ndimh,ofh1,ofh2,h)
  !- Add scalar product (phi1 (vsm*phi2)) to  h
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n1    :number of c1 block of orbitals
  !i   n2    :number of c2 block of orbitals
  !i   ng    :leading dimension of c1,c2
  !i   ncut2 :number of G-vectors
  !i   c1    :product (phi1 vsm) in reciprocal space
  !i   c2    :phi2 in reciprocal space
  !i   ndimh :dimension of h
  !i   ofh1 :row offset to h for first orbital in block
  !i   ofh2 :col offset to h for first orbital in block
  !o Outputs
  !o   h     :scalar product (phi1 vsm) phi2 is added to h
  !r Remarks
  !r   This routine takes M^2 N operations, and is the most
  !r   time consuming step for moderate to large systems.
  !r   (M = number of orbitals, N = number of mesh points)
  !u Updates
  !u   22 May 00 Adapted from nfp pvhsiq
  ! ----------------------------------------------------------------------
  implicit none
  integer :: n1,n2,ofh1,ofh2,ncut2(n2),ndimh,ng
  double complex c1(ng,n1),c2(ng,n2),cf1(n1),cf2(n2),h(n1,n2) !h(ndimh,ndimh)
  integer :: i,i1,i2,ncut
  double complex csum
  call tcn('hsibl2')
  !$$$#if NBAR
  !$$$C     xx = 1
  !$$$C     do   i = 1, n2
  !$$$C       xx = xx*min(ng,ncut2(i))
  !$$$C     enddo
  !$$$C     nbar = xx**(1.d0/dble(n2))
  !$$$      nbar = min(ng,ncut2(1))
  !$$$      do  i = 1, n2
  !$$$        nbar = min(nbar,min(ng,ncut2(i)))
  !$$$      enddo
  !$$$      call zgemm('C','N',n1,n2,nbar,dcmplx(1d0,0d0),c1,ng,c2,ng,
  !$$$     .dcmplx(0d0,0d0),wk,n1)
  !$$$      do  i2 = 1, n2
  !$$$        ncut = min(ng,ncut2(i2))
  !$$$        do  i1 = 1, n1
  !$$$          csum = wk(i1,i2)
  !$$$          do  i = nbar+1, ncut
  !$$$            csum = csum + dconjg(c1(i,i1))*c2(i,i2)
  !$$$          enddo
  !$$$C         This reduces accuracy, but makes exactly compatible
  !$$$          do  i = ncut+1, nbar
  !$$$            csum = csum - dconjg(c1(i,i1))*c2(i,i2)
  !$$$          enddo
  !$$$          csum = csum * dconjg(cf1(i1))*cf2(i2)
  !$$$          h(i1+ofh1,i2+ofh2) = h(i1+ofh1,i2+ofh2) + csum
  !$$$        enddo
  !$$$      enddo
  !$$$#elif CRAY
  !$$$      do  i2 = 1, n2
  !$$$        ncut = min(ng,ncut2(i2))
  !$$$        do  i1 = 1, n1
  !$$$          wk(i1) = 0
  !$$$        enddo
  !$$$        do  i = 1, ncut
  !$$$          do  i1 = 1, n1
  !$$$            wk(i1) = wk(i1) + dconjg(c1(i,i1))*c2(i,i2)
  !$$$          enddo
  !$$$        enddo
  !$$$        do  i1 = 1, n1
  !$$$          csum = wk(i1) * dconjg(cf1(i1))*cf2(i2)
  !$$$          h(i1+ofh1,i2+ofh2) = h(i1+ofh1,i2+ofh2) + csum
  !$$$        enddo
  !$$$      enddo
  !$$$#else
  !      print *,' bbbb:',ndimh,ofh1+n1,ofh2+n2
  do  i2 = 1, n2
     ncut = min(ng,ncut2(i2))
     do  i1 = 1, n1
        csum = 0
        do  i = 1, ncut
           csum = csum + dconjg(c1(i,i1))*c2(i,i2)
        enddo
        csum = csum * dconjg(cf1(i1))*cf2(i2)
        h(i1+ofh1,i2+ofh2) = h(i1+ofh1,i2+ofh2) + csum
     enddo
  enddo
  !$$$#endif
  call tcx('hsibl2')
end subroutine hsibl2

subroutine hsibl3(ie,ir,etab,rtab,vol,nlm1,nlm2,offlm,ng,yl,he,hr, &
     cosgp,singp,wk,c,cfac)
  !- FT of smooth Hankels, without constant factors
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ie    :index to appropriate entry in energy factor table he
  !i   ir    :index to appropriate entry in sm radius factor table hr
  !i   etab  :table of all inequivalent energies
  !i   rtab  :table of all inequivalent smoothing radii
  !i   vol   :cell volume
  !i   nlm1  :starting L for which to accumulate FT
  !i   nlm2  :final L for which to accumulate FT
  !i   offlm :offset to ilm for storing c
  !i   ng    :number of G-vectors
  !i   yl    :spherical harmonics for ng vectors
  !i   he    :table of energy factors
  !i   hr    :table of smoothing radius factors
  !i   cosgp :table of phases (real part)
  !i   singp :table of phases (imaginary part)
  !i   wk    :real work array of length ng
  !o Outputs
  !o   c     :c(ilm+offlm) = phase he(ie) hr(ir) yl(ilm)
  !r Remarks
  !r   This routine requires M^2 N operations, and is one of the most
  !r   time consuming for large systems.
  !u Updates
  !u   22 May 00 Adapted from nfp su_hkft
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ie,ir,nlm1,nlm2,offlm,ng
  double precision :: yl(ng,1),he(ng,1),hr(ng,1),cosgp(1),singp(1), &
       wk(ng),etab(ie),rtab(ir),vol
  double complex c(ng,offlm+1+nlm2-nlm1),cfac(offlm+1+nlm2-nlm1) !c(ng,nlm2),cfac(nlm2)
  ! ... Local parameters
  integer :: i,ilm,offi,lmax,ll,l,m
  double precision :: xxx,fac1,pi
  double complex cf
  parameter (pi = 3.1415926535897931d0)

  if (nlm2 == 0) return
  call tcn('hsibl3')
  do  i = 1, ng
     wk(i) = he(i,ie)*hr(i,ir)
  enddo
  offi = offlm-nlm1+1
  do  ilm = nlm1, nlm2
     do  i = 1, ng
        xxx = wk(i)*yl(i,ilm)
        c(i,ilm+offi) = dcmplx(xxx*cosgp(i),xxx*singp(i))
     enddo
  enddo
  !     Constant factor
  fac1 = -4d0*pi*dexp(etab(ie)*rtab(ir)**2/4)/dsqrt(vol)
  lmax = ll(nlm2)
  ilm = 0
  cf = (0d0,1d0)
  do  l = 0, lmax
     cf = cf*(0d0,-1d0)
     do  m = -l, l
        ilm = ilm+1
        if (ilm >= nlm1) cfac(ilm+offi) = cf * fac1
     enddo
  enddo
  call tcx('hsibl3')
end subroutine hsibl3

subroutine hsibl4(n1,n2,n3,k1,k2,k3,vsm,f,ng,kv,nc,c)
  !- FFT to real space, multiply by potential, FTT back
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n1..3 :FT mesh
  !i   k1..3 :dimensions vsm,f,c
  !i   vsm   :potential
  !i   f     :work array to hold intermediate FFT
  !i   ng    :number of G-vectors
  !i   kv    :indices for gather/scatter operations (gvlist.f)
  !i   nf    :number of wave functions
  ! o Inputs/Outputs
  ! o  c     :On input, holds FT of wave function
  ! o        :On output, holds FT of wave function * potential
  !r Remarks
  !r   This routine takes M N logN operations, and is often the most
  !r   time consuming for small to moderate systems.
  !r   (M = number of orbitals, N = number of mesh points)
  !u Updates
  !u   22 May 00 Adapted from nfp shkpot
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: n1,n2,n3,k1,k2,k3,ng,nc,kv(ng,3)
  double complex c(ng,nc),f(k1,k2,k3),vsm(k1,k2,k3)
  ! ... Local parameters
  integer :: i,i1,i2,i3
  call tcn('hsibl4')
  do  i = 1, nc
     call gvputf(ng,1,kv,k1,k2,k3,c(1,i),f)
     call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,1)
     !       call zprm3('psir',0,f,n1,n2,n3)
     do  i3 = 1, n3
        do  i2 = 1, n2
           do  i1 = 1, n1
              f(i1,i2,i3) = f(i1,i2,i3)*vsm(i1,i2,i3)
           enddo
        enddo
     enddo
     !       call zprm3('v*psir',0,f,n1,n2,n3)
     call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,-1)
     call gvgetf(ng,1,kv,k1,k2,k3,f,c(1,i))
  enddo
  call tcx('hsibl4')
end subroutine hsibl4

subroutine hsibl5(norb1,blks1,offl1,ndim1,norb2,blks2,offl2,ndim2,ndimh,hwk,h)
  !- Adds a subblock of matrix elements into the hamiltonian
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   norb1 :number of row orbital blocks
  !i   blks1 :size of each row orbital block
  !i   offl1 :offset to start of each row orbital block
  !i   ndim1 :leading row dimension to hwk
  !i   norb2 :number of column orbital blocks
  !i   blks2 :size of each column orbital block
  !i   offl2 :offset to start of each column orbital block
  !i   ndim2 :leading column dimension to hwk
  !i   ndimh :dimensions hamiltonian
  !i   hwk   :matrix elements of this block to be added to h
  !o Outputs
  !o   h     :matrix elements added to hamiltonian for this block
  !l Local variables
  !l         :
  !r Remarks
  !r
  !u Updates
  !u   16 Aug 04 First created
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ndimh,ndim1,ndim2
  integer :: norb1,blks1(norb1),offl1(norb1)
  integer :: norb2,blks2(norb2),offl2(norb2)
  double complex h(ndimh,ndimh),hwk(ndim1,ndim2)
  ! ... Local parameters
  integer :: io1,nlm1,ofh1,i1,ofw1
  integer :: io2,nlm2,ofh2,i2,ofw2
  ofw1 = 0
  do  io1 = 1, norb1
     if (blks1(io1) /= 0) then
        ofh1 = offl1(io1)
        nlm1 = blks1(io1)
        ofw2 = 0
        do  io2 = 1, norb2
           if (blks2(io2) /= 0) then
              ofh2 = offl2(io2)
              nlm2 = blks2(io2)
              do  i1 = 1, nlm1
                 do  i2 = 1, nlm2
                    h(ofh1+i1,ofh2+i2) = h(ofh1+i1,ofh2+i2) + &
                         hwk(ofw1+i1,ofw2+i2)
                 enddo
              enddo
              ofw2 = ofw2 + blks2(io2)
           endif
        enddo
        ofw1 = ofw1 + blks1(io1)
     endif
  enddo
end subroutine hsibl5

subroutine hsibl6(ndimh,nlmto,norb1,blks1,offl1,ng,igv,napw,igapw,c1,cf1,h)
  !- Make matrix elements <Hsm | V | PW>
  ! ----------------------------------------------------------------------
  !i   ndimh :dimension of hamiltonian
  !i   nlmto :number of lmto's
  !i   ng    :number of G-vectors to represent lmto's
  !i   napw   :number of PW's
  !i   igv   :list of reciprocal lattice vectors G (from gvlist)
  !o Outputs
  !o   h     : <Hsm | V | PW> is added to h
  !l Local variables
  !l         :
  !r Remarks
  !b Bugs
  !b   ifindiv should be replaced with index passed through
  !u Updates
  !u   04 Jul 08 (T Kotani) first created
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ng,napw,ndimh,nlmto,igv(ng,3),igapw(3,napw)
  integer :: offl1(*),blks1(*)
  double complex h(ndimh,ndimh),c1(ng,nlmto),cf1(nlmto)
  ! ... Local parameters
  integer :: ig,i2,i2x,ofw1,io1,norb1,ofh1,nlm1,i1

  do  ig = 1, napw
     i2  = nlmto+ig
     !       index matching igv,igapw
     i2x = ifindiv(igapw(1,ig),igv,ng)
     ofw1 = 0
     do  io1 = 1, norb1
        if (blks1(io1) == 0) cycle
        ofh1 = offl1(io1)
        nlm1 = blks1(io1)
        do  i1 = 1, nlm1
           h(ofh1+i1,i2) = h(ofh1+i1,i2) &
                + dconjg( cf1(ofw1+i1)*c1(i2x, ofw1+i1) )
        enddo
        ofw1 = ofw1 + blks1(io1)
     enddo
  enddo
end subroutine hsibl6

integer function ifindiv(igapw,igv,ng)
  !- Find index in igv that corresponds to igapw
  ! ----------------------------------------------------------------------
  !i   igapw :vector of APWs, in units of reciprocal lattice vectors
  !i   igv   :List of G vectors
  !i   ng    :number of group operations
  !o Outputs
  !o   ifindiv:index to igv that matches igapw
  !u Updates
  !u   19 Jan 09 Original cleaned up, made more efficient
  ! ----------------------------------------------------------------------
  implicit none
  integer :: ng,igapw(3),igv(ng,3)
  integer :: ig
  ifindiv=-999999
  do  ig = 1, ng
     if (igapw(1) /= igv(ig,1)) cycle
     if (igapw(2) /= igv(ig,2)) cycle
     if (igapw(3) /= igv(ig,3)) cycle
     ifindiv = ig
     return
  enddo
  call rx('ifindiv: igapw not found in igv')
end function ifindiv

end module m_hsibl
