subroutine rsibl(ssite,sspec,lfrce,nbas,isp,q,iq,ndimh,nspc,& 
  napw,igapw,iprmb,nevec,evec,ewgt,k1,k2,k3,smpot,smrho,f)
  use m_struc_def  !Cgetarg
  use m_lmfinit,only: lat_alat,nspec
  use m_lattic,only: lat_qlat
  use m_lattic,only: lat_vol
  use m_supot,only: lat_nabc
  use m_supot,only: lat_ng
  use m_supot,only: lat_gmax
  use m_uspecb,only:uspecb
  use m_lattic,only:lat_plat

  !- Add smooth part of output density into smrho and forces.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   lfrce :if nonzero, accumulate contribution to force
  !i   nbas  :size of basis
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i     Stored:    *
  !i     Passed to: rsibl1
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: ngcut
  !i     Stored:    *
  !i     Passed to: tbhsi rsibl1 uspecb
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: alat plat qlat gmax nabc ng ogv okv vol
  !i     Stored:    *
  !i     Passed to: *
  !i   lfrce :1 calculate contribution to forces
  !i   nbas  :size of basis
  !i   q     :Bloch vector
  !i   iq    :index to current k-point
  !i   ndimh :dimension of hamiltonian
  !i   nspc  :2 for coupled spins; otherwise 1
  !i   napw  :number of augmented PWs in basis
  !i   igapw :vector of APWs, in units of reciprocal lattice vectors
  !i   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
  !i   nevec :number of eigenvectors with nonzero weights
  !i   evec  :eigenvectors
  !i   ewgt  :eigenvector weights
  !i   k1..3 :dimensions smpot,smrho
  !i   smpot :smooth potential on uniform mesh, needed for forces
  !o Outputs
  !o   smrho :smooth density accumulated for this qp
  !o   f     :force contribution accumulated for this qp
  !r Remarks
  !m MPI
  !m   Parallelise over the eigenvector loop. The vector block size is
  !m   chosen (in the range 6-16, by dstrbp.f) so as to distribute the
  !m   work optimally across processes. Two work arrays of the size of
  !m   smrho are allocated from the heap as buffers. Only one will be
  !m   needed under MPI-2. See comments in hsibl.
  !b Bugs
  !b    replace call to gvgvcomp and pass ipv as input
  !b    The non-F90 version should work, but it is no longer tested
  !u Updates
  !u   29 Dec 08 Unsuccessful attempt to make work with openmp
  !u   05 Jul 08 (T. Kotani) output density for new PW part
  !u   10 Sep 06 Added MPI parallelization in the spin-coupled case
  !u   23 Dec 04 Extended to spin-coupled case
  !u   25 Aug 04 Adapted to extended local orbitals
  !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
  !u   15 Feb 02 (ATP) Added MPI parallelization
  !u   27 Aug 01 Extended to local orbitals.
  !u   12 Oct 00 Use q-dependent list of G vectors
  !u    6 Jul 00 attempt to vectorize by grouping eigenvectors in blocks
  !u   17 Jun 00 Spin polarized
  !u   23 May 00 Adapted from nfp rsif_q.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: procid, master, nproc, mpipid
  integer :: lfrce,isp,k1,k2,k3,ndimh,nevec,iprmb(1),iq,nspc
  integer :: napw,igapw(3,napw),nbas
  real(8):: q(3) , ewgt(nevec) , f(3,nbas)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  double complex evec(ndimh,nspc,nevec),smrho(k1,k2,k3,isp), &
       smpot(k1,k2,k3,isp)
  ! ... Local parameters
  integer :: n0,nkap0,nermx,npmx,nblk,nlmto
  parameter (n0=10,nkap0=3,nermx=100,npmx=128)
  integer:: ngabc(3) , n1 , n2 , n3 , nrt , net , ng , &
       nglob , ltop , nlmtop , ogq , og2 , ohe , ohr , oyl , oylw , &
       oiv , iprint
  integer,allocatable :: iv_a_okv(:)
  real(8),allocatable :: rv_a_ogv(:)

  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
  integer :: iprt(n0,nkap0,nermx),ipet(n0,nkap0,nermx),i_copy_size
  double precision :: alat,qlat(3,3),plat(3,3),q0(3),gmax,xx,w
  double precision :: vol
  double precision :: etab(nermx),rtab(nermx)
  integer :: ivec,nvec
  integer,allocatable:: ivp(:)
  complex(8),allocatable::psi(:,:,:),psir(:,:,:),vpsi(:,:,:), &
       wk(:,:,:)
  real(8),allocatable:: cosi(:),sini(:),wk2(:)
  integer:: ivecini,ivecend
  integer,allocatable:: w_oiv(:)
  real(8),allocatable:: w_ogq(:),w_oyl(:),w_oylw(:),w_og2(:),w_ohe(:),w_ohr(:)
  complex(8),allocatable:: w_osmbuf(:)
  real(8),allocatable:: w_ofrbuf(:)
  nproc  = mpipid(0)
  procid = mpipid(1)
  if (nevec <= 0) return
  call tcn('rsibl')
  nlmto = ndimh-napw
  ! ... First setup
  alat=lat_alat
  plat = lat_plat
  qlat = lat_qlat
  gmax=lat_gmax
  ngabc=lat_nabc
  ng=lat_ng
  vol=lat_vol
  ! ... Setup for q-dependent gv ... also makes kv, gv+q and iv
  !     NB: gv generated by gvlst2 has q already added to it!
  call tcn('gvlst2')
  call pshpr(iprint()-30)
  call gvlst2(alat,plat,q,n1,n2,n3,0d0,gmax,0,500,0,ng,xx,xx,xx,xx)
  allocate(rv_a_ogv(abs(ng*3)))
  rv_a_ogv(:)=0.0d0
  allocate(iv_a_okv(abs(ng*3)))
  iv_a_okv(:)=0
  allocate(w_oiv(ng*3))
  call gvlst2(alat,plat, q, n1, n2, n3, 0d0,gmax,0,509, ng, ng, iv_a_okv, rv_a_ogv, w_oiv, w_oiv )
  call poppr
  call tcx('gvlst2')
  !     For PW basis ... for now.
  if (napw > 0) then
     allocate(ivp(napw))
     call gvgvcomp(ng,w_oiv,napw,igapw,ivp)
  else
     allocate(ivp(1))
  endif
  ! --- Tables of energies, rsm, indices to them ---
  call tbhsi(sspec,nspec,nermx,net,etab,ipet,nrt,rtab,iprt,ltop)
  ! --- Allocate and occupy arrays for yl, energy factors, rsm factors ---
  nlmtop = (ltop+1)**2
  allocate(w_ogq(ng*3), w_oyl(ng*nlmtop), w_oylw(ng*nlmtop), w_og2(ng), w_ohe(ng*net), w_ohr(ng*nrt))

  ! ino H_L(G)= \frac{-4 pi}{e-G^2} {cal Y}_L(-iG) exp(gamma(e-G^2))
  ! ino hsibl1 calculaets he=1/(e-G^2) and hr=exp(-gamma G^2)
  ! ino the other parts are calculated in rsibl5.
  call dpzero(q0,3)
  if (nlmto > 0) then
     call hsibl1 ( net , etab , nrt , rtab , ltop , alat , q0 , ng &
          , rv_a_ogv , w_ogq , w_og2 , w_oyl , w_ohe ,  w_ohr )
  endif
  deallocate(w_og2)

  nblk = 16
  !  --- Loop over eigenstates ---
  allocate(psi(ng,nspc,nblk),vpsi(ng,nspc,nblk),wk(ng,nspc,nblk))
  allocate(psir(k1,k2,k3),cosi(ng),sini(ng),wk2(ng))
  ivecini= 1
  ivecend= nevec
  do  ivec = ivecini,ivecend, nblk !blocked calculation for future
     nvec = min(nblk, nevec-ivec+1)
     call rsibl1(0,ssite,sspec,q,nbas,iprmb,ng,w_ogq,w_oiv,n1,n2, &
          n3,qlat,cosi,sini,w_oyl,w_oylw,w_ohe,w_ohr,wk, &
          wk2,vol,iprt,ipet,etab,rtab,ndimh,nlmto,nspc, &
          ewgt,ivec,nvec,evec,w,psi,w)
     ! ino    rsiblp adds PW(G) to psi
     call rsiblp(ng,ndimh,nlmto,nspc,napw,ivp,nvec,dsqrt(vol), &
          evec(1,1,ivec),psi)
     ! ino now psi= H(G) + PW(G)
     !   ... Add to real-space mesh, optionally make smpot*psi for forces
     ! ino rsibl2 executes FFT to get psi(r), which is F0
     ! ino and also calculates <psi|psi>(=F0F0) to get real space charge density.
     call rsibl2 ( ng , nspc , nvec , psi , n1 , n2 , n3 , k1 , k2 &
          , k3 , iv_a_okv ,  ewgt ( ivec ) , lfrce , smpot ( &
          1 , 1 , 1 , isp ) , psir , smrho ( 1 , 1 , 1 , isp ) , vpsi &
          )
     !    --- Add to forces ---
     if (lfrce /= 0) then
        call rsibl1(1,ssite,sspec,q,nbas,iprmb,ng,w_ogq,w_oiv,n1,n2, &
             n3,qlat,cosi,sini,w_oyl,w_oylw,w_ohe,w_ohr, &
             wk,wk2,vol,iprt,ipet,etab,rtab,ndimh,nlmto,nspc, &
             ewgt,ivec,nvec,evec,vpsi,psi,f)
     endif
  enddo
  deallocate(psi,vpsi,wk,psir,cosi,sini,wk2)
  if (allocated(rv_a_ogv)) deallocate(rv_a_ogv)
  if (allocated(iv_a_okv)) deallocate(iv_a_okv)
  deallocate(ivp)
  call tcx('rsibl')
end subroutine rsibl


subroutine rsibl1(mode,ssite,sspec,q,nbas,iprmb,ng,gq,iv,n1,n2,n3, &
     qlat,cosgp,singp,yl,ylw,he,hr,psi0,wk2,vol,iprt,ipet,etab,rtab, &
     ndimh,nlmto,nspc,ewgt,ivec,nvec,evec,vpsi,psi,f)
  use m_uspecb,only:uspecb
  use m_struc_def  !Cgetarg
  use m_orbl,only: Orblib,ktab,ltab,offl,norb
  !- Make wave function for a block of evecs, or add contr. to forces
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 make wave function
  !i         :1 Add 2*Re( (v psi+) grad(psi) ) to f
  !i   ssite :struct containing site-specific information
  !i   sspec :struct containing species-specific information
  !i   q     :Bloch wave number
  !i   nbas  :size of basis
  !i   ng    :number of G-vectors
  !i   gq    :2*pi/alat * (q+G) for all G-vectors
  !i   iv    :g-vectors as integer multiples of qlat (suphs0)
  !i   n1..3 :size uniform mesh for smooth density and potential
  !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
  !i   cosgp :cos(phase) for each g-vector
  !i   singp :sin(phase) for each g-vector
  !i   yl    :spherical harmonics for ng vectors
  !i   ylw   :work array of same dimension as yl
  !i   he    :table of energy factors
  !i   hr    :table of smoothing radius factors
  !i   psi0  :work array (dim ng*2*nspc*nev): psi sans phase factors
  !i   wk2   :work array of dimension ng
  !i   vol   :cell volume
  !o   iprt  :index to which entry in rt a given orbital belongs
  !i   ipet  :index to which entry in etab a given orbital belongs
  !i   etab  :table of all inequivalent energies
  !i   rtab  :table of all inequivalent smoothing radii
  !i   ndimh :dimensions evec
  !i   nspc  :2 for coupled spins; otherwise 1
  !i   ewgt  :weights for each of the trial fermi levels
  !i   ivec  :first of current block of eigenvectors
  !i   nvec  :number of eigenstates to generate
  !i   evec  :eigenvectors
  !i   vspi  :potential * wave function, needed only for mode=1
  !o Outputs
  !o   psi   :wave function (mode=0); work area (mode=1)
  !o   f     :term added to forces (mode=1)
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  real(8):: q(3)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)

  integer :: mode,nbas,ng,ndimh,nlmto,nspc,ivec,nvec,iv(ng,3), &
       n1,n2,n3,n0,nkap0,iprmb(*),i_copy_size
  parameter (n0=10,nkap0=3)
  integer :: iprt(n0,nkap0,*),ipet(n0,nkap0,*)
  double precision :: vol,yl(ng,*),ylw(ng,*),he(ng,*),hr(ng,*), &
       psi0(ng,2,nspc,nvec),wk2(ng),cosgp(ng),singp(ng),etab(*), &
       rtab(*),gq(ng,3),f(3,nbas),ewgt(nvec+ivec-1),qlat(3,3)
  double complex &
       psi(ng,nspc,nvec),evec(ndimh,nspc,ivec),vpsi(ng,nspc,nvec)
  ! ... Local parameters
  !      integer norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
  integer:: blks(n0*nkap0),ntab(n0*nkap0),ncut(n0,nkap0),lh(nkap0),nkapi
  double precision :: e,rsm,eh(n0,nkap0),rsmh(n0,nkap0),f0(3)
  double precision :: xx(n0),wt,p(3)
  integer :: ib,is,io,jo,l2,kp,ie,ir,ioff,nlm1,nlm2,iq,kb,lt,i
  ! takao
  integer::ncutt
  call dpzero(psi, 2*ng*nspc*nvec)

  if (nlmto == 0) return

  do  ib = 1, nbas

     is=ssite(ib)%spec
     i_copy_size=size(ssite(ib)%pos)
     call dcopy(i_copy_size,ssite(ib)%pos,1,p,1)


     i_copy_size=size(sspec(is)%ngcut)
     call icopy(i_copy_size,sspec(is)%ngcut,1,ncut,1)

     call suphas(q,p,ng,iv,n1,n2,n3,qlat,cosgp,singp)
     !       List of orbitals, their l- and k- indices, and ham offsets
     call orblib(ib) !,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
     !       Block into groups with consecutive l and common (e,rsm)
     call uspecb(is,rsmh,eh)
     call gtbsl1(7+16,norb,ltab,ktab,rsmh,eh,ntab,blks)

     call dpzero(psi0,ng*2*nspc*nvec)
     if (mode == 1) call dpzero(psi, 2*ng*nspc*nvec)
     do  io = 1, norb
        if (blks(io) /= 0) then
           jo = ntab(io)
           l2 = ltab(io)
           lt = ltab(jo)
           kp = ktab(io)
           ie = ipet(l2+1,kp,is)
           ir = iprt(l2+1,kp,is)
           ioff = offl(io)
           nlm1 = l2**2+1
           nlm2 = nlm1 + blks(io)-1
           rsm = rtab(ir)
           e   = etab(ie)
           ! takao Apr2009
           ncutt=ncut(lt+1,kp)
           call ncutcorrect(ncutt,1,gq,ng)

           ! ino    rsibl5 calculates 4 pi exp(e gamma) and {cal Y}_L(-iG)
           ! ino    and make psi0(G), which is H(-iG), with hr and he.
           call rsibl5(ie,ir,e,rsm,vol,nlm1,nlm2,ng,min(ng,ncutt) &
                ,yl,ylw,he,hr,wk2,ioff,evec(1,1,ivec),ndimh,nspc,nvec,psi0)
        endif
     enddo
     ! ino    multiply exp(i G R_i) * psi0 to make psi()
     call rsibl6(ng,nspc,nvec,cosgp,singp,psi0,psi)

     if (mode == 1) then
        do  i = 1, nvec
           call rsibl4(vol,ng,nspc,gq,vpsi(1,1,i),psi(1,1,i),f0)
           wt = ewgt(i+ivec-1)
           f(1,ib) = f(1,ib) + wt*f0(1)
           f(2,ib) = f(2,ib) + wt*f0(2)
           f(3,ib) = f(3,ib) + wt*f0(3)
           !             This shouldn't be necessary
           do  kb = 1, nbas
              f(1,kb) = f(1,kb) - wt*f0(1)/nbas
              f(2,kb) = f(2,kb) - wt*f0(2)/nbas
              f(3,kb) = f(3,kb) - wt*f0(3)/nbas
           enddo

        enddo
     endif
  enddo
end subroutine rsibl1


subroutine rsibl2(ng,nspc,nev,psi,n1,n2,n3,k1,k2,k3,kv,ewgt, &
     lfrce,smpot,f,smrho,vpsi)

  !- FT wave function to real space and add square into mesh density
  !  and optionally make smpot * psi
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ng    :number of G-vectors
  !i   nspc  :2 for coupled spins; otherwise 1
  !i   nev   :number of wave functions
  !i   psi   :wave function in reciprocal space
  !i   n1..3 :size of FT mesh
  !i   k1..3 :dimensions smpot,smrho
  !i   kv    :indices for gather/scatter operations (gvlist.f)
  !i   ewgt  :weights for each of the trial fermi levels
  !i   lfrce :if nonzero, make vpsi
  !i   smpot :potential, needed if lfrce is nonzero
  !o Outputs
  !o   f     :psi in real space
  !o   smrho :ewgt (f+)f added to smooth density
  !o   vpsi  :FT (smpot * r.s. wave function) if lfrce is nonzero
  !r Remarks
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: k1,k2,k3,n1,n2,n3,ng,nspc,nev,kv(ng,3),lfrce
  double precision :: ewgt(nev)
  double complex psi(ng,nspc,nev),vpsi(ng,nspc,nev),f(k1,k2,k3)
  double complex smrho(k1,k2,k3,nspc),smpot(k1,k2,k3,nspc)
  ! ... Local parameters
  integer :: i1,i2,i3,iq,i,ispc
  double precision :: wgt1

  call tcn('rsibl2')
  do  ispc = 1, nspc
     do  i = 1, nev
        call gvputf(ng,1,kv,k1,k2,k3,psi(1,ispc,i),f)
        call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,1)
        !          do  iq = 1, numq
        wgt1 = ewgt(i)
        do  i3 = 1, n3
           do  i2 = 1, n2
              do  i1 = 1, n1
                 smrho(i1,i2,i3,ispc) = smrho(i1,i2,i3,ispc) + &
                      wgt1*dconjg(f(i1,i2,i3))*f(i1,i2,i3)
              enddo
           enddo
        enddo
        !          enddo

        if (lfrce /= 0) then
           do  i3 = 1, n3
              do  i2 = 1, n2
                 do  i1 = 1, n1
                    f(i1,i2,i3) = f(i1,i2,i3)*smpot(i1,i2,i3,ispc)
                 enddo
              enddo
           enddo
           call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,-1)
           call gvgetf(ng,1,kv,k1,k2,k3,f,vpsi(1,ispc,i))
        endif
     enddo
  enddo
  call tcx('rsibl2')
end subroutine rsibl2


!      subroutine rsibl3(ng,n1,n2,n3,k1,k2,k3,kv,smpot,f,vpsi)
!C- Make f*smpot and transform to reciprocal space
!      implicit none
!C ... Passed parameters
!      integer k1,k2,k3,n1,n2,n3,ng,kv(ng,3)
!      double complex vpsi(ng),f(k1,k2,k3),smpot(k1,k2,k3)
!C ... Local parameters
!      integer i1,i2,i3

!      call tcn('rsibl3')
!      do  i3 = 1, n3
!        do  i2 = 1, n2
!          do  i1 = 1, n1
!            f(i1,i2,i3) = f(i1,i2,i3)*smpot(i1,i2,i3)
!          enddo
!        enddo
!      enddo
!      call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,-1)
!      call gvgetf(ng,1,kv,k1,k2,k3,f,vpsi)
!      call tcx('rsibl3')
!      end

subroutine rsibl4(vol,ng,nspc,gq,vpsi,psi,f0)
  !- Force term 2*Re( (psi_nu+) vsm grad(psi_nu) )
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   vol   :cell volume
  !i   ng    :number of G-vectors
  !i   nspc  :2 for coupled spins; otherwise 1
  !i   gq    :2*pi/alat * (q+G) for all G-vectors
  !i   vpsi  :(psi vsm) in reciprocal space
  !i   psi   :portion of wave function associated with one site ib
  !o Outputs
  !o   f0    :2*Re( (psi_nu+) vsm grad(psi_ib,nu) ) put in f0
  !r Remarks
  !r   gradient operator is i*G
  !u Updates
  !u   23 Dec 04 Extended to noncollinear case
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ng,nspc
  double precision :: vol,gq(ng,3),f0(3)
  double precision :: vpsi(2,nspc,ng),psi(2,nspc,ng)
  ! ... Local parameters
  integer :: ispc,ig
  double precision :: sum1,sum2,sum3,xx
  !     double complex cc,ovl

  call tcn('rsibl4')
  sum1 = 0
  sum2 = 0
  sum3 = 0
  !     ovl = 0
  do  ispc = 1, nspc
     do  ig = 1, ng
        !        cc  = dcmplx(psi(1,ispc,ig),-psi(2,ispc,ig))
        !     .       *dcmplx(vpsi(1,ispc,ig),vpsi(2,ispc,ig)) * vol
        !        ovl = ovl + cc
        !     .      * dcmplx(0d0,1d0)*gq(ig,1)
        xx = vpsi(2,ispc,ig)*psi(1,ispc,ig) &
             - vpsi(1,ispc,ig)*psi(2,ispc,ig)
        sum1 = sum1 + gq(ig,1)*xx
        sum2 = sum2 + gq(ig,2)*xx
        sum3 = sum3 + gq(ig,3)*xx
     enddo
  enddo

  f0(1) = 2*vol*sum1
  f0(2) = 2*vol*sum2
  f0(3) = 2*vol*sum3

  !     print *, ovl

  call tcx('rsibl4')

end subroutine rsibl4


subroutine rsibl5(ie,ir,e,rsm,vol,nlm1,nlm2,ng,ncut,yl,ylw,he,hr, &
     wk,ioff,evec,ndimh,nspc,nvec,psi)

  !- Add contribution to wave function from one block of orbitals
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ie    :index to appropriate entry in energy factor table he
  !i   ir    :index to appropriate entry in sm radius factor table hr
  !i   e     :hankel energy
  !i   rsm   :smoothing radius
  !i   vol   :cell volume
  !i   nlm1  :starting orbital L for which to accumulate wave function
  !i   nlm2  :final orbital L for which to accumulate wave function
  !i   ng    :number of G-vectors
  !i   ncut  :G-cutoff for wave function
  !i   yl    :spherical harmonics for ng vectors
  !i   ylw   :work array dimensioned same as yl
  !i   he    :table of energy factors
  !i   hr    :table of smoothing radius factors
  !i   wk    :work array of dimension at least ncut
  !i   ioff  :offset to hamiltonian (eigenvector) for this orbital block
  !i   ndimh :dimension of hamiltonian
  !i   nspc  :2 for coupled spins; otherwise 1
  !i   evec  :eigenvectors
  !i   nvec  :number of eigenvectors
  !o Outputs
  !o   psi   :contribution to psi from this block accumulated
  !r Remarks
  !u Updates
  !u   23 Dec 04 Extended to noncollinear case
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ie,ioff,ir,ncut,ng,nlm1,nlm2,ndimh,nspc,nvec
  double precision :: e,rsm,vol,yl(ng,*),ylw(ng,*),he(ng,*),hr(ng,*), &
       wk(ncut)
  double complex evec(ndimh,nspc,nvec),psi(ng,nspc,nvec)
  ! ... Local parameters
  integer :: i,ii,ilm,l,ll,lmax,m,iv,nlmx,ispc
  parameter (nlmx=100)
  double complex cfac(nlmx),cc,evp(nlmx,nspc,nvec)
  double precision :: pi,fac
  parameter (pi=3.1415926535897931d0)

  if (nlm2 == 0) return
  call tcn('rsibl5')

  ! ... Phase and other factors
  lmax = ll(nlm2)
  ! ino -4 pi exp(e gamma)
  fac = -4d0*pi*dexp(e*rsm*rsm*0.25d0)/vol
  ! ino phase factor of {cal Y}_L(-iG)
  cc = (0d0,1d0)*fac
  ilm = 0
  do  l = 0, lmax
     cc = cc*(0d0,-1d0)
     do m = -l,l
        ilm = ilm+1
        cfac(ilm) = cc
     enddo
  enddo

  ! ... Combine G-dependent energy, rsm and YL factors
  do  i = 1, ncut
     wk(i) = he(i,ie)*hr(i,ir)
  enddo
  do  ilm = nlm1, nlm2
     do  i = 1, ncut
        ylw(i,ilm) = wk(i)*yl(i,ilm)
     enddo
  enddo

  ! ... Make vector evec*phase
  do  ispc = 1, nspc
     do  ilm = nlm1, nlm2
        ii = ilm-nlm1+ioff+1
        do  iv = 1, nvec
           !          cc = evec(ii,ispc,iv)*cfac(ilm)
           !          evpr(ilm,ispc,iv) = dble(cc)
           !          evpi(ilm,ispc,iv) = dimag(cc)
           evp(ilm,ispc,iv) = evec(ii,ispc,iv)*cfac(ilm)
        enddo
     enddo

     ! ... For each orbital and evecs 1..nvec, accumulate psi
     !     ii = ilm-nlm1+ioff+1
     do  ilm = nlm1, nlm2
        do  iv = 1, nvec
           do  i = 1, ncut
              psi(i,ispc,iv) = psi(i,ispc,iv)+ylw(i,ilm)*evp(ilm,ispc,iv)
           enddo
        enddo
     enddo
  enddo

  call tcx('rsibl5')
end subroutine rsibl5


subroutine rsibl6(ng,nspc,nvec,cosgp,singp,psi0,psi)

  !- Multiply by phase to make final FT of partial wave function
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ng    :number of G-vectors
  !i   nspc  :2 for coupled spins; otherwise 1
  !i   nvec  :number of eigenvectors
  !i   cosgp :cos(phase)
  !i   singp :sin(phase)
  !i   wr    :real part of psi, unscaled by phase
  !i   wi    :imaginary part of psi, unscaled by phase
  !o Outputs
  !o   psi   :(wr,si)*(cosgp,singp) is added into psi
  !r Remarks
  !u Updates
  !u   23 Dec 04 Extended to noncollinear case
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ng,nspc,nvec
  double precision :: cosgp(ng),singp(ng)
  double complex psi0(ng,nspc,nvec),psi(ng,nspc,nvec)
  ! ... Local parameters
  integer :: i,iv,ispc

  call tcn('rsibl6')

  ! ino phase factor exp(i G R_spc)
  do  iv = 1, nvec
     do  ispc = 1, nspc
        do  i = 1, ng
           psi(i,ispc,iv) = psi(i,ispc,iv) &
                + psi0(i,ispc,iv)*dcmplx(cosgp(i),singp(i))
        enddo
     enddo
  enddo

  call tcx('rsibl6')
end subroutine rsibl6


subroutine rsibl7(xsmrho,k1,k2,k3,smrho)

  !- Combine rho from separate parallel threads
  implicit none
  ! ... Passed parameters
  integer :: k1,k2,k3
  double complex smrho(k1,k2,k3),xsmrho(k1,k2,k3)
  ! ... Local parameters
  integer :: ik1,ik2,ik3,iq
  !      do  iq = 1, numq
  do  ik3 = 1, k3
     do  ik2 = 1, k2
        do  ik1 = 1, k1
           smrho(ik1,ik2,ik3) = smrho(ik1,ik2,ik3) + &
                xsmrho(ik1,ik2,ik3)
        enddo
     enddo
  enddo
  !      enddo

end subroutine rsibl7


subroutine gvgvcomp(ng,igv,napw,igapw,ivp)

  !- Find pointer ivp that maps igv to igapw.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ng    :number of G-vectors
  !o   igv   :list of reciprocal lattice vectors G (gvlist.f)
  !i   napw   :number of R.L.V for PW basis
  !i   igapw :reciprocal lattice vectors for PW basis.
  !o Outputs
  !o   ivp   :if ig = ivp(jg), igv(jg) and nvec(ig) are same vector
  !r Remarks
  !r  This routine should be be cleaned up and ivp
  !r  used by rest of program in place of igapw
  !u Updates
  !u   05 Jul 08 (T. Kotani) first created
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ng,igv(ng,3),napw,igapw(3,napw),ivp(napw)
  ! ... Local parameters
  integer :: jg,ig

  !     Redesign ... inefficient.
  ! FCPP#if F90
  do  ig = 1, napw
     do  jg = 1, ng
        if(sum(abs(igv(jg,:)-igapw(:,ig))) == 0) then
           ivp(ig) = jg
           goto 333
        endif
     enddo
     ivp(ig) = -9999
     call rx('gvgvcomp wrong 111! maybe enlarge GMAX or so')
333  continue
     if( sum(abs( igapw(:,ig)-igv(ivp(ig),:) )) /=0) &
          call rx('bug in gvgvcomp.  Cannot find ivp')
  enddo
  ! FCPP#endif
end subroutine gvgvcomp


subroutine rsiblp(ng,ndimh,nlmto,nspc,napw,ivp,nvec,sqv,evec,psi)

  !- Plane wave part of evec
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ng    :number of G-vectors
  !i   nvec  :number of eigenstates to generate
  !i   evec  :eigenvectors
  !i   vspi  :potential * wave function, needed only for mode=1
  !i   sqv   :square root of volume
  !o Outputs
  !o   psi   :wave function
  !r Remarks
  !u Updates
  !u   05 Jul 08 (T. Kotani) first created
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: ng,ndimh,nlmto,nspc,nvec
  integer :: napw,ivp(napw)
  double precision :: sqv
  double complex psi(ng,nspc,nvec),evec(ndimh,nspc,nvec)
  ! ... Local parameters
  integer :: i,ispc,igv

  if (napw <= 0) return
  !     sqvol = dsqrt()
  do  ispc = 1, nspc
     do  i = 1, nvec
        do  igv = 1, napw
           psi(ivp(igv),ispc,i) = psi(ivp(igv),ispc,i) &
                + evec(nlmto+igv,ispc,i)/sqv
        enddo
     enddo
  enddo

end subroutine rsiblp

