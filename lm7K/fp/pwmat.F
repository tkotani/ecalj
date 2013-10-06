      subroutine pwmat(slat,ssite,sspec,nbas,ndimh,napw,igapw,iprmb,q,
     .ngp,nlmax,igv,GcutH,inn,ppovl,pwhovl)
      use m_struc_def  !Cgetarg
      implicit none
!! inn:  takao added 27 Apr2009
!!   We have  q+G(igvx; internal in pwmat) = qp + G(igapw).
!!  (because qp(shortedned q)  is used when w(oigv2)=igapw is generated. igapw is used in  hambls).
!!   Thus, we need to uss
!!           igv(internally in pwmat) = igapw + inn, where inn = qlatinv*(qp-q).

C- Matrix elements (IPW,IPW) and (IPW,envelope function)
C ----------------------------------------------------------------------
Ci Inputs
Ci   slat  :struct containing information about the lattice
Ci     Elts read: alat plat qlat vol
Ci     Stored:
Ci   ssite :struct containing site-specific information
Ci     Elts read: spec
Ci     Stored:    pos spec
Ci   sspec :struct containing species-specific information
Ci     Elts read: rmt
Ci     Stored:
Ci   nbas  :size of basis
Ci   ndimh :dimension of hamiltonian
Ci   napw  :number of augmented PWs in basis
Ci   igapw :list of G vectors for APW part of basis as multiples of
Ci         :primitive reciprocal lattice vectors
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   q     :Bloch wave vector
Ci   nlmax :maximum value of (1+augmentation l)**2 within any sphere
Ci   ngp   :no. G vectors for eigenfunction expansion (depends on q)
Ci         :Cutoff specified by QpGcut_psi in GWinput
Ci   igv   :list of ngp G vectors for PW expansion of basis as multiples
Ci         :of primitive reciprocal lattice vectors.
Ci   GcutH :G cutoff for smoothed Hankel basis, LDA expansion
Co Outputs
Co   ppovl : <IPW_G1 | IPW_G2 >
Co   pwhovl: <IPW_G1 | basis function = smooth Hankel or APW >
Co         : matrix element; see Remarks
Cl Local variables
Cl   igvx  :same function as igv : a list of G-vectors that is used
Cl         :to expand eigenfunctions.  Difference is cutoff:
Cl         :igvx : cutoff is fixed by LMTO input HAM->GMAX
Cl         :igv  : cutoff is fixed by GWinput QpGcut_psi
Cl         :This inconsistency probably should be removed.
Cl   pwh   :coefficients to PW expansion of basis functions, and also
Cl         :therefore of IPW expansion of interstitial part of basis
Cl         :pwh(ig,j) = Fourier coff ig for envelope function j
Cr Remarks
Cr   IPW(G) denotes projected plane wave, i.e. PW with parts within MT
Cr   spheres projected out.
Cr   The interstitial part of the basis function phi_j is:
Cr      phi_j (istl) = sum_G  pwh(G,j) IPW(G)
Cr   Matrix elements between IPWs and basis functions are then
Cr      pwhovl(G,j) = int [IPW(G) sum_G1 pwh(G1,j) IPW(G1)]
Cr                = sum_G1 int [IPW(G) IPW(G1)] * pwh(G1,j) =
Cr                = sum_G1 ppovl(G,G1) pwh(G1,j) = ppovl * pwh
Cu Updates
Cu   29 Jan 09 Incorporate APW basis
Cu   25 Aug 04 Adapted to extended local orbitals
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   09 Apr 01 Adapted from Kotani's pplmat2
C ----------------------------------------------------------------------
      integer ngp,nlmax,igv(3,ngp),nbas,ndimh,iprmb(1),
     .napw,igapw(3,napw)
      real(8):: q(3) , GcutH
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex ppovl(ngp,ngp), pwhovl(ngp,ndimh)
C ... Local parameters
      integer ips(nbas),ib,is,igetss,ngmx,ig,lmxax,ll,iwk(3),nlmto,iga,
     .ifindiv2
      double precision alat,plat(3,3),qlat(3,3),vol,pi,pi4,tpiba,xx,
     .tripl,dgetss,bas(3,nbas),rmax(nbas),qpg(3),qpg2,denom,gam,srvol
      integer n0,nkap0
      parameter (n0=10, nkap0=3)
      integer lh(nkap0),nkapi
      integer norb,io,l,ik,offh,ilm,ltab(n0*nkap0),ktab(n0*nkap0),
     .offl(n0*nkap0),blks(n0*nkap0),ntab(n0*nkap0)

      double precision eh(n0,nkap0),rsmh(n0,nkap0)
      double complex phase,img,fach,mimgl(0:n0)
      integer,allocatable:: igvx(:,:),kv(:,:)
      double precision,allocatable:: yl(:)
      double complex,allocatable:: pwh(:,:),ppovlx(:,:)

      integer :: inn(3),matmul_pwhovl
      logical:: debug=.false.

      integer:: i_copy_size, i_spackv, i_spacks
C --- Setup ---

      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      i_copy_size=size(slat%qlat)
      call dcopy(i_copy_size,slat%qlat,1,qlat,1)
      vol=slat%vol

      pi  = 4d0*datan(1d0)
      pi4 = 4d0*pi
      tpiba = 2*pi/alat
      vol = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      srvol = dsqrt(vol)
      img = dcmplx(0d0,1d0)
      mimgl(0) = 1
      do  l = 1, n0
        mimgl(l) = (-img)**l
      enddo
C     Create bas,ips,rmax
      i_copy_size=size(ssite(1)%pos)
      do i_spackv=1,nbas
        call spackv_array_copy_r8_r8('u',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,bas)
      enddo

      i_copy_size=1;
      do i_spackv=1,nbas
        call spackv_array_copy_i8_i('u',ssite(i_spackv)%spec,i_copy_size,i_spackv+1-1,ips)
      enddo

      do  ib = 1, nbas
        is = int(ssite(ib)%spec)

        rmax ( ib ) = (sspec(is)%rmt)

      enddo
      nlmto = ndimh-napw

C --- Overlaps between IPW's ---
      call ipwovl(alat,plat,qlat,ngp,igv,ngp,igv,nbas,rmax,bas,ppovl)
C     call zprm('ppovl',2,ppovl,ngp,ngp,ngp)

C ... G vectors for the envelope (smooth hankel) functions
c      call iinit(iwk,3)
      call pshpr(0)
C     Find ngmx = number of G's in PW expansion of basis for this qp:
C     ngmx dimensions igvx,pwh.
c      print *,'pwmat: goto gvlst2 1111 gcutH=',gcutH
      call gvlst2(alat,plat,q,0,0,0,0d0,gcutH,0,100,1,
     .ngmx,xx,xx,xx,xx)
      allocate(igvx(3,ngmx),kv(3,ngmx))
c      call iinit(iwk,3)
c      print *,'pwmat: goto gvlst2 2222'
      call gvlst2(alat,plat,q,0,0,0,0d0,gcutH,0,102,ngmx, !1,
     .ngmx,kv,xx,xx,igvx)
      call poppr

c      print *,'pwmat: end of gvlst2'

C ... Patch: use igv for igvx
C     deallocate(igvx)
C     ngmx = ngp; allocate(igvx(3,ngmx)); igvx = igv

C --- Expansion of envelope basis functions in PWs ---
C     pwh(ig,j) = Fourier coff ig for envelope function j, where
C     ig is index to igvx
      allocate(pwh(ngmx,ndimh))
      pwh = 0d0
      allocate(yl(nlmax))
      lmxax = ll(nlmax)
C ... Fourier coefficients of smoothed hankels for all LMTOs
C     Could be optimized, ala hsibl, but not a critical step for GW.
      do  ig = 1, ngmx
        qpg(1:3) = tpiba * (q(1:3) + matmul(qlat, igvx(1:3,ig)))
        call ropyln(1,qpg(1),qpg(2),qpg(3),lmxax,1,yl,qpg2)
        do  ib = 1, nbas
          phase = exp( -img * sum( qpg*bas(:,ib)*alat )  )
          is = int(ssite(ib)%spec)

          call uspecb(0,2,sspec,is,is,lh,rsmh,eh,nkapi)
          call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
          call gtbsl1(8+16,norb,ltab,ktab,rsmh,eh,ntab,blks)

          do  io = 1, norb
            if (blks(io) .ne. 0) then
C           l,ik = l and kaph indices, needed to address eh,rsmh
              l  = ltab(io)
              ik = ktab(io)
C           offh = hamiltonian offset to this block
              offh  = offl(io)
              denom = eh(l+1,ik) - qpg2
              gam   = 1d0/4d0*rsmh(l+1,ik)**2
              offh  = offl(io)
              fach  = -pi4/vol/denom * phase * mimgl(l) * exp(gam*denom)
              do  ilm = l**2+1, (l+1)**2
                offh = offh+1
                pwh(ig,offh) = fach * yl(ilm)
              enddo
            endif
          enddo
        enddo
      enddo

C ... Fourier coefficients to APWs
C     APWs are normalized:  |G> = 1/sqrt(vol) exp[i G.r]
      do  iga = 1, napw
        pwh(:,iga+nlmto) = 0d0
C       Index to igvx that corresponds to igapw
c        ig = ifindiv2(igapw(1,iga),igvx,ngmx)
        ig = ifindiv2(igapw(1:3,iga)+inn(1:3),igvx,ngmx)
        pwh(ig,iga+nlmto) = 1d0/srvol
      enddo

C     call zprm('pwh',2,pwh,ngmx,ngmx,ndimh)

C --- Matrix elements between each (IPW,envelope function) pair ---
C     See Remarks
      allocate(ppovlx(ngp,ngmx))
      call ipwovl(alat,plat,qlat,ngp,igv,ngmx,igvx,nbas,rmax,bas,ppovlx)

      if(debug) print *,'sss:  pwmat ppovlx=',sum(abs(ppovlx))
      if(debug) print *,'sss:  pwmat pwh   =',sum(abs(pwh))

c takao find acml4.2.0 for gfortran gives wrong results.
      if(matmul_pwhovl()==1) then
        pwhovl= matmul(ppovlx,pwh)
      elseif(matmul_pwhovl()==2) then
        call matm(ppovlx,pwh,pwhovl,ngp,ngmx,ndimh)
      elseif(matmul_pwhovl()==3) then
        call zgemm('N','N',ngp,ndimh,ngmx,dcmplx(1d0,0d0),
     .  ppovlx,ngp,pwh,ngmx,dcmplx(0d0,0d0),pwhovl,ngp)
      else
        stop 'pwmat: matmul_pwovl is wrong'
      endif


C     call zprm('pwhovl',2,pwhovl,ngp,ngp,ndimh)
      if(debug) print *,'sss:  pwmat pwhovl=',sum(abs(pwhovl))

C --- Cleanup ---
      deallocate(yl,igvx,pwh,kv,ppovlx)

      end subroutine pwmat









      subroutine pwmat2(slat,ssite,sspec,nbas,ndimh,napw,igapw,iprmb,q,
     .ngp,nlmax,igv,ovlp,pwh)

      use m_struc_def  !Cgetarg

C- Matrix elements (IPW,IPW); expansion coefficients of basis fns in IPWs
C ----------------------------------------------------------------------
Ci Inputs
Ci   slat  :struct containing information about the lattice
Ci     Elts read: alat plat qlat vol
Ci     Stored:
Ci   ssite :struct containing site-specific information
Ci     Elts read: spec
Ci     Stored:    pos spec
Ci   sspec :struct containing species-specific information
Ci     Elts read: rmt
Ci     Stored:
Ci   nbas  :size of basis
Ci   ndimh :dimension of hamiltonian
Ci   napw  :number of augmented PWs in basis
Ci   igapw :list of G vectors for APW part of basis as multiples of
Ci         :primitive reciprocal lattice vectors
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   q     :Bloch wave vector
Ci   nlmax :maximum value of (1+augmentation l)**2 within any sphere
Ci   ngp   :no. G vectors for eigenfunction expansion (depends on q)
Ci         :Cutoff specified by QpGcut_psi in GWinput
Ci   igv   :list of ngp G vectors for PW expansion of basis as multiples
Ci         :of primitive reciprocal lattice vectors.
Co Outputs
Co   ovlp  :<IPW_G1 | IPW_G2 >
Co   pwh   :coefficients to PW expansion of basis functions, and also
Co         :therefore of IPW expansion of interstitial part of basis
Co         :pwh(ig,j) = Fourier coff ig for envelope function j
Co         :where ig is an index to igv
Co         :The interstitial part of the basis function phi_j is:
Co         :phi_j (istl) = sum_G  pwh(G,j) IPW(G)
Cl Local variables
Cr Remarks
Cr   pwmat2 is similar to pwmat, but returns pwh rather than pwhovl.
Cr   The reasons for this are discussed in sugw.f
Cu Updates
Cu   26 Jan 09 Adapted from pwmat
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ngp,nlmax,igv(3,ngp),nbas,ndimh,iprmb(1),
     .napw,igapw(3,napw)
      real(8):: q(3)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex ovlp(ngp,ngp), pwh(ngp,ndimh)
C ... Local parameters
      integer ips(nbas),ib,is,igetss,ig,lmxax,ll,nlmto,iga,ifindiv2
      double precision alat,plat(3,3),qlat(3,3),vol,pi,pi4,tpiba,xx,
     .tripl,dgetss,bas(3,nbas),rmax(nbas),qpg(3),qpg2,denom,gam,srvol
      integer n0,nkap0
      parameter (n0=10, nkap0=3)
      integer lh(nkap0),nkapi
      integer norb,io,l,ik,offh,ilm,ltab(n0*nkap0),ktab(n0*nkap0),
     .offl(n0*nkap0),blks(n0*nkap0),ntab(n0*nkap0)

      double precision eh(n0,nkap0),rsmh(n0,nkap0)
      double complex phase,img,fach,mimgl(0:n0)
      double precision,allocatable:: yl(:)

C --- Setup ---

      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      i_copy_size=size(slat%qlat)
      call dcopy(i_copy_size,slat%qlat,1,qlat,1)
      vol=slat%vol

      pi  = 4d0*datan(1d0)
      pi4 = 4d0*pi
      tpiba = 2*pi/alat
      vol = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      srvol = dsqrt(vol)
      img = dcmplx(0d0,1d0)
      mimgl(0) = 1
      do  l = 1, n0
        mimgl(l) = (-img)**l
      enddo
C     Create bas,ips,rmax
      i_copy_size=size(ssite(1)%pos)
      do i_spackv=1,nbas
        call spackv_array_copy_r8_r8('u',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,bas)
      enddo

      i_copy_size=1;
      do i_spackv=1,nbas
        call spackv_array_copy_i8_i('u',ssite(i_spackv)%spec,i_copy_size,i_spackv+1-1,ips)
      enddo

      do  ib = 1, nbas
        is = int(ssite(ib)%spec)

        rmax ( ib ) = (sspec(is)%rmt)

      enddo
      nlmto = ndimh-napw

C --- Overlaps between IPW's ---
      call ipwovl(alat,plat,qlat,ngp,igv,ngp,igv,nbas,rmax,bas,ovlp)
C     call zprm('ovlp',2,ovlp,ngp,ngp,ngp)

C --- Expansion pwh of envelope basis functions in PWs ---
      call dpzero(pwh,ngp*ndimh*2)
      allocate(yl(nlmax))
      lmxax = ll(nlmax)

C ... Fourier coefficients of smoothed hankels for all LMTOs
      do  ig = 1, ngp
        qpg(1:3) = tpiba * (q(1:3) + matmul(qlat, igv(1:3,ig)))
        call ropyln(1,qpg(1),qpg(2),qpg(3),lmxax,1,yl,qpg2)
        do  ib = 1, nbas
          phase = exp( -img * sum( qpg*bas(:,ib)*alat )  )
          is = int(ssite(ib)%spec)

          call uspecb(0,2,sspec,is,is,lh,rsmh,eh,nkapi)
          call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
          call gtbsl1(8+16,norb,ltab,ktab,rsmh,eh,ntab,blks)

          do  io = 1, norb
            if (blks(io) .ne. 0) then
C           l,ik = l and kaph indices, needed to address eh,rsmh
              l  = ltab(io)
              ik = ktab(io)
C           offh = hamiltonian offset to this block
              offh  = offl(io)
              denom = eh(l+1,ik) - qpg2
              gam   = 1d0/4d0*rsmh(l+1,ik)**2
              offh  = offl(io)
              fach  = -pi4/vol/denom * phase * mimgl(l) * exp(gam*denom)
              do  ilm = l**2+1, (l+1)**2
                offh = offh+1
                pwh(ig,offh) = fach * yl(ilm)
              enddo
            endif
          enddo
        enddo
      enddo

C ... Fourier coefficients to APWs
C     APWs are normalized:  |G> = 1/sqrt(vol) exp[i G.r]
      do  iga = 1, napw
C       Index to igv that corresponds to igapw
        ig = ifindiv2(igapw(1,iga),igv,ngp)
        pwh(ig,iga+nlmto) = 1d0/srvol
      enddo

C     call zprm('pwh',2,pwh,ngp,ngp,ndimh)

C --- Cleanup ---
      deallocate(yl)

      end subroutine pwmat2


      subroutine ipwovl(alat,plat,qlat,ng1,igv1,ng2,igv2,nbas,
     .rmax,bas,ppovl)

C- Overlap matrix elements between interstitial plane waves
C ----------------------------------------------------------------------
Ci Inputs
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   ng1   :number of G vectors G1
Ci   igv1  :list of G vectors G1 as multiples of lattice vectors
Ci   ng2   :number of G vectors G2
Ci   igv2  :list of G vectors G2 as multiples of lattice vectors
Ci   nbas  :size of basis
Ci   rmax  :augmentation radius, in a.u.
Ci   bas   :basis vectors, in units of alat
Co Outputs
Co   ppovl :<IPW1|IPW2> overlap matrix where IPW1 and IPW2 denote IPWs,
Co         :e.g. int exp(i (G2-G1).r) - int_(spheres) exp(i (G2-G1).r)
Cr Remarks
Cr   IPW(G1)+ * IPW(G2) = exp[i(G2-G1).r] -  (spheres) [i(G2-G1).r]
Cr                      = IPW(G2-G1)
Cr   <IPW1|IPW2> = int IPW(G2-G1)
Cr               = int exp[i(G2-G1).r] - int (spheres) [i(G2-G1).r]
Cu Updates
Cu   30 Mar 01 Taken from Kotani's mkppovl2
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,ng1,ng2,igv1(3,ng1),igv2(3,ng2)
      double precision alat,plat(3,3),qlat(3,3),rmax(nbas),bas(3,nbas)
      double complex ppovl(ng1,ng2)
C ... Local parameters
      integer nx(3),ig1,ig2,n1x,n2x,n3x,n1m,n2m,n3m
      double precision tripl,pi,tpibaqlat(3,3),vol
      double complex,allocatable :: ppox(:,:,:)

      pi  = 4d0*datan(1d0)
      vol = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      tpibaqlat = 2*pi/alat * qlat

C ... Find the range of G = G1*G2
      n1x = maxval( igv2(1,:)) - minval( igv1(1,:))
      n1m = minval( igv2(1,:)) - maxval( igv1(1,:))
      n2x = maxval( igv2(2,:)) - minval( igv1(2,:))
      n2m = minval( igv2(2,:)) - maxval( igv1(2,:))
      n3x = maxval( igv2(3,:)) - minval( igv1(3,:))
      n3m = minval( igv2(3,:)) - maxval( igv1(3,:))

      allocate( ppox(n1m:n1x,n2m:n2x,n3m:n3x) )
      ppox = 99999999999d0

C ... For each G in set of G2*G1, integral(cell-MT spheres) G2*G1
      do  ig1 = 1, ng1
        do  ig2 = 1, ng2
C       nx = G2 * G1
          nx(1:3) = igv2(1:3,ig2) - igv1(1:3,ig1)
          if (ppox(nx(1),nx(2),nx(3)) .eq. 99999999999d0) then
            call matgg2(alat,bas,rmax,nbas,vol,tpibaqlat,nx,
     .      ppox(nx(1),nx(2),nx(3)))
          endif
        enddo
      enddo

C ... For each G2*G1, poke integral into ppovl
      do  ig1 = 1,ng1
        do  ig2 = 1,ng2
          nx(1:3) = igv2(1:3,ig2) - igv1(1:3,ig1)
          ppovl(ig1,ig2) = ppox(nx(1),nx(2),nx(3) )
        enddo
      enddo

      deallocate(ppox)
      end subroutine ipwovl


      subroutine matgg2(alat,bas,rmax,nbas,vol,tpibaqlat,igv,ppovl)

C ----------------------------------------------------------------------
Ci Inputs
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   bas   :basis vectors, in units of alat
Ci   rmax  :augmentation radius, in a.u.
Ci   nbas  :size of basis
Ci   vol   :cell volume
Ci   tpibaqlat
Ci   igv   :List of G vectors, as integer multiples of primitive r.l.v.
Co Outputs
Co   ppovl :integral_cell IPW(G) d^3r = integral( PW - sum_ib PW_ib)
Cl Local variables
Cl   ggvec : G in atomic units
Cb Bugs
Cb   Sites with floating orbitals should not be looped over.
Cb   However, they have no effect since rmax=0 for those sites.
Cr Remarks
Cr   PW defined as exp[i G.r].  Thus  <G | G> = vol
Cr   IPW subtracts out parts inside sphere
Cr   The PW has expansion
Cr     exp iG.r = 4 pi sum_l i^l j_l(Gr) sum m_(-l,l) Y_L(hat r) Y_L(hat G)
Cr   The integral of the PW inside a sphere entails the l=0 part:
Cr     int_S iG.r d^3r = int d^3 4 pi Y_00^2 j_0(Gr) = 4 pi int dr r^2 j_0(Gr)
Cr   The l=0 Bessel function is j_0(x) = sin(x)/x.  For a sphere at the origin
Cr     int_S iG.r d^3r = 4 pi / G^3 int_(0,G rmax) dx x sin x
Cr                     = 4 pi / G^3 [sin(G rmax) - G rmax cos(G rmax)]
Cu Updates
Cu   30 Mar 01 adapted from Kotani
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,igv(3)
      double precision alat,bas(3,nbas),rmax(nbas),vol,tpibaqlat(3,3)
      double complex ppovl
C ... Local parameters
      integer ibas
      double precision pi4,absg,ggvec(3),grmx
      double complex img
      parameter (pi4=12.566370614359172d0)

      img = (0d0,1d0)
      ggvec = matmul(tpibaqlat,igv)
      absg  =  sqrt(sum(ggvec(1:3)**2))
C     Integral of exp(i G.r) = v0l * delta_(G,0)
      ppovl = 0d0
      if (absg .eq. 0d0) ppovl = vol

      do  ibas = 1, nbas
        if (absg .eq. 0d0) then
          ppovl = ppovl - pi4*rmax(ibas)**3/3d0
        else
          grmx  = absg*rmax(ibas)
          ppovl = ppovl -
     .    exp( img * sum(ggvec*bas(1:3,ibas))*alat )
     .    * pi4/absg**3 * ( -grmx * cos(grmx) + sin(grmx) )
        endif
      enddo
      end subroutine matgg2


