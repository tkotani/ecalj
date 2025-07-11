module m_pwmat
  use m_ll,only:ll
  use m_nvfortran,only:findloc
  public pwmat,mkppovl2
  private
  complex(8),allocatable:: ppovl_save(:,:,:)
  logical,allocatable:: has_ppovl(:,:,:)
  real(8):: vol_ppovl
  contains
! Matrix elements (IPW,IPW) and (IPW,envelope function)
subroutine pwmat(nbas,ndimh,napw,igapw,q,ngp,nlmax,igv,GcutH,ppovl,pwhovl)
  use m_struc_def     
  use m_lmfinit,only: alat=>lat_alat,ispec,n0,nkap0,pi,pi4,rmt_i=>rmt
  use m_lattic,only:  qlat=>lat_qlat, vol=>lat_vol,plat=>lat_plat,rv_a_opos
  use m_uspecb,only:uspecb
  use m_orbl,only: Orblib,ktab,ltab,offl,norb
  use m_ropyln,only: ropyln
  use m_blas,only: zmm => zmm_h, zmv => zmv_h, m_op_T
!  integer,parameter:: n0=10,nkap0=3
!  real(8),parameter:: pi=4d0*datan(1d0), pi4=4d0*pi
  implicit none
  !!   We have  q+G(igvx; internal in pwmat) = qp + G(igapw).
  !!  (because qp may be shortedned q (when w(oigv2)=igapw is generated) igapw is used in hambl.
  !!   Thus, we need to use
  !!           igv(internally in pwmat) = igapw 
  ! ----------------------------------------------------------------------
  !i   nbas  :size of basis
  !i   ndimh :dimension of hamiltonian
  !i   napw  :number of augmented PWs in basis
  !i   igapw :list of G vectors for APW part of basis as multiples of
  !i         :primitive reciprocal lattice vectors
  !i   q     :Bloch wave vector
  !i   nlmax :maximum value of (1+augmentation l)**2 within any sphere
  !i   ngp   :number of G vectors for eigenfunction expansion (depends on q)
  !i         :Cutoff specified by QpGcut_psi in GWinput
  !i   igv   :list of ngp G vectors for PW expansion of basis as multiples
  !i         :of primitive reciprocal lattice vectors.
  !i   GcutH :G cutoff for smoothed Hankel basis, LDA expansion
  !o   ppovl : <IPW_G1 | IPW_G2 >
  !o   pwhovl: <IPW_G1 | basis function = smooth Hankel or APW >
  !o         : matrix element; see Remarks
  !l Local variables
  !l   igvx  :same function as igv : a list of G-vectors that is used
  !l         :to expand eigenfunctions.  Difference is cutoff:
  !l         :igvx : cutoff is fixed by LMTO input HAM->GMAX
  !l         :igv  : cutoff is fixed by GWinput QpGcut_psi
  !l         :This inconsistency probably should be removed.
  !l   pwh   :coefficients to PW expansion of basis functions, and also
  !l         :therefore of IPW expansion of interstitial part of basis
  !l         :pwh(ig,j) = Fourier coff ig for envelope function j
  !r Remarks
  !r   IPW(G) denotes projected plane wave, i.e. PW with parts within MT
  !r   spheres projected out.
  !r   The interstitial part of the basis function phi_j is:
  !r      phi_j (istl) = sum_G  pwh(G,j) IPW(G)
  !r   Matrix elements between IPWs and basis functions are then
  !r      pwhovl(G,j) = int [IPW(G) sum_G1 pwh(G1,j) IPW(G1)]
  !r                = sum_G1 int [IPW(G) IPW(G1)] * pwh(G1,j) =
  !r                = sum_G1 ppovl(G,G1) pwh(G1,j) = ppovl * pwh
  !u Updates
  !u   29 Jan 09 Incorporate APW basis
  !u   25 Aug 04 Adapted to extended local orbitals
  !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
  !u   09 Apr 01 Adapted from Kotani's pplmat2
  ! ----------------------------------------------------------------------
  integer:: ngp,nlmax,igv(3,ngp),nbas,ndimh, napw,igapw(3,napw),&
       ips(nbas),ib,is,igetss,ngmx,ig,lmxax,iwk(3),nlmto,iga, ifindiv2,&
       matmul_pwhovl,&
       lh(nkap0),nkapi,io,l,ik,offh,ilm,blks(n0*nkap0),ntab(n0*nkap0),oi,ol,igx
  real(8):: q(3),GcutH,tpiba,xx,tripl,dgetss,bas(3,nbas),rmax(nbas),qpg(3),qpg2(1),denom,gam,srvol,&
       eh(n0,nkap0),rsmh(n0,nkap0)
  complex(8):: ppovl(ngp,ngp), pwhovl(ngp,ndimh),phase,img,fach,mimgl(0:n0)
  integer,allocatable:: igvx(:,:),kv(:,:)
  real(8),allocatable:: yl(:)
  complex(8),allocatable:: pwh(:,:),ppovlx(:,:)
  logical:: debug=.false.
  integer :: istat
  tpiba = 2*pi/alat !vol = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  srvol = dsqrt(vol)
  img = dcmplx(0d0,1d0)
  mimgl(0) = 1
  do l = 1, n0
     mimgl(l) = (-img)**l
  enddo
  do  ib = 1, nbas
     bas(:,ib)=rv_a_opos(:,ib)
     ips(ib) = ispec(ib)
     is = ispec(ib) 
     rmax(ib) = rmt_i(is)
  enddo
  nlmto = ndimh-napw
  call ipwovl(alat,plat,qlat,ngp,igv,ngp,igv,nbas,rmax,bas,ppovl)! --- Overlaps between IPW's ---
  ! ... G vectors for the envelope (smooth hankel) functions
  !     Find ngmx = number of G's in PW expansion of basis for this qp:
  call pshpr(0)
!  call gvlst2(alat,plat,q,0,0,0,0d0,gcutH,0,0,1,    ngmx,xx,xx,xx)
  call getgv2(alat,plat,qlat,q, gcutH,1, ngmx,xx) 
  allocate(igvx(3,ngmx))!,kv(3,ngmx))
!  block
!    integer:: igvxx(ngmx,3)
!    igvxx=0
    !call gvlst2(alat,plat,q,0,0,0,0d0,gcutH,0,1,ngmx, ngmx,kv,xx,igvxx)
    call getgv2(alat,plat,qlat,q, gcutH,2, ngmx,igvx) 
!    igvx=transpose(igvxx)
!  endblock
  call poppr
  ! --- Expansion of envelope basis functions in PWs ---
  !     pwh(ig,j) = Fourier coff ig for envelope function j, where ig is index to igvx
  ! allocate(pwh(ngmx,ndimh),yl(nlmax))
  allocate(yl(nlmax))
  ! pwh = 0d0
  lmxax = ll(nlmax)
  ! ... Fourier coefficients of smoothed hankels for all LMTOs
  !     Could be optimized, ala hsibl, but not a critical step for GW.
  block
    use m_stopwatch
    type(stopwatch) :: sw1, sw2,sw3
    integer:: ig_start, ig_end, igp
    integer, parameter:: ngblock = 1024
    call stopwatch_init(sw1, 'set ppovlx,pwh')
    call stopwatch_init(sw2, 'zmm')
    call stopwatch_init(sw3, 'zmm2')
    call get_ppovl_init(ngp, igv, ngmx, igvx)
    pwhovl(:,:) = (0d0, 0d0)
    gblock_loop: do ig_start = 1, ngmx, ngblock
      ig_end = min(ngmx, ig_start + ngblock -1)
      allocate(pwh(ndimh,ig_start:ig_end), source = (0d0,0d0))
      allocate(ppovlx(ngp,ig_start:ig_end), source = (0d0,0d0))
      call stopwatch_start(sw1)
      do ig = ig_start, ig_end
        qpg(1:3) = tpiba * (q(1:3) + matmul(qlat, igvx(1:3,ig)))
        call ropyln(1,qpg(1),qpg(2),qpg(3),lmxax,1,yl,qpg2)
        do ib = 1, nbas
          phase = exp( -img * sum( qpg*bas(:,ib)*alat )  )
          is = ispec(ib) 
          call uspecb(is,rsmh,eh)
          call orblib(ib) !return norb,ltab,ktab,offl
          call gtbsl8(norb,ltab,ktab,rsmh,eh,ntab,blks)
          do  io = 1, norb
            l  = ltab(io) ! l,ik = l and kaph indices, needed to address eh,rsmh
            ik = ktab(io)
            ol = ltab(io)**2
            oi = offl(io) ! offh = hamiltonian offset to this block
            denom = eh(l+1,ik) - qpg2(1)
            gam   = 1d0/4d0*rsmh(l+1,ik)**2
            if(abs(denom)<1d-10) cycle  !2023feb ok?
            fach  = -pi4/vol/denom * phase * mimgl(l) * exp(gam*denom)
            pwh(oi+1:oi+blks(io),ig) = fach * yl(ol+1:ol+blks(io))
          enddo
        enddo
        do igp=1, ngp
          ppovlx(igp,ig) = get_ppovl(alat, plat, qlat, rmax, bas, nbas, igv(1:3,igp), igvx(1:3,ig))
        enddo
      enddo
      call stopwatch_pause(sw1)
      call stopwatch_start(sw3)
      pwh(nlmto+1:,:) = 0d0
      do iga= 1,napw 
         ig = findloc([(all(igapw(:,iga)==igvx(:,igx)),igx=ig_start, ig_end)],value=.true.,dim=1)
         !if current igx's (ig_start:ig_end) has overlap with iga
         if(1 <= ig .and. ig <= ig_end - ig_start + 1) pwh(iga+nlmto,ig + ig_start - 1) = 1d0/srvol
      enddo
      call stopwatch_pause(sw3)
      call stopwatch_start(sw2)
      istat = zmm(ppovlx, pwh, pwhovl, m=ngp, n=ndimh, k=(ig_end-ig_start+1), opB=m_op_T, beta = (1D0, 0d0))
      deallocate(ppovlx, pwh)
      call stopwatch_pause(sw2)
    enddo gblock_loop
    call get_ppovl_finalize()
    call stopwatch_show(sw1)
    call stopwatch_show(sw2)
    call stopwatch_show(sw3)
  endblock
  ! ... Fourier coefficients to APWs. APWs are normalized:  |G> = 1/sqrt(vol) exp[i G.r]
  ! --- Matrix elements between each (IPW,envelope function) pair ---
  ! allocate(ppovlx(ngp,ngmx))
  ! call ipwovl(alat,plat,qlat,ngp,igv,ngmx,igvx,nbas,rmax,bas,ppovlx)
  ! pwhovl= matmul(ppovlx,pwh)
  ! 2024-11-09 MO replaced matmul by a BLAS call because matmul in mic(intel) was very slow
  ! istat = zmm(ppovlx, pwh, pwhovl, m=ngp, n=ndimh, k=ngmx)
  ! deallocate(yl,igvx,pwh,ppovlx)
  deallocate(yl,igvx)
end subroutine pwmat
subroutine ipwovl(alat,plat,qlat,ng1,igv1,ng2,igv2,nbas, rmax,bas,ppovl)
  !- Overlap matrix elements between interstitial plane waves
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   plat  :primitive lattice vectors, in units of alat
  !i   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
  !i   ng1   :number of G vectors G1
  !i   igv1  :list of G vectors G1 as multiples of lattice vectors
  !i   ng2   :number of G vectors G2
  !i   igv2  :list of G vectors G2 as multiples of lattice vectors
  !i   nbas  :size of basis
  !i   rmax  :augmentation radius, in a.u.
  !i   bas   :basis vectors, in units of alat
  !o Outputs
  !o   ppovl :<IPW1|IPW2> overlap matrix where IPW1 and IPW2 denote IPWs,
  !o         :e.g. int exp(i (G2-G1).r) - int_(spheres) exp(i (G2-G1).r)
  !r Remarks
  !r   IPW(G1)+ * IPW(G2) = exp[i(G2-G1).r] -  (spheres) [i(G2-G1).r]
  !r                      = IPW(G2-G1)
  !r   <IPW1|IPW2> = int IPW(G2-G1)
  !r               = int exp[i(G2-G1).r] - int (spheres) [i(G2-G1).r]
  !u Updates
  !u   30 Mar 01 Taken from Kotani's mkppovl2
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nbas,ng1,ng2,igv1(3,ng1),igv2(3,ng2)
  real(8):: alat,plat(3,3),qlat(3,3),rmax(nbas),bas(3,nbas),tripl,pi,tpibaqlat(3,3),vol
  complex(8):: ppovl(ng1,ng2)
  integer :: nx(3),ig1,ig2,n1x,n2x,n3x,n1m,n2m,n3m
  complex(8),allocatable :: ppox(:,:,:)
  pi  = 4d0*datan(1d0)
  vol = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  ! ... Find the range of G = G1*G2
  n1x = maxval( igv2(1,:)) - minval( igv1(1,:))
  n1m = minval( igv2(1,:)) - maxval( igv1(1,:))
  n2x = maxval( igv2(2,:)) - minval( igv1(2,:))
  n2m = minval( igv2(2,:)) - maxval( igv1(2,:))
  n3x = maxval( igv2(3,:)) - minval( igv1(3,:))
  n3m = minval( igv2(3,:)) - maxval( igv1(3,:))
  allocate( ppox(n1m:n1x,n2m:n2x,n3m:n3x) )
  ppox = 99999999999d0
  ! ... For each G in set of G2*G1, integral(cell-MT spheres) G2*G1
  do  ig1 = 1, ng1
     do  ig2 = 1, ng2
        nx(1:3) = igv2(1:3,ig2) - igv1(1:3,ig1)
        if (ppox(nx(1),nx(2),nx(3)) == 99999999999d0) then
           call matgg2(alat,bas,rmax,nbas,vol,2*pi/alat*qlat,nx, ppox(nx(1),nx(2),nx(3)))
        endif
     enddo
  enddo
  ! ... For each G2*G1, poke integral into ppovl
  do  ig1 = 1,ng1
     do  ig2 = 1,ng2
        nx(1:3) = igv2(1:3,ig2) - igv1(1:3,ig1)
        ppovl(ig1,ig2) = ppox(nx(1),nx(2),nx(3) )
     enddo
  enddo
  deallocate(ppox)
end subroutine ipwovl
subroutine matgg2(alat,bas,rmax,nbas,vol,tpibaqlat,igv,ppovl)
  use m_lmfinit,only:pi4
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   bas   :basis vectors, in units of alat
  !i   rmax  :augmentation radius, in a.u.
  !i   nbas  :size of basis
  !i   vol   :cell volume
  !i   tpibaqlat
  !i   igv   :List of G vectors, as integer multiples of primitive r.l.v.
  !o Outputs
  !o   ppovl :integral_cell IPW(G) d^3r = integral( PW - sum_ib PW_ib)
  !l Local variables
  !l   ggvec : G in atomic units
  !b Bugs
  !b   Sites with floating orbitals should not be looped over.
  !b   However, they have no effect since rmax=0 for those sites.
  !r Remarks
  !r   PW defined as exp[i G.r].  Thus  <G | G> = vol
  !r   IPW subtracts out parts inside sphere
  !r   The PW has expansion
  !r     exp iG.r = 4 pi sum_l i^l j_l(Gr) sum m_(-l,l) Y_L(hat r) Y_L(hat G)
  !r   The integral of the PW inside a sphere entails the l=0 part:
  !r     int_S iG.r d^3r = int d^3 4 pi Y_00^2 j_0(Gr) = 4 pi int dr r^2 j_0(Gr)
  !r   The l=0 Bessel function is j_0(x) = sin(x)/x.  For a sphere at the origin
  !r     int_S iG.r d^3r = 4 pi / G^3 int_(0,G rmax) dx x sin x
  !r                     = 4 pi / G^3 [sin(G rmax) - G rmax cos(G rmax)]
  implicit none
  integer:: nbas,igv(3),ibas
  real(8):: alat,bas(3,nbas),rmax(nbas),vol,tpibaqlat(3,3),absg,ggvec(3),grmx
  complex(8):: ppovl,  img = (0d0,1d0)
  ggvec = matmul(tpibaqlat,igv)
  absg  =  sqrt(sum(ggvec(1:3)**2))
  !     Integral of exp(i G.r) = v0l * delta_(G,0)
  ppovl = 0d0
  if(absg == 0d0) ppovl = vol
  do ibas = 1, nbas
     if (absg==0d0) then
        ppovl= ppovl - pi4*rmax(ibas)**3/3d0
     else
        grmx = absg*rmax(ibas)
        ppovl= ppovl - exp(img*sum(ggvec*bas(1:3,ibas))*alat)*pi4/absg**3 *(-grmx*cos(grmx)+sin(grmx))
     endif
  enddo
end subroutine matgg2
subroutine mkppovl2(alat,plat,qlat,ng1,ngvec1,ng2,ngvec2,nbas,rmax,bas, ppovl)
  !- < G1 | G2 > matrix where G1 denotes IPW, zero within MT sphere.
  implicit none
  integer::  nbas, ng1,ng2,nx(3),ig1,ig2, ngvec1(3,ng1),ngvec2(3,ng2), &
       n1x,n2x,n3x,n1m,n2m,n3m
  real(8) :: tripl,rmax(nbas),pi = 4d0*datan(1d0)
  real(8) :: plat(3,3),qlat(3,3), &
       alat,tpibaqlat(3,3),voltot, bas(3,nbas)
  complex(8) :: ppovl(ng1,ng2)
  complex(8),allocatable :: ppox(:,:,:)
  logical:: debug=.false.
  voltot    = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  tpibaqlat =  2*pi/alat *qlat
  ! < G1 | G2 >
  n1x = maxval( ngvec2(1,:)) - minval( ngvec1(1,:))
  n1m = minval( ngvec2(1,:)) - maxval( ngvec1(1,:))
  n2x = maxval( ngvec2(2,:)) - minval( ngvec1(2,:))
  n2m = minval( ngvec2(2,:)) - maxval( ngvec1(2,:))
  n3x = maxval( ngvec2(3,:)) - minval( ngvec1(3,:))
  n3m = minval( ngvec2(3,:)) - maxval( ngvec1(3,:))
  if(debug) print *,' mkppovl2: 1 ',n1x,n1m,n2x,n2m,n3x,n3m
  allocate( ppox(n1m:n1x,n2m:n2x,n3m:n3x) )
  ppox = 1d99
  do ig1  = 1, ng1
     do ig2  = 1, ng2
        nx(1:3) = ngvec2(1:3,ig2) - ngvec1(1:3,ig1) ! G2-G1
        if( ppox(nx(1),nx(2),nx(3)) == 1d99 ) then
           call matgg2(alat,bas,rmax,nbas,voltot, tpibaqlat, &
                nx(1:3),ppox( nx(1),nx(2),nx(3)))
        endif
     enddo
  enddo
  if(debug) print *,' mkppovl2: 2 ',n1x,n1m,n2x,n2m,n3x,n3m
  do ig1 = 1,ng1
     do ig2 = 1,ng2
        nx(1:3) = ngvec2(1:3,ig2) -ngvec1(1:3,ig1) ! G2-G1
        ppovl(ig1,ig2) = ppox( nx(1),nx(2),nx(3) )
     enddo
  enddo
  deallocate(ppox)
end subroutine mkppovl2
!2024-11-28 MO added the set of get_ppovl interfaces
subroutine get_ppovl_init(ng1, igv1, ng2, igv2)
  implicit none
  integer, intent(in) :: ng1, ng2, igv1(3,ng1), igv2(3,ng2)
  integer :: nx(3), n1x, n1m, n2x, n2m, n3x, n3m
  n1x = maxval( igv2(1,:)) - minval( igv1(1,:))
  n1m = minval( igv2(1,:)) - maxval( igv1(1,:))
  n2x = maxval( igv2(2,:)) - minval( igv1(2,:))
  n2m = minval( igv2(2,:)) - maxval( igv1(2,:))
  n3x = maxval( igv2(3,:)) - minval( igv1(3,:))
  n3m = minval( igv2(3,:)) - maxval( igv1(3,:))
  allocate(ppovl_save(n1m:n1x,n2m:n2x,n3m:n3x), source = (0d0, 0d0))
  allocate( has_ppovl(n1m:n1x,n2m:n2x,n3m:n3x), source = .false.)
  ! write(06,*) 'get_ppovl_init: ', n1x, n1m, n2x, n2m, n3x, n3m
  ! write(06,*) 'get_ppovl_init: size ppovl', size(ppovl_save)
end subroutine get_ppovl_init
function get_ppovl(alat, plat, qlat, rmax, bas, nbas, igv1 , igv2) result(ppovl)
  use m_lmfinit, only: pi
  implicit none
  integer, intent(in) :: igv1(3), igv2(3), nbas
  integer :: nx(3)
  complex(8) :: ppovl
  real(8) :: alat, plat(3,3), qlat(3,3), rmax(nbas), bas(3,nbas), vol,tripl
  vol = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  nx(:) = igv2(:) - igv1(:)
  if (.not.has_ppovl(nx(1),nx(2),nx(3))) then
    call matgg2(alat, bas, rmax, nbas, vol, 2*pi/alat*qlat, nx, ppovl_save(nx(1),nx(2),nx(3)))
    has_ppovl(nx(1),nx(2),nx(3)) = .true.
  endif
  ppovl = ppovl_save(nx(1),nx(2),nx(3))
end function get_ppovl
subroutine get_ppovl_finalize()
  deallocate(has_ppovl, ppovl_save)
end subroutine get_ppovl_finalize
end module m_pwmat

