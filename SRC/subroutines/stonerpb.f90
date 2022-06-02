subroutine stonerpb(nq,nw,nmbas,nbloch, qp,momsite,mmnorm &
     ,emesh,zxq)

  !=============================================================
  ! ---  <~e| X | e~>  matrix calculation. by using X^{-1}=X0^{-1}+U
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nq       : total number of q points in calculation of X(q,w)
  !i   nw       : total number of w points in calculation of X(q,w)
  !i   nmbas    : number of magnetic sites  = nmag
  !i   qp       : qp(3,nq)  vector of each q
  !i   momsite  : magnetic moment m
  !i   mmnorm   : <m|m>
  !i   eiqrm    : eiqrm_i = <~e_i|e^{iqr}> =  <M_i|eiqr>/sqrt(<M_i|M_i>)
  !i   emesh    : energy w mesh
  !i   x0et      : <~e|X0|~e>
  !o Outputs:
  !o
  !o
  !o
  !o
  !l Local variables
  !l   rho0    : radial part of true electron density on each mesh (xyz) = rho1/r**2
  !l   rho0xyz : true electron density on each mesh point specified by (ir,ip)
  !l             rho0xyz(ir,ip,isp)=Sum_{ilm}{ rho0(ir,ilm,isp)*yl(ip,ilm) }
  !l   vxc     : vxc(ir,ip,1:2) for spin_1 and spin_2 ;
  !l             vxc(:,:,4) = ( vxc(:,:,1)-vxc(:,:,2) ) / 2
  !l             vxc(:,:,3) = idol
  !l   bom     : bom(ir,ip) = b(r)/m(r)
  !l   sum_bibj: <B_i | B_j> integrate inside sphere
  !l   sum_bibombj: <B_i | b(r)/m(r) | B_j>
  !r Remarks
  !r    In ASA, m(r_) is a radial function, but not in FP.
  !r    we need to build a spherical mesh, and sum the rho on each mesh
  !r    point.
  !r    Vxc generated from LMTO package is in units of Ryberg.
  !r    Numerical spherical mesh is tabulated in wxp.chk file, which is generated in LMTO
  !r    package ./lmfsph
  !u Updates
  !u    20 May 08 First created

  !     implicit none
  ! ... Passed parameters
  integer :: nq,nw,nmbas,nbloch
  real(8) momsite(nmbas), mmnorm(nmbas),emesh(nw),qp(3,nq)
  complex(8) zxq(nbloch,nbloch,nq,nw)

  ! ... Local parameters
  integer :: jb
  integer :: iw,i,i1,i2,i3,j,j1,j2,j3, nwx,ix,iy,iq, nw_intp &
       ,imax,imin,jmax,jmin, iwpole,iwpolf,iwpola
  complex(8),allocatable :: &
       x0mean(:,:,:), x0inv(:,:,:,:), xinv(:,:,:,:), x2et(:,:,:,:) &
       , qxq_intp(:,:), xinvh(:,:,:,:),dxidw(:,:,:,:) &
       ,dxidw_eb(:,:,:,:)
  real(8), allocatable :: mxevl_xinvh(:,:), freq2(:) &
       ,mxevl2_xinvh(:,:)
  real(8) uub(nmbas,nw), uu(nmbas,nw), eval(nmbas)
  real(8) emin_intp,dw_intp,emax_intp
  real(8) freq(nw),freqm(nw),freq_mev(nw), mmnorm2(nmbas)
  complex(8):: cxtmp(nmbas,nmbas),img=(0d0,1d0), qxq(nw,nq) &
       ,meffi(nmbas,nmbas,nq),   xinvh_w0(nmbas,nmbas,nq) &
       ,meffi_eb(nmbas,nmbas,nq),xinvh_w0eb(nmbas,nmbas,nq) &
       ,meffi2(nmbas,nmbas,nq),meffi2_eb(nmbas,nmbas,nq)
  real(8) qxq_r(nw,nq),qxq_i(nw,nq), rydberg,polinta
  real(8) qxq_inv_r(nw,nq),qxq_inv_i(nw,nq)
  real(8) rtmp,rtmp1, rtmp2,rtmp3,rtmp4, rymev,omg
  real(8) eout,rrrx,sumx,cccxmx,elimit,iiix,epole &
       ,epole_fm(nq),epole_af(nq),vpole_fm(nq),vpole_af(nq)
  real(8) , allocatable:: freq_intp(:),e_intp(:)
  parameter (rydberg=13.6058d0)
  external :: polinta

  ! for jmat calculation
  complex(8)  oo(nmbas,nmbas),evc(nmbas,nmbas) &
       ,jjmat(nmbas,nmbas),sqm(nmbas,nmbas)
  integer:: nmx,nev

  nwx = 400
  rymev = rydberg*1d3

  !      allocate( x0inv(nmbas,nmbas,nq,nwx) )
  !      allocate( xinv(nmbas,nmbas,nq,nwx), x2et(nmbas,nmbas,nq,nwx) )
  !      allocate( xinvh(nmbas,nmbas,nq,nwx),mxevl_xinvh(nq,nwx) )
  !      allocate( dxidw(nmbas,nmbas,nq,nwx))
  !      allocate( dxidw_eb(nmbas,nmbas,nq,nwx))
  !      allocate( mxevl2_xinvh(nq,nwx) )

  ! ... sanity check
  if ( emesh(1) /= 0d0 ) call rx('w(1) /= 0')
  freq=emesh

  !./lmf lsmo56 '--chimedit~new 2~read tkrs'

  jb=1
  open(111, file='X0pbqw.allqb')
  write(111,"( 4i4 )") nq,nmbas,nw
  do iw=1,nw
     write(111,301)  (zxq( jb , jb,iq,iw), iq=1,nq )
  enddo
301 format(1000d23.15)
end subroutine stonerpb

