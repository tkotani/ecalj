module m_stonerrsa
  private
  contains
subroutine stonerrsa(nq,nw,nmbas,qp,momsite,mmnorm,eiqrm,freq, x0et)
  use m_lgunit,only: stdo
  use m_ftox
  !- Transverse susceptibility matrix <~e| X |e~> from X^{-1}=X0^{-1}+U
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nq     :number of q points for which to calculate X(q,w)
  !i   nw     :number of frequencies where input x0et is given
  !i   nmbas  :number of magnetic sites
  !i   qp     :qp(3,nq)  vector of each q
  !i   momsite:magnetic moment m
  !i   mmnorm :<m|m>
  !i   eiqrm  :eiqrm_i = <~e_i|e^{iqr}> =  <M_i|eiqr>/sqrt(<M_i|M_i>)
  !i   freq   :frequency mesh
  !i   x0et   :<~e|X0|~e>
  !o Outputs:
  !o   ... finish this
  !l Local variables
  !r Remarks
  !u Updates
  !u   07 Feb 09 (L. Ke) adapted from T. Kotani
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nq,nw,nmbas
  real(8) momsite(nmbas), mmnorm(nmbas),freq(nw),qp(3,nq)
  complex(8) eiqrm(nmbas,nq),x0et(nmbas,nmbas,nq,nw)
  ! ... Local parameters
  integer :: i,ifi,ifi2,imax,imin,iq,iw,iwpola,iwpolf,ipr, &
       ix,iy,j,jmax,jmin,nglob,nw_intp,nwx
  complex(8),allocatable :: qxq_intp(:,:),wk(:,:), &
       x0inv(:,:,:,:),xinv(:,:,:,:),x2et(:,:,:,:), &
       xinvh(:,:,:,:),dxidw(:,:,:,:),dxidw_eb(:,:,:,:)
  real(8),allocatable:: mxevl_xinvh(:,:),mxevl2_xinvh(:,:),e_intp(:)
  real(8) uub(nmbas,nw),uu(nmbas,nw),eval(nmbas)
  real(8) freqm(nw),mmnorm2(nmbas)
  real(8) emin_intp,dw_intp,emax_intp
  complex(8):: cxtmp(nmbas,nmbas),img=(0d0,1d0),qxq(nw,nq), &
       meffi(nmbas,nmbas,nq),xinvh_w0(nmbas,nmbas,nq), &
       meffi_eb(nmbas,nmbas,nq),xinvh_w0eb(nmbas,nmbas,nq), &
       meffi2(nmbas,nmbas,nq),meffi2_eb(nmbas,nmbas,nq)
  real(8) qxq_r(nw,nq),qxq_i(nw,nq),rydberg,polinta
  real(8) qxq_inv_r(nw,nq),qxq_inv_i(nw,nq), &
       epole_fm(nq),epole_af(nq),vpole_fm(nq),vpole_af(nq)
  real(8) rtmp,rtmp1,rtmp2,rymev,epole
  external :: polinta
  !     For calculating JMAT:
  complex(8) oo(nmbas,nmbas),evc(nmbas,nmbas),sqm(nmbas,nmbas), &
       jjmat(nmbas,nmbas),jjmat2(nmbas,nmbas),mxi(nmbas,nmbas)
  integer:: nmx,nev
  parameter (rydberg=13.6058d0)

  call getpr(ipr)
  nwx = 400
  rymev = rydberg*1d3
  if (freq(1) /= 0d0) &
       call rx1('nonzero 1st frequency, w(1) = %g',freq(1))
  allocate(x0inv(nmbas,nmbas,nq,nwx))
  allocate(xinv(nmbas,nmbas,nq,nwx),x2et(nmbas,nmbas,nq,nwx) )
  allocate(xinvh(nmbas,nmbas,nq,nwx),mxevl_xinvh(nq,nwx) )
  allocate(dxidw(nmbas,nmbas,nq,nwx))
  allocate(dxidw_eb(nmbas,nmbas,nq,nwx))
  allocate(mxevl2_xinvh(nq,nwx) )

! #if DEBUG
!   open(ifi,file='x0qw.allqb',form='unformatted')
!   write (ifi) nq,nmbas,nw,nwx,qp
!   write (ifi) freq(1:nw)
!   do  iq = 1, nq
!      do  iw = 1, nw
!         write(ifi) dble(x0et(1:nmbas,1:nmbas,iq,iw))
!         write(ifi) dimag(x0et(1:nmbas,1:nmbas,iq,iw))
!      enddo
!   enddo
!   close(ifi)
! #endif

  ! --- Determine <~e|U(w)|e~> by only using boundary condition at q=0 ---
  x0inv(:,:,1,:) = x0et(:,:,1,:) ! only need x0(q=0,w)
  mmnorm2 = mmnorm**2
  do  iw = 1, nwx
     call matcinv(nmbas, x0inv(:,:,1,iw)) ! x0inv=<~e|x0^-1|~e>
     uub(1:nmbas,iw) = &
          matmul(dble(x0inv(1:nmbas,1:nmbas,1,iw)),mmnorm(1:nmbas))
     do  i = 1, nmbas
        if (abs(momsite(i)) > 1d-3) then
           uu(i,iw) = freq(iw)*momsite(i)/mmnorm2(i) &
                - uub(i,iw)/mmnorm(i)
        endif
     enddo
     !      write(stdo,"(i4, f13.5, d15.7)") iw,rymev*freq(iw),uu(1,iw)
  enddo

! #if DEBUG
!   !     open(112,file='etUet.allw')
!   !     ifi = fopnx('etUet.allw',2,2,-1)
!   open(newunit=ifi,file='etUet')
!   do  iw = 1, nwx
!      write(ifi,'(f21.13,3x,255d18.10)') &
!           freq(iw),(uu(ix,iw),ix=1,nmbas)
!   enddo
!   close(ifi)
! #endif

  ! --- x2et=<~e|X(w)|e~>  &&  <eiqr|~e><~e|X|~e><eiqr|~e> ---
  write(stdo,*)' STONERRSA:  calculate full Xi for each q'// &
       ' Magnetic moments, by (magnetic) site:'
  if (ipr >= 20) write(stdo,"(1x,12f8.4)") momsite
  !     write(stdo,"( 5x,'q',10x,'qvec' )")
  x0inv = x0et
  do  iq = 1, nq
     !     write(stdo,'( 3x,i3,3f7.3 )') iq, qp(1:3,iq)
     do  iw = 1, nwx
        call matcinv(nmbas,x0inv(:,:,iq,iw))
        xinv(:,:,iq,iw) = x0inv(:,:,iq,iw)
        do  ix = 1, nmbas
           xinv(ix,ix,iq,iw) = x0inv(ix,ix,iq,iw) + uu(ix,iw)
        enddo
        cxtmp(:,:) = xinv(:,:,iq,iw)
        !       Liqin:  why is this in a do loop?
        do  ix = 1, nmbas
           cxtmp(:,:) = cxtmp(:,:) + img*1d-30
        enddo
        !       Liqin: this is poor programming technique
        call matcinv(nmbas,cxtmp(:,:)) !this is full x_+-  !Matrix inversion.
        x2et(:,:,iq,iw) = cxtmp(:,:) ! x2et=<~e|X(w)|e~>
        qxq(iw,iq) = sum( eiqrm(:,iq)& !  <eiqr|~e><~e|X|~e><~|eiqr>
             *matmul(x2et(:,:,iq,iw),dconjg(eiqrm(:,iq) )))
        qxq_r(iw,iq) = dble(qxq(iw,iq) )
        qxq_i(iw,iq) = dimag(qxq(iw,iq) )

        qxq_inv_r(iw,iq) = dble( 1d0/qxq(iw,iq) ) !for interpolation
        qxq_inv_i(iw,iq) = dimag( 1d0/qxq(iw,iq) )

        !      write(stdo,"( i4, f13.5, 2d15.7 )") iw,rymev*freq(iw)
        !     .             , qxq_r(iw,iq),qxq_i(iw,iq)[
     enddo
  enddo

  ! --- Interpolate <eiqr|X|eiqr> ---
  !     Liqin: these should be passed as arguments
  emin_intp = 0d0; dw_intp = 1d-2; emax_intp = 1000d0
  nw_intp = int((emax_intp-emin_intp)/dw_intp) + 1
  write(stdo,ftox)' Make <q|X|q>: for ',nwx,' energy points; '// &
       'emax = ',ftof(rydberg*freq(nwx)),'eV'
  write(stdo,ftox)' Interpolate energy window with '// &
       'emin:dw:emax = ',ftof([emin_intp,dw_intp,emax_intp]),'meV (',nw_intp,' points)'
  allocate(qxq_intp(nw_intp,nq),e_intp(nw_intp) )
  do  iq = 1, nq
     do  iw = 1, nw_intp
        if (iq == 1) then
           e_intp(iw) = (emin_intp+(iw-1)*dw_intp )/(rymev)
        endif
        rtmp1 = polinta(e_intp(iw),freq(1:nwx),qxq_inv_r(1:nwx,iq),nwx)
        rtmp2 = polinta(e_intp(iw),freq(1:nwx),qxq_inv_i(1:nwx,iq),nwx)
        qxq_intp(iw,iq) = 1d0/(rtmp1 + img*rtmp2) !for interpolation
     enddo
  enddo

  if (ipr >= 30) then
     write(stdo,*)' Data for pole search of evl(<q|X|q>)'// &
          'q qxq_r_max qxq_r_min qxq_i_max qxq_i_min'
     do  iq = 1, nq
        jmax=maxloc(-dble(qxq_intp(1:nw_intp,iq)),dim=1)
        jmin=minloc(-dble(qxq_intp(1:nw_intp,iq)),dim=1)
        imax=maxloc(-dimag(qxq_intp(1:nw_intp,iq)),dim=1)
        imin=minloc(-dimag(qxq_intp(1:nw_intp,iq)),dim=1)
        write(stdo,101) iq, &
             rymev*e_intp(jmax), -dble(qxq_intp(jmax,iq)), &
             rymev*e_intp(jmin), -dble(qxq_intp(jmin,iq)), &
             rymev*e_intp(imax), -dimag(qxq_intp(imax,iq)), &
             rymev*e_intp(imin), -dimag(qxq_intp(imin,iq))
101     format(3x,i3,1x, f7.2,'(',d15.7,')',1x &
             ,f7.2,'(' ,d15.7,')',1x,f7.2,'(' ,d15.7,')' &
             ,1x,f7.2,'(' ,d15.7,')' )
     enddo
  endif

! #if DEBUG
!   !      open(106,file='qxqi.allq')
!   !      open(107,file='qxqr.allq')
!   !      ifi  = fopna('qxqi',-1,0)
!   !     ifi2 = fopna('qxqr',-1,0)
!   open(ifi,file='qxqi')
!   open(ifi2,file='qxqr')
!   do  iw = 1, nw_intp
!      write(ifi,"( f13.5, 100d15.7)") rymev*e_intp(iw),(-dimag(qxq_intp(iw,iq)),iq=1,nq)
!      write(ifi2,"( f13.5, 100d15.7)") rymev*e_intp(iw),(-dble(qxq_intp(iw,iq)),iq=1,nq)
!   enddo
!   close(ifi)
!   close(ifi2)
! #endif
  ! ... Finished interpolation of <eiqr|X|eiqr>

  ! --- Eigenvalue of hermitian xinvh ---
  do  iq = 1, nq
     do  iw = 1, nwx
        xinvh(:,:,iq,iw) = .5d0*( xinv(:,:,iq,iw) &
             + transpose(dconjg(xinv(:,:,iq,iw))) )
        call zevl(nmbas,xinvh(1,1,iq,iw),eval)
        mxevl_xinvh(iq,iw) = maxval(eval)
     enddo
  enddo

! #if DEBUG
!   open(ifi,file='evl_xh.allq')
!   do  iw = 1, nwx
!      write(ifi,"( f13.5, 100d15.7)")rymev*freq(iw),(mxevl_xinvh(iq,iw),iq=1,nq)
!   enddo
!   close(ifi)
! #endif
  ! ... Finished finding eigenvalues

  ! --- Pole search of xinvh ---
  ! Liqin: there doesn't seem to be a check whether FM or AFM
  epole = 0
  ! ... Ferromagnetic case
  do  iq = 1, nq
     iwpolf = 1
     !       Coarse search: bracket frequency where evl crosses zero
     do  iw = 1, nwx
        if (freq(iw) >= 0d0 ) then
           if (mxevl_xinvh(iq,iw) < 0d0 .AND. &
                mxevl_xinvh(iq,iw+1) > 0d0 ) then
              iwpolf = iw
              epole = freq(iw)
              exit
           endif
        endif
     enddo

     if (iq == 1) then
        rtmp = mxevl_xinvh(iq,iwpolf)
     elseif (iwpolf /= 1) then ! fine search
        do  ! Liqin ... infinite loop is dangerous
           epole = epole + 1d-7/rydberg
           rtmp  = polinta(epole,freq(iwpolf-1:iwpolf+2), &
                mxevl_xinvh(iq,iwpolf-1:iwpolf+2),4)
           if (rtmp > 0) exit
        enddo
     endif
     epole_fm(iq) = epole
     vpole_fm(iq) = rtmp
  enddo

  ! ... Antiferromagnetic case
  do  iq = 1, nq
     iwpola = 1
     epole = freq(iwpola)
     vpole_af(iq) = mxevl_xinvh(iq,iwpola)
     !       Coarse search: bracket freq where evl crosses zero
     do  iw = 1, nwx
        if (freq(iw) >= 0d0 ) then
           if (mxevl_xinvh(iq,iw) > 0d0 .AND. &
                mxevl_xinvh(iq,iw+1) < 0d0 ) then
              iwpola = iw
              epole = freq(iw)
              exit
           endif
        endif
     enddo

     if (iq == 1) then
        rtmp = mxevl_xinvh(iq,iwpola)
     elseif (iwpola /= 1) then ! fine search
        do
           epole = epole + 1d-7/rydberg
           rtmp  = polinta(epole,freq(iwpola-1:iwpola+2), &
                mxevl_xinvh(iq,iwpola-1:iwpola+2),4)
           if (rtmp < 0) exit
        enddo
     elseif (iwpola == 1) then
        epole = freq(iwpola)
        rtmp = mxevl_xinvh(iq,iwpola)
     endif
     epole_af(iq) = -epole
     vpole_af(iq) = rtmp
  enddo

  write(stdo,*)' Results for pole search of evl(<q|X|q>); q FM pole AFM pole'
  do  iq = 1, nq
     write(stdo,102) iq, rymev*epole_fm(iq), vpole_fm(iq) &
          ,rymev*epole_af(iq), vpole_af(iq)
102  format(3x,i3,3x,f8.2,3x, d15.7, 5x,f8.2,3x,d15.7)
  enddo
  ! ... End of pole search

  ! --- Determine meffi ---
  do  iq = 1, nq
     rtmp1 = epole_fm(iq)
     rtmp2 = epole_af(iq)
     do  i = 1, nmbas
        do  j = 1, nmbas
           do  iw = 1, nwx
              dxidw(i,j,iq,iw)  = &
                                !     .   (xinv(i,j,iq,iw+1) - xinv(i,j,iq,1))/(freq(iw+1)-freq(1))
                   (xinvh(i,j,iq,iw+1) - xinvh(i,j,iq,1))/(freq(iw+1)-freq(1))
              freqm(iw) = 0.5d0*(freq(iw+1) + freq(iw))
           enddo
           meffi(i,j,iq) = &
                polinta(rtmp1,freqm(1:8),dble(dxidw(i,j,iq,1:8)),8) &
                + img* polinta(rtmp1,freqm(1:8),dimag(dxidw(i,j,iq,1:8)),8)
           meffi2(i,j,iq) = &
                polinta(rtmp2,freqm(1:8),dble(dxidw(i,j,iq,1:8)),8) &
                + img* polinta(rtmp2,freqm(1:8),dimag(dxidw(i,j,iq,1:8)),8)
        enddo
     enddo
  enddo

  ! --- Xinvh_w0  dxidw_et -> dxidw_ebar ---
  do  iq = 1, nq
     xinvh_w0(:,:,iq) = xinvh(:,:,iq,1)
     do  i = 1, nmbas ! projected on |e-> instead of |e~>
        do  j = 1, nmbas
           xinvh_w0eb(i,j,iq) = xinvh_w0(i,j,iq) &
                * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
           meffi_eb(i,j,iq) = meffi(i,j,iq) &
                * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
           meffi2_eb(i,j,iq) = meffi2(i,j,iq) &
                * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
           do  iw = 1, nwx
              dxidw_eb(i,j,iq,iw) = dxidw(i,j,iq,iw) &
                   * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
           enddo
        enddo
     enddo
  enddo

  write(stdo,*)' Inverse of effective mag. moment meffi(1,1)'// &
       ' = [xinvh(1,1,iw)-xinvh(1,1,iw=0)]/w'// &
       ' q iw FM pole AFM pole intp FM intp AF'
  do  iq = 1, nq
     write(stdo,  "(3x,i3, 5(2d12.4,1x) )" ) iq,dxidw_eb(1,1,iq,1) &
          ,dxidw_eb(1,1,iq,iwpolf),dxidw_eb(1,1,iq,iwpola) &
          , meffi_eb(1,1,iq), meffi2_eb(1,1,iq)
  enddo

! #if DEBUG
!   open(ifi,file='dxidw.allq')
!   do  iw = 1, nwx
!      write(ifi,"(f8.2,100d15.7)") rymev*freqm(iw),(dxidw_eb(1,1,iq,iw), iq=1,nq)
!   enddo
!   close(ifi)
! #endif

  ! --- Make Jmat file for FM and AFM cases ---
  write(stdo,*)'Writing files to disk: Jmat.allq = J(q,w=intp FM pole), and  '// &
       'Jmat_X2w0 = J(q,w=0)'
  nmx = nmbas
  nev = nmbas
  open(ifi,file='Jmat.allq')
  do  iq = 1, nq
     mxi = 0d0
     do  i = 1, nmbas
        mxi(i,i) = 1d0
     enddo
     allocate(wk(11,nmbas))
     call zhev(nmbas,meffi_eb(1,1,iq),wk, &
          .false.,.true.,nmx,1d99,nev,wk,.false.,-1,eval,evc)
     deallocate(wk)
     oo = 0d0
     do  i = 1, nmbas
        if (eval(i) >= 0) then
           oo(i,i) = 1d0/sqrt(eval(i))
        else
           oo(i,i) = 1d0/csqrt( cmplx(eval(i)) )
        endif
     enddo
     sqm = matmul(evc, matmul(oo, transpose(dconjg(evc))) )
     jjmat = matmul(sqm, matmul(xinvh_w0eb(:,:,iq),sqm))
     jjmat2 = matmul(sqm, matmul(mxi, sqm))
     do  ix = 1, nmbas
        do  iy = 1, nmbas
           jjmat(ix,iy) = jjmat(ix,iy)/sqrt(momsite(ix)*momsite(iy))
        enddo
     enddo
     write(ifi,"('JJMAT: ',3d18.10,3x,255d18.10)") &
          qp(:,iq), ((jjmat(ix,iy), ix=1,nmbas),iy=1,nmbas)
     call zevl(nmbas,jjmat,eval)
  enddo
  close(ifi)
  ! --- Save files: JJMAT: J=X0^(-1)  or J=X^(-1) at (w=0) ---
  open(ifi,file='Jmat_X0w0.allq')
  open(ifi2,file='Jmat_X2w0.allq')
  do  iq = 1, nq
     write(ifi,"('JJMAT: ',3d18.10, 3x, 255d18.10)") &
          qp(:,iq), (( x0inv(ix,iy,iq,1) ,ix=1,nmbas),iy=1,nmbas)
     write(ifi2,"('xinvh_w0eb: ',3d18.10, 3x, 255d18.10)") &
          qp(:,iq), ((  xinvh_w0eb(:,:,iq) ,ix=1,nmbas),iy=1,nmbas)
  enddo
  close(ifi)
  close(ifi2)
end subroutine stonerrsa

subroutine zevl(n,h,eval)
  !- Return eigenvalues of a hermitian matrix h, leaving h unaltered
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n     :dimension of h
  !i   h     :small h = cnu-enu+sqrdel*S^beta*sqrdel
  !o Outputs
  !o   eval  :eigenvalues
  !l Local variables
  !r Remarks
  !u Updates
  !u   08 Feb 09 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: n
  double complex h(n,n)
  double precision :: eval(n)
  ! ... Local parameters
  integer :: nev
  double precision :: xx
  complex(8),allocatable:: hloc(:,:),z(:,:),wk(:,:)

  allocate(wk(n,n),z(n,n),hloc(n,n))
  hloc = h
  call zhevx(n,n,hloc,xx,0,.true.,n,1d99,nev,wk,.false.,eval,n,z)
  deallocate(z,wk,hloc)

end subroutine zevl

!=======================================================================
subroutine stonerpb(nq,nw,nmbas,nbloch, qp,momsite,mmnorm ,emesh,zxq)
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
  implicit none
  integer :: nq,nw,nmbas,nbloch
  real(8) momsite(nmbas), mmnorm(nmbas),emesh(nw),qp(3,nq)
  complex(8) zxq(nbloch,nbloch,nq,nw)
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
  complex(8)  oo(nmbas,nmbas),evc(nmbas,nmbas) &
       ,jjmat(nmbas,nmbas),sqm(nmbas,nmbas)
  integer:: nmx,nev
  nwx = 400
  rymev = rydberg*1d3
  if ( emesh(1) /= 0d0 ) call rx('w(1) /= 0')
  freq=emesh
  jb=1
  open(111, file='X0pbqw.allqb')
  write(111,"( 4i4 )") nq,nmbas,nw
  do iw=1,nw
     write(111,301)  (zxq( jb , jb,iq,iw), iq=1,nq )
  enddo
301 format(1000d23.15)
end subroutine stonerpb

endmodule m_stonerrsa
