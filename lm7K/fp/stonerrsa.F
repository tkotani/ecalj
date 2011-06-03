      subroutine stonerrsa(nq,nw,nmbas,qp,momsite,mmnorm,eiqrm,freq,
     .x0et)
      use m_globalvariables
C- Transverse susceptibility matrix <~e| X |e~> from X^{-1}=X0^{-1}+U
C ----------------------------------------------------------------------
Ci Inputs
Ci   nq     :number of q points for which to calculate X(q,w)
Ci   nw     :number of frequencies where input x0et is given
Ci   nmbas  :number of magnetic sites
Ci   qp     :qp(3,nq)  vector of each q
Ci   momsite:magnetic moment m
Ci   mmnorm :<m|m>
Ci   eiqrm  :eiqrm_i = <~e_i|e^{iqr}> =  <M_i|eiqr>/sqrt(<M_i|M_i>)
Ci   freq   :frequency mesh
Ci   x0et   :<~e|X0|~e>
Co Outputs:
Co   ... finish this
Cl Local variables
Cr Remarks
Cu Updates
Cu   07 Feb 09 (L. Ke) adapted from T. Kotani
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nq,nw,nmbas
      real(8) momsite(nmbas), mmnorm(nmbas),freq(nw),qp(3,nq)
      complex(8) eiqrm(nmbas,nq),x0et(nmbas,nmbas,nq,nw)
C ... Local parameters
      integer fopna,fopnx,i,ifi,ifi2,imax,imin,iq,iw,iwpola,iwpolf,ipr,
     .ix,iy,j,jmax,jmin,nglob,nw_intp,nwx,stdo
      complex(8),allocatable :: qxq_intp(:,:),wk(:,:),
     .x0inv(:,:,:,:),xinv(:,:,:,:),x2et(:,:,:,:),
     .xinvh(:,:,:,:),dxidw(:,:,:,:),dxidw_eb(:,:,:,:)
      real(8),allocatable:: mxevl_xinvh(:,:),mxevl2_xinvh(:,:),e_intp(:)
      real(8) uub(nmbas,nw),uu(nmbas,nw),eval(nmbas)
      real(8) freqm(nw),mmnorm2(nmbas)
      real(8) emin_intp,dw_intp,emax_intp
      complex(8):: cxtmp(nmbas,nmbas),img=(0d0,1d0),qxq(nw,nq),
     .meffi(nmbas,nmbas,nq),xinvh_w0(nmbas,nmbas,nq),
     .meffi_eb(nmbas,nmbas,nq),xinvh_w0eb(nmbas,nmbas,nq),
     .meffi2(nmbas,nmbas,nq),meffi2_eb(nmbas,nmbas,nq)
      real(8) qxq_r(nw,nq),qxq_i(nw,nq),rydberg,polinta
      real(8) qxq_inv_r(nw,nq),qxq_inv_i(nw,nq),
     .epole_fm(nq),epole_af(nq),vpole_fm(nq),vpole_af(nq)
      real(8) rtmp,rtmp1,rtmp2,rymev,epole
Ckino 2011.05 del external,       external :: polinta
C     For calculating JMAT:
      complex(8) oo(nmbas,nmbas),evc(nmbas,nmbas),sqm(nmbas,nmbas),
     .jjmat(nmbas,nmbas),jjmat2(nmbas,nmbas),mxi(nmbas,nmbas)
      integer:: nmx,nev
      parameter (rydberg=13.6058d0)

C --- Setup ---
Changenglob      stdo = nglob('stdo')
      stdo = globalvariables%stdo
      call getpr(ipr)
      nwx = 400
      rymev = rydberg*1d3
      if (freq(1) .ne. 0d0)
     .call rx1('nonzero 1st frequency, w(1) = %g',freq(1))
      allocate(x0inv(nmbas,nmbas,nq,nwx))
      allocate(xinv(nmbas,nmbas,nq,nwx),x2et(nmbas,nmbas,nq,nwx) )
      allocate(xinvh(nmbas,nmbas,nq,nwx),mxevl_xinvh(nq,nwx) )
      allocate(dxidw(nmbas,nmbas,nq,nwx))
      allocate(dxidw_eb(nmbas,nmbas,nq,nwx))
      allocate(mxevl2_xinvh(nq,nwx) )

#if DEBUG
C     open(111, file='x0qw.allqb',form='unformatted')
      ifi = fopnx('x0qw.allqb',2,2+4,-1)
      write (ifi) nq,nmbas,nw,nwx,qp
      write (ifi) freq(1:nw)
      do  iq = 1, nq
        do  iw = 1, nw
          write(ifi) dble(x0et(1:nmbas,1:nmbas,iq,iw))
          write(ifi) dimag(x0et(1:nmbas,1:nmbas,iq,iw))
        enddo
      enddo
#endif

C --- Determine <~e|U(w)|e~> by only using boundary condition at q=0 ---
      x0inv(:,:,1,:) = x0et(:,:,1,:) ! only need x0(q=0,w)
      mmnorm2 = mmnorm**2
      do  iw = 1, nwx
        call matcinv(nmbas, x0inv(:,:,1,iw)) ! x0inv=<~e|x0^-1|~e>
        uub(1:nmbas,iw) =
     .  matmul(dble(x0inv(1:nmbas,1:nmbas,1,iw)),mmnorm(1:nmbas))
        do  i = 1, nmbas
          if (abs(momsite(i)) .gt. 1d-3) then
            uu(i,iw) = freq(iw)*momsite(i)/mmnorm2(i)
     .      - uub(i,iw)/mmnorm(i)
          endif
        enddo
C      write(stdo,"(i4, f13.5, d15.7)") iw,rymev*freq(iw),uu(1,iw)
      enddo

#if DEBUG
C     open(112,file='etUet.allw')
C     ifi = fopnx('etUet.allw',2,2,-1)
      ifi = fopna('etUet',-1,0)
      do  iw = 1, nwx
        write(ifi,'(f21.13,3x,255d18.10)')
     .  freq(iw),(uu(ix,iw),ix=1,nmbas)
      enddo
      call fclose(ifi)
#endif

C --- x2et=<~e|X(w)|e~>  &&  <eiqr|~e><~e|X|~e><eiqr|~e> ---
      call info0(20,1,0,' STONERRSA:  calculate full Xi for each q'//
     .'%N Magnetic moments, by (magnetic) site:')
      if (ipr.ge.20) write(stdo,"(1x,12f8.4)") momsite
C     write(stdo,"( 5x,'q',10x,'qvec' )")
      x0inv = x0et
      do  iq = 1, nq
C     write(stdo,'( 3x,i3,3f7.3 )') iq, qp(1:3,iq)
        do  iw = 1, nwx
          call matcinv(nmbas,x0inv(:,:,iq,iw))
          xinv(:,:,iq,iw) = x0inv(:,:,iq,iw)
          do  ix = 1, nmbas
            xinv(ix,ix,iq,iw) = x0inv(ix,ix,iq,iw) + uu(ix,iw)
          enddo
          cxtmp(:,:) = xinv(:,:,iq,iw)
C       Liqin:  why is this in a do loop?
          do  ix = 1, nmbas
            cxtmp(:,:) = cxtmp(:,:) + img*1d-30
          enddo
C       Liqin: this is poor programming technique
          call matcinv(nmbas,cxtmp(:,:)) !this is full x_+-  !Matrix inversion.
          x2et(:,:,iq,iw) = cxtmp(:,:) ! x2et=<~e|X(w)|e~>
          qxq(iw,iq) = sum( eiqrm(:,iq) ! <eiqr|~e><~e|X|~e><~|eiqr>
     .    *matmul(x2et(:,:,iq,iw),dconjg(eiqrm(:,iq) )))
          qxq_r(iw,iq) = dble(qxq(iw,iq) )
          qxq_i(iw,iq) = dimag(qxq(iw,iq) )

          qxq_inv_r(iw,iq) = dble( 1d0/qxq(iw,iq) ) !for interpolation
          qxq_inv_i(iw,iq) = dimag( 1d0/qxq(iw,iq) )

c      write(stdo,"( i4, f13.5, 2d15.7 )") iw,rymev*freq(iw)
c     .             , qxq_r(iw,iq),qxq_i(iw,iq)
        enddo
      enddo

C --- Interpolate <eiqr|X|eiqr> ---
C     Liqin: these should be passed as arguments
      emin_intp = 0d0; dw_intp = 1d-2; emax_intp = 1000d0
      nw_intp = int((emax_intp-emin_intp)/dw_intp) + 1
      call info(20,0,0,' Make <q|X|q>: for %i energy points; '//
     .'emax = %;5,5d eV',nwx,rydberg*freq(nwx))
      call info5(20,0,0,' Interpolate energy window with '//
     .'emin:dw:emax = %;4d:%;4d:%;4d meV (%i points)',emin_intp,
     .dw_intp,emax_intp,nw_intp,0)

      allocate(qxq_intp(nw_intp,nq),e_intp(nw_intp) )
      do  iq = 1, nq
        do  iw = 1, nw_intp
          if (iq .eq. 1) then
            e_intp(iw) = (emin_intp+(iw-1)*dw_intp )/(rymev)
          endif
          rtmp1 = polinta(e_intp(iw),freq(1:nwx),qxq_inv_r(1:nwx,iq),nwx)
          rtmp2 = polinta(e_intp(iw),freq(1:nwx),qxq_inv_i(1:nwx,iq),nwx)
          qxq_intp(iw,iq) = 1d0/(rtmp1 + img*rtmp2) !for interpolation
        enddo
      enddo

      if (ipr .ge. 30) then
        call info0(30,1,0,' Data for pole search of evl(<q|X|q>)%N'//
     .  '%5fq %13fqxq_r_max %14fqxq_r_min %14fqxq_i_max %14fqxq_i_min')
        do  iq = 1, nq
          jmax=maxloc(-dble(qxq_intp(1:nw_intp,iq)),dim=1)
          jmin=minloc(-dble(qxq_intp(1:nw_intp,iq)),dim=1)
          imax=maxloc(-dimag(qxq_intp(1:nw_intp,iq)),dim=1)
          imin=minloc(-dimag(qxq_intp(1:nw_intp,iq)),dim=1)
          write(stdo,101) iq,
     .    rymev*e_intp(jmax), -dble(qxq_intp(jmax,iq)),
     .    rymev*e_intp(jmin), -dble(qxq_intp(jmin,iq)),
     .    rymev*e_intp(imax), -dimag(qxq_intp(imax,iq)),
     .    rymev*e_intp(imin), -dimag(qxq_intp(imin,iq))
  101     format(3x,i3,1x, f7.2,'(',d15.7,')',1x
     .    ,f7.2,'(' ,d15.7,')',1x,f7.2,'(' ,d15.7,')'
     .    ,1x,f7.2,'(' ,d15.7,')' )
        enddo
      endif

#if DEBUG
C      open(106,file='qxqi.allq')
C      open(107,file='qxqr.allq')
      ifi  = fopna('qxqi',-1,0)
      ifi2 = fopna('qxqr',-1,0)
      do  iw = 1, nw_intp
        write(ifi,"( f13.5, 100d15.7)") rymev*e_intp(iw),
     .  (-dimag(qxq_intp(iw,iq)),iq=1,nq)
        write(ifi2,"( f13.5, 100d15.7)") rymev*e_intp(iw),
     .  (-dble(qxq_intp(iw,iq)),iq=1,nq)
      enddo
      call fclose(ifi); call fclose(ifi2)
#endif
C ... Finished interpolation of <eiqr|X|eiqr>

C --- Eigenvalue of hermitian xinvh ---
      do  iq = 1, nq
        do  iw = 1, nwx
          xinvh(:,:,iq,iw) = .5d0*( xinv(:,:,iq,iw)
     .    + transpose(dconjg(xinv(:,:,iq,iw))) )
          call zevl(nmbas,xinvh(1,1,iq,iw),eval)
          mxevl_xinvh(iq,iw) = maxval(eval)
        enddo
      enddo

#if DEBUG
      ifi = fopnx('evl_xh.allq',2,2,-1)
      do  iw = 1, nwx
        write(ifi,"( f13.5, 100d15.7)") rymev*freq(iw),
     .  (mxevl_xinvh(iq,iw),iq=1,nq)
      enddo
      call fclose(ifi)
#endif
C ... Finished finding eigenvalues

C --- Pole search of xinvh ---
C Liqin: there doesn't seem to be a check whether FM or AFM
      epole = 0
C ... Ferromagnetic case
      do  iq = 1, nq
        iwpolf = 1
C       Coarse search: bracket frequency where evl crosses zero
        do  iw = 1, nwx
          if (freq(iw) .ge. 0d0 ) then
            if (mxevl_xinvh(iq,iw) .lt. 0d0 .and.
     .      mxevl_xinvh(iq,iw+1) .gt. 0d0 ) then
              iwpolf = iw
              epole = freq(iw)
              exit
            endif
          endif
        enddo

        if (iq .eq. 1) then
          rtmp = mxevl_xinvh(iq,iwpolf)
        elseif (iwpolf .ne. 1) then ! fine search
          do  ! Liqin ... infinite loop is dangerous
            epole = epole + 1d-7/rydberg
            rtmp  = polinta(epole,freq(iwpolf-1:iwpolf+2),
     .      mxevl_xinvh(iq,iwpolf-1:iwpolf+2),4)
            if (rtmp .gt. 0) exit
          enddo
        endif
        epole_fm(iq) = epole
        vpole_fm(iq) = rtmp
      enddo

C ... Antiferromagnetic case
      do  iq = 1, nq
        iwpola = 1
        epole = freq(iwpola)
        vpole_af(iq) = mxevl_xinvh(iq,iwpola)
C       Coarse search: bracket freq where evl crosses zero
        do  iw = 1, nwx
          if (freq(iw) .ge. 0d0 ) then
            if (mxevl_xinvh(iq,iw) .gt. 0d0 .and.
     .      mxevl_xinvh(iq,iw+1) .lt. 0d0 ) then
              iwpola = iw
              epole = freq(iw)
              exit
            endif
          endif
        enddo

        if (iq .eq. 1) then
          rtmp = mxevl_xinvh(iq,iwpola)
        elseif (iwpola .ne. 1) then ! fine search
          do
            epole = epole + 1d-7/rydberg
            rtmp  = polinta(epole,freq(iwpola-1:iwpola+2),
     .      mxevl_xinvh(iq,iwpola-1:iwpola+2),4)
            if (rtmp .lt. 0) exit
          enddo
        elseif (iwpola .eq. 1) then
          epole = freq(iwpola)
          rtmp = mxevl_xinvh(iq,iwpola)
        endif
        epole_af(iq) = -epole
        vpole_af(iq) = rtmp
      enddo

      call info0(20,1,0,' Results for pole search of evl(<q|X|q>)%N'//
     .'%5fq %9fFM pole %23fAFM pole')
      do  iq = 1, nq
        write(stdo,102) iq, rymev*epole_fm(iq), vpole_fm(iq)
     .  ,rymev*epole_af(iq), vpole_af(iq)
 102    format(3x,i3,3x,f8.2,3x, d15.7, 5x,f8.2,3x,d15.7)
      enddo
C ... End of pole search

C --- Determine meffi ---
      do  iq = 1, nq
        rtmp1 = epole_fm(iq)
        rtmp2 = epole_af(iq)
        do  i = 1, nmbas
          do  j = 1, nmbas
            do  iw = 1, nwx
              dxidw(i,j,iq,iw)  =
c     .   (xinv(i,j,iq,iw+1) - xinv(i,j,iq,1))/(freq(iw+1)-freq(1))
     .        (xinvh(i,j,iq,iw+1) - xinvh(i,j,iq,1))/(freq(iw+1)-freq(1))
              freqm(iw) = 0.5d0*(freq(iw+1) + freq(iw))
            enddo
            meffi(i,j,iq) =
     .      polinta(rtmp1,freqm(1:8),dble(dxidw(i,j,iq,1:8)),8)
     .      + img* polinta(rtmp1,freqm(1:8),dimag(dxidw(i,j,iq,1:8)),8)
            meffi2(i,j,iq) =
     .      polinta(rtmp2,freqm(1:8),dble(dxidw(i,j,iq,1:8)),8)
     .      + img* polinta(rtmp2,freqm(1:8),dimag(dxidw(i,j,iq,1:8)),8)
          enddo
        enddo
      enddo

C --- Xinvh_w0  dxidw_et -> dxidw_ebar ---
      do  iq = 1, nq
        xinvh_w0(:,:,iq) = xinvh(:,:,iq,1)
        do  i = 1, nmbas ! projected on |e-> instead of |e~>
          do  j = 1, nmbas
            xinvh_w0eb(i,j,iq) = xinvh_w0(i,j,iq)
     .      * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
            meffi_eb(i,j,iq) = meffi(i,j,iq)
     .      * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
            meffi2_eb(i,j,iq) = meffi2(i,j,iq)
     .      * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
            do  iw = 1, nwx
              dxidw_eb(i,j,iq,iw) = dxidw(i,j,iq,iw)
     .        * mmnorm(i)/momsite(i) * mmnorm(j)/momsite(j)
            enddo
          enddo
        enddo
      enddo

      call info0(20,1,0,' Inverse of effective mag. moment meffi(1,1)'//
     .' = [xinvh(1,1,iw)-xinvh(1,1,iw=0)]/w%N'//
     .'%5fq%11fiw=1%21fFM pole%18fAFM pole%17fintp FM%18fintp AF')
      do  iq = 1, nq
        write(stdo,  "(3x,i3, 5(2d12.4,1x) )" ) iq,dxidw_eb(1,1,iq,1)
     .  ,dxidw_eb(1,1,iq,iwpolf),dxidw_eb(1,1,iq,iwpola)
     .  , meffi_eb(1,1,iq), meffi2_eb(1,1,iq)
      enddo

#if DEBUG
C     open(109,file='dxidw.allq')
C     ifi = fopna('dxidw',-1,0)
      ifi = fopnx('dxidw.allq',2,2,-1)
      do  iw = 1, nwx
        write(ifi,"(f8.2,100d15.7)") rymev*freqm(iw),
     .  (dxidw_eb(1,1,iq,iw), iq=1,nq)
      enddo
      call fclose(ifi)
#endif

C --- Make Jmat file for FM and AFM cases ---
      call info0(20,1,0,'Writing files to disk: '//
     .'Jmat.allq = J(q,w=intp FM pole), and  '//
     .'Jmat_X2w0 = J(q,w=0)')
      nmx = nmbas
      nev = nmbas
C     open(110,file='Jmat.allq')
      ifi = fopnx('Jmat.allq',2,2,-1)
      do  iq = 1, nq
C       oo = 0d0
        mxi = 0d0
        do  i = 1, nmbas
C         oo(i,i) = 1d0
          mxi(i,i) = 1d0
        enddo
        allocate(wk(11,nmbas))
        call zhev(nmbas,meffi_eb(1,1,iq),wk,
     .  .false.,.true.,nmx,1d99,nev,wk,.false.,-1,eval,evc)
        deallocate(wk)
C       call diagcv(oo,meffi_eb(1,1,iq),evc, nmbas, eval,nmx,1d99,nev)

        oo = 0d0
        do  i = 1, nmbas
          if (eval(i) .ge. 0) then
            oo(i,i) = 1d0/sqrt(eval(i))
          else
            oo(i,i) = 1d0/csqrt( cmplx(eval(i)) )
          endif
C         oo(i,i)= 1d0/sqrt(eval(i))
        enddo
C       write(stdo,"( 'oo', 100f13.7)") (oo(i,i),i=1,nmbas)
        sqm = matmul(evc, matmul(oo, transpose(dconjg(evc))) )
        jjmat = matmul(sqm, matmul(xinvh_w0eb(:,:,iq),sqm))
        jjmat2 = matmul(sqm, matmul(mxi, sqm))
c       jjmat = xinvh_w0eb(:,:,iq)
        do  ix = 1, nmbas
          do  iy = 1, nmbas
            jjmat(ix,iy) = jjmat(ix,iy)/sqrt(momsite(ix)*momsite(iy))
          enddo
        enddo

C        jjmat2 = xinvh_w0eb(:,:,iq)
C        do  ix = 1, nmbas
C          do  iy = 1, nmbas
C            jjmat2(ix,iy)= jjmat2(ix,iy)/sqrt(momsite(ix)*momsite(iy))
C          enddo
C        enddo

C         open(103,file='Jmat_X3w0.allq')
C         write(103,"('sqmIsqm: ',3d18.10, 3x, 255d18.10)")
C     .     qp(:,iq), (( jjmat2(ix,iy) ,ix=1,nmbas),iy=1,nmbas)


        write(ifi,"('JJMAT: ',3d18.10,3x,255d18.10)")
     .  qp(:,iq), ((jjmat(ix,iy), ix=1,nmbas),iy=1,nmbas)
        call zevl(nmbas,jjmat,eval)

C         write(stdo,"('e sw  ', 3f8.4,2x,255d15.7)")
C     .     qp(:,iq),(-momsite(1)*1d3*rydberg*eval(ix),ix=1,nmbas)
      enddo

C --- Save files: JJMAT: J=X0^(-1)  or J=X^(-1) at (w=0) ---
C     open(104,file='Jmat_X0w0.allq')
C     open(105,file='Jmat_X2w0.allq')
      ifi  = fopnx('Jmat_X0w0.allq',2,2,-1)
      ifi2 = fopnx('Jmat_X2w0.allq',2,2,-1)
      do  iq = 1, nq
        write(ifi,"('JJMAT: ',3d18.10, 3x, 255d18.10)")
     .  qp(:,iq), (( x0inv(ix,iy,iq,1) ,ix=1,nmbas),iy=1,nmbas)
        write(ifi2,"('xinvh_w0eb: ',3d18.10, 3x, 255d18.10)")
     .  qp(:,iq), ((  xinvh_w0eb(:,:,iq) ,ix=1,nmbas),iy=1,nmbas)
      enddo

      end subroutine stonerrsa

      subroutine zevl(n,h,eval)
C- Return eigenvalues of a hermitian matrix h, leaving h unaltered
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :dimension of h
Ci   h     :small h = cnu-enu+sqrdel*S^beta*sqrdel
Co Outputs
Co   eval  :eigenvalues
Cl Local variables
Cr Remarks
Cu Updates
Cu   08 Feb 09 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n
      double complex h(n,n)
      double precision eval(n)
C ... Local parameters
      integer nev
      double precision xx
      complex(8),allocatable:: hloc(:,:),z(:,:),wk(:,:)

      allocate(wk(n,n),z(n,n),hloc(n,n))
      hloc = h
      call zhevx(n,n,hloc,xx,0,.true.,n,1d99,nev,wk,.false.,eval,n,z)
      deallocate(z,wk,hloc)

      end

