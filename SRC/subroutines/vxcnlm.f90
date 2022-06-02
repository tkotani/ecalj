!!= Gradient correction to smoothed rho(q) tabulated on a mesh =
!!*Kotani's version newmode with xcpbe.F in abinit Aug2010
subroutine vxcnlm(lxcg,nsp,k1,k2,k3,smrho,repnl,rmunl,vavgnl,vxnl,vcnl,vxcnl)
  use m_supot,only: iv_a_okv,rv_a_ogv
  !      use m_struc_def, only:s_lat
  use m_xcpbe,  only: xcpbe
  use m_lmfinit,only:    lat_alat
  use m_lattic,only: lat_vol
  use m_supot,only: lat_nabc
  use m_supot,only: lat_ng


  !! ----------------------------------------------------------------------
  ! i Inputs
  ! i   lxcg  : dummy now.  (need to set option in xcpbe)
  ! i   slat,smrho(k1,k2,k3,nsp)
  ! o Outputs (for newmode=T).
  ! o   repnl : integral smrho * eps
  ! o   rmunl : integral smrho * vxc
  ! o   vavgnl:average NL XC potential
  ! o   vxcnl : XC potential on uniform mesh.
  !!   vcnl  : dummy (it was correlation part of vxcnl)
  !!   vxnl  : dummy (it was exchange part of vxcnl)
  !! ----------------------------------------------------------------------

  ! cccccccccccccccccccccccccc
  !  old document below. Kink can exist for (grad |grad rho|) (imagine a case with rho=x^2+1)

  ! cccccccccccccSpecifies GGA for old case
  !i         :  0    LSDA
  !i         :  1    Langreth-Mehl
  !i         :  2    PW91
  !i         :  3    PBE
  !i         :  4    PBE with Becke exchange
  !i   nsp   : 2 for spin-polarized case, otherwise 1
  !i   k1..k3: dimensions of smrho,vnl for smooth mesh density
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: nabc ng ogv okv alat vol
  !i     Stored:
  !i     Passed to: vxcgga vxnlcc vxnloc
  !i   smrho :smooth density on uniform mesh
  !l Local variables :
  !l   agr(*,1)  : |grad rhop| or |grad rho| if nsp=1
  !l   agr(*,2)  : |grad rhom| (nsp=2)
  !l   agr(*,k)  : |grad total rho|. k=3 for nsp=2; else k=1
  !l   agr(*,4)  : grad rho+ . grad rho- (only for Langreth-Mehl-Hu)
  !l   ggr(*,1)  : Laplacian of rhop (total rho if nsp=1)
  !l   ggr(*,2)  : Laplacian of rhom (nsp=2)
  !l   gagr(*,k) : (grad rho).(grad |grad rho|)
  !l   gagr(*,1) : (grad rhop).(grad |grad rhop|) (total rho if nsp=1)
  !l   gagr(*,2) : (grad rhom).(grad |grad rhom|) (nsp=2)
  !l   gagr(*,k) : (grad rho).(grad |grad rho|). k=3 for nsp=2; else k=1
  !r Remarks
  !r
  !u Updates
  !u   06 Apr 09 Adapted from vxcnlp.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: lxcg,k1,k2,k3,nsp
  real(8):: repnl(2),rmunl(2),vavgnl(2)
  complex(8):: smrho(k1,k2,k3,nsp),vxcnl(k1,k2,k3,nsp)
  complex(8):: vxnl(k1,k2,k3,nsp), vcnl(k1,k2,k3,nsp),tpibai

  !      type(s_lat)::  slat
  integer :: ip,i,i1,i2,i3,lcut,ng,np,ipr,n1,n2,n3,ngabc(3),nnn
  real(8),allocatable :: ggr(:,:),agr(:,:),gagr(:,:),rho(:,:)
  real(8),allocatable :: enl(:,:,:),vnl(:,:,:,:),enlbk(:,:,:)
  real(8),allocatable:: dvxcdgr(:,:,:,:),grho2_updn(:,:,:,:),grho(:,:,:,:,:),gv(:,:)
  real(8),allocatable::  grho2_updn_forcall(:,:,:,:)
  complex(8),allocatable:: zgrho(:,:,:),gzgrho(:,:,:),zggrho(:,:)
  complex(8),allocatable:: fgrd(:,:,:),fn(:,:,:),fg(:,:),fgg(:)

  real(8):: alat,vol,xx,fac,tpiba,pi,smmin,sss
  integer ::ig,dummy,j,isp
  ! newmode switch is also in smvxcm
  !      logical:: newmode=.false. , debug=.false., plottest=.false.
  logical:: newmode=.true. , debug=.false., plottest=.false.
  !      logical:: newmode=.true. , debug=.false., plottest=.true.
  !      logical:: newmode=.true. , debug=.false., plottest=.true.

  real(8),allocatable:: r_smrho(:,:,:,:)

  !!== Setup ==
  call tcn('vxcnlm')
  if(debug) print *,'smvxcm:calling vxcnlm sum=',sum(abs(smrho))
  call getpr(ipr)
  ngabc =lat_nabc
  n1= ngabc(1)
  n2= ngabc(2)
  n3= ngabc(3)
  if(abs(n1-k1)+abs(n2-k2)+abs(n3-k3)/=0) call rx('vxcnlm: ni/=ki')
  ng    = lat_ng
  allocate(gv(ng,3))
  call dcopy(3*ng,rv_a_ogv,1,gv,1)
  alat  =lat_alat
  vol   =lat_vol
  np = n1*n2*n3 !k1*k2*k3
  if(debug) print *,'vxcnlm: sum check smrho=',sum(abs(smrho))
  !!== New mode start here ==
  !!== Obtain grho= \nabla smrho (on real mesh) ==
  if(newmode) then
     if(debug) print * ,'goto newmode 111',smrho(1,1,1,1)
     allocate(zgrho(np,3,nsp))
     do  i = 1, nsp
        call grfmsh ( 201 , alat , ng , rv_a_ogv , iv_a_okv , k1 , k2 &
             , k3 , n1 , n2 , n3 , smrho ( 1 , 1 , 1 , i ) , zgrho ( 1 , 1 &
             , i ) , dummy ) !dummy = zzgrho is  not touched for grfmsh(isw=201,...
     enddo
     allocate(grho(k1,k2,k3,3,nsp))
     call dcopy(3*np*nsp,zgrho,2,grho,1) ! grho contains $\nabla$smrho in real space
     deallocate(zgrho)
     !! == grho2_updn = (\nabla smrho) **2 ==
     if(debug) print * ,'newmode 222 nsp=',nsp
     allocate(grho2_updn(n1,n2,n3,2*nsp-1) )
     if(debug) print * ,'newmode 222aaa1'
     if(debug) print *,'newmode333 init sum=',sum(abs(smrho))
     do isp=1,nsp
        do i1=1,n1
           do i2=1,n2
              do i3=1,n3
                 grho2_updn(i1,i2,i3,isp) = sum( grho(i1,i2,i3,:,isp)**2 )
                 if(nsp==2) grho2_updn(i1,i2,i3,3) = &
                      sum(grho(i1,i2,i3,1,:))**2 + sum(grho(i1,i2,i3,2,:))**2 + sum(grho(i1,i2,i3,3,:))**2
              enddo
           enddo
        enddo
     enddo
     if(debug) print * ,'newmode 222aaa2',k1,k2,k3,n1,n2,n3,nsp
     !!== call xcpbe in abinit ==
     if(debug) print * ,'newmode 333 a',k1,k2,k3,nsp
     allocate( vnl(k1,k2,k3,nsp) )
     if(debug) print * ,'newmode 333 c'
     allocate( enl(k1,k2,k3))
     if(debug) print * ,'newmode 333xxx'
     allocate( dvxcdgr(k1,k2,k3,3))
     if(debug) print * ,'newmode 333 b'
     fac=1d0/2d0  !This fac is required since rp (:,:,isp=1) contains total density in the case of nsp=1.
     if(nsp==2) fac=1d0
     if(debug) print * ,'newmode goto xcpbe'
     ! i
     allocate(r_smrho(k1,k2,k3,nsp))
     r_smrho=fac*dreal(smrho)
     ! allocate for calling a subroutine
     allocate(grho2_updn_forcall(n1,n2,n3,2*nsp-1))
     do isp=1,2*nsp-1
        do i3=1,n3
           do i2=1,n2
              do i1=1,n1
                 grho2_updn_forcall(i1,i2,i3,isp)=fac**2*grho2_updn(i1,i2,i3,isp)
              enddo
           enddo
        enddo
     enddo
     call xcpbe(exci=enl,npts=n1*n2*n3,nspden=nsp, &
          option=2,&!  Choice of the functional =2:PBE-GGA
          order=1, &!  order=1 means we only calculate first derivative of rho*exc(rho,\nable rho).
                !  rho_updn=fac*dreal(smrho),vxci=vnl,ndvxci=0,ngr2=2*nsp-1,nd2vxci=0,  !Mandatory Arguments
          rho_updn=r_smrho,vxci=vnl,ndvxci=0,ngr2=2*nsp-1,nd2vxci=0,&  ! & Mandatory Arguments
                !  dvxcdgr=dvxcdgr, grho2_updn=fac**2*grho2_updn)   !Optional Arguments
          dvxcdgr=dvxcdgr, grho2_updn=grho2_updn_forcall)   !Optional Arguments
     deallocate(grho2_updn_forcall) ! deallocate temporary data
     deallocate(r_smrho)

     !!=== Output: converted to Ry.===
     enl = 2d0*enl !in Ry.
     vnl = 2d0*vnl !in Ry.
     dvxcdgr= 2d0*dvxcdgr !in Ry.
     if(debug) print * ,'newmode 333 111 end of xcpbe'


     !!== vxcnl is given ==
     fac=1d0/2d0
     if(nsp==2) fac=1d0
     allocate(fg(ng,3),fn(k1,k2,k3), fgg(ng) )
     allocate(fgrd(k1,k2,k3))
     do isp=1,nsp
        do j=1,3 !x,y,z components of grho.
           fn(:,:,:) = fac*grho(:,:,:,j,isp)*dvxcdgr(:,:,:,isp)      !exchange part for spin density
           do i1=1,n1
              do i2=1,n2
                 do i3=1,n3
                    fn(i1,i2,i3) = fn(i1,i2,i3) &
                         + sum(grho(i1,i2,i3,j,1:nsp)) * dvxcdgr(i1,i2,i3,3) !correlation part for total density
                    ! sum(grho(:,:,:,j,1:nsp),dim=5) means that a sum only for 1:nsp.-->this cause compile error in gfortran
                    ! n (complex)
                    !          print *,'sumcheck fn j =',sum(abs(fn))
                 enddo
              enddo
           enddo
           call fftz3(fn,n1,n2,n3,k1,k2,k3,1,0,-1)  ! fn (complex) is converted from real space to reciprocal space
           call gvgetf(ng,1,iv_a_okv,k1,k2,k3,fn,fg(1,j))
           !          print *,'sumcheck fg j =',sum(abs(fg(:,j)))
        enddo
        !!== make i G fg ---> FFT back ==
        if(debug) print * ,'newmode 444'
        !         print *, 'check xxx fg(1,1:3) sum =',sum(abs(gv(1,1:3) *fg(1,1:3)))
        pi = 4d0*datan(1d0)
        tpiba = 2d0*pi/alat
        tpibai = dcmplx(0d0,1d0)*tpiba
        do ig=1,ng
           fgg(ig) =  tpibai*sum( gv(ig,1:3) * fg(ig,1:3))
           !           print *, ig,abs(fgg(ig)),sum(gv(ig,1:3)**2),sum(fg(ig,1:3)**2)
        enddo
        !         print *, 'check xxx fgg sum =',sum(abs(fgg)),ng,sum(abs(gv)),sum(abs(fg))
        call gvputf(ng,1,iv_a_okv,k1,k2,k3,fgg,fgrd)
        call fftz3(fgrd,n1,n2,n3,k1,k2,k3,1,0,1)
        vxcnl(:,:,:,isp) = vnl(:,:,:,isp) - dreal(fgrd)
     enddo

     !!=== plottest check write for debug ===
     if(plottest) then
        isp=1
        do i1=1,1
           do i2=1,n2
              do i3=1,n3
                 write(8006,"(3i4,10e12.4)") i1,i2,i3,vxcnl(i1,i2,i3,isp) ,fgrd(i1,i2,i3)
                 write(9006,"(3i4,10e12.4)") i1,i2,i3,enl(i1,i2,i3)
              enddo
              write(8006,*)
              write(9006,*)
           enddo
        enddo
        !        stop ' test end of plot 8006'
     endif
     deallocate(fgrd)
     deallocate(fn,fgg,fg)

     if(debug) print * ,'newmode 555'

     !!=== vxnl and vcnl are dummy now ===
     vxnl=0d0 !dummy now
     vcnl=0d0 !dummy now

     !!== Make reps, rmu ==
     repnl=0d0
     rmunl=0d0
     vavgnl=0d0
     do  i = 1, nsp
        do i3 = 1, n3
           do i2 = 1, n2
              do i1 = 1, n1
                 repnl(i) = repnl(i) + dble(smrho(i1,i2,i3,i))*enl(i1,i2,i3)
                 rmunl(i) = rmunl(i) + dble(smrho(i1,i2,i3,i))*vxcnl(i1,i2,i3,i) !all total
                 vavgnl(i) = vavgnl(i) + vxcnl(i1,i2,i3,i)                       !all total
              enddo
           enddo
        enddo
        repnl(i)  = repnl(i)*vol/(n1*n2*n3)
        rmunl(i)  = rmunl(i)*vol/(n1*n2*n3)
        vavgnl(i) = vavgnl(i)/(n1*n2*n3)
     enddo
     if(plottest) then
        allocate(enlbk(k1,k2,k3))
        enlbk=enl
     endif
     if(debug) print * ,'newmode 666'

     deallocate(grho,grho2_updn,dvxcdgr,vnl,enl)

     call tcx('vxcnlm')
     if( .NOT. plottest) return
  endif
  !!* This is the end of newmode.
  !!--------------------------------------------------------------

  !! == From here on, orinal mode by Mark ==
  allocate(agr(np,3*nsp-2),gagr(np,2*nsp-1),ggr(np,nsp))
  allocate(zgrho(np,3,nsp),gzgrho(np,3,2*nsp-1),zggrho(np,nsp))
  ! --- Grad rho_i and Laplacian rho_i (complex) ---
  do  i = 1, nsp
     call grfmsh ( 601 , alat , ng , rv_a_ogv , iv_a_okv , k1 , k2 &
          , k3 , n1 , n2 , n3 , smrho ( 1 , 1 , 1 , i ) , zgrho ( 1 , 1 &
          , i ) , zggrho ( 1 , i ) )
  enddo

  ! --- agr_i : |grad rho_i|, i=1,2 and agr_i(3) : |grad rho| ---
  !     and ggr_i = lap rho_i.  Also agr(4) : grad rho+ . grad rho-
  do  i = 1, nsp
     do  ip = 1, np
        agr(ip,i) = dsqrt(dble(zgrho(ip,1,i))**2 + &
             dble(zgrho(ip,2,i))**2 + &
             dble(zgrho(ip,3,i))**2)
        ggr(ip,i) = dble(zggrho(ip,i))
     enddo
  enddo
  if (nsp == 2) then
     do  ip = 1, np
        agr(ip,3) = dsqrt(dble(zgrho(ip,1,1)+zgrho(ip,1,2))**2 + &
             dble(zgrho(ip,2,1)+zgrho(ip,2,2))**2 + &
             dble(zgrho(ip,3,1)+zgrho(ip,3,2))**2)
        agr(ip,4) =       dble(zgrho(ip,1,1)*zgrho(ip,1,2)) + &
             dble(zgrho(ip,2,1)*zgrho(ip,2,2)) + &
             dble(zgrho(ip,3,1)*zgrho(ip,3,2))
     enddo
  endif
  !     do  i = 1, 3*nsp-2
  !       call prm3('|grad rho(isp=%i)|',i,agr(1,i),n1,n2,n3)
  !     enddo

  ! --- gzgrho (complex) : grad |grad rho_i|, i=1,2,3 (see above for i=3) ---
  !     Use zggrho as complex work array
  do  i = 1, 2*nsp-1
     call dpzero(zggrho,np*2)
     call dcopy(np,agr(1,i),1,zggrho,2)
     !       call zprm3('|grad rho_i|',0,zggrho(1,i),n1,n2,n3)
     call grfmsh ( 201 , alat , ng , rv_a_ogv , iv_a_okv , k1 , k2 &
          , k3 , n1 , n2 , n3 , zggrho , gzgrho ( 1 , 1 , i ) , xx )
  enddo

  deallocate(zggrho)

  ! --- gagr : grad rho_i . grad |grad rho_i|, i=1,2,3 (see above for i=3) ---
  do  i = 1, nsp
     do  ip = 1, np
        gagr(ip,i) = &
             dble(zgrho(ip,1,i))*dble(gzgrho(ip,1,i)) + &
             dble(zgrho(ip,2,i))*dble(gzgrho(ip,2,i)) + &
             dble(zgrho(ip,3,i))*dble(gzgrho(ip,3,i))
     enddo
     !       call prm3('grad rho . grad |grad rho_%i|',i,gagr(1,i),n1,n2,n3)
  enddo
  if (nsp == 2) then
     do  ip = 1, np
        gagr(ip,3) = &
             dble(zgrho(ip,1,1)+zgrho(ip,1,2))*dble(gzgrho(ip,1,3)) + &
             dble(zgrho(ip,2,1)+zgrho(ip,2,2))*dble(gzgrho(ip,2,3)) + &
             dble(zgrho(ip,3,1)+zgrho(ip,3,2))*dble(gzgrho(ip,3,3))
     enddo
     !       call prm3('grad rho . grad |grad rho_%i|',3,gagr(1,3),n1,n2,n3)
  endif

  deallocate(zgrho,gzgrho)

  ! --- Nonlocal potential for all points  ---
  allocate(vnl(k1,k2,k3,nsp),enl(k1,k2,k3),rho(np,nsp))
  call dpzero(vnl,np*nsp)
  call dpzero(enl,np)
  do  i = 1, nsp
     call dcopy(np,smrho(1,1,1,i),2,rho(1,i),1)
     !       call zprm3('smrho_%i',i,smrho(1,1,1,i),n1,n2,n3)
     !       call prm3 ('rho_%i',i,rho(1,i),n1,n2,n3)
     !       call prm3 ('lap-rho_%i',i,ggr(1,i),n1,n2,n3)
  enddo

  print *,'vxcnlm lxcg =',lxcg
  if (lxcg > 2) then
     call vxcgga(lxcg,np,nsp,rho,rho(1,nsp),agr(1,1),agr(1,nsp), &
          ggr(1,1),ggr(1,nsp),agr(1,2*nsp-1),agr(1,4), &
          gagr(1,2*nsp-1),gagr(1,1),gagr(1,nsp),vnl,vnl(1,1,1,nsp),enl)
  else
     lcut = 1
     if (lcut == 1) then
        call vxnlcc(np,nsp,rho,rho(1,nsp),agr(1,1),agr(1,nsp), &
             ggr(1,1),ggr(1,nsp),agr(1,2*nsp-1),agr(1,4),gagr(1,2*nsp-1), &
             gagr(1,1),gagr(1,nsp),vnl,vnl(1,1,1,nsp),enl)
     else
        call vxnloc(np,nsp,rho,rho(1,nsp),agr(1,1),agr(1,nsp), &
             ggr(1,1),ggr(1,nsp),agr(1,2*nsp-1),agr(1,4),gagr(1,2*nsp-1), &
             gagr(1,1),gagr(1,nsp),vnl,vnl(1,1,1,nsp),enl)
     endif
  endif

  !      call prm3('enl',i,enl,n1,n2,n3)
  !      do  i = 1, nsp
  !        call prm3('vnl(isp=%i)',i,vnl(1,1,1,i),n1,n2,n3)
  !      enddo


  ! --- Make nonlocal reps, rmu ---
  repnl=0d0
  rmunl=0d0
  vavgnl=0d0
  do  i = 1, nsp
     do  i3 = 1, n3
        do  i2 = 1, n2
           do  i1 = 1, n1
              repnl(i) = repnl(i) + dble(smrho(i1,i2,i3,i))*enl(i1,i2,i3)
              rmunl(i) = rmunl(i) + dble(smrho(i1,i2,i3,i))*vnl(i1,i2,i3,i)
              vxcnl(i1,i2,i3,i) = vxcnl(i1,i2,i3,i) + vnl(i1,i2,i3,i)
              vavgnl(i) = vavgnl(i) + vnl(i1,i2,i3,i)
           enddo
        enddo
     enddo
     repnl(i)  = repnl(i)*vol/(n1*n2*n3)
     rmunl(i)  = rmunl(i)*vol/(n1*n2*n3)
     vavgnl(i) = vavgnl(i)/(n1*n2*n3)
  enddo

  ! ccccccccccccccccccccccccccccc
  if(plottest) then
     isp=1
     !        smmin=0d0
     do i1=1,n1
        do i2=1,n2
           do i3=1,n3
              write(8007,"(3i4,10e14.6)") i1,i2,i3,vxcnl(i1,i2,i3,isp)
              write(9007,"(3i4,10e14.6)") i1,i2,i3,enl(i1,i2,i3) ,enlbk(i1,i2,i3) ,enl(i1,i2,i3)-enlbk(i1,i2,i3)
              !          sss=abs(enl(i1,i2,i3)-enlbk(i1,i2,i3) )
              !          if(sss>smmin)smmin=sss
           enddo
           write(8007,*)
           write(9007,*)
        enddo
     enddo
     !        print *,' max error=',smmin
     !        stop ' test end of plot 8007'
  endif
  ! ccccccccccccccccccccccccccccc
  deallocate(rho,agr,gagr,ggr,enl,vnl)
  call tcx('vxcnlm')
  stop 'vvvvvvvvvvvvv test end xxxxxxxxxxxxxxxxxxxxx qqqqqqq'
end subroutine vxcnlm
