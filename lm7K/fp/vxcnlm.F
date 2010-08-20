!!= Gradient correction to smoothed rho(q) tabulated on a mesh =
!!*Kotani's version newmode with xcpbe.F in abinit Aug2010
      subroutine vxcnlm(lxcg,nsp,k1,k2,k3,slat,smrho,
     .repnl,rmunl,vavgnl,vxnl,vcnl,vxcnl)
      use m_struc_def, only:s_lat
c      use m_xcpbe_interface
!! ----------------------------------------------------------------------
!!i Inputs
!!i   lxcg  : dummy now.  (need to set option in xcpbe)
!!i   slat,smrho(k1,k2,k3,nsp)
!!o Outputs (for newmode=T).
!!o   repnl : integral smrho * eps
!!o   rmunl : integral smrho * vxc
!!o   vavgnl:average NL XC potential
!!o   vxcnl : XC potential on uniform mesh.
!!   vcnl  : dummy (it was correlation part of vxcnl)
!!   vxnl  : dummy (it was exchange part of vxcnl)
!! ----------------------------------------------------------------------

cccccccccccccccccccccccccccc
c  old document below. Kink can exist for (grad |grad rho|) (imagine a case with rho=x^2+1)
c  
cccccccccccccccSpecifies GGA for old case
Ci         :  0    LSDA
Ci         :  1    Langreth-Mehl
Ci         :  2    PW91
Ci         :  3    PBE
Ci         :  4    PBE with Becke exchange
Ci   nsp   : 2 for spin-polarized case, otherwise 1
Ci   k1..k3: dimensions of smrho,vnl for smooth mesh density
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: nabc ng ogv okv alat vol
Ci     Stored:
Ci     Passed to: vxcgga vxnlcc vxnloc
Ci   smrho :smooth density on uniform mesh
Cl Local variables :
Cl   agr(*,1)  : |grad rhop| or |grad rho| if nsp=1
Cl   agr(*,2)  : |grad rhom| (nsp=2)
Cl   agr(*,k)  : |grad total rho|. k=3 for nsp=2; else k=1
Cl   agr(*,4)  : grad rho+ . grad rho- (only for Langreth-Mehl-Hu)
Cl   ggr(*,1)  : Laplacian of rhop (total rho if nsp=1)
Cl   ggr(*,2)  : Laplacian of rhom (nsp=2)
Cl   gagr(*,k) : (grad rho).(grad |grad rho|)
Cl   gagr(*,1) : (grad rhop).(grad |grad rhop|) (total rho if nsp=1)
Cl   gagr(*,2) : (grad rhom).(grad |grad rhom|) (nsp=2)
Cl   gagr(*,k) : (grad rho).(grad |grad rho|). k=3 for nsp=2; else k=1
Cr Remarks
Cr
Cu Updates
Cu   06 Apr 09 Adapted from vxcnlp.f
C ----------------------------------------------------------------------
      implicit none
      integer lxcg,k1,k2,k3,nsp
      real(8):: repnl(2),rmunl(2),vavgnl(2)
      complex(8):: smrho(k1,k2,k3,nsp),vxcnl(k1,k2,k3,nsp)
      complex(8):: vxnl(k1,k2,k3,nsp), vcnl(k1,k2,k3,nsp),tpibai

      type(s_lat)::  slat
      integer ip,i,i1,i2,i3,lcut,ng,np,ipr,n1,n2,n3,ngabc(3),nnn
      integer,pointer :: iv_p_okv(:)
      real(8),pointer :: rv_p_ogv(:)

      real(8),allocatable :: ggr(:,:),agr(:,:),gagr(:,:),rho(:,:)
      real(8),allocatable :: enl(:,:,:),vnl(:,:,:,:),enlbk(:,:,:)
      real(8),allocatable:: dvxcdgr(:,:,:,:),grho2_updn(:,:,:,:),grho(:,:,:,:,:),gv(:,:)
      complex(8),allocatable:: zgrho(:,:,:),gzgrho(:,:,:),zggrho(:,:)
      complex(8),allocatable:: fgrd(:,:,:),fn(:,:,:),fg(:,:),fgg(:)

      real(8):: alat,vol,xx,fac,tpiba,pi,smmin,sss
      integer ::ig,dummy,j,isp
c newmode switch is also in smvxcm
c      logical:: newmode=.false. , debug=.false., plottest=.false.
      logical:: newmode=.true. , debug=.false., plottest=.false.
c      logical:: newmode=.true. , debug=.false., plottest=.true.
c      logical:: newmode=.true. , debug=.false., plottest=.true.

      interface
         subroutine xcpbe(exci,npts,nspden,option,order,rho_updn,vxci,ndvxci,ngr2,nd2vxci, !Mandatory Arguments
     &   d2vxci,dvxcdgr,dvxci,exexch,grho2_updn) !Optional Arguments
         use defs_basis
         implicit none
         integer,intent(in),optional :: exexch
         integer,intent(in) :: nd2vxci
         integer,intent(in) :: ndvxci
         integer,intent(in) :: ngr2
         integer,intent(in) :: npts
         integer,intent(in) :: nspden
         integer,intent(in) :: option
         integer,intent(in) :: order
         real(dp),intent(out),optional :: d2vxci(npts,nd2vxci)
         real(dp),intent(out),optional :: dvxcdgr(npts,3)
         real(dp),intent(out),optional :: dvxci(npts,ndvxci)
         real(dp),intent(out) :: exci(npts)
         real(dp),intent(in),optional :: grho2_updn(npts,ngr2)
         real(dp),intent(in) :: rho_updn(npts,nspden)
         real(dp),intent(out) :: vxci(npts,nspden)
         end subroutine xcpbe
      end interface

!!== Setup ==
      call tcn('vxcnlm')
      call getpr(ipr)
      ngabc =slat%nabc
      n1= ngabc(1)
      n2= ngabc(2)
      n3= ngabc(3)
      if(abs(n1-k1)+abs(n2-k2)+abs(n3-k3)/=0) call rx('vxcnlm: ni/=ki')
      ng    = slat%ng
      rv_p_ogv => slat%rv_p_ogv 
      allocate(gv(ng,3))
      call dcopy(3*ng,rv_p_ogv,1,gv,1)
      iv_p_okv => slat%iv_p_okv 
      alat  =slat%alat
      vol   =slat%vol
      np = n1*n2*n3 !k1*k2*k3

!!== New mode start here ==
!!== Obtain grho= \nabla smrho (on real mesh) ==
      if(newmode) then
       if(debug) print * ,'goto newmode 111'
       allocate(zgrho(np,3,nsp))
       do  i = 1, nsp
          call grfmsh ( 201 , alat , ng , rv_p_ogv , iv_p_okv , k1 , k2 
     .    , k3 , n1 , n2 , n3 , smrho ( 1 , 1 , 1 , i ) , zgrho ( 1 , 1
     .    , i ) , dummy ) !dummy = zzgrho is  not touched for grfmsh(isw=201,...
       enddo
       allocate(grho(k1,k2,k3,3,nsp)) 
       call dcopy(3*np*nsp,zgrho,2,grho,1) ! grho contains $\nabla$smrho in real space
       deallocate(zgrho)

!! == grho2_updn = (\nabla smrho) **2 ==
       if(debug) print * ,'newmode 222'
       allocate(grho2_updn(n1,n2,n3,nsp) )
       do isp=1,nsp
        do i1=1,n1
        do i2=1,n2
        do i3=1,n3
          grho2_updn(i1,i2,i3,isp) = sum( grho(i1,i2,i3,:,isp)**2 )
          if(nsp==2) grho2_updn(i1,i2,i3,3) = 
     &    sum(grho(i1,i2,i3,1,:))**2 + sum(grho(i1,i2,i3,2,:))**2 + sum(grho(i1,i2,i3,3,:))**2
        enddo  
        enddo  
        enddo  
       enddo

ccccccccccccccccccccccccccccc
c      grho2_updn=0d0
cccccccccccccccccccccccccccccc

!!== negative smrho check==
!! We need full mesh method(all G numner is needed) or so,
!! to assure positive definite, xcpbe seems to interpret abs(rho) when rho is negative.
!!  ---> see xcpbe.F90 (especially, invcb).
       do isp=1,nsp
        nnn=0
        smmin=1d0
        do i1=1,n1
        do i2=1,n2
        do i3=1,n3
          sss=dreal(smrho(i1,i2,i3,isp))
          if(sss<0d0) then
            nnn=nnn+1
            if(sss<smmin) then
              smmin=sss
            endif
          endif
        enddo
        enddo
        enddo
        if(nnn>0) then
          write(6,*) 'warning:vxcnlm negative smrho; isp number min=',isp,nnn,smmin
        endif
       enddo

!!== call xcpbe in abinit ==
       if(debug) print * ,'newmode 333'
       allocate(dvxcdgr(k1,k2,k3,3), vnl(k1,k2,k3,nsp), enl(k1,k2,k3))
       fac=1d0/2d0  !This fac is required since rp (:,:,isp=1) contains total density in the case of nsp=1.
       if(nsp==2) fac=1d0
       call xcpbe(exci=enl,npts=n1*n2*n3,nspden=nsp,
     &           option=2,!Choice of the functional =2:PBE-GGA
     &           order=1, !order=1 means we only calculate first derivative of rho*exc(rho,\nable rho).
     &           rho_updn=fac*dreal(smrho),vxci=vnl,ndvxci=0,ngr2=2*nsp-1,nd2vxci=0,  !Mandatory Arguments
     &           dvxcdgr=dvxcdgr, grho2_updn=fac**2*grho2_updn)   !Optional Arguments
!!=== Output: converted to Ry.===
       enl = 2d0*enl !in Ry.
       vnl = 2d0*vnl !in Ry.
       dvxcdgr= 2d0*dvxcdgr !in Ry.
       if(debug) print * ,'newmode 333 111 end of xcpbe'

       fac=1d0/2d0
       if(nsp==2) fac=1d0
       allocate(fg(ng,3),fn(k1,k2,k3))

       do isp=1,nsp
         do j=1,3 !x,y,z components of grho.
          fn(:,:,:) = fac*grho(:,:,:,j,isp)*dvxcdgr(:,:,:,isp)      !exchange part for spin density
     &          + sum(grho(:,:,:,j,1:nsp),dim=4) * dvxcdgr(:,:,:,3) !correlation part for total density
          ! sum(grho(:,:,:,j,1:nsp),dim=4) means that a sum only for 1:nsp.
          !fn (complex)
c          print *,'sumcheck fn j =',sum(abs(fn))
          call fftz3(fn,n1,n2,n3,k1,k2,k3,1,0,-1)  ! fn (complex) is converted from real space to reciprocal space
          call gvgetf(ng,1,iv_p_okv,k1,k2,k3,fn,fg(1,j))
c          print *,'sumcheck fg j =',sum(abs(fg(:,j)))
         enddo  
         deallocate(fn)

!!== make i G fg ---> FFT back ==
         if(debug) print * ,'newmode 444'
c         print *, 'check xxx fg(1,1:3) sum =',sum(abs(gv(1,1:3) *fg(1,1:3)))
         pi = 4d0*datan(1d0)
         tpiba = 2d0*pi/alat
         tpibai = dcmplx(0d0,1d0)*tpiba
         allocate( fgg(ng) )
         do ig=1,ng
           fgg(ig) =  tpibai*sum( gv(ig,1:3) * fg(ig,1:3))
c           print *, ig,abs(fgg(ig)),sum(gv(ig,1:3)**2),sum(fg(ig,1:3)**2)
         enddo
c         print *, 'check xxx fgg sum =',sum(abs(fgg)),ng,sum(abs(gv)),sum(abs(fg))
         deallocate(fg)
         allocate(fgrd(k1,k2,k3))
         call gvputf(ng,1,iv_p_okv,k1,k2,k3,fgg,fgrd)
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
c        stop ' test end of plot 8006'
        endif

        deallocate(fgg,fgrd)


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
        deallocate(grho,grho2_updn,dvxcdgr,vnl,enl)

        call tcx('vxcnlm')
        if(.not. plottest) return
      endif
!!* This is the end of newmode.
!!--------------------------------------------------------------





!! == From here on, orinal mode by Mark ==
      allocate(agr(np,3*nsp-2),gagr(np,2*nsp-1),ggr(np,nsp))
      allocate(zgrho(np,3,nsp),gzgrho(np,3,2*nsp-1),zggrho(np,nsp))
C --- Grad rho_i and Laplacian rho_i (complex) ---
c      print *,'vxcnlm xxx 000000',slat%nabc
c      print *,'vxcnlm xxx 000111',n1,n2,n3,k1,k2,k3
      do  i = 1, nsp
C       call zprm3('smrho(isp=%i)',i,smrho(1,1,1,i),n1,n2,n3)
Cchp1         call grfmsh(601,alat,ng,w(ogv),w(okv),k1,k2,k3,n1,n2,n3,
Cchp1      .  smrho(1,1,1,i),zgrho(1,1,i),zggrho(1,i))
Cchp1          call grfmsh ( 601 , alat , ng , rv_p_ogv , w ( okv ) , k1 , k2 
Cchp1      .   , k3 , n1 , n2 , n3 , smrho ( 1 , 1 , 1 , i ) , zgrho ( 1 , 1 
Cchp1      .   , i ) , zggrho ( 1 , i ) ) 
          call grfmsh ( 601 , alat , ng , rv_p_ogv , iv_p_okv , k1 , k2 
     .    , k3 , n1 , n2 , n3 , smrho ( 1 , 1 , 1 , i ) , zgrho ( 1 , 1
     .    , i ) , zggrho ( 1 , i ) ) 
C       call zprm3('gradx smrho(isp=%i)',i,zgrho(1,1,i),n1,n2,n3)
C       call zprm3('grady smrho(isp=%i)',i,zgrho(1,2,i),n1,n2,n3)
C       call zprm3('gradz smrho(isp=%i)',i,zgrho(1,3,i),n1,n2,n3)
C       call zprm3('lap smrho(isp=%i)',i,zggrho(1,i),n1,n2,n3)
      enddo


C --- agr_i : |grad rho_i|, i=1,2 and agr_i(3) : |grad rho| ---
C     and ggr_i = lap rho_i.  Also agr(4) : grad rho+ . grad rho-
      do  i = 1, nsp
        do  ip = 1, np
          agr(ip,i) = dsqrt(dble(zgrho(ip,1,i))**2 +
     .    dble(zgrho(ip,2,i))**2 +
     .    dble(zgrho(ip,3,i))**2)
          ggr(ip,i) = dble(zggrho(ip,i))
        enddo
      enddo
      if (nsp .eq. 2) then
        do  ip = 1, np
          agr(ip,3) = dsqrt(dble(zgrho(ip,1,1)+zgrho(ip,1,2))**2 +
     .    dble(zgrho(ip,2,1)+zgrho(ip,2,2))**2 +
     .    dble(zgrho(ip,3,1)+zgrho(ip,3,2))**2)
          agr(ip,4) =       dble(zgrho(ip,1,1)*zgrho(ip,1,2)) +
     .    dble(zgrho(ip,2,1)*zgrho(ip,2,2)) +
     .    dble(zgrho(ip,3,1)*zgrho(ip,3,2))
        enddo
      endif
C     do  i = 1, 3*nsp-2
C       call prm3('|grad rho(isp=%i)|',i,agr(1,i),n1,n2,n3)
C     enddo

C --- gzgrho (complex) : grad |grad rho_i|, i=1,2,3 (see above for i=3) ---
C     Use zggrho as complex work array
      do  i = 1, 2*nsp-1
        call dpzero(zggrho,np*2)
        call dcopy(np,agr(1,i),1,zggrho,2)
C       call zprm3('|grad rho_i|',0,zggrho(1,i),n1,n2,n3)
Cchp1         call grfmsh(201,alat,ng,w(ogv),w(okv),k1,k2,k3,n1,n2,n3,
Cchp1      .  zggrho,gzgrho(1,1,i),xx)
Cchp1          call grfmsh ( 201 , alat , ng , rv_p_ogv , w ( okv ) , k1 , k2 
Cchp1      .   , k3 , n1 , n2 , n3 , zggrho , gzgrho ( 1 , 1 , i ) , xx ) 
          call grfmsh ( 201 , alat , ng , rv_p_ogv , iv_p_okv , k1 , k2 
     .    , k3 , n1 , n2 , n3 , zggrho , gzgrho ( 1 , 1 , i ) , xx ) 



C        call zprm3('gradx |grad rho_%i|',i,gzgrho(1,1,i),n1,n2,n3)
C        call zprm3('grady |grad rho_%i|',i,gzgrho(1,2,i),n1,n2,n3)
C        call zprm3('gradz |grad rho_%i|',i,gzgrho(1,3,i),n1,n2,n3)
      enddo

      deallocate(zggrho)

C --- gagr : grad rho_i . grad |grad rho_i|, i=1,2,3 (see above for i=3) ---
      do  i = 1, nsp
        do  ip = 1, np
          gagr(ip,i) =
     .    dble(zgrho(ip,1,i))*dble(gzgrho(ip,1,i)) +
     .    dble(zgrho(ip,2,i))*dble(gzgrho(ip,2,i)) +
     .    dble(zgrho(ip,3,i))*dble(gzgrho(ip,3,i))
        enddo
C       call prm3('grad rho . grad |grad rho_%i|',i,gagr(1,i),n1,n2,n3)
      enddo
      if (nsp .eq. 2) then
        do  ip = 1, np
          gagr(ip,3) =
     .    dble(zgrho(ip,1,1)+zgrho(ip,1,2))*dble(gzgrho(ip,1,3)) +
     .    dble(zgrho(ip,2,1)+zgrho(ip,2,2))*dble(gzgrho(ip,2,3)) +
     .    dble(zgrho(ip,3,1)+zgrho(ip,3,2))*dble(gzgrho(ip,3,3))
        enddo
C       call prm3('grad rho . grad |grad rho_%i|',3,gagr(1,3),n1,n2,n3)
      endif

      deallocate(zgrho,gzgrho)

C --- Nonlocal potential for all points  ---
      allocate(vnl(k1,k2,k3,nsp),enl(k1,k2,k3),rho(np,nsp))
      call dpzero(vnl,np*nsp)
      call dpzero(enl,np)
      do  i = 1, nsp
        call dcopy(np,smrho(1,1,1,i),2,rho(1,i),1)
C       call zprm3('smrho_%i',i,smrho(1,1,1,i),n1,n2,n3)
C       call prm3 ('rho_%i',i,rho(1,i),n1,n2,n3)
C       call prm3 ('lap-rho_%i',i,ggr(1,i),n1,n2,n3)
      enddo

      print *,'vxcnlm lxcg =',lxcg
      if (lxcg .gt. 2) then
        call vxcgga(lxcg,np,nsp,rho,rho(1,nsp),agr(1,1),agr(1,nsp),
     .  ggr(1,1),ggr(1,nsp),agr(1,2*nsp-1),agr(1,4),
     .  gagr(1,2*nsp-1),gagr(1,1),gagr(1,nsp),vnl,vnl(1,1,1,nsp),enl)
      else
        lcut = 1
        if (lcut .eq. 1) then
          call vxnlcc(np,nsp,rho,rho(1,nsp),agr(1,1),agr(1,nsp),
     .    ggr(1,1),ggr(1,nsp),agr(1,2*nsp-1),agr(1,4),gagr(1,2*nsp-1),
     .    gagr(1,1),gagr(1,nsp),vnl,vnl(1,1,1,nsp),enl)
        else
          call vxnloc(np,nsp,rho,rho(1,nsp),agr(1,1),agr(1,nsp),
     .    ggr(1,1),ggr(1,nsp),agr(1,2*nsp-1),agr(1,4),gagr(1,2*nsp-1),
     .    gagr(1,1),gagr(1,nsp),vnl,vnl(1,1,1,nsp),enl)
        endif
      endif

C      call prm3('enl',i,enl,n1,n2,n3)
c      do  i = 1, nsp
c        call prm3('vnl(isp=%i)',i,vnl(1,1,1,i),n1,n2,n3)
c      enddo


C --- Make nonlocal reps, rmu ---
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

ccccccccccccccccccccccccccccccc
        if(plottest) then
        isp=1
c        smmin=0d0
        do i1=1,n1
        do i2=1,n2
        do i3=1,n3
          write(8007,"(3i4,10e14.6)") i1,i2,i3,vxcnl(i1,i2,i3,isp)
          write(9007,"(3i4,10e14.6)") i1,i2,i3,enl(i1,i2,i3) ,enlbk(i1,i2,i3) ,enl(i1,i2,i3)-enlbk(i1,i2,i3) 
c          sss=abs(enl(i1,i2,i3)-enlbk(i1,i2,i3) )
c          if(sss>smmin)smmin=sss
        enddo
          write(8007,*)
          write(9007,*)
        enddo
        enddo
c        print *,' max error=',smmin
c        stop ' test end of plot 8007'
        endif
ccccccccccccccccccccccccccccccc

      deallocate(rho,agr,gagr,ggr,enl,vnl)
      call tcx('vxcnlm')
      end

