      subroutine vxcnlm(lxcg,nsp,k1,k2,k3,slat,smrho,
     .repnl,rmunl,vavgnl,vxnl,vcnl,vxcnl)
      use m_struc_def, only:s_lat
C- Gradient correction to smoothed rho(q) tabulated on a mesh
C ----------------------------------------------------------------------
Ci Inputs
Ci   lxcg  : Specifies GGA
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
Co Outputs
Co   repnl :integral smrho * epsnl[smrho]
Co   rmunl :integral smrho * vxcnl[smrho]
Ci   vavgnl:average NL XC potential
Ci   vcnl  :correlation part of vxcnl
Ci   vxnl  :exchange part of vxcnl
Co   vxcnl :nonlocal XC potential on uniform mesh added to vxcnl
Cl Local variables
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
C     implicit none
C ... Passed parameters
      integer lxcg,k1,k2,k3,nsp
      double precision repnl(2),rmunl(2),vavgnl(2) !slat(1),
      double complex   smrho(k1,k2,k3,nsp),vxcnl(k1,k2,k3,nsp)
      double complex   vxnl(k1,k2,k3,nsp),vcnl(k1,k2,k3,nsp)
C ... Local parameters
      integer ip,i,i1,i2,i3,lcut,ng,np,ipr,n1,n2,n3,ngabc(3)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      integer ogv,okv
      double precision alat,vol,xx
      real(8),allocatable :: ggr(:,:),agr(:,:),gagr(:,:),rho(:,:)
      real(8),allocatable :: enl(:,:,:),vnl(:,:,:,:)
      complex(8),allocatable:: zgrho(:,:,:),gzgrho(:,:,:),zggrho(:,:)
c
      type(s_lat)::  slat


C     Heap
      integer w(1)
      common /w/ w

C --- Setup and memory allocation ---
      call tcn('vxcnlm')
      call getpr(ipr)
c      call upack('lat nabc ng ogv okv',slat,ngabc,ng,ogv,okv,0)
c      call upack('lat alat vol',slat,alat,vol,0,0,0)
      ngabc =slat%nabc
      ng    =slat%ng
      ogv   =slat%ogv
      okv   =slat%okv
      alat  =slat%alat
      vol   =slat%vol

      np = k1*k2*k3
      allocate(agr(np,3*nsp-2),gagr(np,2*nsp-1),ggr(np,nsp))
      allocate(zgrho(np,3,nsp),gzgrho(np,3,2*nsp-1),zggrho(np,nsp))

C      call defrr(oagr,   np*(3*nsp-2))
C      call defrr(ogagr,  np*3*nsp)
C      call defrr(oggr,   np*nsp)
C      call defrr(oexc,   np*nsp)

C --- Grad rho_i and Laplacian rho_i (complex) ---
      do  i = 1, nsp
C       call zprm3('smrho(isp=%i)',i,smrho(1,1,1,i),n1,n2,n3)
        call grfmsh(601,alat,ng,w(ogv),w(okv),k1,k2,k3,n1,n2,n3,
     .  smrho(1,1,1,i),zgrho(1,1,i),zggrho(1,i))
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
        call grfmsh(201,alat,ng,w(ogv),w(okv),k1,k2,k3,n1,n2,n3,
     .  zggrho,gzgrho(1,1,i),xx)

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
C      do  i = 1, nsp
C        call prm3('vnl(isp=%i)',i,vnl(1,1,1,i),n1,n2,n3)
C      enddo

C --- Make nonlocal reps, rmu ---
      do  i = 1, nsp
        repnl(i) = 0
        rmunl(i) = 0
        vavgnl(i) = 0
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
      deallocate(rho,agr,gagr,ggr,enl,vnl)

      call tcx('vxcnlm')
      end

