      subroutine gwinput_v2x(ifi,
     &  nsp,   nbas,   nclass, nrmx,   ncoremx,
     &  lmxamx,ngpmx,  nband,  ldim2,   nqtt ,
     &  nphi,nphimx,nindx,lindx,ibasindx,
     &  iclass,lmxa,nr,konf,ncore,
     &  zz,aa,bb,bas,alat,plat,
     &  qtt,efermi,ec,evl,vxclda,gx,gcore,cphi,geig,nocore)
c
c     &  ifomit ) !If ifomit=T, omit to store cphi and geig
c- Read and Write the informations required for GW calculations.
c
cr When you write these data, you first have to open a file DATA4GW with the file handle ifi>0.
cr ,e.g,  open(ifi ,file='DATA4GW',form='unformatted')
cr
cr Radial mesh is specified from, nr, aa, and bb; mesh is r(i) = bb*(exp(aa*(i-1)) -1) ,i=1,nr
cr
cr When you read these data, you first need to read the sizes by the subbroutine gwinput0
cr in order to allocate required arrays.
cr
cr NLAindx file should be written by echo 0|lmfgw -------
cr   It should be
cr----NLAindx start---------------
cr     ldim2
cr     n l ibas for each line --- Repeat this line for the number of ldim2.
cr                            --- nphimx=maximum value of n.
cr
c----------------------------------------------------------------------------------------------
C     implicit none
      integer(4) :: ifi, !File handle. Write for ifi<0. Read for ifi>0.
     &  nsp,      ! =1 or 2, corresponding to para or ferro.
     &  nbas,     ! Number of atom in the primitive cell
     &  nclass,   ! Number of atomic class (or type) for in the primitive cell
     &  nrmx,     ! = maxval(nr(1:nclass))  Maximum number of nr
     &  ncoremx,  ! = maxval(ncore(1:nclass))
c
     &  lmxamx,   ! = maxval(lmxa(1:nclass))
     &  ngpmx,    ! Maximum number of G vector.
     &  nband,    ! Number of bands given by GWIN0
cv2     &  ldim,     ! = sum ( (lmxa(1:nbas)+1)**2 )
     &  nqtt,    ! = nqbz*(1+nq0i). Number of q points given by qg4gw
c----
     &  iclass(nbas),   ! class is starting from 1.
     &  lmxa  (nclass), ! Maximum l number for each atom for augmentation.
     &  nr(nclass),     ! Size of radial mesh.
     &  konf(0:lmxamx,nclass),! Principle quantum numbers of valence electron.
      ! For example, 4s4p3d means konf(0)=4, konf(1)=4, konf(2)=3.
      ! Core orbitals are specified by
      !   1, 2,... konf(0)-1 for s
      !   2, 3,... konf(1)-1 for p
      !   3, 4,... konf(2)-1 for d, and so on.
      !
     &  ncore(nclass)   ! ncore = \sum_l  (konf(l)-1) - l
      ! Number of different core orbitals for each atom.

cv2
      integer(4):: nphi,          ! number of augmentation nRl channels, as distinct from:
     &             ldim2,         ! number of nRLm channels = 2*ldim for phi+phidot case.
     &             nphimx,        ! Maxmum number of phi for all l ---  2 for phi+phidot case.
     &             nindx(ldim2),  ! n    index
     &             lindx(ldim2),  ! l    index
     &             ibasindx(ldim2)! ibas index
c
      real(8) ::
     &  zz(nclass),     ! Atomic number.
     &  aa(nclass),bb(nclass),! Radial mesh are specified by these parameters with nr.
     &  bas(3,nbas),    ! Atomic posion in the Cartesian coordinate (alat unit),.
     &  alat,           ! Lattice constant in a.u.
     &  plat(3,3),      ! Primitive translation vector in the unit of alat. plat(1:3,1) is 1st vector, ans so on.
     &  qtt(3,nqtt),  ! q points given by qg4gw
     &  efermi,         ! Fermi energy. It should be calculated for n1 n2 n3 given by GIWN0.
c----
     &  ec(ncoremx, nclass, nsp),   ! Eigenvalues for core
     &  evl    (nband, nqtt, nsp), ! Eigenvalues
     &  vxclda (nband, nqtt, nsp), ! Lda XC potential <psi|Vxc(n_total)  |psi>
      !  or   <psi|Vxc(n_valence)|psi>. See nocore.

cv2     &  gx (nrmx, 0:lmxamx, 2,  nclass,nsp), ! Radial function. gx (nrmx, lmxamx, 1, nclass,nsp) is phi
cv2                                          !                  gx (nrmx, lmxaxm, 2, nclass,nsp) is phidot
     &  gx (nrmx, 0:lmxamx, nphimx, nclass,nsp), ! Radial function. gx (nrmx, lmxamx, 1, nclass,nsp) is phi
      !                  gx (nrmx, lmxaxm, 2, nclass,nsp) is phidot
     &  gcore(nrmx, ncoremx, nclass,nsp)  ! Core radial function.
c
c     These radial functions are only the major part given by the scalar relativistic calculations.
c     gx and gcore = r \phi(r) = u(r), where \phi is the major part of the true radial functions.
c     gx is normalized as 1 = \int dr gx**2
c     gcore is the major part of the true radial functions.
c

      complex(8)::
cv2     &  cphi(ldim, 2, nband, nqtt,nsp), ! Coefficients of eigenfunction
cv2                                         ! for phi    = cphi(ldim, 1, nband, nqtt), and
cv2                                         ! for phidot = cphi(ldim, 2, nband, nqtt).
     &  cphi(ldim2,  nband, nqtt,nsp),     ! Coefficients of eigenfunction
c
     &  geig(ngpmx,  nband, nqtt,nsp)  ! Coefficienets of eigenfunction for IPW.


      logical ::
     &   nocore  ! Flag file NoCore exist or not.
      ! nocore is obtained by inquire(file='NoCore',exist=nocore) in the upper routine.
      ! If nocore exist. you have to supply
      !  <psi|Vxc(n_valence)|psi>  to  vxclda (nband, nqtt).
      ! If not, <psi|Vxc(n_total)|psi>  to  vxclda.
c      logical :: ifomit
c----
      if(ifi>0) then
c        if(ifomit) then
        read(ifi)
     &              iclass,lmxa,nr,konf,ncore,
     &              zz,aa,bb,bas,alat,plat,
     &              qtt,efermi,ec,evl,vxclda,gx,gcore,nocore !omit cphi and geig
cv2
     &              ,nindx,lindx,ibasindx
c        else
c          read(ifi)
c     &              iclass,lmxa,nr,konf,ncore,
c     &              zz,aa,bb,bas,alat,plat,
c     &              qtt,efermi,ec,evl,vxclda,gx,gcore,cphi,geig,nocore
ccv2
c     &              ,nindx,lindx,ibasindx
c        endif
ccccccccccccccccc
        write(6,*)'gwinput_v2x: DATA4GW_V2 readin size=',
     &  nsp,nbas, nclass, nrmx,ncoremx
        write(6,*)'                 =', lmxamx,ngpmx,nband,ldim2, nqtt
        write(6,*)'ncoremx=',ncoremx
        write(6,*)'ecup=',ec(1:ncoremx, 1, 1)
        if(nsp==2) write(6,*)'ecdn=',ec(1:ncoremx, 1, 2)
ccccccccccccccccc

ccccccccccccccccccccccc
c--- This is a swap test ---
c      call swapc(cphi(5:9,    1:nband,1:nqtt,1:nsp)
c     &          ,cphi(201:205,1:nband,1:nqtt,1:nsp),5*nband*nqtt*nsp)
c      call swapr(gx(1:nrmx,2,1,1,1:nsp),gx(1:nrmx,2,3,1,1:nsp),nrmx*nsp)
c      print *,' ----swap test------'
ccccccccccccccccccccccc


ccccccccccccccccccccccc
c       cphi(201:205,1:nband,1:nqtt,1:nsp) = 0d0
c       print *,' --- skip 4d test ---'
ccccccccccccccccccccccc


c---------------------------
      elseif(ifi<0) then
        write(-ifi)
     &              nsp,   nbas,   nclass, nrmx,   ncoremx,
     &              lmxamx,ngpmx,  nband,  ldim2,   nqtt,
cv2
     &              nphi,nphimx
c
        write(6,*)
        write(6,*)' *** DATA4GW_V2 size=', nsp,nbas,nclass,nrmx,ncoremx
        write(6,*)'                 =', lmxamx,ngpmx, nband,ldim2, nqtt
        write(6,*)'                 =', nphi, nphimx

c        if(ifomit) then
        write(-ifi)
     &              iclass,lmxa,nr,konf,ncore,
     &              zz,aa,bb,bas,alat,plat,
     &              qtt,efermi,ec,evl,vxclda,gx,gcore,nocore !omit cphi and geig
cv2
     &              ,nindx,lindx,ibasindx
c        else
c          write(-ifi)
c     &              iclass,lmxa,nr,konf,ncore,
c     &              zz,aa,bb,bas,alat,plat,
c     &              qtt,efermi,ec,evl,vxclda,gx,gcore,cphi,geig,nocore
ccv2
c     &              ,nindx,lindx,ibasindx
c        endif
      else
        stop 'gwinput_v2x: wrong ifi'
      endif
      end
c----------------------------------------------------------------------
      subroutine gwinput0_v2(ifi,
     &              nsp,   nbas,   nclass, nrmx,   ncoremx,
     &              lmxamx,ngpmx,  nband,  ldim2,   nqtt,
     &              nphi, nphimx)
c  Read dimensions only.
C     implicit none
      integer(4)::  ifi,
     &              nsp,   nbas,   nclass, nrmx,   ncoremx,
     &              lmxamx,ngpmx,  nband,  ldim2,   nqtt ,
     &              nphi, nphimx
      read(ifi)
     &              nsp,   nbas,   nclass, nrmx,   ncoremx,
     &              lmxamx,ngpmx,  nband,  ldim2,   nqtt,
     &              nphi, nphimx
      end
c--------------------
      subroutine rdibasindx(ldim2,nphi,nphimx,nindx,lindx,ibasindx,ifi)
C- Reads indexing data for augmentation functions
C ----------------------------------------------------------------------
Ci Inputs
Ci   ldim2 :number of augmentation functions nRlm in unit cell
Ci         :(used here as a sanity check)
Ci   ifi   :file handle ifi>0 for file read
Co Outputs
Co   nphi  :number of augmentation functions nRl  in unit cell
Co   nphimx:maximum number of aug. functions in a single Rl channel
Co   nindx :nindx(i) type of augm. function for the ith Rl channel
Co         :e.g. nindx=1 : phi, nindx=2 : phidot,  nindx=3 : loc. orb
Co   lindx :l quantum number for the ith Rl channel
Co ibasindx:site index for the ith Rl channel
Cr Remarks
Cu Updates
Cu   05 Mar 02 First created
C ----------------------------------------------------------------------
C     implicit none
      integer ldim2,nphi,nphimx,ifi,
     .  nindx(ldim2),           ! n    index
     .  lindx(ldim2),           ! l    index
     .  ibasindx(ldim2)         ! ibas index

      character*80 outs
      integer i,ipqn,l,ib,ii

C --- File read ---
      if (ifi .gt. 0) then
        rewind ifi
        read(ifi,'(a)') outs
        read(ifi,*) i
C       If passed ldim2>0, check that it matches file value
        if (ldim2 .gt. 0 .and. i .ne. ldim2)
     .    call rx('rdibasindx: file mismatch ldim2')
        ldim2 = i
        nphimx = 0
        do  i = 1, ldim2
          read(ifi,'(a)',err=101,end=101) outs
          read(outs,*) ipqn,l,ib,ii
          nindx(i) = ipqn
          lindx(i) = l
          ibasindx(i) = ib
          nphimx = max(nphimx,ipqn)
          nphi = i
        enddo
  101   continue
C --- File read ---
      else
        call rx('ioaindx: file write not implemented')
      endif

      end

c--------------------
      subroutine rdibasindx2(ldim2,nphi,nphimx,nindx,lindx,ibasindx,ifi)
C- Reads indexing data for augmentation functions
C ----------------------------------------------------------------------
Ci Inputs
Ci   ldim2 :number of augmentation functions nRlm in unit cell
Ci         :(used here as a sanity check)
Ci   ifi   :file handle ifi>0 for file read
Co Outputs
Co   nphi  :number of augmentation functions nRl  in unit cell
Co   nphimx:maximum number of aug. functions in a single Rl channel
Co   nindx :nindx(i) type of augm. function for the ith Rl channel
Co         :e.g. nindx=1 : phi, nindx=2 : phidot,  nindx=3 : loc. orb
Co   lindx :l quantum number for the ith Rl channel
Co ibasindx:site index for the ith Rl channel
Cr Remarks
Cu Updates
Cu   05 Mar 02 First created
C ----------------------------------------------------------------------
C     implicit none
      integer ldim2,nphi,nphimx,ifi,
     .  nindx(ldim2),           ! n    index
     .  lindx(ldim2),           ! l    index
     .  ibasindx(ldim2)         ! ibas index

      character*80 outs
      integer i,ipqn,l,ib,ii,nn,m

C --- File read ---
      if (ifi .gt. 0) then
        rewind ifi
        read(ifi,'(a)') outs
        read(ifi,*) i
C       If passed ldim2>0, check that it matches file value
        if (ldim2 .gt. 0 .and. i .ne. ldim2)
     .    call rx('rdibasindx: file mismatch ldim2')
        ldim2 = i
        nphimx = 0
c        do  i = 1, ldim2
        nn = 0
        do
          read(ifi,'(a)',err=101,end=101) outs
          read(outs,*) ipqn,l,ib,ii
          do m=1,2*l+1
            nn = nn+1
            nindx   (nn) = ipqn
            lindx   (nn) = l
            ibasindx(nn) = ib
          enddo
          nphimx = max(nphimx,ipqn)
          nphi = i
        enddo
  101   continue

C --- File read ---
      else
        call rx('ioaindx: file write not implemented')
      endif

      end

c--------------------
      subroutine swapc(a,b,n)
      complex(8),allocatable:: c(:)
      complex(8):: a(n),b(n)
      allocate(c(n))
      c=a
      a=b
      b=c
      deallocate(c)
      end
c--------------------
      subroutine swapr(a,b,n)
      real(8),allocatable:: c(:)
      real(8):: a(n),b(n)
      allocate(c(n))
      c=a
      a=b
      b=c
      deallocate(c)
      end















