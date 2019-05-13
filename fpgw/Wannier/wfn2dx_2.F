      subroutine wfn2dx_2(alat,plat,nsp,nq_wfn,nband_wfn,q_wfn,bindx_wfn,
     &     mesh,rini,rfin,phipw,phiaug,phitot)
      implicit none
c input
      double precision :: alat,plat(3,3),rini(3),rfin(3)
      integer :: nsp,nq_wfn,nband_wfn,bindx_wfn(nband_wfn),mesh(3)
      double precision :: q_wfn(3,nq_wfn)
      double complex :: 
     &     phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp),
     &     phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

c local
      integer :: isp,iq,ib,iband,ifile,i1,i2,i3,ifile_handle
      double precision :: rtmp(3),r(3)
      character*(23) :: fn

      ifile=ifile_handle()
      do isp=1,nsp
      do iq=1,nq_wfn
      do ib=1,nband_wfn
        iband=bindx_wfn(ib)
        write(fn,"(a4,i1,a1,i4.4,a1,i4.4,a8)")
     &       'phis',isp,'q',iq,'b',iband,'.general'

        write(*,*) 'open ',fn
        open(ifile,file=fn)

        write(ifile,"(a,i2)") '#isp=',isp
        write(ifile,"(a,i5,3f10.4)") '#iq_wfn,q=',iq,q_wfn(1:3,iq)
        write(ifile,"(a,2i5)") '#ib,bindx=',ib,iband
        write(ifile,"(a)") 'header = marker "DATA\n"'
        write(ifile,"(a,i4,a,i4,a,i4)")
     &       'grid = ', mesh(1)+1,' x ',mesh(2)+1,' x ',mesh(3)+1
        write(ifile,"(a)") 'format = ascii'      
        write(ifile,"(a)") 'interleaving = field'
        write(ifile,"(a)") 'majority = row'
        write(ifile,"(a,a,a)") 'field = locations, ',
     &       'Re_phipw, Im_phipw, Re_phiaug, ',
     &       'Im_phiaug, Re_phitot, Im_phitot, |phitot|'
        write(ifile,"(a,a)") 'structure = 3-vector,',
     &       ' scalar, scalar, scalar, scalar, scalar, scalar, scalar'
        write(ifile,"(a,a)") 'type = float,',
     &       'float, float, float, float, float, float, float'
        write(ifile,"(a,a,a)") 'dependency = positions,',
     &       ' positions, positions,',
     &       ' positions, positions, positions, positions, positions'
        write(ifile,"(a)") 'end'
        write(ifile,"(a)") 'DATA'
        do i1=1,mesh(1)+1
          do i2=1,mesh(2)+1
            do i3=1,mesh(3)+1
              rtmp(1)=(i1-1)/dble(mesh(1))
              rtmp(2)=(i2-1)/dble(mesh(2))
              rtmp(3)=(i3-1)/dble(mesh(3))
!              call mymatvec(plat,rtmp,r,3,3)
              r(:) = rini(:) + (rfin(:)-rini(:))*rtmp(:)
              r(1:3)=alat*(r(1:3))

              write(ifile,"(10f20.14)") r(1:3),
     &             phipw(i1,i2,i3,ib,iq,isp),
     &             phiaug(i1,i2,i3,ib,iq,isp),
     &             phitot(i1,i2,i3,ib,iq,isp),
     &             abs(phitot(i1,i2,i3,ib,iq,isp))
            enddo
          enddo
        enddo
        close(ifile)
      enddo
      enddo
      enddo
      end subroutine wfn2dx_2
ccccccccccccccccccccccccccccccccccc
      subroutine crystal2dx_2(alat,plat,rini,rfin,
     &   natom,apos,nclass,iclass,zz)
      implicit none
      double precision :: alat,plat(3,3),rini(3),rfin(3)
      integer :: natom,nclass,iclass(natom)
      double precision :: apos(3,natom),zz(nclass)

c local
      integer :: nrange(3)
      double precision :: tmpvec(3),Z(natom)
      integer :: i,j,ic,ia

      do ia=1,natom
        ic=iclass(ia)
        Z(ia)=zz(ic)
      enddo
c 
      nrange(1:3)=2

      call dump_ucell(alat,plat,'ucell.general')
      call dump_box_2(alat,rini,rfin,plat,'box.general')
      call dump_apos_2(alat,plat,rini,rfin,natom,apos,
     &   Z,nrange,'apos.general')
      call dump_apos_xyz(alat,plat,rini,rfin,natom,apos,
     &   Z,nrange,'apos.xyz')
      end subroutine crystal2dx_2
cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dump_ucell(alat,plat,fn)
      implicit none
      double precision,intent(in) :: alat,plat(3,3)
      character*(*) :: fn
c local
      integer :: ifile,i1,i2,i3,ifile_handle
      double precision :: v1(3),v2(3)
      ifile=ifile_handle()
      open(ifile,file=fn)
      write(ifile,"(a)") 'header = marker "DATA\n"'
      write(ifile,"(a)") 'grid = 2 x 2 x 2'
      write(ifile,"(a)") 'format = ascii'
      write(ifile,"(a)") 'interleaving = field'
      write(ifile,"(a)") 'majority = row'
      write(ifile,"(a)") 'field = locations, dummy'
      write(ifile,"(a)") 'structure = 3-vector, scalar'
      write(ifile,"(a)") 'type = float, float'
      write(ifile,"(a)") 'dependency = positions, positions'
      write(ifile,"(a)") 'end'
      write(ifile,"(a)") 'DATA'

      do i1=0,1
      do i2=0,1
      do i3=0,1
        v1(1)=dble(i1)
        v1(2)=dble(i2)
        v1(3)=dble(i3)
        call mymatvec(plat,v1,v2,3,3)
        v2(1:3)=alat*v2(1:3)
        write(ifile,"(3f20.12,f8.2)") v2(1:3),0.0d0
      enddo
      enddo
      enddo
      close(ifile)
      end subroutine dump_ucell
cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dump_apos_2(alat,plat,rini,rfin,natom,apos,Z,nrange,fn)
      implicit none
      double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3),
     &     apos(3,natom),Z(natom)
      integer,intent(in) :: natom,nrange(3)
      character*(*) :: fn
c local
      integer :: natomall,ifile_handle
      integer :: ifile,i,i1,i2,i3
      double precision :: v1(3),v2(3)

      natomall=natom*(2*nrange(1)+1)*(2*nrange(2)+1)*(2*nrange(3)+1)
      ifile=ifile_handle()
      open(ifile,file=fn)
      write(ifile,"(a)") 'header = marker "DATA\n"'
      write(ifile,"(a,i6)") 'points = ', natomall
      write(ifile,"(a)") 'format = ascii'
      write(ifile,"(a)") 'interleaving = field'
      write(ifile,"(a)") 'majority = row'
      write(ifile,"(a)") 'field = locations, Z'
      write(ifile,"(a)") 'structure = 3-vector, scalar'
      write(ifile,"(a)") 'type = float, float'
      write(ifile,"(a)") 'dependency = positions, positions'
      write(ifile,"(a)") 'end'
      write(ifile,"(a)") 'DATA'

      do i=1,natom
        do i1=-nrange(1),nrange(1)
        do i2=-nrange(2),nrange(2)
        do i3=-nrange(3),nrange(3)
          v1(1)=dble(i1)
          v1(2)=dble(i2)
          v1(3)=dble(i3)
          call mymatvec(plat,v1,v2,3,3)
          v2(1:3)=alat*(v2(1:3)+apos(1:3,i))
          write(ifile,"(3f20.12,f8.2)") v2(1:3),Z(i)
        enddo
        enddo
        enddo
      enddo
      close(ifile)
      end subroutine dump_apos_2
cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dump_apos_xyz(alat,plat,rini,rfin,natom,apos,
     & Z,nrange,fn)
      implicit none
      double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3),
     &     apos(3,natom),Z(natom)
      integer,intent(in) :: natom,nrange(3)
      character*(*) :: fn
c local
      integer:: ifile_handle
      integer :: natomall,natomx
      integer :: ifile,i,i1,i2,i3
      double precision :: v1(3),v2(3),aini(3),afin(3),eps
      double precision,allocatable :: rall(:,:)
      character*2,allocatable :: zall(:)
      character*2 :: z2element

      eps = 0.05d0
      aini = alat*(rini-eps)
      afin = alat*(rfin+eps)

      natomx=natom*(2*nrange(1)+1)*(2*nrange(2)+1)*(2*nrange(3)+1)
      allocate(rall(3,natomx),zall(natomx))
      natomall = 0
      do i=1,natom
        do i1=-nrange(1),nrange(1)
        do i2=-nrange(2),nrange(2)
        do i3=-nrange(3),nrange(3)
          v1(1)=dble(i1)
          v1(2)=dble(i2)
          v1(3)=dble(i3)
          call mymatvec(plat,v1,v2,3,3)
          v2(1:3)=alat*(v2(1:3)+apos(1:3,i))
          if ( (v2(1).ge.aini(1).and.v2(1).le.afin(1))
     &    .and.(v2(2).ge.aini(2).and.v2(2).le.afin(2))
     &    .and.(v2(3).ge.aini(3).and.v2(3).le.afin(3)) ) then
             natomall = natomall + 1
             rall(1:3,natomall) = v2(1:3)
             zall(natomall) = z2element(Z(i))
          endif   
        enddo
        enddo
        enddo
      enddo

      ifile=ifile_handle()
      open(ifile,file=fn)
      write(ifile,*)natomall
      write(ifile,*)
      do i=1,natomall
         write(ifile,"(a2,3f20.12)") zall(i),rall(1:3,i)
      enddo
      close(ifile)

      deallocate(rall,zall)
      end subroutine dump_apos_xyz
cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dump_box_2(alat,rini,rfin,plat,fn)
      implicit none
      integer:: ifile_handle
      double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3)
      double precision :: eps
      character*(*) :: fn
c local
      integer :: ifile,i1,i2,i3
      double precision :: v1(3),v2(3)
      ifile=ifile_handle()
      open(ifile,file=fn)
      write(ifile,"(a)") 'header = marker "DATA\n"'
      write(ifile,"(a)") 'grid = 2 x 2 x 2'
      write(ifile,"(a)") 'format = ascii'
      write(ifile,"(a)") 'interleaving = field'
      write(ifile,"(a)") 'majority = row'
      write(ifile,"(a)") 'field = locations, dummy'
      write(ifile,"(a)") 'structure = 3-vector, scalar'
      write(ifile,"(a)") 'type = float, float'
      write(ifile,"(a)") 'dependency = positions, positions'
      write(ifile,"(a)") 'end'
      write(ifile,"(a)") 'DATA'

      eps = 0.1d0
      do i1=0,1
        if (i1.eq.0) then
          v1(1)= -eps
        else
          v1(1)= 1.0d0 + eps
        endif   
      do i2=0,1
        if (i2.eq.0) then
          v1(2)= -eps
        else
          v1(2)= 1.0d0 + eps
        endif   
      do i3=0,1
        if (i3.eq.0) then
          v1(3)= -eps
        else
          v1(3)= 1.0d0 + eps
        endif   
!        call mymatvec(plat,v1,v2,3,3)
        v2(:)=rini(:)+(rfin(:)-rini(:))*v1(:)
        v2(1:3)=alat*v2(1:3)
        write(ifile,"(3f20.12,f8.2)") v2(1:3),0.0d0
      enddo
      enddo
      enddo
      close(ifile)
      end subroutine dump_box_2
cccccccccccccccccccccccccccccccccccccccccccccccc
      character*2 function z2element(Z)
      implicit none
c input
      double precision :: Z
      integer :: iz

      iz=nint(Z)

      if (iz.eq.1) then
         z2element=' H'
      elseif (iz.eq.2) then
         z2element='He'
      elseif (iz.eq.3) then
         z2element='Li'
      elseif (iz.eq.4) then
         z2element='Be'
      elseif (iz.eq.5) then
         z2element=' B'
      elseif (iz.eq.6) then
         z2element=' C'
      elseif (iz.eq.7) then
         z2element=' N'
      elseif (iz.eq.8) then
         z2element=' O'
      elseif (iz.eq.9) then
         z2element=' F'
      elseif (iz.eq.10) then
         z2element='Ne'
      elseif (iz.eq.11) then
         z2element='Na'
      elseif (iz.eq.12) then
         z2element='Mg'
      elseif (iz.eq.13) then
         z2element='Al'
      elseif (iz.eq.14) then
         z2element='Si'
      elseif (iz.eq.15) then
         z2element=' P'
      elseif (iz.eq.16) then
         z2element=' S'
      elseif (iz.eq.17) then
         z2element='Cl'
      elseif (iz.eq.18) then
         z2element='Al'
      elseif (iz.eq.19) then
         z2element=' K'
      elseif (iz.eq.20) then
         z2element='Ca'
      elseif (iz.eq.21) then
         z2element='Sc'
      elseif (iz.eq.22) then
         z2element='Ti'
      elseif (iz.eq.23) then
         z2element=' V'
      elseif (iz.eq.24) then
         z2element='Cr'
      elseif (iz.eq.25) then
         z2element='Mn'
      elseif (iz.eq.26) then
         z2element='Fe'
      elseif (iz.eq.27) then
         z2element='Co'
      elseif (iz.eq.28) then
         z2element='Ni'
      elseif (iz.eq.29) then
         z2element='Cu'
      elseif (iz.eq.30) then
         z2element='Zn'
      elseif (iz.eq.31) then
         z2element='Ga'
      elseif (iz.eq.32) then
         z2element='Ge'
      elseif (iz.eq.33) then
         z2element='As'
      elseif (iz.eq.34) then
         z2element='Se'
      elseif (iz.eq.35) then
         z2element='Br'
      elseif (iz.eq.36) then
         z2element='Kr'
      elseif (iz.eq.37) then
         z2element='Rb'
      elseif (iz.eq.38) then
         z2element='Sr'
      elseif (iz.eq.39) then
         z2element=' Y'
      elseif (iz.eq.40) then
         z2element='Zr'
      elseif (iz.eq.41) then
         z2element='Nb'
      elseif (iz.eq.42) then
         z2element='Mo'
      elseif (iz.eq.43) then
         z2element='Tc'
      elseif (iz.eq.44) then
         z2element='Ru'
      elseif (iz.eq.45) then
         z2element='Rh'
      elseif (iz.eq.46) then
         z2element='Pd'
      elseif (iz.eq.47) then
         z2element='Ag'
      elseif (iz.eq.48) then
         z2element='Cd'
      elseif (iz.eq.49) then
         z2element='In'
      elseif (iz.eq.50) then
         z2element='Sn'
      elseif (iz.eq.51) then
         z2element='Sb'
      elseif (iz.eq.52) then
         z2element='Te'
      elseif (iz.eq.53) then
         z2element=' I'
      elseif (iz.eq.54) then
         z2element='Xe'
      elseif (iz.eq.55) then
         z2element='Cs'
      elseif (iz.eq.56) then
         z2element='Ba'
      elseif (iz.eq.57) then
         z2element='La'
      elseif (iz.eq.58) then
         z2element='Ce'
      elseif (iz.eq.59) then
         z2element='Pr'
      elseif (iz.eq.60) then
         z2element='Nd'
      elseif (iz.eq.61) then
         z2element='Pm'
      elseif (iz.eq.62) then
         z2element='Sm'
      elseif (iz.eq.63) then
         z2element='Eu'
      elseif (iz.eq.64) then
         z2element='Gd'
      elseif (iz.eq.65) then
         z2element='Tb'
      elseif (iz.eq.66) then
         z2element='Dy'
      elseif (iz.eq.67) then
         z2element='Ho'
      elseif (iz.eq.68) then
         z2element='Er'
      elseif (iz.eq.69) then
         z2element='Tm'
      elseif (iz.eq.70) then
         z2element='Yb'
      elseif (iz.eq.71) then
         z2element='Lu'
      elseif (iz.eq.72) then
         z2element='Hf'
      elseif (iz.eq.73) then
         z2element='Ta'
      elseif (iz.eq.74) then
         z2element=' W'
      elseif (iz.eq.75) then
         z2element='Re'
      elseif (iz.eq.76) then
         z2element='Os'
      elseif (iz.eq.77) then
         z2element='Ir'
      elseif (iz.eq.78) then
         z2element='Pt'
      elseif (iz.eq.79) then
         z2element='Au'
      elseif (iz.eq.80) then
         z2element='Hg'
      elseif (iz.eq.81) then
         z2element='Tl'
      elseif (iz.eq.82) then
         z2element='Pb'
      elseif (iz.eq.83) then
         z2element='Bi'
      elseif (iz.eq.84) then
         z2element='Po'
      elseif (iz.eq.85) then
         z2element='At'
      elseif (iz.eq.86) then
         z2element='Rn'
      elseif (iz.eq.87) then
         z2element='Fr'
      elseif (iz.eq.88) then
         z2element='Ra'
      elseif (iz.eq.89) then
         z2element='Ac'
      elseif (iz.eq.90) then
         z2element='Th'
      elseif (iz.eq.91) then
         z2element='Pa'
      elseif (iz.eq.92) then
         z2element=' U'
      elseif (iz.eq.93) then
         z2element='Np'
      elseif (iz.eq.94) then
         z2element='Pu'
      elseif (iz.eq.95) then
         z2element='Am'
      elseif (iz.eq.96) then
         z2element='Cm'
      elseif (iz.eq.97) then
         z2element='Bk'
      elseif (iz.eq.98) then
         z2element='Cf'
      elseif (iz.eq.99) then
         z2element='Es'
      elseif (iz.eq.100) then
         z2element='Fm'
      elseif (iz.eq.101) then
         z2element='Md'
      elseif (iz.eq.102) then
         z2element='No'
      elseif (iz.eq.103) then
         z2element='Lr'
      endif

c debug:
! for TiFe2
!      if (iz.eq.26) z2element=' B' ! Fe
!      if (iz.eq.22) z2element=' C' ! Ti

! for TiFe2
      if (iz.eq.57) z2element=' N' ! La
      if (iz.eq.38) z2element=' H' ! Sr
      if (iz.eq.22) z2element=' C' ! Ti
      if (iz.eq.13) z2element=' O' ! Al
      if (iz.eq.8)  z2element=' B' ! O

      end function z2element

