subroutine wfn2dx_2(alat,plat,nsp,nq_wfn,nband_wfn,q_wfn,bindx_wfn, &
     mesh,rini,rfin,phipw,phiaug,phitot)
  implicit none
  ! input
  double precision ::  alat,plat(3,3),rini(3),rfin(3)
  integer :: nsp,nq_wfn,nband_wfn,bindx_wfn(nband_wfn),mesh(3)
  double precision ::  q_wfn(3,nq_wfn)
  double complex :: &
       phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
       phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp), &
       phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nq_wfn,nsp)

  ! local
  integer :: isp,iq,ib,iband,ifile,i1,i2,i3,ifile_handle
  double precision ::  rtmp(3),r(3)
  character*(23) :: fn

  ifile=ifile_handle()
  do isp=1,nsp
     do iq=1,nq_wfn
        do ib=1,nband_wfn
           iband=bindx_wfn(ib)
           write(fn,"(a4,i1,a1,i4.4,a1,i4.4,a8)") &
                'phis',isp,'q',iq,'b',iband,'.general'

           write(*,*) 'open ',fn
           open(ifile,file=fn)

           write(ifile,"(a,i2)") '#isp=',isp
           write(ifile,"(a,i5,3f10.4)") '#iq_wfn,q=',iq,q_wfn(1:3,iq)
           write(ifile,"(a,2i5)") '#ib,bindx=',ib,iband
           write(ifile,"(a)") 'header = marker "DATA\n"'
           write(ifile,"(a,i4,a,i4,a,i4)") &
                'grid = ', mesh(1)+1,' x ',mesh(2)+1,' x ',mesh(3)+1
           write(ifile,"(a)") 'format = ascii'
           write(ifile,"(a)") 'interleaving = field'
           write(ifile,"(a)") 'majority = row'
           write(ifile,"(a,a,a)") 'field = locations, ', &
                'Re_phipw, Im_phipw, Re_phiaug, ', &
                'Im_phiaug, Re_phitot, Im_phitot, |phitot|'
           write(ifile,"(a,a)") 'structure = 3-vector,', &
                ' scalar, scalar, scalar, scalar, scalar, scalar, scalar'
           write(ifile,"(a,a)") 'type = float,', &
                'float, float, float, float, float, float, float'
           write(ifile,"(a,a,a)") 'dependency = positions,', &
                ' positions, positions,', &
                ' positions, positions, positions, positions, positions'
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

                    write(ifile,"(10f20.14)") r(1:3), &
                         phipw(i1,i2,i3,ib,iq,isp), &
                         phiaug(i1,i2,i3,ib,iq,isp), &
                         phitot(i1,i2,i3,ib,iq,isp), &
                         abs(phitot(i1,i2,i3,ib,iq,isp))
                 enddo
              enddo
           enddo
           close(ifile)
        enddo
     enddo
  enddo
end subroutine wfn2dx_2
! ccccccccccccccccccccccccccccccccc
subroutine crystal2dx_2(alat,plat,rini,rfin, &
     natom,apos,nclass,iclass,zz)
  implicit none
  double precision ::  alat,plat(3,3),rini(3),rfin(3)
  integer :: natom,nclass,iclass(natom)
  double precision ::  apos(3,natom),zz(nclass)

  ! local
  integer :: nrange(3)
  double precision ::  tmpvec(3),Z(natom)
  integer :: i,j,ic,ia

  do ia=1,natom
     ic=iclass(ia)
     Z(ia)=zz(ic)
  enddo

  nrange(1:3)=2

  call dump_ucell(alat,plat,'ucell.general')
  call dump_box_2(alat,rini,rfin,plat,'box.general')
  call dump_apos_2(alat,plat,rini,rfin,natom,apos, &
       Z,nrange,'apos.general')
  call dump_apos_xyz(alat,plat,rini,rfin,natom,apos, &
       Z,nrange,'apos.xyz')
end subroutine crystal2dx_2
! cccccccccccccccccccccccccccccccccccccccccccccc
subroutine dump_ucell(alat,plat,fn)
  implicit none
  double precision,intent(in) :: alat,plat(3,3)
  character*(*) :: fn
  ! local
  integer :: ifile,i1,i2,i3,ifile_handle
  double precision ::  v1(3),v2(3)
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
! cccccccccccccccccccccccccccccccccccccccccccccc
subroutine dump_apos_2(alat,plat,rini,rfin,natom,apos,Z,nrange,fn)
  implicit none
  double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3), &
       apos(3,natom),Z(natom)
  integer,intent(in) :: natom,nrange(3)
  character*(*) :: fn
  ! local
  integer :: natomall,ifile_handle
  integer :: ifile,i,i1,i2,i3
  double precision ::  v1(3),v2(3)

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
! cccccccccccccccccccccccccccccccccccccccccccccc
subroutine dump_apos_xyz(alat,plat,rini,rfin,natom,apos, &
     Z,nrange,fn)
  implicit none
  double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3), &
       apos(3,natom),Z(natom)
  integer,intent(in) :: natom,nrange(3)
  character*(*) :: fn
  ! local
  integer:: ifile_handle
  integer :: natomall,natomx
  integer :: ifile,i,i1,i2,i3
  double precision ::  v1(3),v2(3),aini(3),afin(3),eps
  double precision,allocatable :: rall(:,:)
  character*2,allocatable :: zall(:)
  character(2) ::  z2element

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
              if ( (v2(1) >= aini(1) .AND. v2(1) <= afin(1)) &
                   .AND. (v2(2) >= aini(2) .AND. v2(2) <= afin(2)) &
                   .AND. (v2(3) >= aini(3) .AND. v2(3) <= afin(3)) ) then
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
! cccccccccccccccccccccccccccccccccccccccccccccc
subroutine dump_box_2(alat,rini,rfin,plat,fn)
  implicit none
  integer:: ifile_handle
  double precision,intent(in) :: alat,plat(3,3),rini(3),rfin(3)
  double precision ::  eps
  character*(*) :: fn
  ! local
  integer :: ifile,i1,i2,i3
  double precision ::  v1(3),v2(3)
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
     if (i1 == 0) then
        v1(1)= -eps
     else
        v1(1)= 1.0d0 + eps
     endif
     do i2=0,1
        if (i2 == 0) then
           v1(2)= -eps
        else
           v1(2)= 1.0d0 + eps
        endif
        do i3=0,1
           if (i3 == 0) then
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
! cccccccccccccccccccccccccccccccccccccccccccccc
character(2) function z2element(Z)
  implicit none
  ! input
  double precision ::  Z
  integer :: iz

  iz=nint(Z)

  if (iz == 1) then
     z2element=' H'
  elseif (iz == 2) then
     z2element='He'
  elseif (iz == 3) then
     z2element='Li'
  elseif (iz == 4) then
     z2element='Be'
  elseif (iz == 5) then
     z2element=' B'
  elseif (iz == 6) then
     z2element=' C'
  elseif (iz == 7) then
     z2element=' N'
  elseif (iz == 8) then
     z2element=' O'
  elseif (iz == 9) then
     z2element=' F'
  elseif (iz == 10) then
     z2element='Ne'
  elseif (iz == 11) then
     z2element='Na'
  elseif (iz == 12) then
     z2element='Mg'
  elseif (iz == 13) then
     z2element='Al'
  elseif (iz == 14) then
     z2element='Si'
  elseif (iz == 15) then
     z2element=' P'
  elseif (iz == 16) then
     z2element=' S'
  elseif (iz == 17) then
     z2element='Cl'
  elseif (iz == 18) then
     z2element='Al'
  elseif (iz == 19) then
     z2element=' K'
  elseif (iz == 20) then
     z2element='Ca'
  elseif (iz == 21) then
     z2element='Sc'
  elseif (iz == 22) then
     z2element='Ti'
  elseif (iz == 23) then
     z2element=' V'
  elseif (iz == 24) then
     z2element='Cr'
  elseif (iz == 25) then
     z2element='Mn'
  elseif (iz == 26) then
     z2element='Fe'
  elseif (iz == 27) then
     z2element='Co'
  elseif (iz == 28) then
     z2element='Ni'
  elseif (iz == 29) then
     z2element='Cu'
  elseif (iz == 30) then
     z2element='Zn'
  elseif (iz == 31) then
     z2element='Ga'
  elseif (iz == 32) then
     z2element='Ge'
  elseif (iz == 33) then
     z2element='As'
  elseif (iz == 34) then
     z2element='Se'
  elseif (iz == 35) then
     z2element='Br'
  elseif (iz == 36) then
     z2element='Kr'
  elseif (iz == 37) then
     z2element='Rb'
  elseif (iz == 38) then
     z2element='Sr'
  elseif (iz == 39) then
     z2element=' Y'
  elseif (iz == 40) then
     z2element='Zr'
  elseif (iz == 41) then
     z2element='Nb'
  elseif (iz == 42) then
     z2element='Mo'
  elseif (iz == 43) then
     z2element='Tc'
  elseif (iz == 44) then
     z2element='Ru'
  elseif (iz == 45) then
     z2element='Rh'
  elseif (iz == 46) then
     z2element='Pd'
  elseif (iz == 47) then
     z2element='Ag'
  elseif (iz == 48) then
     z2element='Cd'
  elseif (iz == 49) then
     z2element='In'
  elseif (iz == 50) then
     z2element='Sn'
  elseif (iz == 51) then
     z2element='Sb'
  elseif (iz == 52) then
     z2element='Te'
  elseif (iz == 53) then
     z2element=' I'
  elseif (iz == 54) then
     z2element='Xe'
  elseif (iz == 55) then
     z2element='Cs'
  elseif (iz == 56) then
     z2element='Ba'
  elseif (iz == 57) then
     z2element='La'
  elseif (iz == 58) then
     z2element='Ce'
  elseif (iz == 59) then
     z2element='Pr'
  elseif (iz == 60) then
     z2element='Nd'
  elseif (iz == 61) then
     z2element='Pm'
  elseif (iz == 62) then
     z2element='Sm'
  elseif (iz == 63) then
     z2element='Eu'
  elseif (iz == 64) then
     z2element='Gd'
  elseif (iz == 65) then
     z2element='Tb'
  elseif (iz == 66) then
     z2element='Dy'
  elseif (iz == 67) then
     z2element='Ho'
  elseif (iz == 68) then
     z2element='Er'
  elseif (iz == 69) then
     z2element='Tm'
  elseif (iz == 70) then
     z2element='Yb'
  elseif (iz == 71) then
     z2element='Lu'
  elseif (iz == 72) then
     z2element='Hf'
  elseif (iz == 73) then
     z2element='Ta'
  elseif (iz == 74) then
     z2element=' W'
  elseif (iz == 75) then
     z2element='Re'
  elseif (iz == 76) then
     z2element='Os'
  elseif (iz == 77) then
     z2element='Ir'
  elseif (iz == 78) then
     z2element='Pt'
  elseif (iz == 79) then
     z2element='Au'
  elseif (iz == 80) then
     z2element='Hg'
  elseif (iz == 81) then
     z2element='Tl'
  elseif (iz == 82) then
     z2element='Pb'
  elseif (iz == 83) then
     z2element='Bi'
  elseif (iz == 84) then
     z2element='Po'
  elseif (iz == 85) then
     z2element='At'
  elseif (iz == 86) then
     z2element='Rn'
  elseif (iz == 87) then
     z2element='Fr'
  elseif (iz == 88) then
     z2element='Ra'
  elseif (iz == 89) then
     z2element='Ac'
  elseif (iz == 90) then
     z2element='Th'
  elseif (iz == 91) then
     z2element='Pa'
  elseif (iz == 92) then
     z2element=' U'
  elseif (iz == 93) then
     z2element='Np'
  elseif (iz == 94) then
     z2element='Pu'
  elseif (iz == 95) then
     z2element='Am'
  elseif (iz == 96) then
     z2element='Cm'
  elseif (iz == 97) then
     z2element='Bk'
  elseif (iz == 98) then
     z2element='Cf'
  elseif (iz == 99) then
     z2element='Es'
  elseif (iz == 100) then
     z2element='Fm'
  elseif (iz == 101) then
     z2element='Md'
  elseif (iz == 102) then
     z2element='No'
  elseif (iz == 103) then
     z2element='Lr'
  endif

  ! debug:
  ! for TiFe2
  !      if (iz.eq.26) z2element=' B' ! Fe
  !      if (iz.eq.22) z2element=' C' ! Ti

  ! for TiFe2
  if (iz == 57) z2element=' N' ! La
  if (iz == 38) z2element=' H' ! Sr
  if (iz == 22) z2element=' C' ! Ti
  if (iz == 13) z2element=' O' ! Al
  if (iz == 8)  z2element=' B' ! O

end function z2element

