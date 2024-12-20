module m_hqpe
  contains
subroutine hqpe() bind(C)  ! Jul,2000 t.kotani started from hqpe by Ferdi.Aryasetiawan
  ! calculates quasiparticle energies
  ! E(k,t) = e(k,t) + Z [SEx(k,t) + SEc(k,t) - xcLDA(k,t)]
  ! e(k,t) = LDA eigenvalue
  ! Z      = [1 - dSEc(e(k,t))/dw]^(-1)
  ! SEx(k,t)   = <psi(k,t)| SEx |psi(k,t)>
  ! SEc(k,t)   = <psi(k,t)| SEc |psi(k,t)>, SEc = GWc
  ! xcLDA(k,t) = <psi(k,t)| vxc |psi(k,t)>
  ! SEx and xcLDA are in file SEX
  ! SEc is in file SEC
  use m_genallcf_v3,only: genallcf_v3,laf
  use m_keyvalue,only: getkeyvalue
  implicit real*8 (a-h,o-z)
  implicit integer(i-n)
  dimension ifsex(2),ifsexcore(2),ifxc(2),ifsec(2),ifqpe(2) &
       ,iftote(2),iftote2(2)
  integer,allocatable :: itxc(:),itc(:),itx(:)
  real(8),allocatable :: qxc(:,:,:),eldaxc(:,:),vxc(:,:), &
       qc(:,:,:),eldac(:,:),sex(:,:),sexcore(:,:), &
       qx(:,:,:),eldax(:,:),rsec(:,:,:),csec(:,:,:),zfac(:,:)
  integer:: ret,iix
  logical :: nozmode=.false.
  call Genallcf_v3(0) 
  jin=0
  open(newunit=ifsex(1)   ,file='SEXU')
  open(newunit=ifsexcore(1) ,file='SEXcoreU')
  open(newunit=ifsec(1)   ,file='SECU')
  open(newunit=ifxc(1)    ,file='XCU')
  call readx   (ifsex(1),50)
  read (ifsex(1),*) nspin,nq,ntq
  if (nspin == 2 .AND. .NOT. laf) then
     open(newunit=ifsex(2)   ,file='SEXD')
     open(newunit=ifsexcore(2),file='SEXcoreD')
     open(newunit=ifsec(2)   ,file='SECD')
     open(newunit=ifxc(2)    ,file='XCD')
  endif
  rewind (ifsex(1))
  !> output file
  open(newunit=ifqpe(1)   ,file='QPU')
  open(newunit=iftote(1)  ,file='TOTE.UP')
  open(newunit=iftote2(1) ,file='TOTE2.UP')
  if (nspin == 2) then
     open(newunit=ifqpe(2)   ,file='QPD')
     open(newunit=iftote(2)  ,file='TOTE.DN')
     open(newunit=iftote2(2) ,file='TOTE2.DN')
  endif
  ! loop over spin
  do      is = 1,nspin
     write(6,*) ' --- is=',is
     ! read dimensions
     call readx   (ifsex(is),50)
     read (ifsex(is),*) nspinx,nqx,ntqx
     !      call readx   (ifsex(is),50)
     read (ifsex(is),*)
     read (ifsex(is),*) deltaw
     read (ifsex(is),*) alat
     read (ifsex(is),*) ef
     call readx(ifsec(is),50)
     read (ifsec(is),*) nspinc,nqc,ntqc
     call readx   (ifxc(is),50)
     read (ifxc(is),*) nspinxc,nqxc,ntqxc
     if (nspin /= nspinx)  call rx( 'hqpe: wrong nspin SEx')
     if (nspin /= nspinc)  call rx( 'hqpe: wrong nspin SEc')
     if (nspin /= nspinxc) call rx( 'hqpe: wrong nspin vxc')
     if (nq /= nqx)        call rx( 'hqpe: wrong nq SEx')
     if (nq /= nqc)        call rx( 'hqpe: wrong nq SEc')
     if (nq /= nqxc)       call rx( 'hqpe: wrong nq vxc')
     if (ntq /= ntqx)      call rx( 'hqpe: wrong ntq SEx')
     if (ntq /= ntqc)      call rx( 'hqpe: wrong ntq SEc')
     if (ntq /= ntqxc)     call rx( 'hqpe: wrong ntq vxc')

     if(is==1) write(6,*)' ###  readin XCU'
     if(is==2) write(6,*)' ###  readin XCD'
     allocate( itxc(ntq),qxc(3,ntq,nq),eldaxc(ntq,nq),vxc(ntq,nq) )
     call readx (ifxc(is),50)
     read(ifxc(is),*)
     do ip = 1,nq
        do i  = 1,ntq
           read(ifxc(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                itxc(i),ipxx,isxxx, qxc(1:3,i,ip), eldaxc(i,ip), &
                vxc(i,ip)
        enddo
     enddo

     if(is==1) write(6,*)' ###  readin SEXU'
     if(is==2) write(6,*)' ###  readin SEXD'
     allocate( itx(ntq), qx (3,ntq,nq),eldax (ntq,nq),sex(ntq,nq) )
     call readx   (ifsex(is),50)
     read(ifsex(is),*)
     do ip = 1,nq
        do i  = 1,ntq
           read(ifsex(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                itx(i),ipxx,isxxx, qx(1:3,i,ip), eldax(i,ip), &
                sex(i,ip)
        enddo
     enddo

     if(is==1) write(6,*)' ###  readin SEXcoreU'
     if(is==2) write(6,*)' ###  readin SEXcoreD'
     allocate( sexcore(ntq,nq) )
     call readx   (ifsexcore(is),50)
     call readx   (ifsexcore(is),50)
     read(ifsexcore(is),*)
     do ip = 1,nq
        do i  = 1,ntq
           read(ifsexcore(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                ixx1,ixx2,ixx3, qxxx1,qxxx2,qxxx3, exxx, sexcore(i,ip)
        enddo
     enddo

     if(is==1) write(6,*)' ###  readin SECU'
     if(is==2) write(6,*)' ###  readin SECD'
     allocate( itc(ntq), qc (3,ntq,nq),eldac (ntq,nq) &
          ,rsec(3,ntq,nq),csec(3,ntq,nq),zfac(ntq,nq))
     call readx   (ifsec(is),50)
     read(ifsec(is),*)
     rsec=0d0
     csec=0d0
     do ip = 1,nq
        do i  = 1,ntq
           if(nozmode) then
              read(ifsec(is),*) &
                   itc(i),ipxxx,isxxx, qc(1:3,i,ip), eldac(i,ip), &
                   rsec(2,i,ip),csec(2,i,ip)
              zfac(i,ip)=1.00
              write(*,*) i,ip,csec(2,i,ip)
           else
              read(ifsec(is),*) &
                   itc(i),ipxxx,isxxx, qc(1:3,i,ip), eldac(i,ip), &
                   rsec(1:3,i,ip),csec(1:3,i,ip),zfac(i,ip)
           endif
        enddo
     enddo
     if (sum(abs(itx(1:ntq)-itc(1:ntq)))/=0)  call rx('hqpe itx /= itc ')
     if (sum(abs(itx(1:ntq)-itxc(1:ntq)))/=0) call rx('hqpe itx /= itxc ')
     call qpe1     (ifqpe(is),iftote(is),iftote2(is),itc,qc, &
          eldac,vxc,sex,sexcore, &
          rsec,csec,zfac,jin,deltaw,alat,ef, &
          ntq,nq,is, &
          eshift0,eshift02,eshtlda)
     deallocate( itxc,qxc,eldaxc,vxc ,itc, qc ,eldac, &
          sexcore ,rsec,csec,zfac, &
          itx, qx ,eldax,sex)
     if (laf) exit
     if (jin > 0) jin = 999999
  end do
  call rx0s( ' OK! hqpe ')
end subroutine hqpe
end module m_hqpe
