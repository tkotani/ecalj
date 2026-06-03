module m_w_psir_unusednow
contains
  subroutine psir_minmax(n1,n2,n3,k1,k2,k3,f, emin,emax,norm,imagmax,imaxi)
    implicit none
    integer,intent(in):: n1,n2,n3,k1,k2,k3
    complex(8),intent(in):: f(k1,k2,k3)
    real(8),intent(out):: emin(2),emax(2)
    real(8),intent(out),optional:: norm,imagmax
    integer,intent(out),optional:: imaxi(3)
    integer:: i1,i2,i3
    i1=1;i2=1;i3=1
    emin(1)=dreal(f(i1,i2,i3))
    emin(2)=dimag(f(i1,i2,i3))
    emax(1)=dreal(f(i1,i2,i3))
    emax(2)=dimag(f(i1,i2,i3))
    if (present(imagmax)) imagmax=abs(dimag(f(i1,i2,i3)))
    if (present(imaxi)) imaxi=(/i1,i2,i3/)
    if (present(norm)) norm=0.0d0
    do i3=1,n3
       do i2=1,n2
          do i1=1,n1
             if (emin(1)> dreal(f(i1,i2,i3))) emin(1)=dreal(f(i1,i2,i3))
             if (emax(1)< dreal(f(i1,i2,i3))) emax(1)=dreal(f(i1,i2,i3))
             if (emin(2)> dimag(f(i1,i2,i3))) emin(2)=dimag(f(i1,i2,i3))
             if (emax(2)< dimag(f(i1,i2,i3))) emax(2)=dimag(f(i1,i2,i3))
             if (present(imagmax) .AND. present(imaxi) ) then
                if (imagmax < abs(dimag(f(i1,i2,i3))) ) then
                   imaxi=(/i1,i2,i3/)
                   imagmax=abs(dimag(f(i1,i2,i3)))
                endif
             endif
             if (present(norm)) norm=norm+dreal(f(i1,i2,i3))**2+dimag(f(i1,i2,i3))**2
          enddo
       enddo
    enddo
  end subroutine psir_minmax
  subroutine w_psir_xcrysden(ifi,plat,alat,nbas,pos,z,n1,n2,n3,k1,k2,k3,psi)
    implicit none
    real(8),intent(in):: plat(3,3),alat,pos(3,nbas),z(nbas)
    integer,intent(in):: nbas ,ifi
    integer,intent(in):: n1,n2,n3,k1,k2,k3
    complex(8),intent(in):: psi(k1,k2,k3)
    integer:: isp,i,i1,i2,i3
    character(10):: psiform='(6E15.5)'
    write(ifi,'("CRYSTAL")')
    write(ifi,'("PRIMVEC")')
    write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3) &
         ,i2=1,3)
    write(ifi,'("PRIMCOORD")')
    write(ifi,'(2i5)') nbas,1
    do i = 1, nbas
       write(ifi,'(i4,2x,3f10.5)') int(z(i)), &
            (pos(i2,i)*alat*0.529177208,i2=1,3)
    enddo
    do isp=1,3
       write(ifi,'("BEGIN_BLOCK_DATAGRID_3D")')
       write(ifi,'(" wavefunction_complex_",i1)') isp
       write(ifi,'(" BEGIN_DATAGRID_3D_complex_",i1)') isp
       write(ifi,'(3i4)') n1,n2,n3
       write(ifi,'(3F10.5)') 0.,0.,0.
       write(ifi,'(3F10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3) &
            ,i2=1,3)
       if (isp == 1) then
          do i3=1,n3
             do i2=1,n2
                write(ifi,psiform) (dble(psi(i1,i2,i3)),i1=1,n1)
             enddo
          enddo
       elseif (isp == 2) then
          do i3=1,n3
             do i2=1,n2
                write(ifi,psiform) (dimag(psi(i1,i2,i3)),i1=1,n1)
             enddo
          enddo
       elseif (isp == 3) then
          do i3=1,n3
             do i2=1,n2
                write(ifi,psiform)(dimag(psi(i1,i2,i3))**2+dreal(psi(i1,i2,i3))**2,i1=1,n1)
             enddo
          enddo
       endif
       write(ifi,'(" END_DATAGRID_3D")')
       write(ifi,'("END_BLOCK_DATAGRID_3D")')
    enddo
  end subroutine w_psir_xcrysden
  !!
  subroutine w_psir(ng,nspc,nev,psi,n1,n2,n3,k1,k2,k3,kv,isp,q,iq,n_eiglist,eiglist,plat,alat,nbas,pos,z,f) ! FT wave function to real space and add square into mesh density
    !  and optionally make smpot * psi
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   ng    :number of G-vectors
    !i   nspc  :2 for coupled spins; otherwise 1
    !i   nev   :number of wave functions
    !i   psi   :wave function in reciprocal space
    !i   n1..3 :size of FT mesh
    !i   k1..3 :dimensions smpot,smrho
    !i   kv    :indices for gather/scatter operations (gvlist.f)
    !i   isp  : spin index
    !i   n_eiglist: number of  eiglist
    !i   eiglist: eig list to print psir
    !i   plat: lattice  vector
    !i   alat: unit of plat
    !i   nbas: # of atoms
    !i   pos: postion of atoms
    !i   z:   atomic number of atoms
    !o Work
    !o   f     :psi in real space
    !r Remarks
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer,intent(in):: k1,k2,k3,n1,n2,n3,ng,nspc,nev,kv(ng,3)
    integer,intent(in):: isp  ! spin ?
    complex(8),intent(in):: psi(ng,nspc,nev)
    complex(8)::f(k1,k2,k3)
    integer,intent(in)::  n_eiglist
    integer,intent(in):: eiglist(n_eiglist)
    integer,intent(in):: nbas
    real(8),intent(in):: plat(3,3),alat,pos(3,nbas),z(nbas),q(3)
    character(6):: thisfunc='w_psir'
    integer :: i1,i2,i3,i,ispc,ieig,iq
    logical:: l_execute
    character(40):: fnamer,fnamei
    integer:: ifiler,ifilei,stat
    real(8):: emin(2),emax(2),norm,imagmax
    complex(8):: zphase
    integer:: imaxi(3)

    logical:: l_rot=.false.
    real(8),parameter:: zero=0.0d0
    write(*,*) thisfunc,':start '
    if (q(1) == zero .AND. q(2) == zero .AND. q(3) == zero ) then
       l_rot=.true.
       write(*,*)thisfunc,'rotation is done for gamma point'
    endif
    !$$$#if 0
    !$$$      write(*,'(a,9F15.5)')'plat=',plat
    !$$$      write(*,'(a,F15.5)')'alat=',alat
    !$$$      write(*,'(a,i5)')'nbas=',nbas
    !$$$      write(*,*)'pos='
    !$$$      write(*,'(20F15.5)')pos
    !$$$      write(*,*)'z='
    !$$$      write(*,'(20F15.5)')z
    !$$$#endif
    !$$$
    do  ispc = 1, nspc
       do ieig=1,n_eiglist
          i=eiglist(ieig)
          write(fnamer,'(a,I1,a,I1,a,I4.4,a,I4.4,a)') 'psi_s',isp,'_',ispc,'_q',iq,'_e',i,'.xsf'
          write(*,*)thisfunc,':open ',trim(fnamer),' for psi'
          open(newunit=ifiler,file=fnamer,status='unknown',iostat=stat,action='write')
          if (stat /= 0) then
             write(*,*)thisfunc,':Error: failed to open ',fnamer,'. in ',thisfunc, &
                  ' but continue.'
             return
          endif
          call gvputf(ng,1,kv,k1,k2,k3,psi(1,ispc,i),f)
          call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,1)
          call psir_minmax(n1,n2,n3, k1,k2,k3, f, &
               emin,emax,norm,imagmax,imaxi)
          write(*,'(a,a,4i4)') thisfunc,':ispc,ispc,iq,ie=',isp,ispc,iq,i
          write(*,'(a,a,3F15.7)') thisfunc,':real part, min,max and norm=',emin(1),emax(1),norm
          write(*,'(a,a,3F15.7)')thisfunc, ':imag part, min,max,|max|   =',emin(2),emax(2),imagmax
          zphase=f( imaxi(1),imaxi(2),imaxi(3) )
          if ( abs(dimag(zphase)) > abs(dreal(zphase))*1d-3 ) then
             zphase=dconjg(zphase)
             zphase=zphase / abs(zphase)
             do i3=1,n3
                do i2=1,n2
                   do i1=1,n1
                      f(i1,i2,i3)=f(i1,i2,i3)*zphase
                   enddo
                enddo
             enddo
             call psir_minmax(n1,n2,n3, k1,k2,k3, f, &
                  emin,emax,norm,imagmax,imaxi)
             write(*,'(a,a,4i4)') thisfunc,':isp,ispc,iq,ie=',isp,ispc,iq,i
             write(*,'(a,a,3F15.7)') thisfunc,':rot real part, min,max and norm=',emin(1),emax(1),norm
             write(*,'(a,a,3F15.7)') thisfunc,':rot imag part, min,max,|max|   =',emin(2),emax(2),imagmax
          endif
          write(ifiler,'(a,4i5,3F15.5)') '#isp,ispc,iq,ie,q=',isp,ispc,iq,i,q
          call w_psir_xcrysden(ifiler,plat,alat,nbas,pos,z &
               ,n1,n2,n3,k1,k2,k3, f  )
          close(ifiler)
       enddo
    enddo
  end subroutine w_psir
end module m_w_psir_unusednow
