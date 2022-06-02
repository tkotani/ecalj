!      subroutine sopert(mode,ndimh,nspc,wk,hin,hout)
subroutine sopert2(hin,hout)
  use m_igv2x,only:   ndimh
  use m_lmfinit,only: nspc
  !- Manipulates blocks of hamiltonian for noncollinear case
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 for now.  See remarks
  !i   ndimh :hamiltonian dimension
  !i   nspc  :routine does nothing if nspc ne 2
  !i   wk    :double complex work array of dimension (ndimh*nsp)**2
  !i   hin   :blocks of hamiltonian
  !o Outputs
  !o   hout  :(mode 0) hamiltonian blocks ordered in form suitable for
  !o         :         for diagonalization.
  !o         :         hout and hin may occupy the same address space.
  !o         :           part of hin    spin block of h
  !o         :            hin(*,*,1)       (1,1)
  !o         :            hin(*,*,2)       (2,2)
  !o         :            hin(*,*,3)       (1,2)
  !o         :(mode 1) same as mode 0 but (1,2) block assumed to be 0
  !l Local variables
  !l         :
  !r Remarks
  !r   mode 0 orders hamiltonian spin subblocks blocks into full matrix
  !r   mode 1 does the same, but the (1,2) spin block is set to zero.
  !r Checks:
  !r   Combine s11,s22 into large matrix s
  !r   mc s22 s11 -sub 1,nr+nr,1,nc+nc -suba nr/2+1,nr/2+1 -herm -bw s
  !u Updates
  !u   05 Feb 05 Added mode 1
  !u   23 Dec 04 First created (mode 0)
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer:: mode=0 !,ndimh,nspc
  double complex hin(ndimh,ndimh,*),hout(ndimh,2,ndimh,2)
  double complex wk(ndimh,2,ndimh,2)
  ! ... Local parameters
  integer :: ndimhx,i,j,ks,is,js,mapki(3),mapkj(3)
  data mapki /1,2,1/, mapkj /1,2,2/


  if (nspc /= 2) return

  !     n2 = ndimh**2
  ndimhx = ndimh*nspc
  if (mode == 0 .OR. mode == 1) then
     do  ks = 1, 3
        is = mapki(ks)
        js = mapkj(ks)
        if (is == js .OR. mode == 0) then
           do  j = 1, ndimh
              do  i = 1, ndimh
                 wk(i,is,j,js) = hin(i,j,ks)
              enddo
           enddo
           if (ks == 3) then
              do  j = 1, ndimh
                 do  i = 1, ndimh
                    wk(j,js,i,is) = dconjg(hin(i,j,ks))
                 enddo
              enddo
           endif
        else
           do  j = 1, ndimh
              do  i = 1, ndimh
                 wk(i,is,j,js) = (0d0,0d0)
                 wk(i,js,j,is) = (0d0,0d0)
              enddo
           enddo
        endif
     enddo
     call dcopy(ndimhx**2*2,wk,1,hout,1)
     !       call zprm('h(nc)',2,hout,ndimhx,ndimhx,ndimhx)

  else
     call rxi('sopert: not implemented, mode=',mode)
  endif

end subroutine sopert2


subroutine sopert(mode,ndimh,nspc,wk,hin,hout)
  !- Manipulates blocks of hamiltonian for noncollinear case
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 for now.  See remarks
  !i   ndimh :hamiltonian dimension
  !i   nspc  :routine does nothing if nspc ne 2
  !i   wk    :double complex work array of dimension (ndimh*nsp)**2
  !i   hin   :blocks of hamiltonian
  !o Outputs
  !o   hout  :(mode 0) hamiltonian blocks ordered in form suitable for
  !o         :         for diagonalization.
  !o         :         hout and hin may occupy the same address space.
  !o         :           part of hin    spin block of h
  !o         :            hin(*,*,1)       (1,1)
  !o         :            hin(*,*,2)       (2,2)
  !o         :            hin(*,*,3)       (1,2)
  !o         :(mode 1) same as mode 0 but (1,2) block assumed to be 0
  !l Local variables
  !l         :
  !r Remarks
  !r   mode 0 orders hamiltonian spin subblocks blocks into full matrix
  !r   mode 1 does the same, but the (1,2) spin block is set to zero.
  !r Checks:
  !r   Combine s11,s22 into large matrix s
  !r   mc s22 s11 -sub 1,nr+nr,1,nc+nc -suba nr/2+1,nr/2+1 -herm -bw s
  !u Updates
  !u   05 Feb 05 Added mode 1
  !u   23 Dec 04 First created (mode 0)
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,ndimh,nspc
  double complex hin(ndimh,ndimh,*),hout(ndimh,2,ndimh,2)
  double complex wk(ndimh,2,ndimh,2)
  ! ... Local parameters
  integer :: ndimhx,i,j,ks,is,js,mapki(3),mapkj(3)
  data mapki /1,2,1/, mapkj /1,2,2/


  if (nspc /= 2) return

  !     n2 = ndimh**2
  ndimhx = ndimh*nspc
  if (mode == 0 .OR. mode == 1) then
     do  ks = 1, 3
        is = mapki(ks)
        js = mapkj(ks)
        if (is == js .OR. mode == 0) then
           do  j = 1, ndimh
              do  i = 1, ndimh
                 wk(i,is,j,js) = hin(i,j,ks)
              enddo
           enddo
           if (ks == 3) then
              do  j = 1, ndimh
                 do  i = 1, ndimh
                    wk(j,js,i,is) = dconjg(hin(i,j,ks))
                 enddo
              enddo
           endif
        else
           do  j = 1, ndimh
              do  i = 1, ndimh
                 wk(i,is,j,js) = (0d0,0d0)
                 wk(i,js,j,is) = (0d0,0d0)
              enddo
           enddo
        endif
     enddo
     call dcopy(ndimhx**2*2,wk,1,hout,1)
     !       call zprm('h(nc)',2,hout,ndimhx,ndimhx,ndimhx)

  else
     call rxi('sopert: not implemented, mode=',mode)
  endif

end subroutine sopert

