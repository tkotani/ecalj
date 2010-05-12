      subroutine sopert(mode,ndimh,nspc,wk,hin,hout)
C- Manipulates blocks of hamiltonian for noncollinear case
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 for now.  See remarks
Ci   ndimh :hamiltonian dimension
Ci   nspc  :routine does nothing if nspc ne 2
Ci   wk    :double complex work array of dimension (ndimh*nsp)**2
Ci   hin   :blocks of hamiltonian
Co Outputs
Co   hout  :(mode 0) hamiltonian blocks ordered in form suitable for
Co         :         for diagonalization.
Co         :         hout and hin may occupy the same address space.
Co         :           part of hin    spin block of h
Co         :            hin(*,*,1)       (1,1)
Co         :            hin(*,*,2)       (2,2)
Co         :            hin(*,*,3)       (1,2)
Co         :(mode 1) same as mode 0 but (1,2) block assumed to be 0
Cl Local variables
Cl         :
Cr Remarks
Cr   mode 0 orders hamiltonian spin subblocks blocks into full matrix
Cr   mode 1 does the same, but the (1,2) spin block is set to zero.
Cr Checks:
Cr   Combine s11,s22 into large matrix s
Cr   mc s22 s11 -sub 1,nr+nr,1,nc+nc -suba nr/2+1,nr/2+1 -herm -bw s
Cu Updates
Cu   05 Feb 05 Added mode 1
Cu   23 Dec 04 First created (mode 0)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ndimh,nspc
      double complex hin(ndimh,ndimh,*),hout(ndimh,2,ndimh,2)
      double complex wk(ndimh,2,ndimh,2)
C ... Local parameters
      integer ndimhx,i,j,ks,is,js,mapki(3),mapkj(3)
      data mapki /1,2,1/, mapkj /1,2,2/


      if (nspc .ne. 2) return

C     n2 = ndimh**2
      ndimhx = ndimh*nspc
      if (mode .eq. 0 .or. mode .eq. 1) then
        do  ks = 1, 3
          is = mapki(ks)
          js = mapkj(ks)
          if (is .eq. js .or. mode .eq. 0) then
            do  j = 1, ndimh
              do  i = 1, ndimh
                wk(i,is,j,js) = hin(i,j,ks)
              enddo
            enddo
            if (ks .eq. 3) then
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
C       call zprm('h(nc)',2,hout,ndimhx,ndimhx,ndimhx)

      else
        call rxi('sopert: not implemented, mode=',mode)
      endif

      end

