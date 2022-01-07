      integer function maxocc2 (nspin,ef, nband, qbz,nqbz)
c maximum no. occupied states
      use m_readeigen, only: readeval
      implicit none
      integer(4):: nspin,nqbz,nband,noccx,is,iq,noccxt,noccx1
      real(8) :: qbz(3,nqbz),ef
      real(8),allocatable :: ekt(:,:)
      allocate( ekt(nband, nqbz ) )
      noccx      = 0
      do  is = 1,nspin
        do iq = 1,nqbz
          ekt(:,iq)= readeval(qbz(:,iq),is)
        enddo
        noccxt  = noccx1 (ekt,nqbz,nband,ef)
        if (noccxt .gt. noccx) noccx = noccxt
      enddo
      maxocc2 = noccx
      deallocate( ekt )
      end
      
c$$$!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      integer function noccx1 (ekt,nk,nt,ef)
c 92.02.27
c find the highest number of occupied states
c ekt   = eigenvalues for all k-points and states
c nk    = number of k-points in the 1st BZ
c nt    = number of states
c ef    = fermi level
c noccx1= highest number of occupied states
      implicit real*8(a-h,o-z)
      dimension ekt(nt,nk)
      logical:: allocc
      integer:: noccx,k,it,nk,nt
!! Find the highest number of occupied states ---
!!   Previous version works,But this is necessary when all bands occuied (for model calculaitons). -->I had bug -->see 
!! aug2017 below (okumura test)
!!   bugfix (probably affects little only for metallic case). For single band model, it stops since previous version givs noccx1=0.
      noccx = 0
      allocc=.true.
      do   k  = 1,nk !k index
        do it = 1,nt !band index
          if(ekt(it,k) > ef) then
             allocc=.false.
             if(it-1> noccx) noccx = it-1
!     bugfix aug2017takao => previous verion check noccx only at it=nt because of bug.
            exit
          endif
        enddo
      enddo
      if(allocc) noccx=nt   !completely occupied case. aug2017takao
      noccx1 = noccx
      return
      end



c$$$  integer function nocc (ek,ef,n)
c$$$c ek   = eigenvalues
c$$$c ef   = fermi level
c$$$c locc = .true.  ==> occupied
c$$$c      = .false. ==> unoccupied
c$$$c n    = no. states
c$$$      implicit none
c$$$c     logical, intent(in) :: locc
c$$$      integer,intent(in):: n
c$$$      real(8), intent(in) :: ek(n), ef
c$$$      integer :: i
c$$$      nocc       = 0
c$$$      do       i = 1,n
c$$$         if (ek(i) .le. ef) nocc = nocc + 1
c$$$      end do
c$$$      end function nocc
c
c$$$c 92.03.14
c$$$c 92.10.13 from idxlmto (see above) modified for all atoms
c$$$c indexing of LMTO basis functions
c$$$c follows that in TB-LMTO program
c$$$c
c$$$c l    0     1                 2
c$$$c n    1  2  1        2        1               2
c$$$c m    0  0 -1  0 -1 -1  0  1 -2 -1  0  1  2  -2  -1   0   1   2
c$$$c ind  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17
c$$$c
c$$$c il(ind),in(ind),im(ind) = l,n,m
c$$$c ilmto(n,lm,atom) = index of lmto
c$$$c lm = l*l + l + m + 1
