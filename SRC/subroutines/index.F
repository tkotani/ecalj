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
      integer function maxocc2 (nspin,ef, nband, qbz,nqbz) ! maximum no. occupied states
      use m_readeigen, only: readeval
      implicit none
      integer(4):: nspin,nqbz,nband,noccx,is,iq,noccxt !,noccx1
      real(8) :: qbz(3,nqbz),ef,ekt(nband, nqbz,nspin )
      noccx      = 0
      do is = 1,nspin
      do iq = 1,nqbz
          ekt(:,iq,is)= readeval(qbz(:,iq),is)
      enddo
      enddo
      maxocc2 = maxval(count(ekt(1:nband,1:nqbz,1:nspin)<ef,1)) !noccx1 (ekt,nqbz,nband,ef)
      end
      
c$$$c$$$!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
c$$$      integer function noccx1 (ekt,nk,nt,ef)
c$$$c 92.02.27
c$$$c find the highest number of occupied states
c$$$c ekt   = eigenvalues for all k-points and states
c$$$c nk    = number of k-points in the 1st BZ
c$$$c nt    = number of states
c$$$c ef    = fermi level
c$$$c noccx1= highest number of occupied states
c$$$      implicit real*8(a-h,o-z)
c$$$      dimension ekt(nt,nk)
c$$$      logical:: allocc
c$$$      integer:: noccx,k,it,nk,nt
c$$$!! Find the highest number of occupied states ---
c$$$!!   Previous version works,But this is necessary when all bands occuied (for model calculaitons). -->I had bug -->see 
c$$$!! aug2017 below (okumura test)
c$$$!!   bugfix (probably affects little only for metallic case). For single band model, it stops since previous version givs noccx1=0.
c$$$      noccx = 0
c$$$      allocc=.true.
c$$$      do   k  = 1,nk !k index
c$$$        do it = 1,nt !band index
c$$$          if(ekt(it,k) > ef) then
c$$$             allocc=.false.
c$$$             if(it-1> noccx) noccx = it-1
c$$$!     bugfix aug2017takao => previous verion check noccx only at it=nt because of bug.
c$$$            exit
c$$$          endif
c$$$        enddo
c$$$      enddo
c$$$      if(allocc) noccx=nt   !completely occupied case. aug2017takao
c$$$      noccx1 = noccx
c$$$      return
c$$$      end



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
