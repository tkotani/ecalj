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
      maxocc2 = maxval(count(ekt(1:nband,1:nqbz,1:nspin)<ef,1)) 
      end
      
