! Bloch transform of real-space matrix  sll(j,i) =\sum_T hrr(T,i,j,isp,T)*phase
! , where phase is the avarage for Tbar (shortest equivalent to T).  exp(ikT)
      subroutine bloch2(qp,isp, sll)
      use m_rdsigm2,only: hrr,ndimsig,nk1,nk2,nk3
      use m_hamindex, only: ib_table
      use m_gennlat_sig,only: nlatS,nlatE,nlat
      use m_lattic,only: plat=>lat_plat
c      use m_mpitk,only: master_mpi,procid
      implicit none
      integer:: isp,iset
      real(8):: qp(3)
      real(8),parameter:: twopi = 8*datan(1d0)
      complex(8):: sll(ndimsig,ndimsig)
      complex(8):: phase,img=(0d0,1d0)
      integer:: ik1,ik2,ik3,i,j,ib1,ib2,ix,nS,nE,nnn
      call tcn('bloch2')
      sll=0d0
      nnn=nk1*nk2*nk3
      do ik1=1,nk1
      do ik2=1,nk2
      do ik3=1,nk3
        do i=1,ndimsig 
        do j=1,ndimsig 
         ib1 = ib_table(i) !atomic-site index in the primitive cell
         ib2 = ib_table(j)
         nS=nlatS(ik1-1,ik2-1,ik3-1,ib1,ib2)
         nE=nlatE(ik1-1,ik2-1,ik3-1,ib1,ib2)
         phase=0d0
         do ix=nS,nE
         phase=phase+exp(-img*twopi*sum(qp*matmul(plat,nlat(:,ix,ib1,ib2))))
         enddo
         phase= phase/(nE-nS+1) 
         sll(i,j) = sll(i,j) + hrr(ik1,ik2,ik3,i,j,isp)/nnn*phase
      enddo
      enddo
      enddo
      enddo
      enddo
      call tcx('bloch2')
      end subroutine bloch2
