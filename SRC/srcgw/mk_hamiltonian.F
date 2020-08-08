      subroutine mk_hamiltonian(nbas,s_atom,s_spec,s_lat,q,
     .   k1,k2,k3,smpot,vconst,osig,otau,oppi,alfa,ndimh,h,s,
     &    kmaxx,nlmax,bmat)

C =================================================================
C- Packages the steps to make the hamiltonian given the potential.
C  Input:
C    q               Bloch vector
C    smpot           smooth potential on mesh
C    vconst          constant potential term
C    osig,otau,oppi  augmentation matrices
C  Output:
C    h               Hamiltonian matrix
C    s               overlap matrix
C =================================================================
      implicit real*8 (a-h,p-z), integer(o)
      dimension s_atom(1),s_spec(1),s_lat(1),q(3),
     .   osig(3,nbas),otau(3,nbas),oppi(3,nbas)
      complex*16 smpot(k1,k2,k3),h(ndimh,ndimh),s(ndimh,ndimh)

      integer(4):: kmaxx,nlmax
      complex(8):: bmat(*)

      call tcn('mk_hamiltonian')
      call u_lat_FT(s_lat,n1,n2,n3)
      call dpzero(h,ndimh*ndimh*2)
      call dpzero(s,ndimh*ndimh*2)

C --- Augmentation part of hamiltonian and overlap
      call augm_q(s_atom,s_spec,s_lat,q,osig,otau,oppi,ndimh,h,s,
     &    kmaxx,nlmax,nbas,bmat)

      if (dabs(alfa).gt.1d-10) call pvmkhm(ndimh,alfa,s,h)

C --- Smoothed hamiltonian and overlap, S and T done analytically
      vavg=vconst
      call smhs_q(s_atom,s_spec,s_lat,vavg,q,ndimh,h,s)
      call hsi_q(s_atom,s_spec,s_lat,n1,n2,n3,k1,k2,k3,smpot,
     .   q,ndimh,h)

C --- Alternatively, smoothed H and S done completely numerically
c|      call hsmi_q(s_atom,s_spec,s_lat,n1,n2,n3,k1,k2,k3,smpot,
c|     .   vconst,q,ndimh,h,s)


c --- zero out matrices exactly for distant pairs
      call setzmat(s_atom,s_spec,s_lat,ndimh,h,s)

c --- fill second half of matrix
      do i=1,ndimh
        do j=i,ndimh
          h(j,i)=dconjg(h(i,j))
          s(j,i)=dconjg(s(i,j))
        enddo
c ...   make sure diagonal elements are exactly zero
        h(i,i)=dcmplx(dreal(h(i,i)),0d0)
        s(i,i)=dcmplx(dreal(s(i,i)),0d0)

      enddo

      call tcx('mk_hamiltonian')
      end


c --- pv routine: add alfa times Sloc-squared to H
      subroutine pvmkhm(ndimh,alfa,s,h)
      implicit real*8 (a-h,p-z), integer(o)
      complex*16 s(ndimh,ndimh),h(ndimh,ndimh),csum
      call tcn('pvmkhm')

      do j=1,ndimh
        do i=1,j
          csum=0d0
          do k=1,ndimh
            csum=csum+dconjg(s(k,i))*s(k,j)
          enddo
          h(i,j)=h(i,j)+alfa*csum
          h(j,i)=dconjg(h(i,j))
        enddo
      enddo

      call tcx('pvmkhm')

      end
