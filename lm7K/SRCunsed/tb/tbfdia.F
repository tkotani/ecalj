      subroutine tbfdia(nbas,nl,nsp,nclass,nsites,npr,ip2,ipc,dh)
C- Accumulate parameter derivatives for diagonal hamiltonian MEs
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas,nl,nsp,nclass,ipc
Ci   nsites, total number of neighbors in all clusters
Ci   npr(1,i): offset in accumulated list of all neighbors to the
Ci             ith cluster associated with the ith atom
Ci   ip2: pointer to locations in full list of variables
Co Outputs
Co   dh: deriv. wrt TB parameters for all sites, diagonal parts set up
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,nl,nsp,nclass,nsites
      integer npr(0:1,nbas),ip2(nl,nsp,nclass,5:6),ipc(nbas)
      double precision dh(nl**2,nl**2,nsites*nsp**2,1)
C Local parameters
      integer isp,i,j,l,iv,ll

      do  10  isp = 1, nsp
        do  10  i = 1, nbas
          j = npr(1,i) + (isp-1)*3*nsites + 1
          do  10  l = 1, nl**2
            iv = ip2(ll(l)+1,isp,ipc(i),5)
            if (iv .eq. 0) goto 10
            dh(l,l,j,iv) = 1d0
   10 continue

      end
