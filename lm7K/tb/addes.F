      subroutine addes(nsites,nl,nsp,nbas,ltb,npr,iax,ov,dov,ediag,
     .  hrs,dh)
C- Add to hamiltonian and derivative MEs: hLL' -> hLL' + ebarLL' * sLL'
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl,nsp,nbas,ltb
Ci   nsites: total number of neighbors in all clusters
Ci   npr(1,i): offset in accumulated list of all neighbors to the
Ci             ith cluster associated with the ith atom
Ci   iax(1,j) = cental atom of cluster
Ci   iax(2,j) = second atom of the pair
Ci   ov: overlap matrix elements
Ci   dov: overlap derivative matrix elements
Ci   ediag: work array for M-averaged site energies
Co Outputs
Co   hrs: ebarLL' * sLL' added to hamiltonian MEs hLL'
Co   dh: add corresponding terms to hamiltonian derivatives if needed
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nsites,nl,nsp,nbas,ltb,niax
      parameter (niax=10)
      integer npr(0:1,1),iax(niax,nsites)
      double precision ov(nl**2,nl**2,nsites*nsp**2),
     .  dov(nl**2,nl**2,nsites*nsp**2,4),ediag(0:nl-1,nbas),
     .  hrs(nl**2,nl**2,nsites*nsp**2),dh(nl**2,nl**2,nsites*nsp**2,4)
C Local parameters
      double precision e0,ebar
      integer isp,ib1,ib2,j,l1,l2,ll,i,ii
      logical hdiffs,cryf,bittst

      hdiffs = bittst(ltb,16) .or. bittst(ltb,128)
      cryf   = bittst(ltb,2)
      if (bittst(ltb,2**18)) call rx('ADDES: not ready for TB+U')
      if (hdiffs .and. cryf) print 100
  100 format(/' (warning) forces do not include crystal field ',
     .  'terms in off-diagonal MEs')

      do  60  isp = 1, nsp
C --- Accumulate M-averaged site energies ---
        call dpzero(ediag,nl*nbas)
        do  20  ib1 = 1, nbas
          j = npr(1,ib1) + (isp-1)*3*nsites + 1
          do  10  l1 = 1, nl**2
            ediag(ll(l1),ib1) = ediag(ll(l1),ib1)
     .                        + hrs(l1,l1,j) / (2*ll(l1)+1)
   10     continue
   20   continue

C --- Add to hamiltonian and (if needed) derivative matrix elements ---
        do  50  i = 1, nsites
          ii = i + (isp-1)*3*nsites
          ib1 = iax(1,i)
          ib2 = iax(2,i)
          j = npr(1,ib1) + (isp-1)*3*nsites + 1
          if (ii .eq. j) goto 50
          do  40  l1 = 1, nl**2
            e0 = ediag(ll(l1),ib1)
            do  30  l2 = 1, nl**2
              ebar = (e0 + ediag(ll(l2),ib2)) / 2
              hrs(l1,l2,ii) = hrs(l1,l2,ii) + ebar*ov(l1,l2,ii)
              if (hdiffs) then
                dh(l1,l2,ii,1) = dh(l1,l2,ii,1) + ebar*dov(l1,l2,ii,1)
                dh(l1,l2,ii,2) = dh(l1,l2,ii,2) + ebar*dov(l1,l2,ii,2)
                dh(l1,l2,ii,3) = dh(l1,l2,ii,3) + ebar*dov(l1,l2,ii,3)
                dh(l1,l2,ii,4) = dh(l1,l2,ii,4) + ebar*dov(l1,l2,ii,4)
              endif
   30       continue
   40     continue
   50   continue
   60 continue

      end
