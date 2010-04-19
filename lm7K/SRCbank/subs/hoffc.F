      subroutine hoffc(nbas,ntab,iax,iclass,ioffb,lmx,iprs,iofc,
     .  iofh,ndimh)
C- Makes a table of offsets in cluster hamiltonian to isite
Ci Inputs:
Ci   nbas,ntab,iax,iclass,ioffb,lmx
Ci   iprs,iofc: work arrays
Co Outputs:
Co   iofh: offset in h of the block for isite
Co   ndimh: total dimension needed to store hamiltonian
Cr Remarks
Cr   should be replaced with hoffcl
C     implicit none
      integer nbas,ntab(1),niax,iclass(56),lmx(2),iprs(15),
     .  iofc(15),ioffb(nbas),iofh(1),ndimh
      parameter (niax=10)
      integer iax(niax,20)
      integer it,ib,jb,is,js,nlmi,nlmj,nhs,nt,kb,ic,nclus

      ndimh = 0
      nt = ntab(nbas+1)

C --- Offsets for mode=0  ---
      nhs = ioffb(nbas+1)
      do  10  kb = 1, nbas
        ic  = 1+ntab(kb)
        do  11  it = ic, ntab(kb+1)
   11   iprs(it-ntab(kb)) = iclass(iax(2,it))
        nclus = ntab(kb+1) - ntab(kb) 
        call hoffs(nclus,lmx,iprs,iofc)
        do  10  it = ic, ntab(kb+1)
          ib = iax(1,it)
          jb = iax(2,it)
C          is = iclass(ib)
C          js = iclass(jb)
C          nlmi = iofc(2)-iofc(1)
C          nlmj = iofc(it-ic+2)-iofc(it-ic+1)
          iofh(it) = ioffb(ib) + nhs*ioffb(jb)
C          print 333, it,ib,jb,iofh(it),nlmi*nlmj
C  333     format(' it=',i4,'  ib=',i4,'  jb=',i4,
C     .      '  ndimh=',i9,i6)
   10 continue
      ndimh = nhs**2
      
      end

