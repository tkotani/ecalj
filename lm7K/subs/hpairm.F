      subroutine hpair(nbas,ips,alat,plat,pos,rham,
     .  nttab,ntab,iax,rtab,mxcsiz,ordrd)
C- Neighbor table for Hamiltonian matrix (molecules version)
C  Inputs:  nbas,ips,alat,plat,pos,rham
C  Outputs:
C    nttab: total number of pairs within range dr <= rham(i)+rham(j)
C    iax(1:2,*) sites ib,jb of each pair
C    iax(3:5,*) multiples of plat connecting sites (zero for molecules)
C    iax(6,*)   points to conjugate (jb,ib) pair
C    iax(7,*)   permutation table ordering cluster by increasing ib
C    ntab(1:nbas+1):  offset to first element in iax and rtab for ib
C    rtab(1:3,*): connecting vector defined as pos(jb)-pos(ib)
C    mxcsiz: number of atoms in largest cluster (for dimensioning)
C    ordrd: 0 order each cluster by increasing jb; 1 by increasing dist.
C     implicit none
      integer nbas,nttab,iax(10,1),ips(nbas),ntab(nbas+1),mxcsiz,ordrd
      double precision alat,plat(3,3),pos(3,nbas),rham(nbas),rtab(3,1)
      integer ib,is,jb,js,ipr,nt,mtab,i,owk1,owk2,owk3,it,jt
      double precision r1,r2,rr
      real w(1)
      common /w/ w

C      print *, 'comment out setpr'
C      call setpr(41)

      call getpr(ipr)
      if (ipr .gt. 40) print 332
  332 format(/' hpair:     ib   jb',9x,
     .  '------- dr --------',11x,'d      r1+r2')

C --- Accumulate all pairs for which dr < ri + rj ---
      nttab = 1
      ntab(1) = 0
      mtab = 1
      mxcsiz = 0
      do  10  ib = 1, nbas
      is = ips(ib)
      r1 = rham(is)
      do  20  jb = 1, nbas
        js = ips(jb)
        r2 = rham(js)
C ...   Add to list if connecting vector < r1+r2
        rtab(1,nttab) = alat*(pos(1,jb)-pos(1,ib))
        rtab(2,nttab) = alat*(pos(2,jb)-pos(2,ib))
        rtab(3,nttab) = alat*(pos(3,jb)-pos(3,ib))
        rr = dsqrt(rtab(1,nttab)**2+rtab(2,nttab)**2+rtab(3,nttab)**2)
        if (rr .lt. r1+r2) then
          iax(1,nttab) = ib
          iax(2,nttab) = jb
          iax(3,nttab) = 0
          iax(4,nttab) = 0
          iax(5,nttab) = 0
          nttab = nttab+1
        endif
   20 continue

C --- Sort table by increasing length ---
      if (ordrd .eq. 1) then
        call defrr(owk1, (nttab-mtab)*6)
        call defrr(owk2, (nttab-mtab))
        call defrr(owk3, (nttab-mtab)*3)
        call xhpair(nttab-mtab,iax(1,mtab),rtab(1,mtab),
     .    w(owk1),w(owk2),w(owk3))
        call rlse(owk1)
      endif

C --- Printout for this ib ---
      if (ipr .ge. 30) print 345, ib, nttab-mtab
  345 format(' hpair, ib=',i3,':',i4,' neighbors')
      if (ipr .gt. 40) then
        do  12  nt = mtab, nttab-1
        rr = dsqrt(rtab(1,nt)**2 + rtab(2,nt)**2 + rtab(3,nt)**2)
        r1 = rham(ips(iax(1,nt)))
        r2 = rham(ips(iax(2,nt)))
        print 333, iax(1,nt),iax(2,nt),(rtab(i,nt),i=1,3), rr, r1+r2
C    .    ,iax(3,nt),iax(4,nt),iax(5,nt)
     .    ,iax(7,nt)
  333   format(i14,i5,3f11.6,2f9.4,2x,3i3)
   12   continue
      endif

C --- Cleanup for this ib ---
      mtab = nttab
      ntab(ib+1) = nttab-1
      mxcsiz = max(mxcsiz,ntab(ib+1)-ntab(ib))
   10 continue
      nttab = nttab-1
C     call awrit2('xx ntab %n:1i',' ',80,6,nbas+1,ntab)
      if (ipr.ge.20 .and. ipr .le. 40) print
     .  '('' hpair:'',i5,'' pairs total'',i5,'' is max cluster size'')',
     .  nttab, mxcsiz
C --- For each pair, find matching pair ---
      do  74  it = 1,  nttab
   74 iax(6,it) = 0
      do  70  ib = 1, nbas
      do  70  it = ntab(ib)+1, ntab(ib+1)
        if (iax(6,it) .ne. 0) goto 70
        jb = iax(2,it)
        do  72  jt = ntab(jb)+1, ntab(jb+1)
          if (iax(2,jt) .eq. ib .and. 
     .        iax(3,it) .eq. -iax(3,jt) .and.
     .        iax(4,it) .eq. -iax(4,jt) .and.
     .        iax(5,it) .eq. -iax(5,jt))  then
            iax(6,it) = jt
            iax(6,jt) = it
            goto 73
          endif
   72   continue
        call rx('bug in hpair')
   73   continue
   70 continue

C     do  80  ib = 1, nbas
C  80 print '(99i3)', (iax(6,i), i=ntab(ib)+1, ntab(ib+1))
c     call rx('hpairm')
      end
      subroutine xhpair(nttab,iax,rtab,iwk,iwk2,rwk)
C- Sort neighbor table by distance
C     implicit none
      integer nttab,iax(10,nttab),iwk(6,nttab),iwk2(nttab),i,j,k
      double precision rtab(3,nttab),rwk(3,nttab)

      do  10  i = 1, nttab
        rwk(1,i) = rtab(1,i)
        rwk(2,i) = rtab(2,i)
        rwk(3,i) = rtab(3,i)
        do  12  k = 1, 6
   12   iwk(k,i) = iax(k,i)
   10 continue
      call dvshel(3,nttab,rtab,iwk2,11)
      do  20  i = 1, nttab
        j = iwk2(i)+1
        rtab(1,i) = rwk(1,j)
        rtab(2,i) = rwk(2,j)
        rtab(3,i) = rwk(3,j)
        do  22  k = 1, 6
   22   iax(k,i) = iwk(k,j)
        iax(7,j) = i
   20 continue
      end
      subroutine pairmc(kb,nclus3,iax,iax3,ntab,rtab,ipsite)
C- Table of pointers to pair table for cluster around a third site
C  Output: ipsite(ic,jc)
C  Three-center terms arising from augmentation site kc contribute to
C  the hamiltonian for each pair (ic,jc) that connects to kc.
C  is=ipsite(ic,jc) locates the iax(*,is) corresponding to (ic,jc).
C
C  For efficient implementation, this routine relies on a sorted iax
C  table.  iax(7,*) is a table of permutations ip(*) that sorts each
C  local cluster by basis atom ib (see hpairm).  A sorted
C  iax table guarantees that ib(ip(i+1)) > ib(ip(i)); therefore
C  ipsite(ib(ip(ic)),jb(ip(jc+1))) > ipsite(ib,jb(ip(jc))).
C     implicit none
      integer kb,nclus3,iax(10,9),iax3(nclus3),ntab(1),
     .  ipsite(nclus3,nclus3)
      double precision rtab(3,1)
      integer kc,i,j,ii,jj,ip,jp,ib,jb,nc,ic,icp,ic0,is,nclus,ipr
      integer ipx(200),ipxi(200)

      call getpr(ipr)
      if (nclus3 .gt. 200) call rx('pairmc: increase dim of ipx')
      do  10  ii = 1, nclus3**2
   10 ipsite(ii,1) = 0
      kc = 1+ntab(kb)

C --- Make permutation table for all elts in 3C range ---
      nclus = ntab(kb+1)-ntab(kb)
      i = 0
      do  14  ii = 1, nclus
        ip = iax(7,kc-1+ii)
C   ... if ip is a member of iax(3), keep it; else go to next
        do  16  j = 1, nclus3
          if (ip .eq. iax3(j)) then
            i = i+1
            ipx(i) = ip
            ipxi(i) = j
            goto 17
          endif
   16   continue
   17   continue
   14 continue

C      print 345, (iax3(i), i=1,nclus3)
C      print 345, (iax(7,kc-1+i), i=1,nclus)
C      print 345, (ipx(i), i=1,nclus3)
C      print 345, (ipxi(i), i=1,nclus3)
C      print 345, (iax(2,kc-1+ipx(i)), i=1,nclus3)
C  345 format(34i4)

C --- Loop over 3C pairs i,j in permuted order ---
      do  30  ii = 1, nclus3
        ip = ipx(ii)
        ib = iax(2,kc-1+ip)
C   ... ic is offset corresponding to iax(1,ic) = 1st pair for ib
        ic0 = ntab(ib)
        ic  = ic0
        nc  = ntab(ib+1)
        do  40  jj = 1, nclus3
          jp = ipx(jj)
C     ... Do upper triangle only
          if (jp .lt. ip) goto 40
          jb = iax(2,kc-1+jp)
C     ... Increment ic until this (ib,jb) pair is found
   42     ic = ic+1
          icp = ic0+iax(7,ic)
C     ... Skip if missing pair corresponding to this 3C (ib,jb)
          if (ic .gt. nc) goto 40
C          call awrit5('want ib,jb=%i %i, site %i has %i %i',
C     .      ' ',80,6,ib,jb,icp,ib,iax(2,icp))
C          if (ic .gt. nc .or. icp .lt. 0) then
C            call awrit2('no pair ib,jb=%i %i ...',' ',80,6,ib,jb)
C            pause
C            goto 40
C          endif
          if (jb .ne. iax(2,icp)) goto 42
C          call awrit7('found ip,jp=%i %i  ib,jb=%i %i ic,jc= %i %i'//
C     .      ' site %i',' ',80,6,ip,jp,ib,jb,ipxi(ii),ipxi(jj),icp)
          ipsite(ipxi(ii),ipxi(jj)) = icp
   40   continue
   30 continue

C --- Construct lower triangle of table ---
      do  50  i = 1, nclus3
      do  50  j = i+1, nclus3
        is = ipsite(i,j)
        if (is .gt. 0) ipsite(j,i) = iax(6,is)
   50 continue

C --- Check ---
      if (ipr .lt. 70) return
C      print '(''% rows'',i4,'' cols'',i4)', nclus3, nclus3
C      do  60  i = 1, nclus3
C   60 print '(35i4)', (ipsite(i,j), j=1, nclus3)
      do  70  ii = 1, nclus3
      do  70  jj = 1, nclus3
        i = iax3(ii)
        j = iax3(jj)
        ib = iax(2,kc-1+i)
        jb = iax(2,kc-1+j)
        is = ipsite(ii,jj)
        if (is .eq. 0) then
          print 356, kb,ib,jb
  356     format(' pairmc, kb=',i3,':  missing 2C pair for ib,jb=',2i4)
        elseif (iax(1,is) .ne. ib .or. iax(2,is) .ne. jb) then
          print *, 'problem with i,j,is=',i,j,ib,jb,is
          call rx('bug in pairmc')
        endif
   70 continue
      end
C Testing ...
C      subroutine fmain
C      implicit none
C      integer nbas
C      parameter (nbas=3)
C      double precision alat,plat(3,3),pos(3,nbas),rham(nbas),rtab(3,1)
C      integer ips(nbas)
C      integer mxcsiz,nttab,oiax,ontab,o,ortab
C      double precision w(10000)
C
C      data alat/1d0/ plat/10d0,0d0,0d0,0d0,10d0,0d0,0d0,0d0,10d0/,
C     .  pos/0d0,0d0,0d0,1.3804348d0,0d0,1.05924488d0,
C     .  -1.3804348d0,0d0,1.05924488d0/,
C     .  rham/2d0,3d0,0d0/, ips/2,1,1/
C      common /w/ w
C
C      call wkinit(10000)
C      call pshpr(100)
C
CC --- Obtain neighbor table for hamiltonian ---
C      call defrr(ontab, nbas+1)
C      call defrr(oiax, 200*nbas*6)
C      call defrr(o,     200*nbas*3)
C      call hpair(nbas,ips,alat,plat,pos,rham,nttab,w(ontab),
C     .  w(oiax),w(o),mxcsiz,1)
CC ... reallocate exact amount of space needed for iax,rtab
C      call rlse(oiax)
C      call defrr(oiax, nttab*6)
C      call defrr(ortab, nttab*3)
C      call dpcopy(w(o),w(ortab),1,nttab*3,1d0)
C
C      end

