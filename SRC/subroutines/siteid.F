      subroutine siteid(iax,ncl,ntab,plat,bas,pos,iwk,nid)
C- Assign a unique id for every different site in the neighbor table
C ----------------------------------------------------------------------
Ci Inputs
Ci   iax   :array of parameters containing info about each pair
Ci   ncl   :number of clusters for which neighbor table is made
Ci   ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairs.f)
Ci   plat  :primitive lattice vectors, in units of alat (input)
Ci   bas   :basis vectors (input)
Ci   pos   :work array, of dimension 3*nttab, with nttab=ntab(ncl+1)
Ci   iwk   :work array, of dimension nttab, with nttab=ntab(ncl+1)
Co   iax   :iax(7) is needed, which orders cluster by increasing
Cr          (x,y,z) from the origin (ppair4)
Co Outputs
Co   iax   :iax(10) is assigned a unique site ID for each entry
Co   nid   :number of unique sites
Cr Remarks
Cr   A unique ID (in iax(10)) is assigned to each entry in the
Cr   neighbor table.  Sites with the same ID are identical sites.
Cr   The ordering is the same as in ppair4.
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer niax,ncl,ntab(ncl+1),iwk(*),nid
      double precision bas(3,ncl),pos(3,*),plat(3,3)
      parameter (niax=10)
      integer iax(niax,1)
C Local variables
      integer ic,jb,ix,nttab,id,ip,imin,iminp,i1,icp
      double precision tol,dmatch(3)
C     tolerance should be same as ppair4
      parameter (tol=1d-5)

C --- Table of absolute positions ---
      nttab = ntab(ncl+1)
      do  10  ic = 1, ncl
        do  12  ip = ntab(ic)+1, ntab(ic+1)
C       ib = iax(i,ip)
          jb = iax(2,ip)
          do  14  ix = 1, 3
            pos(ix,ip) = bas(ix,jb) + plat(ix,1)*iax(3,ip) +
     .      plat(ix,2)*iax(4,ip) + plat(ix,3)*iax(5,ip)
   14     continue
   12   continue
   10 continue

C --- Assign a unique id ---
C     For cluster ic, iwk(ic)=index to next site in ic needing id
      do  30  ic = 1, ncl
        iwk(ic) = 1
   30 continue
      dmatch(1) = 9d9
      dmatch(2) = 9d9
      dmatch(3) = 9d9
      id = 0
C --- Loop through all neighbors in the entire cluster table ---
      do  20  ip = 1, nttab
C   ... Get first site not marked in cluster table
        imin = 0
   33   imin = imin+1
C       Skip this cluster if every member already assigned id
        if (iwk(imin) .gt. ntab(imin+1)-ntab(imin)) goto 33
C   ... iminp points to next element in sorted list of pairs at imin
        iminp = ntab(imin) + iax(7,ntab(imin)+iwk(imin))
        i1 = imin+1
        do  34  ic = i1, ncl
          if (iwk(ic) .gt. ntab(ic+1)-ntab(ic)) goto 34
          icp = ntab(ic) + iax(7,ntab(ic)+iwk(ic))
C          print 333, pos(1,iminp),pos(2,iminp),pos(3,iminp)
C          print 333, pos(1,icp),pos(2,icp),pos(3,icp)
C  333     format(3f12.6)
C    ...  Exclude ic if ic(1)>imin(1)
          if (pos(1,iminp)+tol .lt. pos(1,icp)) goto 34
C    ...  ic becomes new imin if ic(1)<imin(1)
          if (abs(pos(1,iminp)-pos(1,icp)) .gt. tol) goto 35
C    ...  Exclude ic if ic(2)>imin(2)
          if (pos(2,iminp)+tol .lt. pos(2,icp)) goto 34
C    ...  ic becomes new imin if ic(2)<imin(2)
          if (abs(pos(2,iminp)-pos(2,icp)) .gt. tol) goto 35
C    ...  Exclude ic if ic(3)>imin(3)
          if (pos(3,iminp)+tol .lt. pos(3,icp)) goto 34
C    ...  ic becomes new imin if ic(3)<imin(3)
          if (abs(pos(3,iminp)-pos(3,icp)) .le. tol) goto 34
   35     continue
          imin = ic
          iminp = icp
   34   continue
C    .. imin holds among (1..ib2) next lowest element
        iwk(imin) = iwk(imin) + 1
C       call awrit2('iwk %n:1i',' ',180,6,ncl,iwk(1))

C   ... If no match with previous, increment new id
        if (abs(pos(1,iminp)-dmatch(1)) .gt. tol .or.
     .  abs(pos(2,iminp)-dmatch(2)) .gt. tol .or.
     .  abs(pos(3,iminp)-dmatch(3)) .gt. tol) then
          id = id+1
          dmatch(1) = pos(1,iminp)
          dmatch(2) = pos(2,iminp)
          dmatch(3) = pos(3,iminp)
        endif
        iax(10,iminp) = id
   20 continue
      nid = id

C    .. Debugging check
      do  38  ic = 1, ncl
        if (iwk(ic)-1 .ne. ntab(ic+1)-ntab(ic)) call rx('bug in siteid')
   38 continue

C ... Debugging printout
C      call dvheap(3,nttab,pos,iwk(1),tol,1)
C      do  50  ip = 1, nttab
C        jc = iwk(ip)
C        print 333, ip, jc, iax(10,jc), pos(1,jc),pos(2,jc),pos(3,jc)
C  333   format(3i6,3f12.6)
C   50 continue
C      call rx('done')
      end
      subroutine nsitsh(mode,ia,ib,iax,ntab,nlst,lsta,lstb)
C- Count the number of sites two clusters share in common
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 count the sites only
Ci          1 also make a list of the sites
Ci   ia,ib :pair of sites for which to seek common elements
Ci   nds   :leading dimension of sid
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ntab  :ntab(ib)=# pairs in iax table preceding ib (pairc.f)
Co Outputs
Co   nlst  :number of sites in common
Co   lsta  :list of
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C     Passed parameters
      integer mode,ia,ib,ntab(1),niax,lsta(1),lstb(1),nlst
      parameter (niax=10)
      integer iax(niax,1)
C     Local variables
      integer ic,ica,icb,n,na,nb,ipa,ipb,sida,sidb

      ica = ntab(ia)+1
      icb = ntab(ib)+1
      ipa = iax(7,ica)
      ipb = iax(7,icb)
      sida = iax(10,ipa+ntab(ia))
      sidb = iax(10,ipb+ntab(ib))
      na = ntab(ia+1)
      nb = ntab(ib+1)
      n = na-ica + nb-icb + 2
      nlst = 0
      do  10  ic = 1, n
        if (ica .gt. na .or. icb .gt. nb) return

C   ... A match; increment nlst
        if (sida .eq. sidb) then
          print *, sida
          nlst = nlst+1
          ica = ica+1
          icb = icb+1
          ipa = iax(7,ica)
          ipb = iax(7,icb)
          sida = iax(10,ipa+ntab(ia))
          sidb = iax(10,ipb+ntab(ib))
          if (mode .ne. 0) then
            lsta(nlst) = ica
            lstb(nlst) = icb
          endif
        elseif (sida .gt. sidb) then
          icb = icb+1
          ipb = iax(7,icb)
          sidb = iax(10,ipb+ntab(ib))
        else
          ica = ica+1
          ipa = iax(7,ica)
          sida = iax(10,ipa+ntab(ia))
        endif
   10 continue
      end

