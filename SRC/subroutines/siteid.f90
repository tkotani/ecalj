subroutine siteid(iax,ncl,ntab,plat,bas,pos,iwk,nid)
  !- Assign a unique id for every different site in the neighbor table
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   iax   :array of parameters containing info about each pair
  !i   ncl   :number of clusters for which neighbor table is made
  !i   ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairs.f)
  !i   plat  :primitive lattice vectors, in units of alat (input)
  !i   bas   :basis vectors (input)
  !i   pos   :work array, of dimension 3*nttab, with nttab=ntab(ncl+1)
  !i   iwk   :work array, of dimension nttab, with nttab=ntab(ncl+1)
  !o   iax   :iax(7) is needed, which orders cluster by increasing
  !r          (x,y,z) from the origin (ppair4)
  !o Outputs
  !o   iax   :iax(10) is assigned a unique site ID for each entry
  !o   nid   :number of unique sites
  !r Remarks
  !r   A unique ID (in iax(10)) is assigned to each entry in the
  !r   neighbor table.  Sites with the same ID are identical sites.
  !r   The ordering is the same as in ppair4.
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: niax,ncl,ntab(ncl+1),iwk(*),nid
  double precision :: bas(3,ncl),pos(3,*),plat(3,3)
  parameter (niax=10)
  integer :: iax(niax,1)
  ! Local variables
  integer :: ic,jb,ix,nttab,id,ip,imin,iminp,i1,icp
  double precision :: tol,dmatch(3)
  !     tolerance should be same as ppair4
  parameter (tol=1d-5)

  ! --- Table of absolute positions ---
  nttab = ntab(ncl+1)
  do  10  ic = 1, ncl
     do  12  ip = ntab(ic)+1, ntab(ic+1)
        !       ib = iax(i,ip)
        jb = iax(2,ip)
        do  14  ix = 1, 3
           pos(ix,ip) = bas(ix,jb) + plat(ix,1)*iax(3,ip) + &
                plat(ix,2)*iax(4,ip) + plat(ix,3)*iax(5,ip)
14      enddo
12   enddo
10 enddo

  ! --- Assign a unique id ---
  !     For cluster ic, iwk(ic)=index to next site in ic needing id
  do  30  ic = 1, ncl
     iwk(ic) = 1
30 enddo
  dmatch(1) = 9d9
  dmatch(2) = 9d9
  dmatch(3) = 9d9
  id = 0
  ! --- Loop through all neighbors in the entire cluster table ---
  do  20  ip = 1, nttab
     !   ... Get first site not marked in cluster table
     imin = 0
33   imin = imin+1
     !       Skip this cluster if every member already assigned id
     if (iwk(imin) > ntab(imin+1)-ntab(imin)) goto 33
     !   ... iminp points to next element in sorted list of pairs at imin
     iminp = ntab(imin) + iax(7,ntab(imin)+iwk(imin))
     i1 = imin+1
     do  34  ic = i1, ncl
        if (iwk(ic) > ntab(ic+1)-ntab(ic)) goto 34
        icp = ntab(ic) + iax(7,ntab(ic)+iwk(ic))
        !          print 333, pos(1,iminp),pos(2,iminp),pos(3,iminp)
        !          print 333, pos(1,icp),pos(2,icp),pos(3,icp)
        !  333     format(3f12.6)
        !    ...  Exclude ic if ic(1)>imin(1)
        if (pos(1,iminp)+tol < pos(1,icp)) goto 34
        !    ...  ic becomes new imin if ic(1)<imin(1)
        if (abs(pos(1,iminp)-pos(1,icp)) > tol) goto 35
        !    ...  Exclude ic if ic(2)>imin(2)
        if (pos(2,iminp)+tol < pos(2,icp)) goto 34
        !    ...  ic becomes new imin if ic(2)<imin(2)
        if (abs(pos(2,iminp)-pos(2,icp)) > tol) goto 35
        !    ...  Exclude ic if ic(3)>imin(3)
        if (pos(3,iminp)+tol < pos(3,icp)) goto 34
        !    ...  ic becomes new imin if ic(3)<imin(3)
        if (abs(pos(3,iminp)-pos(3,icp)) <= tol) goto 34
35      continue
        imin = ic
        iminp = icp
34   enddo
     !    .. imin holds among (1..ib2) next lowest element
     iwk(imin) = iwk(imin) + 1
     !       call awrit2('iwk %n:1i',' ',180,6,ncl,iwk(1))

     !   ... If no match with previous, increment new id
     if (abs(pos(1,iminp)-dmatch(1)) > tol .OR. &
          abs(pos(2,iminp)-dmatch(2)) > tol .OR. &
          abs(pos(3,iminp)-dmatch(3)) > tol) then
        id = id+1
        dmatch(1) = pos(1,iminp)
        dmatch(2) = pos(2,iminp)
        dmatch(3) = pos(3,iminp)
     endif
     iax(10,iminp) = id
20 enddo
  nid = id

  !    .. Debugging check
  do  38  ic = 1, ncl
     if (iwk(ic)-1 /= ntab(ic+1)-ntab(ic)) call rx('bug in siteid')
38 enddo

  ! ... Debugging printout
  !      call dvheap(3,nttab,pos,iwk(1),tol,1)
  !      do  50  ip = 1, nttab
  !        jc = iwk(ip)
  !        print 333, ip, jc, iax(10,jc), pos(1,jc),pos(2,jc),pos(3,jc)
  !  333   format(3i6,3f12.6)
  !   50 continue
  !      call rx('done')
end subroutine siteid
subroutine nsitsh(mode,ia,ib,iax,ntab,nlst,lsta,lstb)
  !- Count the number of sites two clusters share in common
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 count the sites only
  !i          1 also make a list of the sites
  !i   ia,ib :pair of sites for which to seek common elements
  !i   nds   :leading dimension of sid
  !i   iax   :neighbor table containing pair information (pairc.f)
  !i   ntab  :ntab(ib)=# pairs in iax table preceding ib (pairc.f)
  !o Outputs
  !o   nlst  :number of sites in common
  !o   lsta  :list of
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  !     Passed parameters
  integer :: mode,ia,ib,ntab(1),niax,lsta(1),lstb(1),nlst
  parameter (niax=10)
  integer :: iax(niax,1)
  !     Local variables
  integer :: ic,ica,icb,n,na,nb,ipa,ipb,sida,sidb

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
     if (ica > na .OR. icb > nb) return

     !   ... A match; increment nlst
     if (sida == sidb) then
        print *, sida
        nlst = nlst+1
        ica = ica+1
        icb = icb+1
        ipa = iax(7,ica)
        ipb = iax(7,icb)
        sida = iax(10,ipa+ntab(ia))
        sidb = iax(10,ipb+ntab(ib))
        if (mode /= 0) then
           lsta(nlst) = ica
           lstb(nlst) = icb
        endif
     elseif (sida > sidb) then
        icb = icb+1
        ipb = iax(7,icb)
        sidb = iax(10,ipb+ntab(ib))
     else
        ica = ica+1
        ipa = iax(7,ica)
        sida = iax(10,ipa+ntab(ia))
     endif
10 enddo
end subroutine nsitsh

