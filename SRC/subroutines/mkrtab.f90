subroutine mkrtab(mode,alat,plat,pos,iaxb,nttab,ctr,rtab)
  !- Make a table of connecting vectors for each cluster pair
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1s digit
  !i          0 make rtab in unpermuted order
  !i          1 make rtab in iax(7) order
  !i         :10s digit
  !i          1 put vector length into rtab
  !i         :100s digit
  !i          1 cluster index is always 1
  !i            (there is a single cluster; iaxb(1,ip) is not used)
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   plat  :primitive lattice vectors, in units of alat
  !i   pos   :basis vectors
  !i   iaxb  :neighbor table containing pair information (pairc.f)
  !i          May be a cluster or a hamiltonian neighbor table;
  !i          see description for ctr.
  !i   nttab :size of neighbor table
  !i   ctr   :table of cluster centers.   Must correspond to the index
  !i          iaxb(1,*).  Thus if the latter is a site index, ctr(ic)
  !i          must correspond to a site, i.e. ctr(ic)=pos(ic), while
  !i          if a cluster index, ctr(ic) must correspond to a cluster.
  !o Outputs
  !o   rtab  :connecting vectors rtab(1..3,ip) = pos(jb)-ctr(ic)
  !o          for pair jb=iaxb(2,ip) and ic=iaxb(1,ip).
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,niax
  parameter (niax=10)
  integer :: iaxb(niax,1),nttab
  double precision :: alat,plat(3,3),pos(3,1)
  double precision :: rtab(3,nttab),ctr(3,1)
  ! ... Local parameters
  logical :: lprm,llen,lonec
  integer :: ip,ipp,jb,ic,ix,i,offi
  double precision :: rtabi(3)

  lprm = mod(mode,10) .ne. 0
  llen = mod(mode/10,10) .ne. 0
  lonec= mod(mode/100,10) .ne. 0
  i = 0
  offi = 0
  do  40  ip = 1, nttab
     ic = iaxb(1,ip)
     if (lonec) ic = 1
     ipp = ip
     if (lprm) then
        !         New cluster
        if (ic /= i) offi = ip-1
        ipp = offi + iaxb(7,ip)
        if (ipp > nttab) call rx('mkrtab: bad iax table')
     endif
     i = ic
     jb = iaxb(2,ipp)
     do  42  ix = 1, 3
        rtabi(ix) = alat * &
             ( pos(ix,jb) - ctr(ix,ic) &
             + plat(ix,1)*iaxb(3,ipp) &
             + plat(ix,2)*iaxb(4,ipp) &
             + plat(ix,3)*iaxb(5,ipp))
        if ( .NOT. llen) rtab(ix,ip) = rtabi(ix)
42   enddo
     if (llen) then
        rtab(ip,1) = dsqrt(rtabi(1)**2+rtabi(2)**2+rtabi(3)**2)
     endif
40 enddo
  !      if (llen) call prmx('rtab',rtab,nttab,nttab,1)
  !      if (.not. llen) call prmx('rtab',rtab,3,3,nttab)
end subroutine mkrtab

