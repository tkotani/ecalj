      subroutine mkrtab(mode,alat,plat,pos,iaxb,nttab,ctr,rtab)
C- Make a table of connecting vectors for each cluster pair
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci          0 make rtab in unpermuted order
Ci          1 make rtab in iax(7) order
Ci         :10s digit
Ci          1 put vector length into rtab
Ci         :100s digit
Ci          1 cluster index is always 1
Ci            (there is a single cluster; iaxb(1,ip) is not used)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   iaxb  :neighbor table containing pair information (pairc.f)
Ci          May be a cluster or a hamiltonian neighbor table;
Ci          see description for ctr.
Ci   nttab :size of neighbor table
Ci   ctr   :table of cluster centers.   Must correspond to the index
Ci          iaxb(1,*).  Thus if the latter is a site index, ctr(ic)
Ci          must correspond to a site, i.e. ctr(ic)=pos(ic), while
Ci          if a cluster index, ctr(ic) must correspond to a cluster.
Co Outputs
Co   rtab  :connecting vectors rtab(1..3,ip) = pos(jb)-ctr(ic)
Co          for pair jb=iaxb(2,ip) and ic=iaxb(1,ip).
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,niax
      parameter (niax=10)
      integer iaxb(niax,1),nttab
      double precision alat,plat(3,3),pos(3,1)
      double precision rtab(3,nttab),ctr(3,1)
C ... Local parameters
      logical lprm,llen,lonec
      integer ip,ipp,jb,ic,ix,i,offi
      double precision rtabi(3)

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
C         New cluster
          if (ic .ne. i) offi = ip-1
          ipp = offi + iaxb(7,ip)
          if (ipp .gt. nttab) call rx('mkrtab: bad iax table')
        endif
        i = ic
        jb = iaxb(2,ipp)
        do  42  ix = 1, 3
          rtabi(ix) = alat *
     .    ( pos(ix,jb) - ctr(ix,ic)
     .    + plat(ix,1)*iaxb(3,ipp)
     .    + plat(ix,2)*iaxb(4,ipp)
     .    + plat(ix,3)*iaxb(5,ipp))
          if (.not. llen) rtab(ix,ip) = rtabi(ix)
   42   continue
        if (llen) then
          rtab(ip,1) = dsqrt(rtabi(1)**2+rtabi(2)**2+rtabi(3)**2)
        endif
   40 continue
C      if (llen) call prmx('rtab',rtab,nttab,nttab,1)
C      if (.not. llen) call prmx('rtab',rtab,3,3,nttab)
      end

