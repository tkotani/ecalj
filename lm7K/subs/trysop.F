      subroutine trysop(ib,jb,g,nbas,ipclas,tau,rb,qb,lok,trans)
C- Tests whether crystal as seen from site ib same as that seen from jb
C  after point operation g.  If so, returns trans and lok=1.
C     implicit none
      integer ipclas(1)
      integer lok,ib,jb,nbas
      double precision rb(9),qb(9),g(3,3),tau(3,1),d(3),dif(3),trans(3)
      integer kb,m,k,kbp,kk,iprint,i1mach,iclas1,iclas2

      lok = 0
      do 10 kb = 1,nbas
        do 2 m = 1,3
          d(m) = 0
          do 2 k = 1,3
            d(m) = d(m)+g(m,k)*(tau(k,kb)-tau(k,ib))
    2   continue
        kbp = 0
        do 11 kk = 1,nbas
          do 12 m = 1,3
            dif(m) = d(m)-(tau(m,kk)-tau(m,jb))
   12     continue
          call shorbz(dif,dif,rb,qb)
          if(dif(1)**2+dif(2)**2+dif(3)**2.lt.1.d-10) then
            kbp = kk
            goto 13
          endif
   11   continue
        if (iprint() .ge. 60) then
          call awrit3(' trysop ib=%i:  no site analogous to'//
     .    ' site %i from site %i',' ',80,i1mach(2),ib,kb,jb)
          if (iprint() .ge. 80) then
            call awrit1(' missing atom at: %3:1,6;6d',' ',80,i1mach(2),dif)
            do 32 m = 1,3
              dif(m) = tau(m,kb)-tau(m,ib)
   32       continue
            call awrit1('    vector kb-ib: %3:1,6;6d',' ',80,i1mach(2),dif)
            call awrit1('  after rotation: %3:1,6;6d',' ',80,i1mach(2),d)
          endif
        endif
        return
   13   iclas1 = ipclas(kb)
        iclas2 = ipclas(kbp)
        if (iprint().ge.60 .and.iclas1.ne.iclas2 .or. iprint().ge.70) then
          call awrit6(' trysop ib=%i jb=%i:  site %i class %i'//
     .    ' maps to site %i class %i',' ',80,i1mach(2),ib,jb,
     .    kb,iclas1,kbp,iclas2)
          if (iprint() .ge. 80) then
            do 42 m = 1,3
              dif(m) = tau(m,kb)-tau(m,ib)
   42       continue
            call awrit1('    vector kb-ib: %3:1,6;6d',' ',80,i1mach(2),dif)
            call awrit1('  after rotation: %3:1,6;6d',' ',80,i1mach(2),d)
          endif
        endif

        if (iclas1.ne.iclas2) return
   10 continue

      lok = 1
      do 3 m = 1,3
        trans(m) = tau(m,jb)
        do 3 k = 1,3
          trans(m) = trans(m)-g(m,k)*tau(k,ib)
    3 continue
      call shorbz(trans,trans,rb,qb)
      end

