      integer function streqv(ntab,ldot,nla,alpha,adot,ib,ib0,plat,pos,
     .  iax,nitab)
C- Determines wether a local site has a structure equivalent to another
C     implicit none
C Passed parameters
      logical ldot
      integer ntab(1),niax,ib,ib0,nitab,nla
      parameter (niax=10)
      integer iax(niax,45)
      double precision alpha(nla,1),adot(nla,1),pos(3,18),plat(3,3)
C Local parameters
      integer ia,iclusb,nclus,iclusa,i,j,k,ic,jcb,jca,ix,
     .  ia1,ia2,jb,ja
      double precision tolb,tola,xx,drb(3),dra(3),dx
      parameter (tolb=1d-6, tola=1d-8)
      dx(ia1,ia2,i,j,k) = pos(ix,ia2) - pos(ix,ia1) +
     .                    plat(ix,1)*i + plat(ix,2)*j + plat(ix,3)*k

      streqv = -1
      nclus = ntab(ib+1)-ntab(ib)
      iclusb = ntab(ib)

C --- Loop over all prior clusters ---
      do  30  ia  = ib0+1, ib-1

C   ... Match in cluster size? ...
        if (nclus .ne. ntab(ia+1)-ntab(ia)) goto 30
        iclusa = ntab(ia)

C   --- Loop over all pairs in iax(7) order ---
        do  40  ic = 1, nclus

          jcb = iclusb + iax(7,ic+iclusb)
          jca = iclusa + iax(7,ic+iclusa)

C     ... Connecting vectors within tol?
          jb = iax(2,jcb)
          ja = iax(2,jca)
          do  42  ix = 1, 3
            drb(ix) = dx(ib,jb,iax(3,jcb),iax(4,jcb),iax(5,jcb))
            dra(ix) = dx(ia,ja,iax(3,jca),iax(4,jca),iax(5,jca))
            if (abs(drb(ix)-dra(ix)) .gt. tolb) goto 30
   42     continue
C          print 333, ib,jb,jcb,drb
C          print 333, ia,ja,jca,dra
C  333     format(3i4,3f11.6,2f9.4,3x,3i3,i5,i4)
            
C     ... alpha, alpha-dot within tol?
          do  44  ix = 1, nla
            xx = alpha(ix,jb) - alpha(ix,ja)
            if (abs(xx) .gt. tola) goto 30
            if (ldot) then
              xx = adot(ix,jb) - adot(ix,ja)
              if (abs(xx) .gt. tola) goto 30
            endif
   44     continue

C     ... In case ia equiv to ib, set iax(8) now
          iax(8,jcb) = jca
          if (iax(8,jca) .ne. 0) iax(8,jcb) = iax(8,jca)
   40   continue

C   --- ia passes the equivalence test: exeunt ---
        streqv = ia
        return
   30 continue

C --- No preceding strx equivalent to this cluster; set iax(8)=0 ---
      do  50  ic = ntab(ib)+1, ntab(ib+1)
      nitab = nitab+1
   50 iax(8,ic) = nitab

      end


