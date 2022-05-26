      subroutine shoist(istab,nbas,ag,g,ng)
      use m_lgunit,only:stdo
C- Show istab
C     implicit none
      integer ng,nbas
      double precision g(3,3,ng),ag(3,ng)
      integer istab(nbas,ng)
      integer i,ig
      write(stdo,*)'  ib  istab ...'
      do  30  i = 1, nbas
        write(stdo,333) i, (istab(i,ig), ig=1,ng)
  333   format(i4,':',48i3)
   30 continue
      end

