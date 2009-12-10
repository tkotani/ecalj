      subroutine shoist(istab,nbas,ag,g,ng)
C- Show istab
C     implicit none
      integer ng,nbas
      double precision g(3,3,ng),ag(3,ng)
      integer istab(nbas,ng)
      integer i,ig,lgunit,stdo

C      do  20  ig = 1, ng
C      do  20  ig = 2,2
C      print 331, (ag(i,ig), i=1,3)
C   20 print 331, ((g(i,k,ig), k=1,3), i=1,3)
C  331 format(3(3f12.7/))
C
      call info0(0,0,0,'  ib  istab ...')
      stdo = lgunit(1)

      do  30  i = 1, nbas
        write(stdo,333) i, (istab(i,ig), ig=1,ng)
  333   format(i4,':',48i3)
   30 continue
      end

