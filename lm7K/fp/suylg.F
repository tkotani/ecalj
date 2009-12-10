      subroutine suylg(ltop,alat,ng,gv,g,g2,yl)
C- Set up vectors g, g2, yl from list of vectors gv
C ----------------------------------------------------------------------
Ci Inputs
Ci   ltop  :l-cutoff for YL
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Co Outputs
Co   g     :gv scaled by (2 pi / alat)
Co   g2    :square of g
Co   yl    :YL(g)
Cu Updates
Cu   30 May 00 adapted from nfp su_ylg
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ltop,ng
      double precision alat,gv(ng,3),g(ng,3),yl(ng,1),g2(ng)
C ... Local parameters
      integer i
      double precision pi,tpiba

C ... Make (2*pi/alat)*gv in g
      pi = 4d0*datan(1d0)
      tpiba = 2d0*pi/alat
      do  i = 1, ng
        g(i,1) = tpiba*gv(i,1)
        g(i,2) = tpiba*gv(i,2)
        g(i,3) = tpiba*gv(i,3)
      enddo
      
C ... Make the yl's and g2
      call ropyln(ng,g(1,1),g(1,2),g(1,3),ltop,ng,yl,g2)

      end

