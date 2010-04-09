Cgetarg...info...           structure ['hgugbl', 'lat', 'slat']
Cgetarg...info...           use_to_add ['hgugbl', 'm_struc_def']
      subroutine hgugbl(p1,p2,rsm1,rsm2,e1,nlm1,nlm2,ndim1,ndim2,slat,
     .  s,ds)
       
       use m_struc_def  !Cgetarg

C- Estatic energy integrals between Bloch Hankels and gaussians, and grad
C ----------------------------------------------------------------------
Ci Inputs
Ci   slat  :struct containing information about the lattice
Ci   p1    :first center
Ci   p2    :second center
Ci   rsm1  :smoothing radius of Hankels at p1
Ci   rsm2  :smoothing radius of gaussians at p2
Ci   e1    :energy  of Hankels at p1
Ci   nlm1  :L-max for  Hankels at p1
Ci   nlm2  :L-max for  gaussians at p2
Ci   ndim1 :leading dimensions of s,ds
Ci   ndim2 :second dimensions of s,ds
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   s     :integrals between Bloch Hankels and gaussians
Co   ds    :gradient of s; see Remarks
Cr Remarks
Cr   Gradient is wrt p1; use -ds for grad wrt p2.
Cu Updates
Cu   22 Apr 00 Adapted from nfp hgug_bl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlm1,nlm2,ndim1,ndim2
Cgetarg       double precision slat(1),rsm1,rsm2,p1(3),p2(3),e1
       real(8):: rsm1 , rsm2 , p1(3) , p2(3) , e1 
       type(s_lat)::slat

      double complex s(ndim1,ndim2),ds(ndim1,ndim2,3)
C ... Local parameters
      integer ocg,ojcg,oidxcg,ocy,kmax,kdim,ilm2,ilm1
      double precision q(3),e2,vol
C ... Heap
      integer w(1)
      common /w/ w
      data q /0d0,0d0,0d0/

Cgetarg...info...           integer :: ocg , integer(8) :: slat%ocg 1
Cgetarg...info...           integer :: ojcg , integer(8) :: slat%ojcg 1
Cgetarg...info...           integer :: oidxcg , integer(8) :: slat%oidxcg 1
Cgetarg...info...           integer :: ocy , integer(8) :: slat%ocy 1
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg       call upack('lat ocg ojcg oidxcg ocy vol',slat,ocg,ojcg,oidxcg,ocy,
Cgetarg      .  vol)
       
       ocg=slat%ocg
       ojcg=slat%ojcg
       oidxcg=slat%oidxcg
       ocy=slat%ocy
       vol=slat%vol

      kmax = 0
      kdim = 0
      e2 = 0d0

      call hhigbl(0,p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax,ndim1,ndim2,
     .   kdim,w(ocg),w(oidxcg),w(ojcg),w(ocy),slat,s,ds)

      do  ilm2 = 1, nlm2
        do  ilm1 = 1, nlm1
          s(ilm1,ilm2) = 2d0*s(ilm1,ilm2)
          ds(ilm1,ilm2,1) = 2d0*ds(ilm1,ilm2,1)
          ds(ilm1,ilm2,2) = 2d0*ds(ilm1,ilm2,2)
          ds(ilm1,ilm2,3) = 2d0*ds(ilm1,ilm2,3)
        enddo
      enddo

Cgetarg       end
       end subroutine hgugbl 


