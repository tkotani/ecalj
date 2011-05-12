      subroutine hgugbl(p1,p2,rsm1,rsm2,e1,nlm1,nlm2,ndim1,ndim2,slat,
     .s,ds)

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
      real(8):: rsm1 , rsm2 , p1(3) , p2(3) , e1
      type(s_lat)::slat

      double complex s(ndim1,ndim2),ds(ndim1,ndim2,3)
C ... Local parameters
          integer:: kmax , kdim , ilm2 , ilm1 
          integer,pointer :: iv_p_ojcg(:) =>NULL()

         integer,pointer :: iv_p_oidxcg(:) =>NULL()

        real(8),pointer :: rv_p_ocy(:) =>NULL()

       real(8),pointer :: rv_p_ocg(:) =>NULL()

      double precision q(3),e2,vol
C ... Heap
      data q /0d0,0d0,0d0/


       rv_p_ocg => slat%rv_p_ocg 

       iv_p_ojcg => slat%iv_p_ojcg 

       iv_p_oidxcg => slat%iv_p_oidxcg 

       rv_p_ocy => slat%rv_p_ocy 

      vol=slat%vol

      kmax = 0
      kdim = 0
      e2 = 0d0

          call hhigbl ( 0 , p1 , p2 , q , rsm1 , rsm2 , e1 , e2 , nlm1 
     .    , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_p_ocg , iv_p_oidxcg 
     .    , iv_p_ojcg , rv_p_ocy , slat , s , ds ) 





      do  ilm2 = 1, nlm2
        do  ilm1 = 1, nlm1
          s(ilm1,ilm2) = 2d0*s(ilm1,ilm2)
          ds(ilm1,ilm2,1) = 2d0*ds(ilm1,ilm2,1)
          ds(ilm1,ilm2,2) = 2d0*ds(ilm1,ilm2,2)
          ds(ilm1,ilm2,3) = 2d0*ds(ilm1,ilm2,3)
        enddo
      enddo

      end subroutine hgugbl


