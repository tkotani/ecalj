Cchp1 cmdvar= integer,iv_p_,ojcg
Cchp1 cmdvar= integer,iv_p_,oidxcg
Cchp1 cmdvar= real(8),rv_p_,ocy
Cchp1 cmdvar= real(8),rv_p_,ocg
      subroutine ggugbl(p1,p2,rsm1,rsm2,nlm1,nlm2,ndim1,ndim2,slat,s,ds)

      use m_struc_def  !Cgetarg

C- Estatic energy integrals between Bloch gaussians, and gradients.
C ----------------------------------------------------------------------
Ci Inputs
Ci   slat  :struct containing information about the lattice
Ci   p1    :first center
Ci   p2    :second center
Ci   rsm1  :smoothing radius of Gaussians at p1
Ci   rsm2  :smoothing radius of Gaussians at p2
Ci   e1    :energy  of Gaussians at p1
Ci   e2    :energy  of Gaussians at p2
Ci   nlm1  :L-max for  Gaussians at p1
Ci   nlm2  :L-max for  Gaussians at p2
Ci   ndim1 :leading dimensions of s,ds
Ci   ndim2 :second dimensions of s,ds
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   s     :integrals between Bloch Gaussians
Co   ds    :gradient of s; see Remarks
Cr Remarks
Cr   Gradient is wrt p1; use -ds for grad wrt p2.
Cu Updates
Cu   22 Apr 00 Adapted from nfp ggug_bl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlm1,nlm2,ndim1,ndim2
      real(8):: rsm1 , rsm2 , p1(3) , p2(3)
      type(s_lat)::slat

      double complex s(ndim1,ndim2),ds(ndim1,ndim2,3)
C ... Local parameters
Cchp1       integer ocg,ojcg,oidxcg,ocy,kmax,kdim,ilm2,ilm1
Cchp1        integer:: ojcg , oidxcg , ocy , kmax , kdim , ilm2 , ilm1 
Cchp1         integer:: ojcg , oidxcg , kmax , kdim , ilm2 , ilm1 
Cchp1          integer:: ojcg , kmax , kdim , ilm2 , ilm1 
          integer:: kmax , kdim , ilm2 , ilm1 
          integer,pointer :: iv_p_ojcg(:)

         integer,pointer :: iv_p_oidxcg(:)

        real(8),pointer :: rv_p_ocy(:)

       real(8),pointer :: rv_p_ocg(:)

C ... Heap
      integer w(1)
      common /w/ w


Cchp1       ocg=slat%ocg
       rv_p_ocg => slat%rv_p_ocg 

Cchp1       ojcg=slat%ojcg
       iv_p_ojcg => slat%iv_p_ojcg 

Cchp1       oidxcg=slat%oidxcg
       iv_p_oidxcg => slat%iv_p_oidxcg 

Cchp1       ocy=slat%ocy
       rv_p_ocy => slat%rv_p_ocy 


      kmax = 0
      kdim = 0
Cchp1       call gfigbl(p1,p2,rsm1,rsm2,nlm1,nlm2,kmax,ndim1,ndim2,kdim,
Cchp1      .w(ocg),w(oidxcg),w(ojcg),w(ocy),slat,s,ds)
Cchp1        call gfigbl ( p1 , p2 , rsm1 , rsm2 , nlm1 , nlm2 , kmax , ndim1 
Cchp1      . , ndim2 , kdim , rv_p_ocg , w ( oidxcg ) , w ( ojcg ) , w ( ocy 
Cchp1      . ) , slat , s , ds ) 
Cchp1         call gfigbl ( p1 , p2 , rsm1 , rsm2 , nlm1 , nlm2 , kmax , ndim1 
Cchp1      .  , ndim2 , kdim , rv_p_ocg , w ( oidxcg ) , w ( ojcg ) , rv_p_ocy 
Cchp1      .  , slat , s , ds ) 
Cchp1          call gfigbl ( p1 , p2 , rsm1 , rsm2 , nlm1 , nlm2 , kmax , ndim1 
Cchp1      .   , ndim2 , kdim , rv_p_ocg , iv_p_oidxcg , w ( ojcg ) , rv_p_ocy 
Cchp1      .   , slat , s , ds ) 
          call gfigbl ( p1 , p2 , rsm1 , rsm2 , nlm1 , nlm2 , kmax , ndim1 
     .    , ndim2 , kdim , rv_p_ocg , iv_p_oidxcg , iv_p_ojcg , rv_p_ocy 
     .    , slat , s , ds ) 





      do  ilm2 = 1, nlm2
        do  ilm1 = 1, nlm1
          s(ilm1,ilm2) = 2d0*s(ilm1,ilm2)
          ds(ilm1,ilm2,1) = 2d0*ds(ilm1,ilm2,1)
          ds(ilm1,ilm2,2) = 2d0*ds(ilm1,ilm2,2)
          ds(ilm1,ilm2,3) = 2d0*ds(ilm1,ilm2,3)
        enddo
      enddo

      end subroutine ggugbl


