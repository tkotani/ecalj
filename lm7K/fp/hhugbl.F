Cchp1 cmdvar= integer,iv_p_,ojcg
Cchp1 cmdvar= integer,iv_p_,oidxcg
Cchp1 cmdvar= real(8),rv_p_,ocy
Cchp1 cmdvar= real(8),rv_p_,ocg
      subroutine hhugbl(mode,p1,p2,rsm1,rsm2,e1,e2,nlm1,nlm2,ndim1,
     .ndim2,slat,wk,dwk,s,ds)

      use m_struc_def  !Cgetarg

C- Estatic energy integrals between Bloch Hankels, and gradients.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 rsm1,rsm2,e1,e2 are scalars
Ci         :1 rsm1,rsm2,e1,e2 are l-dependent
Ci   p1    :first center
Ci   p2    :second center
Ci   rsm1  :smoothing radius of Hankels at p1 (l-dependent)
Ci   rsm2  :smoothing radius of Hankels at p2 (l-dependent)
Ci   e1    :energy  of Hankels at p1 (l-dependent)
Ci   e2    :energy  of Hankels at p2 (l-dependent)
Ci   nlm1  :L-max for  Hankels at p1
Ci   nlm2  :L-max for  Hankels at p2
Ci   ndim1 :leading dimensions of s,ds
Ci   ndim2 :second dimensions of s,ds
Ci   slat  :struct containing information about the lattice
Ci   wk    :work space of same size as s
Ci   dwk   :work space of same size as ds
Cr Remarks
Cr   Gradient is wrt p1; use -ds for grad wrt p2.
Cu Updates
Cu   23 May 00 Made rsm1,e1,rsm2,e2 l-dependent
Cu   22 Apr 00 Adapted from nfp hhug_bl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nlm1,nlm2,ndim1,ndim2
      real(8):: rsm1(0:*) , rsm2(0:*) , e1(0:*) , e2(0:*) , p1(3) , 
     .p2(3)
      type(s_lat)::slat

      double complex s(ndim1,ndim2),ds(ndim1,ndim2,3),
     .wk(ndim1,ndim2),dwk(ndim1,ndim2,3)
C ... Local parameters
Cchp1       integer ocg,ojcg,oidxcg,ocy,kmax,kdim,i2,i1,lmxx,nlmx,l,ll
Cchp1        integer:: ojcg , oidxcg , ocy , kmax , kdim , i2 , i1 , lmxx 
Cchp1      . , nlmx , l , ll 
Cchp1         integer:: ojcg , oidxcg , kmax , kdim , i2 , i1 , lmxx , nlmx 
Cchp1      .  , l , ll 
Cchp1          integer:: ojcg , kmax , kdim , i2 , i1 , lmxx , nlmx , l , ll 
          integer:: kmax , kdim , i2 , i1 , lmxx , nlmx , l , ll 
          integer,pointer :: iv_p_ojcg(:)

         integer,pointer :: iv_p_oidxcg(:)

        real(8),pointer :: rv_p_ocy(:)

       real(8),pointer :: rv_p_ocg(:)

      parameter (lmxx=6,nlmx=(lmxx+1)**2)
      double precision add,pi,q(3),fpi,vol,gam1,gam2,xx1,xx2
      double precision zer(0:lmxx),bet1(nlmx),fac(nlmx)
C ... Heap
      integer w(1)
      common /w/ w
      data q /0d0,0d0,0d0/
      data zer /lmxx*0d0,0d0/

      pi = 4d0*datan(1d0)

Cchp1       ocg=slat%ocg
       rv_p_ocg => slat%rv_p_ocg 

Cchp1       ojcg=slat%ojcg
       iv_p_ojcg => slat%iv_p_ojcg 

Cchp1       oidxcg=slat%oidxcg
       iv_p_oidxcg => slat%iv_p_oidxcg 

Cchp1       ocy=slat%ocy
       rv_p_ocy => slat%rv_p_ocy 

      vol=slat%vol

      kmax = 0
      kdim = 0
      if (nlm1 .gt. nlmx) call rx('hhugbl: increase lmxx')

Cchp1       call hhigbl(mode,p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax,ndim1,
Cchp1      .ndim2,kdim,w(ocg),w(oidxcg),w(ojcg),w(ocy),slat,s,ds)
Cchp1        call hhigbl ( mode , p1 , p2 , q , rsm1 , rsm2 , e1 , e2 , nlm1 
Cchp1      . , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_p_ocg , w ( oidxcg 
Cchp1      . ) , w ( ojcg ) , w ( ocy ) , slat , s , ds ) 
Cchp1         call hhigbl ( mode , p1 , p2 , q , rsm1 , rsm2 , e1 , e2 , nlm1 
Cchp1      .  , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_p_ocg , w ( oidxcg 
Cchp1      .  ) , w ( ojcg ) , rv_p_ocy , slat , s , ds ) 
Cchp1          call hhigbl ( mode , p1 , p2 , q , rsm1 , rsm2 , e1 , e2 , nlm1 
Cchp1      .   , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_p_ocg , iv_p_oidxcg 
Cchp1      .   , w ( ojcg ) , rv_p_ocy , slat , s , ds ) 
          call hhigbl ( mode , p1 , p2 , q , rsm1 , rsm2 , e1 , e2 , nlm1 
     .    , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_p_ocg , iv_p_oidxcg 
     .    , iv_p_ojcg , rv_p_ocy , slat , s , ds ) 





Cchp1       call hhigbl(mode,p1,p2,q,rsm1,rsm2,e1,zer,nlm1,nlm2,kmax,ndim1,
Cchp1      .ndim2,kdim,w(ocg),w(oidxcg),w(ojcg),w(ocy),slat,wk,dwk)
Cchp1        call hhigbl ( mode , p1 , p2 , q , rsm1 , rsm2 , e1 , zer , nlm1 
Cchp1      . , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_p_ocg , w ( oidxcg 
Cchp1      . ) , w ( ojcg ) , w ( ocy ) , slat , wk , dwk ) 
Cchp1         call hhigbl ( mode , p1 , p2 , q , rsm1 , rsm2 , e1 , zer , nlm1 
Cchp1      .  , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_p_ocg , w ( oidxcg 
Cchp1      .  ) , w ( ojcg ) , rv_p_ocy , slat , wk , dwk ) 
Cchp1          call hhigbl ( mode , p1 , p2 , q , rsm1 , rsm2 , e1 , zer , nlm1 
Cchp1      .   , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_p_ocg , iv_p_oidxcg 
Cchp1      .   , w ( ojcg ) , rv_p_ocy , slat , wk , dwk ) 
          call hhigbl ( mode , p1 , p2 , q , rsm1 , rsm2 , e1 , zer , nlm1 
     .    , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_p_ocg , iv_p_oidxcg 
     .    , iv_p_ojcg , rv_p_ocy , slat , wk , dwk ) 





      do  i2 = 1, nlm2
        if (mode .eq. 0) then
          l = 0
        else
          l = ll(i2)
        endif
        bet1(i2) = dexp(e2(l)*rsm2(l)*rsm2(l)/4d0)
        fac(i2) = 8d0*pi/e2(l)
      enddo
      do  i2 = 1, nlm2
        do  i1 = 1, nlm1
          s(i1,i2)    = fac(i2)*(s(i1,i2)    - bet1(i2)*wk(i1,i2))
          ds(i1,i2,1) = fac(i2)*(ds(i1,i2,1) - bet1(i2)*dwk(i1,i2,1))
          ds(i1,i2,2) = fac(i2)*(ds(i1,i2,2) - bet1(i2)*dwk(i1,i2,2))
          ds(i1,i2,3) = fac(i2)*(ds(i1,i2,3) - bet1(i2)*dwk(i1,i2,3))
        enddo
      enddo

C ... Extra term for l1=l2=0
      fpi = 4d0*pi
      gam1 = 0.25d0*rsm1(0)*rsm1(0)
      gam2 = 0.25d0*rsm2(0)*rsm2(0)
      xx1 = fpi*dexp(gam1*e1(0))/(vol*e1(0))
      xx2 = fpi*dexp(gam2*e2(0))/(vol*e2(0))
      add = -2*xx1*xx2*vol/e2(0)
      s(1,1) = s(1,1) + add

      end subroutine hhugbl


