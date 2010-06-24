Cchp1 cmdvar= real(8),rv_p_,ocy
Cchp1 cmdvar= real(8),rv_p_,ocg
      subroutine setcg(slat,lmxcg,lmxcy)

      use m_struc_def  !Cgetarg

C- Allocate space for, and make Clebsch-Gordan coeffs
C     implicit none
      integer lmxcg,lmxcy
      type(s_lat)::slat

Cchp1       integer lnjcg,lnxcg,nlm,ocy,ocg,ojcg,oidxcg
Cchp1        integer:: lnjcg , lnxcg , nlm , ocy , ojcg , oidxcg 
        integer:: lnjcg , lnxcg , nlm , ojcg , oidxcg 
        real(8),pointer :: rv_p_ocy(:)

       real(8),pointer :: rv_p_ocg(:)

      integer w(1)
      common /w/ w

C ... Choose dimensions for arrays
      if (lmxcg .le. 6) then
        lnjcg = 6500
        lnxcg = 1300
      else if (lmxcg .le. 8) then
        lnjcg = 22700
        lnxcg = 3400
      else if (lmxcg .le. 10) then
        lnjcg = 62200
        lnxcg = 7400
      else
        call rxi('setcg: cannot handle lmxcg=',lmxcg)
      endif

C ... Allocate and occupy the arrays
      nlm=(lmxcy+1)**2
Cchp1       call defrr (ocy,      nlm)
       allocate(rv_p_ocy(abs(nlm)))
       if (nlm<0) rv_p_ocy(:)=0.0d0

Cchp1       call defrr (ocg,      lnjcg)
       allocate(rv_p_ocg(abs(lnjcg)))
       if (lnjcg<0) rv_p_ocg(:)=0.0d0

      call defi  (ojcg,     lnjcg)
      call defi  (oidxcg,   lnxcg)

Cchp1       call sylmnc(w(ocy),lmxcy)
       call sylmnc ( rv_p_ocy , lmxcy ) 

Cchp1       call scg(lmxcg,w(ocg),w(oidxcg),w(ojcg))
       call scg ( lmxcg , rv_p_ocg , w ( oidxcg ) , w ( ojcg ) ) 



Cchp1       slat%ocg=ocg
       slat%rv_p_ocg => rv_p_ocg 

      slat%ojcg=ojcg
      slat%oidxcg=oidxcg
Cchp1       slat%ocy=ocy
       slat%rv_p_ocy => rv_p_ocy 



      end subroutine setcg


