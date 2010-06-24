Cchp1 cmdvar= integer,iv_p_,ojcg
Cchp1 cmdvar= integer,iv_p_,oidxcg
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
Cchp1         integer:: lnjcg , lnxcg , nlm , ojcg , oidxcg 
Cchp1          integer:: lnjcg , lnxcg , nlm , ojcg 
          integer:: lnjcg , lnxcg , nlm 
          integer,pointer :: iv_p_ojcg(:)

         integer,pointer :: iv_p_oidxcg(:)

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

Cchp1       call defi  (ojcg,     lnjcg)
       allocate(iv_p_ojcg(abs(lnjcg)))
       if (lnjcg<0) iv_p_ojcg(:)=0

Cchp1       call defi  (oidxcg,   lnxcg)
       allocate(iv_p_oidxcg(abs(lnxcg)))
       if (lnxcg<0) iv_p_oidxcg(:)=0


Cchp1       call sylmnc(w(ocy),lmxcy)
       call sylmnc ( rv_p_ocy , lmxcy ) 

Cchp1       call scg(lmxcg,w(ocg),w(oidxcg),w(ojcg))
Cchp1        call scg ( lmxcg , rv_p_ocg , w ( oidxcg ) , w ( ojcg ) ) 
Cchp1         call scg ( lmxcg , rv_p_ocg , iv_p_oidxcg , w ( ojcg ) ) 
         call scg ( lmxcg , rv_p_ocg , iv_p_oidxcg , iv_p_ojcg ) 





Cchp1       slat%ocg=ocg
       slat%rv_p_ocg => rv_p_ocg 

Cchp1       slat%ojcg=ojcg
       slat%iv_p_ojcg => iv_p_ojcg 

Cchp1       slat%oidxcg=oidxcg
       slat%iv_p_oidxcg => iv_p_oidxcg 

Cchp1       slat%ocy=ocy
       slat%rv_p_ocy => rv_p_ocy 



      end subroutine setcg


