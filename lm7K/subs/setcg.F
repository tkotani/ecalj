      subroutine setcg(slat,lmxcg,lmxcy)

      use m_struc_def  !Cgetarg

C- Allocate space for, and make Clebsch-Gordan coeffs
C     implicit none
      integer lmxcg,lmxcy
      type(s_lat)::slat

          integer:: lnjcg , lnxcg , nlm 
          integer,pointer :: iv_p_ojcg(:)

         integer,pointer :: iv_p_oidxcg(:)

        real(8),pointer :: rv_p_ocy(:)

       real(8),pointer :: rv_p_ocg(:)


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
       allocate(rv_p_ocy(abs(nlm)))
       if (nlm<0) rv_p_ocy(:)=0.0d0

       allocate(rv_p_ocg(abs(lnjcg)))
       if (lnjcg<0) rv_p_ocg(:)=0.0d0

       allocate(iv_p_ojcg(abs(lnjcg)))
       if (lnjcg<0) iv_p_ojcg(:)=0

       allocate(iv_p_oidxcg(abs(lnxcg)))
       if (lnxcg<0) iv_p_oidxcg(:)=0


       call sylmnc ( rv_p_ocy , lmxcy ) 

         call scg ( lmxcg , rv_p_ocg , iv_p_oidxcg , iv_p_ojcg ) 





       slat%rv_p_ocg => rv_p_ocg 

       slat%iv_p_ojcg => iv_p_ojcg 

       slat%iv_p_oidxcg => iv_p_oidxcg 

       slat%rv_p_ocy => rv_p_ocy 



      end subroutine setcg


