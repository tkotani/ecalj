      subroutine setcg(slat,lmxcg,lmxcy)
      use m_work,only: rv_a_ocy, rv_a_ocg, iv_a_ojcg, iv_a_oidxcg
      use m_struc_def  
C- Allocate space for, and make Clebsch-Gordan coeffs
C     implicit none
      integer lmxcg,lmxcy
      type(s_lat)::slat
      integer:: lnjcg , lnxcg , nlm
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
      allocate(rv_a_ocy(abs(nlm)))
      allocate(rv_a_ocg(abs(lnjcg)))
      allocate(iv_a_ojcg(abs(lnjcg)))
      allocate(iv_a_oidxcg(abs(lnxcg)))
      call sylmnc ( rv_a_ocy , lmxcy )
      call scg ( lmxcg , rv_a_ocg , iv_a_oidxcg , iv_a_ojcg )
      end subroutine setcg


