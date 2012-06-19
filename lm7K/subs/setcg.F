      subroutine setcg(slat,lmxcg,lmxcy)

      use m_struc_def  !Cgetarg

C- Allocate space for, and make Clebsch-Gordan coeffs
C     implicit none
      integer lmxcg,lmxcy
      type(s_lat)::slat

      integer:: lnjcg , lnxcg , nlm
ckino Dec.14.2011:            integer,pointer :: iv_p_ojcg(:) =>NULL()

ckino Dec.14.2011:           integer,pointer :: iv_p_oidxcg(:) =>NULL()

ckino Dec.14.2011:          real(8),pointer :: rv_p_ocy(:) =>NULL()

ckino Dec.14.2011:         real(8),pointer :: rv_p_ocg(:) =>NULL()


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
ckino Dec.14.2011:         allocate(rv_p_ocy(abs(nlm)))
ckino Dec.27.2011:         allocate(slat%rv_p_ocy(abs(nlm)))
      allocate(slat%rv_a_ocy(abs(nlm)))
ckino Dec.14.2011:         if (nlm<0) rv_p_ocy(:)=0.0d0
ckino Dec.27.2011:         if (nlm<0) slat%rv_p_ocy(:)=0.0d0
      if (nlm<0) slat%rv_a_ocy(:)=0.0d0

ckino Dec.14.2011:         allocate(rv_p_ocg(abs(lnjcg)))
ckino Dec.27.2011:         allocate(slat%rv_p_ocg(abs(lnjcg)))
      allocate(slat%rv_a_ocg(abs(lnjcg)))
ckino Dec.14.2011:         if (lnjcg<0) rv_p_ocg(:)=0.0d0
ckino Dec.27.2011:         if (lnjcg<0) slat%rv_p_ocg(:)=0.0d0
      if (lnjcg<0) slat%rv_a_ocg(:)=0.0d0

ckino Dec.14.2011:         allocate(iv_p_ojcg(abs(lnjcg)))
ckino Dec.28.2011:         allocate(slat%iv_p_ojcg(abs(lnjcg)))
      allocate(slat%iv_a_ojcg(abs(lnjcg)))
ckino Dec.14.2011:         if (lnjcg<0) iv_p_ojcg(:)=0
ckino Dec.28.2011:         if (lnjcg<0) slat%iv_p_ojcg(:)=0
      if (lnjcg<0) slat%iv_a_ojcg(:)=0

ckino Dec.14.2011:         allocate(iv_p_oidxcg(abs(lnxcg)))
ckino Dec.27.2011:         allocate(slat%iv_p_oidxcg(abs(lnxcg)))
      allocate(slat%iv_a_oidxcg(abs(lnxcg)))
ckino Dec.14.2011:         if (lnxcg<0) iv_p_oidxcg(:)=0
ckino Dec.27.2011:         if (lnxcg<0) slat%iv_p_oidxcg(:)=0
      if (lnxcg<0) slat%iv_a_oidxcg(:)=0


ckino Dec.14.2011:         call sylmnc ( rv_p_ocy , lmxcy )
ckino Dec.27.2011:         call sylmnc ( slat%rv_p_ocy , lmxcy )
      call sylmnc ( slat%rv_a_ocy , lmxcy )

ckino Dec.14.2011:           call scg ( lmxcg , rv_p_ocg , iv_p_oidxcg , iv_p_ojcg )
ckino Dec.14.2011:           call scg ( lmxcg , rv_p_ocg , iv_p_oidxcg , slat%iv_p_ojcg )
ckino Dec.14.2011:           call scg ( lmxcg , rv_p_ocg , slat%iv_p_oidxcg , slat%iv_p_ojcg )
ckino Dec.27.2011:           call scg ( lmxcg , slat%rv_p_ocg , slat%iv_p_oidxcg , slat%iv_p_ojcg )
ckino Dec.27.2011:           call scg ( lmxcg , slat%rv_p_ocg , slat%iv_a_oidxcg , slat%iv_p_ojcg )
ckino Dec.28.2011:           call scg ( lmxcg , slat%rv_a_ocg , slat%iv_a_oidxcg , slat%iv_p_ojcg )
      call scg ( lmxcg , slat%rv_a_ocg , slat%iv_a_oidxcg , slat%iv_a_ojcg )





ckino Dec.14.2011:         slat%rv_p_ocg => rv_p_ocg

ckino Dec.14.2011:         slat%iv_p_ojcg => iv_p_ojcg

ckino Dec.14.2011:         slat%iv_p_oidxcg => iv_p_oidxcg

ckino Dec.14.2011:         slat%rv_p_ocy => rv_p_ocy



      end subroutine setcg


