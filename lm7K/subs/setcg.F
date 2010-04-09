Cgetarg...info...           structure ['setcg', 'lat', 'slat']
Cgetarg...info...           use_to_add ['setcg', 'm_struc_def']
      subroutine setcg(slat,lmxcg,lmxcy)
       
       use m_struc_def  !Cgetarg

C- Allocate space for, and make Clebsch-Gordan coeffs
C     implicit none
      integer lmxcg,lmxcy
Cgetarg       double precision slat(1)
       type(s_lat)::slat

      integer lnjcg,lnxcg,nlm,ocy,ocg,ojcg,oidxcg
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
      call defrr (ocy,      nlm)
      call defrr (ocg,      lnjcg)
      call defi  (ojcg,     lnjcg)
      call defi  (oidxcg,   lnxcg)
      
      call sylmnc(w(ocy),lmxcy)
      call scg(lmxcg,w(ocg),w(oidxcg),w(ojcg))
      
Cgetarg...info...           integer :: ocg , integer(8) :: slat%ocg 1
Cgetarg...info...           integer :: ojcg , integer(8) :: slat%ojcg 1
Cgetarg...info...           integer :: oidxcg , integer(8) :: slat%oidxcg 1
Cgetarg...info...           integer :: ocy , integer(8) :: slat%ocy 1
Cgetarg       call pack5('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)
       
       slat%ocg=ocg 
       slat%ojcg=ojcg 
       slat%oidxcg=oidxcg 
       slat%ocy=ocy 


Cgetarg       end
       end subroutine setcg 


