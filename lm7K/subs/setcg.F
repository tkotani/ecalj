Cdelw1 var ['ocy', '0', 'ocy', 'nlm', 'real(8)', 'rv_w_', 'setcg']
Cdelw1 var ['ocg', '0', 'ocg', 'lnjcg', 'real(8)', 'rv_w_', 'setcg']
Cdelw1 var ['ojcg', '0', 'ojcg', 'lnjcg', 'integer', 'iv_w_', 'setcg']
Cdelw1 var ['oidxcg', '0', 'oidxcg', 'lnxcg', 'integer', 'iv_w_', 'setcg']
Cdelw1 undel [['wref'], 'ocg', 'setcg']
Cdelw1 undel [['wref'], 'ojcg', 'setcg']
Cdelw1 undel [['wref'], 'oidxcg', 'setcg']
Cdelw1 undel [['wref'], 'ocy', 'setcg']
Cgetarg...info...           structure ['setcg', 'lat', 'slat']
Cgetarg...info...           use_to_add ['setcg', 'm_struc_def']
      subroutine setcg(slat,lmxcg,lmxcy)

      use m_struc_def  !Cgetarg

C- Allocate space for, and make Clebsch-Gordan coeffs
C     implicit none
      integer lmxcg,lmxcy
Cgetarg       double precision slat(1)
      type(s_lat)::slat

Cdelw1 do not change  [['wref'], 'ocy', 'setcg']
Cdelw1 do not change  [['wref'], 'ocg', 'setcg']
Cdelw1 do not change  [['wref'], 'ojcg', 'setcg']
Cdelw1 do not change  [['wref'], 'oidxcg', 'setcg']
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
Cdelw1 do not change ocy because of ['wref']
      call defrr (ocy,      nlm)
Cdelw1 do not change ocg because of ['wref']
      call defrr (ocg,      lnjcg)
Cdelw1 do not change ojcg because of ['wref']
      call defi  (ojcg,     lnjcg)
Cdelw1 do not change oidxcg because of ['wref']
      call defi  (oidxcg,   lnxcg)

Cdelw1 do not change ocy because of ['wref']
      call sylmnc(w(ocy),lmxcy)
Cdelw1 do not change ocg because of ['wref']
Cdelw1 do not change ojcg because of ['wref']
Cdelw1 do not change oidxcg because of ['wref']
      call scg(lmxcg,w(ocg),w(oidxcg),w(ojcg))

Cgetarg...info...           integer :: ocg , integer(8) :: slat%ocg 1
Cgetarg...info...           integer :: ojcg , integer(8) :: slat%ojcg 1
Cgetarg...info...           integer :: oidxcg , integer(8) :: slat%oidxcg 1
Cgetarg...info...           integer :: ocy , integer(8) :: slat%ocy 1
Cgetarg       call pack5('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)

Cdelw1 do not change ocg because of ['wref']
      slat%ocg=ocg
Cdelw1 do not change ojcg because of ['wref']
      slat%ojcg=ojcg
Cdelw1 do not change oidxcg because of ['wref']
      slat%oidxcg=oidxcg
Cdelw1 do not change ocy because of ['wref']
      slat%ocy=ocy


Cgetarg       end
Cdelw1 w_varlist remains: ocy ocg ojcg oidxcg
Cdelw1 w_varlistundel: [ocg ojcg oidxcg ocy]
Cdelw1 w_varlist (undel), remains: [ocy ocg ojcg oidxcg]
Cdelw1 w_varlist (del), remains: (None)
Cdelw1 not deallocate oidxcg because of [wref]
Cdelw1 not deallocate ojcg because of [wref]
Cdelw1 not deallocate ocg because of [wref]
Cdelw1 not deallocate ocy because of [wref]
      end subroutine setcg


