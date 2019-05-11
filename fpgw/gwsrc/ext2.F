      integer function ifile_handle()
!! find open file handle
      implicit none
      integer:: i
      logical:: nexist
      do i=5001,99999
         inquire(unit=i,opened=nexist)
         if(.not.nexist) then
            ifile_handle=i
            return
         endif
      enddo
      call rx('ifile_hundle: we did not find open file hundle')
      end

C#define IEEE
C#define NEWVERS
C --- First line of derfc ---
C#ifdefC DOUBLE16
C      double precision function derfc (x)
CC- complement of error function, real*16 precision
CC ----------------------------------------------------------------
CCi Inputs
CCi   x
CCo Outputs
CCo   complement of error function
CCr Remarks
CCr   series for erf on the interval  0 to 1
CCr           with weighted error   1.28e-32
CCr                log weighted error  31.89
CCr      significant figures required  31.05
CCr           decimal places required  32.55
CCr   Adapted from July 1977 edition, W. Fullerton, c3,
CCr   Los Alamos Scientific Lab.
CCr   For real*8 double precision see below
CC ----------------------------------------------------------------
C      double precision x, erfcs(21), erfccs(59), erc2cs(49), sqeps,
C     .  sqrtpi, xmax, xsml, y,  d1mach, dcsevl, dexp, dlog, dsqrt
C      external d1mach, dcsevl, initds
C      save
C      data erf cs(  1) / -.4904612123 4691808039 9845440333 76 d-1 /
C      data erf cs(  2) / -.1422612051 0371364237 8247418996 31 d+0 /
C      data erf cs(  3) / +.1003558218 7599795575 7546767129 33 d-1 /
C      data erf cs(  4) / -.5768764699 7674847650 8270255091 67 d-3 /
C      data erf cs(  5) / +.2741993125 2196061034 4221607914 71 d-4 /
C      data erf cs(  6) / -.1104317550 7344507604 1353812959 05 d-5 /
C      data erf cs(  7) / +.3848875542 0345036949 9613114981 74 d-7 /
C      data erf cs(  8) / -.1180858253 3875466969 6317518015 81 d-8 /
C      data erf cs(  9) / +.3233421582 6050909646 4029309533 54 d-10/
C      data erf cs( 10) / -.7991015947 0045487581 6073747085 95 d-12/
C      data erf cs( 11) / +.1799072511 3961455611 9672454866 34 d-13/
C      data erf cs( 12) / -.3718635487 8186926382 3168282094 93 d-15/
C      data erf cs( 13) / +.7103599003 7142529711 6899083946 66 d-17/
C      data erf cs( 14) / -.1261245511 9155225832 4954248533 33 d-18/
C      data erf cs( 15) / +.2091640694 1769294369 1705002666 66 d-20/
C      data erf cs( 16) / -.3253973102 9314072982 3641600000 00 d-22/
C      data erf cs( 17) / +.4766867209 7976748332 3733333333 33 d-24/
C      data erf cs( 18) / -.6598012078 2851343155 1999999999 99 d-26/
C      data erf cs( 19) / +.8655011469 9637626197 3333333333 33 d-28/
C      data erf cs( 20) / -.1078892517 7498064213 3333333333 33 d-29/
C      data erf cs( 21) / +.1281188399 3017002666 6666666666 66 d-31/
Cc
Cc series for erc2       on the interval  2.50000e-01 to  1.00000e+00
Cc                                        with weighted error   2.67e-32
Cc                                         log weighted error  31.57
Cc                               significant figures required  30.31
Cc                                    decimal places required  32.42
Cc
C      data erc2cs(  1) / -.6960134660 2309501127 3915082619 7 d-1 /
C      data erc2cs(  2) / -.4110133936 2620893489 8221208466 6 d-1 /
C      data erc2cs(  3) / +.3914495866 6896268815 6114370524 4 d-2 /
C      data erc2cs(  4) / -.4906395650 5489791612 8093545077 4 d-3 /
C      data erc2cs(  5) / +.7157479001 3770363807 6089414182 5 d-4 /
C      data erc2cs(  6) / -.1153071634 1312328338 0823284791 2 d-4 /
C      data erc2cs(  7) / +.1994670590 2019976350 5231486770 9 d-5 /
C      data erc2cs(  8) / -.3642666471 5992228739 3611843071 1 d-6 /
C      data erc2cs(  9) / +.6944372610 0050125899 3127721463 3 d-7 /
C      data erc2cs( 10) / -.1371220902 1043660195 3460514121 0 d-7 /
C      data erc2cs( 11) / +.2788389661 0071371319 6386034808 7 d-8 /
C      data erc2cs( 12) / -.5814164724 3311615518 6479105031 6 d-9 /
C      data erc2cs( 13) / +.1238920491 7527531811 8016881795 0 d-9 /
C      data erc2cs( 14) / -.2690639145 3067434323 9042493788 9 d-10/
C      data erc2cs( 15) / +.5942614350 8479109824 4470968384 0 d-11/
C      data erc2cs( 16) / -.1332386735 7581195792 8775442057 0 d-11/
C      data erc2cs( 17) / +.3028046806 1771320171 7369724330 4 d-12/
C      data erc2cs( 18) / -.6966648814 9410325887 9586758895 4 d-13/
C      data erc2cs( 19) / +.1620854541 0539229698 1289322762 8 d-13/
C      data erc2cs( 20) / -.3809934465 2504919998 7691305772 9 d-14/
C      data erc2cs( 21) / +.9040487815 9788311493 6897101297 5 d-15/
C      data erc2cs( 22) / -.2164006195 0896073478 0981204700 3 d-15/
C      data erc2cs( 23) / +.5222102233 9958549846 0798024417 2 d-16/
C      data erc2cs( 24) / -.1269729602 3645553363 7241552778 0 d-16/
C      data erc2cs( 25) / +.3109145504 2761975838 3622741295 1 d-17/
C      data erc2cs( 26) / -.7663762920 3203855240 0956671481 1 d-18/
C      data erc2cs( 27) / +.1900819251 3627452025 3692973329 0 d-18/
C      data erc2cs( 28) / -.4742207279 0690395452 2565599996 5 d-19/
C      data erc2cs( 29) / +.1189649200 0765283828 8068307845 1 d-19/
C      data erc2cs( 30) / -.3000035590 3257802568 4527131306 6 d-20/
C      data erc2cs( 31) / +.7602993453 0432461730 1938527709 8 d-21/
C      data erc2cs( 32) / -.1935909447 6068728815 6981104913 0 d-21/
C      data erc2cs( 33) / +.4951399124 7733378810 0004238677 3 d-22/
C      data erc2cs( 34) / -.1271807481 3363718796 0862198988 8 d-22/
C      data erc2cs( 35) / +.3280049600 4695130433 1584165205 3 d-23/
C      data erc2cs( 36) / -.8492320176 8228965689 2479242239 9 d-24/
C      data erc2cs( 37) / +.2206917892 8075602235 1987998719 9 d-24/
C      data erc2cs( 38) / -.5755617245 6965284983 1281950719 9 d-25/
C      data erc2cs( 39) / +.1506191533 6392342503 5414405119 9 d-25/
C      data erc2cs( 40) / -.3954502959 0187969531 0428569599 9 d-26/
C      data erc2cs( 41) / +.1041529704 1515009799 8464505173 3 d-26/
C      data erc2cs( 42) / -.2751487795 2787650794 5017890133 3 d-27/
C      data erc2cs( 43) / +.7290058205 4975574089 9770368000 0 d-28/
C      data erc2cs( 44) / -.1936939645 9159478040 7750109866 6 d-28/
C      data erc2cs( 45) / +.5160357112 0514872983 7005482666 6 d-29/
C      data erc2cs( 46) / -.1378419322 1930940993 8964480000 0 d-29/
C      data erc2cs( 47) / +.3691326793 1070690422 5109333333 3 d-30/
C      data erc2cs( 48) / -.9909389590 6243654206 5322666666 6 d-31/
C      data erc2cs( 49) / +.2666491705 1953884133 2394666666 6 d-31/
Cc
Cc series for erfc       on the interval  0.          to  2.50000e-01
Cc                                        with weighted error   1.53e-31
Cc                                         log weighted error  30.82
Cc                               significant figures required  29.47
Cc                                    decimal places required  31.70
Cc
C      data erfccs(  1) / +.7151793102 0292477450 3697709496 d-1 /
C      data erfccs(  2) / -.2653243433 7606715755 8893386681 d-1 /
C      data erfccs(  3) / +.1711153977 9208558833 2699194606 d-2 /
C      data erfccs(  4) / -.1637516634 5851788416 3746404749 d-3 /
C      data erfccs(  5) / +.1987129350 0552036499 5974806758 d-4 /
C      data erfccs(  6) / -.2843712412 7665550875 0175183152 d-5 /
C      data erfccs(  7) / +.4606161308 9631303696 9379968464 d-6 /
C      data erfccs(  8) / -.8227753025 8792084205 7766536366 d-7 /
C      data erfccs(  9) / +.1592141872 7709011298 9358340826 d-7 /
C      data erfccs( 10) / -.3295071362 2528432148 6631665072 d-8 /
C      data erfccs( 11) / +.7223439760 4005554658 1261153890 d-9 /
C      data erfccs( 12) / -.1664855813 3987295934 4695966886 d-9 /
C      data erfccs( 13) / +.4010392588 2376648207 7671768814 d-10/
C      data erfccs( 14) / -.1004816214 4257311327 2170176283 d-10/
C      data erfccs( 15) / +.2608275913 3003338085 9341009439 d-11/
C      data erfccs( 16) / -.6991110560 4040248655 7697812476 d-12/
C      data erfccs( 17) / +.1929492333 2617070862 4205749803 d-12/
C      data erfccs( 18) / -.5470131188 7543310649 0125085271 d-13/
C      data erfccs( 19) / +.1589663309 7626974483 9084032762 d-13/
C      data erfccs( 20) / -.4726893980 1975548392 0369584290 d-14/
C      data erfccs( 21) / +.1435873376 7849847867 2873997840 d-14/
C      data erfccs( 22) / -.4449510561 8173583941 7250062829 d-15/
C      data erfccs( 23) / +.1404810884 7682334373 7305537466 d-15/
C      data erfccs( 24) / -.4513818387 7642108962 5963281623 d-16/
C      data erfccs( 25) / +.1474521541 0451330778 7018713262 d-16/
C      data erfccs( 26) / -.4892621406 9457761543 6841552532 d-17/
C      data erfccs( 27) / +.1647612141 4106467389 5301522827 d-17/
C      data erfccs( 28) / -.5626817176 3294080929 9928521323 d-18/
C      data erfccs( 29) / +.1947443382 2320785142 9197867821 d-18/
C      data erfccs( 30) / -.6826305642 9484207295 6664144723 d-19/
C      data erfccs( 31) / +.2421988887 2986492401 8301125438 d-19/
C      data erfccs( 32) / -.8693414133 5030704256 3800861857 d-20/
C      data erfccs( 33) / +.3155180346 2280855712 2363401262 d-20/
C      data erfccs( 34) / -.1157372324 0496087426 1239486742 d-20/
C      data erfccs( 35) / +.4288947161 6056539462 3737097442 d-21/
C      data erfccs( 36) / -.1605030742 0576168500 5737770964 d-21/
C      data erfccs( 37) / +.6063298757 4538026449 5069923027 d-22/
C      data erfccs( 38) / -.2311404251 6979584909 8840801367 d-22/
C      data erfccs( 39) / +.8888778540 6618855255 4702955697 d-23/
C      data erfccs( 40) / -.3447260576 6513765223 0718495566 d-23/
C      data erfccs( 41) / +.1347865460 2069650682 7582774181 d-23/
C      data erfccs( 42) / -.5311794071 1250217364 5873201807 d-24/
C      data erfccs( 43) / +.2109341058 6197831682 8954734537 d-24/
C      data erfccs( 44) / -.8438365587 9237891159 8133256738 d-25/
C      data erfccs( 45) / +.3399982524 9452089062 7359576337 d-25/
C      data erfccs( 46) / -.1379452388 0732420900 2238377110 d-25/
C      data erfccs( 47) / +.5634490311 8332526151 3392634811 d-26/
C      data erfccs( 48) / -.2316490434 4770654482 3427752700 d-26/
C      data erfccs( 49) / +.9584462844 6018101526 3158381226 d-27/
C      data erfccs( 50) / -.3990722880 3301097262 4224850193 d-27/
C      data erfccs( 51) / +.1672129225 9444773601 7228709669 d-27/
C      data erfccs( 52) / -.7045991522 7660138563 8803782587 d-28/
C      data erfccs( 53) / +.2979768402 8642063541 2357989444 d-28/
C      data erfccs( 54) / -.1262522466 4606192972 2422632994 d-28/
C      data erfccs( 55) / +.5395438704 5424879398 5299653154 d-29/
C      data erfccs( 56) / -.2380992882 5314591867 5346190062 d-29/
C      data erfccs( 57) / +.1099052830 1027615735 9726683750 d-29/
C      data erfccs( 58) / -.4867713741 6449657273 2518677435 d-30/
C      data erfccs( 59) / +.1525877264 1103575676 3200828211 d-30/
Cc
C      data sqrtpi / 1.772453850 9055160272 9816748334 115d0 /
C      data nterf, nterfc, nterc2, xsml, xmax, sqeps / 3*0, 3*0.d0 /
Cc
C      if (nterf .ne. 0) goto 10
C      eta = 0.1*sngl(d1mach(3))
C      nterf = initds (erfcs, 21, eta)
C      nterfc = initds (erfccs, 59, eta)
C      nterc2 = initds (erc2cs, 49, eta)
Cc
C      xsml = -dsqrt (-dlog(sqrtpi*d1mach(3)))
C      xmax = dsqrt (-dlog(sqrtpi*d1mach(1)) )
C      xmax = xmax - 0.5d0*dlog(xmax)/xmax - 0.01d0
C      sqeps = dsqrt (2.0d0*d1mach(3))
Cc
C   10 if (x .gt. xsml) goto 20
Cc
Cc === erfc(x) = 1.0 - erf(x)  for  x .lt. xsml ===
C      derfc = 2.0d0
C      return
Cc
C   20 if (x .gt. xmax) goto 40
C      y = dabs(x)
C      if (y .gt. 1.0d0) goto 30
Cc
Cc === erfc(x) = 1.0 - erf(x)  for abs(x) .le. 1.0 ===
C      if (y .lt. sqeps) derfc = 1.0d0 - 2.0d0*x/sqrtpi
C      if (y .ge. sqeps)
C     . derfc = 1.0d0 - x*(1.0d0 + dcsevl (2.d0*x*x-1.d0,erfcs, nterf))
C      return
Cc
Cc === erfc(x) = 1.0 - erf(x)  for  1.0 .lt. abs(x) .le. xmax ===
C   30 y = y*y
C      if (y .le. 4.d0) derfc = dexp(-y)/dabs(x) *
C     .  (0.5d0 + dcsevl ((8.d0/y-5.d0)/3.d0, erc2cs, nterc2))
C      if (y .gt. 4.d0) derfc = dexp(-y)/dabs(x) *
C     .   (0.5d0 + dcsevl (8.d0/y-1.d0, erfccs, nterfc))
C      if (x .lt. 0.d0) derfc = 2.0d0 - derfc
C      return
Cc
C   40 call errmsg ('DERFC: underflow', 1)
C      derfc = 0.d0
C      return
Cc
C      end
C#else
      double precision function derfc (x)
C- complement of error function, real*8 precision
C ----------------------------------------------------------------
Ci Inputs
Ci   x
Co Outputs
Co   complement of error function
Cr Remarks
Cr   erfcs: series for erf on the interval  0 to 1
Cr              with weighted error  7.10d-18
Cr                  log weighted error  17.15
Cr        significant figures required  16.31
Cr             decimal places required  17.71
Cr   erc2s: series for erc2 on the interval  .25 to 1
Cr             with weighted error   5.22d-17
Cr                  log weighted error  16.28
Cr         significant figures required  15.0
Cr             decimal places required  16.96
Cr   erfccs: series for erfc on the interval  0 to  .25
Cr              with weighted error  4.81d-17
Cr                  log weighted error  16.32
Cr         significant figures required  15.0
Cr             decimal places required  17.01
Cr   Adapted from July 1977 edition, W. Fullerton, c3,
Cr   Los Alamos Scientific Lab.
Cr   For real*16 double precision see below
C ----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension erfcs(13), erfccs(24), erc2cs(23)
      real eta
      external dcsevl, initds, d1mach
      integer iprint
      save

      data erf cs( 1) /   -.0490461212 34691808d0 /
      data erf cs( 2) /   -.1422612051 0371364d0 /
      data erf cs( 3) /    .0100355821 87599796d0 /
      data erf cs( 4) /   -.0005768764 69976748d0 /
      data erf cs( 5) /    .0000274199 31252196d0 /
      data erf cs( 6) /   -.0000011043 17550734d0 /
      data erf cs( 7) /    .0000000384 88755420d0 /
      data erf cs( 8) /   -.0000000011 80858253d0 /
      data erf cs( 9) /    .0000000000 32334215d0 /
      data erf cs(10) /   -.0000000000 00799101d0 /
      data erf cs(11) /    .0000000000 00017990d0 /
      data erf cs(12) /   -.0000000000 00000371d0 /
      data erf cs(13) /    .0000000000 00000007d0 /
c
      data erc2cs( 1) /   -.0696013466 02309501d0 /
      data erc2cs( 2) /   -.0411013393 62620893d0 /
      data erc2cs( 3) /    .0039144958 66689626d0 /
      data erc2cs( 4) /   -.0004906395 65054897d0 /
      data erc2cs( 5) /    .0000715747 90013770d0 /
      data erc2cs( 6) /   -.0000115307 16341312d0 /
      data erc2cs( 7) /    .0000019946 70590201d0 /
      data erc2cs( 8) /   -.0000003642 66647159d0 /
      data erc2cs( 9) /    .0000000694 43726100d0 /
      data erc2cs(10) /   -.0000000137 12209021d0 /
      data erc2cs(11) /    .0000000027 88389661d0 /
      data erc2cs(12) /   -.0000000005 81416472d0 /
      data erc2cs(13) /    .0000000001 23892049d0 /
      data erc2cs(14) /   -.0000000000 26906391d0 /
      data erc2cs(15) /    .0000000000 05942614d0 /
      data erc2cs(16) /   -.0000000000 01332386d0 /
      data erc2cs(17) /    .0000000000 00302804d0 /
      data erc2cs(18) /   -.0000000000 00069666d0 /
      data erc2cs(19) /    .0000000000 00016208d0 /
      data erc2cs(20) /   -.0000000000 00003809d0 /
      data erc2cs(21) /    .0000000000 00000904d0 /
      data erc2cs(22) /   -.0000000000 00000216d0 /
      data erc2cs(23) /    .0000000000 00000052d0 /
c
      data erfccs( 1) /   0.0715179310 202925d0 /
      data erfccs( 2) /   -.0265324343 37606719d0 /
      data erfccs( 3) /    .0017111539 77920853d0 /
      data erfccs( 4) /   -.0001637516 63458512d0 /
      data erfccs( 5) /    .0000198712 93500549d0 /
      data erfccs( 6) /   -.0000028437 12412769d0 /
      data erfccs( 7) /    .0000004606 16130901d0 /
      data erfccs( 8) /   -.0000000822 77530261d0 /
      data erfccs( 9) /    .0000000159 21418724d0 /
      data erfccs(10) /   -.0000000032 95071356d0 /
      data erfccs(11) /    .0000000007 22343973d0 /
      data erfccs(12) /   -.0000000001 66485584d0 /
      data erfccs(13) /    .0000000000 40103931d0 /
      data erfccs(14) /   -.0000000000 10048164d0 /
      data erfccs(15) /    .0000000000 02608272d0 /
      data erfccs(16) /   -.0000000000 00699105d0 /
      data erfccs(17) /    .0000000000 00192946d0 /
      data erfccs(18) /   -.0000000000 00054704d0 /
      data erfccs(19) /    .0000000000 00015901d0 /
      data erfccs(20) /   -.0000000000 00004729d0 /
      data erfccs(21) /    .0000000000 00001432d0 /
      data erfccs(22) /   -.0000000000 00000439d0 /
      data erfccs(23) /    .0000000000 00000138d0 /
      data erfccs(24) /   -.0000000000 00000048d0 /
c
      data sqrtpi /1.772453850 9055160d0/
      data nterf, nterfc, nterc2, xsml, xmax, sqeps /3*0, 3*0.d0/
c
      if (nterf .ne. 0d0) goto 10
      eta = 0.1*sngl(d1mach(3))
      nterf = initds (erfcs, 13, eta)
      nterfc = initds (erfccs, 24, eta)
      nterc2 = initds (erc2cs, 23, eta)

      xsml = -dsqrt (-dlog(sqrtpi*d1mach(3)))
      xmax = dsqrt (-dlog(sqrtpi*d1mach(1)))
      xmax = xmax - 0.5d0*dlog(xmax)/xmax - 0.01d0
      sqeps = dsqrt (2.0d0*d1mach(3))
c
   10 if (x .gt. xsml) goto 20
c
C --- derfc(x) = 1.0d0 - erf(x) for x .lt. xsml ---
      derfc = 2.d0
      return
c
   20 if (x .gt. xmax) goto 40
      y = dabs(x)
      if (y .gt. 1.0d0) goto 30
c
C --- derfc(x) = 1.0d0 - erf(x) for -1.d0 .le. x .le. 1.d0 ---
      if (y .lt. sqeps) derfc = 1.0d0 - 2.0d0*x/sqrtpi
      if (y .ge. sqeps) derfc = 1.0d0 -
     .  x*(1.0d0 + dcsevl (2.d0*x*x-1.d0, erfcs, nterf) )
      return
c
C --- derfc(x) = 1.0d0 - erf(x) for 1.d0 .lt. dabs(x) .le. xmax ---
   30 y = y*y
      if (y .le. 4.d0) derfc = dexp(-y)/dabs(x) *
     .  (0.5d0 + dcsevl ((8.d0/y-5.d0)/3.d0, erc2cs, nterc2) )
      if (y .gt. 4.d0) derfc = dexp(-y)/dabs(x) *
     .  (0.5d0 + dcsevl (8.d0/y-1.d0, erfccs, nterfc) )
      if (x .lt. 0.d0) derfc = 2.0d0 - derfc
      return
c
   40 if (iprint() .gt. 50) call errmsg ('DERFC: underflow', 1)
      derfc = 0.d0
      return
c
      end
C#endif
      double precision function derf (x)
      double precision derfc,x
      derf = 1d0 - derfc(x)
      end

      subroutine ishell(n,iarray)
      integer n
      integer iarray(1)
      integer lognb2,i,j,k,l,m,nn,it

      lognb2 = int(log(float(n+1))*1.4426950)
      m = n
      do  12  nn = 1, lognb2
        m = m/2
        k = n - m
        do  11  j = 1, k
          i = j
    3     continue
          l = i + m
          if (iarray(l) .lt. iarray(i)) then
            it = iarray(i)
            iarray(i) = iarray(l)
            iarray(l) = it
            i = i - m
            if (i .ge. 1) goto 3
          endif
   11   continue
   12 continue
      return
      end

      SUBROUTINE DMCPY(A,NCA,NRA,B,NCB,NRB,N,M)
C- general matrix copy
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   b,ncb,nrb is the right matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   n,m: the number of columns and rows, respectively, to calculate
Co Outputs:
Co   result matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix copy routine,
Cr   copying a subblock of matrix a to a subblock of matrix b.
Cr   Normally matrix nc{a,b} is the row dimension of matrix {a,b}
Cr   and nr{a,b} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop.
Cr
Cr   Example: Set 3-by-2 block of matrix c to constant z
Cr     call dmcpy(z,0,0,c,nc,1,3,2)
Cr   Note scalar z is represented by an array of 0 dimension
C ----------------------------------------------------------------
C
      INTEGER NCA,NRA,NCB,NRB,N,M
      DOUBLE PRECISION A(0:*), B(0:*)
      INTEGER I,J,IA,IB

      DO  200  I = N-1, 0, -1
        IA = I*NRA+M*NCA
        IB = I*NRB+M*NCB
      DO  200  J = M-1, 0, -1
        IA = IA-NCA
        IB = IB-NCB
        B(IB) = A(IA)
  200 CONTINUE
      RETURN
      END
      logical function bittst(n,bit)
C- returns true when a bit is set in an integer
C ----------------------------------------------------------------
Ci Inputs
Ci   n: integer
Ci   bit: a bit, ie 1,2,4,8,16, etc
Co Outputs
Co   bittst: true when bit in n is set, false otherwise
C ----------------------------------------------------------------
      integer n,bit
      bittst = (mod(n,bit+bit) - mod(n,bit) .eq. bit)
      end
      integer function getdig(n,i,base)
C- Returns the a digit from an integer
C ----------------------------------------------------------------
Ci Inputs
Ci   n,i,base
Co Outputs
Co   getdig = ith digit from n, base "base"; eg 4=getdig(12345,1,10)
C ----------------------------------------------------------------
      integer n,i,base
      getdig = mod(n/base**i,base)
      end
C#ifdef NEWVERS
      subroutine zhmul(ndim,hk,d,nid,ck)
C- Multiplies matrix ck = hk*d*hk, where d is a diagonal matrix.
C- Only the upper triangle of ck is calculated.
C ----------------------------------------------------------------
Ci Inputs
Ci   ndim,hk,d
Ci   nid: skip length between elements of d
Co Outputs
Co   ck
Cr Remarks
Cr   Hk and ck are dimensioned (ndim*ndim*2) (real followed by imag)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndim,nid
      double precision hk(ndim*ndim*2),ck(ndim*ndim*2),d(nid*ndim)
C Local parameters
      integer n2,i,j,k,jndim,j2ndim,kndim,k2ndim
      double precision xikr,xiki,xikim

      n2 = ndim**2
      do  100  j = 1, ndim
        jndim = (j-1)*ndim
        j2ndim = jndim+n2
      do  100  k = 1, ndim
        kndim = (k-1)*ndim
        k2ndim = kndim+n2
        xikr = d(k)*hk(jndim+k)
        xiki = d(k)*hk(j2ndim+k)
        xikim = -xiki
        do  110  i = 1, j
          ck(jndim+i) = ck(jndim+i) + xikr * hk(kndim+i) +
     .                                 xikim * hk(k2ndim+i)
          ck(j2ndim+i) = ck(j2ndim+i) + xiki * hk(kndim+i) +
     .                                   xikr * hk(k2ndim+i)
  110   continue
  100 continue
      end
C#elseC
C      subroutine zhmul(ndim,hk,d,nid,ck)
CC- Multiplies matrix ck = hk*d*hk, where d is a diagonal matrix
CC ----------------------------------------------------------------
CCi Inputs
CCi   ndim,hk,d
CCi   nid: skip length between elements of d
CCo Outputs
CCo   ck
CCr Remarks
CCr   Hk and ck are dimensioned (ndim,ndim,2) (real followed by imag)
CCr   The looping order is chosen so as to vectorize the inner loop
CC ----------------------------------------------------------------
CC     implicit none
CC Passed parameters
C      integer ndim,nid
C      double precision hk(ndim,ndim),ck(ndim,ndim),d(nid,1)
CC Local parameters
C      integer n2,i,j,k
C      double precision xikr,xiki
CC
C      n2 = ndim**2
C      do  100  k = 1, ndim
C        do  100  i = 1, ndim
C        xikr = d(1,k)*hk(i,k)
C        xiki = d(1,k)*hk(n2+i,k)
C#ifdefC BLAS
C        call daxpy(ndim, xikr,hk(   k,1),ndim,ck(i,1),   ndim)
C        call daxpy(ndim,-xiki,hk(n2+k,1),ndim,ck(i,1),   ndim)
C        call daxpy(ndim, xikr,hk(n2+k,1),ndim,ck(n2+i,1),ndim)
C        call daxpy(ndim, xiki,hk(   k,1),ndim,ck(n2+i,1),ndim)
C#elseC
C        do  110  j = 1, ndim
C          ck(i,j)    = ck(i,j)    + xikr*hk(k,j) - xiki*hk(n2+k,j)
C          ck(n2+i,j) = ck(n2+i,j) + xikr*hk(n2+k,j) + xiki*hk(k,j)
C  110 continue
C#endifC
C  100 continue
C      end
C#endif
      subroutine zmul0(ndim,a,b,c)
C- Multiplies complex matrix c = a*b,
C- All matrices have the same dimension
C ----------------------------------------------------------------
Ci Inputs
Ci   ndim,a,b
Co Outputs
Co   c
Cr Remarks
Cr   a,b, and c are dimensioned (ndim*ndim*2) (real followed by imag)
Cr   The looping order is chosen so as to vectorize the inner loop
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndim
      double precision a(ndim*ndim*2),b(ndim*ndim*2),c(ndim*ndim*2)
C Local parameters
      integer n2,i,j,k,jndim,j2ndim,kndim,k2ndim
      double precision bkj,bk2j,bk2jm

      n2 = ndim**2

C#ifdefC CRAY
C      call mxma(a,1,ndim,b,1,ndim,c,1,ndim,ndim,ndim,ndim)
C      call mxma(a,1,ndim,b(n2+1),1,ndim,c(n2+1),1,ndim,ndim,ndim,ndim)
C      call dinit(a,n2)
C      call mxma(a(n2+1),1,ndim,b(n2+1),1,ndim,a,1,ndim,ndim,ndim,ndim)
C      do 10 i = 1, n2
C        c(i) = c(i) - a(i)
C   10 continue
C      call dinit(a,n2)
C      call mxma(a(n2+1),1,ndim,b,1,ndim,a,1,ndim,ndim,ndim,ndim)
C      do 20 i = 1, n2
C        c(n2+i) = c(n2+i) + a(i)
C   20 continue
C#else
      do  100  j = 1, ndim
        jndim = (j-1)*ndim
        j2ndim = jndim+n2
      do  100  k = 1, ndim
        kndim = (k-1)*ndim
        k2ndim = kndim+n2
        bkj = b(jndim+k)
        bk2j = b(j2ndim+k)
        bk2jm = -bk2j
        do  110  i = 1, ndim
          c(jndim+i) = c(jndim+i) + bkj*a(kndim+i)
     .                                + bk2jm*a(k2ndim+i)
          c(j2ndim+i) = c(j2ndim+i) + bk2j*a(kndim+i)
     .                                  + bkj*a(k2ndim+i)
  110   continue
  100 continue
C#endif
      end
      subroutine zmul(a,nleft,b,nmid,c,nright)
C- General complex matrix copy c = a*b, cij = \sum aik*bkj
C ----------------------------------------------------------------
Ci Inputs
Ci   nleft,nmid,nright,a,b
Co Outputs
Co   c
Cr Remarks
Cr   c(nleft,nright) = a(nleft,nmid) * b(nmid,nright)
Cr   The looping order is chosen so as to vectorize the inner loop
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nleft,nmid,nright
      double precision a(nleft*nmid*2),b(nmid*nright*2),
     .                 c(nleft*nright*2)
C Local parameters
      integer nlm,nmr,nlr,i,j,k,jb,j2b,jc,j2c,ka,k2a
      double precision bkj,bk2j,bk2jm

      nlm = nleft * nmid
      nmr = nmid * nright
      nlr = nleft * nright

      do  100  j = 1, nright
        jc = (j-1) * nleft
        j2c = jc + nlr
        jb = (j-1) * nmid
        j2b = jb + nmr
      do  100  k = 1, nmid
        ka = (k-1) * nleft
        k2a = ka + nlm
        bkj = b(jb+k)
        bk2j = b(j2b+k)
        bk2jm = -bk2j
        do  110  i = 1, nleft
          c(jc+i)  = c(jc+i)  + bkj*a(ka+i)  + bk2jm*a(k2a+i)
          c(j2c+i) = c(j2c+i) + bk2j*a(ka+i) + bkj*a(k2a+i)
  110   continue
  100 continue
      end
      subroutine zmult(a,nleft,b,nmid,c,nright)
C- General complex matrix copy c = at*b,
C-                             cij = \sum atik*bkj = \sum aki*bkj
C- at is the transposed of a
C ----------------------------------------------------------------
Ci Inputs
Ci   nleft,nmid,nright,a,b
Co Outputs
Co   c
Cr Remarks
Cr   c(nleft,nright) = at(nleft,nmid) * b(nmid,nright)
Cr   c(nleft,nright) = a (nmid,nleft) * b(nmid,nright)
Cr   The looping order is chosen so as to vectorize the inner loop
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nleft,nmid,nright
      double precision a(nmid*nleft*2),b(nmid*nright*2),
     .                 c(nleft*nright*2)
C Local parameters
      integer nml,nmr,nlr,i,j,k,ia,i2a,jb,j2b,jc,j2c
      double precision cic,ci2c

      nml = nmid * nleft
      nmr = nmid * nright
      nlr = nleft * nright

      do  100  j = 1, nright
        jc = (j-1) * nleft
        j2c = jc + nlr
        jb = (j-1) * nmid
        j2b = jb + nmr
      do  100  i = 1, nleft
        ia = (i-1) * nmid
        i2a = ia + nml
        cic = 0
        ci2c = 0
        do  110  k = 1, nmid
          cic  = cic  + a(ia+k)*b(jb+k) - a(i2a+k)*b(j2b+k)
          ci2c = ci2c + a(i2a+k)*b(jb+k) + a(ia+k)*b(j2b+k)
  110   continue
        c(jc+i) = c(jc+i) + cic
        c(j2c+i) = c(j2c+i) + ci2c
  100 continue
      end
      SUBROUTINE DVCPY(A,NCA,B,NCB,N)
C- general vector copy
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca is the source vector and the number of
Ci      elements separating each element in the vector
Ci   b,ncb,nrb is the destination vector and the number of
Ci      elements separating each element in the vector
Ci   n: the number of elements to calculate
Co Outputs:
Co   result matrix stored in b
Cr Remarks:
Cr   This is a general-purpose vectore copy routine
Cr   Example: Set all elements of 3-by-2 matrix c to -1.d0
Cr     call dvcpy(-1.d0,0,c,1,3*2)
Cr   Example: Set block (n,m) of array a(p,m)=0, and a(i,i)=1, i=1,m
Cr     call dmcpy(0.d0,0,0,a,p,1,n,m)
Cr     call dvcpy(1.d0,0,a,p+1,m)
C ----------------------------------------------------------------
      INTEGER NCA,NCB,N
      DOUBLE PRECISION A(0:*), B(0:*)
      INTEGER I,IA,IB

      IA = N*NCA
      IB = N*NCB
      DO  200  I = N-1, 0, -1
        IA = IA-NCA
        IB = IB-NCB
        B(IB) = A(IA)
  200 CONTINUE

      return
      end
      subroutine ivshel(m,n,iarray,iwk,lopt)
C- shell sort of a array of integer vectors
C ----------------------------------------------------------------
Ci Inputs
Ci   iarray(m,n)
Ci   iwk: a work array of dimension n
Ci   lopt:if T, iwk returned as list of indices to iarray to sort it,
Ci              while iarray is unchanged.
Ci        if F, iarray returned sorted
Co Outputs
Co   iwk a table of indices to array iarray (lopt true)
C ----------------------------------------------------------------
      integer m,n
      logical lopt
      integer iarray(m,0:n-1),iwk(n)
      integer lognb2,i,j,k,l,n2,nn,it
      lognb2 = int(log(float(n+1))*1.4426950)
      n2 = n
      do  2  i = 1, n
        iwk(i) = i-1
    2 continue
      do  12  nn = 1, lognb2
        n2 = n2/2
        k = n - n2
        do  11  j = 1, k
          i = j
    3     continue
          l = i + n2
c      print *, 'test ',i,l,iwk(i),iwk(l)
          do  15  mm = 1, m
            if (iarray(mm,iwk(l)) - iarray(mm,iwk(i))) 16,15,11
   16       continue
            if (lopt) then
              it = iwk(i)
              iwk(i) = iwk(l)
              iwk(l) = it
c      print 800, (iwk(nnn), nnn=1,n)
c800   format(' swap', 11i5)
            else
c      print 800,  (iarray(nnn,i-1), nnn=1,m)
c      print 800,  (iarray(nnn,l-1), nnn=1,m)
              do  14  mmm = 1, m
                it = iarray(mmm,i-1)
                iarray(mmm,i-1) = iarray(mmm,l-1)
                iarray(mmm,l-1) = it
   14         continue
            endif
            i = i - n2
            if (i .ge. 1) goto 3
            goto 11
   15     continue
   11   continue
   12 continue

      return
      end



      SUBROUTINE DMPY(A,NCA,NRA,B,NCB,NRB,C,NCC,NRC,N,M,L)
C- matrix multiplication
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   b,ncb,nrb is the right matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   c,ncc,nrc is the product matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   n,m: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed n*m*l times.  No attempt is made to optimize
Cr   the outer loops, executed n*m times.
Cr     Examples: product of (n,l) subblock of a into (l,m) subblock of b
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,n,m,l)
Cr     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
Cr     To generate the tranpose of that product, use:
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,n,m,l)
C ----------------------------------------------------------------
      INTEGER NCA,NRA,NCB,NRB,NCC,NRC,N,M,L
      DOUBLE PRECISION A(0:*), B(0:*), C(0:*)
      DOUBLE PRECISION SUM
      INTEGER I,J,K,NAKPI,NBJPK
C
C#ifdefC CRAY
C      CALL MXMA(A,NRA,NCA,B,NRB,NCB,C,NRC,NCC,N,L,M)
C#else
      DO  200  I = N-1, 0, -1
      DO  200  J = M-1, 0, -1
        SUM = 0
        NAKPI = NRA*I
        NBJPK = NCB*J
        DO  210  K = L-1, 0, -1
          SUM = SUM + A(NAKPI)*B(NBJPK)
          NAKPI = NAKPI + NCA
          NBJPK = NBJPK + NRB
  210   CONTINUE
        C(I*NRC+J*NCC) = SUM
  200 CONTINUE
C#endif
      END
      double precision function dcsevl (x, a, n)
C- Evaluate the n-term Chebyshev series a at x.
C ----------------------------------------------------------------
Ci Inputs
Ci   x:  dble prec value at which the series is to be evaluated.
Ci   a:  dble prec array of n terms of a chebyshev series.
Ci       In evaluating a, only half the first coef is summed.
Ci   n:  number of terms in array a.
Co Outputs
Co   dcsevl
Cr Remarks
Cr   Adapted from R. Broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
C ----------------------------------------------------------------
      double precision a(n), x, twox, b0, b1, b2
c
      if (dabs(x) .gt. 1.1d0) call errmsg('DCSEVL:  x outside (-1,1)',2)
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do  10  i = 1, n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
   10 continue
      dcsevl = .5d0*(b0-b2)

      return
      end
      subroutine errmsg (messg, iopt)
C- Write error message to standard error device
C ----------------------------------------------------------------
Ci Inputs
Ci   iopt: 0, return without message printed
Ci         1, return with message printed
Ci         2, stop with message printed
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      character*(*) messg

      if (iopt .ne. 0) write(i1mach(4),*) messg
      if (iopt .lt. 2) return
      end
      function initds (dos, nos, eta)
C- Initialize things for Chebychev series
C ----------------------------------------------------------------
Ci Inputs
Ci   dos: dble prec array of nos coefficients in an orthogonal series.
Ci   nos: number of coefficients in dos.
Ci   eta: requested accuracy of series (real).
Co Outputs
Co   Returns number of terms necessary in series for spec'd eta
Cr Remarks
Cr   Initialize the double precision orthogonal series dos so that
Cr   initds is the number of terms needed to insure the error is no
Cr   larger than eta.  ordinarily eta will be chosen to be one-tenth
Cr   machine precision.
Cr   Adapted from June 1977 edition W. Fullerton,
Cr   c3, los alamos scientific lab.
C ----------------------------------------------------------------
      double precision dos(nos)
      real eta
c
      err = 0.
      do  10  ii = 1, nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err .gt. eta) goto 20
   10 continue
c
   20 continue
c     if (i .eq. nos) call errmsg('INITDS: eta may be too small',1)
      initds = i
c
      return
      end
      subroutine zampy(a,nca,nra,nia,b,ncb,nrb,nib,c,ncc,nrc,nic,n,m,l)
C- complex matrix multiplication: matrix incremented, not overwritten
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of a (in real words)
Ci   b,ncb,nrb is the right matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of b (in real words)
Ci   c,ncc,nrc is the product matrix and respectively the number of
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of c (in real words)
Ci   n,m: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored added into c
Cr Remarks:
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed n*m*l times.  No attempt is made to optimize
Cr   the outer loops, executed n*m times.
Cr     Examples: product of complex matrix c = a*b  (arrays b,c
Cr     dimensioned complex*16; a real*8 with imag following real)
Cr     call zampy(a,n,1,ndim**2,b,2*n,2,1,c,2*n,2,1,n,n,n)
Cr     To generate c = a*b
Cr     call zampy(a,2*n,2,1,b,2,2*n,1,c,2*n,2,1,n,n,n)
Cr   This version suitable for Cray
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,nia,ncb,nrb,nib,ncc,nrc,nic,n,m,l
      double precision a(0:*), b(0:*), c(0:*)
      integer i,j,k,nrci,nccj,ncbj,nrcicj
      double precision ar,ai

C --- Do multiplication ---
      do  20  k = l-1, 0, -1
      do  20  i = n-1, 0, -1
        ar = a(      nra*i + nca*k)
        ai = a(nia + nra*i + nca*k)
        nrci = nrc*i
C#ifdefC BLAS
C        nrbk = nrb*k
C        call daxpy(m, ar,b(    nrbk),ncb,c(nrci),    ncc)
C        call daxpy(m,-ai,b(nib+nrbk),ncb,c(nrci),    ncc)
C        call daxpy(m, ar,b(nib+nrbk),ncb,c(nic+nrci),ncc)
C        call daxpy(m, ai,b(    nrbk),ncb,c(nic+nrci),ncc)
C#else
        nccj = -ncc
        ncbj = -ncb + nrb*k
      do  20  j = m-1, 0, -1
        nccj = nccj + ncc
        ncbj = ncbj + ncb
        nrcicj = nrci + nccj
        c(nrcicj)     = c(nrcicj)     + ar*b(ncbj) - ai*b(nib+ncbj)
        c(nic+nrcicj) = c(nic+nrcicj) + ar*b(nib+ncbj) + ai*b(ncbj)
C#endif
   20 continue
      end
      subroutine zmpy(a,nca,nra,nia,b,ncb,nrb,nib,c,ncc,nrc,nic,n,m,l)
C- Double precision complex matrix multiplication
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of a (in real words)
Ci   b,ncb,nrb is the right matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of b (in real words)
Ci   c,ncc,nrc is the product matrix and respectively the number of
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of c (in real words)
Ci   n,m: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed n*m*l times.  No attempt is made to optimize
Cr   the outer loops, executed n*m times.
Cr     Examples: product of complex matrix c = a*b  (arrays b,c
Cr     dimensioned complex*16; a real*8 with imag following real)
Cr     call zmpy(a,n,1,ndim**2,b,2*n,2,1,c,2*n,2,1,n,n,n)
Cr     To generate c = a*b
Cr     call zmpy(a,2*n,2,1,b,2,2*n,1,c,2*n,2,1,n,n,n)
Cr   Warning: this routine has not been thoroughly checked!
Cr   This version suitable for Cray
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,nia,ncb,nrb,nib,ncc,nrc,nic,n,m,l
      double precision a(0:*), b(0:*), c(0:*)
      integer i,j,k,nrci,nccj,ncbj,nrcicj
      double precision ar,ai

C --- Initialize array to zero ---
      do  10  i = n-1, 0, -1
        nrci = nrc*i
        nccj = -ncc
      do  10  j = m-1, 0, -1
        nccj = nccj + ncc
        nrcicj = nrci + nccj
        c(nrcicj)     = 0
        c(nic+nrcicj) = 0
   10 continue

C --- Do multiplication ---
      do  20  k = l-1, 0, -1
      do  20  i = n-1, 0, -1
        ar = a(      nra*i + nca*k)
        ai = a(nia + nra*i + nca*k)
        nrci = nrc*i
C#ifdefC BLAS
C        nrbk = nrb*k
C        call daxpy(m, ar,b(    nrbk),ncb,c(nrci),    ncc)
C        call daxpy(m,-ai,b(nib+nrbk),ncb,c(nrci),    ncc)
C        call daxpy(m, ar,b(nib+nrbk),ncb,c(nic+nrci),ncc)
C        call daxpy(m, ai,b(    nrbk),ncb,c(nic+nrci),ncc)
C#else
        nccj = -ncc
        ncbj = -ncb + nrb*k
      do  20  j = m-1, 0, -1
        nccj = nccj + ncc
        ncbj = ncbj + ncb
        nrcicj = nrci + nccj
        c(nrcicj)     = c(nrcicj)     + ar*b(ncbj) - ai*b(nib+ncbj)
        c(nic+nrcicj) = c(nic+nrcicj) + ar*b(nib+ncbj) + ai*b(ncbj)
C#endif
   20 continue
      end

      function ran1()
C- Return a random deviate between 0.0 and 1.0.
C ----------------------------------------------------------------
Ci Inputs
Co Outputs
Co   ran1
Cr Remarks
Cr   Algorithm from Knuth; adapted here from Numerical Recipes, chapter
Cr   Uses three linear
Cr   congruential generators (two for high and low order, a third
Cr   to shuffle).  Use ran1in to initialize.
C ----------------------------------------------------------------
C     implicit none
C Local
      real ran1
      integer j
      integer m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
      parameter (m1=259200,ia1=7141,ic1=54773)
      parameter (m2=134456,ia2=8121,ic2=28411)
      parameter (m3=243000,ia3=4561,ic3=51349)
C Static
      integer ix1,ix2,ix3
      real r(97),rm1,rm2
      common /dran1/ ix1,ix2,ix3,rm1,rm2,r

C Generate next number for each sequence;
c use third to generate random integer between 1 and 97
      ix1 = mod(ia1*ix1+ic1,m1)
      ix2 = mod(ia2*ix2+ic2,m2)
      ix3 = mod(ia3*ix3+ic3,m3)
      j = 1 + (97*ix3)/m3
C#ifdefC TEST
C      if (j .gt. 97 .or. j .lt. 1) pause
C#endif
C Return the table entry ...
      ran1 = r(j)
C And refill it.
      r(j) = (float(ix1) + float(ix2)*rm2)*rm1
      return
      end
      subroutine ran1in(iseed)
C- A simple one-parameter initializer for ran1
C ----------------------------------------------------------------
Ci Inputs
Ci   iseed
Co Outputs
Co   ran1 is set up
Cr Remarks
C ----------------------------------------------------------------
C     implicit none
C Passed
      integer iseed
C Local
      integer j
      integer m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
      parameter (m1=259200,ia1=7141,ic1=54773)
      parameter (m2=134456,ia2=8121,ic2=28411)
      parameter (m3=243000,ia3=4561,ic3=51349)
C To preserve
      integer ix1,ix2,ix3
      real r(97),rm1,rm2
      common /dran1/ ix1,ix2,ix3,rm1,rm2,r

      rm1 = 1./m1
      rm2 = 1./m2
C Seed the first, second and third sequences
      ix1 = mod(ic1-iseed,m1)
      ix1 = mod(ia1*ix1+ic1,m1)
      ix2 = mod(ix1,m2)
      ix1 = mod(ia1*ix1+ic1,m1)
      ix3 = mod(ix1,m3)
C Fill table with sequential uniform deviates generated by first two
      do  11  j = 1, 97
        ix1 = mod(ia1*ix1+ic1,m1)
        ix2 = mod(ia2*ix2+ic2,m2)
        r(j) = (float(ix1) + float(ix2)*rm2)*rm1
   11 continue
      return
      end
C#ifdefC TEST
C      PROGRAM D7R2
CC     Driver for routine RAN1
CC     Calculates pi statistically using volume of unit n-sphere
C      parameter(pi=3.1415926)
C      dimension iy(3),yprob(3)
C      fnc(x1,x2,x3,x4) = sngl(dsqrt(
C     .  dble(x1)**2+dble(x2)**2+dble(x3)**2+dble(x4)**2))
CC
C      call ran1in(-1)
C      do  11  i = 1, 3
C        iy(i) = 0
C   11 continue
C      write(*,'(1x,/,t15,a)') 'Volume of unit n-sphere, n=2,3,4'
C      write(*,'(1x,/,t3,a,t17,a,t26,a,t37,a)')
C     .  '# points','pi','(4/3)*pi','(1/2)*pi^2'
C      do  14  j = 1, 15
C        do  12  k=2**(j-1),2**j
C          x1 = ran1()
C          x2 = ran1()
C          x3 = ran1()
C          x4 = ran1()
C          if (fnc(x1,x2,0.0,0.0) .lt. 1.0) iy(1) = iy(1)+1
C          if (fnc(x1,x2,x3,0.0) .lt. 1.0) iy(2) = iy(2)+1
C          if (fnc(x1,x2,x3,x4) .lt. 1.0) iy(3) = iy(3)+1
C   12   continue
C        do  13  i = 1, 3
C          yprob(i) = 1.0*(2**(i+1))*iy(i)/(2**j)
C   13   continue
C        write(*,'(1x,i8,3f12.6)') 2**j,(yprob(i),i=1,3)
C   14 continue
C      write(*,'(1x,/,t4,a,3f12.6,/)') 'actual',pi,4.0*pi/3.0,0.5*(pi**
C      end
C#endif


c$$$      integer function amix_notused(wk,nelts,nmix,mmix,beta,norm,t,kpvt,ipr,tm,
c$$$     .                      rmsdel)
c$$$C- Anderson mixing of a vector
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   nmix: number of previous iterations to fold into mix
c$$$Ci         nmix = 0 => linear mixing (x* = x0)
c$$$Ci         Practical experience for LMTO programs shows that nmix>2
c$$$Ci         does not work well.
c$$$Ci   mmix: maximum number of previous iterations to fold into mix
c$$$Ci         (used for dimensioning of work array)
c$$$Ci   nelts:number of elements to mix
c$$$Ci   wk:   array of dimension (nelts,1+mmix,2) where:
c$$$Ci         wk(*,i,1) holds f(xi) (see remarks)
c$$$Ci         wk(*,i,2) holds   xi  (see remarks)
c$$$Ci   beta: new x is beta f(x*) + (1-beta) x*
c$$$Ci   tm:   upper limit to any tj:  if any tj exceeds tm, effective
c$$$Ci         nmix is set to zero.
c$$$Co Outputs
c$$$Co   f(x_i) => f(x_i+1); x_i => x_i+1
c$$$Co   wk(*,1,2): new x = beta f(x*) + (1 - beta) x*
c$$$Co   rmsdel:rms difference between x_0 and f(x_0)
c$$$Co   amix:  returns effective nmix (see input tm)
c$$$Cr Remarks
c$$$Cr   Given a vector function f(x), where x is some
c$$$Cr   vector, we want to find x* such that f(x*) = x*.  We want to find
c$$$Cr   x* with the minimum number of computations of f(x).  Supposing
c$$$Cr   that we have x0,f(x0); x1,f(x1); x2,f(x2); ...; x_n+1,f(x_n+1).
c$$$Cr   (Usually x_j corresponds to x at the jth previous iteration.)
c$$$Cr   We take a linear combination x* of x_0, x_1, x_2, ... x_n that
c$$$Cr   minimizes <(x* - f(x*))^2>.  We then seek t_1, t_2, ... t_n in
c$$$Cr     x* = x_0 -
c$$$Cr   To evaluate f(x*) we linearize d(x) = f(x)-x as
c$$$Cr     f(x*)-x*  =  d*  =  d_0  -
c$$$Cr   Then
c$$$Cr     < (d_0 -
c$$$Cr   constitute a set of n simultaneous equations in the n unknowns t_k.
c$$$Cr   Note that d's enter into these equations, not the f's.
c$$$Cr   The d's are stored on disk, unit n1; the x's on unit n2.
c$$$Cr   Given the t_k's, x* can be estimated from (2).  To dampen
c$$$Cr   instablities, a linear combination of (1) and (2) is taken as
c$$$Cr       beta f(x*)  +  (1-beta) x*  = beta d*  +  x*           (4)
c$$$Cr   beta is an input parameter, which can usually be taken to be 1.
c$$$C ----------------------------------------------------------------
c$$$C     implicit none
c$$$C Passed parameters
c$$$      integer nelts,mmix,nmix,kpvt(1),ipr
c$$$      double precision norm(mmix,mmix),wk(nelts,0:mmix+1,2),t(mmix),
c$$$     .  beta,tm,rmsdel
c$$$C Local parameters
c$$$      integer i,j,nwk,jmix
c$$$      character*1 mixis0
c$$$      double precision dsqrt,ddot
c$$$
c$$$      nwk = nelts*(mmix+2)
c$$$
c$$$C --- d_0 = f-x  =>  wk(*,0,1) ---
c$$$      call daxpy(nelts,-1d0,wk(1,0,2),1,wk,1)
c$$$
c$$$C --- Copy x_0 and d_0 to end of wk:  x*, d* constructed here ---
c$$$      call dmcpy(wk,nwk,1,wk(1,mmix+1,1),nwk,1,nelts,2)
c$$$
c$$$C --- Make < (d_0 - d_j) (d_0 - d_k) > and  < d_0 (d_0 - d_j) >  ---
c$$$      do  20  i = 1, nmix
c$$$        t(i) = ddot(nelts,wk(1,0,1),1,wk(1,0,1),1) -
c$$$     .         ddot(nelts,wk(1,0,1),1,wk(1,i,1),1)
c$$$      do  20  j = 1, nmix
c$$$        norm(i,j) =  ddot(nelts,wk(1,0,1),1,wk(1,0,1),1)
c$$$     .             - ddot(nelts,wk(1,0,1),1,wk(1,j,1),1)
c$$$     .             - ddot(nelts,wk(1,i,1),1,wk(1,0,1),1)
c$$$     .             + ddot(nelts,wk(1,i,1),1,wk(1,j,1),1)
c$$$   20 continue
c$$$
c$$$C --- Solve the simultaneous equations ---
c$$$      call dsifa(norm,mmix,nmix,kpvt,i)
c$$$      if (i .ne. 0) then
c$$$        if (ipr .ge. 20) print *,  'AMIX: normal eqns singular'
c$$$      else
c$$$        call dsisl(norm,mmix,nmix,kpvt,t)
c$$$      endif
c$$$
c$$$C --- Set jmix = effective nmix to zero if any t_j exceeds tm ---
c$$$      jmix = nmix
c$$$      mixis0 = ' '
c$$$      do  30  j = 1, nmix
c$$$        if (dabs(t(j)) .gt. dabs(tm)) jmix = 0
c$$$   30 continue
c$$$      if (jmix .ne. nmix) mixis0 = '*'
c$$$      amix_notused = jmix
c$$$
c$$$C --- Make (d,x)* = (d,x)_0 -
c$$$      do  40  j = 1, jmix
c$$$        call dscal(nelts,1-t(j),wk(1,mmix+1,1),1)
c$$$        call daxpy(nelts,t(j),wk(1,j,1),1,wk(1,mmix+1,1),1)
c$$$        call dscal(nelts,1-t(j),wk(1,mmix+1,2),1)
c$$$        call daxpy(nelts,t(j),wk(1,j,2),1,wk(1,mmix+1,2),1)
c$$$   40 continue
c$$$
c$$$C --- Copy arrays to new positions ---
c$$$      do  50  i = mmix, 1, -1
c$$$        call dmcpy(wk(1,i-1,1),nwk,1,wk(1,i,1),nwk,1,nelts,2)
c$$$   50 continue
c$$$
c$$$C --- x* + beta d* ---
c$$$   10 call daxpy(nelts,beta,wk(1,mmix+1,1),1,wk(1,mmix+1,2),1)
c$$$
c$$$C --- Calculate rms change ---
c$$$      rmsdel = dsqrt(ddot(nelts,wk,1,wk,1)/nelts)
c$$$
c$$$C --- Printout ---
c$$$      if (ipr .lt. 20) goto 60
c$$$ca    write(*,133) nmix,mixis0,mmix,nelts,beta,tm,rmsdel
c$$$      write(*,133) nmix,mixis0,mmix,nelts,beta,rmsdel
c$$$ca133 format(/' AMIX:  nmix=',i1,a1,' mmix=',i1,'  nelts=',i3,
c$$$ca   .        '  beta=',f7.5,'  tm=',f8.5,'  rmsdel=',1pd8.2)
c$$$  133 format(/' AMIX:  nmix=',i1,a1,' mmix=',i1,'  nelts=',i3,
c$$$     .        '  beta=',f7.5,'  rmsdel=',1pd8.2)
c$$$      if (nmix .gt. 0) write(*,134) (t(j), j=1,nmix)
c$$$  134 format(3x,'tj:',7(f8.5,2x))
c$$$
c$$$      if (ipr .le. 30) goto 60
c$$$      write(*,110)
c$$$      do  12  i = 1, nelts
c$$$        if (dabs(wk(i,0,1)) + dabs(wk(i,mmix+1,2)-wk(i,0,2)) .ge. 5.d-7)
c$$$     .  write(*,111) i,wk(i,0,2),wk(i,0,2)+wk(i,0,1),
c$$$     .                 wk(i,0,1),wk(i,mmix+1,2)
c$$$   12 continue
c$$$
c$$$C --- Restore d* and x* + beta d* from end of wk ---
c$$$   60 call dmcpy(wk(1,mmix+1,1),nwk,1,wk,nwk,1,nelts,2)
c$$$
c$$$  104 format(1p,4d18.11)
c$$$  111 format(i5,4f14.6)
c$$$  110 format(14x,'OLD',11X,' NEW',9X,'DIFF',10X,'MIXED')
c$$$      end

!!#############################################################################
      double precision function tripl(a,b,c)
!! == tripl (determinant of 3x3 matrix) ==
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(0:2),b(0:2),c(0:2)
c     tripl=a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)
c     .     -a(3)*b(2)*c(1)-a(2)*b(1)*c(3)-a(1)*b(3)*c(2)
c     ... g77 needs this rewriting
      tmp = 0.d0
      do i = 0,2
         j = mod(i + 1, 3)
         k = mod(i + 2, 3)
         tmp = tmp + a(i) * (b(j)*c(k) - b(k)*c(j))
      enddo
      tripl = tmp
      end

!!#############################################################################
      integer(4) function llxxx(ilm)
      integer(4),parameter :: lmx=50
      integer(4),save:: lla((lmx+1)**2)
      logical:: init=.true.
      if(ilm>(lmx+1)**2) call rx( 'll: ilm too large')
      if(init) then
         do l=0,lmx
            lini= l**2 + 1
            lend=(l+1)**2
            lla(lini:lend)=l
         enddo
      endif
      llxxx = lla(ilm)
      return
      end
