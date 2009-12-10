      double precision function dci (x)
c december 1980 edition, w. fullerton, bell labs.
      double precision x, cics(19), xsml, y, f, g, sinx, dcsevl,
     1  d1mach, dcos, dlog, dsin, dsqrt
      external d1mach, dcos, dcsevl, dlog, dsin, dsqrt, initds
c
c series for ci   on the interval  0.00000e+00 to  1.60000e+01
c                                        with weighted error   9.26e-33
c                                         log weighted error  32.03
c                               significant figures required  32.06
c                                    decimal places required  32.67
c
      data ci  cs(  1) / -0.3400428185 6055363156 2810766331 29873d0/
      data ci  cs(  2) / -1.0330216640 1177456807 1592710401 63751d0/
      data ci  cs(  3) /  0.1938822265 9917082876 7158746060 81709d0/
      data ci  cs(  4) / -0.0191826043 6019865893 9463462701 75301d0/
      data ci  cs(  5) /  0.0011078925 2584784967 1840980992 66118d0/
      data ci  cs(  6) / -0.0000415723 4558247208 8038402318 14601d0/
      data ci  cs(  7) /  0.0000010927 8524300228 7152955789 66285d0/
      data ci  cs(  8) / -0.0000000212 3285954183 4652196012 80329d0/
      data ci  cs(  9) /  0.0000000003 1733482164 3485448651 29873d0/
      data ci  cs( 10) / -0.0000000000 0376141547 9876836993 81798d0/
      data ci  cs( 11) /  0.0000000000 0003622653 4884839643 36956d0/
      data ci  cs( 12) / -0.0000000000 0000028911 5284936518 52433d0/
      data ci  cs( 13) /  0.0000000000 0000000194 3278606764 94420d0/
      data ci  cs( 14) / -0.0000000000 0000000001 1151831826 50184d0/
      data ci  cs( 15) /  0.0000000000 0000000000 0055278588 87706d0/
      data ci  cs( 16) / -0.0000000000 0000000000 0000239070 13943d0/
      data ci  cs( 17) /  0.0000000000 0000000000 0000000910 01612d0/
      data ci  cs( 18) / -0.0000000000 0000000000 0000000003 07233d0/
      data ci  cs( 19) /  0.0000000000 0000000000 0000000000 00926d0/
c
      data nci, xsml /0, 0.0d0/
c
      if (nci.ne.0) go to 10
      nci = initds (cics, 19, 0.1*sngl(d1mach(3)))
      xsml = dsqrt (d1mach(3))
c
 10   if (x.le.0.0d0) call seteru ('17hdci     x is le 0', 17, 1, 2)
c
      if (x.gt.4.0d0) go to 20
      y = -1.0d0
      if (x.gt.xsml) y = (x*x-8.d0)*0.125d0
c
      dci = dlog(x) - 0.5d0 + dcsevl (y, cics, nci)
      return
c
 20   call d9sifg (x, f, g)
      sinx = dsin (x)
      call erroff
      dci = f*sinx - g*dcos(x)
c
      return
      end

