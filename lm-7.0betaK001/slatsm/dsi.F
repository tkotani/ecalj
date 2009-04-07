      double precision function dsi (x)
c december 1980 edition, w. fullerton, bell labs.
      double precision x, sics(18), pi2, xsml, absx, f, g, cosx,
     1  dcsevl, d1mach, dcos, dsin, dsqrt
      external d1mach, dcos, dcsevl, dsin, dsqrt, initds
c
c series for si   on the interval  0.00000e+00 to  1.60000e+01
c                                        with weighted error   8.58e-32
c                                         log weighted error  31.07
c                               significant figures required  30.53
c                                    decimal places required  31.69
c
      data si  cs(  1) / -0.1315646598 1848419289 0427517300 0457d0/
      data si  cs(  2) / -0.2776578526 9736018920 4828766015 7299d0/
      data si  cs(  3) /  0.0354414054 8666591797 4913546471 0086d0/
      data si  cs(  4) / -0.0025631631 4479339776 5875278836 1530d0/
      data si  cs(  5) /  0.0001162365 3904970092 8126492148 2985d0/
      data si  cs(  6) / -0.0000035904 3272416060 4267000434 7148d0/
      data si  cs(  7) /  0.0000000802 3421237057 1016230865 2976d0/
      data si  cs(  8) / -0.0000000013 5629976925 4025064993 1846d0/
      data si  cs(  9) /  0.0000000000 1794407215 9973677556 7759d0/
      data si  cs( 10) / -0.0000000000 0019083873 4308714549 0737d0/
      data si  cs( 11) /  0.0000000000 0000166699 8958682433 0853d0/
      data si  cs( 12) / -0.0000000000 0000001217 3098836850 3042d0/
      data si  cs( 13) /  0.0000000000 0000000007 5418186699 3865d0/
      data si  cs( 14) / -0.0000000000 0000000000 0401417884 2446d0/
      data si  cs( 15) /  0.0000000000 0000000000 0001855369 0716d0/
      data si  cs( 16) / -0.0000000000 0000000000 0000007516 6966d0/
      data si  cs( 17) /  0.0000000000 0000000000 0000000026 9113d0/
      data si  cs( 18) / -0.0000000000 0000000000 0000000000 0858d0/
c
      data pi2 / 1.5707963267 9489661923 1321691639 75 d0 /
      data nsi, xsml /0, 0.0d0/
c
      if (nsi.ne.0) go to 10
      nsi = initds (sics, 18, 0.1*sngl(d1mach(3)))
      xsml = dsqrt(d1mach(3))
c
 10   absx = dabs(x)
      if (absx.gt.4.0d0) go to 20
      dsi = x
      if (absx.lt.xsml) return
c
      dsi = x*(0.75d0 + dcsevl ((x*x-8.0d0)*.125d0, sics, nsi))
      return
c
 20   call d9sifg (absx, f, g)
      cosx = dcos (absx)
      call erroff
      dsi = pi2 - f*cosx - g*dsin(x)
      if (x.lt.0.0d0) dsi = -dsi
c
      return
      end

