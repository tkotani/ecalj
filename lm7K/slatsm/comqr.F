      subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)
c
      integer i,j,l,n,en,ll,nm,igh,itn,its,low,lp1,enm1,ierr
      double precision hr(nm,n),hi(nm,n),wr(n),wi(n)
      double precision si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr,norm,tst1,tst2,
     xpythag
c
c     this subroutine is a translation of a unitary analogue of the
c     algol procedure  comlr, num. math. 12, 369-376(1968) by martin
c     and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
c     the unitary analogue substitutes the qr algorithm of francis
c     (comp. jour. 4, 332-345(1962)) for the lr algorithm.
c
c     this subroutine finds the eigenvalues of a complex
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  cbal.  if  cbal  has not been used,
c          set low=1, igh=n.
c
c        hr and hi contain the real and imaginary parts,
c          respectively, of the complex upper hessenberg matrix.
c          their lower triangles below the subdiagonal contain
c          information about the unitary transformations used in
c          the reduction by  corth, if performed.
c
c     on output
c
c        the upper hessenberg portions of hr and hi have been
c          destroyed.  therefore, they must be saved before
c          calling  comqr  if subsequent calculation of
c          eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  if an error
c          exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     calls cdiv for complex division.
c     calls csroot for complex square root.
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      if (low .eq. igh) go to 180
c     .......... create real subdiagonal elements ..........
      l = low + 1
c
      do 170 i = l, igh
        ll = min0(i+1,igh)
        if (hi(i,i-1) .eq. 0.0d0) go to 170
        norm = pythag(hr(i,i-1),hi(i,i-1))
        yr = hr(i,i-1) / norm
        yi = hi(i,i-1) / norm
        hr(i,i-1) = norm
        hi(i,i-1) = 0.0d0
c
        do 155 j = i, igh
          si = yr * hi(i,j) - yi * hr(i,j)
          hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
          hi(i,j) = si
  155   continue
c
        do 160 j = low, ll
          si = yr * hi(j,i) + yi * hr(j,i)
          hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
          hi(j,i) = si
  160   continue
c
  170 continue
c     .......... store roots isolated by cbal ..........
  180 do 200 i = 1, n
        if (i .ge. low .and. i .le. igh) go to 200
        wr(i) = hr(i,i)
        wi(i) = hi(i,i)
  200 continue
c
      en = igh
      tr = 0.0d0
      ti = 0.0d0
      itn = 30*n
c     .......... search for next eigenvalue ..........
  220 if (en .lt. low) go to 1001
      its = 0
      enm1 = en - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low d0 -- ..........
  240 do 260 ll = low, en
        l = en + low - ll
        if (l .eq. low) go to 300
        tst1 = dabs(hr(l-1,l-1)) + dabs(hi(l-1,l-1))
     x  + dabs(hr(l,l)) + dabs(hi(l,l))
        tst2 = tst1 + dabs(hr(l,l-1))
        if (tst2 .eq. tst1) go to 300
  260 continue
c     .......... form shift ..........
  300 if (l .eq. en) go to 660
      if (itn .eq. 0) go to 1000
      if (its .eq. 10 .or. its .eq. 20) go to 320
      sr = hr(en,en)
      si = hi(en,en)
      xr = hr(enm1,en) * hr(en,enm1)
      xi = hi(enm1,en) * hr(en,enm1)
      if (xr .eq. 0.0d0 .and. xi .eq. 0.0d0) go to 340
      yr = (hr(enm1,enm1) - sr) / 2.0d0
      yi = (hi(enm1,enm1) - si) / 2.0d0
      call csroot(yr**2-yi**2+xr,2.0d0*yr*yi+xi,zzr,zzi)
      if (yr * zzr + yi * zzi .ge. 0.0d0) go to 310
      zzr = -zzr
      zzi = -zzi
  310 call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
      sr = sr - xr
      si = si - xi
      go to 340
c     .......... form exceptional shift ..........
  320 sr = dabs(hr(en,enm1)) + dabs(hr(enm1,en-2))
      si = 0.0d0
c
  340 do 360 i = low, en
        hr(i,i) = hr(i,i) - sr
        hi(i,i) = hi(i,i) - si
  360 continue
c
      tr = tr + sr
      ti = ti + si
      its = its + 1
      itn = itn - 1
c     .......... reduce to triangle (rows) ..........
      lp1 = l + 1
c
      do 500 i = lp1, en
        sr = hr(i,i-1)
        hr(i,i-1) = 0.0d0
        norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
        xr = hr(i-1,i-1) / norm
        wr(i-1) = xr
        xi = hi(i-1,i-1) / norm
        wi(i-1) = xi
        hr(i-1,i-1) = norm
        hi(i-1,i-1) = 0.0d0
        hi(i,i-1) = sr / norm
c
        do 490 j = i, en
          yr = hr(i-1,j)
          yi = hi(i-1,j)
          zzr = hr(i,j)
          zzi = hi(i,j)
          hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
          hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
          hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
          hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
  490   continue
c
  500 continue
c
      si = hi(en,en)
      if (si .eq. 0.0d0) go to 540
      norm = pythag(hr(en,en),si)
      sr = hr(en,en) / norm
      si = si / norm
      hr(en,en) = norm
      hi(en,en) = 0.0d0
c     .......... inverse operation (columns) ..........
  540 do 600 j = lp1, en
        xr = wr(j-1)
        xi = wi(j-1)
c
        do 580 i = l, j
          yr = hr(i,j-1)
          yi = 0.0d0
          zzr = hr(i,j)
          zzi = hi(i,j)
          if (i .eq. j) go to 560
          yi = hi(i,j-1)
          hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
  560     hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
          hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
          hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
  580   continue
c
  600 continue
c
      if (si .eq. 0.0d0) go to 240
c
      do 630 i = l, en
        yr = hr(i,en)
        yi = hi(i,en)
        hr(i,en) = sr * yr - si * yi
        hi(i,en) = sr * yi + si * yr
  630 continue
c
      go to 240
c     .......... a root found ..........
  660 wr(en) = hr(en,en) + tr
      wi(en) = hi(en,en) + ti
      en = enm1
      go to 220
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end

