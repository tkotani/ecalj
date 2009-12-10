      subroutine shopol(nl,nclass,nsp,sspec,ics,pp)
C- Find poles for third order potential functions
C-----------------------------------------------------------------------
Ci Inputs
Ci   nl,nbas,nsp,isp,ipc,pp (in gamma representation)
Co Outputs
Cr Remarks
C-----------------------------------------------------------------------
C     implicit none
C Passed parameters
      double precision sspec(1)
      integer nl,nsp,nclass,ics(nclass)
      double precision pp(6,nl,nsp,*)
C Local parameters
      integer l,ic,isp,ipr,i,lmxa,is
      double precision oalpha,pgam,enu,C,gam,del,pa,alp,xx
      double complex zz,phold(3),zmc,guess,update

      call getpr(ipr)
      ipr = 55
      if (ipr .ge. 55) then
        if (nsp .eq. 1) then
          print '(''  ib   l   Poles P_l ...'')'
        else
          print '(''  ib spin  Poles P_l(ib) ...'')'
        endif
      endif

      do  4  isp = 1, nsp
        do  3  ic = 1, nclass
          is = ics(ic)
          call upack('spec lmxa',sspec,is,lmxa,0,0,0)
          do  2  l = 0, nl-1
              if (l .gt. lmxa) goto 2
              enu = pp(1,l+1,isp,ic)
              C   = pp(2,l+1,isp,ic)
              del = pp(3,l+1,isp,ic)**2
              pa  = pp(4,l+1,isp,ic)
              gam = pp(5,l+1,isp,ic)
              alp = pp(6,l+1,isp,ic)
              pgam = pa - oalpha(enu,C,del,alp,gam)**2

              guess = C
C             guess = (-1.045d0, .002d0)
              i = 1
   10         continue
              zz = guess + pgam * (guess - enu)**3
              zmc = zz - C
              update = zmc/(1+pgam*3*(guess-enu)**2)
              guess = guess - update
              if (cdabs(update) .gt. 1d-10) goto 10
              phold(i) = guess
              if (i .eq. 1) then
                i = 2
                xx = (C-enu)**2*pgam
                guess = dcmplx(enu-(C-enu)*(1-xx)/2,
     .                         sqrt(1/pgam)*(1+3*xx/8))
                goto 10
              endif
              if (i .eq. 2) then
                i = 3
                guess = dconjg(guess)
                goto 10
              endif
              if (ipr .ge. 55 .and. nsp .eq. 1) then
                print 333,  ic, l, (phold(i), i=1,3)
              elseif (ipr .ge. 55) then
                print 334,  ic, isp, l, (phold(i), i=1,3)
              endif
    2     continue
  333       format(2i4,10(1x,2f10.5,2x))
  334       format(3i4,10(1x,2f10.5,1x))
    3   continue
    4 continue
      end

