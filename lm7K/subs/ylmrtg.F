      subroutine ylmrtg(nlm,rotp,rmat)
C- Matrix to rotate cubic harmonics for a given rotation matrix
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nlm: (lmax+1)**2 for rmat
Ci   rotp: 3x3 rotation matrix
Co Outputs:
Co   rmat     matrix that rotates Y_L(r) to Y_L(rotp r)
Cr Remarks:
Cr   Y_L(rotp r)   = sum_M rmat(L,M) Y_M(r)  and also
Cr   Y_L(rotp-1 r) = sum_M Y_M(r) rmat(M,L)
Cr   rmat is block diagonal in l.
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters:
      integer nlm
      double precision rmat(nlm,nlm),rotp(9)
C Local parameters:
      integer lmx,mmx,ndim,ll,lmax,i,l,m,mmax,nm,offs,ierr,nlmi
      parameter ( lmx=8, mmx=2*lmx+1, ndim=(lmx+1)**2 )
      integer iwk(ndim)
      double precision p(3,mmx),rp(mmx,3),rr(mmx),
     .yp(ndim,ndim),ylm(mmx,ndim),y(ndim,ndim)
      integer nlmsav
      save nlmsav,y
      data nlmsav /0/
C ... A collection of random vectors, to make rmat = yrot(p) y^-1
      data p/ 0.020,-.025,-.118,
     .0.419,-.538,0.513,
     .0.245,-.717,-.600,
     .-.056,0.224,-.309,
     .-.034,-.180,0.207,
     .-.351,-.614,0.950,
     .-.782,-.134,-.308,
     .0.568,0.716,-.457,
     .-.528,-.927,-.562,
     .-.856,-.443,0.267,
     .-.111,0.794,0.598,
     .-.985,-.144,-.617,
     .0.678,0.400,-.617,
     .0.730,-.207,-.101,
     .0.540,-.137,-.773,
     .-.758,-.992,-.561,
     .0.321,-.363,-.988/

C     call tcn('ylmrtg')

C --- Initialization  ---
      lmax = ll(nlm)
      mmax = 2*lmax+1
      if (lmax .gt. lmx) call rx('increase lmx in ylmrtg')

C --- Set up and invert matrix ylm(p) ---
      if (nlm .gt. nlmsav) then
        call dpzero(y,ndim**2)

C   ... Normalize p, make ylm(p)
        do  8  i = 1, mmax
          call dscal(3,1/dsqrt(p(i,1)**2+p(2,i)**2+p(3,i)**2),p(1,i),1)
          rp(i,1) =p(1,i)
          rp(i,2) =p(2,i)
          rp(i,3) =p(3,i)
    8   continue
        call ropyln(mmax,rp,rp(1,2),rp(1,3),lmax,mmx,ylm,rr)

C   ... Generate matrix y(p)
        do  10  i = 1, mmax
          do  12  l = 0, lmax
            nm = 2*l+1
            offs = l*l
            if (i .le. nm) then
              do  14   m = 1, nm
                y(m+offs,i+offs) = ylm(i,m+offs)
   14         continue
            endif
   12     continue
   10   continue

C   ... Invert matrix y(p)
        y(1,1) = 1/y(1,1)
        do  16  l = 1, lmax
          offs = l**2
          nlmi = 2*l+1

C         call prmx('y',y(offs+1,offs+1),ndim,nlmi,nlmi)
          call dgetrf(nlmi,nlmi,y(offs+1,offs+1),ndim,iwk,ierr)
          if (ierr .ne. 0) call rx('ylmrtg cannot invert for rmat')
          call dgetri(nlmi,y(offs+1,offs+1),ndim,iwk,yp,ndim,ierr)
C         call prmx('y',y(offs+1,offs+1),ndim,nlmi,nlmi)
   16   continue

        nlmsav = nlm
C       print *, 'generated setup for ylmrtg'

      endif

C --- Set up matrix ylm(rotp*p) ---
      call dpzero(rmat,nlm**2)

C ... Make rp = rotp*p in with rp dimensioned (3,mmax)
C     call dgemm('N','N',3,mmax,3,1d0,rotp,3,p,3,0d0,rp,3)
C     call prmx('rp',rp,3,3,mmax)
C ... Make rp = rotp*p in with rp dimensioned (mmax,3)
      call dgemm('T','T',mmax,3,3,1d0,p,3,rotp,3,0d0,rp,mmx)
C     call prmx('rp',rp,mmx,mmax,3)
      call ropyln(mmax,rp,rp(1,2),rp(1,3),lmax,mmx,ylm,rr)

C ... Make matrix y(rp)
      do  20  i = 1, mmax
        do  22  l = 0, lmax
          nm = 2*l+1
          offs = l*l
          if (i .le. nm) then
            do  24   m = 1, nm
              yp(m+offs,i+offs) = ylm(i,m+offs)
   24       continue
          endif
   22   continue
   20 continue

C --- rmat = yrot * y^-1 ---
      rmat(1,1) = yp(1,1)*y(1,1)
      do  40  l = 1, lmax
        offs = l**2
        nlmi = 2*l+1

        call dgemm('N','N',nlmi,nlmi,nlmi,1d0,yp(offs+1,offs+1),ndim,
     .  y(offs+1,offs+1),ndim,0d0,rmat(offs+1,offs+1),nlm)

   40 continue
C      call prmx('ylmrtg: rmat',rmat,nlm,nlm,nlm)

C     call tcx('ylmrtg')
      end

C ... testing
C      subroutine fmain
C      implicit none
C      integer lmax,nlm,ng,i,wksize
C      parameter (lmax=8, nlm=(lmax+1)**2, wksize=30000)
C      double precision p(3*nlm),xmat(nlm**2),cy(289),
C     .  g(9,10),rold(nlm**2),rnew(nlm**2),diff,ag(3,10)
C      character*40 strn
C
C      integer w(wksize)
C      common /w/ w
C
C      call wkinit(wksize)
C      call sylmnc(cy,lmax+1)
C      call ylmrt0(lmax,nlm,xmat,p,cy)
C    1 print *, 'generator?'
C      read(*,'(a40)') strn
C      i = 0
C      call psymop(strn,w,g,ag,ng)
C      call awrit1(' group op %9:1d',' ',80,6,g)
C
C      call ylmrt1(lmax,nlm,g(1,1),rold,xmat,p,cy)
C      call ylmrtg(nlm,g,rnew)
C      call prmx('ok rotation matrix',rold,nlm,nlm,nlm)
C      call prmx('   rotation matrix',rnew,nlm,nlm,nlm)
C
C      do  10  i = 1, nlm**2
C        diff = rnew(i)-rold(i)
C        if (abs(diff) .gt. 1d-12) then
C          print *, 'error, i=',i,rnew(i),rold(i)
C          stop
C        endif
C   10 continue
C      print *, 'rotation matrix ok, max err=', sngl(diff)
C      goto 1
C      end

