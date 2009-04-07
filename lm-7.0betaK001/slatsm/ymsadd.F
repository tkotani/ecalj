      subroutine ymsadd(nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd,
     .  alpha,beta,src,offsi,dest,offdi)
C- Add a subblock block of a complex matrix to a destination
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlma : row subblock size
Ci   nlmb : col subblock size
Ci   ndas : row dimension of source matrix
Ci   ndad : row dimension of destination matrix
Ci   ofas : offset to first row of source matrix
Ci   ofbs : offset to first col of source matrix
Ci   ofad : offset to first row of destination matrix
Ci   ofbd : offset to first col of destination matrix
Ci   alpha: complex scalar alpha scaling src; see Remarks
Ci   beta : complex scalar beta scaling dest; see Remarks
Ci   src  : source matrix
Ci   offsi: separation between real,imaginary parts of src
Ci   offdi: separation between real,imaginary parts of dest
Co Outputs
Co   dest : is overwritten by alpha * src + beta * dest
Cr Remarks
Cr   ymsadd adds a subblock alpha * src into a subblock of dest:
Cr   dest is overwritten by
Cr       dest <- alpha * src + beta * dest
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlma,nlmb,ndas,ndad,ofas,ofbs,ofad,ofbd,offsi,offdi
      double precision src(ndas,1),dest(ndad,1),alpha(2),beta(2)
C ... Local parameters
      integer ia,ib,ofsi,ofdi

      if (nlma .le. 0 .or. nlmb .le. 0) return

      ofsi = ofas + offsi
      ofdi = ofad + offdi

      if (beta(1) .eq. 0 .and. beta(2) .eq. 0) then
C       dest = src
        if (alpha(1) .eq. 1 .and. alpha(2) .eq. 0) then
          do  10  ib = 1, nlmb
          do  10  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = src(ia+ofas,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = src(ia+ofsi,ib+ofbs)
   10     continue
C       dest = -src
        elseif (alpha(1) .eq. -1 .and. alpha(2) .eq. 0) then
          do  11  ib = 1, nlmb
          do  11  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = -src(ia+ofas,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = -src(ia+ofsi,ib+ofbs)
   11     continue
C       dest = alpha * src
        else
          do  12  ib = 1, nlmb
          do  12  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = alpha(1)*src(ia+ofas,ib+ofbs) -
     .                              alpha(2)*src(ia+ofsi,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = alpha(1)*src(ia+ofsi,ib+ofbs) +
     .                              alpha(2)*src(ia+ofas,ib+ofbs)
   12     continue
        endif
      elseif (beta(1) .eq. 1 .and. beta(2) .eq. 0) then
C       dest = dest + src
        if (alpha(1) .eq. 1 .and. alpha(2) .eq. 0) then
          do  20  ib = 1, nlmb
          do  20  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = dest(ia+ofad,ib+ofbd) +
     .                              src(ia+ofas,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = dest(ia+ofdi,ib+ofbd) +
     .                              src(ia+ofsi,ib+ofbs)
   20     continue
C       dest = dest - src
        elseif (alpha(1) .eq. -1 .and. alpha(2) .eq. 0) then
          do  21  ib = 1, nlmb
          do  21  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = dest(ia+ofad,ib+ofbd) -
     .                              src(ia+ofas,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = dest(ia+ofdi,ib+ofbd) -
     .                              src(ia+ofsi,ib+ofbs)
   21     continue
C       dest = dest + alpha * src
        else
          do  22  ib = 1, nlmb
          do  22  ia = 1, nlma
            dest(ia+ofad,ib+ofbd) = dest(ia+ofad,ib+ofbd) +
     .                              alpha(1)*src(ia+ofas,ib+ofbs) -
     .                              alpha(2)*src(ia+ofsi,ib+ofbs)
            dest(ia+ofdi,ib+ofbd) = dest(ia+ofdi,ib+ofbd) +
     .                              alpha(1)*src(ia+ofsi,ib+ofbs) +
     .                              alpha(2)*src(ia+ofas,ib+ofbs)
   22     continue
        endif
C     dest = beta * dest + alpha * src
      else
        do  30  ib = 1, nlmb
        do  30  ia = 1, nlma
          dest(ia+ofad,ib+ofbd) = beta(1)*dest(ia+ofad,ib+ofbd) -
     .                            beta(2)*dest(ia+ofdi,ib+ofbd) +
     .                            alpha(1)*src(ia+ofas,ib+ofbs) -
     .                            alpha(2)*src(ia+ofsi,ib+ofbs)
          dest(ia+ofdi,ib+ofbd) = beta(1)*dest(ia+ofdi,ib+ofbd) +
     .                            beta(2)*dest(ia+ofad,ib+ofbd) +
     .                            alpha(1)*src(ia+ofsi,ib+ofbs) +
     .                            alpha(2)*src(ia+ofas,ib+ofbs)
   30   continue
      endif
      end
#if TEST
      subroutine fmain
C      implicit none
      integer lda,ldb,nr,nc,ir,ic,i1mach,ld1,ld2,ldr,ldr1,oi,oi1,kcplxi,
     .  kcplxf,lerr,so,lda1,ldb1
      parameter (lda=5,ldb=6,nr=4,nc=3,lda1=7,ldb1=8)
      double precision alpha(2), beta(2)
      double precision yy(lda,ldb,2),y(lda,2,ldb),z(2,lda,ldb)
      double precision yy0(lda,ldb,2),diff
      double precision yy1(lda1,ldb1,2),yy10(lda1,ldb1,2)
      character*8 fm
      equivalence (yy,y)
      equivalence (yy,z)

      so = i1mach(2)
      fm = '(8f8.2)'

      lerr = 0

      do  10  ir = 1, lda
      do  10  ic = 1, ldb
        yy(ir,ic,1) = 100*ir + ic
        yy(ir,ic,2) = -(100*ir + ic)

        yy0(ir,ic,1) = 100*ir + ic
        yy0(ir,ic,2) = -(100*ir + ic)

   10 continue

      do  20  ir = 1, lda1
      do  20  ic = 1, ldb1
        yy1(ir,ic,1) = 100*ir + ic
        yy1(ir,ic,2) = -(100*ir + ic)

        yy10(ir,ic,1) = 100*ir + ic
        yy10(ir,ic,2) = -(100*ir + ic)

   20 continue

      do  30  kcplxi = 0, 2, 2
      call cplxdm(kcplxi,lda,ldb,ld1,ld2,ldr,oi)
      call cplxdm(kcplxi,lda1,ldb1,ld1,ld2,ldr1,oi1)
      call ztoyy(yy,lda,ldb,lda,ldb,0,kcplxi)
      call ztoyy(yy1,lda1,ldb1,lda1,ldb1,0,kcplxi)
      call ztoyy(yy10,lda1,ldb1,lda1,ldb1,0,kcplxi)

C ... straight copy, with an offset
      alpha(1) = 1
      alpha(2) = 0
      beta(1) = 0
      beta(2) = 0
      call dvset(yy1,1,lda1*ldb1*2,-1d0)
      call ymscop(4,2,3,ldr,ldr1,2,1,1,2,yy,oi,yy1,oi1)
      call ywrm(0,'real copy yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
C     call ymsadd(2,3,ldr,ldr1,2,1,1,2,alpha,beta,yy,oi,yy1,oi1)
      call ymscop(0,2,3,ldr,ldr1,2,1,1,2,yy,oi,yy1,oi1)
      call ywrm(0,'complex copy yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)

C ... add with alpha=1-i, beta=1
      alpha(1) = 1
      alpha(2) = -1
      beta(1) = 1
      beta(2) = 0
      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
      call ymsadd(3,2,ldr,ldr1,0,0,0,0,alpha,beta,yy,oi,yy1,oi1)
      call ywrm(0,'yy1+(1-i)yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)

C ... with alpha=-1, beta=0
      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
      alpha(1) = -1
      alpha(2) = 0
      beta(1) = 0
      beta(2) = 0
      call ymsadd(3,2,ldr,ldr1,0,0,0,0,alpha,beta,yy,oi,yy1,oi1)
      call ywrm(0,'- yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)

C ... with alpha=-1, beta=1
      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
      alpha(1) = -1
      alpha(2) = 0
      beta(1) = 1
      beta(2) = 0
      call ymsadd(2,3,ldr,ldr1,2,1,1,2,alpha,beta,yy,oi,yy1,oi1)
      call ywrm(0,'yy1 - yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
C ... ditto, for real matrix using ymscop mode=5
      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
      call dscal(lda*ldb*2,-1d0,yy,1)
      call ymscop(5,2,3,ldr,ldr1,2,1,1,2,yy,oi,yy1,oi1)
      call ywrm(0,'real yy1 - yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)
      call dscal(lda*ldb*2,-1d0,yy,1)

C ... with beta=0
      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
      alpha(1) = 0
      alpha(2) = 1
      beta(1) = 0
      beta(2) = 0
      call ymsadd(3,2,ldr,ldr1,0,0,0,0,alpha,beta,yy,oi,yy1,oi1)
      call ywrm(0,'i * yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)

C ... add with alpha=beta=1
      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
      alpha(1) = 1
      alpha(2) = 0
      beta(1) = 1
      beta(2) = 0
      call ymsadd(3,2,ldr,ldr1,0,0,0,0,alpha,beta,yy,oi,yy1,oi1)
      call ywrm(0,'y + yy',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)

C ... add with general alpha,beta
      alpha(1) = 1
      alpha(2) = -1
      beta(1) = -1d0/2
      beta(2) = 1d0/2
      call dcopy(lda1*ldb1*2,yy10,1,yy1,1)
      call ymsadd(3,2,ldr,ldr1,0,0,0,0,alpha,beta,yy,oi,yy1,oi1)
      call ywrm(0,'(1-i)(yy-yy1/2)',kcplxi+2,so,fm,yy1,oi1,lda1,nr,ldb1)

   30 continue

      end
#endif

