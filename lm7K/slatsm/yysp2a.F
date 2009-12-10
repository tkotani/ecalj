      subroutine yysp2a(mode,i1,i2,j1,j2,apr,api,ldap,ija,offs,ar,ai,
     .  lda)
C- Unpacks a subblock of a complex sparse matrix into standard form
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit not used now
Ci          0  a is real
Ci          1  only the diagonal blocks of a are complex
Ci         >1  a is complex
Ci   i1,i2 :generate a(i,j) for row subblocks i = i1..i2
Ci   j1,j2 :generate a(i,j) for col subblocks j = j1..j2
Ci  apr,api:sparse matrix, stored in block form by rows.
Ci          ap consists of a vector of matrix subblocks:
Ci          apr(*,*,i) = real part of matrix subblock i
Ci          api(*,*,i) = imaginary part of matrix subblock i
Ci   ldap  :leading dimension of ap
Ci   ija   :column index packed array pointer data to array a
Ci         ija(1,*) follows essentially the same conventions
Ci         as for scalar packed arrays (see da2spr)
Ci         except that indices now refer to matrix subblocks.
Ci         ija(1,1)= n+2, where n = max(number of rows, number of cols)
Ci         ija(1,i), i = 1,..n+1 = points to first entry in a for row i
Ci         ija(1,i), i = n+2... column index element a(i).  Thus
Co                   for row i, k ranges from ija(i) <= k < ija(i+1) and
Co                   sum_j a_ij x_j -> sum_k a_(ija(2,k)) x_(ija(1,k))
Ci         ija(2,*)  pointers to the matrix subblocks blocks in a:
Ci         ija(2,i), i=1..n  pointers to blocks on the diagonal of a
Ci         ija(2,i), i=n+2.. pointers to elements of a, grouped by rows
Ci   offs  :offsets to first entries in matrix subblocks
Ci          Thus the dimension of row i = offs(i+1) - offs(i)
Ci          If a consists of scalar subblocks, offs(i) = i-1.
Ci   lda   :leading dimension of ar,ai
Co Outputs
Co   ar,ai :dense (normal) form of matrix subblock a(i1..i2,j1..j2)
Cr Remarks
Cr   The storage convention for packed form of matrices follows that of
Cr   Numerical Recipes (which see for further description), except that
Cr   here entries are themselves matrix subblocks, instead of scalars.
Cr   The matrix subblocks are storage in apr,api.
Cr   offs holds the information about the matrix subblock sizes.
Cr   ija(2,*) hold pointers the appropriate subblock in apr,api
C ----------------------------------------------------------------------
C     implicit none
      integer mode,i1,i2,j1,j2,ldap,lda,ija(2,*),offs(i2)
      double precision apr(ldap,ldap,1),api(ldap,ldap,1),
     .  ar(lda,1),ai(lda,1)
      integer i,j,ir,ic,pa,ofrb,ofcb,k,nra,nca,mode1
C     integer ofr,ofc

#if DEBUG
      character*10 fmt
      data fmt /'(9f12.6)'/
#endif

C     offsets shifting origin of rows and columns
C     ofr = 0
C     ofc = 0

C      do  6  j = offs(j1)-ofc+1, offs(j2+1)-ofc
C      do  6  i = offs(i1)-ofr+1, offs(i2+1)-ofr
      do  6  j = offs(j1)+1, offs(j2+1)
      do  6  i = offs(i1)+1, offs(i2+1)
        ar(i,j) = 0
        ai(i,j) = 0
    6 continue

C ... Fill diagonal elements a(ir,ir). Skip if missing.
      do  10  ir = i1, i2
C       offset to row for this subblock
C       ofrb = offs(ir) - ofr
        ofrb = offs(ir)
C       row dimension of this subblock
        nra = offs(ir+1) - offs(ir)
C       pointer to diagonal subblock in a
        pa  = ija(2,ir)
C   ... a(ir,ir) <- ap(pa).  Skip if missing.
        if (pa .ne. 0 .and. ir .ge. j1 .and. ir .le. j2) then
C         column offset to a for this subblock
C         ofcb = offs(ir) - ofc
          ofcb = offs(ir)
          do  14  j = 1, nra
          do  14  i = 1, nra
            ar(i+ofrb,j+ofcb) = apr(i,j,pa)
            ai(i+ofrb,j+ofcb) = api(i,j,pa)
   14     continue
        endif

        do  11  k = ija(1,ir), ija(1,ir+1)-1
C         column index to a
          ic  = ija(1,k)
          if (ic .lt. j1 .or. ic .gt. j2) goto 11
C         offset to column for this subblock
C         ofcb = offs(ic) - ofc
          ofcb = offs(ic)
C         column dimension of subblock
          nca = offs(ic+1) - offs(ic)
C         pointer to subblock in packed matrix
          pa = ija(2,k)
C         a(ir,ija(k)) <- ap(pa)
          if (mode .lt. 2) then
            do  15  j = 1, nca
            do  15  i = 1, nra
              ar(i+ofrb,j+ofcb) = apr(i,j,pa)
   15       continue
          else
            do  16  j = 1, nca
            do  16  i = 1, nra
              ar(i+ofrb,j+ofcb) = apr(i,j,pa)
              ai(i+ofrb,j+ofcb) = api(i,j,pa)
   16       continue
          endif
   11   continue
   10 continue

c     call yprm(.false.,'a',2,6,fmt,ar,lda,offs(i2+1),lda,offs(j2+1))
      end

