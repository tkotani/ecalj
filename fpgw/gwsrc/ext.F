C --- Simple FORTRAN dynamic memory allocation ---
C Memory is allocated from a single heap.  The heap is declared as
C an integer array in common block W by the main program.  Its size
C is fixed by the size of W as declared by the main program.  It
C is possible to allocated memory from the top of the heap, and to
C free memory already allocated from the top of the heap.  It is
C not possible to free memory allocated from the middle of the heap.
C
C To use these routines, first call WKINIT with NSIZE = size of W.
C
C Character, integer, real, double precision, complex or double
C complex arrays can then be allocated by calling DEFCH, DEFI,
C DEFR, DEFDR, DEFC or DEFDC with an index IPNAME and LENGTH the
C desired length.  IPNAME is returned as the effective offset of
C the integer array W.  If it is desired that the array be
C initialized to zero, pass LENGTH as the negative of the desired
C length.
C
C Allocated memory can be freed by calling RLSE(IPNAME) where
C IPNAME must one of the indices from a previous allocation.  All
C blocks of memory allocated after IPNAME are also released.
C
C WKPRNT turns on or off a trace that prints a message each time
C        memory is allocated or released.  Also when turned on,
C        runs through links each time RLSE is called.
C DEFASK returns the amount of memory remaining in the pool
C WKSAV  sets an internal switch so that a heap allocation request
C        returns the the negative of length of the request in IPNAME
C        when not enough memory is available, rather than aborting.
C        In this case, no memory is allocated and no other internal
C        variables are changed.
C WKINFO prints information about the arrays allocated
C WKCHK  runs through links to see if any array was overwritten
C        call wkchk(string), where string terminated by '$'
C ----------------------------------------------------------------
      SUBROUTINE WKINIT(NSIZE)
C- Initialize conditions for the 'common pool' data W
ca !!!!!!!!!!!
ca    INTEGER W(2)
      INTEGER W(1)
      CHARACTER*1 STRING(60),STR(60)
      logical lopt,lerr
      SAVE
      COMMON /W/ W

      IRND(I) = (I+499)/1000

C ----- DEFINE STORAGE SIZE ------
C  START OF FIRST ARRAY AND MAX NUMBER TO BE DEFINED:
      IP0 = 5
      NDEFMX = 100
      LIMIT = NSIZE
      IPMAX = 0
      NDFDMX = 0
      NDFSUM = 0
      JPR = 0
      IPFREE = 5
      lerr = .false.
      WRITE(*,*) 'WKINIT:  size=',IRND(NSIZE),'K'
      WRITE(*,*) ' '
      RETURN
      ENTRY WKPRNT(JPRINT)
C- Set debug switch for heap management
      if (jprint .eq. 2) then
        jpr = 1-jpr
      else
        JPR = JPRINT
      endif
      RETURN
      entry wksav(lopt)
      lerr = lopt
      return
C ------ SUBROUTINES TO DEFINE ARRAYS OF VARIOUS TYPES -----
      ENTRY DEFCH(IPNAME,LENG)
C- Allocate character array
      LENGTH = (LENG+3)/4
      GOTO 10
      ENTRY DEFI(IPNAME,LENG)
C- Allocate integer array
      LENGTH = LENG
      JOPT = 1
      GOTO 10
      ENTRY DEFR(IPNAME,LENG)
C- Allocate single precision real array
      LENGTH = LENG*I1MACH(17)
      JOPT = 2
      GOTO 10
      ENTRY DEFRR(IPNAME,LENG)
      ENTRY DEFDR(IPNAME,LENG)
C- Allocate double precision real array
      LENGTH = LENG*I1MACH(18)
      JOPT = 3
      GOTO 10
      ENTRY DEFC(IPNAME,LENG)
C- Allocate single precision complex array
      LENGTH = LENG*2*I1MACH(17)
      JOPT = 4
      GOTO 10
      ENTRY DEFCC(IPNAME,LENG)
      ENTRY DEFDC(IPNAME,LENG)
C- Allocate double precision complex array
      LENGTH = LENG*2*I1MACH(18)
      JOPT = 5
   10 IOPT = 0
      IF (LENGTH .LT. 0) THEN
        IOPT = 1
        LENGTH = -LENGTH
      ENDIF
      IF (LENGTH .EQ. 0) LENGTH = 1
      IMOD = 1
      GOTO 83
   84 IPNAME = IPFREE
      if (lerr .and. ipfree+length+2 .gt. limit) then
        ipname = -LENGTH
        if (jpr .gt. 0)  print *,
     .    'ALLOC: heap storage exceeded; returning -LENGTH=',-LENGTH
        return
      endif
      IPFREE = IPFREE + LENGTH + 1
      IPFREE = 4*((IPFREE+2)/4)+1
      IPMAX = MAX0(IPMAX,IPFREE)
      W(IPNAME-1)=IPFREE
      NDEFD = NDEFD + 1
      NDFDMX = MAX0(NDFDMX,NDEFD)
      NDFSUM = NDFSUM + LENGTH
      IF (JPR .GT. 0) WRITE(*,100) NDEFD,LENG,LENGTH,IPNAME,IPFREE-1
  100 FORMAT(' define array',I4,':   els=',I8,'   length=',I8,',',
     .   I8,'  to',I8)
      IF (IPFREE .LE. LIMIT) THEN
c takao apl9 2008
c        IF (IOPT .NE. 0) GOTO (201,202,203,204,205) JOPT
        IF (IOPT .NE. 0) then
          if(JOPT==1)  GOTO 201
          if(JOPT==2)  GOTO 202
          if(JOPT==3)  GOTO 203
          if(JOPT==4)  GOTO 204
          if(JOPT==5)  GOTO 205
        ENDIF

        RETURN
  201   CALL IINIT(W(IPNAME),-LENG)
        RETURN
  202   CALL AINIT(W(IPNAME),-LENG)
        RETURN
  203   CALL DINIT(W(IPNAME),-LENG)
        RETURN
  204   CALL CINIT(W(IPNAME),-LENG)
        RETURN
  205   CALL ZINIT(W(IPNAME),-LENG)
        RETURN
      ENDIF
      WRITE(*,101) IPFREE
  101 FORMAT(' ALLOC: WORKSPACE OVERFLOW, NEED AT LEAST',I8)
      STOP
C- Release data up to pointer
      ENTRY RLSE(IPNAME)
      IF (IPNAME .GT. LIMIT) STOP 'RLSE: release pointer exceeds limit'
      IF (IPNAME .LT. 3) STOP 'RLSE: release pointer less than 3'
      IF (JPR .eq. 0) goto 82
      imod = 3
      goto 83
   87 WRITE(*,*) 'RLSE from: ',IPNAME
   82 IPFREE = IPNAME
      return
      ENTRY DEFASK(LREST)
C- Return number of words left in common pool
      LREST = LIMIT - IPFREE - 2
      IF (JPR .GT. 0) WRITE(*,*) 'SPACE LEFT=',LREST,'  SINGLE WORDS'
      RETURN
      ENTRY WKINFO()
C- Output workspace information
      IMOD = 2
      GOTO 83
  81  WRITE(*,601) IRND(LIMIT),IRND(NDFSUM),IRND(IPMAX-1),
     .             IRND(IPFREE-1),NDFDMX,NDEFD
  601 FORMAT(
     .  /'  total workspace size =',I5,' K',
     .  /'  total space allocated=',I5,' K',
     .  /'  workspace used:    max',I5,' K   now',I4,' K',
     .  /'  arrays defined:    max',I7,  '   now',I6)
      IF (IPFREE .EQ. IP0) RETURN
      if (jpr .gt. 0) WRITE(*,602)
  602 FORMAT(/'  array',6X,'begin',7X,'end',7X,'length')
      IPX = IP0
      DO  30  I = 1, NDEFMX
        IPY = W(IPX-1)
        IF (IPX .EQ. IPFREE) RETURN
        IF (IPY .LT. IP0 .OR. IPY .GT. LIMIT) WRITE(*,*) '   . . . . . '
        IF (IPY .LT. IP0 .OR. IPY .GT. LIMIT) RETURN
        IF (JPR .GT. 0) WRITE(*,603) I,IPX,IPY,IPY-IPX
  603   FORMAT(4(I6,5X))
        IPX = IPY
   30 CONTINUE
      RETURN
      ENTRY WKCHK(STRING)
C- Run through links to see if any dynamic array was overwritten
      IMOD = 0
      DO  88  I = 1, 60
        STR(I)=STRING(I)
        NSTR = I-1
        IF (STRING(I) .EQ. '$') GOTO 89
   88 continue
   89 WRITE(*,*) 'WKCHK: ',(STR(I),I = 1,NSTR)
   83 NDEFD = 0
      IPX = IP0
      IPPLOC = -999
      DO  35  I = 1, NDEFMX
        IF (IPX .LT. IP0 .OR. IPX .GT. LIMIT) THEN
          WRITE(*,888) NDEFD,IPX,IPPLOC
  888     FORMAT(' ALLOC: LINK DESTROYED AT START OF ARRAY',I3,
     .     ',  PTR=',I8,' AT',I8)
          STOP
        ENDIF
        IF (IPX .EQ. IPFREE) GOTO 86
        NDEFD = NDEFD + 1
        IPPLOC = IPX - 1
        IPX = W(IPPLOC)
   35 continue
  86  CONTINUE
c apr9 2008
c      GOTO (84,81,87), imod
      if(imod==1) goto 84
      if(imod==2) goto 81
      if(imod==3) goto 87

      WRITE(*,360) NDEFD,IPFREE-1
  360 FORMAT('     LINKS OK   NDEFD=',I3,'   SPACE USED=',I7)
      RETURN
      END
      SUBROUTINE CINIT(ARRAY,LENG)
C- Initializes complex array to zero
      INTEGER LENG
      REAL ARRAY(2*LENG)
      CALL AINIT(ARRAY,LENG+LENG)
      RETURN
      END
      SUBROUTINE ZINIT(ARRAY,LENG)
C- Initializes complex array to zero
      INTEGER LENG
      DOUBLE PRECISION ARRAY(2*LENG)
      CALL DINIT(ARRAY,LENG+LENG)
      RETURN
      END
      SUBROUTINE IINIT(ARRAY,LENG)
C- Initializes integer array to zero
      INTEGER LENG
      INTEGER ARRAY(LENG)
      DO 10 I=1,LENG
        ARRAY(I)=0
   10 CONTINUE
      RETURN
      END
      SUBROUTINE AINIT(ARRAY,LENG)
C- Initializes real array to zero
      INTEGER LENG
      REAL ARRAY(LENG)
      DO 10 I=1,LENG
        ARRAY(I)=0.
   10 CONTINUE
      RETURN
      END
      SUBROUTINE DINIT(ARRAY,LENG)
C- Initializes double precision array to zero
      INTEGER LENG
      DOUBLE PRECISION ARRAY(LENG)
      DO 10 I=1,LENG
        ARRAY(I)=0.D0
   10 CONTINUE
      RETURN
      END
      FUNCTION RVAL(ARRAY,INDEX)
C- Returns the real value of ARRAY(INDEX)
      REAL ARRAY(INDEX)
      RVAL=ARRAY(INDEX)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DRVAL(ARRAY,INDEX)
C- Returns the double precision value of ARRAY(INDEX)
      DOUBLE PRECISION ARRAY(INDEX)
      DRVAL=ARRAY(INDEX)
      RETURN
      END
      INTEGER FUNCTION IVAL(ARRAY,INDEX)
C- Returns the integer value of ARRAY(INDEX)
      INTEGER ARRAY(INDEX)
      IVAL=ARRAY(INDEX)
      RETURN
      END
      COMPLEX FUNCTION CVAL(ARRAY,INDEX)
C- Returns the complex value of ARRAY(INDEX)
      COMPLEX ARRAY(INDEX)
      CVAL=ARRAY(INDEX)
      RETURN
      END



      SUBROUTINE DMADD(A,NCA,NRA,SCALEA,B,NCB,NRB,SCALEB,C,NCC,NRC,N,M)
C- general matrix addition
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   b,ncb,nrb is the right matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   c,ncc,nrc is the result matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   n,m: the number of rows and columns, respectively, to calculate
Co Outputs:
Co   result matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix linear combination routine,
Cr   adding a subblock of matrix a to a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop.
Cr
Cr   Destination matrix c can coincide with either a or b, provided that
Cr   the transpose of the coincident matrix is not taken.
Cr   Example: Add 3-by-2 block of (transpose of a - .5*b) into c
Cr     call dmadd(a,1,na,1.d0,b,nb,0,-.5d0,c,nc,1,3,2)
Cr     OLD call dmadd(a,na,1,b,nb,0,-.5d0,c,nc,3,2)
C ----------------------------------------------------------------
C
      INTEGER NCA,NRA,NCB,NRB,NCC,NRC,N,M
      DOUBLE PRECISION A(0:1), B(0:1), C(0:1), SCALEA, SCALEB
      INTEGER I,J,IA,IB,IC

      DO  200  I = N-1, 0, -1
        IA = I*NRA+M*NCA
        IB = I*NRB+M*NCB
        IC = I*NRC+M*NCC
      DO  200  J = M-1, 0, -1
        IA = IA-NCA
        IB = IB-NCB
        IC = IC-NCC
        C(IC) = A(IA)*SCALEA + B(IB)*SCALEB
  200 CONTINUE
      RETURN
      END
      subroutine dpmpy(a,b,nscb,nsrb,c,nscc,nsrc,nr,nc,l)
C- matrix multiplication, (packed) (normal) -> (normal)
C ----------------------------------------------------------------
Ci Inputs:
Ci   a is the left matrix (packed)
Ci   b,nscb,nsrb is the right matrix and respectively the spacing
Ci      between column elements and row elements.
Ci   c,nscc,nsrc is the product matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   nr,nc: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed nr*nc*l times.  No attempt is made to optimize
Cr   the outer loops, executed nr*nc times.
Cr     Examples: product of (nr,l) subblock of a into (l,nc) subblock of
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,nr,nc,l)
Cr     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
Cr     To generate the tranpose of that product, use:
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,nr,nc,l)
C ----------------------------------------------------------------
C Passed Parameters
      integer nscb,nsrb,nscc,nsrc,nr,nc,l
      double precision a(0:*), b(0:*), c(0:*)
C Local parameters
      double precision sum
      integer i,j,k,offa,offb

      do  20  i = 0, nr-1
      do  20  j = 0, nc-1
        sum = 0
        offa = (i*(i+1))/2
        offb = nscb*j
        do  21  k = 0, i-1
          sum = sum + a(offa)*b(offb)
          offa = offa + 1
          offb = offb + nsrb
   21   continue
        do  22  k = i, l-1
          sum = sum + a(offa)*b(offb)
          offa = offa + k+1
          offb = offb + nsrb
   22   continue
        c(i*nsrc+j*nscc) = sum
   20 continue
      end
      REAL FUNCTION R1MACH(I)
C
C  SINGLE-PRECISION MACHINE CONSTANTS
C
C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  R1MACH(5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C
C  WHERE POSSIBLE, OCTAL OR HEXADECIMAL CONSTANTS HAVE BEEN USED
C  TO SPECIFY THE CONSTANTS EXACTLY WHICH HAS IN SOME CASES
C  REQUIRED THE USE OF EQUIVALENT INTEGER ARRAYS.
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      REAL RMACH(5)
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
C#ifdef IEEE | RIEEE
      DATA SMALL(1) /     8388608 /
      DATA LARGE(1) /  2139095039 /
      DATA RIGHT(1) /   864026624 /
      DATA DIVER(1) /   872415232 /
      DATA LOG10(1) /  1050288283 /
C#endif
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1) /    1048576 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  990904320 /
C      DATA DIVER(1) / 1007681536 /
C      DATA LOG10(1) / 1091781651 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA RMACH(1) / Z400800000 /
C      DATA RMACH(2) / Z5FFFFFFFF /
C      DATA RMACH(3) / Z4E9800000 /
C      DATA RMACH(4) / Z4EA800000 /
C      DATA RMACH(5) / Z500E730E8 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.
C
C      DATA RMACH(1) / O1771000000000000 /
C      DATA RMACH(2) / O0777777777777777 /
C      DATA RMACH(3) / O1311000000000000 /
C      DATA RMACH(4) / O1301000000000000 /
C      DATA RMACH(5) / O1157163034761675 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C      DATA RMACH(1) / 00014000000000000000B /
C      DATA RMACH(2) / 37767777777777777777B /
C      DATA RMACH(3) / 16404000000000000000B /
C      DATA RMACH(4) / 16414000000000000000B /
C      DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA RMACH(1) / '00800000'X /
C      DATA RMACH(2) / '7FFFFFFF'X /
C      DATA RMACH(3) / '34800000'X /
C      DATA RMACH(4) / '35000000'X /
C      DATA RMACH(5) / '3F9A209B'X /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C#ifdefC CRAY
C      DATA RMACH(1) / 200034000000000000000B /
C      DATA RMACH(2) / 577767777777777777776B /
C      DATA RMACH(3) / 377224000000000000000B /
C      DATA RMACH(4) / 377234000000000000000B /
C      DATA RMACH(5) / 377774642023241175720B /
C#endif
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC RMACH(5)
C
C      DATA SMALL/20K,0/,LARGE/77777K,177777K/
C      DATA RIGHT/35420K,0/,DIVER/36020K,0/
C      DATA LOG10/40423K,42023K/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '00000177 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000353 /
C      DATA LOG10(1),LOG10(2) / '23210115, '00000377 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA RMACH(1) / O402400000000 /
C      DATA RMACH(2) / O376777777777 /
C      DATA RMACH(3) / O714400000000 /
C      DATA RMACH(4) / O716400000000 /
C      DATA RMACH(5) / O776464202324 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C#ifdefC IBM_VM | IBM370 | IBM3080 | IBM3090
C      DATA RMACH(1) / Z00100000 /
C      DATA RMACH(2) / Z7FFFFFFF /
C      DATA RMACH(3) / Z3B100000 /
C      DATA RMACH(4) / Z3C100000 /
C      DATA RMACH(5) / Z41134413 /
C#endif
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA RMACH(1) / Z'00100000' /
C      DATA RMACH(2) / Z'7EFFFFFF' /
C      DATA RMACH(3) / Z'3B100000' /
C      DATA RMACH(4) / Z'3C100000' /
C      DATA RMACH(5) / Z'41134413' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C
C      DATA RMACH(1) / "000400000000 /
C      DATA RMACH(2) / "377777777777 /
C      DATA RMACH(3) / "146400000000 /
C      DATA RMACH(4) / "147400000000 /
C      DATA RMACH(5) / "177464202324 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1) /    8388608 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  880803840 /
C      DATA DIVER(1) /  889192448 /
C      DATA LOG10(1) / 1067065499 /
C
C      DATA RMACH(1) / O00040000000 /
C      DATA RMACH(2) / O17777777777 /
C      DATA RMACH(3) / O06440000000 /
C      DATA RMACH(4) / O06500000000 /
C      DATA RMACH(5) / O07746420233 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /   128,     0 /
C      DATA LARGE(1),LARGE(2) / 32767,    -1 /
C      DATA RIGHT(1),RIGHT(2) / 13440,     0 /
C      DATA DIVER(1),DIVER(2) / 13568,     0 /
C      DATA LOG10(1),LOG10(2) / 16282,  8347 /
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA RIGHT(1),RIGHT(2) / O032200, O000000 /
C      DATA DIVER(1),DIVER(2) / O032400, O000000 /
C      DATA LOG10(1),LOG10(2) / O037632, O020233 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA SMALL(1) / $00800000 /
C      DATA LARGE(1) / $7F7FFFFF /
C      DATA RIGHT(1) / $33800000 /
C      DATA DIVER(1) / $34000000 /
C      DATA LOG10(1) / $3E9A209B /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA RMACH(1) / O000400000000 /
C      DATA RMACH(2) / O377777777777 /
C      DATA RMACH(3) / O146400000000 /
C      DATA RMACH(4) / O147400000000 /
C      DATA RMACH(5) / O177464202324 /
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER.
C
C      DATA SMALL(1) /       128 /
C      DATA LARGE(1) /    -32769 /
C      DATA RIGHT(1) /     13440 /
C      DATA DIVER(1) /     13568 /
C      DATA LOG10(1) / 547045274 /
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER.
C
C      DATA RMACH(1) / Z00000080 /
C      DATA RMACH(2) / ZFFFF7FFF /
C      DATA RMACH(3) / Z00003480 /
C      DATA RMACH(4) / Z00003500 /
C      DATA RMACH(5) / Z209B3F9A /
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2.
C
C#ifdefC VMS
C      DATA RMACH(1) /       '80'X /
C      DATA RMACH(2) / 'FFFF7FFF'X /
C      DATA RMACH(3) /     '3480'X /
C      DATA RMACH(4) /     '3500'X /
C      DATA RMACH(5) / '209B3F9A'X /
C#endif
C
c aek
      write(i1mach(2),123)
 123  format('WARNING: R1MACH IS USED !!')
      IF (I .LT. 1  .OR.  I .GT. 5) GOTO 999
      R1MACH = RMACH(I)
      RETURN
  999 WRITE(I1MACH(2),1999) I
 1999 FORMAT(' R1MACH - I OUT OF BOUNDS',I10)
      STOP
      END
      DOUBLE PRECISION FUNCTION D1MACH(I)
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
c aek this version is only for use of d1mach(1)--d1mach(3)
c aek dmach(1): smallest possible value
c aek dmach(2): largest possible value
c aek dmach(3): smallest eps so that 1.d0+eps > 1.d0
      double precision dmach(3)
      integer i
      if (i .lt. 1 .or. i .gt. 3) go to 999
      dmach(1) = 1.d-99
      dmach(2) = 1.d+99
      dmach(3) = 1.d-15
      d1mach = dmach(i)
      return
 999  write(i1mach(2),123) i
 123  format(' D1MACH: I=',i1,' OUT OF BOUNDS !')
      stop
      end
C#define EXTENDED
C This adaptation of I1MACH, is identical to the public version, except
C that I1MACH has additional elements (17,18) specifying the length of a
C real word and of a double precision word.  BEWARE that these numbers
C have so far only been included for machines this program has been
C tested on.
C
      INTEGER FUNCTION I1MACH(I)
C
C  I/O UNIT NUMBERS.
C
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C
C  WORDS.
C
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
C                 FOR  FORTRAN 77, THIS IS ALWAYS 1.  FOR FORTRAN 66,
C                 CHARACTER STORAGE UNIT = INTEGER STORAGE UNIT.
C
C  INTEGERS.
C
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
C
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
C
C    I1MACH( 7) = A, THE BASE.
C
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C    BASE-B FORM
C
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
C
C    I1MACH(10) = B, THE BASE.
C
C  SINGLE-PRECISION
C
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
C *** Extensions: *** (only implemented for IEEE machines)
C
C    I1MACH(17) = number of integer words that fit into a real word
C
C    I1MACH(18) = number of integer words that fit into a real*8 word
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
C  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
C  WITH THE LOCAL OPERATING SYSTEM.  FOR FORTRAN 77, YOU MAY WISH
C  TO ADJUST THE DATA STATEMENT SO IMACH(6) IS SET TO 1, AND
C  THEN TO COMMENT OUT THE EXECUTABLE TEST ON I .EQ. 6 BELOW.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE, EXCEPT PERHAPS
C  FOR IMACH(1) - IMACH(4).
C
C#ifdef EXTENDED
      INTEGER IMACH(18),OUTPUT,SANITY
C#elseC
C      INTEGER IMACH(16),OUTPUT,SANITY
C#endif
C
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
C#ifdef IEEE | RIEEE
c aek      DATA IMACH( 1) /    5 /
      DATA IMACH( 2) /    6 /
c aek      DATA IMACH( 3) /    6 /
      DATA IMACH( 4) /    6 /
c aek      DATA IMACH( 5) /   32 /
      DATA IMACH( 6) /    4 /
c aek      DATA IMACH( 7) /    2 /
c aek      DATA IMACH( 8) /   31 /
c aek      DATA IMACH( 9) / 2147483647 /
c aek      DATA IMACH(10) /    2 /
c aek      DATA IMACH(11) /   24 /
c aek      DATA IMACH(12) / -125 /
c aek      DATA IMACH(13) /  128 /
c aek      DATA IMACH(14) /   53 /
c aek      DATA IMACH(15) / -1021 /
C#ifdef EXTENDED
c aek      DATA IMACH(16) /  1024 /
      DATA IMACH(17) /     1 /
      DATA IMACH(18) /     2 /, SANITY/987/
C#elseC
C      DATA IMACH(16) /  1024 /, SANITY/987/
C#endif
C#endif
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA IMACH( 1) /    7 /
C      DATA IMACH( 2) /    2 /
C      DATA IMACH( 3) /    2 /
C      DATA IMACH( 4) /    2 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   33 /
C      DATA IMACH( 9) / Z1FFFFFFFF /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -256 /
C      DATA IMACH(13) /  255 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) / -256 /
C      DATA IMACH(16) /  255 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -50 /
C      DATA IMACH(16) /  76 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -32754 /
C      DATA IMACH(16) /  32780 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / 00007777777777777777B /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   48 /
C      DATA IMACH(12) / -974 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   96 /
C      DATA IMACH(15) / -927 /
C      DATA IMACH(16) / 1070 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C#ifdefC CRAY
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /   102 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    63 /
C      DATA IMACH( 9) /  777777777777777777777B /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    47 /
C      DATA IMACH(12) / -8189 /
C      DATA IMACH(13) /  8190 /
C      DATA IMACH(14) /    94 /
C      DATA IMACH(15) / -8099 /
C#ifdefC EXTENDED
C      DATA IMACH(16) /  8190 /
C      DATA IMACH(17) /     1 /
C      DATA IMACH(18) /     1 /, SANITY/987/
C#elseC
C      DATA IMACH(16) /  8190 /, SANITY/987/
C#endifC
C#endif
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C      DATA IMACH( 1) /   11 /
C      DATA IMACH( 2) /   12 /
C      DATA IMACH( 3) /    8 /
C      DATA IMACH( 4) /   10 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) /32767 /
C      DATA IMACH(10) /   16 /
C      DATA IMACH(11) /    6 /
C      DATA IMACH(12) /  -64 /
C      DATA IMACH(13) /   63 /
C      DATA IMACH(14) /   14 /
C      DATA IMACH(15) /  -64 /
C      DATA IMACH(16) /   63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA IMACH( 1) /       5 /
C      DATA IMACH( 2) /       6 /
C      DATA IMACH( 3) /       0 /
C      DATA IMACH( 4) /       6 /
C      DATA IMACH( 5) /      24 /
C      DATA IMACH( 6) /       3 /
C      DATA IMACH( 7) /       2 /
C      DATA IMACH( 8) /      23 /
C      DATA IMACH( 9) / 8388607 /
C      DATA IMACH(10) /       2 /
C      DATA IMACH(11) /      23 /
C      DATA IMACH(12) /    -127 /
C      DATA IMACH(13) /     127 /
C      DATA IMACH(14) /      38 /
C      DATA IMACH(15) /    -127 /
C      DATA IMACH(16) /     127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /   43 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   63 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C#ifdefC IBM_VM | IBM370 | IBM3080 | IBM3090
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z7FFFFFFF /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C#ifdefC EXTENDED
C      DATA IMACH(16) /  63 /
C      DATA IMACH(17) /   1 /
C      DATA IMACH(18) /   2 /, SANITY/987/
C#elseC
C      DATA IMACH(16) /  63 /, SANITY/987/
C#endifC   EXTENDED
C#endif   IBM_VM
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   6 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z'7FFFFFFF' /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  62 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  62 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   54 /
C      DATA IMACH(15) / -101 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   62 /
C      DATA IMACH(15) / -128 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) / 32767 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
C     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
C     SUPPLIED BY IGOR BRAY.
C
C      DATA IMACH( 1) /            1 /
C      DATA IMACH( 2) /            1 /
C      DATA IMACH( 3) /            2 /
C      DATA IMACH( 4) /            1 /
C      DATA IMACH( 5) /           32 /
C      DATA IMACH( 6) /            4 /
C      DATA IMACH( 7) /            2 /
C      DATA IMACH( 8) /           31 /
C      DATA IMACH( 9) / :17777777777 /
C      DATA IMACH(10) /            2 /
C      DATA IMACH(11) /           23 /
C      DATA IMACH(12) /         -127 /
C      DATA IMACH(13) /         +127 /
C      DATA IMACH(14) /           47 /
C      DATA IMACH(15) /       -32895 /
C      DATA IMACH(16) /       +32637 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA IMACH( 1) /     0 /
C      DATA IMACH( 2) /     0 /
C      DATA IMACH( 3) /     7 /
C      DATA IMACH( 4) /     0 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     1 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    31 /
C      DATA IMACH( 9) /  2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    24 /
C      DATA IMACH(12) /  -125 /
C      DATA IMACH(13) /   128 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
C     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
C     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    6 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR VAX.
C
C#ifdefC VAX | VMS
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
CC
C#ifdefC EXTENDED
C      DATA IMACH(16) /  127 /
C      DATA IMACH(17) /    1 /
C      DATA IMACH(18) /    2 /, SANITY/987/
C#elseC
C      DATA IMACH(16) /  127 /, SANITY/987/
C#endifC
C#endif
C  ***  ISSUE STOP 777 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SANITY .NE. 987) STOP 777
C#ifdef EXTENDED
      IF (I .LT. 1  .OR.  I .GT. 18) GO TO 999
c aek
      if (i.ne.2.and.i.ne.4.and.i.ne.6.and.i.ne.17.and.i.ne.18) goto 99
      go to 98
 99   write(output,123) i
 123  format('WARNING: i1mach used with i=',i2)
C#elseC
C      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 999
C#endif
 98   I1MACH=IMACH(I)

C#ifdefC FTN77
C      IF(I.EQ.6) I1MACH=1
C#endif
C#ifdefC GFLOAT_VAX
CC machine constants for g_float optional compilation on the VAX
CC -------------------------------------------------------------
C      double precision dmach(5)
C      data dmach(1) / 1.112536929253601E-308 /
C      data dmach(2) / 4.494232837155789D+307 /
C      data dmach(3) / 1.110223024625157D-016 /
C      data dmach(4) / 2.220446049250313D-016 /
C      data dmach(5) / 0.301029995663981      /
C      d1mach = dmach(i)
C#endif
      RETURN
  999 WRITE(OUTPUT,1999) I
 1999 FORMAT(' I1MACH:  I OUT OF BOUNDS',I10)
      STOP
      END
