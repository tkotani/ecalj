c$$$C --- Simple FORTRAN dynamic memory allocation ---
c$$$C Memory is allocated from a single heap.  The heap is declared as
c$$$C an integer array in common block W by the main program.  Its size
c$$$C is fixed by the size of W as declared by the main program.  It
c$$$C is possible to allocated memory from the top of the heap, and to
c$$$C free memory already allocated from the top of the heap.  It is
c$$$C not possible to free memory allocated from the middle of the heap.
c$$$C
c$$$C To use these routines, first call WKINIT with NSIZE = size of W.
c$$$C
c$$$C Character, integer, real, double precision, complex or double
c$$$C complex arrays can then be allocated by calling DEFCH, DEFI,
c$$$C DEFR, DEFDR, DEFC or DEFDC with an index IPNAME and LENGTH the
c$$$C desired length.  IPNAME is returned as the effective offset of
c$$$C the integer array W.  If it is desired that the array be
c$$$C initialized to zero, pass LENGTH as the negative of the desired
c$$$C length.
c$$$C
c$$$C Allocated memory can be freed by calling RLSE(IPNAME) where
c$$$C IPNAME must one of the indices from a previous allocation.  All
c$$$C blocks of memory allocated after IPNAME are also released.
c$$$C
c$$$C WKPRNT turns on or off a trace that prints a message each time
c$$$C        memory is allocated or released.  Also when turned on,
c$$$C        runs through links each time RLSE is called.
c$$$C DEFASK returns the amount of memory remaining in the pool
c$$$C WKSAV  sets an internal switch so that a heap allocation request
c$$$C        returns the the negative of length of the request in IPNAME
c$$$C        when not enough memory is available, rather than aborting.
c$$$C        In this case, no memory is allocated and no other internal
c$$$C        variables are changed.
c$$$C WKINFO prints information about the arrays allocated
c$$$C WKCHK  runs through links to see if any array was overwritten
c$$$C        call wkchk(string), where string terminated by '$'
c$$$C ----------------------------------------------------------------
c$$$      SUBROUTINE WKINIT(NSIZE)
c$$$C- Initialize conditions for the 'common pool' data W
c$$$ca !!!!!!!!!!!
c$$$ca    INTEGER W(2)
c$$$      INTEGER W(1)
c$$$      CHARACTER*1 STRING(60),STR(60)
c$$$      logical lopt,lerr
c$$$      SAVE
c$$$      COMMON /W/ W
c$$$
c$$$      IRND(I) = (I+499)/1000
c$$$
c$$$C ----- DEFINE STORAGE SIZE ------
c$$$C  START OF FIRST ARRAY AND MAX NUMBER TO BE DEFINED:
c$$$      IP0 = 5
c$$$      NDEFMX = 100
c$$$      LIMIT = NSIZE
c$$$      IPMAX = 0
c$$$      NDFDMX = 0
c$$$      NDFSUM = 0
c$$$      JPR = 0
c$$$      IPFREE = 5
c$$$      lerr = .false.
c$$$      WRITE(*,*) 'WKINIT:  size=',IRND(NSIZE),'K'
c$$$      WRITE(*,*) ' '
c$$$      RETURN
c$$$      ENTRY WKPRNT(JPRINT)
c$$$C- Set debug switch for heap management
c$$$      if (jprint .eq. 2) then
c$$$        jpr = 1-jpr
c$$$      else
c$$$        JPR = JPRINT
c$$$      endif
c$$$      RETURN
c$$$      entry wksav(lopt)
c$$$      lerr = lopt
c$$$      return
c$$$C ------ SUBROUTINES TO DEFINE ARRAYS OF VARIOUS TYPES -----
c$$$      ENTRY DEFCH(IPNAME,LENG)
c$$$C- Allocate character array
c$$$      LENGTH = (LENG+3)/4
c$$$      GOTO 10
c$$$      ENTRY DEFI(IPNAME,LENG)
c$$$C- Allocate integer array
c$$$      LENGTH = LENG
c$$$      JOPT = 1
c$$$      GOTO 10
c$$$      ENTRY DEFR(IPNAME,LENG)
c$$$C- Allocate single precision real array
c$$$      LENGTH = LENG*I1MACH(17)
c$$$      JOPT = 2
c$$$      GOTO 10
c$$$      ENTRY DEFRR(IPNAME,LENG)
c$$$      ENTRY DEFDR(IPNAME,LENG)
c$$$C- Allocate double precision real array
c$$$      LENGTH = LENG*I1MACH(18)
c$$$      JOPT = 3
c$$$      GOTO 10
c$$$      ENTRY DEFC(IPNAME,LENG)
c$$$C- Allocate single precision complex array
c$$$      LENGTH = LENG*2*I1MACH(17)
c$$$      JOPT = 4
c$$$      GOTO 10
c$$$      ENTRY DEFCC(IPNAME,LENG)
c$$$      ENTRY DEFDC(IPNAME,LENG)
c$$$C- Allocate double precision complex array
c$$$      LENGTH = LENG*2*I1MACH(18)
c$$$      JOPT = 5
c$$$   10 IOPT = 0
c$$$      IF (LENGTH .LT. 0) THEN
c$$$        IOPT = 1
c$$$        LENGTH = -LENGTH
c$$$      ENDIF
c$$$      IF (LENGTH .EQ. 0) LENGTH = 1
c$$$      IMOD = 1
c$$$      GOTO 83
c$$$   84 IPNAME = IPFREE
c$$$      if (lerr .and. ipfree+length+2 .gt. limit) then
c$$$        ipname = -LENGTH
c$$$        if (jpr .gt. 0)  print *,
c$$$     .    'ALLOC: heap storage exceeded; returning -LENGTH=',-LENGTH
c$$$        return
c$$$      endif
c$$$      IPFREE = IPFREE + LENGTH + 1
c$$$      IPFREE = 4*((IPFREE+2)/4)+1
c$$$      IPMAX = MAX0(IPMAX,IPFREE)
c$$$      W(IPNAME-1)=IPFREE
c$$$      NDEFD = NDEFD + 1
c$$$      NDFDMX = MAX0(NDFDMX,NDEFD)
c$$$      NDFSUM = NDFSUM + LENGTH
c$$$      IF (JPR .GT. 0) WRITE(*,100) NDEFD,LENG,LENGTH,IPNAME,IPFREE-1
c$$$  100 FORMAT(' define array',I4,':   els=',I8,'   length=',I8,',',
c$$$     .   I8,'  to',I8)
c$$$      IF (IPFREE .LE. LIMIT) THEN
c$$$c takao apl9 2008
c$$$c        IF (IOPT .NE. 0) GOTO (201,202,203,204,205) JOPT
c$$$        IF (IOPT .NE. 0) then
c$$$          if(JOPT==1)  GOTO 201
c$$$          if(JOPT==2)  GOTO 202
c$$$          if(JOPT==3)  GOTO 203
c$$$          if(JOPT==4)  GOTO 204
c$$$          if(JOPT==5)  GOTO 205
c$$$        ENDIF
c$$$
c$$$        RETURN
c$$$  201   CALL IINIT(W(IPNAME),-LENG)
c$$$        RETURN
c$$$  202   CALL AINIT(W(IPNAME),-LENG)
c$$$        RETURN
c$$$  203   CALL DINIT(W(IPNAME),-LENG)
c$$$        RETURN
c$$$  204   CALL CINIT(W(IPNAME),-LENG)
c$$$        RETURN
c$$$  205   CALL ZINIT(W(IPNAME),-LENG)
c$$$        RETURN
c$$$      ENDIF
c$$$      WRITE(*,101) IPFREE
c$$$  101 FORMAT(' ALLOC: WORKSPACE OVERFLOW, NEED AT LEAST',I8)
c$$$      STOP
c$$$C- Release data up to pointer
c$$$      ENTRY RLSE(IPNAME)
c$$$      IF (IPNAME .GT. LIMIT) STOP 'RLSE: release pointer exceeds limit'
c$$$      IF (IPNAME .LT. 3) STOP 'RLSE: release pointer less than 3'
c$$$      IF (JPR .eq. 0) goto 82
c$$$      imod = 3
c$$$      goto 83
c$$$   87 WRITE(*,*) 'RLSE from: ',IPNAME
c$$$   82 IPFREE = IPNAME
c$$$      return
c$$$      ENTRY DEFASK(LREST)
c$$$C- Return number of words left in common pool
c$$$      LREST = LIMIT - IPFREE - 2
c$$$      IF (JPR .GT. 0) WRITE(*,*) 'SPACE LEFT=',LREST,'  SINGLE WORDS'
c$$$      RETURN
c$$$      ENTRY WKINFO()
c$$$C- Output workspace information
c$$$      IMOD = 2
c$$$      GOTO 83
c$$$  81  WRITE(*,601) IRND(LIMIT),IRND(NDFSUM),IRND(IPMAX-1),
c$$$     .             IRND(IPFREE-1),NDFDMX,NDEFD
c$$$  601 FORMAT(
c$$$     .  /'  total workspace size =',I5,' K',
c$$$     .  /'  total space allocated=',I5,' K',
c$$$     .  /'  workspace used:    max',I5,' K   now',I4,' K',
c$$$     .  /'  arrays defined:    max',I7,  '   now',I6)
c$$$      IF (IPFREE .EQ. IP0) RETURN
c$$$      if (jpr .gt. 0) WRITE(*,602)
c$$$  602 FORMAT(/'  array',6X,'begin',7X,'end',7X,'length')
c$$$      IPX = IP0
c$$$      DO  30  I = 1, NDEFMX
c$$$        IPY = W(IPX-1)
c$$$        IF (IPX .EQ. IPFREE) RETURN
c$$$        IF (IPY .LT. IP0 .OR. IPY .GT. LIMIT) WRITE(*,*) '   . . . . . '
c$$$        IF (IPY .LT. IP0 .OR. IPY .GT. LIMIT) RETURN
c$$$        IF (JPR .GT. 0) WRITE(*,603) I,IPX,IPY,IPY-IPX
c$$$  603   FORMAT(4(I6,5X))
c$$$        IPX = IPY
c$$$   30 CONTINUE
c$$$      RETURN
c$$$      ENTRY WKCHK(STRING)
c$$$C- Run through links to see if any dynamic array was overwritten
c$$$      IMOD = 0
c$$$      DO  88  I = 1, 60
c$$$        STR(I)=STRING(I)
c$$$        NSTR = I-1
c$$$        IF (STRING(I) .EQ. '$') GOTO 89
c$$$   88 continue
c$$$   89 WRITE(*,*) 'WKCHK: ',(STR(I),I = 1,NSTR)
c$$$   83 NDEFD = 0
c$$$      IPX = IP0
c$$$      IPPLOC = -999
c$$$      DO  35  I = 1, NDEFMX
c$$$        IF (IPX .LT. IP0 .OR. IPX .GT. LIMIT) THEN
c$$$          WRITE(*,888) NDEFD,IPX,IPPLOC
c$$$  888     FORMAT(' ALLOC: LINK DESTROYED AT START OF ARRAY',I3,
c$$$     .     ',  PTR=',I8,' AT',I8)
c$$$          STOP
c$$$        ENDIF
c$$$        IF (IPX .EQ. IPFREE) GOTO 86
c$$$        NDEFD = NDEFD + 1
c$$$        IPPLOC = IPX - 1
c$$$        IPX = W(IPPLOC)
c$$$   35 continue
c$$$  86  CONTINUE
c$$$c apr9 2008
c$$$c      GOTO (84,81,87), imod
c$$$      if(imod==1) goto 84
c$$$      if(imod==2) goto 81
c$$$      if(imod==3) goto 87
c$$$
c$$$      WRITE(*,360) NDEFD,IPFREE-1
c$$$  360 FORMAT('     LINKS OK   NDEFD=',I3,'   SPACE USED=',I7)
c$$$      RETURN
c$$$      END
c$$$      SUBROUTINE CINIT(ARRAY,LENG)
c$$$C- Initializes complex array to zero
c$$$      INTEGER LENG
c$$$      REAL ARRAY(2*LENG)
c$$$      CALL AINIT(ARRAY,LENG+LENG)
c$$$      RETURN
c$$$      END
c$$$      SUBROUTINE ZINIT(ARRAY,LENG)
c$$$C- Initializes complex array to zero
c$$$      INTEGER LENG
c$$$      DOUBLE PRECISION ARRAY(2*LENG)
c$$$      CALL DINIT(ARRAY,LENG+LENG)
c$$$      RETURN
c$$$      END
c$$$      SUBROUTINE IINIT(ARRAY,LENG)
c$$$C- Initializes integer array to zero
c$$$      INTEGER LENG
c$$$      INTEGER ARRAY(LENG)
c$$$      DO 10 I=1,LENG
c$$$        ARRAY(I)=0
c$$$   10 CONTINUE
c$$$      RETURN
c$$$      END
c$$$      SUBROUTINE AINIT(ARRAY,LENG)
c$$$C- Initializes real array to zero
c$$$      INTEGER LENG
c$$$      REAL ARRAY(LENG)
c$$$      DO 10 I=1,LENG
c$$$        ARRAY(I)=0.
c$$$   10 CONTINUE
c$$$      RETURN
c$$$      END
c$$$      SUBROUTINE DINIT(ARRAY,LENG)
c$$$C- Initializes double precision array to zero
c$$$      INTEGER LENG
c$$$      DOUBLE PRECISION ARRAY(LENG)
c$$$      DO 10 I=1,LENG
c$$$        ARRAY(I)=0.D0
c$$$   10 CONTINUE
c$$$      RETURN
c$$$      END
c$$$      FUNCTION RVAL(ARRAY,INDEX)
c$$$C- Returns the real value of ARRAY(INDEX)
c$$$      REAL ARRAY(INDEX)
c$$$      RVAL=ARRAY(INDEX)
c$$$      RETURN
c$$$      END
c$$$      DOUBLE PRECISION FUNCTION DRVAL(ARRAY,INDEX)
c$$$C- Returns the double precision value of ARRAY(INDEX)
c$$$      DOUBLE PRECISION ARRAY(INDEX)
c$$$      DRVAL=ARRAY(INDEX)
c$$$      RETURN
c$$$      END
c$$$      INTEGER FUNCTION IVAL(ARRAY,INDEX)
c$$$C- Returns the integer value of ARRAY(INDEX)
c$$$      INTEGER ARRAY(INDEX)
c$$$      IVAL=ARRAY(INDEX)
c$$$      RETURN
c$$$      END
c$$$      COMPLEX FUNCTION CVAL(ARRAY,INDEX)
c$$$C- Returns the complex value of ARRAY(INDEX)
c$$$      COMPLEX ARRAY(INDEX)
c$$$      CVAL=ARRAY(INDEX)
c$$$      RETURN
c$$$      END
c$$$
c$$$
c$$$
c$$$      SUBROUTINE DMADD(A,NCA,NRA,SCALEA,B,NCB,NRB,SCALEB,C,NCC,NRC,N,M)
c$$$C- general matrix addition
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs:
c$$$Ci   a,nca,nra is the left matrix and respectively the number of
c$$$Ci      elements separating columns and rows.
c$$$Ci   b,ncb,nrb is the right matrix and respectively the number of
c$$$Ci      elements separating columns and rows.
c$$$Ci   c,ncc,nrc is the result matrix and respectively the number of
c$$$Ci      elements separating columns and rows.
c$$$Ci   n,m: the number of rows and columns, respectively, to calculate
c$$$Co Outputs:
c$$$Co   result matrix stored in c
c$$$Cr Remarks:
c$$$Cr   This is a general-purpose matrix linear combination routine,
c$$$Cr   adding a subblock of matrix a to a subblock of matrix b.
c$$$Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
c$$$Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
c$$$Cr   Arrays are locally one-dimensional so as to optimize inner loop.
c$$$Cr
c$$$Cr   Destination matrix c can coincide with either a or b, provided that
c$$$Cr   the transpose of the coincident matrix is not taken.
c$$$Cr   Example: Add 3-by-2 block of (transpose of a - .5*b) into c
c$$$Cr     call dmadd(a,1,na,1.d0,b,nb,0,-.5d0,c,nc,1,3,2)
c$$$Cr     OLD call dmadd(a,na,1,b,nb,0,-.5d0,c,nc,3,2)
c$$$C ----------------------------------------------------------------
c$$$C
c$$$      INTEGER NCA,NRA,NCB,NRB,NCC,NRC,N,M
c$$$      DOUBLE PRECISION A(0:1), B(0:1), C(0:1), SCALEA, SCALEB
c$$$      INTEGER I,J,IA,IB,IC
c$$$
c$$$      DO  200  I = N-1, 0, -1
c$$$        IA = I*NRA+M*NCA
c$$$        IB = I*NRB+M*NCB
c$$$        IC = I*NRC+M*NCC
c$$$      DO  200  J = M-1, 0, -1
c$$$        IA = IA-NCA
c$$$        IB = IB-NCB
c$$$        IC = IC-NCC
c$$$        C(IC) = A(IA)*SCALEA + B(IB)*SCALEB
c$$$  200 CONTINUE
c$$$      RETURN
c$$$      END
c$$$      subroutine dpmpy(a,b,nscb,nsrb,c,nscc,nsrc,nr,nc,l)
c$$$C- matrix multiplication, (packed) (normal) -> (normal)
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs:
c$$$Ci   a is the left matrix (packed)
c$$$Ci   b,nscb,nsrb is the right matrix and respectively the spacing
c$$$Ci      between column elements and row elements.
c$$$Ci   c,nscc,nsrc is the product matrix and respectively the number of
c$$$Ci      elements separating columns and rows.
c$$$Ci   nr,nc: the number of rows and columns, respectively, to calculate
c$$$Ci   l:   length of vector for matrix multiply
c$$$Co Outputs:
c$$$Co   product matrix stored in c
c$$$Cr Remarks:
c$$$Cr   This is a general-purpose matrix multiplication routine,
c$$$Cr   multiplying a subblock of matrix a by a subblock of matrix b.
c$$$Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
c$$$Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
c$$$Cr   Arrays are locally one-dimensional so as to optimize inner loop,
c$$$Cr   which is executed nr*nc*l times.  No attempt is made to optimize
c$$$Cr   the outer loops, executed nr*nc times.
c$$$Cr     Examples: product of (nr,l) subblock of a into (l,nc) subblock of
c$$$Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,nr,nc,l)
c$$$Cr     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
c$$$Cr     To generate the tranpose of that product, use:
c$$$Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,nr,nc,l)
c$$$C ----------------------------------------------------------------
c$$$C Passed Parameters
c$$$      integer nscb,nsrb,nscc,nsrc,nr,nc,l
c$$$      double precision a(0:*), b(0:*), c(0:*)
c$$$C Local parameters
c$$$      double precision sum
c$$$      integer i,j,k,offa,offb
c$$$
c$$$      do  20  i = 0, nr-1
c$$$      do  20  j = 0, nc-1
c$$$        sum = 0
c$$$        offa = (i*(i+1))/2
c$$$        offb = nscb*j
c$$$        do  21  k = 0, i-1
c$$$          sum = sum + a(offa)*b(offb)
c$$$          offa = offa + 1
c$$$          offb = offb + nsrb
c$$$   21   continue
c$$$        do  22  k = i, l-1
c$$$          sum = sum + a(offa)*b(offb)
c$$$          offa = offa + k+1
c$$$          offb = offb + nsrb
c$$$   22   continue
c$$$        c(i*nsrc+j*nscc) = sum
c$$$   20 continue
c$$$      end

C       REAL FUNCTION R1MACH(I)
C C
C C  SINGLE-PRECISION MACHINE CONSTANTS
C C
C C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C C
C C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C C
C C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C C
C C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C C
C C  R1MACH(5) = LOG10(B)
C C
C C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C C  REMOVING THE C FROM COLUMN 1.
C C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C C
C C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C C
C C  WHERE POSSIBLE, OCTAL OR HEXADECIMAL CONSTANTS HAVE BEEN USED
C C  TO SPECIFY THE CONSTANTS EXACTLY WHICH HAS IN SOME CASES
C C  REQUIRED THE USE OF EQUIVALENT INTEGER ARRAYS.
C C
C       INTEGER SMALL(2)
C       INTEGER LARGE(2)
C       INTEGER RIGHT(2)
C       INTEGER DIVER(2)
C       INTEGER LOG10(2)
C C
C       REAL RMACH(5)
C C
C       EQUIVALENCE (RMACH(1),SMALL(1))
C       EQUIVALENCE (RMACH(2),LARGE(1))
C       EQUIVALENCE (RMACH(3),RIGHT(1))
C       EQUIVALENCE (RMACH(4),DIVER(1))
C       EQUIVALENCE (RMACH(5),LOG10(1))
C C
C C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C C
C C#ifdef IEEE | RIEEE
C       DATA SMALL(1) /     8388608 /
C       DATA LARGE(1) /  2139095039 /
C       DATA RIGHT(1) /   864026624 /
C       DATA DIVER(1) /   872415232 /
C       DATA LOG10(1) /  1050288283 /
C C#endif
C C
C C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C C
C C      DATA SMALL(1) /    1048576 /
C C      DATA LARGE(1) / 2147483647 /
C C      DATA RIGHT(1) /  990904320 /
C C      DATA DIVER(1) / 1007681536 /
C C      DATA LOG10(1) / 1091781651 /
C C
C C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C C
C C      DATA RMACH(1) / Z400800000 /
C C      DATA RMACH(2) / Z5FFFFFFFF /
C C      DATA RMACH(3) / Z4E9800000 /
C C      DATA RMACH(4) / Z4EA800000 /
C C      DATA RMACH(5) / Z500E730E8 /
C C
C C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.
C C
C C      DATA RMACH(1) / O1771000000000000 /
C C      DATA RMACH(2) / O0777777777777777 /
C C      DATA RMACH(3) / O1311000000000000 /
C C      DATA RMACH(4) / O1301000000000000 /
C C      DATA RMACH(5) / O1157163034761675 /
C C
C C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C C
C C      DATA RMACH(1) / 00014000000000000000B /
C C      DATA RMACH(2) / 37767777777777777777B /
C C      DATA RMACH(3) / 16404000000000000000B /
C C      DATA RMACH(4) / 16414000000000000000B /
C C      DATA RMACH(5) / 17164642023241175720B /
C C
C C     MACHINE CONSTANTS FOR CONVEX C-1.
C C
C C      DATA RMACH(1) / '00800000'X /
C C      DATA RMACH(2) / '7FFFFFFF'X /
C C      DATA RMACH(3) / '34800000'X /
C C      DATA RMACH(4) / '35000000'X /
C C      DATA RMACH(5) / '3F9A209B'X /
C C
C C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C C
C C#ifdefC CRAY
C C      DATA RMACH(1) / 200034000000000000000B /
C C      DATA RMACH(2) / 577767777777777777776B /
C C      DATA RMACH(3) / 377224000000000000000B /
C C      DATA RMACH(4) / 377234000000000000000B /
C C      DATA RMACH(5) / 377774642023241175720B /
C C#endif
C C
C C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C C
C C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C C     STATIC RMACH(5)
C C
C C      DATA SMALL/20K,0/,LARGE/77777K,177777K/
C C      DATA RIGHT/35420K,0/,DIVER/36020K,0/
C C      DATA LOG10/40423K,42023K/
C C
C C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C C
C C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C C      DATA LARGE(1),LARGE(2) / '37777777, '00000177 /
C C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /
C C      DATA DIVER(1),DIVER(2) / '20000000, '00000353 /
C C      DATA LOG10(1),LOG10(2) / '23210115, '00000377 /
C C
C C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C C
C C      DATA RMACH(1) / O402400000000 /
C C      DATA RMACH(2) / O376777777777 /
C C      DATA RMACH(3) / O714400000000 /
C C      DATA RMACH(4) / O716400000000 /
C C      DATA RMACH(5) / O776464202324 /
C C
C C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C C
C C#ifdefC IBM_VM | IBM370 | IBM3080 | IBM3090
C C      DATA RMACH(1) / Z00100000 /
C C      DATA RMACH(2) / Z7FFFFFFF /
C C      DATA RMACH(3) / Z3B100000 /
C C      DATA RMACH(4) / Z3C100000 /
C C      DATA RMACH(5) / Z41134413 /
C C#endif
C C
C C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C C
C C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C C
C C      DATA RMACH(1) / Z'00100000' /
C C      DATA RMACH(2) / Z'7EFFFFFF' /
C C      DATA RMACH(3) / Z'3B100000' /
C C      DATA RMACH(4) / Z'3C100000' /
C C      DATA RMACH(5) / Z'41134413' /
C C
C C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C C
C C      DATA RMACH(1) / "000400000000 /
C C      DATA RMACH(2) / "377777777777 /
C C      DATA RMACH(3) / "146400000000 /
C C      DATA RMACH(4) / "147400000000 /
C C      DATA RMACH(5) / "177464202324 /
C C
C C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C C
C C      DATA SMALL(1) /    8388608 /
C C      DATA LARGE(1) / 2147483647 /
C C      DATA RIGHT(1) /  880803840 /
C C      DATA DIVER(1) /  889192448 /
C C      DATA LOG10(1) / 1067065499 /
C C
C C      DATA RMACH(1) / O00040000000 /
C C      DATA RMACH(2) / O17777777777 /
C C      DATA RMACH(3) / O06440000000 /
C C      DATA RMACH(4) / O06500000000 /
C C      DATA RMACH(5) / O07746420233 /
C C
C C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C C
C C      DATA SMALL(1),SMALL(2) /   128,     0 /
C C      DATA LARGE(1),LARGE(2) / 32767,    -1 /
C C      DATA RIGHT(1),RIGHT(2) / 13440,     0 /
C C      DATA DIVER(1),DIVER(2) / 13568,     0 /
C C      DATA LOG10(1),LOG10(2) / 16282,  8347 /
C C
C C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C C      DATA RIGHT(1),RIGHT(2) / O032200, O000000 /
C C      DATA DIVER(1),DIVER(2) / O032400, O000000 /
C C      DATA LOG10(1),LOG10(2) / O037632, O020233 /
C C
C C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C C
C C      DATA SMALL(1) / $00800000 /
C C      DATA LARGE(1) / $7F7FFFFF /
C C      DATA RIGHT(1) / $33800000 /
C C      DATA DIVER(1) / $34000000 /
C C      DATA LOG10(1) / $3E9A209B /
C C
C C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C C
C C      DATA RMACH(1) / O000400000000 /
C C      DATA RMACH(2) / O377777777777 /
C C      DATA RMACH(3) / O146400000000 /
C C      DATA RMACH(4) / O147400000000 /
C C      DATA RMACH(5) / O177464202324 /
C C
C C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER.
C C
C C      DATA SMALL(1) /       128 /
C C      DATA LARGE(1) /    -32769 /
C C      DATA RIGHT(1) /     13440 /
C C      DATA DIVER(1) /     13568 /
C C      DATA LOG10(1) / 547045274 /
C C
C C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C C     FORTRAN IV-PLUS COMPILER.
C C
C C      DATA RMACH(1) / Z00000080 /
C C      DATA RMACH(2) / ZFFFF7FFF /
C C      DATA RMACH(3) / Z00003480 /
C C      DATA RMACH(4) / Z00003500 /
C C      DATA RMACH(5) / Z209B3F9A /
C C
C C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2.
C C
C C#ifdefC VMS
C C      DATA RMACH(1) /       '80'X /
C C      DATA RMACH(2) / 'FFFF7FFF'X /
C C      DATA RMACH(3) /     '3480'X /
C C      DATA RMACH(4) /     '3500'X /
C C      DATA RMACH(5) / '209B3F9A'X /
C C#endif
C C
C c aek
C       write(i1mach(2),123)
C  123  format('WARNING: R1MACH IS USED !!')
C       IF (I .LT. 1  .OR.  I .GT. 5) GOTO 999
C       R1MACH = RMACH(I)
C       RETURN
C   999 WRITE(I1MACH(2),1999) I
C  1999 FORMAT(' R1MACH - I OUT OF BOUNDS',I10)
C       STOP
C       END

