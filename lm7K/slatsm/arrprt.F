      subroutine arrprt(hfmt,fmt,casts,ntab,nrr,ncol,nskk,sep,
     .a1,a2,a3,a4,a5,a6,a7,a8)
C- Print out a table in multi-column format
C ----------------------------------------------------------------------
Ci Inputs
Ci   hfmt  :format for header; see Remarks
Ci   fmt   :(string) awrite format for table; see Remarks
Ci   casts :(string) list of casts in table; see Remarks
Ci   ntab  :size of the table
Ci   nrr   :not used now
Ci   ncol  :print the table ncol columns
Ci   nskk  :not used now
Ci   sep   :(string) separator between columns
Ci   a1    :a1(1..ntab) hold elements for first member in table
Ci   a2    :a2(1..ntab) hold elements for second member in table
Ci   a3    :etc
Ci   a4    :etc
Ci   a5    :etc
Ci   a6    :etc
Ci   a7    :etc
Ci   a8    :etc
Co Outputs
Co   A table is written to stdo, in multicolumn format.
Cl Local variables
Cl   nmemb :number of members each entry has
Cl   is    :current group of quantities (1..ntab) being printed
Cr Remarks
Cr   arrprt prints a table of objects in multicolumn format, for example:
Cr
Cr   Class   dq     amom      Class   dq     amom      Class   dq     amom
Cr     1  -1.2898 -0.0147  |    5   0.1902 -0.0001  |    9  -0.0628  0.0035
Cr     2  -0.9817 -2.7848  |    6  -1.2856 -0.0005  |   10   0.6502  0.0003
Cr     3  -0.1006 -0.0240  |    7  -0.0625 -0.0029  |
Cr     4   0.6511 -0.0008  |    8   0.6506 -0.0003  |
Cr
Cr   Each object is a group of numbers (in the example, the group consists
Cr   of three numbers.  There is a maximum of 8 numbers in a group).
Cr
Cr   Each number is specified as one of:
Cr     - the index to the current entry in the table, i.e. 1 for the first
Cr       object, 2 for the second, etc (in the example, see the first entry)
Cr     - an integer (no entries in the example are integers)
Cr     - a double-precision number (2nd and 3rd entries in the example)
Cr
Cr   For the second and third specifications, you pass this information in
Cr   the arrays a1,a2,a3, etc.  Thus a1 holds a vector of numbers for the
Cr   first entry (it's not used in the example, since first entry is an
Cr   index); a2 holds a vector of numbers for the second entry (in the
Cr   example, a double precision array), etc.
Cr
Cr   Argument 'casts' is a string that specifies both how many numbers are
Cr   in a group, and which type each member of the group is.  The first
Cr   character specifies the type of the first member, the second the type
Cr   of the second, etc.  The size of the 'casts' specifies the size of the
Cr   group.  Each character in 'casts' must be one of 'I', 'i', or 'd',
Cr   specifying that is respectively an index, an, integer, or a double
Cr   precision number).  Thus in the above example, casts='Idd'.
Cr
Cr   'hfmt' specifies the header appearing at the top of each column.
Cr   In the example hfmt=' Class   dq     amom'
Cr
Cr   How each object appears is specified by the 'fmt' argument, which is
Cr   an awrite format.  In the example, fmt='%,4i%:-3,4;4d%:-2,4;4d'
Cr
Cr   Columns are separated by string 'sep'.  In the example, sep='  | '
Cu Updates
Cu   01 Jun 01 Input simplified: nrr and nskk no longer used.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) hfmt,fmt,casts,sep
      integer ncol,nrr,nskk,ntab
      double precision a1(1),a2(1),a3(1),a4(1),a5(1),a6(1),a7(1),a8(1)
C ... Local parameters
      integer i1,i2,ia,ip,ic,is,ir,lgunit,nrow,awrite,nmemb
      character*(150) sout,head
      double precision xx(8)
C ... External calls
Ckino 2011.05 del external,       external awrit0,vwrt,word

      nrow = (ntab+ncol-1)/ncol
      call word(casts,1,i1,i2)
      nmemb = i2-i1+1
      if (ncol .le. 0 .or. nmemb .le. 0) return

      do  10  ir = 1, nrow
        ip = 1
        sout = ' '
        is = ir-nrow
        do  20  ic = 1, ncol
          is = is + nrow
          if (is .gt. ntab) goto 20

          do  30  ia = 1, nmemb
C  ... 'I' prints out column index
            if (casts(ia:ia) .eq. 'I')  then
              call vwrt(1,1,is,a2,a3,a4,a5,a6,a7,a8,2,xx(ia),xx(ia))
C  ...  this argment is of cast integer
            elseif (casts(ia:ia) .eq. 'i')  then
              call vwrt(ia,is,a1,a2,a3,a4,a5,a6,a7,a8,2,xx(ia),xx(ia))
C  ...  this argment is of cast double
            else
              call vwrt(ia,is,a1,a2,a3,a4,a5,a6,a7,a8,4,xx(ia),xx(ia))
            endif
   30     continue
          head(ip:) = hfmt
          ip = ip + awrite(fmt,sout(ip:),len(sout),0,xx(1),xx(2),xx(3),
     .    xx(4),xx(5),xx(6),xx(7),xx(8))
          if (ic .lt. ncol) then
            sout(ip:ip-1+len(sep)) = sep
            ip = ip + len(sep)
          endif
   20   continue

        if (ir .eq. 1) call awrit0(head,' ',-len(head),lgunit(1))
        call awrit0(sout,' ',-len(sout),lgunit(1))
   10 continue
      end

