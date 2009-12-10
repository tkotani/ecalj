      integer function lstshr(mode,n1,lst1,n2,lst2,lstc)
C- Find elements common to two sorted vectors of integers
C ----------------------------------------------------------------
Ci Inputs
Ci   mode:    governs how output information is provided; see Remarks
Ci   lst1,n1: first list, and number
Ci   lst2,n2: second list, and number
Co Outputs
Co   lst: information about the elements common to both lists,
Co        information according to mode (see Remarks)
Co   lstshr: number of entries in lst
Cr Remarks
Cr   mode 0: lstc = list of integers common to both
Cr        1: lstc = vector of dimension (2,1..lstshr)
Cr                  lstc(1,*) points to entry in lst1
Cr                  lstc(2,*) points to entry in lst2
Cr        2: lstc = vector of dimension (3,1..lstshr)
Cr                  entries are grouped according to start, number:
Cr                  lstc(1,*) points to 1st entry in lst1
Cr                  lstc(2,*) points to 1st entry in lst2
Cr                  lstc(3,*) holds number of consecutive entries that
Cr                            are common to lst1,lst2
Cr        3: lstc as in mode 2, but integers must be sequential in order
Cr                            to be grouped together.
C ----------------------------------------------------------------
C     implicit none
      integer mode,n1,n2,lst1(n1),lst2(n2),lstc(5)
C Local
      integer i1l,i2l,i1h,i2h,nc
      logical lequ

      nc = 0
      i1l = 1
      i2l = 1
C --- Find the next common point ---
   10 continue
      if (i1l .gt. n1 .or. i2l .gt. n2) goto 99
      if (lst2(i2l)-lst1(i1l)) 11,13,12
C ... lst1 is larger
   11 continue
      i2l = i2l+1
      goto 10
C ... lst2 is larger
   12 continue
      i1l = i1l+1
      goto 10
C ... They are equal ... add to list
   13 continue
      nc = nc + 1

C --- Add to list ---
      if (mode .eq. 0) then
        lstc(nc) = lst1(i1l)
        i1l = i1l+1
        i2l = i2l+1
        goto 10
      elseif (mode .eq. 1) then
        lstc(2*nc-1) = i1l
        lstc(2*nc) = i2l
        i1l = i1l+1
        i2l = i2l+1
        goto 10
      elseif (mode .eq. 2 .or. mode .eq. 3) then
        lstc(3*nc-2) = i1l
        lstc(3*nc-1) = i2l
        i1h = i1l
        i2h = i2l
C   ... find top
   21   continue
        i1h = i1h+1
        i2h = i2h+1
        if (i1h .le. n1 .and. i2h .le. n2) then
          lequ = lst2(i2h) .eq. lst1(i1h)
          if (mode .eq. 3) lequ = lequ .and.
     .                            lst1(i1h) .eq. lst1(i1h-1)+1
          if (lequ) goto 21
          lstc(3*nc) = i1h-i1l
          i1l = i1h
          i2l = i2h
          goto 10
        else
          lstc(3*nc) = i1h-i1l
          goto 99
        endif
      else
        call rxi('lstshr: bad mode',mode)
      endif

C ... Increment i1l until at least as big as i2l
      if (lst1(i1l) .lt. lst2(i2l)) then
        i1l = i1l+1
        goto 10
      endif
C ... Increment i2l until at least as big as i1l
      if (lst1(i1l) .gt. lst2(i2l)) then
        i2l = i2l+1
        goto 10
      endif


   99 lstshr = nc
      return

      end
C test lstshr
C      subroutine fmain
C      implicit none
C      integer nclmx,ncnmx
C      parameter (nclmx=5,ncnmx=5)
C      integer cln(0:ncnmx,2,nclmx),lstc(ncnmx*3),i,mode,ic,lstshr
C      integer lst(ncnmx)
C
C      data cln /
C     . 3, 11,22,44,0,0,  5, 11,22,33,44,55,
C     . 4, 11,22,44,55,0, 5, 1,22,33,44,55,
C     . 4, 2,3,4,5,6,  5, 1,2,3,4,5,
C     . 5, 1,2,3,4,6,  4, 3,4,5,6,7,
C     . 3, 2,4,5,0,0,  4, 1,2,3,4,5/
C
C      do  10  ic = 1, nclmx
C        call awrit2('lst1  %n:1i',' ',255,6,cln(0,1,ic),cln(1,1,ic))
C        call awrit2('lst2  %n:1i',' ',255,6,cln(0,2,ic),cln(1,2,ic))
C        do  12  mode = 0, 3
C          call iinit(lstc,ncnmx*3)
C          i = lstshr(mode,cln(0,1,ic),cln(1,1,ic),cln(0,2,ic),
C     .      cln(1,2,ic),lstc)
C          if (mode .eq. 0) then
C            call awrit4('mode=%i returned %i : %n:1i',' ',255,6,mode,i,
C     .        i,lstc)
C          elseif (mode .eq. 1) then
C           call icopy(i,lstc(1),2,lst,1)
C           call awrit4('mode=%i returned %i : %n:1i',' ',255,6,mode,i,
C     .       i,lst)
C           call icopy(i,lstc(2),2,lst,1)
C           call awrit2('                    %n:1i',' ',255,6,i,lst)
C          else
C           call icopy(i,lstc(1),3,lst,1)
C           call awrit4('mode=%i returned %i : %n:1i',' ',255,6,mode,i,
C     .       i,lst)
C           call icopy(i,lstc(2),3,lst,1)
C           call awrit2('                    %n:1i',' ',255,6,i,lst)
C           call icopy(i,lstc(3),3,lst,1)
C           call awrit2('                    %n:1i',' ',255,6,i,lst)
C          endif
C   12   continue
C   10 continue
C      end

