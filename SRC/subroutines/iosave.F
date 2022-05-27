      subroutine iosave(flg,vstrn,vars,ifi,nvario)
      use m_ftox
C- Write line to save file
C ----------------------------------------------------------------
Ci Inputs
Ci   flg: 1-character label
Ci        optional second character can be a digit specifying
Ci                 precision to save output variables.
Ci                 If character is missing or not between
Ci                 '0' and '9', defaults to '6'
Ci        optional third character can be a digit specifying
Ci                 precision to save vstrn variables
Ci                 If character is missing or not between
Ci                 '0' and '9', defaults to '6'
Ci        optional fourth character, set to 'f' or 'F'
Ci                 suppresses rewinding file.
Ci   vstrn: list of variable names, separated by commas, eg
Ci          "time,mmom,etot"
Ci   vars:  list of variables corresponding to string
Ci   ifi:   file unit (<0 means output)
Ci          ifi>0 not implemented
Ci   nvario: number of variables to i/o
Co Outputs
Co   variables and ehk,ehf written to -ifi
Cr Remarks
Cu Updates
Cu   17 Aug 01 Added (optional) fourth character to flg
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*(*) vstrn
      character*(*) flg
      double precision vars(1)
      integer ifi,nvario
C Local parameters
      logical leof
      integer nvar,i,ib
      double precision val
      integer lstr,k,k0,j,lgunit
      parameter (lstr=1024)
      character outstr*(lstr), nam*16, p1*1, p2*1
      integer::p1x,p2x
      call rxx(ifi.gt.0,'iosave:  attempt to read from save file')
C ... Print fractions w/out leading 0, to save space
      ib = -1
c      call bin2a0(ib)
c      call bin2a0(10)
      call numsyv(nvar)
      outstr = flg(1:1)
      p1 = '6'
      p2 = '6'
C --- Output variables up to nvario ---
      do  10  i = 4, min(nvario,nvar)
        nam = ' '
        call watsyv(nam,val,i)
        outstr=trim(outstr)//' '//trim(nam)//'='//ftom(val)
   10 continue
C --- Output all variables in string ---
      k0 = 0
      call skipbl(vstrn,len(vstrn),k0)
      i = 0
   20 k = k0
      i = i+1
      call chrps2(vstrn,', ',2,len(vstrn),k,j)
      nam = vstrn(k0+1:k)
      outstr=trim(outstr)//' '//trim(nam)//'='//ftof(vars(i),6)
      k0 = k+1
      if (j .eq. 1 .and. k .lt. len(vstrn)) goto 20
      if (leof) call poseof(-ifi)
      write(-ifi,"(a)") trim(outstr)
c      call bin2a0(ib)
      end

