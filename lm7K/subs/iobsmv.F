      subroutine iobsmv(mode,tim,T,V,nbas,alat,pos,ifi)
C- File I/O for new atom positions, xbs format
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 not used now.  Must be zero
Ci   tim   :time
Ci   T     :kinetic energy
Ci   V     :potential energy
Ci   nbas  :size of basis
Ci   pos   :position vectors
Ci   ifi   :file logical unit
Ci         :ifi<0 file WRITE
Ci         :ifi>0 file READ
Co Outputs :position vectors are written to disk
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   18 Nov 02
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,ifi
      double precision tim,T,V,alat,pos(3*nbas)
C ... Local parameters
      integer i,jfi
C     double precision

      call isanrg(mode,0,0,' iobsmv','mode',.true.)
      if (ifi .lt. 0) then
        jfi = -ifi
        call awrit4('frame t= %g T= %g V= %g T+V= %g',' ',80,jfi,
     .  tim,T,V,T+V)
        write(jfi,100) (alat*pos(i), i=1, 3*nbas)
  100   format(3000f8.3)
        return
      else
        call rx('iobsmv: file READ not implemented')
      endif

      end

      subroutine iobsm0(mode,tim,T,V,nbas,alat,pos,ifi)
C- Initialization for iobsmv
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0
Ci         :position file to EOF
Ci         :Attempt to read time from last write.
Ci         :  If successful, return
Ci         :  Else, set time=0; write initial positions
Ci         :No other modes implemented at present
Ci   ifi   :file handle
Co Outputs
Co   tim   :time at last step, if file
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   18 Nov 02
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,ifi
      double precision tim,T,V,alat,pos(3*nbas)
C ... Local parameters
      logical a2bin
      integer recln,nw,j1,j2,i
      parameter (recln=80)
      character*(recln) s

      call isanrg(mode,0,0,' iobsm0','mode',.true.)
      call poseof(ifi)
      backspace ifi
      backspace ifi
C     Read second-to-last line; jump to 999 if missing
      read(ifi,'(a80)',err=999,end=999) s
      call words(s,nw)
C     Various checks before reading tim
C     This line should have structure 'frame t= # ...'
      if (nw .lt. 3) goto 998
      call word(s,1,j1,j2)
      if (s(j1:j2) .ne. 'frame') goto 998
      call word(s,2,j1,j2)
      if (s(j1:j1+1) .ne. 't=') goto 998
      call word(s,3,i,j2)
      i = j1+1
      if (.not. a2bin(s,tim,4,0,' ',i,80)) goto 998
      call poseof(ifi)
      call info(30,0,0,
     .' iobsm0: found initial time from mv file: t = %g',tim,0)
      return

C     Improper format: position ifi to EOF and write initial pos.
  998 continue
      call info(30,0,0,
     .' iobsm0: improper file format; no initial time read',0,0)

C     Reset time; write initial positions to ifi
  999 continue
      tim = 0
      call poseof(ifi)
      call iobsmv(0,tim,T,V,nbas,alat,pos,-ifi)
      call info(30,0,0,
     .' iobsm0: writing initial positions to mv file',0,0)

      end

