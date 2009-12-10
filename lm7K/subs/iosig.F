      integer function iosig(ismod,ifis,iq,nq,ldim,sigma)
C- File I/O for sigma
C ----------------------------------------------------------------------
Cio Inputs/Outputs
Ci    ifis >0 to read, <0 to write
Ci    iq   -1 (ifis>0 only) read header into 
Ci          0, read/write header; input checks match with
Ci             passed ismod,ldim,nq
Ci         >0, read/write sigm(iq)
Ci    nq    number of qp.  On input, used to check; on output, in header
Cio   sigma,ldim   sigma read or written for this iq, depending on mode;
Cio         see Remarks
Cio   iosig returns -1 if failed to read sigma (iq>0) or
Cio         if a file mismatch (iq=0)
Cr  Remarks
Cr    File consists of a header containing ismod, ldim and nq
Cr    followed by Sigma.  For iq>0 iosig returns sigma for 1 qp only:
Cr    *ismod=0  Sigma for every qp is written to disk; on writing iosig
Cr     is called nq+1 times, one for the header and one for each iq.
Cr    *ismod=1  Sigma for 1st qp only is written to disk; on reading
Cr     iosig rewinds the file and reads the same Sigma each time.
Cr
Cr    On reading, iq=0, iosig checks for a match in ldim,
Cr    and if ismod=0, requires file has ismod=0 and checks nq.
C ----------------------------------------------------------------------
C     implicit none
      integer ismod,ifis,iq,nq,ldim
      double precision sigma(ldim,ldim,2)
      integer ismod0,ldim0,nq0

      iosig = 0

C --- File read ---
      if (ifis .gt. 0) then
C   ... Header
        if (iq .eq. -1) then
          rewind ifis
          read (ifis,err=999,end=999) ismod,ldim,nq
        elseif (iq .eq. 0) then
          rewind ifis
          read (ifis,err=999,end=999) ismod0,ldim0,nq0
          if (ldim .ne. ldim0) goto 999
          if (ismod .eq. 0) then
            if (nq .ne. nq0 .or. ismod0 .ne. 0) goto 999
          endif
          return
C   ... Sigma
        else
          if (ismod .eq. 0) then
            read (ifis,err=999,end=999) sigma
          else
            rewind ifis
            read (ifis,err=999,end=999) ismod0
            read (ifis,err=999,end=999) sigma            
          endif
        endif
C --- File write ---
      else
C   ... header
        if (iq .eq. 0) then
          rewind (-ifis)
          write (-ifis) ismod,ldim,nq
C   ... Sigma
        else
          write (-ifis) sigma
        endif
      endif
      return
  999 iosig = -1
      return

      end

