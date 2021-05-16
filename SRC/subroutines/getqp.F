      subroutine getqp(mode,nfilqp,nqp,nkabc,lshft,ntet,q,wts,idtet)
      use m_lmfinit,only: stdo
C- I/O for q-points on disc
C-----------------------------------------------------------------------
Ci Input
Ci   mode  : 0  read header information nqp,nkabc,ntet.  File rewound
Ci         :    in write mode only writes header
Ci         : 1  read or write q; wts not written
Ci         :    in write mode, header is always written
Ci         : 2  read or write q,wts, and tet if ntet>0
Ci         :    in write mode, header is always written
Ci   nfilqp, file-handle for QPTS: positive for read, negative for write
Cio Inputs/Outputs
Cio  nqp   : (output, mode 0; input mode 1) number of q-points
Cio  nkabc : (output, mode 0; input mode 1) no. divisions in each
Cio                                         reciprocal lattice vector
Cio  ntet  : (output, mode 0; input mode 1) number of tetrahedra
Co Output
Co   q     : q-points
Co   wts   : degeneracies
Co   idtet : tetrahedron weights
Cr Remarks
Cr   Read branch is separated into modes 0 and 1, so memory can be
Cr   allocated for q,wts,idtet
Cu Updates
Cu   23 Jun 08 New mode; reads and writes header including nkabc
C-----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer mode,nfilqp,nqp,ntet,nkabc(3),lshft(3),idtet(5,*)
      double precision q(3,*),wts(*)
C Local variables
      integer iqp,i,is,iprint,it,iv(5),parg,a2vec
      double precision xv(5)
      external iprint
      character outs*80

      if (nfilqp .lt. 0) goto 100
      outs = ' '
Changenglob      stdo = nglob('stdo')
c      stdo = globalvariables%stdo
C --- Read branch ---
C     Header string: first line not beginning with a '#'
      if (mode .eq. 0) then
        rewind nfilqp
   80   read(nfilqp,'(a80)',end=99,err=99) outs
        if (outs(1:1) .eq. '#') goto 80
C       Read nkabc, if it is present ; else error
        i = 0
        i = parg('nkp=',2,outs,i,len(outs),', ;',3,1,iv,nqp)
        if (i .le. 0)
     .  call rx('GETQP: failed to parse header for "nkp=expr"')
C       Read nkabc, if it is present ; else set to zero
        i = 0
        i = parg('nkabc=',2,outs,i,len(outs),', ;',3,3,iv,nkabc)
        if (i .le. 0) call iinit(nkabc,3)
C       Read lshft, if it is present ; else set to zero
        i = 0
        i = parg('lshft=',2,outs,i,len(outs),', ;',3,3,iv,lshft)
        if (i .le. 0) call iinit(lshft,3)
C       Read ntet, if it is present else set to zero
        i = 0
        i = parg('ntet=',2,outs,i,len(outs),', ;',3,1,iv,ntet)
        if (i .le. 0) ntet = 0
        call info5(30,0,0,
     .  ' GETQP:  disc file contains %i k-points'//
     .  '%?#(n>0)#%-1j  nkabc =%3:1i  lshft =%3:1i#%j#'//
     .  '%?#(n>0)#%-1j  ntet = %i##'//
     .  ' ',nqp,nkabc,lshft,ntet,0)
        return

C --- Error handling ---
   99   continue
        call rx('GETQP: missing or empty file')

      endif

      if (mode .eq. 1) then
        do  iqp = 1, nqp

C         Straight fortran input
C         read (nfilqp,*) i, (q(i,iqp),i=1,3)
C         Enables alegbraic expresions
   90     read(nfilqp,'(a80)',end=98,err=98) outs
          if (outs(1:1) .eq. '#') goto 90
          is = 0
          i = a2vec(outs,len(outs),is,4,', ',2,3,4,iv,xv)
          if (i .ne. 4) goto 90
          do  i = 1, 3
            q(i,iqp) = xv(1+i)
          enddo
        enddo

        if (iprint() .gt. 50) then
          write (stdo,10)
          do  iqp = 1, nqp
            write (stdo,20) iqp, (q(i,iqp),i=1,3)
          enddo
        endif
      elseif (mode .eq. 2) then
        do  iqp = 1, nqp
          read (nfilqp,*) i, (q(i,iqp),i=1,3), wts(iqp)
        enddo
        if (iprint() .gt. 50) then
          write (stdo,10)
          do  iqp = 1, nqp
            write (stdo,20) iqp, (q(i,iqp),i=1,3), wts(iqp)
          enddo
        endif
        if (ntet .gt. 0) then
          do   it = 1, ntet
            read (nfilqp,*) i,(idtet(i,it),i=1,5)
          enddo
        endif
   10   format(19x,'K',16x,'Weight')
   20   format(i4,4f10.6)
      else
        call rx('GETQP: bad mode')
      endif
      return

   98 continue
      call rxi('GETQP: failed to read qp # ', i)

C --- Write branch ---
  100 continue
      i = ntet
      if (mode .lt. 2) i = 0
      call info5(30,0,0,
     .' PUTQP:  writing %i k-points to disc file'//
     .'%?#(n>0)#%-1j  nkabc =%3:1i  lshft =%3:1i#%j#'//
     .'%?#(n>0)#%-1j  ntet = %i##'//
     .' ',nqp,nkabc,lshft,i,0)

      call awrit4('%x   nkp=%i '//
     .'%?#(n>0)#%-1j  nkabc=%3:1i  lshft=%3:1i#%j#'//
     .'%?#(n>0)#%-1j  ntet= %i##'//
     .' ',outs,len(outs),-nfilqp,nqp,nkabc,lshft,i)

      if (mode .eq. 1) then
        do  iqp = 1, nqp
          write(-nfilqp,110) iqp,(q(i,iqp),i=1,3)
        enddo
      elseif (mode .eq. 2) then
        do  iqp = 1, nqp
          write(-nfilqp,110) iqp,(q(i,iqp),i=1,3),wts(iqp)
        enddo
        if (ntet .gt. 0) then
          do  it = 1, ntet
            write(-nfilqp,120) it,(idtet(i,it),i=1,5)
          enddo
        endif
      endif
  110 format(i5,1p,4d20.12)
  120 format(i8,5i5)

      end

