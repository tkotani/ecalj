      logical function parmxp(iter,strn,lstrn,broy,nmix,wgt,beta,
     .  elind,mixnam,wc,nkill,betv,rmserr)

C- Parse strng to get mixing parameters for current iteration
C --------------------------------------------------
Ci Inputs
Ci iter: current iteration
Ci       iter=-1 => strn parsed to check integrity of string,
Ci       and optionally to display (iprint().ge.10)
Ci       parmxp does not set parameters broy...betv in this case
Ci strn: string, parsed to extract mixing parameters for this iteration.
Cio Inputs/Outputs
Cio rmserr used in a conditional expression to determine whether
Cio       mixing mode should shift (see "r<expr", Remarks)
Cio       On output, rmserr is changed to -rmserr if a new mixing block is
Cio       started caused by the condition rmserr<rmsc. This acts as a signal
Cio       in the event that one may wish to kill the mix files.
Co Outputs
Co   Each of mixing parameters broy,nmix,wgt,beta,mixnam,wc,nkill may be
Co   updated depending on the current iteration and strn.  Parameters
Co   not found withing strn take default values (see Remarks).
Co   broy   0 for linear or Anderson mixing
Co          1 for Broyden mixing
Co          2 for conjugate gradients mixing
Co   nmix   number of prior iterations to include in mix
Co   wgt    relative weights to assign; see Remarks
Co   beta   mixing beta (Anderson and CG mixing)
Co   elind  Lindhard screening parameter (when dielectric F can be est)
Co   mixnam file name to read prior iterations
Co   wc     mixing weights for Broyden mixing (see Remarks)
Co   nkill  for periodic mixing file deletion (see Remarks)
Co   betv   for independent potential mixing
Co   parmxp returns false if the parse fails, otherwise true.
Cr Remarks
Cr The general syntax for strn is a sequence of groups of mixing
Cr parameters.  The syntax of one group looks like the following.
Cr   A[nmix][some-parameters--see below][;another-sequence]  or
Cr   B[nmix][some-parameters--see below][;another-sequence]
Cr   thus a ';' indicates the start of a new block in the sequence
Cr   The mixing parameters are as follows:
Cr   ,b=#:   set mixing beta=#.
Cr           (NB: for Broyden mixing, only meaningful to get started)
Cr   ,bv=#[,#2] set extra potential mixing parameter betv to #.  If
Cr           last in this block, set to #2.
Cr   ,n=#    do this mixing for # iterations
Cr   ,k=#    kill the mixing file after # iterations
Cr   ,fn=nam set the mixing file name to 'nam'
Cr   ,wc=#   set Broyden wc to #.  Smaller wc more heavily weights
Cr           most recent iterations.  wc<0 sets wc to abs(wc)*rms error.
Cr   ,w=w1,w2 (spin pol only).  Spin-pol calculations mix up+down and
Cr           up-down.  w1 and w2 are the relative weights to assign
Cr           to these two channels.  w2=0=> magnetic moments frozen;
Cr           w1=0 => total charge is frozen.
Cr   ,wa=#   weight for additonal parameters to be mixed
Cr   ,elind=# Lindhard screening parameter
Cr   ,r<expr continue this block of mixing sequence until rmserr<expr.
Cr           NB: parmxp temporarily loads into the variables table
Cr           the value of errmin, which may be used in parsing expr.
Cr           The latter is set by calling parmx0(0,0,errmin).
Cr           If r<expr is to be used then n=1 must also be set.
Cr Groups are separated by a ";".  parmxp determines which group belongs
Cr   to the current iteration (iter) by adding the number of iterations
Cr   nit for the first group, second group, etc. until their sum exceeds
Cr   the current iteration number (iter).  Thus if strn = "B,n=2;A,n=3",
Cr   the current group would be "B,n=2" for iter = 1,2,6,7,11,12,...
Cr   and would be "A,n=3" for iter = 3,4,5,8,9,10.
Cr Example value of strn : B30,n=8,w=2,1,fn=mxm,wc=11,k=3;A2,b=1
Cr   does 8 iterations of Broyden mixing, followed by Anderson mixing
Cr   The Broyden iterations weight the (up+down) double that of
Cr   (up-down) for the spin pol case, and iterations are saved in a file
Cr   which is deleted at the end of every third iteration.  WC is 11.
Cr   beta assumes the default value.
Cr   The Anderson iterations mix two prior iterations with beta of 1.
Cr Periodic file deletion (nkill):  parmxp returns nkill as -nkill when
Cr   mod(iter-nitj,nkill) is zero, as a signal that the current mixing
Cr   file is to be deleted after it is used.  Here nitj is the sum of
Cr   the number of iterations prior to the current group.
Cb Bugs:
Cb   parmxp cannot tell if this iterations is the last one in a block
Cb   if the constraint rmsc<rmserr is not satisfied.  If it is not,
Cb   betv always returns bv(1).
Cr Defaults:
Cr   nkill has hardwired defaults of -1 and 0.
Cr   nit defaults to infinity
Cr   broy,nmix,wgt,beta,mixnam,wc,betv default to their input values.
Cr   Setting input wgt(3) to -9 signals that wgt(3) is not used.
Cl Local variables
Cl   iblk:  index to current sequence of mixing parameters
Cl   nitj:  effective iter corresponding to last iteration of last
Cl          mixing block, to determine position in current mixing block
Cl   lstblk:index to block of mixing parameters of last call.  Normally
Cl          set and saved internally.  Caller may set lstblk (entry
Cl           parmx0) to fix the current place in the mixing block group.
Cl   lstitj:performs a dual function, depending on sign.
Cl          lstitj<=0: corresponds to -nitj of prior call.  parmxp sets
Cl                     lstitj internally for this mode.
Cl          lstitj>0:  (set by caller), which iteration within the
Cl                     current mixing block the next iteration will
Cl                     correspond to.
Cu Updates
Cu    20 Jul 08 (ATP) returns rmserr as -rmserr to signal new block
Cu    1  Jun 00 added argument elind
C --------------------------------------------------
C     implicit none
C Passed parameters
      integer broy
C     character strn*(*),mixnam*8
      character strn*(70),mixnam*8
      integer iter,lstrn,nkill,nmix
      double precision wgt(3),beta,elind,wc,betv,rmserr
C Local variables
      integer i,j,np,it(5),parg,nit,nmixj,jp,kp,killj,nitj,
     .  iprint,i1mach,a2vec,lstblk,lstitj,iblk,nbump,k,ia0,lbroy
      logical lpr,lagain,cmdopt
      character outs*100,fnam*8
      double precision bet,elin,wt(3),wcj,rmsc,bv(2),errmin,xx
C ... this is the only way to create static variables in fortran
      common /parmx1/ lstblk,lstitj,errmin

      parmxp = .true.
      ia0 = -1
      call bin2a0(ia0)
C ... Internal defaults
      if (iter .gt. 0) then
        nkill = 0
        nit = -1
      endif
      if (strn .eq. ' ' .or. lstrn .le. 0) goto 9999
C ... Passed defaults
      call bin2a0(10)
      fnam = mixnam
      lbroy = broy
      wcj = wc
      nmixj = nmix
      bet = beta
      elin = elind
      bv(1) = betv
      bv(2) = betv
      wt(1) = wgt(1)
      wt(2) = wgt(2)
      wt(3) = wgt(3)
      nit = -1
      if (wgt(3) .eq. -9) wt(3) = 0
C ... lpr: switch to determine whether to print or not
C     lpr = .false.
C ... iblk,iit: current mixing block and iter within block
      iblk = 0
      nitj = 0
      np = 0
      outs = ' mixing: mode=A'

C --- Entry point for parsing a new set of switches ---
   10 continue
      call skipbl(strn,lstrn,np)
      if (np .ge. lstrn) then
        if (iter .lt. 0) goto 9999
        if (cmdopt('--nomixcycle',12,0,outs)) then
          nkill = -999
          goto 9999
        else
          np = 0
          goto 10
        endif
      endif
C ... Switch for Broyden or Anderson mixing
      jp = np
      call chrps2(strn,'AaBbCc',6,np,jp,it)
      if (it(1) .eq. 0) goto 999
      iblk = iblk+1
C ... If iblk=lstblk, override nitj with lstitj
      if (iter .gt. 0 .and. iblk .eq. lstblk) then
        nitj = -lstitj
        if (lstitj .gt. 0) nitj = iter - lstitj
        lstitj = -nitj
C  ...  A bug if starting iteration bigger than iter
        if (nitj .gt. iter) call rx('parmxp: bad lstitj')
      endif
      lbroy = 0
      if (it(1) .ge. 3) lbroy = 1
      if (it(1) .ge. 5) lbroy = 2
      if (lbroy .eq. 1) call awrit0('%a%bB',outs,len(outs),0)
      if (lbroy .eq. 2) call awrit0('%a%bC',outs,len(outs),0)
C ... Pick up nmix
      jp = np+1
      call chrps2(strn,',; ',3,np+1,jp,it)
      if (it(1) .eq. 0) then
        if (a2vec(strn,lstrn,jp,2,',; ',3,1,1,it,nmixj) .lt. 0) goto 999
        call awrit1('%a  nmix=%i',outs,len(outs),0,nmixj)
      endif
C ... Pick up rmsc
      rmsc = -1
      jp = np
C ... set variable errmin to current value of errmin
      call getsyv('errmin',xx,j)
      call lodsyv('errmin',1,errmin,k)
      i = parg(',r<',4,strn,jp,lstrn,',; ',2,1,it,rmsc)
C ... Put back the original one, or remove newly created one
      if (j .eq. k) then
        call lodsyv('errmin',1,xx,k)
      else
        call clrsyv(k-1)
      endif
      if (i .lt. 0) goto 999
C     if (i .gt. 0) lpr = .true.
      if (rmsc .ge. 0 .and. iter .lt. 0)
     .  call awrit1('%a  err<%1;3g',outs,len(outs),0,rmsc)
      if (rmsc .ge. 0 .and. iter .gt. 0)
     .  call awrit2('%a  err(%1;3g)<%1;3g',outs,len(outs),0,rmserr,rmsc)
C ... Pick up nit
      jp = np
      i = parg(',n=',2,strn,jp,lstrn,',; ',2,1,it,nit)
      if (i .lt. 0) goto 999
C ... increment nit if rmserr>rmsc and iter>=nit+nitj
      nbump = 0
      if (nit .ne. -1) then
        if (iter .lt. 0) call awrit1('%a  nit=%i',outs,len(outs),0,nit)
        if (iter .gt. 0) then
          call awrit2('%a  it %i of %i',outs,len(outs),0,iter,nit+nitj)
          if (iblk .ge. lstblk .and. iter .ge. nit+nitj .and.
     .      rmserr .gt. rmsc .and. rmsc .gt. 0) then
            nbump = iter - (nit+nitj) + 1
            nit = nit + nbump
            call awrit0('%a(*)',outs,len(outs),0)
          endif
        endif
      endif

C ... ATP added this:
      if (rmserr .lt. rmsc .and. iblk .eq. lstblk .and. rmserr .ne. 0)
     .  rmserr = -rmserr

C      if (nit .eq. -1 .and. iter .gt. 0)
C     .  call awrit1('%a  it %i of *',outs,len(outs),0,iter)

C     if (i .gt. 0) lpr = .true.
C ... Pick up file name
      jp = np
      i = parg(',fn=',0,strn,jp,lstrn,',; ',2,0,it,0)
      if (i .gt. 0) then
        kp = jp+1
        call chrps2(strn,',; ',3,jp+5,kp,it)
        fnam = strn(jp+1:kp)
        call awrit0('%a  fnam='//fnam,outs,len(outs),0)
C       lpr = .true.
      endif
C ... Pick up mixing wc
      jp = np
      if (lbroy .eq. 1) then
        i = parg(',wc=',4,strn,jp,lstrn,',; ',2,1,it,wcj)
        if (i .lt. 0) goto 999
C       if (i .gt. 0) lpr = .true.
        if (i .gt. 0) call awrit1('%a  wc=%d',outs,len(outs),0,wcj)
      endif
C ... Pick up mixing beta
      jp = np
      i = parg(',b=',4,strn,jp,lstrn,',; ',2,1,it,bet)
      if (i .lt. 0) goto 999
C     if (i .gt. 0) lpr = .true.
      call awrit1('%a  beta=%d',outs,len(outs),0,bet)
C ... Pick up elind
      jp = np
      i = parg(',elind=',4,strn,jp,lstrn,',; ',2,1,it,elin)
      if (i .lt. 0) goto 999
C     if (i .gt. 0) lpr = .true.
      if (elin.ne.0) call awrit1('%a  elind=%;3d',outs,len(outs),0,elin)
C ... Pick up weights
      jp = np
      i = parg(',w=',4,strn,jp,lstrn,',; ',2,2,it,wt)
      if (i .lt. 0) goto 999
      jp = np
      j = parg(',wa=',4,strn,jp,lstrn,',; ',2,1,it,wt(3))
      if (j .lt. 0) goto 999
      if (i .gt. 0 .or. j .gt. 0) then
        call awrit2('%a  wgt=%d,%d',outs,len(outs),0,wt(1),wt(2))
        if (j .gt. 0 .and. wgt(3) .ne. -9) then
          call awrit1('%a(%d)',outs,len(outs),0,wt(3))
        elseif (j .gt. 0) then
          call awrit0('%a(-)',outs,len(outs),0)
        endif
      endif
C...  Pick up iteration number for file kill
      killj = -1
      jp = np
      i = parg(',k=',2,strn,jp,lstrn,',; ',2,1,it,killj)
      if (i .lt. 0) goto 999
C     if (i .gt. 0) lpr = .true.
      if (killj .ne. -1) call awrit1('%a  kill=%i',outs,len(outs),
     .  0,killj)
C...  Pick up betv
      jp = np
      i = parg(',bv=',4,strn,jp,lstrn,',; ',2,2,it,bv)
C     if only one element found, copy 1st element to second:
      if (i .eq. -2) then
        bv(2) = bv(1)
        i = 1
      endif
      if (i .lt. 0) goto 999
      if (i .gt. 0) lpr = .true.
      if (iter .eq. nitj+nit .and. iblk .ge. lstblk) bv(1) = bv(2)
      if (bv(1).ne.1) call awrit1('%a  betv=%1;3g',outs,len(outs),0,bv)
C ... If iter < 0, printout and parse through all strings
      if (iter .lt. 0) then
        if (iprint().ge.10) call awrit0('%a',outs,-len(outs),-i1mach(2))
        lagain = nit .ne. -1
        outs = '         mode=A'
      endif

C --- If this is last pass, eg nitj <= iter <nitj+nit ---
      if (iter .gt. 0 .and. iblk .ge. lstblk .and. nitj .lt. iter
     .  .and. (iter .le. nitj+nit .or. nit .eq. -1)) then
        if (iprint().ge.20) call awrit0('%a',outs,-len(outs),-i1mach(2))
        broy = lbroy
        nmix = nmixj
        wgt(1) = wt(1)
        wgt(2) = wt(2)
        if (wgt(3) .eq. -9) wt(3) = 0
        wgt(3) = wt(3)
        beta = bet
        elind = elin
        mixnam = fnam
        wc = wcj
        nkill = max(killj,0)
        if (nkill .gt. 1) then
          if (mod(iter-nitj,nkill) .eq. 0) nkill=-nkill
        endif
        if (nkill .eq. 1) nkill=-nkill
        betv = bv(1)
        lstblk = iblk
        lstitj = -nitj
        if (nbump .gt. 0) lstitj = lstitj - (nbump-1)
C       print *, 'exiting', lstitj,nbump
        goto 9999
      else
        nitj = nitj+nit
        lagain = .true.
        outs = ' mixing: mode=A'
      endif

C  99 continue
      if (lagain) then
        call chrps2(strn,'; ',2,lstrn,np,it)
        np = np+1
        goto 10
      else
        goto 9999
      endif

C --- Error exit ---
  999 outs = 'parmxp: parse failed:'//strn(1:lstrn)
      if (iprint() .ge. 10) call awrit0('%a',outs,-len(outs),-i1mach(2))
      parmxp = .false.
C --- Normal exit ---
 9999 continue
      call bin2a0(ia0)
Cgetarg       end
      end function parmxp

      subroutine parmx0(i1,i2,errxx)

C- sets lstblk,lstitj and errmin
C     implicit none
      integer i1,i2,lstblk,lstitj,mode
      double precision errxx,errmin
C ... this is the only way to create static variables in fortran
      common /parmx1/ lstblk,lstitj,errmin
      if (i1 .gt. 0) lstblk = i1
      if (i2 .gt. 0) lstitj = i2
      if (errxx .ge. 0d0) then
        if (errmin .gt. 0) errmin = min(errmin,errxx)
        if (errmin .eq. 0) errmin = errxx
      endif
C     print *, 'errmin=',errmin
      return
      entry parms0(i1,i2,errxx,mode)
      if (mode .gt. 0) then
        i1 = lstblk
        i2 = lstitj
        errxx = errmin
      else
        lstblk = i1
        lstitj = i2
        errmin = errxx
      endif
Cgetarg       end
      end subroutine parmx0

      block data dparmx
      integer lstblk,lstitj
      double precision errmin
      common /parmx1/ lstblk,lstitj,errmin
      data lstblk /0/ lstitj /0/ errmin /0d0/
Cgetarg       end
      end block data dparmx

c ... testing parmxp
C      subroutine fmain
C      implicit none
C      character*90 strn,mixnam*8
C      integer i,ip,it(100),iter,nkill,nmix,i1,i2
C      double precision wgt(3),beta,betv,wc,rmserr
C      real ran1
C      logical lbroy,parmxp
C
C      nmix = -1
C      iter = -1
C      wc = -1
C      beta = .7d0
C      betv = .01d0
C      wgt(1) = 1
C      wgt(2) = 1
C      wgt(3) = 1
C      rmserr = -1
C      mixnam = 'MiXm'
C      strn = ' B30,n=4,w=2,1,wa=9,fn=mxm,wc=11,b=1; A,k=3,bv=.11,.22'
C      if (.not. parmxp(iter,strn,len(strn),lbroy,nmix,wgt,beta,mixnam,
C     .    wc,nkill,betv,rmserr)) stop
C      print *, '-----'
C      do  10  iter = 3, 12
C      if (.not. parmxp(iter,strn,len(strn),lbroy,nmix,wgt,beta,mixnam,
C     .    wc,nkill,betv,rmserr)) stop
C      print 333, iter,lbroy,nmix,wgt,beta,mixnam,wc,nkill,betv,rmserr
C  333 format('iter=',i3,' lbroy=',l1,' nmix=',i2,' wgt=',3f8.4,' beta=',
C     . f8.4,' mixnam=',a,' wc=',f8.4,' nkill=',i3,' betv,rmserr=',2f8.4)
C      print *, ' '
C   10 continue
C
C      call parmx0(2,2,0d0)
C
C      strn = ' B30,n=4,w=2,1,wa=9,fn=mxm,wc=11,b=1;'//
C     .  ' A,k=3,bv=.11,.22,n=3;A,b=.88,r<.2,n=2,bv=.11,.33'
C      strn = ' B30,n=4,w=2,1,wa=9,fn=mxm,wc=11,b=1;'//
C     .  ' A,k=3,bv=.11,.22,n=3;A,b=.88,r<errmin+.2,n=2,bv=.11,.33'
C      if (.not. parmxp(-1,strn,len(strn),lbroy,nmix,wgt,beta,mixnam,
C     .    wc,nkill,betv,rmserr)) stop
C      call ran1in(1)
C      print *, '-----'
C      do  20  iter = 3, 18
C      rmserr = ran1()
C      if (.not. parmxp(iter,strn,len(strn),lbroy,nmix,wgt,beta,mixnam,
C     .    wc,nkill,betv,rmserr)) stop
C      call parmx0(0,0,rmserr)
C      print 333, iter,lbroy,nmix,wgt,beta,mixnam,wc,nkill,betv,rmserr
C      print *, ' '
C   20 continue
C
CC      strn = ' '
CC      read(*,'(a70)') strn
CC      if (.not. parmxp(iter,strn,len(strn),lbroy,nmix,wgt,beta,mixnam,
CC     .  wc,nkill,betv,rmserr)) stop
C
CC      call parms0(i1,i2,rmserr,1)
CC      print *, i1,i2,rmserr
CC      i1 = 8
CC      i2 = -15
CC      rmserr = .1234d0
CC      call parms0(i1,i2,rmserr,-1)
CC      i1 = 0
CC      i2 = 0
CC      rmserr = 0
CC      call parms0(i1,i2,rmserr,1)
CC      print *, i1,i2,rmserr
C
C      end

