      subroutine pqedit(nclass,nl,lmx,nsp,nx,iter,mixmod,dmxp,
     .  pold,qold,xold,pnu,qnu,xnew)
C- Edits mixing file
C ------------------------------------------------------------------
Ci Inputs: see pqmix
Cr Remarks
Cr   switches: -fn=filenam overrides default filenam (from ctrl file)
Cr             -a2bin converts ascii file to binary
Cr             -bin2a converts binary file to ascii
C ----------------------------------------------------------------
C     implicit none
C Passed variables:
      logical broy
      character*(*) mixmod
      integer nclass,nl,nsp,nx,lmx(nclass),iter
      double precision dmxp(33),
     .  pnu(nl,nsp,nclass), qnu(3,nl,nsp,nclass), xnew(nx),
     .  pold(nl,nsp,nclass),qold(3,nl,nsp,nclass),xold(nx)
C Local variables
C mmix: max no. prior iterations available from disk.
      logical parmxp,cmdopt,lbin
      integer mmix,namx,ic,l,nelts,idum,fopna,mxsav,nmix,
     .  nkill,lgunit,iprint,ifi,oa,oaa,i0,npq
      character outs*80,fnam*8
      double precision wt(3),rmsdel,rms2,wc,beta,elind
C heap:
      integer w(1)
      common /w/ w

C     call wkprnt(1)
      if (iprint() .ge. 20) print *, ' '

C --- Iteration-dependent mixing parameters ---
      broy  = dmxp(1) .ne. 0
      beta  = dmxp(2)
      wc    = dmxp(3)
      wt(1) = dmxp(4)
      wt(2) = dmxp(5)
      wt(3) = 0
      mxsav = nint(dmxp(6))
      nmix  = nint(dmxp(7))
      nkill = nint(dmxp(8))
      fnam  = 'mixm'
      rmsdel = dmxp(11)
      elind = dmxp(33)
      if (.not. parmxp(iter,mixmod,len(mixmod),broy,nmix,wt,beta,elind,
     .  fnam,wc,nkill,dmxp(9),dmxp(10))) call rx(
     .  'PQEDIT: parse in parmxp failed')
      if (wt(1)**2+wt(2)**2+wt(3)**2 .eq. 0)
     .  call fexit(-1,111,' Exit -1 PQEDIT: '//
     .  'bad mixing weights w =%3:1;6d',wt)
      if (cmdopt('-fn=',4,0,outs)) fnam = outs(5:80)
      call dscal(3,1/dsqrt(wt(1)**2+wt(2)**2+wt(3)**2),wt,1)

C --- Count number of mixing elements ---
      npq = 0
      rms2 = 0
      do  10  ic = 1, nclass
      do  10  l = 0, lmx(ic)
   10 npq = npq+1
      npq = npq*4
      namx = npq*nsp
      namx = namx + nx

C --- Determine number of previous iterations available from disk ---
      mmix = 0
      if (cmdopt('-a2bin',6,0,outs)) then
        ifi = fopna('mixa',-1,0)
        rewind ifi
        read(ifi,*,err=99,end=99) idum, nelts
        lbin = .false.
      else
        ifi = fopna(fnam,-1,4)
        rewind ifi
        read(ifi,err=99,end=99) idum, nelts
        lbin = .true.
      endif
      if (nelts .ne. namx) then
        call awrit2(' PQEDIT: expecting %i elements but found %i ... '//
     .    'discarding file',' ',80,lgunit(1),namx,nelts)
        goto 30
      endif
      mmix = idum
   30 continue

      call awrit2(' pqedit:  read %i iter, %i elements',outs,80,0,
     .  mmix,namx)
      if (nx .gt. 0) call awrit2('%a (%i+%i)',outs,80,0,namx-nx,nx)
      call awrit0('%a from file '//fnam//'%a',outs,80,-lgunit(1))

C --- Copy new P,Q into holding array, read prior iter from disk ---
      call defdr(oa,-namx*(mxsav+2)*2)
      call defdr(oaa,-namx*(mxsav+2)*2)
      stop 'update this call'
C      call pqmxio(namx,mmix,mxsav,ifi,lbin,nclass,nl,nsp,nx,lmx,
C     .  pnu,pold,qnu,qold,xnew,xold,w(oa),rms2)

C --- Copy a subset of the iterations found into holding array ---
      i0 = 1
   22 nmix = 0
      call pqed2(namx,mxsav,i0,mmix,w(oa),nmix,w(oaa))
      if (nmix .eq. -1) then
        call dcopy(namx*(mxsav+2)*2,w(oa),1,w(oaa),1)
        nmix = mmix
      endif
      if (nmix .ne. mmix) then
        mmix = nmix
        call dcopy(namx*(mxsav+2)*2,w(oaa),1,w(oa),1)
        i0 = 0
        goto 22
      endif

C --- Save updated P,Q on disk ---
      if (cmdopt('-bin2a',6,0,outs)) then
        call fclose(ifi)
        ifi = fopna('mixa',-1,0)
        rewind ifi
        write(ifi,'(2i6)') min(mmix+1,mxsav), namx
        lbin = .false.
      else
        if (cmdopt('-a2bin',6,0,outs)) then
          call fclose(ifi)
          ifi = fopna(fnam,-1,4)
        endif
        rewind ifi
        write(ifi) min(mmix+1,mxsav), namx
        lbin = .true.
      endif
      stop 'update this call'
C      call pqmxio(namx,mmix+1,mxsav,-ifi,lbin,nclass,nl,nsp,nx,lmx,
C     .  pnu,pold,qnu,qold,xnew,xold,w(oa),rms2)
      call fclose(ifi)

C --- Copy mixed P,Q from holding array ---
      call pqmxup(namx,mxsav,nclass,nl,nsp,nx,lmx,pnu,qnu,xnew,
     .  pold,qold,xold,0,namx,w(oa),rms2)

      return

   99 call rx('PQEDIT: read error in mixing file ' // fnam)
      end

      subroutine pqed2(namx,mxsav,i0,mmix,a,nmix,aa)
C- Retain a subset of iterations from mixing file
Ci  i0: starting iteration in a
Cio nmix: starting iteration in aa; nmix total number in aa
C     implicit none
      integer namx,mxsav,mmix,i0,nmix
      double precision a(namx,0:mxsav+1,2),aa(namx,0:mxsav+1,2)
      integer i1mach,i,nlist,list(100),j
      double precision rms2,ddot
C Heap
      character*80 s

      do  10  i = i0, mmix
        rms2 =  dsqrt(dabs(ddot(namx,a(1,i,1),1,a(1,i,1),1)
     .    -2*ddot(namx,a(1,i,1),1,a(1,i,2),1)
     .    + ddot(namx,a(1,i,2),1,a(1,i,2),1))/namx)
        call awrit2('  iter %i   RMS DQ = %1,4;4e',' ',80,
     .    i1mach(2),i,rms2)
   10 continue

      call cwrite(' Prior iterations to keep: ',0,34,0)
      read(*,'(a80)') s
      if (s .eq. ' ') call awrit2('%i:%i',s,80,0,i0,mmix)
      call mkilst(s,nlist,list)

      do  20  i = 1, nlist
        j = list(i)
        if (j .lt. i0 .or. j .gt. mmix) then
          call awrit1('  iter %i out of range ...',' ',80,i1mach(2),j)
        else
          call dcopy(namx,a(1,j,1),1,aa(1,nmix,1),1)
          call dcopy(namx,a(1,j,2),1,aa(1,nmix,2),1)
          nmix = nmix+1
        endif
   20 continue
      nmix = nmix-1
      end

