c$$$      module m_pvioeu
c$$$      contains
c$$$      subroutine pvioeu(mode,ssite,eula,nbas,neul)
c$$$      use m_struc_def  
c$$$C- Writes the Euler angles to a file or to stdout
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   mode  :0 write in format suitable for rdm
c$$$Ci         :1 write in printout format
c$$$Ci         :10s digit
c$$$Ci         :0 printout is for Euler angles
c$$$Ci         :1 printout is for Bfield
c$$$Ci   ssite :struct for site-specific information; see routine usite
c$$$Ci     Elts read: *
c$$$Ci     Stored:    clabel
c$$$Ci   eula  :Euler angles for noncollinear spins
c$$$Ci   nbas  :size of basis
c$$$Ci   neul  :1 if Euler angles are l-independent, nl otherwise
c$$$Co Outputs
c$$$Co   Writes to stdout the Euler angles
c$$$Cr Remarks
c$$$Cr
c$$$Cu Updates
c$$$Cu   14 Feb 03 output can be for euler angles or for B-field
c$$$Cu   22 May 02 Writes label information if passed through in ssite
c$$$Cu             Altered argument list
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      integer mode,nbas,neul
c$$$      real(8):: eula(nbas,neul,3)
c$$$      type(s_site),optional::ssite(*)
c$$$C ... Local parameters
c$$$      integer i,j,k,lgunit,stdo,mode0,mode1,lbl
c$$$      double precision rotm(3,3)
c$$$      character*72 outs, clabl*8
c$$$      character*255:: formatc(344:346),lbll
c$$$      mode0 = mod(mode,10)
c$$$      mode1 = mod(mode/10,10)
c$$$      stdo = lgunit(1)
c$$$      clabl = ' '
c$$$      if (mode0 .ne. 0) then
c$$$        lbl=344
c$$$        if (mode1 .ne. 0) lbl=345
c$$$        if (mode1 .eq. 0 .and. mode0 .ne. 0) lbl=346
c$$$        formatc(344)="('   ib     alpha        beta       gamma':6x,a)"
c$$$        formatc(345)="('   ib      bx          by          bz  ':6x,a)"
c$$$        formatc(346)="('   ib   alpha      beta     gamma':5x,a,"//
c$$$     .  "11x,'z axis rotated to:')"
c$$$        lbll=formatc(lbl)
c$$$
c$$$        if ( neul .gt. 1)  write(outs,trim(lbll))
c$$$c        if (ssite(1)%size .eq. 0 .or. neul .gt. 1)  write(outs,trim(lbll))
c$$$c        if (ssite(1)%size .ne. 0 .and. neul .eq. 1) write(outs,trim(lbll)) 'class'
c$$$      endif
c$$$c  344 format('   ib     alpha        beta       gamma':6x,a)
c$$$c  345 format('   ib      bx          by          bz  ':6x,a)
c$$$c  346 format('   ib   alpha      beta     gamma':5x,a,
c$$$c     .  11x,'z axis rotated to:')
c$$$      if (neul .gt. 4) then
c$$$        call awrit0('%a     lm',outs,len(outs),0)
c$$$      elseif (neul .gt. 1) then
c$$$        call awrit0('%a      l',outs,len(outs),0)
c$$$      endif
c$$$      call awrit0(outs,' ',-len(outs),stdo)
c$$$      do  21  i = 1, nbas
c$$$c        if (ssite(1)%size .ne. 0) then
c$$$c          do i_spacks=i,i
c$$$c            call spacks_copy('u',ssite(i_spacks)%clabel,i,i,clabl,i_spacks)
c$$$c          enddo
c$$$c        endif
c$$$        if (neul .gt. 1 .and. mode0 .ne. 0) then
c$$$          if (clabl .eq. ' ') call awrit1('# site %,4i',' ',80,stdo,i)
c$$$          if (clabl .ne. ' ') call awrit1('# site %,4i  class '//clabl,
c$$$     .    ' ',80,stdo,i)
c$$$          write(stdo,'(i5,3f12.6,i4)')
c$$$     .    (i,(eula(i,j,k),k=1,3),j, j=1,neul)
c$$$        elseif (neul.eq.1 .and. mode0.ne.0 .and. clabl.eq.' ') then
c$$$          write(stdo,'(i5,3f12.6)')
c$$$     .    (i,(eula(i,j,k),k=1,3), j=1,neul)
c$$$        elseif (neul .eq. 1 .and. mode0 .ne. 0) then
c$$$          call eua2rm(eula(i,1,1),eula(i,1,2),eula(i,1,3),rotm)
c$$$          write(stdo,'(i5,3f10.6,3x,a,1x,3f10.6)')
c$$$     .    i,(eula(i,1,k),k=1,3), clabl, (rotm(3,k),k=1,3)
c$$$        else
c$$$          write(stdo,'(3f16.12)') ((eula(i,j,k),k=1,3),j=1,neul)
c$$$        endif
c$$$   21 continue
c$$$
c$$$      end subroutine pvioeu
c$$$      end module m_pvioeu
c$$$
c$$$      subroutine ioeula(nbas,nl,eula,neul,xsi,ifi)
c$$$C- I/O of Euler angles
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   ifi: <0 for write, >0 for read
c$$$Ci   nbas: number of basis atoms
c$$$Co Outputs
c$$$Co   Euler angles are read in or written out
c$$$Co   neul  is read in or written out (either 1, nl or nl**2)
c$$$Co   xsi is read in if xsi=# is available on the first line
c$$$Co   xsi is written out if xsi is nonzero
c$$$Cr Remarks
c$$$Cr   Aborts on read when nbas does not match file
c$$$Cu Updates
c$$$Cu   24 May 08 Only read/write master (MPI)
c$$$Cu             Note: it is the caller's responsibility to broadcast
c$$$Cu             This is because neul and the dimensions of eula
c$$$Cu             are only known after reading.  Typical Broadcast:
c$$$Cu             call mpibc1(neul,1,2,mlog,'susite','neul')
c$$$Cu             call mpibc1(xsi,5,4,mlog,'susite','xsi')
c$$$Cu             call mpibc1(eula,nbas*neul*3,4,mlog,'susite','eula')
c$$$Cu   19 Nov 97 use rdfiln for reading first line
c$$$C ----------------------------------------------------------------
c$$$      use m_pvioeu
c$$$C     implicit none
c$$$      integer nbas,neul,nl,ifi
c$$$      double precision eula(nbas,neul,3),tmp(125),xsi(3)
c$$$      integer i,j,k,ipr,lgunit,nc,rdm,nxsi,a2vec,ix(10)
c$$$      logical parstr,a2bin
c$$$      character*80 ss,ss2*4
c$$$      integer master,mpipid
c$$$C ... for rdfiln
c$$$      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
c$$$      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
c$$$      character recrd*(recl),ctbl(mxchr,2)*(ctlen),a*(recl),
c$$$     .vnam(mxlev)*16,rdarg*6
c$$$      logical loop0(0:mxlev)
c$$$      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
c$$$     .nlist(0:mxlev)
c$$$      equivalence (a,ss)
c$$$      data rdarg /'#{}% c'/
c$$$
c$$$
c$$$C --- I/O through master node only (MPI) ---
c$$$      master = 0
c$$$      if (mpipid(1) .eq. master) then
c$$$
c$$$C --- Write ---
c$$$        if (ifi .lt. 0) then
c$$$C       rewind (-ifi)
c$$$          if (neul .ne. 1 .and. neul .ne. nl .and. neul .ne. nl*nl)
c$$$     .    call rxi('ioeula: bad dim for eula: neula=',neul)
c$$$          call awrit3('%% rows %i cols %i neula %i',ss,80,0,nbas*neul,3,
c$$$     .    neul)
c$$$          nxsi = 0
c$$$          do  11  i = 1, 3
c$$$            if (xsi(i) .eq. 0) goto 12
c$$$            nxsi = nxsi+1
c$$$   11     continue
c$$$   12     continue
c$$$          if (nxsi .gt. 0)
c$$$     .    call awrit3('%a nxsi=%i %n:1;6,6d',ss,80,0,nxsi,nxsi,xsi)
c$$$          call awrit0('%a',ss,80,ifi)
c$$$          do  10  i = 1, nbas
c$$$            if (neul .gt. 1) call awrit1('# ib %,4i',' ',80,-ifi,i)
c$$$            write(-ifi,'(3f16.12)') ((eula(i,j,k),k=1,3),j=1,neul)
c$$$   10     continue
c$$$C --- Read ---
c$$$        else
c$$$          call getpr(ipr)
c$$$          rewind ifi
c$$$          nr = 0
c$$$          call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
c$$$     .    ilist,nlist,vnam,ctbl,mxchr,ss,recrd,recl,nr)
c$$$C   ... No error if file is empty, but don't bother reading
c$$$          if (nr .eq. 0) then
c$$$            if (ipr .gt. 30) print '('' IOEULA:   empty file'')'
c$$$            return
c$$$          endif
c$$$          if (ss(1:1) .eq. '%') then
c$$$C   ... Read xsi if it is there
c$$$            i = 0
c$$$            if (parstr(ss,'nxsi=',len(ss)-5,5,' ',i,j)) then
c$$$              j = j-1
c$$$              i = j
c$$$              if (a2bin(ss,nxsi,2,0,' ',j,len(ss)-5)) then
c$$$                j = a2vec(ss,len(ss)-5,i,4,' #',2,2,nxsi+1,ix,tmp)
c$$$                if (j .ne. nxsi+1) call rx('ioeula:  failed to read xsi')
c$$$                call dcopy(nxsi,tmp(2),1,xsi,1)
c$$$              endif
c$$$            endif
c$$$            neul = 1
c$$$            i = 0
c$$$            if (parstr(ss,'neula ',len(ss)-6,6,' ',i,j)) then
c$$$              if (.not. a2bin(ss,neul,2,0,' ',j,len(ss)))
c$$$     .        call rx('ioeula: failed to parse for neula')
c$$$            endif
c$$$          else
c$$$            neul = 1
c$$$            rewind ifi
c$$$          endif
c$$$
c$$$C   ... Expect nbas*neul*3 rows in file ... skip reading if not
c$$$          if (neul .ne. 1 .and. neul .ne. nl .and. neul .ne. nl*nl)
c$$$     .    call rxi('ioeula: illegal dimension: neula =',neul)
c$$$          nc = 3
c$$$          ss2 = 'in'
c$$$          if (neul .eq. nl) ss2 = ' '
c$$$          if (neul .eq. nl*nl) ss2 = '%bm-'
c$$$          if (ipr .gt. 30) call awrit1(' IOEULA:  reading l-'//ss2//
c$$$     .    '%adependent Euler angles (file %i)',' ',80,lgunit(1),ifi)
c$$$C        if (nr .ne. neul*nbas .or. nc .ne. 3) then
c$$$C          if (ipr .ge. 20) call awrit3(' ioeula: input skipped --'//
c$$$C     .      ' expected %i rows (3 cols) but found %i (%i cols)',
c$$$C     .      ' ',80,lgunit(1),neul*nbas,nr,nc)
c$$$
c$$$C   ... Read the data, 1 atom at a time ...
c$$$          call dpzero(eula,nbas*neul*3)
c$$$          do  20  i = 1, nbas
c$$$            call rxx(rdm(ifi,0,neul*nc,' ',tmp,neul,nc).ne.1,
c$$$     .      'ioeula:  bad file')
c$$$            call dmcpy(tmp,neul,1,eula(i,1,1),nbas*neul,nbas,neul,nc)
c$$$   20     continue
c$$$          if (ipr .ge. 50) call pvioeu(mode=1,eula=eula,nbas=nbas,neul=neul)
c$$$        endif
c$$$      endif
c$$$      end subroutine ioeula
