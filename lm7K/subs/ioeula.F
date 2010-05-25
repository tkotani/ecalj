      subroutine ioeula(nbas,nl,eula,neul,xsi,ifi)

C- I/O of Euler angles
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi: <0 for write, >0 for read
Ci   nbas: number of basis atoms
Co Outputs
Co   Euler angles are read in or written out
Co   neul  is read in or written out (either 1, nl or nl**2)
Co   xsi is read in if xsi=# is available on the first line
Co   xsi is written out if xsi is nonzero
Cr Remarks
Cr   Aborts on read when nbas does not match file
Cu Updates
Cu   24 May 08 Only read/write master (MPI)
Cu             Note: it is the caller's responsibility to broadcast
Cu             This is because neul and the dimensions of eula
Cu             are only known after reading.  Typical Broadcast:
Cu             call mpibc1(neul,1,2,mlog,'susite','neul')
Cu             call mpibc1(xsi,5,4,mlog,'susite','xsi')
Cu             call mpibc1(eula,nbas*neul*3,4,mlog,'susite','eula')
Cu   19 Nov 97 use rdfiln for reading first line
C ----------------------------------------------------------------
C     implicit none
      integer nbas,neul,nl,ifi
      double precision eula(nbas,neul,3),tmp(125),xsi(3)
      integer i,j,k,ipr,lgunit,nc,rdm,nxsi,a2vec,ix(10)
      logical parstr,a2bin
      character*80 ss,ss2*4
      integer master,mpipid
C ... for rdfiln
      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character recrd*(recl),ctbl(mxchr,2)*(ctlen),a*(recl),
     .vnam(mxlev)*16,rdarg*6
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .nlist(0:mxlev)
      equivalence (a,ss)
      data rdarg /'#{}% c'/


C --- I/O through master node only (MPI) ---
      master = 0
      if (mpipid(1) .eq. master) then

C --- Write ---
        if (ifi .lt. 0) then
C       rewind (-ifi)
          if (neul .ne. 1 .and. neul .ne. nl .and. neul .ne. nl*nl)
     .    call rxi('ioeula: bad dim for eula: neula=',neul)
          call awrit3('%% rows %i cols %i neula %i',ss,80,0,nbas*neul,3,
     .    neul)
          nxsi = 0
          do  11  i = 1, 3
            if (xsi(i) .eq. 0) goto 12
            nxsi = nxsi+1
   11     continue
   12     continue
          if (nxsi .gt. 0)
     .    call awrit3('%a nxsi=%i %n:1;6,6d',ss,80,0,nxsi,nxsi,xsi)
          call awrit0('%a',ss,80,ifi)
          do  10  i = 1, nbas
            if (neul .gt. 1) call awrit1('# ib %,4i',' ',80,-ifi,i)
            write(-ifi,'(3f16.12)') ((eula(i,j,k),k=1,3),j=1,neul)
   10     continue
C --- Read ---
        else
          call getpr(ipr)
          rewind ifi
          nr = 0
          call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .    ilist,nlist,vnam,ctbl,mxchr,ss,recrd,recl,nr)
C   ... No error if file is empty, but don't bother reading
          if (nr .eq. 0) then
            if (ipr .gt. 30) print '('' IOEULA:   empty file'')'
            return
          endif
          if (ss(1:1) .eq. '%') then
C   ... Read xsi if it is there
            i = 0
            if (parstr(ss,'nxsi=',len(ss)-5,5,' ',i,j)) then
              j = j-1
              i = j
              if (a2bin(ss,nxsi,2,0,' ',j,len(ss)-5)) then
                j = a2vec(ss,len(ss)-5,i,4,' #',2,2,nxsi+1,ix,tmp)
                if (j .ne. nxsi+1) call rx('ioeula:  failed to read xsi')
                call dcopy(nxsi,tmp(2),1,xsi,1)
              endif
            endif
            neul = 1
            i = 0
            if (parstr(ss,'neula ',len(ss)-6,6,' ',i,j)) then
              if (.not. a2bin(ss,neul,2,0,' ',j,len(ss)))
     .        call rx('ioeula: failed to parse for neula')
            endif
          else
            neul = 1
            rewind ifi
          endif

C   ... Expect nbas*neul*3 rows in file ... skip reading if not
          if (neul .ne. 1 .and. neul .ne. nl .and. neul .ne. nl*nl)
     .    call rxi('ioeula: illegal dimension: neula =',neul)
          nc = 3
          ss2 = 'in'
          if (neul .eq. nl) ss2 = ' '
          if (neul .eq. nl*nl) ss2 = '%bm-'
          if (ipr .gt. 30) call awrit1(' IOEULA:  reading l-'//ss2//
     .    '%adependent Euler angles (file %i)',' ',80,lgunit(1),ifi)
C        if (nr .ne. neul*nbas .or. nc .ne. 3) then
C          if (ipr .ge. 20) call awrit3(' ioeula: input skipped --'//
C     .      ' expected %i rows (3 cols) but found %i (%i cols)',
C     .      ' ',80,lgunit(1),neul*nbas,nr,nc)

C   ... Read the data, 1 atom at a time ...
          call dpzero(eula,nbas*neul*3)
          do  20  i = 1, nbas
            call rxx(rdm(ifi,0,neul*nc,' ',tmp,neul,nc).ne.1,
     .      'ioeula:  bad file')
            call dmcpy(tmp,neul,1,eula(i,1,1),nbas*neul,nbas,neul,nc)
   20     continue
          if (ipr .ge. 50) call pvioeu(1,0d0,eula,nbas,neul)
        endif
      endif
      end subroutine ioeula

      subroutine pvioeu(mode,ssite,eula,nbas,neul)

      use m_struc_def  !Cgetarg

C- Writes the Euler angles to a file or to stdout
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 write in format suitable for rdm
Ci         :1 write in printout format
Ci         :10s digit
Ci         :0 printout is for Euler angles
Ci         :1 printout is for Bfield
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    clabel
Ci   eula  :Euler angles for noncollinear spins
Ci   nbas  :size of basis
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Co Outputs
Co   Writes to stdout the Euler angles
Cr Remarks
Cr
Cu Updates
Cu   14 Feb 03 output can be for euler angles or for B-field
Cu   22 May 02 Writes label information if passed through in ssite
Cu             Altered argument list
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,neul
      real(8):: eula(nbas,neul,3)
      type(s_site)::ssite(*)

C ... Local parameters
      integer i,j,k,lgunit,stdo,mode0,mode1,lbl
      double precision rotm(3,3)
      character*72 outs, clabl*8
      character*255:: formatc(344:346),lbll

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      stdo = lgunit(1)
      clabl = ' '
      if (mode0 .ne. 0) then
c        assign 344 to lbl
c        if (mode1 .ne. 0) assign 345 to lbl
c        if (mode1 .eq. 0 .and. mode0 .ne. 0) assign 346 to lbl
        lbl=344
        if (mode1 .ne. 0) lbl=345
        if (mode1 .eq. 0 .and. mode0 .ne. 0) lbl=346
        formatc(344)="('   ib     alpha        beta       gamma':6x,a)"
        formatc(345)="('   ib      bx          by          bz  ':6x,a)"
        formatc(346)="('   ib   alpha      beta     gamma':5x,a,"//
     .  "11x,'z axis rotated to:')"
        lbll=formatc(lbl)

        if ( neul .gt. 1)  write(outs,trim(lbll))
c        if (ssite(1)%size .eq. 0 .or. neul .gt. 1)  write(outs,trim(lbll))
c        if (ssite(1)%size .ne. 0 .and. neul .eq. 1) write(outs,trim(lbll)) 'class'
      endif
c  344 format('   ib     alpha        beta       gamma':6x,a)
c  345 format('   ib      bx          by          bz  ':6x,a)
c  346 format('   ib   alpha      beta     gamma':5x,a,
c     .  11x,'z axis rotated to:')
      if (neul .gt. 4) then
        call awrit0('%a     lm',outs,len(outs),0)
      elseif (neul .gt. 1) then
        call awrit0('%a      l',outs,len(outs),0)
      endif
      call awrit0(outs,' ',-len(outs),stdo)
      do  21  i = 1, nbas
c        if (ssite(1)%size .ne. 0) then
c          do i_spacks=i,i
c            call spacks_copy('u',ssite(i_spacks)%clabel,i,i,clabl,i_spacks)
c          enddo
c        endif
        if (neul .gt. 1 .and. mode0 .ne. 0) then
          if (clabl .eq. ' ') call awrit1('# site %,4i',' ',80,stdo,i)
          if (clabl .ne. ' ') call awrit1('# site %,4i  class '//clabl,
     .    ' ',80,stdo,i)
          write(stdo,'(i5,3f12.6,i4)')
     .    (i,(eula(i,j,k),k=1,3),j, j=1,neul)
        elseif (neul.eq.1 .and. mode0.ne.0 .and. clabl.eq.' ') then
          write(stdo,'(i5,3f12.6)')
     .    (i,(eula(i,j,k),k=1,3), j=1,neul)
        elseif (neul .eq. 1 .and. mode0 .ne. 0) then
          call eua2rm(eula(i,1,1),eula(i,1,2),eula(i,1,3),rotm)
          write(stdo,'(i5,3f10.6,3x,a,1x,3f10.6)')
     .    i,(eula(i,1,k),k=1,3), clabl, (rotm(3,k),k=1,3)
        else
          write(stdo,'(3f16.12)') ((eula(i,j,k),k=1,3),j=1,neul)
        endif
   21 continue

      end subroutine pvioeu


