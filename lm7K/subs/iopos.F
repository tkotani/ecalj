Cdelw1 var ['owk', '3', 'owk', '3*nbas', 'real(8)', 'rv_w_', 'iopos']
Cdelw1 var ['opos', '0', 'opos', '3*nbas', 'real(8)', 'rv_w_', 'ioposs']
Cgetarg...info...           structure ['ioposs', 'site', 'ssite']
Cgetarg...info...           use_to_add ['ioposs', 'm_struc_def']
      subroutine iopos(lio,mode,filnam,nbas,bas)

C- File i/o of site positions
C ----------------------------------------------------------------
Ci Inputs
Ci  lio:    true for write, false for read
Ci  mode:   0 read all sites as a matrix.
Ci
Ci          1 read each line in the form
Ci            ib x y z
Ci          Here ib must be between 1 and nbas.
Ci          Not all sites need be specified; sites not
Ci          specified are not changed.
Ci
Ci         -1 mode unknown; must be specified by input file;
Ci          see Remarks.
Ci filnam:  file name containing data
Ci  nbas:   size of basis
Cio bas:    is read in (lio=F), or written out (lio=T)
Cr Remarks
Cr   File mode can be specified by input file, which it does
Cr   by first line beginning with % and containing mode=#, e.g.
Cr     % mode=1
Cr
Cr   See also iosits.
Cm MPI
Cm   Master does I/O. If read then bas is broadcast.
C ----------------------------------------------------------------
C     implicit none
#if MPI | MPIK
      include "mpif.h"
      integer procid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
      logical mlog, cmdopt
      character*120 strn
#endif
      logical lio
      integer mode,nbas
      character*(*) filnam
      double precision bas(3,nbas)
C Local variables
      logical parstr,a2bin
Cdelw1       integer j1,j2,ifi,fopna,owk,j,rdm,ipr,lgunit,i,lmode,a2vec,ix(4)
       integer:: j1 , j2 , ifi , fopna , j , rdm , ipr , lgunit , i 
     . , lmode , a2vec , ix(4) 
       real(8) ,allocatable :: rv_w_owk(:)

      double precision xv(4)
      character*100 fnam
C ... for rdfiln
      integer recl,nr,mxchr,mxlev,lstsiz,ctlen
      parameter (mxchr=20,mxlev=4,lstsiz=200,recl=500,ctlen=120)
      character recrd*(recl),ctbl(mxchr,2)*(ctlen),a*(recl),
     .  vnam(mxlev)*16,rdarg*6
      logical loop0(0:mxlev)
      integer nlin(0:mxlev),list(lstsiz,mxlev),ilist(mxlev),
     .  nlist(0:mxlev)
C heap:
      integer w(1)
      common /w/ w

      data rdarg /'#{}% c'/

#if MPI | MPIK
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
      if (procid .eq. master) then
#endif
C ... Open file
      j1 = 1
      fnam = filnam
      call nword(fnam,1,j1,j2)
      i = 1
      if (lio) i = 2
      ifi = fopna(fnam(j1:j2),-1,i)
      rewind ifi
      call getpr(ipr)

C --- File write ---
      if (lio) then
C   ... Write, mode=1
        if (mode .eq. 1) then
          write(ifi,'(''% mode=1'')')
          do  i = 1, nbas
            write(ifi,'(i4,3f13.7)') i,bas(1,i),bas(2,i),bas(3,i)
          enddo
          call fclose(ifi)

C   ... Write, mode=0
        elseif (mode .eq. 0) then
Cdelw1           call defrr(owk, 3*nbas)
           allocate(rv_w_owk(3*nbas))
           if (3*nbas<0) rv_w_owk(:)=0.0d0

Cdelw1           call dmcpy(bas,3,1,w(owk),1,nbas,3,nbas)
           call dmcpy ( bas , 3 , 1 , rv_w_owk , 1 , nbas , 3 , nbas ) 

Cdelw1           call ywrm(0,' ',1,ifi,'(3f14.8)',w(owk),1,nbas,nbas,3)
           call ywrm ( 0 , ' ' , 1 , ifi , '(3f14.8)' , rv_w_owk , 1 , nbas 
     .     , nbas , 3 ) 

          call fclose(ifi)
        else
          call rxi('IOPOS: not implemented, mode=',mode)
        endif

        call info2(10,1,1,' iopos : write %i sites into file '//
     .      fnam(j1:j2)//', mode %i',nbas,mode)

C --- File read ---
      else
C   ... Look for mode on first line
        nr = 0
        call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .    ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
        backspace ifi
        lmode = 0
        if (recrd(1:1) .eq. '%') then
          i = 0
          j = 0
          if (parstr(recrd,'mode=',len(recrd),5,'=',i,j)) then
            i = j
            if (.not. a2bin(recrd,lmode,2,0,' ',i,-1))
     .        call rxs('IOPOS failed to parse ...',recrd(j-5:))
          endif
        endif
        if (mode .ge. 0 .and. lmode .ne. mode)
     .    call rxi('IOPOS: mode mismatch: need mode',mode)
C   ... Read, mode=0
        if (lmode .eq. 0) then
Cdelw1           call defrr(owk, 3*nbas*2)
           allocate(rv_w_owk(3*nbas*2))
           if (3*nbas*2<0) rv_w_owk(:)=0.0d0

Cdelw1           j = rdm(ifi,0,3*nbas,' ',w(owk),nbas,3)
           j = rdm ( ifi , 0 , 3 * nbas , ' ' , rv_w_owk , nbas , 3 ) 

          if (j .lt. 0) call rxi(
     .      'IOPOS: file incompatible with nbas=',nbas)
Cdelw1           call dmcpy(w(owk),nbas,1,bas,1,3,nbas,3)
           call dmcpy ( rv_w_owk , nbas , 1 , bas , 1 , 3 , nbas , 3 ) 

Cdelw1 rlse name= owk old_list= owk 
Cdelw1 rlse name= owk new_list= (None)
Cdelw1           call rlse(owk)
           if (allocated(rv_w_owk)) deallocate(rv_w_owk)

          nr = nbas
C   ... Read, mode=1
        elseif (lmode .eq. 1) then
          nr = 0
   41     call rdfiln(ifi,rdarg,mxlev,loop0,nlin,list,lstsiz,
     .      ilist,nlist,vnam,ctbl,mxchr,a,recrd,recl,nr)
          if (nr .ge. 0) then
            if (recrd(1:1) .eq. '%') goto 41
            i = 0
            j = a2vec(recrd,len(recrd),i,4,', ',2,3,4,ix,xv)
            if (j .ne. 4) call rxs(
     .        'IOPOS: failed to parse line:',recrd)
            i = nint(xv(1))
            call dcopy(3,xv(2),1,bas(1,i),1)
            goto 41
          endif
          nr = -nr-1
        else
          call rxi('IOPOS: unkown mode',lmode)
        endif
C     ... Show input
        if (ipr .ge. 0) then
          call awrit3(' iopos : read %i sites from file '//fnam(j1:j2)
     .      //'%a, mode=%i%?#n>=40#:'//
     .      '%N  ib       x           y           z',
     .      ' ',100,lgunit(1),nr,lmode,ipr)
          if (ipr .ge. 40) then
            do  10  i = 1, nbas
   10       print 310, i, bas(1,i), bas(2,i), bas(3,i)
  310       format(i4,3f12.6)
          endif
        endif
      endif
      call fclose(ifi)
#if MPI | MPIK
      endif
      if (.not. lio) then
        call MPI_BCAST(bas,3*nbas,MPI_DOUBLE_PRECISION,master,
     .                 MPI_COMM_WORLD,ierr)
        if (mlog) then
          call gettime(datim)
          call awrit3(' iopos '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' bcast bas, %i d.p. numbers',' ',
     .      256,lgunit(3),procid,numprocs,3*nbas)
        endif
      endif
#endif

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine iopos 

      subroutine ioposs(lio,mode,filnam,nbas,ssite)
       
       use m_struc_def  !Cgetarg

C- File i/o of site positions from site structure
C     implicit none
      logical lio
      integer mode,nbas
Cgetarg       double precision ssite(1)
       type(s_site)::ssite(*)

      character*(*) filnam
Cdelw1       integer opos
       real(8) ,allocatable :: rv_w_opos(:)

C heap:
      integer w(1)
      common /w/ w

C ... Unpack site positions
Cdelw1       call defrr(opos,3*nbas)
       allocate(rv_w_opos(3*nbas))
       if (3*nbas<0) rv_w_opos(:)=0.0d0

Cgetarg       call spackv(10,'site pos',ssite,1,nbas,w(opos))
       i_copy_size=size(ssite(1)%pos)
       do i_spackv=1,nbas 
Cdelw1        call spackv_array_copy_r8_r8('u',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
        call spackv_array_copy_r8_r8 ( 'u' , ssite ( i_spackv ) %pos 
     .  , i_copy_size , i_spackv + 1 - 1 , rv_w_opos ) 

       enddo

C ... File I/O
Cdelw1       call iopos(lio,mode,filnam,nbas,w(opos))
       call iopos ( lio , mode , filnam , nbas , rv_w_opos ) 

C ... Repack site positions
Cgetarg       call spackv(11,'site pos',ssite,1,nbas,w(opos))
       i_copy_size=size(ssite(1)%pos)
       do i_spackv=1,nbas 
Cdelw1        call spackv_array_copy_r8_r8('p',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
        call spackv_array_copy_r8_r8 ( 'p' , ssite ( i_spackv ) %pos 
     .  , i_copy_size , i_spackv + 1 - 1 , rv_w_opos ) 

       enddo

Cdelw1 rlse name= opos old_list= opos 
Cdelw1 rlse name= opos new_list= (None)
Cdelw1       call rlse(opos)
       if (allocated(rv_w_opos)) deallocate(rv_w_opos)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine ioposs 


