#define F90 1
      integer function iorsj(mode,ifi,nbas,alat,plat,aamom,iax,rtab,
     .  ntab,rsJ,rmax)
C- File read of coefficients to micromagnetics (pairwise) hamiltonian
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :mode controls what iorsj does on file READ or WRITE
Ci         :... for file READ
Ci         :1s digit mode (controls program flow):
Ci         :0 Setup mode. iorsj does:
Ci         :  1. read file nbas and alat
Ci         :  2. check file nbas against passed value
Ci         :     or copy file value to nbas if input nbas = 0
Ci         :  3. read file plat
Ci         :  4. check file plat against passed array
Ci         :     or copy file value to plat if input plat = 0
Ci         :  5. read file local moments
Ci         :  6. Find rmax = length of largest connecting vector
Ci
Ci         :2 read mode. iorsj does:
Ci         :  1. read file local moments
Ci         :  2. Match file entries against iax table,
Ci         :     and read matching entries into rsJ
Ci         :  By default, a match is with pair it found when:
Ci         :  * file ib matches iax(1,it)
Ci         :  * connecting vector rtab(it) matches file value
Ci         :  If 10's digit mode is nonzero, also require for match:
Ci         :  * file jb matches iax(2,it)
Ci
Ci         :1 match mode. Its usual purpose is to prepare the iax table
Ci         :  for purging of pairs missing from the (pair) hamiltonian
Ci         :  read from the file.  iorsj does:
Ci         :  1. Match file entries against iax table, and mark missing
Ci         :     entries by setting iax(1) to zero
Ci         :  The match criterion applies to this mode
Ci         :  in the same way as it does to the read mode.
Ci
Ci         :3 match mode same as mode=1, but uses rsJ from file
Ci
Ci         :NB: match and read modes may be done together.
Ci         :However, generally not advisable because usual purpose of
Ci         :match mode is to setup shortening of iax table, which
Ci         :should be done after match and before read.
Ci
Ci         :10s digit mode (??):
Ci
Ci         :... for file WRITE
Ci
Ci         :0  write header for rsj file
Ci         :1  write exchange parameters in standard form
Ci         :6  write exchange parameters in form suitable for cvm input
Ci
Ci   ifi   :file logical unit, >0 for read and <0 for write
Ci   nbas  :size of basis
Ci   plat  :primitive lattice vectors, in units of alat
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   rtab  :site positions corresponding to entries in a neighbor table
Ci         :relative to some center.  Generate from iax using routine mkrtab
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Co Outputs depend on whether setup, match, or read modes (0,1,2)
Co  iorsj  : (setup read mode) -1 if plat mismatch
Co         : (setup read mode) Otherwise, number of entries in file
Co         : (match read mode) number of file entries matching iax table
Co   rsJ   : (read mode)  coefficients to pairwise hamiltonian
Co   rmax  : (setup read mode) largest connecting vector in table
Co   iax   : (match read mode) iax(1) set to 0 for missing file entries
Cr Remarks
Cr
Cu Updates
Cu   02 Oct 03 Implement write of exch. parameters for CVM input
Cu   24 Nov 02 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ifi,nbas,ntab(nbas+1),niax
      parameter (niax=10)
      integer iax(niax,*)
      double precision alat,rmax
      double precision aamom(nbas),plat(3,3),rsJ(*),rtab(3,*)
C ... Local parameters
      logical lmatch,lread
      integer nread,i,ib,jb,it,nmatch,mode0,mode1,jfi,nbasl,nttab,ib0,
     .  jt,kt
      double precision platl(3,3),ddot,tol,xi(3),cJ,alatl
      parameter (tol=1d-5)
#if F90
      real(8),allocatable :: prham(:)
#else
      double precision prham(1)
#endif

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)

C --- Setup mode ---
C     1. check file plat against passed array
C     2. count number of pairs in table
C     3. Find rmax = length of largest connecting vector
      if (mode0 .eq. 0) then

C       ifi<0: write header data
        if (ifi .lt. 0) then
          jfi = -ifi
          rewind jfi
          write(jfi,'('' nbas, alat:'')')
          write(jfi,'(i4,f12.6)') nbas, alat
          write(jfi,'('' Lattice vectors:''/9f13.7)') plat
          write(jfi,'('' Local moments:'')')
          write(jfi,'(9f12.6)') aamom
          write(jfi,'(''  ib  jb'',16x,''vector'',23x,''J'')')
          iorsj = 0
          return
        endif

        rewind ifi
        iorsj = -1
        rmax = 0
        nread = -1

C       Read file nbas, alat
        read(ifi,*,err=10,end=10)
        read(ifi,*,err=10,end=10) nbasl, alatl
        if (nbas .eq. 0) nbas = nbasl
        call isanrg(nbasl,nbas,nbas,'iorsj:','file''s nbas',.true.)
        call fsanrg(alatl,alat,alat,tol,'iorsj:','file''s alat',.true.)



C       Read file plat and compare or copy to passed array
        read(ifi,*,err=10,end=10)
        read(ifi,*,err=10,end=10) platl
        if (ddot(9,plat,1,plat,1) .eq. 0) then
          call dcopy(9,platl,1,plat,1)
        else
          do  i = 1, 9
            if (dabs(plat(i,1)-platl(i,1)) .gt. tol) return
          enddo
        endif

C       Read file local moments
        read(ifi,*,err=10,end=10)
        read(ifi,*,err=10,end=10) aamom

C       Position file pointer past end of header
        read(ifi,*,err=10,end=10)

        nread = 0
        do  i = 1, 9999999
          read(ifi,*,err=10,end=10) ib,jb,xi,cJ
          rmax = max(rmax,ddot(3,xi,1,xi,1))
          nread = nread+1
        enddo
        call rx('too many lines in micromagnetics file')
   10   continue
        rmax = sqrt(rmax)
        iorsj = nread
        return
      endif

C --- Match and read modes ---
C     1. Match file entries against iax table, and mark matching entries
C        Mark: iax(1) = 0 for entries not matched in file
C        Caller may subsequently purge iax table using
      if (mode0 .ge. 1 .and. mode0 .le. 3) then

C       ifi<0: write entire file
        if (ifi .lt. 0) then
          jfi = -ifi
          rewind jfi
          write(jfi,'('' nbas, alat:'')')
          write(jfi,'(i4,f12.6)') nbas, alat
          write(jfi,'('' Lattice vectors:''/9f13.7)') plat
          write(jfi,'('' Local moments:'')')
          write(jfi,'(9f12.6)') aamom
          write(jfi,'(''    ib    jb'',16x,''vector'',18x,''J'')')

          nttab = ntab(nbas+1)
          do  it = 1, nttab
            ib = iax(1,it)
            jb = iax(2,it)
            call awrit6('%,6i%,6i%;12,7D%;12,7D%;12,7D %;6g',' ',80,jfi,
     .        ib,jb,rtab(1,it),rtab(2,it),rtab(3,it),rsJ(it))
          enddo

          iorsj = nttab
          return
        endif


        lmatch = (mod(mode0,2) .eq. 1)
        lread  = (mode0 .ge. 2)
        iorsj = -1

C   ... Read file local moments, position pointer past end of header
C       Position file past end of header
        rewind ifi
        read(ifi,*)
        read(ifi,*) nbasl, xi(1)
        read(ifi,*)
        read(ifi,*) platl
        read(ifi,*)
        read(ifi,*) aamom
        read(ifi,*)

C       Read each entry, mark in iax by setting iax(1)<0
        nmatch = 0
        nread = 0
        do  i = 1, 9999999
          read(ifi,*,err=110,end=110) ib,jb,xi,cJ
          nread = nread+1
          if (ib .le. nbas) then
C       ... Loop until a match is found for this entry
            do  it = ntab(ib)+1, ntab(ib+1)
              if (iax(1,it) .lt. 0) goto 101
              if (mode1 .ne. 0) then
                if (iax(2,it) .ne. jb) goto 101
              endif
              if (dabs(xi(1)-rtab(1,it)) .gt. tol) goto 101
              if (dabs(xi(2)-rtab(2,it)) .gt. tol) goto 101
              if (dabs(xi(3)-rtab(3,it)) .gt. tol) goto 101
C             A match is found
              nmatch = nmatch+1
              if (lmatch) iax(1,it) = -ib
              if (lread) rsJ(it) = cJ
              goto 102
  101         continue
            enddo
C           No match was found for this entry
C           call info2(30,0,0,' iorsj: no match found for site %i, '//
C    .        'connecting vector %3:1,6;12D',ib,xi)
  102       continue
          endif
        enddo
        call rx('too many lines in micromagnetics file')

C   ... End of file read
  110   continue
        call info(30,0,0,' iorsj: matched %i pair entries out of %i'//
     .    ' file values',nmatch,nread)

C   ... Match mode: put iax(1)=0 for missing entries
        if (lmatch) then
          do  ib = 1, nbas
            do  it = ntab(ib)+1, ntab(ib+1)
              if (iax(1,it) .lt. 0) then
                iax(1,it) = ib
              else
                iax(1,it) = 0
              endif
            enddo
          enddo
        endif

C       Purge iax table of unused entries
C      call symiax(1,plat,nbas,xi,xi,xi,0,ntab,iax,nttab,mxcsz)
C      if (nttab .ne. nmatch) call rx('iorsj: bug in symiax')

        iorsj = nmatch
        return
      endif

C --- German's CVM format ---
      if (mode0 .ge. 6 .and. mode0 .le. 7) then


C       ifi<0: write entire file
        if (ifi .lt. 0) then
          jfi = -ifi
          rewind jfi

          if (mode0 .eq. 7) then
#if F90
            allocate(prham(nbas))
#else
            call rx('iorsj:  mode 7 requires f90 compiler')
#endif
          endif

C     ... this header stuff needs to change
          if (mode0 .eq. 6) then
            write(jfi,'('' nbas, alat:'')')
            write(jfi,'(i4,f12.6)') nbas, alat
            write(jfi,'('' Lattice vectors:''/9f13.7)') plat
            write(jfi,'('' Local moments:'')')
            write(jfi,'(9f12.6)') aamom
            write(jfi,'(6x,''  ib    jb'',6x,''J'')')
          endif

          ib0 = 0
          nttab = ntab(nbas+1)
C     ... For each pair, do
          do  it = 1, nttab

            ib = iax(1,it)

C       ... Case new cluster
            if (ib .ne. ib0) then
              ib0 = ib

C             Find last pair in ib cluster
              do  jt = it, nttab
                if (ib .ne. iax(1,jt)) goto 611
                kt = jt
              enddo
  611         continue
C             print *, 'ib,it,kt=',ib,it,kt


C         ... For each atom jb, count number of pairs and cum J
              if (mode0 .eq. 7) call dpzero(prham,nbas)

              do  jb = 1, nbas
C               ncJ = ncJ+1
C               cJ = 0
                do  jt = it, kt
                  if (jb .eq. iax(2,jt)) then
                    if (mode0 .eq. 7) then
                      prham(jb) = prham(jb) + rsJ(jt)
                    else
                      call awrit4('%,4i%,6i%,6i  %;6g',' ',80,jfi,
     .                  1,ib,jb,rsJ(jt))
                    endif
                  endif
                enddo
              enddo

              if (mode0 .eq. 7) then
                write(jfi,'(1p,12e14.5)') (prham(jt),jt=1,nbas)
              endif


C       ... This cluster is done ... loop until new one starts
            else
            endif

          enddo

          if (mode0 .eq. 7) then
#if F90
            deallocate(prham)
#endif
          endif

          iorsj = nttab
          return
        endif

      endif
      end

