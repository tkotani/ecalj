      logical function ioorbp(mode,nkapi,is1,is2,sspec,nkmax,ifi)

      use m_struc_def  !Cgetarg

C- File i/o of parameters defining orbitals in basis envelope functions
C ----------------------------------------------------------------------
Ci inputs
Ci   mode  :1 parameters consist of rsmh,eh in sspec->orbp
Ci            and if nkapi=2 rsmh2,eh2
Ci            and if PZ is present, values are put into orbp(*,3)
Ci          Add 10 if to determine idxdn locally
Ci          Add 100 if basis is allowed to be only partially complete:
Ci                  basis for missing species are unchanged
Ci   nkapi :number of envelope fns of one l-quantum number in the basis
Ci   is1..2:range of species to read/write basis parameters
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Co Inputs/Outputs
Co  sspec   :struct containing species-specific information
Co          :Case mode = 1:
Co          :sspec->rtab table of sm-H smoothing radii is I/O
Co          :sspec->etab table of sm-H energies is I/O
Co  nkmax   :(used in file read only): largest nkapi encountered in is1..is2
Co          :10's digit set to 1 if any PZ is read
Co  ioorbp  :false: failed to read basis from file
Co          :true: successfully read basis from file
Cr Remarks
Cb Bugs
Cb   This routine doesn't check for local orbitals
Cu Updates
Cu   02 Sep 03 ioorbp can read local orbitals.
Cu             Miscellaneous other additions to read mode
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   24 Aug 01 Extended to local orbitals.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nkapi,is1,is2,ifi,nkmax,i_copy_size
      type(s_spec)::sspec(*)

C ... Local parameters
      logical scat,rdstrn
      character spid*8,sout*256
      integer jfi,n0,parg,norp,ntorb,nkap0,nkape
      parameter (n0=10,nkap0=3,nkape=2)
      integer is,ks,lp1,jp,ikap,lmxb,lmxi,j1,j2,ls,i,j,lpz,
     .it(n0),idxdn(n0,nkap0)
      double precision rtab(n0,nkap0),etab(n0,nkap0),orbp(n0,2,nkap0),
     .dasum
      logical:: l_dummy_isanrg, isanrg

      if (mod(mode,10) .ne. 1 .or. mod(mode/10,10) .gt. 1 .or.
     .mod(mode/100,10) .gt. 1)
     .call rxi('unexpected value for mode, mode=',mode)

      lpz = 0
      ioorbp = .true.
      ls = len(sout)
      if (is2 .lt. is1) return
      if (ifi .gt. 0) then
        nkmax = 1
        ioorbp = .false.
        if (mod(mode,10) .eq. 1) then
          if (.not. scat(ifi,'BASIS:',':',.true.)) return
          if (mod(mode/100,10) .eq. 1) ioorbp = .true.
          ks = is1-1
          do  10  is = is1, is2
            ks = ks+1
   18       if (.not. rdstrn(ifi,sout,ls,.false.)) return
            if (sout(1:1) .eq. '#') goto 18
            call word(sout,1,j1,j2)
   16       continue
c            do i_spacks=ks,ks
c              call spacks_copy('u',sspec(i_spacks)%name,ks,ks,spid,i_spacks)
c            enddo
            spid=sspec(ks)%name
            if (spid .ne. sout(j1:j2)) then
              ks = ks+1
              if (ks .gt. is2) then
                call rxs('ioorbp: species unknown or out of order: ',
     .          sout(j1:j2))
              endif
              goto 16
            endif
            call dpzero(orbp,n0*2*nkap0)

            sspec(ks)%norp=2
            sspec(ks)%ntorb=n0*nkape

            call iinit(idxdn,n0*nkap0)

            lmxb=sspec(ks)%lmxb
            i_copy_size=size(sspec(ks)%idxdn)
            call icopy(i_copy_size,sspec(ks)%idxdn,1,idxdn,1)

C         Case determine idxdn locally
            if (mod(mode/10,10) .ne. 0) call iinit(idxdn,n0*nkape)
            if (lmxb+1 .gt. n0) call rx('increase size(idxdn) in sspec')

C     --- Read first set of parameters RSMH,EH ---
            jp = j2
            i = parg(' RSMH=',4,sout,jp,ls,', ',3,lmxb+1,it,orbp(1,1,1))
C         no smoothing radii read ... skip this species
            if (i .eq. 0 .or. i .eq. -1) goto 10
C         after next block, lmxi = true lmxb for this rsmh,eh
            lmxi = lmxb
C         Fewer rsm read than lmxb : l's beyond last rsm -> high
            if (i .lt. 0) then
              do  j  = -i, lmxb+1
C             print *, j, orbp(j,1,1)
                orbp(j,1,1) = 0
                idxdn(j,1) = 4
              enddo
              lmxi=-i-2
            endif
C         For l's > lmxb, set idxdn as follows:
            do  lp1 = 1, lmxi+1
C           if rsm > 0 and idxdn "undetermined", idxdn => low
              if (orbp(lp1,1,1) .gt. 0 .and. idxdn(lp1,1) .eq. 0)
     .        idxdn(lp1,1) = 1
C           if rsm < 0, idxdn -> high
              if (orbp(lp1,1,1) .lt. 0) idxdn(lp1,1) = 4
C           if rsm = 0, idxdn -> intermediate.  ? Why is this here?
C           if (orbp(lp1,1,1) .eq. 0) idxdn(lp1,1) = max(idxdn(lp1,1),2)
C           if rsm = 0, idxdn -> high
              if (orbp(lp1,1,2) .eq. 0) idxdn(lp1,2) = 4
            enddo
            jp = j2
            i = parg(' EH=',4,sout,jp,ls,', ',3,lmxi+1,it,orbp(1,2,1))
            if (i .eq. 0 .or. i .eq. -1) goto 10
C         if (i .lt. lmxi+1) goto 10
C         Fewer eh read than lmxb : l's beyond last eh -> high
            if (i .lt. 0) then
              do  j  = -i, lmxi+1
c             print *, j, orbp(j,1,1)
                orbp(j,1,1) = 0
                idxdn(j,1) = 4
              enddo
              lmxi = -i-2
            endif

C     --- Read second set of parameters RSMH2,EH2 ---
            if (nkapi .eq. 2) then
              jp = j2
              i=parg(' RSMH2=',4,sout,jp,ls,', ',3,lmxb+1,it,orbp(1,1,2))
              if (i .ne. 0) then
C           lmxi = true lmxb for this rsmh,eh
                lmxi = lmxb
                if (i .lt. 0) lmxi=-i-2
                call ivset(idxdn(1,2),lmxi+2,n0,4)
                do  lp1 = 1, lmxi+1
C             if rsm > 0 and idxdn "undetermined", idxdn => low
                  if (orbp(lp1,1,2) .gt. 0 .and. idxdn(lp1,2) .eq. 0)
     .            idxdn(lp1,2) = 1
C             if rsm < 0, idxdn -> high
                  if (orbp(lp1,1,2) .lt. 0) idxdn(lp1,2) = 4
C             if rsm = 0, idxdn -> high
C             if (orbp(lp1,1,2) .eq. 0) idxdn(lp1,2) = max(idxdn(lp1,2),2)
                  if (orbp(lp1,1,2) .eq. 0) idxdn(lp1,2) = 4
                enddo
                jp = j2
                i =parg(' EH2=',4,sout,jp,ls,', ',3,lmxi+1,it,orbp(1,2,2))
                if (i .lt. lmxi+1) goto 10
                if (dasum(n0,orbp(1,1,2),1) .ne. 0) nkmax = 2
              else
                call ivset(idxdn(1,2),1,n0,4)
              endif
            endif
C         Reduce lmxb
            do  lp1 = lmxb+1, 1, -1
              if (idxdn(lp1,1).ne.4 .or. idxdn(lp1,nkapi).ne.4) goto 12
              lmxb = lmxb-1
            enddo
   12       continue
            jp = j2
            i=parg(' PZ=',4,sout,jp,ls,', ',3,lmxb+1,it,orbp(1,1,3))
            if (i .ne. 0) then
C           lmxi = lmxb
C           if (i .lt. 0) lmxi=-i-2
              i_copy_size=size(sspec(ks)%pz)
              call dcopy(i_copy_size,orbp(1,1,3),1,sspec(ks)%pz,1)
              lpz = 1
            endif
            i_copy_size=size(sspec(ks)%orbp)
            call dcopy(i_copy_size,orbp,1,sspec(ks)%orbp,1)
            i_copy_size=size(sspec(ks)%idxdn)
            call icopy(i_copy_size,idxdn,1,sspec(ks)%idxdn,1)
cc            sspec(ks)%orbp = orbp
cc            sspec(ks)%idxdn= idxdn
            sspec(ks)%lmxb = lmxb
   10     continue
        endif
        ioorbp = .true.
        if (lpz .eq. 1) nkmax = nkmax + 10
      else
        jfi = -ifi
        write(jfi,'(''BASIS:'')')
        if (mode .eq. 1) then
          do  110  is = is1, is2
c            do i_spacks=is,is
c              call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
c            enddo
            spid=sspec(is)%name

            sspec(is)%norp=2
            sspec(is)%ntorb=n0*nkape

            norp=sspec(is)%norp
            ntorb=sspec(is)%ntorb
            i_copy_size=size(sspec(is)%orbp)
            call dcopy(i_copy_size,sspec(is)%orbp,1,orbp,1)
            lmxb=sspec(is)%lmxb

C         Sanity checks
            l_dummy_isanrg=isanrg(norp,2,2,'ioorbp:','norp',.true.)
            l_dummy_isanrg=isanrg(ntorb,n0*nkapi,n0*nkape,'ioorbp:','ntorb',.true.)
            etab=0d0 !call dpzero(etab,n0*2)
            rtab=0d0 !call dpzero(rtab,n0*2)
            sout = ' ' // spid
            do  112  ikap = 1, nkapi
              call dcopy(n0,orbp(1,1,ikap),1,rtab(1,ikap),1)
              call dcopy(n0,orbp(1,2,ikap),1,etab(1,ikap),1)
              if (dasum(n0,orbp(1,1,ikap),1) .ne. 0) then
                call awrit6('%a RSMH%?#n==2#2##=%n:1;3d EH%?#n==2#2##=%n:1;3d',
     .          sout,ls,0,
     .          ikap,lmxb+1,rtab(1,ikap),
     .          ikap,lmxb+1,etab(1,ikap))
              endif
  112       continue
            lpz = 0
            do  lp1 = 1, lmxb+1
              if (orbp(lp1,1,3) .gt. 0) lpz = lp1
            enddo
            if (lpz .gt. 0)
     .      call awrit2('%a PZ=%n:1;3d',sout,ls,0,lpz,orbp(1,1,3))
            call awrit0('%a',sout,-ls,-jfi)
  110     continue
        endif
      endif
      end function ioorbp


