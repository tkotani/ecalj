c$$$      logical function ioorbp(mode,nkapi,is1,is2,sspec,nkmax,ifi)
c$$$
c$$$      use m_struc_def  !Cgetarg
c$$$
c$$$C- File i/o of parameters defining orbitals in basis envelope functions
c$$$C ----------------------------------------------------------------------
c$$$Ci inputs
c$$$Ci   mode  :1 parameters consist of rsmh,eh in sspec->orbp
c$$$Ci            and if nkapi=2 rsmh2,eh2
c$$$Ci            and if PZ is present, values are put into orbp(*,3)
c$$$Ci          Add 10 if to determine idxdn locally
c$$$Ci          Add 100 if basis is allowed to be only partially complete:
c$$$Ci                  basis for missing species are unchanged
c$$$Ci   nkapi :number of envelope fns of one l-quantum number in the basis
c$$$Ci   is1..2:range of species to read/write basis parameters
c$$$Ci   ifi   :file logical unit, but >0 for read, <0 for write
c$$$Co Inputs/Outputs
c$$$Co  sspec   :struct containing species-specific information
c$$$Co          :Case mode = 1:
c$$$Co          :sspec->rtab table of sm-H smoothing radii is I/O
c$$$Co          :sspec->etab table of sm-H energies is I/O
c$$$Co  nkmax   :(used in file read only): largest nkapi encountered in is1..is2
c$$$Co          :10's digit set to 1 if any PZ is read
c$$$Co  ioorbp  :false: failed to read basis from file
c$$$Co          :true: successfully read basis from file
c$$$Cr Remarks
c$$$Cb Bugs
c$$$Cb   This routine doesn't check for local orbitals
c$$$Cu Updates
c$$$Cu   02 Sep 03 ioorbp can read local orbitals.
c$$$Cu             Miscellaneous other additions to read mode
c$$$Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
c$$$Cu   24 Aug 01 Extended to local orbitals.
c$$$C ----------------------------------------------------------------------
c$$$      implicit none
c$$$C ... Passed parameters
c$$$      integer mode,nkapi,is1,is2,ifi,nkmax,i_copy_size
c$$$      type(s_spec)::sspec(*)
c$$$
c$$$C ... Local parameters
c$$$      logical scat,rdstrn
c$$$      character spid*8,sout*256
c$$$      integer jfi,n0,parg,norp,ntorb,nkap0,nkape
c$$$      parameter (n0=10,nkap0=3,nkape=2)
c$$$      integer is,ks,lp1,jp,ikap,lmxb,lmxi,j1,j2,ls,i,j,lpz,
c$$$     .it(n0),idxdn(n0,nkap0)
c$$$      double precision rtab(n0,nkap0),etab(n0,nkap0),orbp(n0,2,nkap0),
c$$$     .dasum
c$$$      logical:: l_dummy_isanrg, isanrg
c$$$
c$$$      if (mod(mode,10) .ne. 1 .or. mod(mode/10,10) .gt. 1 .or.
c$$$     .mod(mode/100,10) .gt. 1)
c$$$     .call rxi('unexpected value for mode, mode=',mode)
c$$$
c$$$      lpz = 0
c$$$      ioorbp = .true.
c$$$      ls = len(sout)
c$$$      if (is2 .lt. is1) return
c$$$      if (ifi .gt. 0) then
c$$$        nkmax = 1
c$$$        ioorbp = .false.
c$$$        if (mod(mode,10) .eq. 1) then
c$$$          if (.not. scat(ifi,'BASIS:',':',.true.)) return
c$$$          if (mod(mode/100,10) .eq. 1) ioorbp = .true.
c$$$          ks = is1-1
c$$$          do  10  is = is1, is2
c$$$            ks = ks+1
c$$$   18       if (.not. rdstrn(ifi,sout,ls,.false.)) return
c$$$            if (sout(1:1) .eq. '#') goto 18
c$$$            call word(sout,1,j1,j2)
c$$$   16       continue
c$$$c            do i_spacks=ks,ks
c$$$c              call spacks_copy('u',sspec(i_spacks)%name,ks,ks,spid,i_spacks)
c$$$c            enddo
c$$$            spid=sspec(ks)%name
c$$$            if (spid .ne. sout(j1:j2)) then
c$$$              ks = ks+1
c$$$              if (ks .gt. is2) then
c$$$                call rxs('ioorbp: species unknown or out of order: ',
c$$$     .          sout(j1:j2))
c$$$              endif
c$$$              goto 16
c$$$            endif
c$$$            call dpzero(orbp,n0*2*nkap0)
c$$$
c$$$            sspec(ks)%norp=2
c$$$            sspec(ks)%ntorb=n0*nkape
c$$$
c$$$            call iinit(idxdn,n0*nkap0)
c$$$
c$$$            lmxb=sspec(ks)%lmxb
c$$$            i_copy_size=size(sspec(ks)%idxdn)
c$$$            call icopy(i_copy_size,sspec(ks)%idxdn,1,idxdn,1)
c$$$
c$$$C         Case determine idxdn locally
c$$$            if (mod(mode/10,10) .ne. 0) call iinit(idxdn,n0*nkape)
c$$$            if (lmxb+1 .gt. n0) call rx('increase size(idxdn) in sspec')
c$$$
c$$$C     --- Read first set of parameters RSMH,EH ---
c$$$            jp = j2
c$$$            i = parg(' RSMH=',4,sout,jp,ls,', ',3,lmxb+1,it,orbp(1,1,1))
c$$$C         no smoothing radii read ... skip this species
c$$$            if (i .eq. 0 .or. i .eq. -1) goto 10
c$$$C         after next block, lmxi = true lmxb for this rsmh,eh
c$$$            lmxi = lmxb
c$$$C         Fewer rsm read than lmxb : l's beyond last rsm -> high
c$$$            if (i .lt. 0) then
c$$$              do  j  = -i, lmxb+1
c$$$C             print *, j, orbp(j,1,1)
c$$$                orbp(j,1,1) = 0
c$$$                idxdn(j,1) = 4
c$$$              enddo
c$$$              lmxi=-i-2
c$$$            endif
c$$$C         For l's > lmxb, set idxdn as follows:
c$$$            do  lp1 = 1, lmxi+1
c$$$C           if rsm > 0 and idxdn "undetermined", idxdn => low
c$$$              if (orbp(lp1,1,1) .gt. 0 .and. idxdn(lp1,1) .eq. 0)
c$$$     .        idxdn(lp1,1) = 1
c$$$C           if rsm < 0, idxdn -> high
c$$$              if (orbp(lp1,1,1) .lt. 0) idxdn(lp1,1) = 4
c$$$C           if rsm = 0, idxdn -> intermediate.  ? Why is this here?
c$$$C           if (orbp(lp1,1,1) .eq. 0) idxdn(lp1,1) = max(idxdn(lp1,1),2)
c$$$C           if rsm = 0, idxdn -> high
c$$$              if (orbp(lp1,1,2) .eq. 0) idxdn(lp1,2) = 4
c$$$            enddo
c$$$            jp = j2
c$$$            i = parg(' EH=',4,sout,jp,ls,', ',3,lmxi+1,it,orbp(1,2,1))
c$$$            if (i .eq. 0 .or. i .eq. -1) goto 10
c$$$C         if (i .lt. lmxi+1) goto 10
c$$$C         Fewer eh read than lmxb : l's beyond last eh -> high
c$$$            if (i .lt. 0) then
c$$$              do  j  = -i, lmxi+1
c$$$c             print *, j, orbp(j,1,1)
c$$$                orbp(j,1,1) = 0
c$$$                idxdn(j,1) = 4
c$$$              enddo
c$$$              lmxi = -i-2
c$$$            endif
c$$$
c$$$C     --- Read second set of parameters RSMH2,EH2 ---
c$$$            if (nkapi .eq. 2) then
c$$$              jp = j2
c$$$              i=parg(' RSMH2=',4,sout,jp,ls,', ',3,lmxb+1,it,orbp(1,1,2))
c$$$              if (i .ne. 0) then
c$$$C           lmxi = true lmxb for this rsmh,eh
c$$$                lmxi = lmxb
c$$$                if (i .lt. 0) lmxi=-i-2
c$$$                call ivset(idxdn(1,2),lmxi+2,n0,4)
c$$$                do  lp1 = 1, lmxi+1
c$$$C             if rsm > 0 and idxdn "undetermined", idxdn => low
c$$$                  if (orbp(lp1,1,2) .gt. 0 .and. idxdn(lp1,2) .eq. 0)
c$$$     .            idxdn(lp1,2) = 1
c$$$C             if rsm < 0, idxdn -> high
c$$$                  if (orbp(lp1,1,2) .lt. 0) idxdn(lp1,2) = 4
c$$$C             if rsm = 0, idxdn -> high
c$$$C             if (orbp(lp1,1,2) .eq. 0) idxdn(lp1,2) = max(idxdn(lp1,2),2)
c$$$                  if (orbp(lp1,1,2) .eq. 0) idxdn(lp1,2) = 4
c$$$                enddo
c$$$                jp = j2
c$$$                i =parg(' EH2=',4,sout,jp,ls,', ',3,lmxi+1,it,orbp(1,2,2))
c$$$                if (i .lt. lmxi+1) goto 10
c$$$                if (dasum(n0,orbp(1,1,2),1) .ne. 0) nkmax = 2
c$$$              else
c$$$                call ivset(idxdn(1,2),1,n0,4)
c$$$              endif
c$$$            endif
c$$$C         Reduce lmxb
c$$$            do  lp1 = lmxb+1, 1, -1
c$$$              if (idxdn(lp1,1).ne.4 .or. idxdn(lp1,nkapi).ne.4) goto 12
c$$$              lmxb = lmxb-1
c$$$            enddo
c$$$   12       continue
c$$$            jp = j2
c$$$            i=parg(' PZ=',4,sout,jp,ls,', ',3,lmxb+1,it,orbp(1,1,3))
c$$$            if (i .ne. 0) then
c$$$C           lmxi = lmxb
c$$$C           if (i .lt. 0) lmxi=-i-2
c$$$              i_copy_size=size(sspec(ks)%pz)
c$$$              call dcopy(i_copy_size,orbp(1,1,3),1,sspec(ks)%pz,1)
c$$$              lpz = 1
c$$$            endif
c$$$            i_copy_size=size(sspec(ks)%orbp)
c$$$            call dcopy(i_copy_size,orbp,1,sspec(ks)%orbp,1)
c$$$            i_copy_size=size(sspec(ks)%idxdn)
c$$$            call icopy(i_copy_size,idxdn,1,sspec(ks)%idxdn,1)
c$$$cc            sspec(ks)%orbp = orbp
c$$$cc            sspec(ks)%idxdn= idxdn
c$$$            sspec(ks)%lmxb = lmxb
c$$$   10     continue
c$$$        endif
c$$$        ioorbp = .true.
c$$$        if (lpz .eq. 1) nkmax = nkmax + 10
c$$$      else
c$$$        jfi = -ifi
c$$$        write(jfi,'(''BASIS:'')')
c$$$        if (mode .eq. 1) then
c$$$          do  110  is = is1, is2
c$$$c            do i_spacks=is,is
c$$$c              call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
c$$$c            enddo
c$$$            spid=sspec(is)%name
c$$$
c$$$            sspec(is)%norp=2
c$$$            sspec(is)%ntorb=n0*nkape
c$$$
c$$$            norp=sspec(is)%norp
c$$$            ntorb=sspec(is)%ntorb
c$$$            i_copy_size=size(sspec(is)%orbp)
c$$$            call dcopy(i_copy_size,sspec(is)%orbp,1,orbp,1)
c$$$            lmxb=sspec(is)%lmxb
c$$$
c$$$C         Sanity checks
c$$$            l_dummy_isanrg=isanrg(norp,2,2,'ioorbp:','norp',.true.)
c$$$            l_dummy_isanrg=isanrg(ntorb,n0*nkapi,n0*nkape,'ioorbp:','ntorb',.true.)
c$$$            etab=0d0 !call dpzero(etab,n0*2)
c$$$            rtab=0d0 !call dpzero(rtab,n0*2)
c$$$            sout = ' ' // spid
c$$$            do  112  ikap = 1, nkapi
c$$$              call dcopy(n0,orbp(1,1,ikap),1,rtab(1,ikap),1)
c$$$              call dcopy(n0,orbp(1,2,ikap),1,etab(1,ikap),1)
c$$$              if (dasum(n0,orbp(1,1,ikap),1) .ne. 0) then
c$$$                call awrit6('%a RSMH%?#n==2#2##=%n:1;3d EH%?#n==2#2##=%n:1;3d',
c$$$     .          sout,ls,0,
c$$$     .          ikap,lmxb+1,rtab(1,ikap),
c$$$     .          ikap,lmxb+1,etab(1,ikap))
c$$$              endif
c$$$  112       continue
c$$$            lpz = 0
c$$$            do  lp1 = 1, lmxb+1
c$$$              if (orbp(lp1,1,3) .gt. 0) lpz = lp1
c$$$            enddo
c$$$            if (lpz .gt. 0)
c$$$     .      call awrit2('%a PZ=%n:1;3d',sout,ls,0,lpz,orbp(1,1,3))
c$$$            call awrit0('%a',sout,-ls,-jfi)
c$$$  110     continue
c$$$        endif
c$$$      endif
c$$$      end function ioorbp
c$$$
c$$$
