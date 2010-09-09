      subroutine uspecb(lpack,mode,sspec,is1,is2,lh,rsmh,eh,nkape)

      use m_struc_def  !Cgetarg
      use m_globalvariables
C- Pack/unpack parameters related to basis
C ----------------------------------------------------------------------
Ci Inputs
Ci   lpack :0 for unpack, 1 for pack
Ci         :NB packing case, it is the caller's responsibility
Ci         :to ensure rsmh and eh are nonzero only for l and energies
Ci         :1..nkape which have envelope functions.
Ci   mode  :At present, this is the only implemented mode
Ci         :-1 Only forces uspecb to read nkaph from global variable
Ci         :0 No reference to envelope function parameters.
Ci         :  Only nkape is returned; see below
Ci         :1,2,3: parameters consist of rsmh,eh, where:
Ci              rsmh,eh         first orbital in basis
Ci              rsmh2,eh2       if nkap=2
Ci              rsmhl,ehl       if extended local orbitals
Ci            Data taken from spec->orbp (see Remarks)
Ci            modes 1..4 are similar except that:
Ci            mode=1 does not return envelope information about
Ci                   local orbitals
Ci            mode=2 returns rsmh,eh in last channel corresponding to
Ci                   loc. orb. (rsmh>0 for extended loc. orbs; else 0)
Ci                   and sets nkape to nkaph
Ci            mode=3 intended for local orbitals of the second type
Ci                   NOT IMPLEMENTED
Ci            mode=4 only returns envelope information about
Ci                   local orbitals in rsmh(*,1),eh(*,1)
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxb pz orbp
Ci     Stored:    name orbp
Ci     Passed to: spacks
Ci  is1,is2:range of species for which to extract parameters
Cio Inputs/Outputs
Co   ... All of the following are generated for lpack=0
Co   lh    :lh(1..nkape) = max l for which envelope basis fn defined
Co         :Also lh(nkaph) = max l for which a local orbital defined
Co         :lh(1+nkape..nkap0) are initialized to -1
Cio  rsmh  :rsmh(1..l+1,1..nkape,is1..is2), rsmh(1..l+1,nkaph,is1..is2)
Cio        :smoothing radii for envelope functions, extended loc. orbitals
Cio        :(rsmh is input for lpack=1)
Cio  eh    :eh(1..l+1,1..nkape,is1..is2), eh(1..l+1,nkaph,is1..is2)
Cio        :energies of envelope functions and extended local orbitals
Cio        :(eh is input for lpack=1)
Cio  nkape :number of envelope function types per l q.n. for is1..is2
Cio        :1s digit contains number of energy channels in rsmh,eh
Cio        :which contain valence orbitals with envelope functions.
Cio        :nkape may be modified depending on mode:
Cio        :mode=1: nkape not affect by local orbitals
Cio        :mode=2: 1 added to nkape if local orbitals are included.
Co         :When called with mode=0, nkape returns the following:
Co         :1s    digit always zero
Co         :10s   digit 1 if local orbital of first type exists
Co         :100s  digit 1 if local orbital of 2nd type exists
Cl Local variables
Cl   npzi  :1 if there are extended local orbitals AND mode>1
Cl         :0 otherwise
Cl   nkapi :number of envelope functions EXCLUDING
Cl         :extended local orbitals
Cr Remarks
Cr   Case parameters contained in sspec->orbp
Cr     orbp is dimensioned orbp(n0,2,nkap0)
Cr     1st index for orbital types ?
Cr     2nd index :1 for rsmh, 2 for eh
Cr     3rd index :1 for (rsmh,eh) 2 for (rsmh2,eh2), 3 loc orb.
Cr     Or, if no rsmh2,eh2, loc. orb occupy 2nd spot
Cu Updates
Cu    8 Jul 05 New mode=-1
Cu   25 Jul 04 Redesigned for extended local orbitals
Cu   10 Apr 02 Redimensioned eh,rsmh to accomodate larger lmax
Cu   19 Feb 02 mode 0 implemented.  Convention for npqni changed.
Cu   12 Feb 02 npqn,nkape are returned for range is1..is2
Cu   25 Aug 01 Extended to return npqn,lh for local orbitals
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nkap0,n0,nkape
      parameter (nkap0=3,n0=10)
      integer lpack,mode,is1,is2,lh(nkap0,is1:is2)
      real(8):: rsmh(n0,nkap0,is1:is2) , eh(n0,nkap0,is1:is2)
      type(s_spec)::sspec(*)

C ... Local parameters
      character spid*8
      integer is,lmxb,ikap,kkap,nkaph,nglob,nkapi,i,j,k,l,npzi,lpzi,
     .lpzl,npze,lmxa
      double precision orbp(n0,2,nkap0),dasum,pnz(n0,2)
      save nkaph
      data nkaph /-1/

Changenglob      if (nkaph .eq. -1) nkaph = nglob('nkaph')
      if (nkaph .eq. -1) nkaph = globalvariables%nkaph

      if (mode .eq. -1) then
Changenglob        nkaph = nglob('nkaph')
        nkaph = globalvariables%nkaph
        return
      else if (mode .eq. 0) then
        i = 0
        j = 0
        l = 0
        do  is = is1, is2

          lmxa=sspec(is)%lmxa
          lmxb=sspec(is)%lmxb
          i_copy_size=size(sspec(is)%pz)
          call dcopy(i_copy_size,sspec(is)%pz,1,pnz,1)

          if (lmxa .gt. -1) then
            do  k = 1, lmxb+1
              if (pnz(k,1) .gt.  0) i = 1
              if (pnz(k,1) .ge. 10) j = 1
              if (pnz(k,1) .ge. 20) l = 1
            enddo
          endif
        enddo
        nkape = 10*i+100*j+1000*l

      elseif (mode .le. 4) then


        if (mode .eq. 3) call rx('uspecb: mode=3 not implemented')

        if (lpack .eq. 0) then
          nkape = 0
          npze = 0
        endif
        do  is = is1, is2
          do i_spacks=is,is
            call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
          enddo

          call dpzero(orbp,n0*nkap0*2)

          i_copy_size=size(sspec(is)%orbp)
          call dcopy(i_copy_size,sspec(is)%orbp,1,orbp,1)
          lmxb=sspec(is)%lmxb
          i_copy_size=size(sspec(is)%pz)
          call dcopy(i_copy_size,sspec(is)%pz,1,pnz,1)


C   --- Unpack ---
          if (lpack .eq. 0) then

            do  ikap = 1, nkap0
              kkap = ikap
              if (mode .eq. 4) kkap = 1
              lh(kkap,is) = -1
            enddo

C     ... make nkapi,npzi
            lpzi = 0
            do  k = 1, lmxb+1
              lpzl = 0
              if (pnz(k,1) .gt.  0) lpzl = 1
              if (pnz(k,1) .ge. 10) lpzl = 2
              if (pnz(k,1) .ge. 20) lpzl = 3
              lpzi = max(lpzi,lpzl)
            enddo
            npzi = 0
            if (lpzi .gt. 1 .and. (mode.ge.2 .and. mode.le.4)) npzi = 1

C     ... Count number of envelope energy channels nkapi (<=nkaph).
C         There are the following possibilities:
C         1. energy channel nkaph is a local orbital
C            1a.  local orbital for this species.  Then lpzi>0
C            1b.  local orbital for another species.  Then all rsmh=0
C                 for channel nkaph.
C         2. This species has no envelopes in energy channel nkapi.
C            Then all rsmh=0
C         Handle case 1a first:
            nkapi = nkaph
            if (lpzi .ne. 0) nkapi = nkaph-1
C         Cases 1b and 2 both handled by checking for rsmh=0 in channel.
   10       continue
            if (dasum(n0,orbp(1,1,nkapi),1) .eq. 0) then
              nkapi = nkapi - 1
              if (nkapi .gt. 0) goto 10
            endif
C          print *, '!!'
C          call dpzero(orbp(1,1,1),n0*2)

C     ... unpack rsmh,eh; make lh
            if (mode .eq. 4) then
              call dpzero(eh(1,1,is),n0)
              call dpzero(rsmh(1,1,is),n0)
              nkapi = 0
            else
              call dpzero(eh(1,1,is),n0*nkap0)
              call dpzero(rsmh(1,1,is),n0*nkap0)
            endif
            do  ikap = 1, nkaph
              if (ikap.le.nkapi .or. ikap.eq.nkaph .and. npzi.ne.0) then
                kkap = ikap
                if (mode .eq. 4) kkap = 1
                call dcopy(n0,orbp(1,1,ikap),1,rsmh(1,kkap,is),1)
                call dcopy(n0,orbp(1,2,ikap),1,eh(1,kkap,is),1)
                do  l = lmxb, 0, -1
                  if (rsmh(l+1,kkap,is) .gt. 0) then
                    lh(kkap,is) = l
                    if (eh(l+1,kkap,is) .ge. 0) call rxi('species '//spid//
     .              '%a has positive Hankel energy for l=',l)
                    goto 116
                  endif
                enddo
  116           continue
              endif
            enddo
C         If no channel has orbital, eliminate kappa
            do  ikap = nkapi, 1, -1
              if (lh(ikap,is) .eq. -1) nkapi = nkapi-1
            enddo

C     ... lh for local orbital
            kkap = nkaph
            if (mode .eq. 4) kkap = 1
            do  l = 0, lmxb
              if (pnz(l+1,1) .ne. 0) then
                lh(kkap,is) = l
              endif
            enddo

            nkape = max(nkape,nkapi)
            npze = max(npze,npzi)

C   --- Pack ---
          else
            if (mode .ge. 2) call rx('mode>1 not allowed for write')
            nkapi = mod(nkape,10)
            npzi = 0
            if (nkape .gt. 10) npzi = 1
            kkap = 0
            do  ikap = 1, nkaph
              if (ikap.le.nkapi .or. ikap.eq.nkaph .and. npzi.ne.0) then
                kkap = kkap+1
                call dcopy(n0,rsmh(1,kkap,is),1,orbp(1,1,ikap),1)
                call dcopy(n0,eh(1,kkap,is),  1,orbp(1,2,ikap),1)
              endif
            enddo

C          do  ikap = 1, nkaph
C            write(*,'(i3,''  RSMH='',4f9.4,''  EH='',4f9.4)')
C     .        ikap,(orbp(l,1,ikap),l=1,4),(orbp(l,2,ikap),l=1,4)
C          enddo


            i_copy_size=size(sspec(is)%orbp)
            call dcopy(i_copy_size,orbp,1,sspec(is)%orbp,1)


          endif

        enddo

        if (lpack .eq. 0) then
          if (npze.ne.0 .and. mode .eq. 2) nkape=nkaph
          if (npze.ne.0 .and. mode .eq. 4) nkape=1
        endif

      else
        call rxi('uspecb: mode not implemented, mode',mode)
      endif

      end subroutine uspecb


