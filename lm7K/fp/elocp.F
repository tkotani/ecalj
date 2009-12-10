      subroutine elocp(nbas,nsp,ssite,sspec,job)
C- Make envlope parameters for extended local orbitals
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :number of spin channels
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pnu ov0 pz
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: a nr rmt z lmxa lmxb rs3 eh3 vmtz
Ci     Stored:    orbp
Ci     Passed to: uspecb
Ci   job   :1s  digit
Ci         : 0 do nothing ; just return
Ci         : 1 make core and augmentation matrices
Co Outputs
Co   sspec->orbp:
Co         : smoothing radius and energy set for extended local orbitals
Cr Remarks
Cu Updates
Cu   06 Jul 05 first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nsp,job,n0
      parameter (n0=10)
      double precision ssite(1),sspec(1)
C ... Local parameters
      character spid*8
      logical eloc
      integer ib,ibs,iclbsj,ipr,iprint,is,k,l,lmxa,lmxb,nglob,nkap0,
     .  nr,nrmx,nrspec,nspec,stdo
      parameter (nrmx=1501, nkap0=3)
      integer lh(nkap0),nkape,nkaph,idamax
      integer oips,ov0
      double precision z,a,rmt,rs3,eh3,vmtz,xx
      double precision rofi(nrmx),vseli(4,n0),vsel(4,n0,nbas),
     .  pnu(n0,2),pnz(n0,2),eh(n0,nkap0),rsmh(n0,nkap0),
     .  pnui(n0,2),pnzi(n0,2),ehl(n0,nkap0),rsml(n0,nkap0)
C ... Heap
      integer w(1)
      common /w/ w

C --- Setup ---
      if (mod(job,10) .eq. 0) return
      call tcn('elocp')
      stdo = nglob('stdo')
      nspec = nglob('nspec')
      ipr = iprint()
      nkaph = nglob('nkaph')
      eloc = .false.
      call dpzero(vsel,4*n0*nbas)

C --- Find val, slo, K.E. for all sites ---
      do  ib = 1, nbas
        call upack('site spec pnu ov0',ssite,ib,is,pnu,ov0,0)
        call upack2('site pz',ssite,ib,pnz)
        call spacks(0,'spec name',sspec,spid,is,is)
        call upack('spec a nr rmt',sspec,is,a,nr,rmt,0)
        call upack('spec z lmxa lmxb',sspec,is,z,lmxa,lmxb,0)
        if (lmxa .eq. -1) goto 10
        if (pnz(idamax(lmxb+1,pnz,1),1) .lt. 10) goto 10
        eloc = .true.

        call radmsh(rmt,a,nr,rofi)
        call loctsh(1101,spid,z,a,nr,nr,nsp,lmxa,rofi,w(ov0),
     .    pnu,pnz,xx,xx,vmtz,vsel(1,1,ib),rsml,ehl)
   10   continue
      enddo

      if (.not. eloc) goto 999
      if (ipr .ge. 30) write(stdo,199)
  199 format(/' elocp:')

C --- Determine shape of smooth Hankel tails for local orbitals ---
      call defi(oips,nbas)
      call spackv(10,'site spec',ssite,1,nbas,w(oips))
C ... Loop over species containing extended local orbitals
      do  is = 1, nspec
        call spacks(0,'spec name',sspec,spid,is,is)
        call upack('spec z lmxa lmxb',sspec,is,z,lmxa,lmxb,0)
        if (lmxa .eq. -1) goto 20
        nrspec = iabs(iclbsj(is,w(oips),-nbas,nbas))
        if (nrspec .eq. 0) goto 20
        ib = iclbsj(is,w(oips),nbas,1)
        call upack('site pnu pz',ssite,ib,pnui,pnzi,0,0)
        if (pnzi(idamax(lmxb+1,pnzi,1),1) .lt. 10) goto 20

C   ... Average over sites within this species
        call dpzero(vseli,4*n0)
        do  ibs = 1, nrspec
          ib = iclbsj(is,w(oips),nbas,ibs)
          call upack2('site pz',ssite,ib,pnz)
          if (pnzi(idamax(lmxb+1,pnzi,1),1) .lt. 10) goto 22
          call dpadd(vseli,vsel(1,1,ib),1,4*n0,1/dble(nrspec))

   22     continue
        enddo

C   ... Printout of input for parameters
        if (ipr .ge. 90/1) then
        write(stdo,261) spid
  261   format(/'  l  site    Eval        Val         Slo         K.E.',
     .    5x,'species ',a)
        do  l = 0, lmxb
          if (pnz(l+1,1) .lt. 10) goto 24
          do  ibs = 1, nrspec
            ib = iclbsj(is,w(oips),nbas,ibs)
            write (stdo,260)
     .        l,ib,vsel(4,l+1,ib),(vsel(k,l+1,ib),k=2,4)
  260       format(i3,i4,4f12.6:a)
  262       format(i3,' avg',4f12.6)
          enddo
          write (stdo,262) l,vseli(4,l+1),(vseli(k,l+1),k=2,4)

   24     continue
        enddo
        endif

C   ... Make parameters for this species
        call dpzero(rsmh,n0*nkap0)
        call dpzero(eh,n0*nkap0)
        call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkape)
        call upack('spec rs3 eh3 vmtz',sspec,is,rs3,eh3,vmtz,0)
        call upack('spec a nr rmt',sspec,is,a,nr,rmt,0)
        call dpzero(rsml,n0*2)
        call dpzero(ehl,n0*2)
        call radmsh(rmt,a,nr,rofi)
        call loctsh(1102,spid,xx,a,nr,nr,nsp,lmxa,rofi,w,
     .    pnui,pnzi,rs3,eh3,vmtz,vseli,rsml,ehl)
        call dcopy(n0,rsml,1,rsmh(1,nkaph),1)
        call dcopy(n0,ehl, 1,eh(1,nkaph),1)
        call uspecb(1,1,sspec,is,is,lh,rsmh,eh,nkaph)

   20   continue
      enddo

  999 continue
      call tcx('elocp')

      end


