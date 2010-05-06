Cgetarg...info...           structure ['nmpot', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['nmpot', 'ham', 'sham']
Cgetarg...info...           structure ['nmpot', 'site', 'ssite']
Cgetarg...info...           structure ['nmpot', 'pot', 'spot']
Cgetarg...info...           structure ['nmpot', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['nmpot', 'm_struc_def']
      subroutine nmpot(job,sctrl,sspec,ssite,sham,spot,lidim,lihdim,
     .  iprmb,dclabl,ppn)
       
       use m_struc_def  !Cgetarg

C- Generates NMTO potential parameters
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :1s digit
Ci         :0 Make ppar for each of ham->nmto energies
Ci         :1 Make ppar for one energy, b.c.'s supplied by pnu
Ci         :  (not implemented)
Ci         :10s  digit
Ci         :0 Read potential from class file (requires dclabl)
Ci         :1 Use potential from site->ov0
Ci         :4 add 4 to 10s digit if pot->ves(ic) should be added
Ci         :  to spherical potential
Ci         :100s  digit
Ci         :1 Use hcr = rmt, rather than spec->hcr
Ci         :  NB: 2nd gen LMTO uses this switch
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbasp nclass nspec nspin
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: a nr z rmt lmxa hcr
Ci     Stored:    *
Ci     Passed to: *
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: pnu spec class ov0
Ci     Stored:    *
Ci     Passed to: *
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: nmto kmto
Ci     Stored:    *
Ci     Passed to: *
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: vmtz oves
Ci     Stored:    *
Ci     Passed to: *
Ci
Ci   lidim :number of lower+intermediate orbitals
Ci
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci
Ci   dclabl:class name, packed as a real number
Ci         :dclabl(0) => classes are not named
Ci
Co Outputs
Co   ppn   :3rd generation potential parameters, in downfolding order
Cu Updates
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters:
      integer job,lidim,lihdim,iprmb(1)
Cgetarg       double precision sctrl(1),sspec(1),ssite(1),sham(1),spot(1)
       type(s_ctrl)::sctrl
       type(s_spec)::sspec(*)
       type(s_site)::ssite(*)
       type(s_ham)::sham
       type(s_pot)::spot

      integer nppn
      parameter (nppn=12)
      double precision ppn(nppn,lihdim,*),dclabl(*)
C Local parameters:
      logical sw,aiopot
      integer fopna,ib,ic,ifi,igetss,ii,ik,iks,is,isp,jb,jj,l,lmr,lmx,
     .  m,n0,nbasp,nclass,nglob,nkaph,nl,nmto,nr,nsp,nspec,offji,job0,
     .  job1,job2
      integer orofi,ov0,oibc,iprint,lgunit,oves,jobx
      parameter(n0=10)
      double precision z,avw,a,b,hcr(n0),rmt,kmto(10),vmtz,dglob
      double precision ppnl(nppn,n0,n0,2),ves,pnu(n0,2)
C     double complex kmtoz(10)
      character*8 clabl
C Heap:
      integer w(1)
      common /w/ w

C --- Unpack variables and other setup ---
      job0 = mod(job,10)
      job1 = mod(job/10,10)
      job2 = mod(job/100,10)
Cgetarg...info...           integer :: nbasp , integer(8) :: sctrl%nbasp 1
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg...info...           integer :: nspec , integer(8) :: sctrl%nspec 1
Cgetarg...info...           integer :: nsp , integer(8) :: sctrl%nspin 1
Cgetarg       call upack('ctrl nbasp nclass nspec nspin',sctrl,nbasp,nclass,
Cgetarg      .  nspec,nsp,0)
       
       nbasp=sctrl%nbasp
       nclass=sctrl%nclass
       nspec=sctrl%nspec
       nsp=sctrl%nspin


C     Unpack other job-specific parameters
      if (job0 .eq. 0) then
Cgetarg...info...           integer :: nmto , integer(8) :: sham%nmto 1
Cgetarg...info...           real(8) :: kmto(10) , real(8) :: sham%kmto 6
Cgetarg         call upack('ham nmto kmto',sham,nmto,kmto,0,0,0)
         
         nmto=sham%nmto
         i_copy_size=size(sham%kmto) 
         call dcopy(i_copy_size,sham%kmto,1,kmto,1) 

        jobx = 2
      elseif (job0 .eq. 1) then
        call rxi('nmpot: ill defined job',job)
Cgetarg...info...           real(8) :: pnu(n0,2) , real(8) :: ssite%pnu 20
Cgetarg         call upack2('site pnu',ssite,ib,pnu)
         
         i_copy_size=size(ssite(ib)%pnu) 
         call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1) 

        jobx = 0
        nmto = 1
      else
        call rxi('nmpot: ill defined job',job)
      endif
Cgetarg...info...           real(8) :: vmtz , real(8) :: spot%vmtz 1
Cgetarg...info...           integer :: oves , integer(8) :: spot%oves 1
Cgetarg       call upack('pot vmtz oves',spot,vmtz,oves,0,0,0)
       
       vmtz=spot%vmtz
       oves=spot%oves

      ves = 0

      nkaph = nglob('nkaph')
      avw = dglob('avw',0d0,0)
      nl = nglob('nl')
      nbasp = nglob('nbasp')
      call setcc(nglob('lrel'))
      if (lidim.ne.lihdim) call rx('nmpot2 not set up for downfolding')
      call dpzero(ppn, nppn*lihdim*nmto*nsp)
      call dpzero(ppnl,nppn*n0*n0*2)
      call defi(oibc,-nbasp)
      if (iprint() .ge. 20) then
        if (job0 .eq. 0) call awrit2('%N NMPOT: parameters for %i%-1j'//
     .    ' kinetic energies :%n:2,1d',' ',80,lgunit(1),nmto,kmto)
        if (job0 .eq. 1) call awrit0('%N NMPOT: parameters for given'//
     .    ' boundary conditions ',' ',80,lgunit(1))
        if (iprint() .ge. 30) then
          write(lgunit(1),301)
  301     format(' l,E',6x,
     .      ' hcr     <phi phi>       a*D',
     .      '        a*Ddot      phi(a)      phip(a)')
        else
          write(lgunit(1),'(1x)')
        endif
      endif

C --- For each class, get potential parameters ---
      lmr = 0
      do  ib = 1, nbasp
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg         ic = igetss('site class',ib,ssite)
         ic = int(ssite(ib)%class) 


C   ... Make pp's for new class ic
        if (w(oibc+ic-1) .eq. 0) then
Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg           call upack('spec a nr z rmt',sspec,is,a,nr,z,rmt)
           
           a=sspec(is)%a
           nr=sspec(is)%nr
           z=sspec(is)%z
           rmt=sspec(is)%rmt

Cgetarg...info...           integer :: lmx , integer(8) :: sspec%lmxa 1
Cgetarg...info...           real(8) :: hcr(n0) , real(8) :: sspec%hcr 10
Cgetarg           call upack('spec lmxa hcr',sspec,is,lmx,hcr,0,0)
           
           lmx=sspec(is)%lmxa
           i_copy_size=size(sspec(is)%hcr) 
           call dcopy(i_copy_size,sspec(is)%hcr,1,hcr,1) 

          if (job2 .eq. 1) then
            call dvset(hcr,1,n0,rmt)
          endif
          call defdr(orofi,nr*2)
          call radmsh(rmt,a,nr,w(orofi))
          b = rmt/(dexp(a*(nr-1)) - 1)

C         Get potential for parameters, depending on 1s digit job
          clabl = ' '
          if (mod(job1,4) .eq. 0) then
            call r8tos8(dclabl(ic),clabl)
            ifi = fopna(clabl,-1,0)
            call defdr(ov0,nr*nsp)
            sw = aiopot(nr,nsp,a,rmt,-99d0,w(ov0),ifi)
            call fclose(ifi)
            if (.not. sw) call rxi('nmpot: no potential for class',ic)
          elseif (mod(job1,4) .eq. 1) then
Cgetarg...info...           integer :: ov0 , integer(8) :: ssite%ov0 1
Cgetarg             call upack2('site ov0',ssite,ib,ov0)
             
             ov0=ssite(ib)%ov0

          else
            call rxi('nmpot: not implemented for job',job)
          endif

          if (job1 .ge. 4) then
            call dpscop(w(oves),ves,1,ic,1,1d0)
          else
          endif

          call nmpot2(jobx,a,avw,b,ib,ic,is,clabl,kmto,pnu,lmx,nmto,nl,
     .      nr,nsp,w(orofi),hcr,w(ov0),vmtz,ves,z,ppnl)
C          kmtoz = kmto
C          call nmpot3(jobx,a,avw,b,ib,ic,is,clabl,kmtoz,pnu,lmx,nmto,nl,
C     .      nr,nsp,w(orofi),hcr,w(ov0),vmtz,ves,z,ppnl)

          do  l = 0, nl-1
          if (iprmb(lmr+1) .gt. lihdim) then
            lmr = lmr + 2*l+1
          else
            do   m = -l, l
              lmr = lmr+1
              ii = iprmb(lmr)
              do  ik = 1, nmto
                do  isp = 1, nsp
                  iks = ik + nmto*(isp-1)
                  call dcopy(nppn,ppnl(1,l+1,ik,isp),1,ppn(1,ii,iks),1)
                enddo
              enddo
            enddo
          endif
          enddo
          if (mod(job1,4) .eq. 0) call rlse(ov0)
          call rlse(orofi)
          w(oibc+ic-1) = ib

C   ... Copy from prior site jb with equivalent class
        else
          jb = w(oibc+ic-1)
          offji = nl*nl*nkaph*(jb-ib)
          do  l = 0, nl-1
          if (iprmb(lmr+1) .gt. lihdim) then
            lmr = lmr + 2*l+1
          else
            do   m = -l, l
              lmr = lmr+1
              ii = iprmb(lmr)
              jj = iprmb(lmr+offji)
              do  ik = 1, nmto
                do  isp = 1, nsp
                  iks = ik + nmto*(isp-1)
                  call dcopy(nppn,ppn(1,jj,iks),1,ppn(1,ii,iks),1)
                enddo
              enddo
            enddo
          endif
          enddo

        endif

      enddo

      call rlse(oibc)

      if (iprint() .ge. 90)
     .  call prm('ppn',0,ppn,nppn*lihdim,nppn*lihdim,nmto*nsp)

Cgetarg       end
       end subroutine nmpot 

      subroutine nmpot2(job,a,avw,b,ib,ic,is,clabl,kmto,pnu,lmx,ne,nl,
     .  nr,nsp,rofi,hcr,v,vmtz,ves,z,ppnl)

C- Generates potential parameters for pot. functions and downfolding
C ----------------------------------------------------------------------
Ci Inputs:
Ci   job   :passed to phidx.  Use
Ci         :0, boundary conditions specified val,slo
Ci         :2, boundary conditions specified by energy e (see Remarks)
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   ib    :site index, for printout only
Ci   ic    :class index, for printout only
Ci   is    :species index, for printout only
Ci   clabl :class name, used for printout only
Ci   kmto  :set of kinetic energies for which to make ppars
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci         :pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci         :(not used now)
Ci   lmx   :l-cutoff
Ci   ne    :number of kmto kinetic energies
Ci   nl    :number of l's (for dimension of arrays)
Ci   nr    :number of mesh points
Ci   nsp   :=1 spin degenerate, =2 non-degenerate
Ci   rofi  :radial mesh points
Ci   hcr   :hard sphere radii
Ci   v     :spherical potential = v_true
Ci   vmtz  :muffin-tin zero
Ci   ves   :constant potential shift
Ci   z     :nuclear charge
Co Outputs:
Co   ppnl  :for each energy,l,isp:
Co         :(1) = inverse potential function
Co         :(2) = normalization of phi
Co         :(3) = a * log derivative of phi at a = hcr
Co         :(4) = a * log derivative of phidot at a = hcr
Co         :(5) = value of wave function at a = hcr
Co         :(6) = value of phidot at a = hcr
Co         :(7) = normalization of phidot (not computed now)
Co         :(8) = not used now
Cr Remarks:
Cr   This was adapted from the Stuttgart LMTO package.
C  ---------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer job,ib,ic,is,ne,nl,nsp,lmx,nr,nppn,n0
      parameter (nppn=12,n0=10)
      double precision a,avw,b,ppnl(nppn,n0,n0,nsp),kmto(ne),pnu(n0,2),
     .  rofi(nr),hcr(nl),v(nr,nsp),vmtz,ves,z
      character*8 clabl
C Local variables:
      logical liophi
      integer ie,k,iprint,isp,il,ir,l,nn,lgunit,og,ogp,awrite
      double precision alpha,dla,dlap(4),dum,ptf,s,vali(5),sloi(5),rmax,
     .  phia,phiap(4),phio,ei,tol,phi,dphi,phip,dphip,p,hcrl
      parameter (tol=1d-12)
      character*80 outs
C Heap:
      integer w(1)
      common /w/ w
C External calls:
      external awrit0,defdr,gintsr,makptf,norm2g,phidx

C     Add constant potential to satisfy boundary conditions
      do  isp = 1, nsp
      do  ir = 1, nr
        v(ir,isp) = v(ir,isp) + ves
      enddo
      enddo

      liophi = .false.
      rmax = rofi(nr)
      call defdr(og,2*nr)
      call defdr(ogp,8*nr)

C     outs = description of this class, depending on info available
      if (iprint() .ge. 30) then
        outs = ' '
        il = awrite(' %?#n#%-1jsite %i, ##'//
     .              '%?#n#%-1jclass %i, ##'//
     .              '%?#n#%-1jspecies %i, ##'//
     .              '%b%b :',outs,len(outs),0,ib,ic,is,0,0,0,0,0)
        if (clabl .ne. ' ') outs(il+2:) = clabl
        call awrit0(outs,' ',-len(outs),lgunit(1))
      endif
      do  isp = 1, nsp
      do  ie = 1, ne
      if (job .eq. 2) then
        ei = kmto(ie)+vmtz
      elseif (job .eq. 2) then
        call rx('not ready')
C       nn = int(pnu(l+1,isp)) - l - 1
      else
        call rxi('nmpot2: nonsensical job',job)
      endif
      do  l = 0, lmx
        il = l+1
        hcrl = hcr(il)

        call phidx(job,z,l,v(1,isp),hcrl,vmtz,rofi,nr,4,tol,ei,vali,
     .    sloi,nn,w(og),w(ogp),phi,dphi,phip,dphip,p,phia,phiap,dla,
     .    dlap)
C       dphip = (sloi(2)-phip)/rmax
        call makptf(avw,kmto(ie),l,hcrl,sloi,vali,rmax,ptf,alpha)
C       Scale wave function to make phi^0(a)=1
        call norm2g(a,b,kmto(ie),w(og),l,liophi,nr,hcrl,sloi,
     .    vali,rmax,rofi,dum,phio)
C       s is resulting overlap
        call gintsr(w(og),w(og),a,b,nr,z,ei,l,v(1,isp),rofi,s)
        ppnl(1,il,ie,isp) = 1d0/ptf
        ppnl(2,il,ie,isp) = s
        ppnl(3,il,ie,isp) = hcrl*dla
        ppnl(4,il,ie,isp) = hcrl*dlap(1)
        ppnl(5,il,ie,isp) = phia
        ppnl(6,il,ie,isp) = phiap(1)

C       if (iprint() .ge. 30) write(lgunit(1),301)
C    .    l,ei,hcrl,(ppnl(k,il,ie,isp),k=1,4)
        if (iprint() .ge. 30) write(lgunit(1),301)
     .    l,ie,hcrl,(ppnl(k,il,ie,isp),k=2,6)
  301   format(2i2,2f12.6,f13.6,10f12.6)
C       call awrit5('%,2i%,2i    %,6#12d%4,6;7#12g',' ',80,
C    .    lgunit(1),l,ie,hcrl,ppnl(2,il,ie,isp),0)

C       call prrmsh('g',rofi,w(og),nr,nr,1)
        enddo
      enddo
      enddo
      call rlse(og)

C     Undo constant potential shift
      do  isp = 1, nsp
      do  ir = 1, nr
        v(ir,isp) = v(ir,isp) - ves
      enddo
      enddo

Cgetarg       end
       end subroutine nmpot2 

C      subroutine nmpot3(job,a,avw,b,ib,ic,is,clabl,kmto,pnu,lmx,ne,nl,
C     .  nr,nsp,rofi,hcr,v,vmtz,ves,z,ppnl)
CC- Potential parameters for pot. functions and downfolding, complex e
CC ----------------------------------------------------------------------
CCi Inputs:
CCi   job   :passed to phidz.  Use
CCi         :0, boundary conditions specified val,slo
CCi         :2, boundary conditions specified by energy e (see Remarks)
CCi   b     :                 -//-
CCi   ib    :site index, for printout only
CCi   ic    :class index, for printout only
CCi   is    :species index, for printout only
CCi   clabl :atom name, for printout only
CCi   lmx   :l-cutoff
CCi   nl    :number of l's (for dimension of arrays)
CCi   nr    :number of mesh points
CCi   nsp   :=1 spin degenerate, =2 non-degenerate
CCi   rofi  :radial mesh points
CCi   hcr   :hard sphere radii
CCi   v     :spherical potential = v_true
CCi   vmtz  :muffin-tin zero
CCi   ves   :constant potential shift
CCi   z     :nuclear charge
CCo Outputs:
CCo   ppnl  :for each energy,l,isp:
CCo         :(1) = inverse potential function
CCo         :(2) = normalization of phi
CCo         :(3) = a * log derivative of phi at a = hcr
CCo         :(4) = a * log derivative of phidot at a = hcr
CCo         :(5) = value of wave function at a = hcr
CCo         :(6) = value of phidot at a = hcr
CCo         :(7) = normalization of phidot (not computed now)
CCo         :(8) = not used now
CCr Remarks:
CCr   This is an analog of nmpot2 for complex energies.
CC  ---------------------------------------------------------------------
C      implicit none
CC Passed variables:
C      integer job,ib,ic,is,ne,nl,nsp,lmx,nr,nppn,n0
C      parameter (nppn=12,n0=10)
C      double precision a,avw,b,pnu(n0,2),
C     .  rofi(nr),hcr(nl),v(nr,nsp),vmtz,ves,z
C      double complex kmto(ne),ppnl(nppn,n0,n0,nsp)
C      character*8 clabl
CC Local variables:
C      logical liophi
C      integer ie,k,iprint,isp,il,ir,l,nn,lgunit,og,ogp,awrite
C      double precision dum,rmax,tol,hcrl
C      double complex ei,sloi,vali,phia,phiap(4),phio,phi,dphi,phip,
C     .  dphip,p,dla,dlap(4),ptf,alpha,s
C      parameter (tol=1d-12)
C      character*80 outs
CC Heap:
C      integer w(1)
C      common /w/ w
CC External calls:
C      external awrit0,defcc,gintz,makpzf,nrm2gz,phidz
C
CC     Add constant potential to satisfy boundary conditions
C      do  isp = 1, nsp
C      do  ir = 1, nr
C        v(ir,isp) = v(ir,isp) + ves
C      enddo
C      enddo
C
C      liophi = .false.
C      rmax = rofi(nr)
C      call defcc(og,2*nr)
C      call defcc(ogp,8*nr)
C
CC     outs = description of this class, depending on info available
C      if (iprint() .ge. 30) then
C        outs = ' '
C        il = awrite(' %?#n#%-1jsite %i, ##'//
C     .              '%?#n#%-1jclass %i, ##'//
C     .              '%?#n#%-1jspecies %i, ##'//
C     .              '%b%b :',outs,len(outs),0,ib,ic,is,0,0,0,0,0)
C        if (clabl .ne. ' ') outs(il+2:) = clabl
C        call awrit0(outs,' ',-len(outs),lgunit(1))
C      endif
C      do  isp = 1, nsp
C      do  ie = 1, ne
C      if (job .eq. 2) then
C        ei = kmto(ie)+vmtz
C      elseif (job .eq. 2) then
C        call rx('not ready')
CC       nn = int(pnu(l+1,isp)) - l - 1
C      else
C        call rxi('nmpot2: nonsensical job',job)
C      endif
C      do  l = 0, lmx
C        il = l+1
C        hcrl = hcr(il)
C        call phidz(job,z,l,v(1,isp),hcrl,vmtz,rofi,nr,4,tol,ei,vali,
C     .    sloi,nn,w(og),w(ogp),phi,dphi,phip,dphip,p,phia,phiap,dla,
C     .    dlap)
C        call makpzf(avw,kmto(ie),l,hcrl,sloi,vali,rmax,ptf,alpha)
CC       Scale wave function to make phi^0(a)=1
C        call nrm2gz(a,b,kmto(ie),w(og),l,liophi,nr,hcrl,sloi,
C     .    vali,rmax,rofi,dum,phio)
CC       s is resulting overlap
C        call gintz(w(og),w(og),a,b,nr,z,ei,l,v(1,isp),rofi,s)
C        ppnl(1,il,ie,isp) = 1d0/ptf
C        ppnl(2,il,ie,isp) = s
C        ppnl(3,il,ie,isp) = hcrl*dla
C        ppnl(4,il,ie,isp) = hcrl*dlap(1)
C        ppnl(5,il,ie,isp) = phia
C        ppnl(6,il,ie,isp) = phiap(1)
C
CC       if (iprint() .ge. 30) write(lgunit(1),301)
CC    .    l,ei,hcrl,(ppnl(k,il,ie,isp),k=1,4)
C        if (iprint() .ge. 30) write(lgunit(1),301)
C     .    l,ie,hcrl,(ppnl(k,il,ie,isp),k=2,6)
C  301   format(2i2,2f12.6,f13.6,10f12.6)
CC       call awrit5('%,2i%,2i    %,6#12d%4,6;7#12g',' ',80,
CC    .    lgunit(1),l,ie,hcrl,ppnl(2,il,ie,isp),0)
C
CC       call prrmsh('g',rofi,w(og),nr,nr,1)
C        enddo
C      enddo
C      enddo
C      call rlse(og)
C
CC     Undo constant potential shift
C      do  isp = 1, nsp
C      do  ir = 1, nr
C        v(ir,isp) = v(ir,isp) - ves
C      enddo
C      enddo
C
C      end

