      module m_elocp
      use m_lmfinit,only: nspec,n0
      real(8),allocatable,protected,public:: ehl(:,:),rsml(:,:)
      public:: elocp
      
      private
      contains
      subroutine elocp() 
      use m_struc_def  
      use m_lmfinit,only: nkaph,stdo,nspec,nbas,nsp,ssite=>v_ssite,sspec=>v_sspec,n0,nkapii!,ex=>ehl,rx=>rsml
c      use m_uspecb,only: uspecb
C- Make envlope parameters for extended local orbitals
C ----------------------------------------------------------------------
Ci Inputs
Ci         : 0 do nothing ; just return
Ci         : 1 make core and augmentation matrices
Co Outputs
Co   ehl,rsml
Co         : smoothing radius and energy set for extended local orbitals
Cr Remarks
Cu Updates
Cu   06 Jul 05 first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
c      integer n0
c      parameter (n0=10)hl
cki      type(s_site)::ssite(*)
c      type(s_site)::ssite(nbas)
c      type(s_spec)::sspec(*)

C ... Local parameters
      character spid*8
      logical eloc
      integer ib,ibs,iclbsj,ipr,iprint,is,k,l,lmxa,lmxb,nglob,nkap0,
     .nr,nrmx,nrspec
      parameter (nrmx=1501, nkap0=3)
      integer:: nkape,idamax,nkapex !lh(nkap0),
      integer ,allocatable :: ips_iv(:)
      double precision z,a,rmt,rs3,eh3,vmtz,xx
      double precision rofi(nrmx),vseli(4,n0),vsel(4,n0,nbas),
     .pnu(n0,2),pnz(n0,2),eh(n0,nkap0),rsmh(n0,nkap0),
     .     pnui(n0,2),pnzi(n0,2)
      real(8):: ehls(n0*2),rsmls(n0*2)
      integer ::iwdummy,i_copy_size,i_spackv
      real(8):: wdummy(1)
C --- Setup ---
      call tcn('elocp')
      if(allocated(ehl)) deallocate(ehl,rsml)
      allocate( ehl(n0,nspec),rsml(n0,nspec))
!!      
c      if(allocated(ex)) deallocate(ex,rx)
c      allocate( ex(n0,nspec),rx(n0,nspec))
      
      ipr = iprint()
      eloc = .false.
      vsel=0d0
C --- Find val, slo, K.E. for all sites ---
      do  ib = 1, nbas
        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pnu)
        call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
        i_copy_size=size(ssite(ib)%pz)
        call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)
        spid = sspec(is)%name
        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt
        z=sspec(is)%z
        lmxa=sspec(is)%lmxa
        lmxb=sspec(is)%lmxb
        if (lmxa .eq. -1) goto 10
        if (pnz(idamax(lmxb+1,pnz,1),1) .lt. 10) goto 10
        eloc = .true.
        call radmsh(rmt,a,nr,rofi)
        call loctsh ( 1101 , spid , z , a , nr , nr , nsp , lmxa , rofi
     .       , ssite(ib)%rv_a_ov0 , pnu , pnz , xx , xx , vmtz , vsel ( 1 , 1 , ib ) 
     .       , rsml , ehl )
   10   continue
      enddo
      if (.not. eloc) goto 999
      if (ipr .ge. 30) write(stdo,199)
  199 format(/' elocp:')
C --- Determine shape of smooth Hankel tails for local orbitals ---
      allocate(ips_iv(nbas))
      i_copy_size=1;
      do i_spackv=1,nbas
        ips_iv(i_spackv)=ssite(i_spackv)%spec
      enddo
      
C ... Loop over species containing extended local orbitals
      do  is = 1, nspec
        spid=sspec(is)%name
        z=sspec(is)%z
        lmxa=sspec(is)%lmxa
        lmxb=sspec(is)%lmxb
        if (lmxa .eq. -1) goto 20
        nrspec = iabs ( iclbsj ( is , ips_iv , - nbas , nbas ) )
        if (nrspec .eq. 0) goto 20
        ib = iclbsj ( is , ips_iv , nbas , 1 )
        i_copy_size=size(ssite(ib)%pnu)
        call dcopy(i_copy_size,ssite(ib)%pnu,1,pnui,1)
        i_copy_size=size(ssite(ib)%pz)
        call dcopy(i_copy_size,ssite(ib)%pz,1,pnzi,1)
        if (pnzi(idamax(lmxb+1,pnzi,1),1) .lt. 10) goto 20
C   ... Average over sites within this species
        vseli=0d0 !call dpzero(vseli,4*n0)
        do  ibs = 1, nrspec
          ib = iclbsj ( is , ips_iv , nbas , ibs )
          i_copy_size=size(ssite(ib)%pz)
          call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)
          if (pnzi(idamax(lmxb+1,pnzi,1),1) .lt. 10) goto 22
          vseli = vseli + vsel(:,:,ib)/dble(nrspec)
   22     continue
        enddo
C   ... Printout of input for parameters
        if (ipr .ge. 90/1) then
          write(stdo,261) spid
  261     format(/'  l  site    Eval        Val         Slo         K.E.',
     .    5x,'species ',a)
          do  l = 0, lmxb
            if (pnz(l+1,1) .lt. 10) goto 24
            do  ibs = 1, nrspec
              ib = iclbsj ( is , ips_iv , nbas , ibs )
              write (stdo,260)
     .        l,ib,vsel(4,l+1,ib),(vsel(k,l+1,ib),k=2,4)
  260         format(i3,i4,4f12.6:a)
  262         format(i3,' avg',4f12.6)
            enddo
            write (stdo,262) l,vseli(4,l+1),(vseli(k,l+1),k=2,4)
   24       continue
          enddo
        endif
!!        
        rs3=sspec(is)%rs3
        eh3=sspec(is)%eh3
        vmtz=sspec(is)%vmtz
        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt
        rsmls=0d0 
        ehls=0d0  
        call radmsh(rmt,a,nr,rofi)
!1102 means spin-averaged ehl,rsml
        call loctsh ( 1102 , spid , xx , a , nr , nr , nsp , lmxa , rofi
     .       , wdummy , pnui , pnzi , rs3 , eh3 , vmtz , vseli , rsmls , ehls )
        rsml(1:n0,is) = rsmls(1:n0) !Rsmooth for PZ
        ehl(1:n0,is) = ehls(1:n0)   !Eh for PZ
   20   continue
      enddo
  999 continue
      call tcx('elocp')
      if (allocated(ips_iv)) deallocate(ips_iv)
      end subroutine elocp

      subroutine loctsh(mode,spid,z,a,nr,nrmt,nsp,lmxb,rofi,v,pnu,pnz,
     .rs3,eh3,vmtz,vsel,rsml,ehl)
      use m_lmfinit,only: stdo
C- Fit value and slope of local orbitals to smoothed Hankel
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit concerns fit to low-lying local orbitals
Ci         : 0 do not fit low-lying local orbitals
Ci         : 1 Generate val,slo, and K.E. from potential
Ci         : 2 Fit low-lying local orbitals attempting to match rs,eh
Ci         :   to val, slope and K.E.
Ci         : 3 combination of 1+2
Ci         : 4  -----
Ci         :10s digit concerns fit to high-lying local orbitals
Ci         :   NB: this branch is not implemented
Ci         : 0 do not fit high-lying local orbitals
Ci         : 1 Generate val,slo, and K.E. from potential
Ci         : 2 -----
Ci         : 4 Vary rsm to fit log deriv of high-lying local orbitals
Ci         :   to phi using specified eh3
Ci         : 5 combination of 1+4
Ci         :100s digit concerns constraint on rsm
Ci         : 0 no constraint on maximum rsm
Ci         : 1 constrain rsm to be <= rmt
Ci         :1000s digit deals with spin-polarized case:
Ci         : 0 Treat spins separately
Ci         : 1 Compute ehl,rsml for average potential
Ci         :10000s digit 1 plot wave functions (not implemented)
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(ir) = b [e^(a(ir-1)) -1]
Ci   nr    :number of radial mesh points in potential sphere
Ci   nrmt  :number of points between 0..rmt
Ci   lmxb  :l-cutoff for basis
Ci   pl    :boundary conditions for valence wavefunctions.
Ci   pnz   :boundary conditions for local orbital. pnz=0 -> no loc. orb.
Ci         :10s digit controls how local orbital included in hamiltonian
Ci         :10s digit nonzero -> smooth Hankel tail is attached.
Ci   rs3   :minimum allowed smoothing radius in attaching Hankel tails
Ci         :to local orbitals
Ci   eh3   :Hankel energy when attaching Hankel tails to high-lying
Ci         :local orbitals
Ci   vmtz  :parameter used in attaching Hankel tails to local orbitals
Ci         :It is used as a constant shift to Hankel energies for the
Ci         :fitting of local orbitals to Hankel tails. Thus vmtz
Ci         :is an estimate for the potential at the MT radius.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   v     :spherical potential (atomsr.f)
Ci   rofi  :radial mesh points
Ci   spid  :species label
Cio Inputs/Outputs
Cio  vsel  :value, slope, K.E. of radial w.f.
Cio        :output if 1s digit (or 10s digit) has 1's bit set
Cio        :input  if 1s digit (or 10s digit) has 1's bit clear
Co Outputs
Co   ehl   :energies of smoothed Hankel tail for local orbital
Co   rsml  :smoothing radii for smoothed Hankel tail of local orbital
Cl Local variables
Cl   rmt   :muffin-tin radius, in a.u.
Cr Remarks
Cu Updates
Cu   06 Jul 05 generation of (val,slo,k.e.) and fitting split into
Cu             independent operations
Cu   27 Jun 04 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
c      use m_globalvariables
      integer lmxb,nr,nrmt,nsp,n0,mode
      double precision a,z,rs3,eh3,vmtz
      parameter (n0=10)
      double precision rofi(1),v(nr,nsp),pnz(n0,nsp),pnu(n0,nsp)
      double precision rsml(n0,2),ehl(n0,2),vsel(4,n0)
C ... Local parameters
      character spid*8,orbit(2)*4,flg(2)*1
      integer:: ipr , iprint , i , l , info , mode0 
     ., mode1 , mode2 , mode3 , mode4 , modei , loclo , nfit , iprt 
     ., isw
      real(8) ,allocatable :: g_rv(:)
      real(8) ,allocatable :: gp_rv(:)
      real(8) ,allocatable :: h_rv(:)
      real(8):: vv(nr,nsp)
      integer nrx,nrbig,ir,lp1
      double precision dphi,dphip,eh,eval,p,phi,phip,
     .pnul,rsm,rsmin,rsmax,rmt,ekin,h,sum1,sum2,qrmt
C     emin and emax are the maximum allowed ranges in Hankel energies
C     for the fitting of local orbitals to Hankel tails.
      double precision emin,emax,vmine
      parameter (nrx=1501)
      parameter (emin=-10d0,emax=-.02d0,IPRT=20)
      double precision rbig,rofix(nrx),rwgtx(nrx),xi(0:20)
      data orbit/'low','high'/
      data flg/'*',' '/
C --- Setup ---
      ipr = iprint()
      if (lmxb .gt. 8) call rx('loctsh:  lmax too large')
      rmt = rofi(nrmt)
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)
      mode4 = mod(mode/10000,10)
      nfit = 0
      do  i = 1, nsp
        allocate(h_rv(nr))
        allocate(g_rv(2*nr))
        allocate(gp_rv(2*nr*4))
C --- Loop over local orbitals ---
        do  l = 0, lmxb
C       itab(l+1,i) = 0
          rsml(l+1,i) = 0
          ehl(l+1,i) = 0
C       pl(l+1,i) = pnu(l+1,i)
C       konfig = mod(pnz(l+1,i),10d0)
C       Skip all but local orbitals with tails attached
          if (pnz(l+1,i) .lt. 10) goto 10
C       Case local orbital deeper than valence
          if (int(pnu(l+1,i)-1) .eq. int(mod(pnz(l+1,i),10d0))) then
            loclo = 1
            if (mode0 .eq. 0) goto 10
            modei = mode0
C       Case local orbital higher than the valence state
          elseif (int(pnu(l+1,i)+1) .eq. int(mod(pnz(l+1,i),10d0))) then
            loclo = 0
            if (mode1 .eq. 0) goto 10
            modei = mode1
C       Local orbital neither low nor high: error
          else
            call fexit3(-1,111,' Exit -1 loctsh, l=%i:  sc PZ=%d '//
     .      'incompatible with valence P=%;3d',l,pnz(l+1,i),pnu(l+1,i))
          endif
          pnul = mod(pnz(l+1,i),10d0)
C   ... Initial printout
          nfit = nfit + 1
          if (nfit .eq. 1 .and. modei .gt. 1) then
            call info2(20,1,0,
     .      ' Fit local orbitals to sm hankels, species '//spid//
     .      '%a, rmt=%;7g',rmt,0)
            if (ipr .ge. IPRT) write (stdo,261)
          endif
C   ... Make value, slope, kinetic energy
          if (mod(modei,2) .eq. 1) then
C     ... Overwrite V+, V- with (V+ + V-)/2, (V+ - V-)/2,
            if (nsp .eq. 2 .and. mode3 .eq. 1) then
c               call dsumdf(nr,0.5d0,v,0,1,v(1,2),0,1)
               vv(:,1)= .5d0*(v(:,1)+v(:,2))
               vv(:,2)= .5d0*(v(:,1)-v(:,2))
            else
               vv= v
            endif
C     ... Wave function and potential parameters at MT sphere
            call makrwf ( 0 , z , rmt , l , vv( 1 , i ) , a , nrmt , rofi 
     .      , pnz ( 1 , i ) , 4 , g_rv , gp_rv , eval , phi , dphi 
     .      , phip , dphip , p )
cC     ... Restore V+, V-
c            if (nsp .eq. 2 .and. mode3 .eq. 1) then
c              call dsumdf(nr,1d0,v,0,1,v(1,2),0,1)
c            endif
            vsel(1,l+1) = phi
            vsel(2,l+1) = dphi
            vmine = eval - (v(nrmt,i)-2*z/rmt)
            vsel(3,l+1) = vmine
            vsel(4,l+1) = eval
          else
            phi   = vsel(1,l+1)
            dphi  = vsel(2,l+1)
            vmine = vsel(3,l+1)
            eval  = vsel(4,l+1)
          endif
C   ... Set conditions on envelope functions
          if (modei .gt. 1) then
            rsmin = rs3
            rsmax = 5
            if (mode2 .eq. 1) rsmax = rmt
            if (eval .lt. emin)
     .      call rx1('increase emin in loctsh: eval=%;4g',eval)
C   ... Match val,slo, and K.E. of Hankel to phi,dphi,vmine
            if (modei .eq. 2 .or. modei .eq. 3) then
              rsm = 0
              eh = min(eval-vmtz,emax)
              call mtchre(103,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,
     .        vmine,dphi,rsm,eh,ekin,info)
C   ... Vary rsm to match sm Hankel to phi,dphi
            elseif (modei .eq. 4 .or. modei .eq. 5) then
              call rx('this branch not checked')
              rsm = rsmin
              call mtchre(100,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,phi,
     .        dphi,rsm,eh,ekin,info)
            else
              call rxi('loctsh: not implemented fitting mode=',modei)
            endif
C   ... Printout of fit functions
            if (ipr .ge. IPRT) then
              call radext(11,nr,nrx,2d0,a,rmt,nrbig,rbig,rofix,rwgtx)
C         make sphere charge for r>rmt for printout
              lp1 = l+1
              sum1 = 0
              do  ir = 1, nr
C           h = r*radial part of sm. Hankel
                call hansmr(rofix(ir),eh,1/rsm,xi,l)
                h = xi(l)*(rofix(ir)**lp1)
                sum1 = sum1 + rwgtx(ir)*h**2
              enddo
              sum2 = 0
              do  ir = nr, nrbig
C           h = r*radial part of sm. Hankel
                call hansmr(rofix(ir),eh,1/rsm,xi,l)
                h = xi(l)*(rofix(ir)**lp1)
                sum2 = sum2 + rwgtx(ir)*h**2
              enddo
              qrmt = sum2/(sum1+sum2)

              write (stdo,260)
     .        l,orbit(2-loclo),pnul,eval,vmine,rsm,eh,qrmt,ekin,
     .        flg(2-isw(dabs(ekin-vmine).gt.1d-5))
  261         format('  l  type    Pnu      Eval        K.E.',
     .        7x,'Rsm       Eh      Q(r>rmt)    Fit K.E.')
  260         format(i3,2x,a,f8.3,2f12.6,3f10.5,f12.6,a1)
            endif
            rsml(l+1,i) = rsm
            ehl(l+1,i) = eh
          endif
   10     continue
        enddo
        if (mode3 .ne. 0) goto 80
      enddo
   80 continue
      if (mode4 .eq. 1) call rx('loctsh: plotting not ready')
      if (allocated(gp_rv)) deallocate(gp_rv)
      if (allocated(g_rv)) deallocate(g_rv)
      if (allocated(h_rv)) deallocate(h_rv)
      end
     
      end module m_elocp
