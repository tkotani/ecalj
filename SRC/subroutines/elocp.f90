module m_elocp ! envlope parameters for extended local orbitals
  use m_lmfinit,only: nspec,n0
  real(8),allocatable,protected,public:: ehl(:,:),rsml(:,:)
  public:: elocp
  private
contains
  subroutine elocp()! Make envlope parameters for extended local orbitals
    use m_lmfinit,only: stdo,nspec,nbas,nsp,ispec,sspec=>v_sspec,n0,nkapii,slabl,vmtz,rs3,eh3,lpzex
    use m_density,only: v0pot,pnuall,pnzall
    ! ----------------------------------------------------------------------
    !i Inputs
    !i         : 0 do nothing ; just return
    !i         : 1 make core and augmentation matrices
    !o Outputs
    !o   ehl,rsml
    !o         : smoothing radius and energy set for extended local orbitals
    !r Remarks
    !u Updates
    !u   06 Jul 05 first created
    ! ----------------------------------------------------------------------
    implicit none
    character spid*8
    logical :: eloc
    integer :: ib,ibs,iclbsj,ipr,iprint,is,k,l,lmxa,lmxb,nglob,nkap0,nr,nrmx,nrspec
    parameter (nrmx=1501, nkap0=3)
    integer:: nkape,idamax,nkapex
    integer,allocatable :: ips_iv(:)
    real(8):: z,a,rmt,xx, rofi(nrmx),vseli(4,n0),vsel(4,n0,nbas), &
         eh(n0,nkap0),rsmh(n0,nkap0) !
    real(8),pointer:: pnu(:,:),pnz(:,:),pnui(:,:),pnzi(:,:)
    real(8):: ehls(n0*2),rsmls(n0*2),wdummy(1)
    ! --- Setup ---
    call tcn('elocp')
    if(allocated(ehl)) deallocate(ehl,rsml)
    allocate( ehl(n0,nspec),rsml(n0,nspec),source=0d0)
!    if(sum(lpzex)/=0) return
    ipr = iprint()
    eloc = .false.
    vsel=0d0
    ! --- Find val, slo, K.E. for all sites ---
    do  ib = 1, nbas
       is=ispec(ib) 
       pnu=>pnuall(:,:,ib)
       pnz=>pnzall(:,:,ib)
       spid = slabl(is) !sspec(is)%name
       a=sspec(is)%a
       nr=sspec(is)%nr
       rmt=sspec(is)%rmt
       z=sspec(is)%z
       lmxa=sspec(is)%lmxa
       lmxb=sspec(is)%lmxb
       if (lmxa == -1) goto 10
       if (pnz(idamax(lmxb+1,pnz,1),1) < 10) goto 10
       eloc = .true.
       call radmsh(rmt,a,nr,rofi)
       call loctsh ( 1 , spid , z , a , nr , nr , nsp , lmxa , rofi & !1101
            , v0pot(ib)%v , pnu , pnz , xx , xx , vmtz(is) , vsel ( 1 , 1 , ib ) &
            , rsml , ehl )
10     continue
    enddo
    if ( .NOT. eloc) goto 999
    if (ipr >= 30) write(stdo,199)
199 format(/' elocp:')
    ! --- Determine shape of smooth Hankel tails for local orbitals ---
    allocate(ips_iv(nbas))
    ips_iv=ispec
    ! ... Loop over species containing extended local orbitals
    do  is = 1, nspec
       spid=slabl(is) !sspec(is)%name
       z=sspec(is)%z
       lmxa=sspec(is)%lmxa
       lmxb=sspec(is)%lmxb
       if (lmxa == -1) goto 20
       nrspec = iabs ( iclbsj ( is , ips_iv , - nbas , nbas ) )
       if (nrspec == 0) goto 20
       ib = iclbsj ( is , ips_iv , nbas , 1 )
       pnui=>pnuall(:,:,ib)
       pnzi=>pnzall(:,:,ib)
       
       if (pnzi(idamax(lmxb+1,pnzi,1),1) < 10) goto 20
       !   ... Average over sites within this species
       vseli=0d0
       do  ibs = 1, nrspec
          ib = iclbsj ( is , ips_iv , nbas , ibs )
          if (pnzi(idamax(lmxb+1,pnzi,1),1) < 10) goto 22
          vseli = vseli + vsel(:,:,ib)/dble(nrspec)
22        continue
       enddo
!        !   ... Printout of input for parameters
!        if (ipr >= 90/1) then
!           write(stdo,261) spid
! 261       format(/'  l  site    Eval        Val         Slo         K.E.', &
!                5x,'species ',a)
!           do  l = 0, lmxb
!              if (pnz(l+1,1) < 10) goto 24
!              do  ibs = 1, nrspec
!                 ib = iclbsj ( is , ips_iv , nbas , ibs )
!                 write (stdo,260) l,ib,vsel(4,l+1,ib),(vsel(k,l+1,ib),k=2,4)
! 260             format(i3,i4,4f12.6:a)
! 262             format(i3,' avg',4f12.6)
!              enddo
!              write (stdo,262) l,vseli(4,l+1),(vseli(k,l+1),k=2,4)
! 24           continue
!           enddo
!        endif
       !rs3= sspec(is)%rs3
       !eh3= sspec(is)%eh3
       !vmtz=sspec(is)%vmtz
       a=sspec(is)%a
       nr=sspec(is)%nr
       rmt=sspec(is)%rmt
       rsmls=0d0
       ehls=0d0
       call radmsh(rmt,a,nr,rofi)
       ! 102 means spin-averaged ehl,rsml
       call loctsh ( 2 , spid , xx , a , nr , nr , nsp , lmxa , rofi & !1102
            , wdummy , pnui , pnzi , rs3(is) , eh3(is) , vmtz(is) , vseli , rsmls , ehls )
       rsml(1:n0,is) = rsmls(1:n0) !Rsmooth for PZ
       ehl(1:n0,is) = ehls(1:n0)   !Eh for PZ
20     continue
    enddo
999 continue
    call tcx('elocp')
    if (allocated(ips_iv)) deallocate(ips_iv)
  end subroutine elocp
  subroutine loctsh(mode0,spid,z,a,nr,nrmt,nsp,lmxb,rofi,v,pnu,pnz,rs3,eh3,vmtz, vsel,rsml,ehl)
    use m_mtchae,only:mtchre
    use m_lmfinit,only: stdo
    use m_hansr,only:hansmr
    use m_atwf,only: makrwf
    use m_ftox
    !- Fit value and slope of local orbitals to smoothed Hankel
    ! ----------------------------------------------------------------------
    !i Inputs mode=110x allowed now
    !i   mode  :1s digit concerns fit to low-lying local orbitals
    !i         : 1 Generate val,slo, and K.E. from potential
    !i         : 2 Fit low-lying local orbitals attempting to match rs,eh
    !i         :   to val, slope and K.E.
    !i         :xxx10s digit concerns fit to high-lying local orbitals
    !i         :   NB: this branch is not implemented
    !i         : 0 do not fit high-lying local orbitals
    !i         :xxx100s digit concerns constraint on rsm
    !i         : 1 constrain rsm to be <= rmt
    !i         :xxx1000s digit deals with spin-polarized case:
    !i         : 1 Compute ehl,rsml for average potential
    !i   z     :nuclear charge
    !i   a     :the mesh points are given by rofi(ir) = b [e^(a(ir-1)) -1]
    !i   nr    :number of radial mesh points in potential sphere
    !i   nrmt  :number of points between 0..rmt
    !i   lmxb  :l-cutoff for basis
    !i   pl    :boundary conditions for valence wavefunctions.
    !i   pnz   :boundary conditions for local orbital. pnz=0 -> no loc. orb.
    !i         :10s digit controls how local orbital included in hamiltonian
    !i         :10s digit nonzero -> smooth Hankel tail is attached.
    !i   rs3   :minimum allowed smoothing radius in attaching Hankel tails
    !i         :to local orbitals
    !i   eh3   :Hankel energy when attaching Hankel tails to high-lying
    !i         :local orbitals
    !i   vmtz  :parameter used in attaching Hankel tails to local orbitals
    !i         :It is used as a constant shift to Hankel energies for the
    !i         :fitting of local orbitals to Hankel tails. Thus vmtz
    !i         :is an estimate for the potential at the MT radius.
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   v     :spherical potential (atomsr.f)
    !i   rofi  :radial mesh points
    !i   spid  :species label
    ! o Inputs/Outputs
    ! o  vsel  :value, slope, K.E. of radial w.f.
    ! o        :output if 1s digit (or 10s digit) has 1's bit set
    ! o        :input  if 1s digit (or 10s digit) has 1's bit clear
    !o Outputs
    !o   ehl   :energies of smoothed Hankel tail for local orbital
    !o   rsml  :smoothing radii for smoothed Hankel tail of local orbital
    !l Local variables
    !l   rmt   :muffin-tin radius, in a.u.
    !r Remarks
    !u Updates
    !u   06 Jul 05 generation of (val,slo,k.e.) and fitting split into
    !u             independent operations
    !u   27 Jun 04 First created
    ! ----------------------------------------------------------------------
    implicit none
    integer :: lmxb,nr,nrmt,nsp,n0,mode
    double precision :: a,z,rs3,eh3,vmtz
    parameter (n0=10)
    double precision :: rofi(1),v(nr,nsp),pnz(n0,nsp),pnu(n0,nsp)
    double precision :: rsml(n0,2),ehl(n0,2),vsel(4,n0)
    character spid*8,orbit(2)*4,flg(2)*1
    integer:: ipr , iprint , i , l , info , mode0, modei , loclo=-99999 , nfit , iprt , isw
    !, mode1 , mode2 , mode3 , mode4
    real(8) ,allocatable :: g_rv(:)
    real(8) ,allocatable :: gp_rv(:)
    real(8) ,allocatable :: h_rv(:)
    real(8):: vv(nr,nsp)
    integer :: nrx,nrbig,ir,lp1
    double precision :: dphi,dphip,eh,eval,p,phi,phip,pnul,rsm,rsmin,rsmax,rmt,ekin,h,sum1,sum2,qrmt
    !     emin and emax are the maximum allowed ranges in Hankel energies
    !     for the fitting of local orbitals to Hankel tails.
    double precision :: emin,emax,vmine
    parameter (nrx=1501)
    parameter (emin=-10d0,emax=-.02d0,IPRT=20)
    double precision :: rbig,rofix(nrx),rwgtx(nrx),xi(0:20)
    character*120:: aaa
    data orbit/'low','high'/
    data flg/'*',' '/
    ipr = iprint()
    if (lmxb > 8) call rx('loctsh:  lmax too large')
    rmt = rofi(nrmt)
    nfit = 0
    allocate(h_rv(nr))
    allocate(g_rv(2*nr))
    allocate(gp_rv(2*nr*4))
    ! --- Loop over local orbitals ---
    i=1 !isp=1 for averaged spin mode
    do  l = 0, lmxb
       rsml(l+1,i) = 0
       ehl(l+1,i) = 0
       !       pl(l+1,i) = pnu(l+1,i)
       !       konfig = mod(pnz(l+1,i),10d0)
       !       Skip all but local orbitals with tails attached
       if (pnz(l+1,i) < 10) goto 10!
       !       Case local orbital deeper than valence
       if (int(pnu(l+1,i)-1) == int(mod(pnz(l+1,i),10d0))) then
          loclo = 1
          if (mode0 == 0) goto 10
          modei = mode0
          !       Case local orbital higher than the valence state
       elseif (int(pnu(l+1,i)+1) == int(mod(pnz(l+1,i),10d0))) then
          loclo = 0
          goto 10 !if (mode1 == 0) goto 10
          modei = 0 !mode1
       else!   Local orbital neither low nor high: error
          write(aaa,ftox)'Exit -1 loctsh',l,ftof(pnz(l+1,i),3),'incompatible with valence P=',ftof(pnu(l+1,i),3)
          call rx(trim(aaa))
       endif
       pnul = mod(pnz(l+1,i),10d0)
       nfit = nfit + 1
       if (nfit == 1 .AND. modei > 1) then
          write(stdo,ftox)'Fit local orbitals to sm hankels, species '//trim(spid),ftof(rmt)
          if(ipr >= IPRT) write (stdo,261)
       endif
       !   ... Make value, slope, kinetic energy
       if (mod(modei,2) == 1) then
          if (nsp == 2) then ! .AND. mode3 == 1) then
             vv(:,1)= .5d0*(v(:,1)+v(:,2)) !spin averaged
             vv(:,2)= .5d0*(v(:,1)-v(:,2)) !up - dn
          else
             vv= v
          endif
          !     ... Wave function and potential parameters at MT sphere
          call makrwf(z,rmt,l,vv(1,i),a,nrmt,rofi,pnz(1,i),4,g_rv,gp_rv,eval,phi,dphi,phip,dphip,p)
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
       !   ... Set conditions on envelope functions
       if (modei > 1) then
          rsmin = rs3
          rsmax = 5
          rsmax = rmt !if (mode2 == 1) 
          if (eval < emin) call rx1('increase emin in loctsh: eval=%;4g',eval)
          !   ... Match val,slo, and K.E. of Hankel to phi,dphi,vmine
          if (modei == 2 .OR. modei == 3) then
             rsm = 0
             eh = min(eval-vmtz,emax)
             call mtchre(103,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,vmine,dphi,rsm,eh,ekin,info)
             !   ... Vary rsm to match sm Hankel to phi,dphi
          elseif (modei == 4 .OR. modei == 5) then
             call rx('this branch not checked')
             rsm = rsmin
             call mtchre(100,l,rsmin,rsmax,emin,emax,rmt,rmt,phi,dphi,phi,dphi,rsm,eh,ekin,info)
          else
             call rxi('loctsh: not implemented fitting mode=',modei)
          endif
          if (ipr >= IPRT) then!   ... Printout of fit functions
             call radext(11,nr,nrx,2d0,a,rmt,nrbig,rbig,rofix,rwgtx)
             !         make sphere charge for r>rmt for printout
             lp1 = l+1
             sum1 = 0
             do  ir = 1, nr
                call hansmr(rofix(ir),eh,1/rsm,xi,l)!h = r*radial part of sm. Hankel
                h = xi(l)*(rofix(ir)**lp1)
                sum1 = sum1 + rwgtx(ir)*h**2
             enddo
             sum2 = 0
             do  ir = nr, nrbig
                call hansmr(rofix(ir),eh,1/rsm,xi,l) !h = r*radial part of sm. Hankel
                h = xi(l)*(rofix(ir)**lp1)
                sum2 = sum2 + rwgtx(ir)*h**2
             enddo
             qrmt = sum2/(sum1+sum2)
             write(stdo,260) l,orbit(2-loclo),pnul,eval,vmine,rsm,eh,qrmt,ekin,flg(2-isw(dabs(ekin-vmine).gt.1d-5))
261          format('  l  type    Pnu      Eval        K.E.',7x,'Rsm       Eh      Q(r>rmt)    Fit K.E.')
260          format(i3,2x,a,f8.3,2f12.6,3f10.5,f12.6,a1)
          endif
          rsml(l+1,i) = rsm
          ehl(l+1,i) = eh
       endif
10     continue
    enddo
    if (allocated(gp_rv)) deallocate(gp_rv)
    if (allocated(g_rv)) deallocate(g_rv)
    if (allocated(h_rv)) deallocate(h_rv)
  end subroutine loctsh
end module m_elocp
