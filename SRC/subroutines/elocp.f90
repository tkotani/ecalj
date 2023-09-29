module m_elocp ! envlope parameters for extended local orbitals
  use m_lmfinit,only: nspec,n0
  real(8),allocatable,protected,public:: ehl(:,:),rsml(:,:)
  public:: elocp
  private
contains
  subroutine elocp()! Make envlope parameters for extended local orbitals
    use m_lmfinit,only: stdo,nspec,nbas,nsp,ispec,sspec=>v_sspec,n0,nkapii,slabl,vmtz,rs3,eh3
    use m_lmfinit,only: z_i=>z,nr_i=>nr,lmxa_i=>lmxa,rmt_i=>rmt,lmxb_i=>lmxb, spec_a
    use m_density,only: v0pot,pnuall,pnzall
    !o Outputs
    !o   ehl,rsml
    !o         : smoothing radius and energy set for extended local orbitals
    implicit none
    character spid*8
    logical :: eloc
    integer :: ib,ibs,iclbsj,ipr,iprint,is,k,l,lmxa,lmxb,nglob,nkap0,nr,nrmx,nrspec
    parameter (nrmx=1501, nkap0=3)
    integer:: nkape,idamax,nkapex
    integer,allocatable :: ips_iv(:)
    real(8):: z,a,rmt,xx, rofi(nrmx),vseli(4,n0),vsel(4,n0,nbas), eh(n0,nkap0),rsmh(n0,nkap0) 
    real(8),pointer:: pnu(:,:),pnz(:,:),pnui(:,:),pnzi(:,:)
    real(8):: ehls(n0*2),rsmls(n0*2),wdummy(1)
    call tcn('elocp')
    if(allocated(ehl)) deallocate(ehl,rsml)
    allocate( ehl(n0,nspec),rsml(n0,nspec),source=0d0)
    ipr = iprint()
    eloc = .false.
    vsel=0d0
    ! --- Find val, slo, K.E. for all sites ---
    do  ib = 1, nbas
       is=ispec(ib) 
       pnu=>pnuall(:,:,ib)
       pnz=>pnzall(:,:,ib)
       spid = slabl(is) 
       a=spec_a(is)
       nr=nr_i(is)
       rmt=rmt_i(is)
       z=z_i(is)
       lmxa=lmxa_i(is)
       lmxb=lmxb_i(is)
       if (lmxa == -1) cycle
       if (pnz(idamax(lmxb+1,pnz,1),1) < 10) cycle
       eloc = .true.
       call radmsh(rmt,a,nr,rofi)
       call loctsh ( 1,spid,z,a,nr,nr,nsp,lmxa,rofi & !mode1101
           ,v0pot(ib)%v,pnu,pnz,xx,xx,vmtz(is),vsel ( 1,1,ib ) ,rsml,ehl )
    enddo
    if ( .NOT. eloc) goto 999
    if (ipr >= 30) write(stdo,"(/' elocp:')")
    ! --- Determine shape of smooth Hankel tails for local orbitals ---
    allocate(ips_iv(nbas))
    ips_iv=ispec
    ! ... Loop over species containing extended local orbitals
    do  is = 1, nspec
       spid=slabl(is)
       z= z_i(is)
       lmxa=lmxa_i(is)
       lmxb=lmxb_i(is)
       if (lmxa == -1) cycle
       nrspec = iabs ( iclbsj ( is,ips_iv,- nbas,nbas ) )
       if (nrspec == 0) cycle
       ib = iclbsj ( is,ips_iv,nbas,1 )
       pnui=>pnuall(:,:,ib)
       pnzi=>pnzall(:,:,ib)
       if (pnzi(idamax(lmxb+1,pnzi,1),1) < 10) cycle
       vseli=0d0 !Average over sites within this species
       do  ibs = 1, nrspec
          ib = iclbsj ( is,ips_iv,nbas,ibs )
          if (pnzi(idamax(lmxb+1,pnzi,1),1) < 10) cycle
          vseli = vseli + vsel(:,:,ib)/dble(nrspec)
       enddo
       a=  spec_a(is)
       nr= nr_i(is)
       rmt=rmt_i(is)
       call radmsh(rmt,a,nr,rofi)
       call loctsh ( 2,spid,xx,a,nr,nr,nsp,lmxa,rofi & !mode1102
           ,wdummy,pnui,pnzi,rs3(is),eh3(is),vmtz(is),vseli,rsmls,ehls )
       rsml(1:n0,is) = rsmls(1:n0) !Rsmooth for PZ
       ehl(1:n0,is) = ehls(1:n0)   !Eh for PZ
    enddo
999 continue
    call tcx('elocp')
    if (allocated(ips_iv)) deallocate(ips_iv)
  end subroutine elocp
  subroutine loctsh(mode0,spid,z,a,nr,nrmt,nsp,lmxb,rofi,v,pnu,pnz,rs3,eh3,vmtz, vsel,rsml,ehl)!Fit value and slope of local orbitals to smoothed Hankel
    use m_mtchae,only:mtchre
    use m_lmfinit,only: stdo
    use m_hansr,only:hansmr
    use m_atwf,only: makrwf
    use m_ftox
    !i Inputs mode=110x allowed now. Only x is supplied
    !i   mode  : For fitting to low-lying local orbitals,
    !i         := 1 Generate val,slo, and K.E. from potential
    !i         := 2 Fit low-lying local orbitals attempting to match rs,eh to val, slope and K.E.
    !i      We do not fit high-lying local orbitals
    !i      We have constrain rsm to be <= rmt
    !i      We Compute ehl,rsml for average potential
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
    integer:: ipr,iprint,i,l,info,mode0, modei,loclo=-99999,nfit,iprt,isw
    !, mode1,mode2,mode3,mode4
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
       if (pnz(l+1,i) < 10) goto 10!
       if (int(pnu(l+1,i)-1) == int(mod(pnz(l+1,i),10d0))) then     ! Case local orbital deeper than valence
          loclo = 1
          if (mode0 == 0) goto 10
          modei = mode0
       elseif (int(pnu(l+1,i)+1) == int(mod(pnz(l+1,i),10d0))) then !  Case local orbital higher than the valence state
          loclo = 0
          goto 10 
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
