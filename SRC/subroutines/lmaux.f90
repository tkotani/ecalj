module m_lmaux
  use m_lgunit,only:stdo
  public:: lmaux
  private
contains

  subroutine lmaux()        !main part of lmchk
    use m_mksym,only: ctrl_nclass,iv_a_oics,iclasst
    use m_lmfinit,only: iv_a_oips,str_mxnbr,str_rmax,ctrl_nbas,ctrl_nspec,ctrl_nspin, &
         ctrl_nl,ctrl_omax1,ctrl_omax2,ctrl_wsrmax,slabl,sspec=>v_sspec, &
         lat_avw,lat_alat,mxcst2
    use m_lattic,only: lat_nkd
    use m_lattic,only: lat_nkq
    use m_struc_def
    !      use m_ovmin , only: ovmin
    use m_lattic,only:lat_plat,rv_a_opos
    !! check crystal structure symmetry and get WSR
    ! ----------------------------------------------------------------------
    !i Inputs
    !u Updates
    !u   12 Aug 08 (L. Ke) empty sphere finder
    ! xxx   04 Nov 04 Upgrade of rsta editor
    ! xxx   26 Jan 03 Call to angtab changed
    !u   17 May 02 Modified MT radii scaling to lower priority for E.S.
    !u   23 Apr 02 Added option (--getwsr) to find MT radii
    !u   01 Mar 02 Updated Import data mode
    ! xxx   05 Oct 01 Adapted mode 2**3 to work with lm v6.11
    !u   24 Nov 97 changed ovmin to run quickly
    ! ----------------------------------------------------------------------
    implicit none
    integer:: mode=1 !,wksize
    character(120) :: outs,fnam(8)
    integer :: NULLI
    logical :: cmdopt,T,F,swtmp
    parameter (T=.true., F=.false., NULLI=-99999)
    integer :: getdig,i,ip,j,k,m,ifi,iprint,lpbc, &
         nbas,nclasp,nclass,nl,nlspc,nsp,modep(3),parg,nbasp, &
         nbaspp,nkd,nkq,nspec,neul,nc,mxcsiz,nttab,igets, & ! & npadl,npadr,
         iosits,cmplat,ngrp,ival,irs(5),nclspp,bitand,igetss, &
         ngmx,nsgrp
    integer:: oeold  , olmx , opold , owk2 , orham , oamsh &
         , onrmsh , oalpha , onpr , os , ormx , oip , opgfsl , mxclas
    integer,allocatable :: iv_a_ontab(:)
    integer,allocatable :: iv_a_oiax(:,:)
    real(8),allocatable :: rv_a_og(:),rv_a_ormax(:)
    real(8) ,allocatable :: pos2_rv(:,:)
    real(8) ,allocatable :: rmt_rv(:)
    integer ,allocatable :: lock_iv(:)
    real(8) ,allocatable :: lockc_rv(:)
    real(8) ,allocatable :: z_rv(:)
    real(8) ,allocatable :: zz_rv(:)
    integer ,allocatable :: ips2_iv(:)
    real(8) ,allocatable :: zc_rv(:)
    real(8) ,allocatable :: rmtc_rv(:)
    double precision :: xv(10),xx,alat,plat(3,3),facrmx,facrng, & ! & ,plat2(9)
         dval,avw,ekap(2),enu,qss(4),ckbas,cksumf,ehterm(4), rmaxs, &
         qlat(9),emad,trumad,vmtz(2),omax1(3),omax2(3),wsrmax
    parameter (ngmx=48,mxclas=1000)
    integer:: i_copy_size, i_spackv, i_spacks
    integer:: ifx,w_dummy(1)=1
    integer,allocatable:: lmxa(:)
    real(8),allocatable:: z(:),rmax(:)
    print *,' lmaux:'
    nbas=ctrl_nbas
    nclass=ctrl_nclass
    nl=ctrl_nl
    nspec=ctrl_nspec
    nsp=ctrl_nspin
    !      modep = ctrl_modep
    modep=99999 !bug fix 2022-6-29 (no initialization before. No problem as long as lmchk works.)
    lpbc = 0
    nclasp=ctrl_nclass !sarray%nclasp
    avw=lat_avw
    alat=lat_alat
    plat=lat_plat
    nkd=lat_nkd
    nkq=lat_nkq
    rmaxs=str_rmax !sstr%rmax
    nclspp = max(nclass,nspec)
    allocate(rv_a_ormax(nclspp))
    rv_a_ormax = sspec(iv_a_oics(1:nclasp))%rmt
    allocate(lmxa(nclasp),z(nclasp))
    lmxa(1:nclasp) = sspec(iv_a_oics(1:nclasp))%lmxa !sarray%
    z   (1:nclasp) = sspec(iv_a_oics(1:nclasp))%z
    !print *,' nclasp,lmxa=',nclasp,lmxa
    !print *,' z   =',z
    nbasp = nbas !+ npadl + npadr
    nbaspp = nbas !2*nbasp - nbas
    !      stdo = lgunit(1)
    j = 10
    if (cmdopt('--shorten',j-1,0,outs)) then
       call shorps ( nbasp , plat , modep , rv_a_opos , rv_a_opos )
    endif
    ! --- Neighbor tables and sphere overlaps ---
    !      if (getdig(mode,0,2) .ne. 0) then
    if (rmaxs <= 0d0) then
       rmaxs = 2.7d0*avw
       call info5(30,0,0,'%1f'// &
            'Use default rmaxs = %;3d a.u. = %;3d*avw = %;3d*alat', &
            rmaxs,rmaxs/avw,rmaxs/alat,0,0)
    endif
    ! ... Get neighbor table iax for each atom in the cluster
!    if (lpbc == 0) then
       i = 3 !3dim
       j = -1
!    elseif (lpbc == 1 .OR. lpbc == 11) then
!       i = 2
!       j = 1
!    else
!       call rx('ASASTR: not implemented for lpbc>1')
!    endif
    mxcsiz = str_mxnbr !int(sstr%mxnbr)

    call pshpr(iprint()-20)
    call pairs ( nbas , nbasp , alat , plat ,(/ rmaxs / 2/) , rv_a_opos &
         , (/- 1/), w_dummy , nttab , iv_a_ontab , iv_a_oiax , mxcsiz ) !, i,j
    call poppr

    ! --- Print out a few superlattice vectors ---
    j = 6
    if (cmdopt('--slat',j-1,0,outs)) then
       if (iprint() >= 10) then
          call info0(10,1,0,' LMCHK:  print multiples of plat%N'// &
               '  i1  i2  i3%7fx%11fy%11fz%11flen')
          do  i = -2, 2
             do  j = -2, 2
                do  k = -2, 2
                   xx = 0
                   do  m = 1, 3
                      xv(m) = i*plat(m,1) + j*plat(m,2) + k*plat(m,3)
                      xx = xx + xv(m)**2
                   enddo
                   xx = dsqrt(xx)
                   print 368, i,j,k, xv(1), xv(2), xv(3), xx
368                format(3i4, 3f12.7, 1x, f12.5)
                enddo
             enddo
          enddo
       endif
    endif

    ! --- Find sphere overlaps ---
    j = 9
    ifx=0
    if (cmdopt('--getwsr',j-1,0,outs)) then
       call info(10,1,0,' ... Make sphere radii',0,0)
       allocate(zz_rv(nspec))
       allocate(rmt_rv(nspec))
       do i_spackv=1,nspec
          zz_rv (i_spackv) = sspec(i_spackv)%z
          rmt_rv(i_spackv) = sspec(i_spackv)%rmt
       enddo
       allocate(lock_iv(nspec))
       lock_iv(:)=0
       do  i = 1, nspec
          lock_iv(i)= mxcst2(i) !bitand ( int ( sspec ( i )%mxcst ) , 2 ) )
       enddo
       if (lpbc == 0) then
          i = 3
       elseif (lpbc == 1 .OR. lpbc == 11) then
          i = 2
       else
          call rx('LMAUX: not implemented for lpbc>1')
       endif
       call makrm0 ( 101 , nspec , nbas , alat , plat , rv_a_opos , &
            slabl , iv_a_oips , modep , lock_iv , zz_rv , rmt_rv ) !sarray%
       !   ... Scale sphere radii satisfying constraints
       i_copy_size=size(ctrl_omax1)
       call dcopy(i_copy_size,ctrl_omax1,1,omax1,1)
       i_copy_size=size(ctrl_omax2)
       call dcopy(i_copy_size,ctrl_omax2,1,omax2,1)
       wsrmax=ctrl_wsrmax
       call sclwsr ( 20 , nbas , nbasp , nspec , alat , plat , rv_a_opos &
            , iv_a_oips , modep , slabl , zz_rv , lock_iv , 1d0 , wsrmax &
            , omax1 , omax2 , rmt_rv )
       i_copy_size=1;
       nclspp = max(2*nclasp-nclass,nspec)
       allocate(rmax(nclspp))
       print *,' zzzz nclspp=',nclspp
       if(allocated(rv_a_ormax)) deallocate(rv_a_ormax) !is this correct???
       do i=1,nclspp
          rmax(i) = rmt_rv(iv_a_oics(i)) !sspec(iv_a_oics(i))%rmt
       enddo
       allocate(rv_a_ormax(nclspp))
       call dcopy ( nclspp , rmax , 1 , rv_a_ormax , 1 )
       ifx=1
    endif
    !-------
    if(ifx==0) then
       allocate(rmax(nclasp))
       call dcopy ( nclasp , rv_a_ormax , 1 , rmax , 1 )
    endif
    ! ... Write positions in Cartesian coordinates and as multiples plat
    if (iprint() >= 50) then
       write(stdo,357)
357    format(/' site spec',8x,'pos (Cartesian coordinates)',9x, &
            'pos (multiples of plat)')
       !     qlat = (plat+)^-1
       call dinv33(plat,1,qlat,xx)
       do  i = 1, nbas
          call dpscop ( rv_a_opos , xv , 3 , 3 * i - 2 , 1 , 1d0 )
          !       posp+ = (plat)^-1 pos+
          call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
          ip = ival ( iv_a_oips , i )
          print 345, i, slabl(ip), (xv(j),j=1,3), (xv(3+j),j=1,3)
345       format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
       enddo
    endif
    ! --- Print overlaps, optionally minimize wrt spec'd sites ---
    outs = ' '
    i = 6
    !swtmp = cmdopt('-mino',5,0,outs)
    !swtmp = cmdopt('--mino',6,0,outs)
    !if (swtmp) i = 7
    j = 1
    if (iprint() < 30) j = 0
    call ovmin ( outs ( i: ) , nbas , nbasp , alat , plat , rmax &
         , rmax , slabl , iclasst , modep , z , iv_a_ontab , iv_a_oiax , &
         rv_a_opos , j )
    !      endif
    ! --- Interpolate core to another mesh ---
    !      if (getdig(mode,4,2) .ne. 0) then
    !        call rx('patch clabl for call to coritp')
    !C       call coritp(nclass,nsp,w(oclabl),nrmsh,amsh,w(ormax))
    !      endif
    deallocate(lmxa,z)
    if (allocated(lockc_rv)) deallocate(lockc_rv)
    if (allocated(rmtc_rv)) deallocate(rmtc_rv)
    if (allocated(zc_rv)) deallocate(zc_rv)
    if (allocated(z_rv)) deallocate(z_rv)
    if (allocated(ips2_iv)) deallocate(ips2_iv)
    if (allocated(pos2_rv)) deallocate(pos2_rv)
    if (allocated(lock_iv)) deallocate(lock_iv)
    if (allocated(rmt_rv)) deallocate(rmt_rv)
    if (allocated(zz_rv)) deallocate(zz_rv)
  end subroutine lmaux

  !ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
  subroutine makrm0(opt,nspec,nbas,alat,plat,pos,slabl,ips,modep, lock,z,rmt)
    use m_lmfinit,only: nsp
    use m_freeat,only:freats
    use m_xclda,only:evxcv
    ! C- Estimate muffin-tin radii overlapping atomic potentials
    ! C ----------------------------------------------------------------------
    ! Ci Inputs
    ! Ci   opt  :specifies options in makrm0.
    ! Ci         :
    ! Ci         :1s digit : prescription for constructing potential
    ! Ci         :  The only option is mode 1 (see Remarks):
    ! Ci         :  set 1's digit opt=1.
    ! Ci         :
    ! Ci         :10s digit : quantity to overlap
    ! Ci         :  0 overlap electrostatic atomic potentials
    ! Ci         :  1 overlap estat atomic potentials + add total xc potential
    ! Ci         :    This option tends to produce slightly more uniform radii.
    ! Ci         :  2 overlap total atomic potentials
    ! Ci         :    This option probably doesn't make sense.
    ! Ci         :
    ! Ci         :100s digit specifies what MT radii are returned.
    ! Ci         :Initially MT radii are generated by site.
    ! Ci         :  0 return MT radii by site:
    ! Ci              rmt(i) corresponds to radius for site i
    ! Ci         :  1 average MT radii by species:
    ! Ci              rmt(i) corresponds to radius for species i
    ! Ci         :  2 smallest MT radii of this species
    ! Ci              rmt(i) corresponds to radius for species i
    ! Ci         :
    ! Ci         :1000s digit specifies whether nspec,slabl,ips refer to
    ! Ci                species or classes (used for printout only)
    ! Ci         :  0 quantities refer to species
    ! Ci         :  1 quantities refer to classes
    ! Ci   nspec :number of species (or classes)
    ! Ci   nbas  :size of basis
    ! Ci   alat  :length scale of lattice and basis vectors, a.u.
    ! Ci   plat  :primitive lattice vectors, in units of alat
    ! Ci   pos   :basis vectors, in units of alat
    ! Ci   slabl :vector of species (or class) labels
    ! Ci   ips   :species (or class) table: site ib belongs to species ips(ib)
    ! Ci   modep :integer vector of length 3 governing how lattice vectors
    ! Ci         :are shortened (shorps).  In particular,
    ! Ci         :0 suppresses shifts along plat(j)
    ! Ci         :2 shifts to minimize length of pos
    ! Ci   z     :table of nuclear charges by species (or class)
    ! Co Outputs
    ! Co   rmt
    ! Cl Local variables
    ! Cl         :
    ! Cr Remarks
    ! Cr   makrm0 attempts to determine an optimal set of muffin-tin radii so
    ! Cr   that the errors due to the muffin-tin approxiation are minimized.
    ! Cr   Inside the MT-spheres the potential should be as spherical as
    ! Cr   posssible, and the potential should be as constant as possible.
    ! Cr
    ! Cr   At present there is only one prescription (mode 1).
    ! Cr   This routine follows the methodology of the Stuttgart code potmax.
    ! Cr
    ! Cr   mode 1 : rmt are chosen as follows:
    ! Cr   1. The free-atom potential is constructed for all species.
    ! Cr   For each atom:
    ! Cr   2. A neighbor table is made for the atom.
    ! Cr   3. The connecting vectors to each neighbor make up a group
    ! Cr         of direction vectors stellating from that atom.
    ! Cr   4. The first occurence of a maximum in the potential along any
    ! Cr         of these vectors determines rmt for the atom
    ! Cr
    ! Cr   Finally, rmt may be averaged over species (or classes).
    ! Cu Updates
    ! Cu   21 Apr 02 First created.
    ! C ----------------------------------------------------------------------
    ! C     implicit none
    ! C Passed variables:
    character*8 slabl(*)
    integer opt,nspec,nbas,ips(nbas),modep(3),lock(*)
    double precision z(nspec),alat,plat(9),pos(3,nbas),rmt(*)
    ! C Local variables:
    character*8 spid,strn*80
    integer nrmx,nxi0,n0,niax,mxcsiz,npmx,pass
    parameter (nrmx=1501,nxi0=10,n0=10,niax=10,mxcsiz=4000,npmx=mxcsiz)
    integer nxi,nrmix(2),lxcfun,is,opt0,opt1,opt2,opt3
    integer idmod(n0),lmxa,nrmt,nr,ib,nttab,k,ir,np,ipr, nrspec
    integer ntab(nbas+1),iax(niax,mxcsiz),ipa(nbas)
    double precision qc,ccof,ceh,rmtl(nspec),rfoca,rsmfa,qcor(2),sumec,sumtc,eref,fac,ddot,fpi,etot
    double precision a(nspec),b(nspec),xx(5)
    double precision pnu(n0,2),pnz(n0,2),qat(n0,2)
    double precision hfc(nxi0,2),exi(nxi0),hfct(nxi0,2)
    double precision rtab(n0,2),etab(n0,2)
    double precision v(nrmx,nspec+1),rho(nrmx,nspec+1), rhoc(nrmx,nspec+1),rofi(nrmx*2),range(nbas), &
         rmti(nbas), vp(npmx,0:2),xp(3,npmx),rp(npmx),vxcp(npmx),excp(npmx)
    real(8) ,allocatable :: excx_rv(:)
    real(8) ,allocatable :: excc_rv(:)
    real(8) ,allocatable :: vxcx_rv(:)
    real(8) ,allocatable :: vxcc_rv(:)

    !C     Sets scale for neighbor table
    double precision facr,facri
    parameter (facr=4d0,facri=1.02d0) !enlarge facr because mp-632250 (H only failed)
    !C     Relative positions
    double precision rpos(3,mxcsiz),ri,rbar,rmin,vrmax,cur,slo,di
    character*1 sym(2)
    !C ... External calls
    !      external dcopy,defpq,defwsr,dpcopy,dpzero,dscal,getpr,
    !     .iinit,info,polint,poppr,pshpr,psymq0,psymr0,rx,rxi
    !C ... Heap
    data sym /' ','*'/

    logical:: isanrg, l_dummy_isanrg
    integer:: ifives

    print *,'makrm0:'
    !C --- Setup ---
    opt0 = mod(opt,10)
    opt1 = mod(opt/10,10)
    opt2 = mod(opt/100,10)
    opt3 = mod(opt/1000,10)
    l_dummy_isanrg=isanrg(opt0,1,1,'makrm0:','1s digit opt',.true.)
    l_dummy_isanrg=isanrg(opt1,0,2,'makrm0:','10s digit opt',.true.)
    l_dummy_isanrg=isanrg(opt2,0,1,'makrm0:','100s digit opt',.true.)
    l_dummy_isanrg=isanrg(opt3,0,1,'makrm0:','1000s digit opt',.true.)
    !c      stdo = lgunit(1)
    call getpr(ipr)
    call dpzero(rpos,3*mxcsiz)
    call dpzero(v,nrmx*nspec)
    call dpzero(rho,nrmx*nspec)
    call dpzero(rhoc,nrmx*nspec)
    call dpzero(vp,npmx*3)
    fpi = 16d0*datan(1d0)

    !C --- Densities and potentials for each species ---
    exi(1) = -1
    exi(2) = -2
    exi(3) = -4
    exi(4) = -6
    exi(5) = -9
    exi(6) = -15
    nxi = 6
    call dpzero(hfct,2*nxi0)
    nrmix(1) = 80
    nrmix(2) = 2
    lxcfun = 1
    do  is = 1, nspec

       !C       nrspec(is) = iabs(iclbsj(is,ips,-nbas,nbas))
       spid = slabl(is)
       rsmfa = 1
       rfoca = 1
       qcor(1) = 0
       qcor(2) = 0
       a(is) = .025d0
       call defwsr(rmtl(is),z(is))
       rsmfa = rmtl(is)/2
       rfoca = rmtl(is)/2
       lmxa = 3
       call dpzero(pnu,2*n0)
       call dpzero(pnz,2*n0)
       call dpzero(qat,2*n0)
       call iinit(idmod,n0)
       call defpq(z(is),lmxa,nsp,pnu,qat)
       eref = 0
       nrmt = 0
       call pshpr(ipr-20)
       call freats(spid,is,nxi0,nxi,exi,rfoca,rsmfa,0,-1,qcor,nrmix(1),0, &
            lxcfun,z(is),rmtl(is),a(is),nrmt,pnu,pnz,qat,0d0,0d0,0d0,[0d0,0d0], &
            idmod,lmxa,eref,rtab,etab,hfc,hfct,nr,rofi,rho(1,is),rhoc(1, &
            is),qc,ccof,ceh,sumec,sumtc,v(1,is),etot, 1, -1,-1) !nmcore=1 july2012 !last -1 -1 means not write ves* files.
       call poppr
       b(is) = rmtl(is)/(dexp(a(is)*(nrmt-1)) - 1d0)
       !C       Restore 4*pi*r**2*(total density) in array rho
       call daxpy(nrmx*nsp,1d0,rhoc(1,is),1,rho(1,is),1)

       !C       Use estat potential: overwrite v with estat
       if (opt1 .le. 1) then
          !C         call prrmsh('vtot',rofi,v(1,is),nr,nr,1)
          vrmax = 0
          call poiss0(z(is),a(is),b(is),rofi,rho(1,is),nr,vrmax,v(1,is), &
               xx(2),xx(4),nsp)
          !C         call prrmsh('ves',rofi,v(1,is),nr,nr,1)
       endif

       if (nsp .eq. 2) then
          call daxpy(nr,1d0,rho(1+nr,is),1,rho(1,is),1)
          call daxpy(nr,1d0,rhoc(1+nr,is),1,rhoc(1,is),1)
          call daxpy(nr,1d0,v(1+nr,is),1,v(1,is),1)
          call dscal(nr,.5d0,v(1,is),1)
          call rx('need check this branch')
       endif

       !C       Overwrite v,rho with proper spherical potential and density
       do  ir = 2, nr
          v(ir,is) = v(ir,is) - 2*z(is)/rofi(ir)
          rho(ir,is) = rho(ir,is)/(fpi*rofi(ir)**2)
       enddo
       !C       call prrmsh('rho',rofi,rho(1,is),nr,nr,1)
       !C       call prrmsh('v',rofi,v(1,is),nr,nr,1)
    enddo
    !C

    !C ... Sets range for neighbor table
    do  ib = 1, nbas
       is = ips(ib)
       range(ib) = facr * rmtl(is)
       print *,'rrrrrrr',ib,range(ib)
    enddo

    !C --- For each site, find rmti = initial estimate for rmtl ---
    do  ib = 1, nbas

       is = ips(ib)
       spid = slabl(is)

       !C   ... Neighbor table connecting sites to this one
       nttab = mxcsiz
       call pairc(ib,ib,nbas,modep,20,[0],alat,plat,pos,pos,range,-1,[1], nttab,ntab,iax,rpos,k)
       if (nttab.gt.mxcsiz)  call rxi('makrm0 : increase mxcsiz: need',nttab)

       !C   ... Maximum in potential ( mode 1)
       !C       Find first maximum along all direction vectors
       !C       which are proportional to rpos, excluding self (first site).
       !C       As a conservative choice, choose initial radius = rmtl/4.
       np = nttab-1
       if (nttab .gt. npmx)  call rxi('makrm0 : increase npmx: need',nttab)
       call dpzero(xp,3*nttab)
       ri = rmtl(is)/4
       if (z(is) .eq. 0) then
          rmti(ib) = rmtl(is)
          goto 31
       endif

       !C       Re-entry for loop until maximum in potential found
       pass = 0
30     continue
       !C       Make nttab-1 connecting vectors of radius ri
       do  k = 2, nttab
          fac = ri/dsqrt(ddot(3,rpos(1,k),1,rpos(1,k),1))
          call dpcopy(rpos(1,k),xp(1,k-1),1,3,fac)
       enddo

       !C       Add superposition of potentials for each vector
       call sumsro(xp,np,ips,a,b,v,nttab,iax,rpos,vp(1,2))

       !C       Add vxc[superposition of densities] for each vector
       if (opt1 .eq. 1) then
          call sumsro(xp,np,ips,a,b,rho,nttab,iax,rpos,rp)
          allocate(excx_rv(np))
          allocate(excc_rv(np))
          allocate(vxcx_rv(np))
          allocate(vxcc_rv(np))
          call evxcv ( rp , rp , np , 1 , lxcfun , excp , excx_rv ,   excc_rv , vxcp , vxcx_rv , vxcc_rv )
          call daxpy(np,1d0,vxcp,1,vp(1,2),1)
          deallocate(vxcc_rv)
          deallocate(vxcx_rv)
          deallocate(excc_rv)
          deallocate(excx_rv)
       endif

       !C       Check against prior point.  Skip if nothing yet to compare
       if (pass .gt. 0) then
          do  k = 1, nttab-1
             !C           This loop executes if a maximum is found
             if (vp(k,1) .gt. vp(k,2)) then
                !C             It shouldn't happen already at the first point
                if (pass .eq. 1) then
                   call   rx('makrm0: starting r is too small.  Check geometry')
                endif
                rmti(ib) = ri/dsqrt(facri)
                if (pass .gt. 1) then
                   cur = vp(k,2) + vp(k,0) - 2*vp(k,1)
                   slo = (vp(k,2) - vp(k,0))/2
                   di = -slo/cur
                   if (dabs(di) .lt. 1) then
                      !C                 print *, pass, di, ri, (ri/facri)*dexp(di*(facri-1))
                      rmti(ib) = (ri/facri)*dexp(di*(facri-1))
                   endif
                endif
                goto 31
             endif
          enddo
       endif

       !C       Copy vp(*,2) to vp(*,1), and increment pass
       call dcopy(npmx,vp(1,1),1,vp(1,0),1)
       call dcopy(npmx,vp(1,2),1,vp(1,1),1)
       pass = pass+1
       ri = ri * facri
       if (ri .gt. 10d0)  call rx('makrm0: failed to find maximum in potential')
       goto 30

       !C       Loop exit: radius rmti(ib) has been found
31     continue

       !C       This is the connecting vectory along which the maximum was found
       !C       iconn = iax(2,k+1)

    enddo

    !C --- Copy appropriate rmti to final result, depending on opt2 ---
    strn = 'estat'
    if (opt1 .eq. 1) strn = 'LDA'
    if (opt1 .eq. 2) strn = 'sum-of-atomic'
    call info(30,1,0,' makrm0: initial MT radii from first '//strn(1:20)//'%a potential maximum',0,0)
    call info(30,0,0,'  site   spec%12frmt'//'%7frmt-%7frmt-%7frold%3flock%N%34f<spec avg>  spec-min',0,0)

    !C ... First loop for printout only, so can hang on to original rmt
    do  ib = 1, nbas
       is = ips(ib)
       spid = slabl(is)
       call psymr0(-2,-is,nbas,ips,xx,xx,ipa,nrspec)
       rbar = 0
       rmin = 9999
       do  k = 1, nrspec
          rbar = rbar + rmti(ipa(k))/nrspec
          rmin = min(rmin,rmti(ipa(k)))
       enddo
       !C       Index where to poke rnew
       k = ib
       if (opt2 .ge. 1) k = is

       if (ipr .ge. 30) then
          write(stdo,344) ib,is,spid,rmti(ib),rmti(ib)-rbar, rmti(ib)-rmin,rmt(k), sym(1+mod(lock(k)/2,2))
344       format(i5,2x,i3,':',a,4f11.4,3x,a)
       endif

    enddo

    !C     print *, rmt

    do  ib = 1, nbas
       is = ips(ib)

       spid = slabl(is)
       call psymr0(-2,-is,nbas,ips,xx,xx,ipa,nrspec)
       rbar = 0
       rmin = 9999
       do  k = 1, nrspec
          rbar = rbar + rmti(ipa(k))/nrspec
          rmin = min(rmin,rmti(ipa(k)))
       enddo
       !C       Index where to poke rnew
       k = ib
       if (opt2 .ge. 1) k = is

       if (lock(k) .ne. 2) then
          if (opt2 .eq. 0) rmt(ib) = rmti(ib)
          if (opt2 .eq. 1) rmt(is) = rbar
          if (opt2 .eq. 2) rmt(is) = rmin
       endif
    enddo
    !C    print *, rmt
  end subroutine makrm0
  subroutine sclwsr(opt,nbas,nbasp,nspec,alat,plat,basp,ips,modep, &
       slabl,z,lock,volfac,wsmax,omax1,omax2,wsr)
    use m_ext,only: sname
    !- Scales sphere radii to specified volume, satisfying constraints.
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   opt  :specifies options in sclwsr.
    !i         :1s digit specifies whether scaling reaches target volfac
    !i         :         is required:
    !i         : 0 meeting volume target is optional
    !i         : 1 meeting volume target is optional but
    !i         :   a warning is printed if not met.
    !i         : 2 meeting volume target is required
    !i         :10s digit concerns treatment of empty spheres
    !i         : 0 ES and other sites are treated symmetrically
    !i         : 1 all sites with z>0 are resized first; then
    !i             all sites are resized.
    !i         : 2 all sites with z>0 are resized first; then
    !i             the ES sites only are resized.
    !i   nbas  :size of basis
    !i   nbasp :size of padded basis (layer programs)
    !i   nspec :number of atom species
    !i   alat  :length scale
    !i   plat  :primitive lattice vectors, in units of alat
    !i   basp  :basis vectors (scaled by alat; padded for layer mode)
    !i   ips   :the jth atom belongs to species ips(j)
    !i   modep :specifies which dimensions have periodic b.c.
    !i         :In particular,
    !i         :0 suppresses shifts along plat(j)
    !i         :2 shifts to minimize length of pos
    !i          are used
    !i   slabl :species labels
    !i   z     :nuclear charge, by species
    !i   lock  :constraints specifying which species are locked and
    !i         :which are free to float.
    !i         :On input, lock should be zero or two for each species.
    !i         :Each species for which lock(i)=2 is constrained not to change.
    !i         :The radii for other species are floated.
    !i         :lock(1..nspec) is OVERWRITTEN on output
    !i   volfac:scale until sum of sphere vol = volfac * cell vol or until
    !i         :all sphere radii are constrained (see Remarks, sclws2)
    !i   wsmax :(wsmax>0) a global maximum on the size of MT spheres.
    !i         :No sphere is allowed to exceed wsmax.
    !i         :(wsmax=0) no constraint is imposed.
    !i   omax1 :max. allowed overlap divided by distance (r1+r2-d)/d<omax1
    !i         :omax1(1) is constraint for atom-atom overlaps
    !i         :omax1(2) is constraint for atom-empty-sphere overlaps
    !i         :omax1(3) is constraint for empty-sphere-empty-sphere overlaps
    !i   omax2 :max. allowed overlap divided by radius  (r1+r2-d)/r1<omax2
    !i         :omax2(1) is constraint for atom-atom overlaps
    !i         :omax2(2) is constraint for atom-empty-sphere overlaps
    !i         :omax2(3) is constraint for empty-sphere-empty-sphere overlaps
    !o Inputs/Outputs:
    ! o  wsr   :Wigner-Seitz sphere radius (in atomic units)
    ! o        :On input, starting values for wsr
    ! o        :On output, final values for wsr
    !l Local:
    !l   dovl1 :maximum overlap divided by distance
    !l   dovl2 :maximum overlap divided by radius
    !l   gamma :(zero) passed to sclws2
    !l   nrspec:number of atoms in each species
    !r Remarks
    !r  Sphere radii are scaled in an iterative procedure.  In any
    !r  iteration, species are divided into those that are 'locked'
    !r  (frozen) and those that are allowed to float.  The largest scaling
    !r  factor is determined for all those species of the latter type, that
    !r  satisfies the constraints (see below).  These species are scaled
    !r  and a new iteration begins.  By construction, each iteration will
    !r  cause at least one new species to be locked; thus, the total number
    !r  of iterations will not exceed the number of species.
    !u Updates
    !u   17 Jan 09  bug fix 10s digit mode=2, no ES
    !u   17 May 02  New 10s digit opt switch gives optional lower priority
    !u              to empty spheres
    !u   22 Apr 02  First created; adapted from Stuttgart LMTO56.
    ! ----------------------------------------------------------------------
    implicit none
    ! Passed variables:
    integer :: opt,nbas,nbasp,nspec,ips(nbas),modep(3),lock(nspec)
    double precision :: alat,basp(3,nbasp),volfac,wsmax, &
         omax1(3),omax2(3),plat(3,3),wsr(nspec),z(nspec)
    character(8) :: slabl(*)
    integer:: i , k , nrspec(nspec) , ib , is , opt1
    real(8) ,allocatable :: wk_rv(:)
    integer :: niax,mxcsiz,mxnbr,nttab,ipr,llock(nspec)
    integer ,allocatable :: ntab_iv(:)
    integer ,allocatable :: iax_iv(:)
    logical :: les
    double precision :: dovl1(3),dovl2(3),gamma,range(nbas),wsrs(nspec)
    double precision :: facr,tiny,avw,vol,avwsr,volnew,volnes,volold
    parameter(niax=10,mxcsiz=2000,facr=2d0,tiny=1d-5)
    integer:: istdo

    !     omax1(1) = -.01
    !      wsr(1) = 3.232247d0
    !      wsr(2) = 3.232247d0
    !      wsr(3) = 2.248243d0
    !      wsr(4) = 2.097284d0
    !      wsr(5) = 1.647651d0
    !      wsr(6) = 1.647651d0

    ! --- Setup ---
    call getpr(ipr)
    !      stdo  = lgunit(1)
    gamma = 0
    !      call maknrs(nbas,ips,ib,nrspec)
    nrspec=0
    do  ib = 1, nbas
       is = ips(ib)
       nrspec(is) = nrspec(is) + 1
    enddo
    !      if (ib .gt. nspec) call rx('sclwsr: wrong number of species')
    opt1 = mod(opt/10,10)
    call dcopy(nspec,wsr,1,wsrs,1)
    avw = avwsr(plat,alat,vol,nbas)

    ! --- Make a neighbor table and adjust llock to freeze ES sites ---
    do  ib = 1, nbas
       is = ips(ib)
       range(ib) = facr * wsr(is)
    enddo
    nbasp = nbas
    mxnbr = mxcsiz*nbas
    allocate(ntab_iv(nbasp+1))

    allocate(iax_iv(abs(-niax*mxnbr)))
    if (-niax*mxnbr<0) iax_iv(:)=0

    allocate(wk_rv(3*mxnbr))

    nttab = mxnbr
    call pshpr(ipr-20)
    call pairc ( 1 , nbas , nbasp , modep , 20 , [0] , alat , plat &
         , basp , basp , range , - 1 , [1] , nttab , ntab_iv , iax_iv &
         , wk_rv , k )

    call poppr

    ! --- Scale the sphere radii, freezing empty spheres ---
    if (opt1 >= 1) then

       !   ... Local copy of lock and wsr, adjusting to freeze ES sites
       call icopy(nspec,lock,1,llock,1)
       les = .false.
       do  is  = 1, nspec
          if (z(is) == 0) then
             les = .true.
             llock(is) = 2
             wsr(is) = 0
          endif
       enddo

       if (les) then

          !     ... Scale wsr with z=0 sites locked at wsr=0
          volold = volsph(nspec,nrspec,wsr)/vol
          volnes = volfac
          call sclws2 ( nbas , nspec , alat , plat , basp , slabl , iax_iv &
               , ips , ntab_iv , z , nrspec , omax1 , omax2 , gamma , wsmax &
               , llock , volnes , dovl1 , dovl2 , wsr )

          if (ipr >= 30) then
             !            call awrit2(' SCLWSR:  initial sphere packing = %;1d%%'//
             !     .      ' scaled to %;1d%% (no empty spheres)',
             !     .      ' ',120,stdo,100*volold,100*volnes)
             write(stdo,"(' SCLWSR:  Initial sphere packing = ',f10.5, &
                  ' scaled to (no empty spheres)',f10.5)")  100*volold,100*volnes
          endif

          !     ... Restore wsr(Z=0)
          call icopy(nspec,lock,1,llock,1)
          do  is  = 1, nspec
             if (z(is) == 0) wsr(is) = wsrs(is)
          enddo
       endif

       !   ... 10s digit opt=2 : freeze wsr(Z>0)
       if (opt1 >= 2 .AND. les) then
          do  is  = 1, nspec
             if (z(is) /= 0) lock(is) = 2
          enddo
       endif

    endif

    ! --- Scale the sphere radii (2nd pass for les) ---
    volold = volsph(nspec,nrspec,wsrs)/vol
    volnew = volfac
    call sclws2 ( nbas , nspec , alat , plat , basp , slabl , iax_iv &
         , ips , ntab_iv , z , nrspec , omax1 , omax2 , gamma , wsmax &
         , lock , volnew , dovl1 , dovl2 , wsr )


    ! --- Printout ---
    ! i#error, have return with len(w_varlist)>0 at line 183
    if ( ipr < 10 ) then
       if (allocated(wk_rv)) deallocate(wk_rv)
       if (allocated(iax_iv)) deallocate(iax_iv)
       if (allocated(ntab_iv)) deallocate(ntab_iv)
       return
    endif

    !      if (ipr .ge. 10) write(stdo,309)
    !     .  vol,100*volold,100*volnew
    !  309 format(/' SCLWSR: vol=',f11.3,
    !     .  '  sphere fraction=',f5.1,
    !     .  '%(initial)  ',f5.1,'%(scaled)')

    if (ipr >= 10) then
       write(stdo,"(' SCLWSR:  vol = ',f10.5,' a.u.**3 ', &
            ' Initial sphere packing = ',f10.5, &
            ' scaled to ',f10.5)")   vol,100*volold,100*volnew
    endif

    if (ipr >= 30) then
       write(stdo,310) &
            (omax1(i)*100,i=1,3),(omax2(i)*100,i=1,3), &
            (dovl1(i)*100,i=1,3),(dovl2(i)*100,i=1,3)
310    format(1x,'constr omax1=',3f6.1,' %    omax2=',3f6.1,' %', &
            /1x, 'actual omax1=',3f6.1,' %    omax2=',3f6.1,' %')
       !        if (ipr .gt. 30) then
       write(stdo,311)
311    format(/' spec  name',8x,'old rmax    new rmax     ratio')
       do  i = 1, nspec
          write(stdo,312) i,slabl(i), wsrs(i), wsr(i), wsr(i)/wsrs(i)
312       format(i4,3x,a,3f12.6)
       enddo
       !        endif
    endif

    open(newunit=istdo,file='rmt.'//trim(sname))
    do  i = 1, nspec
       write(istdo,"(a,f12.6)") slabl(i), wsr(i)
    enddo
    close(istdo)

    if (dabs(volnew-volfac) > tiny) then
       if (mod(opt,10) >= 1) then
          write(stdo,321) int(volfac*100)
321       format(/' SCLWSR (warning): failed to reach target vol (', &
               i3,'% of cell vol)')
       endif
       if (mod(opt,10) == 2) then
          call rx('SCLWSR: failed to reach target VOL.  Increase omax.')
       endif
    endif
    if (allocated(ntab_iv)) deallocate(ntab_iv)
    if (allocated(iax_iv)) deallocate(iax_iv)
    if (allocated(wk_rv)) deallocate(wk_rv)
  end subroutine sclwsr


  subroutine sclws2(nbas,nspec,alat,plat,bas,slabl,iax,ips,ntab,z, &
       nrspec,omax1,omax2,gamma,wsmax,lock,volfac,dovl1,dovl2,wsr)
    !- Enlarges the spheres to reach a specified volume under constraints
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   nbas  :size of basis
    !i   nspec :number of species
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   plat  :primitive lattice vectors, in units of alat
    !i   bas   :basis vectors, in units of alat
    !i   slabl :species labels (for printout)
    !i   iax   :neighbor table containing pair information (pairc.f)
    !i   ips   :species table: site ib belongs to species ips(ib)
    !i   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
    !i   z     :nuclear charge, by species
    !i   nrspec:number of atoms in the ith species
    !i   omax1 :max. allowed overlap divided by distance (s1+s2-d)/d<omax1
    !i   omax2 :max. allowed overlap divided by radius  (s1+s2-d)/s1<omax2
    !i   gamma :a factor that changes scaling wsr from simple multiplicative
    !i         :scaling to a combination of additive + multiplicative scaling
    !i         :That is, in each iteration,
    !i         :scaling is r -> a(r+b) with a*b=gamma*(a-1)*avw
    !i         :gamma>0 tends enlarge small spheres faster than large ones
    !i   wsmax :(wsmax>0) a global maximum on the size of MT spheres.
    !i         :No sphere is allowed to exceed wsmax.
    !i         :(wsmax=0) no constraint is imposed.
    ! o Inputs/Outputs:
    ! o  lock  :constraints specifying which species are locked and
    ! o        :which are free to float.
    ! o        :On input, lock should be zero or two for each species.
    ! o        :Each species for which lock(i)=2 is constrained not to change.
    ! o        :The radii for other species are floated.
    ! o        :lock(1..nspec) is OVERWRITTEN on output
    ! o  volfac:scale until sum of sphere vol = volfac * cell vol
    ! o        :or until all sphere radii are constrained (see Remarks)
    ! o        :On input, volfac=target ratio (sum of sphere vol)/(cell vol)
    ! o        :On output, volfac=actual ratio
    ! o  wsr   :Wigner-Seitz sphere radius (in atomic units)
    ! o        :On input, starting values.
    ! o        :On output, scaled values.
    !o Outputs:
    !o   dovl1 :maximum overlap divided by distance
    !o   dovl2 :maximum overlap divided by radius
    !r Remarks:
    !r  Sphere radii are scaled in an iterative procedure.  In any
    !r  iteration, species are divided into those that are 'locked'
    !r  (frozen) and those that are allowed to float.  The largest scaling
    !r  factor is determined for all those species of the latter type, that
    !r  satisfies the constraints (see below).  These species are scaled
    !r  and a new iteration begins.  By construction, each iteration will
    !r  cause at least one new species to be locked; thus, the total number
    !r  of iterations will not exceed the number of species.
    !r
    !r  Typically sclws2 is called by is a higher-level routine that creates
    !r  the neighbor-table and other necessary arrays.
    !r
    !r  This code was adapted from Stuttgart routine blowup, v LMTO56.
    ! ----------------------------------------------------------------------
    !     implicit none
    ! Passed variables:
    integer :: niax
    parameter (niax=10)
    integer :: nbas,nspec,iax(niax,*),ntab(*),ips(nbas),nrspec(nspec)
    integer :: lock(nspec)
    double precision :: alat,bas(3,nbas),dovl1(3),dovl2(3),wsmax, &
         volfac,gamma,omax1(3),omax2(3),plat(3,3),wsr(nspec),z(nspec)
    character(8) :: slabl(nspec)
    ! Local variables:
    integer :: ib,jb,is,iclbsj,ip,k1,k2,k3, &
         kb,ks,kco,kpr,nloop,stdo,ipr,npcol,locki
    double precision :: a,amax,amax1,amax2,amax3,amax4,avw,b,bmax,d, &
         dm(0:3),dovlap,dscl(nspec), &
         dsclmx,dr(3),fpi3,gw,opo1,omo2,p,q,r,ratio, &
         rik,riko,s,t,tiny,u,v,vol,vola,volb, &
         wsri,wsrk,x,avwsr,Vconst
    logical :: fin
    character(72) :: fmt
    parameter(fpi3=4.18879020478639053d0,tiny=1d-5)

    avw  = avwsr(plat,alat,vol,nbas)
    gw   = avw*gamma
    !      stdo = lgunit(1)
    call getpr(ipr)
    !      print *, '!!'
    !      ipr = 40
    npcol = 7

    !     call dcopy(nspec,1d0,0,dscl(1,0),1)

    fmt = '(''   SPEC:    '',6(1x,a8):/(12x,6(1x,a8)))'
    write(fmt(17:17),'(I1)') npcol
    write(fmt(32:32),'(I1)') npcol
    if (ipr >= 50) then
       write(stdo,'('' '')')
       write(stdo,fmt) (slabl(is),is=1,nspec)
    endif
    fmt = '(1x,''init rmt'',1x,6f9.5:/(10x,6f9.5))'
    write(fmt(19:19),'(I1)') npcol
    write(fmt(31:31),'(I1)') npcol
    !     print *, fmt
    if (ipr >= 50) write(stdo,fmt) (wsr(is),is=1,nspec)

    do  nloop = 1, nspec+1
       amax = 9d9
       do  is = 1, nspec
          if (lock(is) /= 2) lock(is) = 0
       enddo

       !   --- Lock radii of those spheres with maximum allowed overlap ---
       !       and unlock those with radii exceeding maximum allowed.
       do  is = 1, nspec
          wsri = wsr(is)

          !     ... If overlap criterion exactly satisfied for any connecting vector,
          !         set lock for this species
          do  jb = 1, nrspec(is)
             ib = iclbsj(is,ips,nbas,jb)
             do  kpr = ntab(ib)+2, ntab(ib+1)
                kb = iax(2,kpr)
                k1 = iax(3,kpr)
                k2 = iax(4,kpr)
                k3 = iax(5,kpr)
                ks = ips(kb)
                !           ip selects which ommax to use (A-A, A-E, or E-E)
                ip = 2
                if (idnint(z(is)) /= 0 .AND. idnint(z(ks)) /= 0) ip=1
                if (idnint(z(is)) == 0 .AND. idnint(z(ks)) == 0) ip=3
                opo1 = 1+omax1(ip)
                omo2 = 1-omax2(ip)
                wsrk = wsr(ks)
                rik = dsqrt(drr2(plat,bas(1,ib),bas(1,kb),k1,k2,k3,dr))
                rik = rik*alat
                !           Set lock if any of these conditions are met:
                !             wi + wk - rik = o1*rik
                !             wi + wk - rik = o2*wi
                !             wi + wk - rik = o2*wk
                !             wi = wsmax
                locki = lock(is)
                if (dabs(opo1*rik-wsri-wsrk) < tiny) locki=1
                if (dabs(rik-omo2*wsri-wsrk) < tiny) locki=1
                if (dabs(rik-wsri-omo2*wsrk) < tiny) locki=1
                if (wsmax > 0 .AND. dabs(wsri-wsmax) < tiny) locki=1
                lock(is) = max(locki,lock(is))
             enddo
          enddo
          !         print *,'is,lock0=',is,lock(is)

          !     ... If overlap criterion exceeded for any connecting vector,
          !         unset lock for this species
          if (lock(is) /= 2) then
             do  jb = 1, nrspec(is)
                ib = iclbsj(is,ips,nbas,jb)
                do  kpr = ntab(ib)+2, ntab(ib+1)
                   kb = iax(2,kpr)
                   k1 = iax(3,kpr)
                   k2 = iax(4,kpr)
                   k3 = iax(5,kpr)
                   ks = ips(kb)
                   wsrk = wsr(ks)
                   rik = dsqrt(drr2(plat,bas(1,ib),bas(1,kb),k1,k2,k3,dr))
                   rik = rik*alat
                   !           ip selects which ommax to use (A-A, A-E, or E-E)
                   ip = 2
                   if (idnint(z(is)) /= 0 .AND. idnint(z(ks)) /= 0) ip=1
                   if (idnint(z(is)) == 0 .AND. idnint(z(ks)) == 0) ip=3
                   opo1 = 1+omax1(ip)
                   omo2 = 1-omax2(ip)
                   !           Unset lock if any of these conditions are exceeded
                   !             wi + wk - rik < o1*rik
                   !             wi + wk - rik < o2*wi
                   !             wi + wk - rik < o2*wk
                   !             wi > wsmax
                   if (opo1*rik-wsri-wsrk < -tiny) lock(is)=0
                   if (rik-omo2*wsri-wsrk < -tiny) lock(is)=0
                   if (rik-wsri-omo2*wsrk < -tiny) lock(is)=0
                   if (wsmax > 0 .AND. wsri-wsmax > tiny) lock(is)=0
                enddo
             enddo
          endif
       enddo
       !       print *, 'lock', (lock(is), is=1,nspec)

       !  --- Find amax=largest allowed scaling for unlocked species ---
       do  is = 1, nspec
          if (lock(is) == 0) then
             riko = -1
             kco = -1
             wsri = wsr(is)
             do  jb = 1, nrspec(is)
                ib = iclbsj(is,ips,nbas,jb)
                do  kpr = ntab(ib)+2, ntab(ib+1)
                   kb = iax(2,kpr)
                   k1 = iax(3,kpr)
                   k2 = iax(4,kpr)
                   k3 = iax(5,kpr)
                   ks = ips(kb)
                   rik = dsqrt(drr2(plat,bas(1,ib),bas(1,kb),k1,k2,k3,dr))
                   rik = rik*alat
                   ip = 2
                   if (idnint(z(is)) /= 0 .AND. idnint(z(ks)) /= 0) ip=1
                   if (idnint(z(is)) == 0 .AND. idnint(z(ks)) == 0) ip=3
                   opo1 = 1+omax1(ip)
                   omo2 = 1-omax2(ip)
                   if (dabs(rik-riko) > tiny .OR. kco /= ks) then
                      wsrk = wsr(ks)
                      riko = rik
                      kco = ks
                      amax1 = 9d9
                      amax2 = 9d9
                      amax3 = 9d9
                      amax4 = 9d9
                      !               If second site ks is locked, can only scale site i
                      if (lock(ks) /= 0) then
                         amax1 = (opo1*rik-wsrk+gw)/(wsri+gw)
                         if (omo2 > 0d0) &
                              amax2 = (rik-wsrk+omo2*gw)/(omo2*wsri+omo2*gw)
                         amax3 = (rik-omo2*wsrk+gw)/(wsri+gw)
                         !               If neither site is locked, both will scale
                      else
                         amax1 = (opo1*rik+gw+gw)/(wsri+wsrk+gw+gw)
                         if (wsrk+omo2*wsri+(1+omo2)*gw > 0d0) &
                              amax2 = (rik+(1+omo2)*gw)/ &
                              (wsrk+omo2*wsri+(1+omo2)*gw)
                         if (wsri+omo2*wsrk+(1+omo2)*gw > 0d0) &
                              amax3 = (rik+(1+omo2)*gw)/ &
                              (wsri+omo2*wsrk+(1+omo2)*gw)
                      endif
                      if (wsmax > 0) then
                         amax4 = wsmax/wsri
                      endif
                      amax = dmin1(amax,amax1,amax2,amax3,amax4)
                      !               print *, 'kc,kbas,kpr,amax',ks,kb,kpr,amax
                   endif
                enddo
             enddo
          endif
       enddo
       bmax = (1d0-1/amax)*gw

       !   ... Determine what new volume will be after scaling with a,b
       vola = 0d0
       volb = 0d0
       do  ib = 1, nbas
          is = ips(ib)
          if (lock(is) == 0) then
             volb = volb + (amax*wsr(is)+bmax)**3
          else
             vola = vola + wsr(is)**3
          endif
       enddo
       vola = vola * fpi3
       volb = volb * fpi3

       !   --- Case scaling will lead to new volume > final volume ---
       if (vol*volfac < vola+volb) then
          call dpzero(dm,4)
          do  ib = 1, nbas
             is = ips(ib)
             if (lock(is) == 0) then
                a = wsr(is)+gw
                !             For numerical reasons distinguish cases
                if (dabs(gamma) > 1d0) then
                   b = wsr(is)
                else
                   b = -gw
                endif
                !             fpi3*dm(0) = sum of sphere volumes not constrained
                dm(0) = dm(0) +       b*b*b
                dm(1) = dm(1) + 3d0 * a*b*b
                dm(2) = dm(2) + 3d0 * a*a*b
                dm(3) = dm(3) +       a*a*a
             endif
          enddo
          !         Vconst + a**3*Vuncst = Vtarget; Vconst + Vuncst = vola
          !         If Vconst > Vtarget, constraint cannot be satisified
          Vconst = vola - fpi3*dm(3)
          if (Vconst > vol*volfac) then
             call fexit2(-1,111,' Exit -1 : SCLWSR: constrained '// &
                  'sphere vol (%;0d) exceeds target vol=%;0d', &
                  Vconst,vol*volfac)
          endif
          if (dabs(dm(3)) > tiny) then
             r =  dm(2) / dm(3)
             s =  dm(1) / dm(3)
             t = (dm(0) - (vol*volfac-vola)/fpi3) / dm(3)
             p = s - r*r/3d0
             q = 2d0*r*r*r/27d0 - r*s/3d0 + t
             d = p*p*p/27d0 + q*q/4d0
             u = (dsqrt(d)-q/2d0)**(1d0/3d0)
             v = -p/u/3d0
             x = u+v-r/3d0
             if (dabs(gamma) > 1d0) then
                amax = x+1d0
                bmax = x*gw/amax
             else
                amax = x
                bmax = (1d0-1d0/amax)*gw
             endif
             if (ipr >= 100) then
                write(stdo,300)'R S T',r,s,t
                write(stdo,300)'P Q  ',p,q
                write(stdo,300)'  D  ',d
                write(stdo,300)' U V ',u,v
                write(stdo,300)' AMAX',amax
                write(stdo,300)' BMAX',bmax
                write(stdo,300)' -------------------------'
             endif
          endif
       endif

       !   --- Scale unlocked spheres by w <- a*w+b ---
       fin = .true.
       do  is = 1, nspec
          dscl(is) = 1d0
          if (lock(is) == 0) then
             dsclmx  = amax + bmax/wsr(is)
             wsr(is) = dsclmx*wsr(is)
             dscl(is)= dsclmx
             fin = .false.
          endif
       enddo

       !       ratio = new volume / old volume
       ratio = volsph(nspec,nrspec,wsr)/vol
       !       fin=T if volume change is small or all spheres are locked
       fin = fin .or. dabs(ratio-volfac).lt.tiny .or. nloop.eq.nspec+1

       fmt = '(1x,''iter:'',i3,1x,6f9.5:/(10x,6f9.5))'
       !       print *, fmt
       write(fmt(19:19),'(I1)') npcol
       write(fmt(31:31),'(I1)') npcol
       !       print *, fmt
       if (ipr >= 50) write(stdo,fmt) nloop,(dscl(is),is=1,nspec)

       !   ... Last iteration
       if (fin) then
          call dpzero(dovl1,3)
          call dpzero(dovl2,3)
          do  is = 1, nspec
             wsri = wsr(is)
             do  jb = 1, nrspec(is)
                ib = iclbsj(is,ips,nbas,jb)
                do  kpr = ntab(ib)+2, ntab(ib+1)
                   kb = iax(2,kpr)
                   k1 = iax(3,kpr)
                   k2 = iax(4,kpr)
                   k3 = iax(5,kpr)
                   ks = ips(kb)
                   ip = 2
                   if (idnint(z(is)) /= 0 .AND. idnint(z(ks)) /= 0) ip=1
                   if (idnint(z(is)) == 0 .AND. idnint(z(ks)) == 0) ip=3
                   wsrk = wsr(ks)
                   rik = dsqrt(drr2(plat,bas(1,ib),bas(1,kb),k1,k2,k3,dr))
                   rik = rik*alat
                   dovlap = wsri+wsrk-rik
                   dovl1(ip) = dmax1(dovl1(ip),dovlap/rik)
                   dovl2(ip) = dmax1(dovl2(ip),dovlap/wsri,dovlap/wsrk)
                   !              if (dovl1(ip).gt.omax1(ip)+tiny .or.
                   !     .            dovl2(ip).gt.omax2(ip)+tiny)
                   !     .          write(stdo,'('' SCLWS2: warning dovl>omax'')')
                enddo
             enddo
          enddo
          if (ipr >= 50) then
             fmt = '(1x,'' new rmt'',1x,6f9.5:/(10x,6f9.5))'
             write(fmt(19:19),'(I1)') npcol
             write(fmt(31:31),'(I1)') npcol
             write(stdo,fmt) (wsr(is),is=1,nspec)
          endif
          volfac = ratio
          return
       endif
    enddo
    call rx('sclws2: this cannot happen')

300 format(6x,a,3f13.7)
  end subroutine sclws2

  subroutine pairs(nbas,nbasp,alat,plat,rmax,baspp,ipsp, & !nd,iltab, &
       pltab,nttab,iv_a_ontab,iv_a_oiax,mxcsiz)
    !- Allocate memory for and create neighbor table for a crystal
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis (input)
    !i   nbasp :size of padded basis (layer programs)
    !i          nbasp = nbas !+ nbas(left bulk) + nbas(right bulk)
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   plat  :primitive lattice vectors, in units of alat (input)
    !i   rmax  :maximum range for connecting vector, in a.u.
    !i          All connecting vectors with length < rmax(i)+rmax(j)
    !i          are retained.  rmax may be a scalar, a species-dependent
    !i          array, or a site-dependent array, depending on ipsp(1);
    !i          see description of ipsp
    !i   baspp :basis vectors, doubly padded for planar geometry
    !i   ipsp  :index to which species each site belongs, for padded basis;
    !i          identifies which rmax is associated with each site. NB:
    !i          ipsp(1) = -1 => rmax is a global scalar, independent of site
    !i          ipsp(1) =  0 => rmax is site-, not species-dependent.
    !i          In either of these cases, ipsp is not used.
    !i   nd    :number of dimensions for which periodic boundary conditions
    !i          are used
    !i   iltab :iltab<0, has no effect.  Otherwise, see pltabp.
    !i   pltabp:include only pairs for which pltabp(jb)-pltabp(ib) <= iltab
    !i   mxcsiz:if nonzero, use in place of internal formula for mxnbr
    !o Outputs
    !o   nttab   :total number of pairs in neighbor table and iax
    !o   w(ontab):ntab array; see pairc below where it is generated
    !o   w(oiax) :iax array; see pairc, below where it is generated
    !o   mxcsiz  :size of the largest cluster encountered
    ! ----------------------------------------------------------------------
    implicit none
    integer:: nbas , nbasp , ipsp(1) , pltab(1) , nttab , nd , iltab=-1,i_data_size
    integer,allocatable :: iv_a_ontab(:)
    integer,allocatable :: iv_a_oiax(:,:)
    double precision :: alat,plat(9),rmax(1),baspp(3,1)
    double precision :: avw,avwsr,vol
    integer:: modep(3) , nbaspp , mxcsiz , mxnbr , i , niax , isw
    real(8) ,allocatable :: wk_rv(:)
    integer,allocatable:: iv_a_tmp(:,:)
    parameter (niax=10)
    ! ... Set up input for call to pairc
    nbaspp = nbas !2*nbasp - nbas
    ! ... Estimate an upper bound to the size of the neighbor table
    avw = avwsr(plat,alat,vol,nbas)
    mxnbr = 3*(2*rmax(1)/avw)**3*nbasp
    if (mxcsiz > 0) mxnbr = mxcsiz
    if (allocated(iv_a_ontab)) deallocate(iv_a_ontab)
    allocate(iv_a_ontab(abs(nbasp+1)))
    if (allocated(iv_a_oiax)) deallocate(iv_a_oiax)
    allocate(iv_a_oiax(niax,mxnbr))
    iv_a_oiax=0
    allocate(wk_rv(3*mxnbr))
    modep=2
!    do  10  i = 1, 3
!       modep(i) = 2
!       if (i > nd) modep(i) = 0
!10  END DO
    ! ... This makes the neighbor table
    nttab = mxnbr
    isw = 0
    call pairc ( 1 , nbasp , nbaspp , modep , isw , ipsp , alat , &
         plat , baspp , baspp , rmax , iltab , pltab , nttab , iv_a_ontab &
         , iv_a_oiax , wk_rv , mxcsiz )
    ! ... Allocate iax to proper size
    i_data_size=size(iv_a_oiax)
    allocate(iv_a_tmp(niax,i_data_size/niax))
    iv_a_tmp=iv_a_oiax
    deallocate(iv_a_oiax)
    i_data_size=min(i_data_size,niax*nttab)
    allocate(iv_a_oiax(niax,nttab))
    iv_a_oiax(:,:i_data_size/niax)=iv_a_tmp(:,:i_data_size/niax)
    deallocate(iv_a_tmp)
    if (allocated(wk_rv)) deallocate(wk_rv)
  end subroutine pairs

  subroutine pairc(ib1,ib2,nbasp,mode,isw,ips,alat,plat,pos,ctr, &
       range,iltab,pltabp,nttab,ntab,iax,rtab,mxcsiz)
    !- Make a neighbor table (crystal version)
    ! ----------------------------------------------------------------
    !i Inputs:
    !i  ib1,ib2:range of sites for which to generate tables
    !i   nbasp :the size of the basis, plus possible extensions.
    !i          Usually nbasp=nbas, but will differ in special
    !i          cases, such as having padding sites to extend
    !i          to a semi-infinite geometry.
    !i   mode:  vector of length 3 governing how pos shortened (see shorps)
    !i   isw:   1's digit fixes how range is calculated.
    !i           0: vector length must be < range(i)+range(j)
    !i           1: include all connecting vecs w/ r < range(i)
    !i         10's digit sets what part of iax table is not calculated
    !i           1: do not calculate iax(6)
    !i              (may be needed when ctr and pos are different)
    !i           2: calculate only iax(1..5)
    !i   ips   :index to which species each site belongs, for padded basis;
    !i          identifies which rmax is associated with each site. NB:
    !i          ips(1) = -1 => rmax is a global scalar, independent of site
    !i          ips(1) =  0 => rmax is site-, not species-dependent.
    !i          In either of these cases, ips is not used.
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   plat  :primitive lattice vectors, in units of alat (input)
    !i   pos   :site positions (doubly padded for planar geometry)
    !i   ctr   :ctr(1..3,ib) is the effective center of the cluster
    !i          associated with site ib for purposes of computing distance
    !i          pos(jb)-ctr(ib).  May point to the same address space as pos
    !i   range :maximum range for connecting vector, in a.u..
    !i          This quantity may be a scalar, a species-dependent
    !i          array, or a site-dependent array, depending on ips(1);
    !i          see description of ips.  Precisely what meaning range has
    !i          depends on mode and isw.
    !i   iltab :iltab<0, has no effect.  Otherwise, see pltabp.
    !i   pltabp:include only pairs for which pltabp(jb)-pltabp(ib) <= iltab
    !i   nttab :maximum dimension of iax table; used to guard against
    !i          generated table size exceeding dimension of iax.
    !o Outputs:
    !o   nttab    :total number of pairs generated
    !o   iax      :neighbor table containing information about each pair ip
    !o             For each pair ip, information is contained in iax(*,ip).
    !o             as described below.  iax is ordered grouped by the basis
    !o             atoms, so that all pairs connected to site ib are grouped
    !o             together.  For each pair ip, iax(*,ip) contains:
    !o   iax(1)   :site index to basis atoms ib=source;
    !o             all pairs with common ib are contiguous
    !o   iax(2)   :site index to jb=field of each pair
    !o   iax(3..5):multiples of plat added the difference in site positions
    !o             that connect the pair.
    !o   iax(6)   :index to conjugate (jb,ib) pair matching (ib,jb)
    !o             NB: no matching pairs outside (ib1..ib2) can be found.
    !o   iax(7)   :permutation index ordering cluster by increasing
    !o             effective site index; see ppair4.f
    !o   iax(8)   :left untouched by pairc
    !o   iax(9)   :left untouched by pairc
    !o   iax(10)  :effective site index; see siteid.f
    !o   ntab     :ntab(ib)=number of pairs in iax table preceding ib
    !o             ntab is created for ib1:ib2+1.
    !o   rtab     :rtab(1..3,ip) = pos(jb)-ctr(ib) for pair ip
    !o   mxcsiz   :the largest cluster encountered
    !r Remarks
    !r   For each site ib=ib1..ib2, pairc finds all connecting vectors
    !r   for a lattice of points with periodic boundary conditions in
    !r   1, 2, or 3 dimensions, within a specified range of site ib.
    !r   The range can be defined in various ways, depending on isw.
    !u Updates
    !u   23 Apr 02 added option to make only iax(1..5) (isw=20)
    ! ----------------------------------------------------------------
    implicit none
    integer :: ib1,ib2,nbasp,mode(3),isw,nttab,niax,ips(nbasp), &
         ntab(ib1:ib2+1),iltab,pltabp(nbasp),mxcsiz
    parameter (niax=10)
    integer :: iax(niax,1),n,itmp,ixo
    double precision :: alat,plat(3,3),pos(3,nbasp),ctr(3,ib2), &
         range(1),rtab(3,1)
    ! Local variables
    integer:: ib , is , jb , mtab , i , moder , mode2(3) , nlat , &
         mxntab , nsite
    real(8) ,allocatable :: wk1_rv(:)
    real(8) ,allocatable :: wk2_rv(:)
    real(8) ,allocatable :: wk3_rv(:)
    real(8) ,allocatable :: pos_rv(:)
    real(8) ,allocatable :: lat_rv(:)
    real(8) ,allocatable :: ctr_rv(:)

    double precision :: r1,rr,qlat(3,3),p0(3)
    integer ::iwdummy


    !     call tcn('ppair1')

    ! --- Setup ---
    nsite = ib2-ib1+1
    mxntab = nttab
    moder = mod(isw,10)
    do  3  i = 1, 3
       mode2(i) = mode(i)
       if (mode2(i) == 1) mode2(i) = 0
3   END DO
    ! ... Make r1 = 2*maximum range
    r1 = range(1)
    if (ips(1) >= 0) then
       do  5  ib = 1, nbasp
          is = ib
          if (ips(1) > 0) then
             is = ips(ib)
          endif
          r1 = max(r1,range(is))
5      END DO
    endif
    if (moder == 0) r1 = 2*r1
    r1 = 2*r1
    ! ... List of lattice vectors to add to pos(ib)-pos(jb)
    call xlgen ( plat , r1 / alat , 0d0 , 0 , 20 , mode , i , iwdummy )
    allocate(lat_rv(3*i))
    call xlgen ( plat , r1 / alat , 0d0 , i , 0 , mode , nlat , lat_rv)
    ! ... qlat = (plat^-1)^T so that qlat^T . plat = 1
    call mkqlat(plat,qlat,rr)
    ! ... Save true pos in opos
    !     and ctr in octr in case same address space used for ctr
    allocate(pos_rv(3*nbasp))
    call dpcopy ( pos , pos_rv , 1 , 3 * nbasp , 1d0 )
    allocate(ctr_rv(3*nsite))
    call dpcopy ( ctr , ctr_rv , 1 , 3 * nsite , 1d0 )

    ! --- For each ib, find all pairs for which dr < range ---
    nttab = 1
    ntab(ib1) = 0
    mtab = 1
    mxcsiz = 0
    do  10  ib = ib1, ib2
       r1 = range(1)
       if (ips(1) >= 0) then
          is = ib
          if (ips(1) > 0) then
             is = ips(ib)
          endif
          r1 = range(is)
       endif
       ! --- Shorten all pos relative to ctr(ib) ---
       ! ... Make pos-ctr(ib)
       call dpcopy ( pos_rv , pos , 1 , 3 * nbasp , 1d0 )
       call dpcopy ( ctr_rv , ctr , 1 , 3 * nsite , 1d0 )
       do  12  i = 1, 3
          p0(i)  = ctr(i,ib)
          do  14  jb = 1, nbasp
             pos(i,jb) = pos(i,jb) - p0(i)
14        END DO
12     END DO
       ! ... Shorten pos-ctr(ib)
       call shorps(nbasp,plat,mode2,pos,pos)
       ! ... Undo shift -ctr(ib) to restore shortened pos to absolute pos
       do   jb = 1, nbasp
          do   i = 1, 3
             pos(i,jb) = pos(i,jb) + p0(i)
          enddo
       enddo
       ! --- Find all sites in range of ctr ---
       call ppair2 ( nbasp , iltab , pltabp , moder , alat , qlat , &
            pos , p0 , range , ips , rtab , ib , r1 , nlat , lat_rv , &
            pos_rv , mxntab , nttab , iax )
       ! --- Sort table by increasing length ---
       call ppair3 ( nttab - mtab , iax(1, mtab ), rtab(1 , mtab))! , wk1_rv , wk2_rv , wk3_rv )
       ! --- Cleanup for this ib ---
       mtab = nttab
       ntab(ib+1) = nttab-1
       mxcsiz = max(mxcsiz,ntab(ib+1)-ntab(ib))
10  END DO
    nttab = nttab-1
    ! --- Restore original pos,ctr ---
    call dpcopy ( pos_rv , pos , 1 , 3 * nbasp , 1d0 )
    call dpcopy ( ctr_rv , ctr , 1 , 3 * nsite , 1d0 )
    if (allocated(ctr_rv)) deallocate(ctr_rv)
    if (allocated(pos_rv)) deallocate(pos_rv)


    ! --- Fill out iax table ---
    call ppair1(isw,ib1,ib2,nbasp,ips,alat,plat,pos,range, &
         nttab,ntab,iax,mxcsiz)

    if (allocated(lat_rv)) deallocate(lat_rv)

  end subroutine pairc

  subroutine ovmin(sovmin,nbas,nbasp,alat,plat,rmax,rmt,dclabl, &
       ips,mode,z,iv_a_ontab,iv_a_oiax,pos,iprtbl)
    use m_gradzr,only:gradzr
    !- Check volume and sphere overlaps, optionally minimizing them
    ! ----------------------------------------------------------------
    !i Inputs
    !i   sovmin: a set of modifiers, with the syntax
    !i          -mino[:dxmx=#][:xtol=#][:style=#]:site-list
    !i   nbas  :size of basis
    !i   nbasp :size of padded basis (layer programs)
    !i          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   plat  :primitive lattice vectors, in units of alat
    !i   rmax  :potential radius, in a.u.
    !i   rmt   :augmentation radius, in a.u.
    !i   dclabl:class name, packed as a real number
    !i   ips   :species table: site ib belongs to species ips(ib)
    !i   mode:  vector of length 3 governing how pos shortened (see shorps)
    !i   z     :nuclear charge
    !i   w(ontab):ntab(ib)=# pairs in iax table preceding ib (pairc.f)
    !i   w(oiax):neighbor table containing pair information (pairc.f)
    !i   pos   :basis vectors
    !i   iprtbl: nonzero if to call ovlchk and print table of overlaps
    !o Outputs
    !o   Sphere overlaps are printed out
    !r Remarks
    !r   rmt(1)  not used now
    !u Updates
    !u   22 Oct 02  weight ES-ES and atom-ES overlaps differently when
    !u              minimizing sphere overlap positions
    !u    9 Dec 98  replace call to frpmin with call to gradzr.
    !u    8 Sep 98  small patches in minimizing algorithm
    !u   24 Nov 97  changed ovmin to call fovlp for fast execution
    ! ----------------------------------------------------------------
    !     implicit none
    integer:: nbas , nbasp , iprtbl
    integer,allocatable :: iv_a_ontab(:)
    integer,allocatable :: iv_a_oiax(:,:)
    double precision :: plat(3,3),pos(3,nbasp),rmax(1),rmt(1),z(1),alat
    character(8):: dclabl(*) !double precision
    integer :: ips(1),mode(3)
    character sovmin*(*)
    double precision :: alato,plato(9),xx
    integer:: nbaso , nbaspo , mxlst , nlst , modeo(3) , novl
    parameter (mxlst=256)
    integer :: ilst(mxlst)
    real(8) ::wdummy(1,1)
    double precision :: fovl,xtol,gtol,dxmn,dxmx,fovmx
    double precision :: wk(0:27)
    integer :: i1mach,isw,ir,i,j,j1,j2,ls,m,lstyle, &! & op
         iv,parg,nlstc,mxint,nclass,ib,ic,iclbsj,maxit,ipr,n!,lgunit
    character dc*1
    integer,allocatable:: olist(:)
    real(8),allocatable:: w_opos(:,:),w_oz(:),w_ormax(:),w_oips(:),w_op(:)
    ! --- Print out positions and overlaps ---
    call getpr(ipr)
    if (iprtbl > 0) call ovlchk(nbas,nbasp,pos,alat,rmax,[0d0], &
         dclabl,ips,mode,plat,fovmx,xx)
    call fovlp ( 1 , nbas , iv_a_ontab , iv_a_oiax , plat , pos , &
         ips , alat , rmax , z , 6d0 , 1d0 , .75d0 , .5d0 , fovmx , fovl &
         , novl )
    if (novl == 0) novl = 1
    if (ipr >= 10 .OR. iprtbl > 0) &
         write(stdo,"(' OVMIN:  fovl= ',f5.1,' <ovlp> = ',f5.1,'%', &
         '   max ovlp = ',f5.1,'%')") fovl/novl,(fovl/novl)**(1/6d0)*100,fovmx*100
    !     --- Minimize overlaps wrt positions in list ---
    if (sovmin /= ' ') call rx('ovmin: need to recover if necessary.takao See old version')
  end subroutine ovmin
  subroutine ovlchk(nbas,nbasp,pos,alat,rmax,rmt,dclabl, &
       ips,mode,plat,fovl,volsph)
    use m_lmfinit,only: ispec!ssite=>v_ssite
    use m_ftox
    !- Check volume and sphere overlaps
    ! ----------------------------------------------------------------
    !i Inputs
    !i   nbas  :size of basis
    !i   nbasp :size of padded basis (layer programs)
    !i          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
    !i   pos   :basis vectors
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   rmax  :augmentation radius, in a.u.,
    !i   rmt   :augmentation radius, in a.u., by species
    !i   dclabl:class name, packed as a real number
    !i   ips   :species table: site ib belongs to species ips(ib)
    !i   mode  :same as mode in pairc.f
    !i   rmax,rmt (see remarks)
    !o Outputs
    !o   fovl, a function of the positive overlaps, for now set to
    !o         sum (ovl)^6
    !o   volsph sum of sphere volumes
    !r Remarks
    !r   Checks overlaps for both rmax and rmt.
    !r   rmt(1) <= 0 switches off part with rmt.
    !u Updates
    !u   21 Aug 02 Can print out positions as multiples of plat
    ! ----------------------------------------------------------------
    implicit none
    ! Passed parameters
    integer :: nbas,nbasp
    double precision :: plat(3,3),pos(3,nbasp),rmax(1),rmt(1), &
         alat,fovl,volsph
    character(8):: dclabl(*)
    integer :: ips(1),mode(3),is
    ! Local parameters
    double precision :: dd(3),dd1,dd2,sumrs,summt,ovlprs,ovlpmt, &
         ctrs,ctmt,ddot,dlat(3),xx,avwsr,vol
    double precision :: qlat(3,3),volspp
    integer :: ibas,jbas,ic,jc,kc,m,ipr,i1mach,m1,m2,m3,isw,istdo
    character(80) :: a, ch(1)
    logical :: lterse,cmdopt,lrmt
    character(8) :: clabl,clablj
    integer:: ifile_handle,ifp,js
    character(10):: i2char
    call getpr(ipr)
    lrmt = rmt(1) .gt. 0
    !      stdo = lgunit(1)
    call mkqlat(plat,qlat,xx)

    ! --- Determine which linear combination of plat is shortest ---
    dd1 = 9d9
    call dpzero(dlat,3)
    do  102  m1 = 1, 3
       do  101  m2 = 1, 3
          do  10  m3 = 1, 3
             if (mode(m1)*mode(m2)*mode(m3) == 0) goto 10
             do  122  ic = -1, 1
                do  121  jc = -1, 1
                   do  12  kc = -1, 1
                      dd=plat(:,m1)*dble(ic)+plat(:,m2)*dble(jc)+plat(:,m3)*dble(kc)
                      !          call dpcopy(plat(1,m1),dd,1,3,dble(ic))
                      !          call dpadd(dd,plat(1,m2),1,3,dble(jc))
                      !          call dpadd(dd,plat(1,m3),1,3,dble(kc))
                      dd2 = ddot(3,dd,1,dd,1)
                      if (dd1 > dd2 .AND. dd2 > 1d-6) then
                         call dpcopy(dd,dlat,1,3,1d0)
                         dd1 = dd2
                      endif
12                 END DO
121             END DO
122          END DO
10        END DO
101    END DO
102 END DO

    if(ipr>=10) then
       if(lrmt) write(stdo,*)'    Site   ic Spec        Rmax   RMT   Position'
       if(.NOT.lrmt)write(stdo,*)'    Site   ic Spec        Rmax      Position'
    endif
    volsph = 0d0
    ifp=ifile_handle()
    open(ifp,file='SiteInfo.lmchk')

    do  20  ibas = 1, nbasp
       ic = ips(ibas) !class id
       is = ispec(ibas) !ssite(ibas)%spec
       if (ipr <= 10) goto 20
       if (ibas == nbas+1) write(stdo,'(''  ... Padding basis'')')
       !        call r8tos8(dclabl(ic),clabl)
       clabl=dclabl(is)
       if (dclabl(1) == '') clabl= i2char(is)
       if (lrmt) then
          write(stdo,450) ibas,ic,clabl,rmax(ic),rmt(ic),(pos(m,ibas),m=1,3)
          write(ifp,450)  ibas,ic,clabl,rmax(ic),rmt(ic), (pos(m,ibas),m=1,3)
450       format('conf ',i5,3x,i4,2x,a8,2f12.6,3f11.5)
       else
          write(stdo,351) ibas,ic,clabl,rmax(ic),(pos(m,ibas),m=1,3)
          write(ifp,351)ibas,ic,clabl,rmax(ic),(pos(m,ibas),m=1,3)
351       format('conf ',i5,3x,i4,2x,a8,f8.6,2x,3f11.6)
          if (ipr >= 41) then
             call dgemm('T','N',3,1,3,1d0,qlat,3,pos(1,ibas),3,0d0,dd,3)
             !          do  m = 1, 3
             !            dd(m) = pos(1,ibas)*qlat(1,m)+pos(2,ibas)*qlat(2,m)+
             !     .        pos(3,ibas)*qlat(3,m)
             !          enddo
             write(stdo,352) dd
352          format(13x,'as multiples of plat:',3f11.6)
          endif
       endif

       if (ibas <= nbas) volsph = volsph + 4.188790205d0*rmax(ic)**3
       !        volspp = volspp + 4.188790205d0*rmax(ic)**3
20  END DO
    close(ifp)

    xx = avwsr(plat,alat,vol,nbas)
    if (ipr >= 10) then ! .AND. volspp == volsph) then
       write(stdo,ftox)' Cell volume= ',ftof(vol,3),'Sum of sphere volumes=', &
            ftof(volsph,3),'(',ftof(volsph/vol,1),' %)'
    endif
    !$$$      elseif (ipr .ge. 10) then
    !$$$        volspp = 2*volspp-volsph
    !$$$        call info5(0,1,0,
    !$$$     .  ' Cell volume= %1,5;5d'//
    !$$$     .  ' Sum of sphere volumes= %1,5;5d + %1,5;5d(2 x pad) '//
    !$$$     .  ' ratio=%1;5d',vol,volsph,volspp-volsph,volspp/vol,0)
    !$$$      endif

    ! --- Check sphere overlaps ---
    fovl = 0
    lterse = cmdopt('-terse',6,0,a) .or. cmdopt('--terse',7,0,a)
    if (lrmt       .AND. ipr > 10) write(stdo,453)
    if ( .NOT. lrmt .AND. ipr > 10) write(stdo,463)
    do  301  ibas = 1, nbasp
       ic = ips(ibas)
       is = ispec(ibas) !ssite(ibas)%spec
       if (ipr >= 10) then
          !          call r8tos8(dclabl(ic),clabl)
          clabl = trim(dclabl(is))//i2char(ic)
          if (dclabl(1) == '') clabl=i2char(ic) 
          !     all awrit1('%x%,4i',clabl,8,0,ic)
       endif
       do  30  jbas = ibas, nbasp
          jc = ips(jbas)
          js = ispec(jbas) !ssite(jbas)%spec
          if (rmax(ic) == 0 .AND. rmax(jc) == 0) then
             goto 30
          endif
          if (ipr >= 10) then
             !          call r8tos8(dclabl(jc),clablj)
             clablj=trim(dclabl(js))//i2char(jc) !jc)
             if (dclabl(1) == '')clablj=i2char(jc)
             !     all awrit1('%x%,4i',clablj,8,0,jc)
          endif
          if (ibas == jbas) then
             if (ddot(3,dlat,1,dlat,1) == 0) goto 30
             call dcopy(3,dlat,1,dd,1)
          else
             do  33  m = 1, 3
                dd(m) = pos(m,jbas) - pos(m,ibas)
33           END DO
             dd1 = dsqrt(dd(1)**2 + dd(2)**2 + dd(3)**2)
             call shorps(1,plat,mode,dd,dd)
             dd2 = dsqrt(dd(1)**2 + dd(2)**2 + dd(3)**2)
             if (dd2 > dd1+1d-9) call rx('bug in ovlchk')
             ! ...     test against shorbz
             !         call mkqlat(plat,qlat,xx)
             !         call shorbz(dd,dd,plat,qlat)
             !         print *, dd2,dsqrt(dd(1)**2 + dd(2)**2 + dd(3)**2)
          endif
          do  34  m = 1, 3
             dd(m) = dd(m)*alat
34        END DO
          dd1 = dsqrt(dd(1)**2 + dd(2)**2 + dd(3)**2)
          sumrs = rmax(ic) + rmax(jc)
          ovlprs = sumrs - dd1
          if (dd1 < 1d-6) then
             write(stdo,451) ibas,jbas,clabl,clablj,dd,dd1
             call rx('ovlchk: positions coincide')
          endif
          ctrs = nint(1000*ovlprs/dd1)/10d0
          if (lrmt) then
             summt = rmt(ic) + rmt(jc)
             ovlpmt = summt - dd1
             ctmt = nint(1000*ovlpmt/dd1)/10d0
          endif
          fovl = fovl + max(ovlprs/dd1,0d0)**6
          if ((lterse .OR. ipr <= 40) .AND. ctrs <= -10 &
               .OR. ipr <= 10) goto 30
          ch = ' '
          if (ovlprs >= 0d0) ch='*'
          if (lrmt .AND. .FALSE. ) then
             write(stdo,451) ibas,jbas,clabl,clablj,dd,dd1, &
                  sumrs,ovlprs,ctrs,summt,ovlpmt,ctmt,ch
451          format(2i3,2x,a4,1x,a4,3f6.2,2f7.2,f7.2,f5.1,f7.2,f7.2,f5.1,a1)
          else
             write(stdo,461) ibas,jbas,clabl,clablj,dd,dd1, &
                  sumrs,ovlprs,ctrs,ch
461          format(2i3,2x,a8,a8,3f7.3,2f7.3,f7.2,f6.1,a1)
          endif
453       format(/' ib jb',2x,'sp&ic1  sp&ic2',7x,' Pos(jb)-Pos(ib)', &
               4x,'Dist   sumrs   Ovlp   %  summt   Ovlp   %')
463       format(/' ib jb',2x,'sp&ic1  sp&ic2',8x,'Pos(jb)-Pos(ib)', &
               3x,'Dist  sumrs   Ovlp    %')
30     END DO
301 END DO

    !      if (ipr .gt. 0)
    !     .  call awrit1('%N ovlchk: fovl= %;6g',' ',80,i1mach(2),fovl)

  end subroutine ovlchk


  subroutine sumsro(rp,np,ips,a,b,rho,nttab,iax,rpos,rosum)
    !- Add a r.s. superposition of spherical densities at a set of points
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   rp    :set of np points
    !i   np    :number points
    !i   ips   :species table: density for site ib site ib found in ips(ib)
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i         :a(i) is coefficient for species i
    !i   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i         :b(i) is coefficient for species i
    !i   rho   :spherical density for each species
    !i   nttab :total number of pairs in neighbor and iax (pairc.f)
    !i   iax   :neighbor table containing pair information connecting
    !i         :site i to its neighbors
    !i   rpos  :positions of neighbors.
    !i         :rp and rpos must be defined with respect to same origin.
    !o Outputs
    !o   rosum :sum of densities at each of np points
    !l Local variables
    !u Updates
    !u   21 Apr 02 First created
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: np,ips(1),nrmx,niax
    parameter (nrmx=1501, niax=10)
    integer :: nttab,iax(niax,nttab)
    double precision :: rp(3,np),rpos(3,nttab),rho(nrmx,*),rosum(np)
    double precision :: a(*),b(*)
    ! ... Local parameters
    integer :: k,jb,js,ip,jx
    double precision :: dx(3),d,xx,y,e
    double precision :: rp3p,rp5p,rpp3p,rpp32,rpp5p,rppp
    !     double precision rofii(nrmx),dy

    !     x axis for polynomial interpolation.  Use log mesh
    !     Needed only for call to polint
    !      do  k = 1, nrmx
    !        rofii(k) = k
    !      enddo

    ! ... Sum the potential from nttab neighbors
    do  ip = 1, np
       rosum(ip)= 0
       do  k = 1, nttab
          dx(1) = rp(1,ip)-rpos(1,k)
          dx(2) = rp(2,ip)-rpos(2,k)
          dx(3) = rp(3,ip)-rpos(3,k)
          jb = iax(2,k)
          js = ips(jb)
          d = dsqrt(dx(1)**2+dx(2)**2+dx(3)**2)

          xx = 1 + dlog(d/b(js)+1)/a(js)
          if (int(xx) <= nrmx-2) then
             !           General polynomial interpolation.  Works, but slow.
             !           jx = int(xx)-1
             !           call polint(rofii,rho(1,js),nrmx,5,xx,0d0,0,jx,y,dy)

             !           Interpolate in-line by 5-point formula
             jx = nint(xx)-2
             !           first derivative, 3- and 5-point estimate
             rp3p = (rho(jx+3,js)-rho(jx+1,js))/2
             rp5p = (4*rp3p-(rho(jx+4,js)-rho(jx,js))/4)/(4-1)

             !           Second derivative, nearest points and 2nd points
             rpp3p =  rho(jx+1,js) + rho(jx+3,js) - 2*rho(jx+2,js)
             rpp32 = (rho(jx+0,js) + rho(jx+4,js) - 2*rho(jx+2,js))/4
             !           five-point estimate for second derivative
             rpp5p = (rpp3p*4 - rpp32*1) / (4 - 1)

             !           Estimate for third derivative
             rppp = -6*(rp5p - rp3p)

             !           Interpolate function value by polynomial
             e = xx-jx-2
             y = rho(jx+2,js)+e*rp5p+e*e/2*rpp5p+e**3/6*rppp

             rosum(ip) = rosum(ip) + y
          endif

          !       Verbose debugging printout
          !       print 333, k, iax(2,k), js, ip, dx, d, y, rosum(ip)
          ! 333   format(4i4,6f8.3)
       enddo
       !     Debugging printout
       !     print 334, ip, rosum(ip)
       !  334 format(i4,f12.6)
    enddo

  end subroutine sumsro


  subroutine mkqlat(plat,qlat,vol0)
    !- Reciprocal of a lattice vector
    !     implicit none
    double precision :: plat(3,3),qlat(3,3),tripl,vol0
    integer :: i,k
    call cross(plat(1,2),plat(1,3),qlat)
    call cross(plat(1,3),plat(1,1),qlat(1,2))
    call cross(plat(1,1),plat(1,2),qlat(1,3))
    vol0 = tripl(plat,plat(1,2),plat(1,3))
    do  i=1,3
       do  k=1,3
          qlat(i,k)=qlat(i,k)/vol0
       enddo
    enddo
  end subroutine mkqlat

  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine shorps(nbas,plat,mode,pin,pout)
    !- Shift each basis vector by multiples of plat according to mode
    ! ----------------------------------------------------------------
    !i Inputs:  nbas,plat
    !i   nbas<0: sign used as a switch to make shorps return in
    !i         array pout the change in pin, in units of plat
    !i   pin:  position (basis) vectors
    !i   mode  vector of length 3 governing shifts along selected axes.
    !i         0 suppresses shifts along plat(j)
    !i         1 shifts to unit cell at origin (pos in 1st quadrant)
    !i         2 shifts to minimize length of pos
    !o Outputs:
    !o   pout  (may point to the same address space as pin).
    !o   iat:  multiples of plat added
    !r Remarks
    !r   pos = f . plat, with f = integer + fraction for each plat.
    !r   Integer part according to mode.
    !u Updates
    !u   09 Jan 09 Do not change basis vectors if they exactly fall
    !u             on the cell boundary
    !u   10 Apr 02 Patch to handle mixed boundary conditions
    ! ----------------------------------------------------------------
    !     implicit none
    integer :: nbas,mode(3)
    double precision :: plat(3,3),pin(3,nbas),pout(3,nbas)
    ! Local
    double precision :: qlat(3,3),xx,x0,x(3),a2,ap,p0(3),xmin(3),amin,vol
    integer :: ib,i,m,j1,j2,j3,nbas1,j1max,j2max,j3max
    double precision :: tol
    parameter (tol = 1d-12)
    nbas1 = iabs(nbas)
    ! ... qlat = (plat^-1)^T so that qlat^T . plat = 1
    call dinv33(plat,1,qlat,vol)
    do  10  ib = 1, nbas1
       call dpcopy(pin(1,ib),p0,1,3,1d0)
       !   --- Reduce to unit cell centered at or near origin ---
       do  12  i = 1, 3
          !   ... x0 is projection of pin along plat(i)
          x0 = pin(1,ib)*qlat(1,i)+pin(2,ib)*qlat(2,i)+pin(3,ib)*qlat(3,i)
          if (mode(i) <= 0) then
             x(i) = x0
          else
             !   ... leave basis vectors intact if |pin(i,ib)| = 0.5
             if (dabs(x0) <= 0.5d0+tol) then
                xx = 0d0
             else
                xx = idnint(x0)
             endif
             !     ... first octant for mode=1
             if (mode(i) == 1 .AND. x0-xx < -tol) xx = xx-1
             x(i) = x0-xx
          endif
12     END DO
       do  14  m = 1, 3
          pout(m,ib) = x(1)*plat(m,1) + x(2)*plat(m,2) + x(3)*plat(m,3)
14     END DO
       !   --- Try shortening by adding +/- lattice vectors ---
       j1max = 1
       if (mode(1) <= 1) j1max = 0
       j2max = 1
       if (mode(2) <= 1) j2max = 0
       j3max = 1
       if (mode(3) <= 1) j3max = 0
15     continue
       amin = 0
       do  162  j1 = -j1max, j1max
          do  161  j2 = -j2max, j2max
             do  16  j3 = -j3max, j3max
                !     ... (-1,0,1) (plat(1) + (-1,0,1) plat(2)) + (-1,0,1) plat(3))
                do  17  i = 1, 3
                   x(i) = plat(i,1)*j1 + plat(i,2)*j2 + plat(i,3)*j3
17              END DO
                a2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
                ap = pout(1,ib)*x(1) + pout(2,ib)*x(2) + pout(3,ib)*x(3)
                if (a2+2*ap < amin) then
                   xmin(1) = x(1)
                   xmin(2) = x(2)
                   xmin(3) = x(3)
                   amin = a2+2*ap
                endif
16           END DO
161       END DO
162    END DO
       !       if (amin .lt. 0) then
       if (amin < -tol) then
          pout(1,ib) = pout(1,ib) + xmin(1)
          pout(2,ib) = pout(2,ib) + xmin(2)
          pout(3,ib) = pout(3,ib) + xmin(3)
          !         In cases w/ mixed boundary conditions, (-1,0,1) may not be enough
          !         Patched 10 Apr 02
          if (mode(1) == 0 .OR. mode(2) == 0 .OR. mode(3) == 0) goto 15
       endif
       !   --- pout <- pout - pin, units of plat ---
       if (nbas < 0) then
          !           call dpadd(p0,pout(1,ib),1,3,-1d0)
          p0=p0-pout(:,ib)
          do  20  i = 1, 3
             xx = -p0(1)*qlat(1,i) - p0(2)*qlat(2,i) - p0(3)*qlat(3,i)
             if (dabs(xx-nint(xx)) > 1d-10) call rx('bug in shorps')
             pout(i,ib) = xx
20        END DO
       endif
10  END DO
  end subroutine shorps
  !$$$#if TESTforshorps
  !$$$      subroutine fmain
  !$$$C to see that the change in position vectors are multiples of
  !$$$C the lattice vector, copy input,output pos to 'pos','posf'; invoke
  !$$$C mc posf pos -- -t plat -t -i -x
  !$$$      implicit none
  !$$$      integer ix(3),i
  !$$$      double precision pos(3,48),pos2(3,48),plat(3,3),qlat(9),xx
  !$$$
  !$$$      call wkinit(10000)
  !$$$
  !$$$      data plat /
  !$$$     .0.5d0,          .5d0, 0d0,
  !$$$     .0.0d0,          0.d0, 1d0,
  !$$$     .2.570990255d0, -2.570990255d0, 0d0/
  !$$$      data pos /
  !$$$     .-0.697107d0,  1.197107d0,  0.250000d0,
  !$$$     .-0.697107d0,  1.197107d0,  0.750000d0,
  !$$$     .-0.770330d0,  0.770330d0,  0.000000d0,
  !$$$     .-0.770330d0,  0.770330d0,  0.500000d0,
  !$$$     .-0.343553d0,  0.843553d0,  0.250000d0,
  !$$$     .-0.343553d0,  0.843553d0,  0.750000d0,
  !$$$     .-0.416777d0,  0.416777d0,  0.000000d0,
  !$$$     .-0.416777d0,  0.416777d0,  0.500000d0,
  !$$$     .0.010000d0,  0.490000d0,  0.250000d0,
  !$$$     .0.010000d0,  0.490000d0,  0.750000d0,
  !$$$     .0.250000d0,  0.250000d0,  0.500000d0,
  !$$$     .0.500000d0,  0.500000d0,  0.750000d0,
  !$$$     .0.750000d0,  0.750000d0,  1.000000d0,
  !$$$     .1.000000d0,  1.000000d0,  1.250000d0,
  !$$$     .0.250000d0, -0.250000d0,  0.000000d0,
  !$$$     .0.500000d0,  0.000000d0,  0.250000d0,
  !$$$     .0.750000d0,  0.250000d0,  0.500000d0,
  !$$$     .1.000000d0,  0.500000d0,  0.750000d0,
  !$$$     .0.750000d0, -0.250000d0,  0.500000d0,
  !$$$     .1.000000d0,  0.000000d0,  0.750000d0,
  !$$$     .1.250000d0,  0.250000d0,  1.000000d0,
  !$$$     .1.500000d0,  0.500000d0,  1.250000d0,
  !$$$     .0.740000d0, -0.740000d0,  0.000000d0,
  !$$$     .0.740000d0, -0.740000d0,  0.500000d0,
  !$$$     .1.166777d0, -0.666777d0,  0.250000d0,
  !$$$     .1.166777d0, -0.666777d0,  0.750000d0,
  !$$$     .1.093553d0, -1.093553d0,  0.000000d0,
  !$$$     .1.093553d0, -1.093553d0,  0.500000d0,
  !$$$     .1.520330d0, -1.020330d0,  0.250000d0,
  !$$$     .1.520330d0, -1.020330d0,  0.750000d0,
  !$$$     .1.447107d0, -1.447107d0,  0.000000d0,
  !$$$     .1.447107d0, -1.447107d0,  0.500000d0,
  !$$$     .-1.050660d0,  1.550660d0,  0.250000d0,
  !$$$     .-1.050660d0,  1.550660d0,  0.750000d0,
  !$$$     .-1.123883d0,  1.123883d0,  0.000000d0,
  !$$$     .-1.123883d0,  1.123883d0,  0.500000d0,
  !$$$     .1.873883d0, -1.373883d0,  0.250000d0,
  !$$$     .1.873883d0, -1.373883d0,  0.750000d0,
  !$$$     .1.800660d0, -1.800660d0,  0.000000d0,
  !$$$     .1.800660d0, -1.800660d0,  0.500000d0,
  !$$$     .-1.404214d0,  1.904214d0,  0.250000d0,
  !$$$     .-1.404214d0,  1.904214d0,  0.750000d0,
  !$$$     .-1.477437d0,  1.477437d0,  0.000000d0,
  !$$$     .-1.477437d0,  1.477437d0,  0.500000d0,
  !$$$     .2.227437d0, -1.727437d0,  0.250000d0,
  !$$$     .2.227437d0, -1.727437d0,  0.750000d0,
  !$$$     .2.154214d0, -2.154214d0,  0.000000d0,
  !$$$     .2.154214d0, -2.154214d0,  0.500000d0/
  !$$$
  !$$$
  !$$$      call prmx('plat',plat,3,3,3)
  !$$$      call prmx('starting pos',pos,3,3,48)
  !$$$      ix(1) = 2
  !$$$      ix(2) = 2
  !$$$      ix(3) = 2
  !$$$      call shorps(48,plat,ix,pos,pos2)
  !$$$      call prmx('final pos',pos2,3,3,48)
  !$$$
  !$$$      call mkqlat(plat,qlat,xx)
  !$$$      do  10  i = 1, 48
  !$$$        call shorbz(pos(1,i),pos2(1,i),plat,qlat)
  !$$$   10 continue
  !$$$
  !$$$      call prmx('from shorbz',pos2,3,3,48)
  !$$$      end
  !$$$#endif
  !$$$#if TEST2
  !$$$      subroutine fmain
  !$$$C Check special case in which a bug fix, mixed boundary conditions
  !$$$      implicit none
  !$$$      integer ix(3),i
  !$$$      double precision pos(3,1),pos2(3,1),plat(3,3),qlat(9),xx
  !$$$      double precision dd1,dd2
  !$$$
  !$$$      call wkinit(10000)
  !$$$
  !$$$      data plat /-0.5d0,0.5d0,0.0d0,
  !$$$     .0.0d0,0.0d0,1.0d0,
  !$$$     .7.0d0,7.0d0,4.0d0/
  !$$$      data ix /2,2,0/
  !$$$
  !$$$      pos(1,1) = 2.5d0
  !$$$      pos(2,1) = 3.0d0
  !$$$      pos(3,1) = 0.0d0
  !$$$
  !$$$      dd1 = dsqrt(pos(1,1)**2 + pos(2,1)**2 + pos(3,1)**2)
  !$$$      call shorps(1,plat,ix,pos,pos2)
  !$$$      dd2 = dsqrt(pos2(1,1)**2 + pos2(2,1)**2 + pos2(3,1)**2)
  !$$$      print *, dd1, dd2
  !$$$
  !$$$      call mkqlat(plat,qlat,xx)
  !$$$      do  10  i = 1, 1
  !$$$        call shorbz(pos(1,i),pos2(1,i),plat,qlat)
  !$$$   10 continue
  !$$$      dd2 = dsqrt(pos2(1,1)**2 + pos2(2,1)**2 + pos2(3,1)**2)
  !$$$      print *, dd1, dd2
  !$$$
  !$$$      end
  !$$$#endif
  !$$$
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine fovlp(ib1,ib2,ntab,iax,plat,pos,ipc,alat,rmax,z,pwr, &
       facaa,facae,facee,fmax,f,inc)
    !- Analytical function of the sphere overlaps
    !u Updates
    !u   22 Oct 02  weight ES-ES and atom-ES overlaps differently
    !u              New argument list
    !     implicit none
    integer :: ib1,ib2,ntab(ib2+1),niax,ipc(ib2),inc
    parameter (niax=10)
    integer :: iax(niax,1)
    double precision :: alat,plat(3,3),pos(3,ib2),rmax(1),pwr,fmax,f, &
         facaa,facae,facee
    double precision :: dsqr,dr,d,sumrs,z(ib2),zb,zi,fac
    integer :: i,ib,i0,i1,ii,jj,kk,ix,icb,ici

    fmax = -99d0
    f = 0
    inc = 0
    do  10  ib = ib1, ib2
       i0 = ntab(ib)+1
       i1 = ntab(ib+1)
       do  12  i = i0+1, i1
          ii = iax(3,i)-iax(3,i0)
          jj = iax(4,i)-iax(4,i0)
          kk = iax(5,i)-iax(5,i0)
          dsqr = 0
          do  16  ix = 1, 3
             dr = pos(ix,iax(2,i)) - pos(ix,iax(1,i0)) + &
                  plat(ix,1)*ii + plat(ix,2)*jj + plat(ix,3)*kk
             dsqr = dsqr + dr**2
16        END DO
          icb = ipc(ib)
          ici = ipc(iax(2,i))
          zb = z(icb)
          zi = z(ici)
          fac = 1
          if (zb /= 0 .AND. zi /= 0) then
             fac = facaa
          elseif (zb == 0 .AND. zi == 0) then
             fac = facee
          else
             fac = facae
          endif
          sumrs = rmax(icb) + rmax(ici)
          d = alat*dsqrt(dsqr)
          if (sumrs > d) then
             f = f + fac*(sumrs/d-1)**pwr
             inc = inc+1
             !          else
             !            goto 10
          endif
          fmax = max(sumrs/d-1,fmax)
12     END DO
10  END DO
  end subroutine fovlp


  double precision function volsph(nspec,nrspec,wsr)
    !- Sum-of-sphere volumes
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i  nspec :number of species
    !i  nrspec:number of atoms in the ith species
    !i  wsr   :Wigner-Seitz sphere radius
    !o Outputs:
    !o  volsph:sum-of-sphere volumes, in units of wsr**3
    ! ----------------------------------------------------------------------
    !     implicit none
    ! Passed variables:
    integer :: nspec,nrspec(*)
    double precision :: wsr(*)
    ! Local variables:
    integer :: ic
    double precision :: fpi3
    parameter(fpi3=4.18879020478639053d0)

    volsph = 0d0
    do  ic = 1, nspec
       volsph = volsph + fpi3*nrspec(ic)*wsr(ic)**3
    enddo
  end function volsph

  subroutine ppair1(isw,ib1,ib2,nbasp,ips,alat,plat,pos,range,nttab,ntab,iax,mxcsiz)
    use m_lgunit,only:stdo
    ! C- Fill out parts of the aix table
    ! C ----------------------------------------------------------------
    ! Ci  Inputs
    ! Ci   isw   :1's digit fixes how range is calculated.
    ! Ci           0: vector length must be < range(i)+range(j)
    ! Ci           1: include all connecting vecs w/ r < range(i)
    ! Ci         :10's digit sets what part of iax table is calculated
    ! Ci           0: make iax(6),iax(7),iax(10)
    ! Ci           1: make iax(7),iax(10)
    ! Ci           2: no change to iax: printout only
    ! Ci           4: just make iax(6)
    ! Ci   ib1   :fill out iax table for pairs ntab(ib1)+1..ntab(ib2)
    ! Ci   ib2   :fill out iax table for pairs ntab(ib1)+1..ntab(ib2)
    ! Ci   nbasp :size of padded basis (not needed)
    ! Ci   ips   :species table: site ib belongs to species ips(ib)
    ! Ci   alat  :length scale of lattice and basis vectors, a.u.
    ! Ci   plat  :primitive lattice vectors, in units of alat
    ! Ci   pos   :basis vectors
    ! Ci   range :maximum range for connecting vector, in a.u..
    ! Ci          This quantity may be a scalar, a species-dependent
    ! Ci          array, or a site-dependent array, depending on ips(1);
    ! Ci          see description of ips.  See 1s digit of isw for
    ! Ci          how range is used.
    ! Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
    ! Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
    ! Ci   iax   :neighbor table containing pair information (pairc.f)
    ! Ci   mxcsiz:maximum cluster size (for printout only)
    ! Co  Outputs
    ! Co   iax(6)   :index to conjugate (jb,ib) pair matching (ib,jb)
    ! Co             NB: only matching pairs within site list can be found.
    ! Co   iax(7)   :permutation index ordering cluster by increasing
    ! Co             effective site index; see ppair4.f
    ! Co   iax(10)  :effective site index; see siteid.f
    ! C ----------------------------------------------------------------
    implicit none
    integer isw,ib1,ib2,nbasp,nttab,niax,ips(nbasp),ntab(ib1:ib2)
    parameter (niax=10)
    integer iax(niax,1),mxcsiz
    double precision alat,plat(3,3),pos(3,19),range(1)
    integer:: ib , is , jb , js , ipr , i , j , moder , it , jt ,iprint ,  nsite , isw1
    real(8) ,allocatable :: pos_rv(:)
    integer ,allocatable :: iwk_iv(:)
    double precision r1,r2,rr,rcut,vlat(3),tol
    parameter (tol=1d-5)
    isw1 = mod(isw/10,10)
    ipr = iprint()
    moder = mod(isw,10)
    nsite = ib2-ib1+1
    if (isw1 .eq. 2) goto 80
    if (isw1 .eq. 4) goto 71
    ! --- Set iax(7) to sort this cluster ---
    call ppair5(ib1,ib2,plat,pos,tol,ntab,iax)
    ! --- For each pair, find matching pair, store in iax(6) ---
71  continue
    do  74  it = 1,  nttab
       iax(6,it) = 0
74  enddo
    if (mod(isw1,2) .eq. 0) then
       do 170  ib = ib1, ib2
          do  70  it = ntab(ib)+1, ntab(ib+1)
             if (iax(6,it) .ne. 0) goto 70
             jb = iax(2,it)
             !  ... No matching pair for padded sites
             if (jb .lt. ib1 .or. jb .gt. ib2) goto 70
             do  72  jt = ntab(jb)+1, ntab(jb+1)
                if (iax(2,jt) .eq. ib .and. &
                     iax(3,it) .eq. -iax(3,jt) .and. &
                     iax(4,it) .eq. -iax(4,jt) .and. &
                     iax(5,it) .eq. -iax(5,jt))  then
                   iax(6,it) = jt
                   iax(6,jt) = it
                   goto 73
                endif
72           enddo
             call fexit3(-1,1,' Exit -1 pairc: cannot find pair'// &
                  ' matching sites (%i,%i), pair %i',ib,jb,it-ntab(ib))
73           continue
70        enddo
170    enddo
    endif
    if (isw1 .eq. 4) return
    ! ... Assign a unique id for every different site in the cluster table
    allocate(iwk_iv(nttab))
    allocate(pos_rv(3*nttab))
    call siteid ( iax , nsite , ntab , plat , pos , pos_rv , iwk_iv , i )
    if (allocated(pos_rv)) deallocate(pos_rv)
    if (allocated(iwk_iv)) deallocate(iwk_iv)
    ! --- Printout ---
80  if (ipr .lt. 30) goto 91
    if (ipr .le. 40) write(stdo,'(1x)')
    if (ipr .gt. 40) write(stdo,332)
332 format(/'  ib  jb',9x,'--- r(jb)-r(ib) ---',10x,'d       -x-plat-  map ord  id')
    i = 0
    do  90  it = 1, nttab
       ib = iax(1,it)
       jb = iax(2,it)
       rr = dsqrt(drr2(plat,pos(1,ib),pos(1,jb), iax(3,it),iax(4,it),iax(5,it),vlat))
       r1 = range(1)
       r2 = range(1)
       if (ips(1) .ge. 0) then
          is = ib
          if (ips(1) .gt. 0) then
             is = ips(ib)
          endif
          r1 = range(is)
          js = jb
          if (ips(1) .gt. 0) then
             js = ips(jb)
          endif
          r2 = range(js)
       endif
       if (moder .eq. 0) rcut = r1+r2
       if (moder .eq. 1) rcut = r1
       if (ib .ne. i) then
          if (alat .ne. 1) write(stdo,345) ib,ntab(ib+1)-ntab(ib),rcut/alat,rcut
          if (alat .eq. 1) write(stdo,345) ib,ntab(ib+1)-ntab(ib),rcut
345       format(' pairc, ib=',i3,':',i4,' neighbors in range',2f7.3)
       endif
       i = ib
       if (ipr .gt. 40) write(stdo,334) iax(1,it),iax(2,it),(vlat(j),j=1,3), rr, (iax(j,it), j=3,7),iax(10,it)
334    format(i4,i4,3f11.6,f9.4,3x,3i3,i5,2i4)
90  enddo
91  if (ipr .ge. 20) write(stdo,'('' pairc:'',i8,'' pairs total'',i10,'' is max cluster size'')')nttab, mxcsiz
  end subroutine ppair1


  subroutine ppair2(nbas,iltab,pltabp,moder,alat,qlat,pos,ctr,range, &
       ips,rtab,ib,r1,nlat,lat,trupos,mxntab,nttab,iax)
    !- Kernel of pairc to find all sites in range of ctr
    implicit none
    integer :: nbas,ib,iltab,ips(nbas),pltabp(nbas),niax,nlat,moder, &
         mxntab,nttab
    parameter (niax=10)
    integer :: iax(niax,1)
    double precision :: alat,ctr(3),pos(3,nbas),range(nbas),rtab(3,1)
    double precision :: qlat(3,3),trupos(3,nbas),lat(3,*),r1
    ! Local variables
    integer :: i,ilat,jb,js
    double precision :: r2,rr,rcut,vlat(3),xx,rcutba,dpos(3)

    do  20  jb = 1, nbas
       if (iltab >= 0) then
          if (abs(pltabp(jb)-pltabp(ib)) > iltab) goto 20
       endif
       r2 = range(1)
       if (ips(1) >= 0) then
          js = jb
          if (ips(1) > 0) then
             js = ips(jb)
          endif
          r2 = range(js)
       endif
       if (moder == 0) rcut = r1+r2
       if (moder == 1) rcut = r1
       rcutba = (rcut / alat)**2
       dpos(1) = pos(1,jb)-ctr(1)
       dpos(2) = pos(2,jb)-ctr(2)
       dpos(3) = pos(3,jb)-ctr(3)

       !   --- For each (ib,jb,ilat), do ---
       do  30  ilat = 1, nlat

          if (nttab > mxntab) call rxi( &
               'pairc: table exceeds input maximum size,',mxntab)

          ! ...   Add to list if connecting vector within range
          rtab(1,nttab) = dpos(1) + lat(1,ilat)
          rtab(2,nttab) = dpos(2) + lat(2,ilat)
          rtab(3,nttab) = dpos(3) + lat(3,ilat)
          rr = rtab(1,nttab)**2+rtab(2,nttab)**2+rtab(3,nttab)**2

          !*        call awrit5('try ib,jb,ilat= %i %i %i rr=%;4d: %l',' ',80,
          !*     .    6,ib,jb,ilat,rr,rr.lt.rcut)

          !   --- Add to iax table if this pair in range ---
          if (rr < rcutba) then

             !     ... vlat += shortening vector
             do  32  i = 1, 3
                rtab(i,nttab) = alat*rtab(i,nttab)
                !           rtab(i,nttab) = alat*(rtab(i,nttab)+ctr(i)-pos(i,ib))
                vlat(i) = lat(i,ilat) + pos(i,jb) - trupos(i,jb)
32           END DO

             !     ... iax table for this pair
             iax(1,nttab) = ib
             iax(2,nttab) = jb
             do  33  i = 1, 3
                xx = vlat(1)*qlat(1,i)+vlat(2)*qlat(2,i)+vlat(3)*qlat(3,i)
                iax(2+i,nttab) = nint(xx)
33           END DO
             nttab = nttab+1

          endif

30     END DO
20  END DO
  end subroutine ppair2

  subroutine ppair3(nttab,iax,rtab) !,iwk,iwk2,rwk)
    !- Sort neighbor table by distance
    implicit none
    integer :: nttab,niax,iwk2(nttab),i,j,k
    parameter (niax=10)
    integer :: iax(niax,nttab),iwk(niax,nttab)
    double precision :: rtab(3,nttab),rwk(3,nttab)

    do  10  i = 1, nttab
       rwk(1,i) = rtab(1,i)
       rwk(2,i) = rtab(2,i)
       rwk(3,i) = rtab(3,i)
       do  12  k = 1, niax
          iwk(k,i) = iax(k,i)
12     END DO
10  END DO
    call dvshel(3,nttab,rtab,iwk2,11)
    do  20  i = 1, nttab
       j = iwk2(i)+1
       rtab(1,i) = rwk(1,j)
       rtab(2,i) = rwk(2,j)
       rtab(3,i) = rwk(3,j)
       do  22  k = 1, niax
          iax(k,i) = iwk(k,j)
22     END DO
20  END DO
  end subroutine ppair3

  double precision function drr2(plat,tau1,tau2,i,j,k,dr)
    !     - Calculates the vector connecting two sites in a solid
    !     ----------------------------------------------------------------
    !     i Inputs
    !     i   plat: primitive lattice vectors
    !     i   tau1,tau2: basis vectors of the two sites
    !     i   i,j,k:the number of primitive lattice vectors separating sites
    !     o Outputs
    !     o   dr:   connecting vector tau2 - tau1
    !     o   drr2: square of the length of this vector
    !     r Remarks
    !     r   Using the TB package and a table of indices iax, the connecting
    !     r   vector and the square of the distance is obtained by
    !     r      rsqr = drr2(plat,bas(1,iax(1)),bas(1,iax(2)),
    !     r     .            iax(3),iax(4),iax(5),dr)
    !     ----------------------------------------------------------------
    implicit none
    integer :: i,j,k
    double precision :: dr(3)
    double precision :: plat(3,3),tau1(3),tau2(3)
    integer :: ix
    drr2 = 0.d0
    do  10  ix = 1, 3
       dr(ix) = tau2(ix) - tau1(ix) + plat(ix,1)*i + plat(ix,2)*j + plat(ix,3)*k
       drr2 = drr2 + dr(ix)**2
10  END DO
  end function drr2
  subroutine ppair4(iclus,nclus,plat,pos,ctr,iwk,rtab,tol,iax)
    !- Sort cluster by increasing (x,y,z) relative to its center
    ! ----------------------------------------------------------------
    !i Inputs
    !i   iclus,nclus: sort iax(iclus..nclus)
    !i   plat :primitive lattice vectors
    !i    pos :basis vectors
    !i    ctr :cluster origin:does not affect the ordering, but shifts rtab
    !i    iwk :integer work array of length nclus-iclus+1
    !i    tol :tolerance to which positions are considered coincident
    !i         tol<0 => sort iax by iax(1..5)
    !o Outputs
    !o   iax(7,iclus..nclus) orders the cluster by increasing (x,y,z)
    !o         (or increasing iax(1..5) if tol < 0
    !o   rtab  :connecting vectors rtab(1..3,ip) = pos(jb)-ctr
    !o          for pair ip and jb=iax(2,ip)
    !r Remarks
    !r  Each cluster is sorted by increasing (x,y,z),
    !r  sorted by x first, then by y, then by z, thus guaranteeing that
    !r  all sites common to any pair of clusters are ordered the same.
    ! ----------------------------------------------------------------
    !     implicit none
    integer :: iclus,nclus,niax,iwk(15)
    parameter (niax=10)
    integer :: iax(niax,1)
    double precision :: plat(3,3),pos(3,1),ctr(3),rtab(3,29),tol
    integer :: ic,jb,ic0,ix,ia2,i,j,k
    ! Local variables
    double precision :: dx
    !     integer jx
    !     double precision wk2(3,nclus*3)
    dx(ia2,i,j,k) = pos(ix,ia2) + &
         plat(ix,1)*i + plat(ix,2)*j + plat(ix,3)*k

    ic0 = 0
    do  12  ic = iclus, nclus
       jb = iax(2,ic)
       ic0 = ic0+1
       do  14  ix = 1,3
          rtab(ix,ic0) = dx(jb,iax(3,ic),iax(4,ic),iax(5,ic)) - ctr(ix)
14     END DO
12  END DO

    if (tol < 0) then
       call ivheap(niax,nclus-iclus+1,iax(1,iclus),iwk,1)
    else
       call dvheap(3,nclus-iclus+1,rtab,iwk,tol,1)
    endif
    do  20  ic = iclus, nclus
       iax(7,ic) = iwk(ic-iclus+1)
20  END DO
  end subroutine ppair4

  subroutine ppair5(ib1,ib2,plat,pos,tol,ntab,iax)
    !- Sort a range of clusters according to tol
    ! ----------------------------------------------------------------------
    !i Inputs
    !i  ib1,ib2:range of clusters to sort
    !i   plat  :primitive lattice vectors, in units of alat
    !i   pos   :basis vectors
    !i   tol   :tolerance; see ppair4
    !i   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
    !i   iax   :neighbor table containing pair information (pairc.f)
    !o Outputs
    !o   iax   :iax(7) is set to order cluster; see ppair4.
    !r Remarks
    !u Updates
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: niax,ib1,ib2
    parameter (niax=10)
    integer :: iax(niax,1),ntab(ib2)
    double precision :: plat(3,3),pos(3,*),tol
    integer:: ib , nttab
    integer,allocatable :: iwk1_rv(:)
    real(8) ,allocatable :: wk2_rv(:)
    ! --- Set iax(7) to sort this cluster ---
    do  10  ib = ib1, ib2
       nttab = ntab(ib+1)-ntab(ib)
       allocate(iwk1_rv(nttab*2))
       allocate(wk2_rv(nttab*3))
       call ppair4 ( ntab ( ib ) + 1 , ntab ( ib + 1 ) , plat , pos &
            , pos ( 1 , ib ) , iwk1_rv , wk2_rv , tol , iax )
       deallocate(iwk1_rv,wk2_rv)
10  END DO
  end subroutine ppair5


  subroutine siteid(iax,ncl,ntab,plat,bas,pos,iwk,nid)
    !- Assign a unique id for every different site in the neighbor table
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   iax   :array of parameters containing info about each pair
    !i   ncl   :number of clusters for which neighbor table is made
    !i   ntab  :ntab(ib) no. pairs in neighbor table preceding ib (pairs.f)
    !i   plat  :primitive lattice vectors, in units of alat (input)
    !i   bas   :basis vectors (input)
    !i   pos   :work array, of dimension 3*nttab, with nttab=ntab(ncl+1)
    !i   iwk   :work array, of dimension nttab, with nttab=ntab(ncl+1)
    !o   iax   :iax(7) is needed, which orders cluster by increasing
    !r          (x,y,z) from the origin (ppair4)
    !o Outputs
    !o   iax   :iax(10) is assigned a unique site ID for each entry
    !o   nid   :number of unique sites
    !r Remarks
    !r   A unique ID (in iax(10)) is assigned to each entry in the
    !r   neighbor table.  Sites with the same ID are identical sites.
    !r   The ordering is the same as in ppair4.
    !u Updates
    ! ----------------------------------------------------------------------
    !     implicit none
    ! ... Passed parameters
    integer :: niax,ncl,ntab(ncl+1),iwk(*),nid
    double precision :: bas(3,ncl),pos(3,*),plat(3,3)
    parameter (niax=10)
    integer :: iax(niax,1)
    ! Local variables
    integer :: ic,jb,ix,nttab,id,ip,imin,iminp,i1,icp
    double precision :: tol,dmatch(3)
    !     tolerance should be same as ppair4
    parameter (tol=1d-5)

    ! --- Table of absolute positions ---
    nttab = ntab(ncl+1)
    do  10  ic = 1, ncl
       do  12  ip = ntab(ic)+1, ntab(ic+1)
          !       ib = iax(i,ip)
          jb = iax(2,ip)
          do  14  ix = 1, 3
             pos(ix,ip) = bas(ix,jb) + plat(ix,1)*iax(3,ip) + &
                  plat(ix,2)*iax(4,ip) + plat(ix,3)*iax(5,ip)
14        enddo
12     enddo
10  enddo

    ! --- Assign a unique id ---
    !     For cluster ic, iwk(ic)=index to next site in ic needing id
    do  30  ic = 1, ncl
       iwk(ic) = 1
30  enddo
    dmatch(1) = 9d9
    dmatch(2) = 9d9
    dmatch(3) = 9d9
    id = 0
    ! --- Loop through all neighbors in the entire cluster table ---
    do  20  ip = 1, nttab
       !   ... Get first site not marked in cluster table
       imin = 0
33     imin = imin+1
       !       Skip this cluster if every member already assigned id
       if (iwk(imin) > ntab(imin+1)-ntab(imin)) goto 33
       !   ... iminp points to next element in sorted list of pairs at imin
       iminp = ntab(imin) + iax(7,ntab(imin)+iwk(imin))
       i1 = imin+1
       do  34  ic = i1, ncl
          if (iwk(ic) > ntab(ic+1)-ntab(ic)) goto 34
          icp = ntab(ic) + iax(7,ntab(ic)+iwk(ic))
          !          print 333, pos(1,iminp),pos(2,iminp),pos(3,iminp)
          !          print 333, pos(1,icp),pos(2,icp),pos(3,icp)
          !  333     format(3f12.6)
          !    ...  Exclude ic if ic(1)>imin(1)
          if (pos(1,iminp)+tol < pos(1,icp)) goto 34
          !    ...  ic becomes new imin if ic(1)<imin(1)
          if (abs(pos(1,iminp)-pos(1,icp)) > tol) goto 35
          !    ...  Exclude ic if ic(2)>imin(2)
          if (pos(2,iminp)+tol < pos(2,icp)) goto 34
          !    ...  ic becomes new imin if ic(2)<imin(2)
          if (abs(pos(2,iminp)-pos(2,icp)) > tol) goto 35
          !    ...  Exclude ic if ic(3)>imin(3)
          if (pos(3,iminp)+tol < pos(3,icp)) goto 34
          !    ...  ic becomes new imin if ic(3)<imin(3)
          if (abs(pos(3,iminp)-pos(3,icp)) <= tol) goto 34
35        continue
          imin = ic
          iminp = icp
34     enddo
       !    .. imin holds among (1..ib2) next lowest element
       iwk(imin) = iwk(imin) + 1
       !       call awrit2('iwk %n:1i',' ',180,6,ncl,iwk(1))

       !   ... If no match with previous, increment new id
       if (abs(pos(1,iminp)-dmatch(1)) > tol .OR. &
            abs(pos(2,iminp)-dmatch(2)) > tol .OR. &
            abs(pos(3,iminp)-dmatch(3)) > tol) then
          id = id+1
          dmatch(1) = pos(1,iminp)
          dmatch(2) = pos(2,iminp)
          dmatch(3) = pos(3,iminp)
       endif
       iax(10,iminp) = id
20  enddo
    nid = id

    !    .. Debugging check
    do  38  ic = 1, ncl
       if (iwk(ic)-1 /= ntab(ic+1)-ntab(ic)) call rx('bug in siteid')
38  enddo

    ! ... Debugging printout
    !      call dvheap(3,nttab,pos,iwk(1),tol,1)
    !      do  50  ip = 1, nttab
    !        jc = iwk(ip)
    !        print 333, ip, jc, iax(10,jc), pos(1,jc),pos(2,jc),pos(3,jc)
    !  333   format(3i6,3f12.6)
    !   50 continue
    !      call rx('done')
  end subroutine siteid
  subroutine nsitsh(mode,ia,ib,iax,ntab,nlst,lsta,lstb)
    !- Count the number of sites two clusters share in common
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :0 count the sites only
    !i          1 also make a list of the sites
    !i   ia,ib :pair of sites for which to seek common elements
    !i   nds   :leading dimension of sid
    !i   iax   :neighbor table containing pair information (pairc.f)
    !i   ntab  :ntab(ib)=# pairs in iax table preceding ib (pairc.f)
    !o Outputs
    !o   nlst  :number of sites in common
    !o   lsta  :list of
    !r Remarks
    !u Updates
    ! ----------------------------------------------------------------------
    !     implicit none
    !     Passed parameters
    integer :: mode,ia,ib,ntab(1),niax,lsta(1),lstb(1),nlst
    parameter (niax=10)
    integer :: iax(niax,1)
    !     Local variables
    integer :: ic,ica,icb,n,na,nb,ipa,ipb,sida,sidb

    ica = ntab(ia)+1
    icb = ntab(ib)+1
    ipa = iax(7,ica)
    ipb = iax(7,icb)
    sida = iax(10,ipa+ntab(ia))
    sidb = iax(10,ipb+ntab(ib))
    na = ntab(ia+1)
    nb = ntab(ib+1)
    n = na-ica + nb-icb + 2
    nlst = 0
    do  10  ic = 1, n
       if (ica > na .OR. icb > nb) return

       !   ... A match; increment nlst
       if (sida == sidb) then
          print *, sida
          nlst = nlst+1
          ica = ica+1
          icb = icb+1
          ipa = iax(7,ica)
          ipb = iax(7,icb)
          sida = iax(10,ipa+ntab(ia))
          sidb = iax(10,ipb+ntab(ib))
          if (mode /= 0) then
             lsta(nlst) = ica
             lstb(nlst) = icb
          endif
       elseif (sida > sidb) then
          icb = icb+1
          ipb = iax(7,icb)
          sidb = iax(10,ipb+ntab(ib))
       else
          ica = ica+1
          ipa = iax(7,ica)
          sida = iax(10,ipa+ntab(ia))
       endif
10  enddo
  end subroutine nsitsh

end module m_lmaux

