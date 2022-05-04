      module m_lmaux
      public:: lmaux
      private
      contains
      subroutine lmaux()     !main part of lmchk
      use m_mksym,only: ctrl_nclass,iv_a_oics,iv_a_oipc
      use m_lmfinit,only: iv_a_oips,str_mxnbr,str_rmax,ctrl_nbas,ctrl_nspec,ctrl_nspin,
     &     ctrl_nl,ctrl_modep,ctrl_omax1,ctrl_omax2,ctrl_wsrmax,slabl,sspec=>v_sspec,
     &     lat_avw,lat_alat,mxcst2
      use m_lattic,only: lat_nkd
      use m_lattic,only: lat_nkq
      use m_struc_def  
      use m_ovmin , only: ovmin
      use m_lattic,only:lat_plat,rv_a_opos
!! check crystal structure symmetry and get WSR
C ----------------------------------------------------------------------
Ci Inputs
Cu Updates
Cu   12 Aug 08 (L. Ke) empty sphere finder
Cuxxx   04 Nov 04 Upgrade of rsta editor
Cuxxx   26 Jan 03 Call to angtab changed
Cu   17 May 02 Modified MT radii scaling to lower priority for E.S.
Cu   23 Apr 02 Added option (--getwsr) to find MT radii
Cu   01 Mar 02 Updated Import data mode
Cuxxx   05 Oct 01 Adapted mode 2**3 to work with lm v6.11
Cu   24 Nov 97 changed ovmin to run quickly
C ----------------------------------------------------------------------
      implicit none
      integer:: mode=1 !,wksize
      character*120 outs,fnam*8
      integer NULLI
      logical cmdopt,T,F,swtmp
      parameter (T=.true., F=.false., NULLI=-99999)
      integer getdig,i,ip,j,k,m,ifi,iprint,lpbc,lgunit,
     .nbas,nclasp,nclass,nl,nlspc,nsp,modep(3),parg,nbasp,
     .nbaspp,nkd,nkq,nspec,neul,nc,mxcsiz,nttab,igets, !npadl,npadr,
     .iosits,cmplat,ngrp,ival,stdo,irs(5),nclspp,bitand,igetss,
     .ngmx,nsgrp
      integer:: oeold  , olmx , opold , owk2 , orham , oamsh 
     .     , onrmsh , oalpha , onpr , os , ormx , oip , opgfsl , mxclas
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
      double precision xv(10),xx,alat,plat(3,3),facrmx,facrng, !,plat2(9)
     .dval,avw,ekap(2),enu,qss(4),ckbas,cksumf,ehterm(4), rmaxs,
     .qlat(9),emad,trumad,vmtz(2),omax1(3),omax2(3),wsrmax
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
      modep = ctrl_modep
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
      print *,' nclasp,lmxa=',nclasp,lmxa
      print *,' z   =',z
      nbasp = nbas !+ npadl + npadr
      nbaspp = nbas !2*nbasp - nbas
      stdo = lgunit(1)
      j = 10
      if (cmdopt('--shorten',j-1,0,outs)) then
        call shorps ( nbasp , plat , modep , rv_a_opos , rv_a_opos )
      endif
C --- Neighbor tables and sphere overlaps ---
      if (getdig(mode,0,2) .ne. 0) then
        if (rmaxs .le. 0d0) then
          rmaxs = 2.7d0*avw
          call info5(30,0,0,'%1f'//
     .    'Use default rmaxs = %;3d a.u. = %;3d*avw = %;3d*alat',
     .    rmaxs,rmaxs/avw,rmaxs/alat,0,0)
        endif
C ... Get neighbor table iax for each atom in the cluster
        if (lpbc .eq. 0) then
          i = 3
          j = -1
        elseif (lpbc .eq. 1 .or. lpbc .eq. 11) then
          i = 2
          j = 1
        else
          call rx('ASASTR: not implemented for lpbc>1')
        endif
        mxcsiz = str_mxnbr !int(sstr%mxnbr)

        call pshpr(iprint()-20)
        call pairs ( nbas , nbasp , alat , plat ,(/ rmaxs / 2/) , rv_a_opos
     .   , (/- 1/) , i , j , w_dummy , nttab , iv_a_ontab , iv_a_oiax , mxcsiz )
        call poppr

C --- Print out a few superlattice vectors ---
        j = 6
        if (cmdopt('--slat',j-1,0,outs)) then
          if (iprint() .ge. 10) then
            call info0(10,1,0,' LMCHK:  print multiples of plat%N'//
     .      '  i1  i2  i3%7fx%11fy%11fz%11flen')
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
  368             format(3i4, 3f12.7, 1x, f12.5)
                enddo
              enddo
            enddo
          endif
        endif

C --- Find sphere overlaps ---
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
          if (lpbc .eq. 0) then
            i = 3
          elseif (lpbc .eq. 1 .or. lpbc .eq. 11) then
            i = 2
          else
            call rx('LMAUX: not implemented for lpbc>1')
          endif
          call makrm0 ( 101 , nspec , nbas , alat , plat , rv_a_opos ,
     .      slabl , iv_a_oips , modep , lock_iv , zz_rv , rmt_rv ) !sarray%
C   ... Scale sphere radii satisfying constraints
          i_copy_size=size(ctrl_omax1)
          call dcopy(i_copy_size,ctrl_omax1,1,omax1,1)
          i_copy_size=size(ctrl_omax2)
          call dcopy(i_copy_size,ctrl_omax2,1,omax2,1)
          wsrmax=ctrl_wsrmax
          call sclwsr ( 20 , nbas , nbasp , nspec , alat , plat , rv_a_opos
     .      , iv_a_oips , modep , slabl , zz_rv , lock_iv , 1d0 , wsrmax 
     .      , omax1 , omax2 , rmt_rv )
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
c-------
        if(ifx==0) then
          allocate(rmax(nclasp))
          call dcopy ( nclasp , rv_a_ormax , 1 , rmax , 1 )
        endif
C ... Write positions in Cartesian coordinates and as multiples plat
        if (iprint() .ge. 50) then
          write(stdo,357)
  357     format(/' site spec',8x,'pos (Cartesian coordinates)',9x,
     .    'pos (multiples of plat)')
C     qlat = (plat+)^-1
          call dinv33(plat,1,qlat,xx)
          do  i = 1, nbas
            call dpscop ( rv_a_opos , xv , 3 , 3 * i - 2 , 1 , 1d0 )
C       posp+ = (plat)^-1 pos+
            call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
            ip = ival ( iv_a_oips , i )

            print 345, i, slabl(ip), (xv(j),j=1,3), (xv(3+j),j=1,3)
  345       format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
          enddo
        endif
C --- Print overlaps, optionally minimize wrt spec'd sites ---
        outs = ' '
        i = 6
        swtmp = cmdopt('-mino',5,0,outs)
        swtmp = cmdopt('--mino',6,0,outs)
        if (swtmp) i = 7
        j = 1
        if (iprint() .lt. 30) j = 0
        call ovmin ( outs ( i: ) , nbas , nbasp , alat , plat , rmax
     .     , rmax , slabl , iv_a_oipc , modep , z , iv_a_ontab , iv_a_oiax ,
     .     rv_a_opos , j )
      endif
C --- Interpolate core to another mesh ---
      if (getdig(mode,4,2) .ne. 0) then
        call rx('patch clabl for call to coritp')
C       call coritp(nclass,nsp,w(oclabl),nrmsh,amsh,w(ormax))
      endif
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
      subroutine makrm0(opt,nspec,nbas,alat,plat,pos,slabl,ips,modep,
     .lock,z,rmt)
      use m_lmfinit,only: nsp
C- Estimate muffin-tin radii overlapping atomic potentials
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt  :specifies options in makrm0.
Ci         :
Ci         :1s digit : prescription for constructing potential
Ci         :  The only option is mode 1 (see Remarks):
Ci         :  set 1's digit opt=1.
Ci         :
Ci         :10s digit : quantity to overlap
Ci         :  0 overlap electrostatic atomic potentials
Ci         :  1 overlap estat atomic potentials + add total xc potential
Ci         :    This option tends to produce slightly more uniform radii.
Ci         :  2 overlap total atomic potentials
Ci         :    This option probably doesn't make sense.
Ci         :
Ci         :100s digit specifies what MT radii are returned.
Ci         :Initially MT radii are generated by site.
Ci         :  0 return MT radii by site:
Ci              rmt(i) corresponds to radius for site i
Ci         :  1 average MT radii by species:
Ci              rmt(i) corresponds to radius for species i
Ci         :  2 smallest MT radii of this species
Ci              rmt(i) corresponds to radius for species i
Ci         :
Ci         :1000s digit specifies whether nspec,slabl,ips refer to
Ci                species or classes (used for printout only)
Ci         :  0 quantities refer to species
Ci         :  1 quantities refer to classes
Ci   nspec :number of species (or classes)
Ci   nbas  :size of basis
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors, in units of alat
Ci   slabl :vector of species (or class) labels
Ci   ips   :species (or class) table: site ib belongs to species ips(ib)
Ci   modep :integer vector of length 3 governing how lattice vectors
Ci         :are shortened (shorps).  In particular,
Ci         :0 suppresses shifts along plat(j)
Ci         :2 shifts to minimize length of pos
Ci   z     :table of nuclear charges by species (or class)
Co Outputs
Co   rmt
Cl Local variables
Cl         :
Cr Remarks
Cr   makrm0 attempts to determine an optimal set of muffin-tin radii so
Cr   that the errors due to the muffin-tin approxiation are minimized.
Cr   Inside the MT-spheres the potential should be as spherical as
Cr   posssible, and the potential should be as constant as possible.
Cr
Cr   At present there is only one prescription (mode 1).
Cr   This routine follows the methodology of the Stuttgart code potmax.
Cr
Cr   mode 1 : rmt are chosen as follows:
Cr   1. The free-atom potential is constructed for all species.
Cr   For each atom:
Cr   2. A neighbor table is made for the atom.
Cr   3. The connecting vectors to each neighbor make up a group
Cr         of direction vectors stellating from that atom.
Cr   4. The first occurence of a maximum in the potential along any
Cr         of these vectors determines rmt for the atom
Cr
Cr   Finally, rmt may be averaged over species (or classes).
Cu Updates
Cu   21 Apr 02 First created.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      character*8 slabl(*)
      integer opt,nspec,nbas,ips(nbas),modep(3),lock(*)
      double precision z(nspec),alat,plat(9),pos(3,nbas),rmt(*)
C Local variables:
      character*8 spid,strn*80
      integer nrmx,nxi0,n0,niax,mxcsiz,npmx,pass
      parameter (nrmx=1501,nxi0=10,n0=10,niax=10,mxcsiz=400,npmx=mxcsiz)
      integer nxi,nrmix(2),lxcfun,is,opt0,opt1,opt2,opt3
      integer idmod(n0),lmxa,nrmt,nr,ib,nttab,k,ir,np,ipr,stdo,lgunit,
     .nrspec
      integer ntab(nbas+1),iax(niax,mxcsiz),ipa(nbas)
      double precision qc,ccof,ceh,rmtl(nspec),rfoca,rsmfa,qcor(2),
     .sumec,sumtc,eref,fac,ddot,fpi,etot
      double precision a(nspec),b(nspec),xx(5)
      double precision pnu(n0,2),pnz(n0,2),qat(n0,2)
      double precision hfc(nxi0,2),exi(nxi0),hfct(nxi0,2)
      double precision rtab(n0,2),etab(n0,2)
      double precision v(nrmx,nspec+1),rho(nrmx,nspec+1),
     .rhoc(nrmx,nspec+1),rofi(nrmx*2),range(nbas),rmti(nbas)
      double precision vp(npmx,0:2),xp(3,npmx),rp(npmx),
     .vxcp(npmx),excp(npmx)
      real(8) ,allocatable :: excx_rv(:)
      real(8) ,allocatable :: excc_rv(:)
      real(8) ,allocatable :: vxcx_rv(:)
      real(8) ,allocatable :: vxcc_rv(:)

C     Sets scale for neighbor table
      double precision facr,facri
      parameter (facr=2d0,facri=1.02d0)
C     Relative positions
      double precision rpos(3,mxcsiz),ri,rbar,rmin,vrmax,cur,slo,di
      character*1 sym(2)
C ... External calls
      external dcopy,defpq,defwsr,dpcopy,dpzero,dscal,freats,getpr,
     .iinit,info,polint,poppr,pshpr,psymq0,psymr0,rx,rxi,sumsro
C ... Heap
      data sym /' ','*'/

      logical:: isanrg, l_dummy_isanrg
      integer:: ifives

      print *,'makrm0:'
C --- Setup ---
      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
      opt2 = mod(opt/100,10)
      opt3 = mod(opt/1000,10)
Ckino isanrg is logical function,       call isanrg(opt0,1,1,'makrm0:','1s digit opt',.true.)
      l_dummy_isanrg=isanrg(opt0,1,1,'makrm0:','1s digit opt',.true.)
Ckino isanrg is logical function,       call isanrg(opt1,0,2,'makrm0:','10s digit opt',.true.)
      l_dummy_isanrg=isanrg(opt1,0,2,'makrm0:','10s digit opt',.true.)
Ckino isanrg is logical function,       call isanrg(opt2,0,1,'makrm0:','100s digit opt',.true.)
      l_dummy_isanrg=isanrg(opt2,0,1,'makrm0:','100s digit opt',.true.)
Ckino isanrg is logical function,       call isanrg(opt3,0,1,'makrm0:','1000s digit opt',.true.)
      l_dummy_isanrg=isanrg(opt3,0,1,'makrm0:','1000s digit opt',.true.)
Changenglob      nsp = nglob('nsp')
c      nsp = globalvariables%nsp
      stdo = lgunit(1)
      call getpr(ipr)
      call dpzero(rpos,3*mxcsiz)
      call dpzero(v,nrmx*nspec)
      call dpzero(rho,nrmx*nspec)
      call dpzero(rhoc,nrmx*nspec)
      call dpzero(vp,npmx*3)
      fpi = 16d0*datan(1d0)

C --- Densities and potentials for each species ---
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

C       nrspec(is) = iabs(iclbsj(is,ips,-nbas,nbas))
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
        call freats(spid,is,nxi0,nxi,exi,rfoca,rsmfa,0,-1,qcor,nrmix,0,
     .  lxcfun,z(is),rmtl(is),a(is),nrmt,pnu,pnz,qat,0d0,0d0,0d0,0d0,
     .  idmod,lmxa,eref,rtab,etab,hfc,hfct,nr,rofi,rho(1,is),rhoc(1,
     .  is),qc,ccof,ceh,sumec,sumtc,v(1,is),etot, 1, -1,-1) !nmcore=1 july2012 !last -1 -1 means not write ves* files.
        call poppr
        b(is) = rmtl(is)/(dexp(a(is)*(nrmt-1)) - 1d0)
C       Restore 4*pi*r**2*(total density) in array rho
        call daxpy(nrmx*nsp,1d0,rhoc(1,is),1,rho(1,is),1)

C       Use estat potential: overwrite v with estat
        if (opt1 .le. 1) then
C         call prrmsh('vtot',rofi,v(1,is),nr,nr,1)
          vrmax = 0
          call poiss0(z(is),a(is),b(is),rofi,rho(1,is),nr,vrmax,v(1,is),
     .    xx(2),xx(4),nsp)
C         call prrmsh('ves',rofi,v(1,is),nr,nr,1)
        endif

        if (nsp .eq. 2) then
          call daxpy(nr,1d0,rho(1+nr,is),1,rho(1,is),1)
          call daxpy(nr,1d0,rhoc(1+nr,is),1,rhoc(1,is),1)
          call daxpy(nr,1d0,v(1+nr,is),1,v(1,is),1)
          call dscal(nr,.5d0,v(1,is),1)
          call rx('need check this branch')
        endif

C       Overwrite v,rho with proper spherical potential and density
        do  ir = 2, nr
          v(ir,is) = v(ir,is) - 2*z(is)/rofi(ir)
          rho(ir,is) = rho(ir,is)/(fpi*rofi(ir)**2)
        enddo
C       call prrmsh('rho',rofi,rho(1,is),nr,nr,1)
C       call prrmsh('v',rofi,v(1,is),nr,nr,1)
      enddo
C

C ... Sets range for neighbor table
      do  ib = 1, nbas
        is = ips(ib)
        range(ib) = facr * rmtl(is)
      enddo

C --- For each site, find rmti = initial estimate for rmtl ---
      do  ib = 1, nbas

        is = ips(ib)
        spid = slabl(is)

C   ... Neighbor table connecting sites to this one
        nttab = mxcsiz
        call pairc(ib,ib,nbas,modep,20,[0],alat,plat,pos,pos,range,-1,[1],
     .  nttab,ntab,iax,rpos,k)
        if (nttab.gt.mxcsiz)
     .  call rxi('makrm0 : increase mxcsiz: need',nttab)

C   ... Maximum in potential ( mode 1)
C       Find first maximum along all direction vectors
C       which are proportional to rpos, excluding self (first site).
C       As a conservative choice, choose initial radius = rmtl/4.
        np = nttab-1
        if (nttab .gt. npmx)
     .  call rxi('makrm0 : increase npmx: need',nttab)
        call dpzero(xp,3*nttab)
        ri = rmtl(is)/4
        if (z(is) .eq. 0) then
          rmti(ib) = rmtl(is)
          goto 31
        endif

C       Re-entry for loop until maximum in potential found
        pass = 0
   30   continue
C       Make nttab-1 connecting vectors of radius ri
        do  k = 2, nttab
          fac = ri/dsqrt(ddot(3,rpos(1,k),1,rpos(1,k),1))
          call dpcopy(rpos(1,k),xp(1,k-1),1,3,fac)
        enddo

C       Add superposition of potentials for each vector
        call sumsro(xp,np,ips,a,b,v,nttab,iax,rpos,vp(1,2))

C       Add vxc[superposition of densities] for each vector
        if (opt1 .eq. 1) then
          call sumsro(xp,np,ips,a,b,rho,nttab,iax,rpos,rp)
          allocate(excx_rv(np))
          allocate(excc_rv(np))
          allocate(vxcx_rv(np))
          allocate(vxcc_rv(np))
          call evxcv ( rp , rp , np , 1 , lxcfun , excp , excx_rv , 
     .    excc_rv , vxcp , vxcx_rv , vxcc_rv )
          call daxpy(np,1d0,vxcp,1,vp(1,2),1)
          deallocate(vxcc_rv)
          deallocate(vxcx_rv)
          deallocate(excc_rv)
          deallocate(excx_rv)
        endif

C       Check against prior point.  Skip if nothing yet to compare
        if (pass .gt. 0) then
          do  k = 1, nttab-1
C           This loop executes if a maximum is found
            if (vp(k,1) .gt. vp(k,2)) then
C             It shouldn't happen already at the first point
              if (pass .eq. 1) then
                call
     .          rx('makrm0: starting r is too small.  Check geometry')
              endif
              rmti(ib) = ri/dsqrt(facri)
              if (pass .gt. 1) then
                cur = vp(k,2) + vp(k,0) - 2*vp(k,1)
                slo = (vp(k,2) - vp(k,0))/2
                di = -slo/cur
                if (dabs(di) .lt. 1) then
C                 print *, pass, di, ri, (ri/facri)*dexp(di*(facri-1))
                  rmti(ib) = (ri/facri)*dexp(di*(facri-1))
                endif
              endif
              goto 31
            endif
          enddo
        endif

C       Copy vp(*,2) to vp(*,1), and increment pass
        call dcopy(npmx,vp(1,1),1,vp(1,0),1)
        call dcopy(npmx,vp(1,2),1,vp(1,1),1)
        pass = pass+1
        ri = ri * facri
        if (ri .gt. 10d0)
     .  call rx('makrm0: failed to find maximum in potential')
        goto 30

C       Loop exit: radius rmti(ib) has been found
   31   continue

C       This is the connecting vectory along which the maximum was found
C       iconn = iax(2,k+1)

      enddo

C --- Copy appropriate rmti to final result, depending on opt2 ---
      strn = 'estat'
      if (opt1 .eq. 1) strn = 'LDA'
      if (opt1 .eq. 2) strn = 'sum-of-atomic'
      call info(30,1,0,' makrm0: initial MT radii from first '//
     .strn(1:20)//'%a potential maximum',0,0)
      call info(30,0,0,'  site   spec%12frmt'//
     .'%7frmt-%7frmt-%7frold%3flock%N%34f<spec avg>  spec-min',0,0)

C ... First loop for printout only, so can hang on to original rmt
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
C       Index where to poke rnew
        k = ib
        if (opt2 .ge. 1) k = is

        if (ipr .ge. 30) then
          write(stdo,344) ib,is,spid,rmti(ib),rmti(ib)-rbar,
     .    rmti(ib)-rmin,rmt(k), sym(1+mod(lock(k)/2,2))
  344     format(i5,2x,i3,':',a,4f11.4,3x,a)
        endif

      enddo

C     print *, rmt

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
C       Index where to poke rnew
        k = ib
        if (opt2 .ge. 1) k = is

        if (lock(k) .ne. 2) then
          if (opt2 .eq. 0) rmt(ib) = rmti(ib)
          if (opt2 .eq. 1) rmt(is) = rbar
          if (opt2 .eq. 2) rmt(is) = rmin
        endif
      enddo
C    print *, rmt
      end
      
      subroutine sclwsr(opt,nbas,nbasp,nspec,alat,plat,basp,ips,modep,
     .slabl,z,lock,volfac,wsmax,omax1,omax2,wsr)
      use m_ext,only: sname
C- Scales sphere radii to specified volume, satisfying constraints.
C ----------------------------------------------------------------------
Ci Inputs:
Ci   opt  :specifies options in sclwsr.
Ci         :1s digit specifies whether scaling reaches target volfac
Ci         :         is required:
Ci         : 0 meeting volume target is optional
Ci         : 1 meeting volume target is optional but
Ci         :   a warning is printed if not met.
Ci         : 2 meeting volume target is required
Ci         :10s digit concerns treatment of empty spheres
Ci         : 0 ES and other sites are treated symmetrically
Ci         : 1 all sites with z>0 are resized first; then
Ci             all sites are resized.
Ci         : 2 all sites with z>0 are resized first; then
Ci             the ES sites only are resized.
Ci   nbas  :size of basis
Ci   nbasp :size of padded basis (layer programs)
Ci   nspec :number of atom species
Ci   alat  :length scale
Ci   plat  :primitive lattice vectors, in units of alat
Ci   basp  :basis vectors (scaled by alat; padded for layer mode)
Ci   ips   :the jth atom belongs to species ips(j)
Ci   modep :specifies which dimensions have periodic b.c.
Ci         :In particular,
Ci         :0 suppresses shifts along plat(j)
Ci         :2 shifts to minimize length of pos
Ci          are used
Ci   slabl :species labels
Ci   z     :nuclear charge, by species
Ci   lock  :constraints specifying which species are locked and
Ci         :which are free to float.
Ci         :On input, lock should be zero or two for each species.
Ci         :Each species for which lock(i)=2 is constrained not to change.
Ci         :The radii for other species are floated.
Ci         :lock(1..nspec) is OVERWRITTEN on output
Ci   volfac:scale until sum of sphere vol = volfac * cell vol or until
Ci         :all sphere radii are constrained (see Remarks, sclws2)
Ci   wsmax :(wsmax>0) a global maximum on the size of MT spheres.
Ci         :No sphere is allowed to exceed wsmax.
Ci         :(wsmax=0) no constraint is imposed.
Ci   omax1 :max. allowed overlap divided by distance (r1+r2-d)/d<omax1
Ci         :omax1(1) is constraint for atom-atom overlaps
Ci         :omax1(2) is constraint for atom-empty-sphere overlaps
Ci         :omax1(3) is constraint for empty-sphere-empty-sphere overlaps
Ci   omax2 :max. allowed overlap divided by radius  (r1+r2-d)/r1<omax2
Ci         :omax2(1) is constraint for atom-atom overlaps
Ci         :omax2(2) is constraint for atom-empty-sphere overlaps
Ci         :omax2(3) is constraint for empty-sphere-empty-sphere overlaps
Co Inputs/Outputs:
Cio  wsr   :Wigner-Seitz sphere radius (in atomic units)
Cio        :On input, starting values for wsr
Cio        :On output, final values for wsr
Cl Local:
Cl   dovl1 :maximum overlap divided by distance
Cl   dovl2 :maximum overlap divided by radius
Cl   gamma :(zero) passed to sclws2
Cl   nrspec:number of atoms in each species
Cr Remarks
Cr  Sphere radii are scaled in an iterative procedure.  In any
Cr  iteration, species are divided into those that are 'locked'
Cr  (frozen) and those that are allowed to float.  The largest scaling
Cr  factor is determined for all those species of the latter type, that
Cr  satisfies the constraints (see below).  These species are scaled
Cr  and a new iteration begins.  By construction, each iteration will
Cr  cause at least one new species to be locked; thus, the total number
Cr  of iterations will not exceed the number of species.
Cu Updates
Cu   17 Jan 09  bug fix 10s digit mode=2, no ES
Cu   17 May 02  New 10s digit opt switch gives optional lower priority
Cu              to empty spheres
Cu   22 Apr 02  First created; adapted from Stuttgart LMTO56.
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer opt,nbas,nbasp,nspec,ips(nbas),modep(3),lock(nspec)
      double precision alat,basp(3,nbasp),volfac,wsmax,
     .omax1(3),omax2(3),plat(3,3),wsr(nspec),z(nspec)
      character*8 slabl(*)
      integer:: i , k , nrspec(nspec) , ib , is , opt1
      real(8) ,allocatable :: wk_rv(:)
      integer niax,mxcsiz,mxnbr,nttab,ipr,stdo,lgunit,llock(nspec)
      integer ,allocatable :: ntab_iv(:)
      integer ,allocatable :: iax_iv(:)
      logical les
      double precision dovl1(3),dovl2(3),gamma,range(nbas),wsrs(nspec)
      double precision facr,tiny,avw,vol,avwsr,volnew,volnes,volold,volsph
      parameter(niax=10,mxcsiz=200,facr=2d0,tiny=1d-5)
      integer:: istdo

C     omax1(1) = -.01
C      wsr(1) = 3.232247d0
C      wsr(2) = 3.232247d0
C      wsr(3) = 2.248243d0
C      wsr(4) = 2.097284d0
C      wsr(5) = 1.647651d0
C      wsr(6) = 1.647651d0

C --- Setup ---
      call getpr(ipr)
      stdo  = lgunit(1)
      gamma = 0
c      call maknrs(nbas,ips,ib,nrspec)
      nrspec=0
      do  ib = 1, nbas
        is = ips(ib)
        nrspec(is) = nrspec(is) + 1
      enddo
c      if (ib .gt. nspec) call rx('sclwsr: wrong number of species')
      opt1 = mod(opt/10,10)
      call dcopy(nspec,wsr,1,wsrs,1)
      avw = avwsr(plat,alat,vol,nbas)

C --- Make a neighbor table and adjust llock to freeze ES sites ---
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
      call pairc ( 1 , nbas , nbasp , modep , 20 , [0] , alat , plat 
     ., basp , basp , range , - 1 , [1] , nttab , ntab_iv , iax_iv 
     ., wk_rv , k )

      call poppr

C --- Scale the sphere radii, freezing empty spheres ---
      if (opt1 .ge. 1) then

C   ... Local copy of lock and wsr, adjusting to freeze ES sites
        call icopy(nspec,lock,1,llock,1)
        les = .false.
        do  is  = 1, nspec
          if (z(is) .eq. 0) then
            les = .true.
            llock(is) = 2
            wsr(is) = 0
          endif
        enddo

        if (les) then

C     ... Scale wsr with z=0 sites locked at wsr=0
          volold = volsph(nspec,nrspec,wsr)/vol
          volnes = volfac
          call sclws2 ( nbas , nspec , alat , plat , basp , slabl , iax_iv 
     .    , ips , ntab_iv , z , nrspec , omax1 , omax2 , gamma , wsmax 
     .    , llock , volnes , dovl1 , dovl2 , wsr )

          if (ipr .ge. 30) then
c            call awrit2(' SCLWSR:  initial sphere packing = %;1d%%'//
c     .      ' scaled to %;1d%% (no empty spheres)',
c     .      ' ',120,stdo,100*volold,100*volnes)
             write(stdo,"(' SCLWSR:  Initial sphere packing = ',f10.5,
     .       ' scaled to (no empty spheres)',f10.5)")  100*volold,100*volnes
          endif

C     ... Restore wsr(Z=0)
          call icopy(nspec,lock,1,llock,1)
          do  is  = 1, nspec
            if (z(is) .eq. 0) wsr(is) = wsrs(is)
          enddo
        endif

C   ... 10s digit opt=2 : freeze wsr(Z>0)
        if (opt1 .ge. 2 .and. les) then
          do  is  = 1, nspec
            if (z(is) .ne. 0) lock(is) = 2
          enddo
        endif

      endif

C --- Scale the sphere radii (2nd pass for les) ---
      volold = volsph(nspec,nrspec,wsrs)/vol
      volnew = volfac
      call sclws2 ( nbas , nspec , alat , plat , basp , slabl , iax_iv 
     ., ips , ntab_iv , z , nrspec , omax1 , omax2 , gamma , wsmax 
     ., lock , volnew , dovl1 , dovl2 , wsr )


C --- Printout ---
Cki#error, have return with len(w_varlist)>0 at line 183
      if ( ipr .lt. 10 ) then
        if (allocated(wk_rv)) deallocate(wk_rv)
        if (allocated(iax_iv)) deallocate(iax_iv)
        if (allocated(ntab_iv)) deallocate(ntab_iv)
        return
      endif

C      if (ipr .ge. 10) write(stdo,309)
C     .  vol,100*volold,100*volnew
C  309 format(/' SCLWSR: vol=',f11.3,
C     .  '  sphere fraction=',f5.1,
C     .  '%(initial)  ',f5.1,'%(scaled)')

      if (ipr .ge. 10) then
        write(stdo,"(' SCLWSR:  vol = ',f10.5,' a.u.**3 ',
     .  ' Initial sphere packing = ',f10.5,
     .  ' scaled to ',f10.5)")   vol,100*volold,100*volnew
      endif

      if (ipr .ge. 30) then
        write(stdo,310)
     .  (omax1(i)*100,i=1,3),(omax2(i)*100,i=1,3),
     .  (dovl1(i)*100,i=1,3),(dovl2(i)*100,i=1,3)
  310   format(1x,'constr omax1=',3f6.1,' %    omax2=',3f6.1,' %',
     .  /1x, 'actual omax1=',3f6.1,' %    omax2=',3f6.1,' %')
c        if (ipr .gt. 30) then
        write(stdo,311)
  311   format(/' spec  name',8x,'old rmax    new rmax     ratio')
        do  i = 1, nspec
          write(stdo,312) i,slabl(i), wsrs(i), wsr(i), wsr(i)/wsrs(i)
  312     format(i4,3x,a,3f12.6)
        enddo
c        endif
      endif

      open(newunit=istdo,file='rmt.'//trim(sname))
      do  i = 1, nspec
        write(istdo,"(a,f12.6)") slabl(i), wsr(i)
      enddo
      close(istdo)

      if (dabs(volnew-volfac) .gt. tiny) then
        if (mod(opt,10) .ge. 1) then
          write(stdo,321) int(volfac*100)
  321     format(/' SCLWSR (warning): failed to reach target vol (',
     .    i3,'% of cell vol)')
        endif
        if (mod(opt,10) .eq. 2) then
          call rx('SCLWSR: failed to reach target VOL.  Increase omax.')
        endif
      endif
      if (allocated(ntab_iv)) deallocate(ntab_iv)
      if (allocated(iax_iv)) deallocate(iax_iv)
      if (allocated(wk_rv)) deallocate(wk_rv)
      end subroutine sclwsr


      subroutine sclws2(nbas,nspec,alat,plat,bas,slabl,iax,ips,ntab,z,
     .nrspec,omax1,omax2,gamma,wsmax,lock,volfac,dovl1,dovl2,wsr)
C- Enlarges the spheres to reach a specified volume under constraints
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas  :size of basis
Ci   nspec :number of species
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   bas   :basis vectors, in units of alat
Ci   slabl :species labels (for printout)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   z     :nuclear charge, by species
Ci   nrspec:number of atoms in the ith species
Ci   omax1 :max. allowed overlap divided by distance (s1+s2-d)/d<omax1
Ci   omax2 :max. allowed overlap divided by radius  (s1+s2-d)/s1<omax2
Ci   gamma :a factor that changes scaling wsr from simple multiplicative
Ci         :scaling to a combination of additive + multiplicative scaling
Ci         :That is, in each iteration,
Ci         :scaling is r -> a(r+b) with a*b=gamma*(a-1)*avw
Ci         :gamma>0 tends enlarge small spheres faster than large ones
Ci   wsmax :(wsmax>0) a global maximum on the size of MT spheres.
Ci         :No sphere is allowed to exceed wsmax.
Ci         :(wsmax=0) no constraint is imposed.
Cio Inputs/Outputs:
Cio  lock  :constraints specifying which species are locked and
Cio        :which are free to float.
Cio        :On input, lock should be zero or two for each species.
Cio        :Each species for which lock(i)=2 is constrained not to change.
Cio        :The radii for other species are floated.
Cio        :lock(1..nspec) is OVERWRITTEN on output
Cio  volfac:scale until sum of sphere vol = volfac * cell vol
Cio        :or until all sphere radii are constrained (see Remarks)
Cio        :On input, volfac=target ratio (sum of sphere vol)/(cell vol)
Cio        :On output, volfac=actual ratio
Cio  wsr   :Wigner-Seitz sphere radius (in atomic units)
Cio        :On input, starting values.
Cio        :On output, scaled values.
Co Outputs:
Co   dovl1 :maximum overlap divided by distance
Co   dovl2 :maximum overlap divided by radius
Cr Remarks:
Cr  Sphere radii are scaled in an iterative procedure.  In any
Cr  iteration, species are divided into those that are 'locked'
Cr  (frozen) and those that are allowed to float.  The largest scaling
Cr  factor is determined for all those species of the latter type, that
Cr  satisfies the constraints (see below).  These species are scaled
Cr  and a new iteration begins.  By construction, each iteration will
Cr  cause at least one new species to be locked; thus, the total number
Cr  of iterations will not exceed the number of species.
Cr
Cr  Typically sclws2 is called by is a higher-level routine that creates
Cr  the neighbor-table and other necessary arrays.
Cr
Cr  This code was adapted from Stuttgart routine blowup, v LMTO56.
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer niax
      parameter (niax=10)
      integer nbas,nspec,iax(niax,*),ntab(*),ips(nbas),nrspec(nspec)
      integer lock(nspec)
      double precision alat,bas(3,nbas),dovl1(3),dovl2(3),wsmax,
     .volfac,gamma,omax1(3),omax2(3),plat(3,3),wsr(nspec),z(nspec)
      character*8 slabl(nspec)
C Local variables:
      integer ib,jb,is,iclbsj,ip,k1,k2,k3,
     .kb,ks,kco,kpr,lgunit,nloop,stdo,ipr,npcol,locki
      double precision a,amax,amax1,amax2,amax3,amax4,avw,b,bmax,d,
     .dm(0:3),dovlap,dscl(nspec),
     .dsclmx,drr2,dr(3),fpi3,gw,opo1,omo2,p,q,r,ratio,
     .rik,riko,s,t,tiny,u,v,vol,vola,volb,
     .volsph,wsri,wsrk,x,avwsr,Vconst
      logical fin
      character*72 fmt
      parameter(fpi3=4.18879020478639053d0,tiny=1d-5)
C ... External calls
      external dcopy,dpzero,dscal,errmsg,getpr,iinit

      avw  = avwsr(plat,alat,vol,nbas)
      gw   = avw*gamma
      stdo = lgunit(1)
      call getpr(ipr)
C      print *, '!!'
C      ipr = 40
      npcol = 7

C     call dcopy(nspec,1d0,0,dscl(1,0),1)

      fmt = '(''   SPEC:    '',6(1x,a8):/(12x,6(1x,a8)))'
      write(fmt(17:17),'(I1)') npcol
      write(fmt(32:32),'(I1)') npcol
      if (ipr .ge. 50) then
        write(stdo,'('' '')')
        write(stdo,fmt) (slabl(is),is=1,nspec)
      endif
      fmt = '(1x,''init rmt'',1x,6f9.5:/(10x,6f9.5))'
      write(fmt(19:19),'(I1)') npcol
      write(fmt(31:31),'(I1)') npcol
C     print *, fmt
      if (ipr .ge. 50) write(stdo,fmt) (wsr(is),is=1,nspec)

      do  nloop = 1, nspec+1
        amax = 9d9
        do  is = 1, nspec
          if (lock(is) .ne. 2) lock(is) = 0
        enddo

C   --- Lock radii of those spheres with maximum allowed overlap ---
C       and unlock those with radii exceeding maximum allowed.
        do  is = 1, nspec
          wsri = wsr(is)

C     ... If overlap criterion exactly satisfied for any connecting vector,
C         set lock for this species
          do  jb = 1, nrspec(is)
            ib = iclbsj(is,ips,nbas,jb)
            do  kpr = ntab(ib)+2, ntab(ib+1)
              kb = iax(2,kpr)
              k1 = iax(3,kpr)
              k2 = iax(4,kpr)
              k3 = iax(5,kpr)
              ks = ips(kb)
C           ip selects which ommax to use (A-A, A-E, or E-E)
              ip = 2
              if (idnint(z(is)).ne.0 .and. idnint(z(ks)).ne.0) ip=1
              if (idnint(z(is)).eq.0 .and. idnint(z(ks)).eq.0) ip=3
              opo1 = 1+omax1(ip)
              omo2 = 1-omax2(ip)
              wsrk = wsr(ks)
              rik = dsqrt(drr2(plat,bas(1,ib),bas(1,kb),k1,k2,k3,dr))
              rik = rik*alat
C           Set lock if any of these conditions are met:
C             wi + wk - rik = o1*rik
C             wi + wk - rik = o2*wi
C             wi + wk - rik = o2*wk
C             wi = wsmax
              locki = lock(is)
              if (dabs(opo1*rik-wsri-wsrk) .lt. tiny) locki=1
              if (dabs(rik-omo2*wsri-wsrk) .lt. tiny) locki=1
              if (dabs(rik-wsri-omo2*wsrk) .lt. tiny) locki=1
              if (wsmax.gt.0 .and. dabs(wsri-wsmax).lt.tiny) locki=1
              lock(is) = max(locki,lock(is))
            enddo
          enddo
C         print *,'is,lock0=',is,lock(is)

C     ... If overlap criterion exceeded for any connecting vector,
C         unset lock for this species
          if (lock(is) .ne. 2) then
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
C           ip selects which ommax to use (A-A, A-E, or E-E)
                ip = 2
                if (idnint(z(is)).ne.0 .and. idnint(z(ks)).ne.0) ip=1
                if (idnint(z(is)).eq.0 .and. idnint(z(ks)).eq.0) ip=3
                opo1 = 1+omax1(ip)
                omo2 = 1-omax2(ip)
C           Unset lock if any of these conditions are exceeded
C             wi + wk - rik < o1*rik
C             wi + wk - rik < o2*wi
C             wi + wk - rik < o2*wk
C             wi > wsmax
                if (opo1*rik-wsri-wsrk .lt. -tiny) lock(is)=0
                if (rik-omo2*wsri-wsrk .lt. -tiny) lock(is)=0
                if (rik-wsri-omo2*wsrk .lt. -tiny) lock(is)=0
                if (wsmax.gt.0 .and. wsri-wsmax.gt.tiny) lock(is)=0
              enddo
            enddo
          endif
        enddo
C       print *, 'lock', (lock(is), is=1,nspec)

C  --- Find amax=largest allowed scaling for unlocked species ---
        do  is = 1, nspec
          if (lock(is) .eq. 0) then
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
                if (idnint(z(is)).ne.0 .and. idnint(z(ks)).ne.0) ip=1
                if (idnint(z(is)).eq.0 .and. idnint(z(ks)).eq.0) ip=3
                opo1 = 1+omax1(ip)
                omo2 = 1-omax2(ip)
                if (dabs(rik-riko).gt.tiny .or. kco.ne.ks) then
                  wsrk = wsr(ks)
                  riko = rik
                  kco = ks
                  amax1 = 9d9
                  amax2 = 9d9
                  amax3 = 9d9
                  amax4 = 9d9
C               If second site ks is locked, can only scale site i
                  if (lock(ks) .ne. 0) then
                    amax1 = (opo1*rik-wsrk+gw)/(wsri+gw)
                    if (omo2 .gt. 0d0)
     .              amax2 = (rik-wsrk+omo2*gw)/(omo2*wsri+omo2*gw)
                    amax3 = (rik-omo2*wsrk+gw)/(wsri+gw)
C               If neither site is locked, both will scale
                  else
                    amax1 = (opo1*rik+gw+gw)/(wsri+wsrk+gw+gw)
                    if (wsrk+omo2*wsri+(1+omo2)*gw .gt. 0d0)
     .              amax2 = (rik+(1+omo2)*gw)/
     .              (wsrk+omo2*wsri+(1+omo2)*gw)
                    if (wsri+omo2*wsrk+(1+omo2)*gw .gt. 0d0)
     .              amax3 = (rik+(1+omo2)*gw)/
     .              (wsri+omo2*wsrk+(1+omo2)*gw)
                  endif
                  if (wsmax .gt. 0) then
                    amax4 = wsmax/wsri
                  endif
                  amax = dmin1(amax,amax1,amax2,amax3,amax4)
C               print *, 'kc,kbas,kpr,amax',ks,kb,kpr,amax
                endif
              enddo
            enddo
          endif
        enddo
        bmax = (1d0-1/amax)*gw

C   ... Determine what new volume will be after scaling with a,b
        vola = 0d0
        volb = 0d0
        do  ib = 1, nbas
          is = ips(ib)
          if (lock(is) .eq. 0) then
            volb = volb + (amax*wsr(is)+bmax)**3
          else
            vola = vola + wsr(is)**3
          endif
        enddo
        vola = vola * fpi3
        volb = volb * fpi3

C   --- Case scaling will lead to new volume > final volume ---
        if (vol*volfac .lt. vola+volb) then
          call dpzero(dm,4)
          do  ib = 1, nbas
            is = ips(ib)
            if (lock(is) .eq. 0) then
              a = wsr(is)+gw
C             For numerical reasons distinguish cases
              if (dabs(gamma) .gt. 1d0) then
                b = wsr(is)
              else
                b = -gw
              endif
C             fpi3*dm(0) = sum of sphere volumes not constrained
              dm(0) = dm(0) +       b*b*b
              dm(1) = dm(1) + 3d0 * a*b*b
              dm(2) = dm(2) + 3d0 * a*a*b
              dm(3) = dm(3) +       a*a*a
            endif
          enddo
C         Vconst + a**3*Vuncst = Vtarget; Vconst + Vuncst = vola
C         If Vconst > Vtarget, constraint cannot be satisified
          Vconst = vola - fpi3*dm(3)
          if (Vconst .gt. vol*volfac) then
            call fexit2(-1,111,' Exit -1 : SCLWSR: constrained '//
     .      'sphere vol (%;0d) exceeds target vol=%;0d',
     .      Vconst,vol*volfac)
          endif
          if (dabs(dm(3)) .gt. tiny) then
            r =  dm(2) / dm(3)
            s =  dm(1) / dm(3)
            t = (dm(0) - (vol*volfac-vola)/fpi3) / dm(3)
            p = s - r*r/3d0
            q = 2d0*r*r*r/27d0 - r*s/3d0 + t
            d = p*p*p/27d0 + q*q/4d0
            u = (dsqrt(d)-q/2d0)**(1d0/3d0)
            v = -p/u/3d0
            x = u+v-r/3d0
            if (dabs(gamma) .gt. 1d0) then
              amax = x+1d0
              bmax = x*gw/amax
            else
              amax = x
              bmax = (1d0-1d0/amax)*gw
            endif
            if (ipr .ge. 100) then
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

C   --- Scale unlocked spheres by w <- a*w+b ---
        fin = .true.
        do  is = 1, nspec
          dscl(is) = 1d0
          if (lock(is) .eq. 0) then
            dsclmx  = amax + bmax/wsr(is)
            wsr(is) = dsclmx*wsr(is)
            dscl(is)= dsclmx
            fin = .false.
          endif
        enddo

C       ratio = new volume / old volume
        ratio = volsph(nspec,nrspec,wsr)/vol
C       fin=T if volume change is small or all spheres are locked
        fin = fin .or. dabs(ratio-volfac).lt.tiny .or. nloop.eq.nspec+1

        fmt = '(1x,''iter:'',i3,1x,6f9.5:/(10x,6f9.5))'
C       print *, fmt
        write(fmt(19:19),'(I1)') npcol
        write(fmt(31:31),'(I1)') npcol
C       print *, fmt
        if (ipr .ge. 50) write(stdo,fmt) nloop,(dscl(is),is=1,nspec)

C   ... Last iteration
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
                if (idnint(z(is)).ne.0 .and. idnint(z(ks)).ne.0) ip=1
                if (idnint(z(is)).eq.0 .and. idnint(z(ks)).eq.0) ip=3
                wsrk = wsr(ks)
                rik = dsqrt(drr2(plat,bas(1,ib),bas(1,kb),k1,k2,k3,dr))
                rik = rik*alat
                dovlap = wsri+wsrk-rik
                dovl1(ip) = dmax1(dovl1(ip),dovlap/rik)
                dovl2(ip) = dmax1(dovl2(ip),dovlap/wsri,dovlap/wsrk)
C              if (dovl1(ip).gt.omax1(ip)+tiny .or.
C     .            dovl2(ip).gt.omax2(ip)+tiny)
C     .          write(stdo,'('' SCLWS2: warning dovl>omax'')')
              enddo
            enddo
          enddo
          if (ipr .ge. 50) then
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

      subroutine pairs(nbas,nbasp,alat,plat,rmax,baspp,ipsp,nd,iltab,
     .pltab,nttab,iv_a_ontab,iv_a_oiax,mxcsiz)
C- Allocate memory for and create neighbor table for a crystal
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis (input)
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas !+ nbas(left bulk) + nbas(right bulk)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat (input)
Ci   rmax  :maximum range for connecting vector, in a.u.
Ci          All connecting vectors with length < rmax(i)+rmax(j)
Ci          are retained.  rmax may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ipsp(1);
Ci          see description of ipsp
Ci   baspp :basis vectors, doubly padded for planar geometry
Ci   ipsp  :index to which species each site belongs, for padded basis;
Ci          identifies which rmax is associated with each site. NB:
Ci          ipsp(1) = -1 => rmax is a global scalar, independent of site
Ci          ipsp(1) =  0 => rmax is site-, not species-dependent.
Ci          In either of these cases, ipsp is not used.
Ci   nd    :number of dimensions for which periodic boundary conditions
Ci          are used
Ci   iltab :iltab<0, has no effect.  Otherwise, see pltabp.
Ci   pltabp:include only pairs for which pltabp(jb)-pltabp(ib) <= iltab
Ci   mxcsiz:if nonzero, use in place of internal formula for mxnbr
Co Outputs
Co   nttab   :total number of pairs in neighbor table and iax
Co   w(ontab):ntab array; see pairc below where it is generated
Co   w(oiax) :iax array; see pairc, below where it is generated
Co   mxcsiz  :size of the largest cluster encountered
C ----------------------------------------------------------------------
      implicit none
      integer:: nbas , nbasp , ipsp(1) , pltab(1) , nttab , nd , iltab,i_data_size
      integer,allocatable :: iv_a_ontab(:)
      integer,allocatable :: iv_a_oiax(:,:)
      double precision alat,plat(9),rmax(1),baspp(3,1)
      double precision avw,avwsr,vol
      integer:: modep(3) , nbaspp , mxcsiz , mxnbr , i , niax , isw
      real(8) ,allocatable :: wk_rv(:)
      integer,allocatable:: iv_a_tmp(:,:)
      parameter (niax=10)
C ... Set up input for call to pairc
      nbaspp = nbas !2*nbasp - nbas
C ... Estimate an upper bound to the size of the neighbor table
      avw = avwsr(plat,alat,vol,nbas)
      mxnbr = 3*(2*rmax(1)/avw)**3*nbasp
      if (mxcsiz .gt. 0) mxnbr = mxcsiz
      if (allocated(iv_a_ontab)) deallocate(iv_a_ontab)
      allocate(iv_a_ontab(abs(nbasp+1)))
      if (allocated(iv_a_oiax)) deallocate(iv_a_oiax)
      allocate(iv_a_oiax(niax,mxnbr))
      iv_a_oiax=0
      allocate(wk_rv(3*mxnbr))
      do  10  i = 1, 3
        modep(i) = 2
        if (i .gt. nd) modep(i) = 0
   10 continue

C ... This makes the neighbor table
      nttab = mxnbr
      isw = 0
      call pairc ( 1 , nbasp , nbaspp , modep , isw , ipsp , alat ,
     .  plat , baspp , baspp , rmax , iltab , pltab , nttab , iv_a_ontab
     .  , iv_a_oiax , wk_rv , mxcsiz )
C ... Allocate iax to proper size
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

      subroutine pairc(ib1,ib2,nbasp,mode,isw,ips,alat,plat,pos,ctr,
     .range,iltab,pltabp,nttab,ntab,iax,rtab,mxcsiz)
C- Make a neighbor table (crystal version)
C ----------------------------------------------------------------
Ci Inputs:
Ci  ib1,ib2:range of sites for which to generate tables
Ci   nbasp :the size of the basis, plus possible extensions.
Ci          Usually nbasp=nbas, but will differ in special
Ci          cases, such as having padding sites to extend
Ci          to a semi-infinite geometry.
Ci   mode:  vector of length 3 governing how pos shortened (see shorps)
Ci   isw:   1's digit fixes how range is calculated.
Ci           0: vector length must be < range(i)+range(j)
Ci           1: include all connecting vecs w/ r < range(i)
Ci         10's digit sets what part of iax table is not calculated
Ci           1: do not calculate iax(6)
Ci              (may be needed when ctr and pos are different)
Ci           2: calculate only iax(1..5)
Ci   ips   :index to which species each site belongs, for padded basis;
Ci          identifies which rmax is associated with each site. NB:
Ci          ips(1) = -1 => rmax is a global scalar, independent of site
Ci          ips(1) =  0 => rmax is site-, not species-dependent.
Ci          In either of these cases, ips is not used.
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat (input)
Ci   pos   :site positions (doubly padded for planar geometry)
Ci   ctr   :ctr(1..3,ib) is the effective center of the cluster
Ci          associated with site ib for purposes of computing distance
Ci          pos(jb)-ctr(ib).  May point to the same address space as pos
Ci   range :maximum range for connecting vector, in a.u..
Ci          This quantity may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ips(1);
Ci          see description of ips.  Precisely what meaning range has
Ci          depends on mode and isw.
Ci   iltab :iltab<0, has no effect.  Otherwise, see pltabp.
Ci   pltabp:include only pairs for which pltabp(jb)-pltabp(ib) <= iltab
Ci   nttab :maximum dimension of iax table; used to guard against
Ci          generated table size exceeding dimension of iax.
Co Outputs:
Co   nttab    :total number of pairs generated
Co   iax      :neighbor table containing information about each pair ip
Co             For each pair ip, information is contained in iax(*,ip).
Co             as described below.  iax is ordered grouped by the basis
Co             atoms, so that all pairs connected to site ib are grouped
Co             together.  For each pair ip, iax(*,ip) contains:
Co   iax(1)   :site index to basis atoms ib=source;
Co             all pairs with common ib are contiguous
Co   iax(2)   :site index to jb=field of each pair
Co   iax(3..5):multiples of plat added the difference in site positions
Co             that connect the pair.
Co   iax(6)   :index to conjugate (jb,ib) pair matching (ib,jb)
Co             NB: no matching pairs outside (ib1..ib2) can be found.
Co   iax(7)   :permutation index ordering cluster by increasing
Co             effective site index; see ppair4.f
Co   iax(8)   :left untouched by pairc
Co   iax(9)   :left untouched by pairc
Co   iax(10)  :effective site index; see siteid.f
Co   ntab     :ntab(ib)=number of pairs in iax table preceding ib
Co             ntab is created for ib1:ib2+1.
Co   rtab     :rtab(1..3,ip) = pos(jb)-ctr(ib) for pair ip
Co   mxcsiz   :the largest cluster encountered
Cr Remarks
Cr   For each site ib=ib1..ib2, pairc finds all connecting vectors
Cr   for a lattice of points with periodic boundary conditions in
Cr   1, 2, or 3 dimensions, within a specified range of site ib.
Cr   The range can be defined in various ways, depending on isw.
Cu Updates
Cu   23 Apr 02 added option to make only iax(1..5) (isw=20)
C ----------------------------------------------------------------
      implicit none
      integer ib1,ib2,nbasp,mode(3),isw,nttab,niax,ips(nbasp),
     .ntab(ib1:ib2+1),iltab,pltabp(nbasp),mxcsiz
      parameter (niax=10)
      integer iax(niax,1),n,itmp,ixo
      double precision alat,plat(3,3),pos(3,nbasp),ctr(3,ib2),
     .range(1),rtab(3,1)
C Local variables
      integer:: ib , is , jb , mtab , i , moder , mode2(3) , nlat , 
     .mxntab , nsite
      real(8) ,allocatable :: wk1_rv(:)
      real(8) ,allocatable :: wk2_rv(:)
      real(8) ,allocatable :: wk3_rv(:)
      real(8) ,allocatable :: pos_rv(:)
      real(8) ,allocatable :: lat_rv(:)
      real(8) ,allocatable :: ctr_rv(:)

      double precision r1,rr,qlat(3,3),p0(3)
      integer ::iwdummy


C     call tcn('ppair1')

C --- Setup ---
      nsite = ib2-ib1+1
      mxntab = nttab
      moder = mod(isw,10)
      do  3  i = 1, 3
        mode2(i) = mode(i)
        if (mode2(i) .eq. 1) mode2(i) = 0
    3 continue
C ... Make r1 = 2*maximum range
      r1 = range(1)
      if (ips(1) .ge. 0) then
        do  5  ib = 1, nbasp
          is = ib
          if (ips(1) .gt. 0) then
            is = ips(ib)
          endif
          r1 = max(r1,range(is))
    5   continue
      endif
      if (moder .eq. 0) r1 = 2*r1
      r1 = 2*r1
C ... List of lattice vectors to add to pos(ib)-pos(jb)
      call xlgen ( plat , r1 / alat , 0d0 , 0 , 20 , mode , i , iwdummy 
     . )

      allocate(lat_rv(3*i))

      call xlgen ( plat , r1 / alat , 0d0 , i , 0 , mode , nlat , lat_rv 
     .)

C ... qlat = (plat^-1)^T so that qlat^T . plat = 1
      call mkqlat(plat,qlat,rr)
C ... Save true pos in opos
C     and ctr in octr in case same address space used for ctr
      allocate(pos_rv(3*nbasp))

      call dpcopy ( pos , pos_rv , 1 , 3 * nbasp , 1d0 )

      allocate(ctr_rv(3*nsite))

      call dpcopy ( ctr , ctr_rv , 1 , 3 * nsite , 1d0 )


C --- For each ib, find all pairs for which dr < range ---
      nttab = 1
      ntab(ib1) = 0
      mtab = 1
      mxcsiz = 0
      do  10  ib = ib1, ib2
        r1 = range(1)
        if (ips(1) .ge. 0) then
          is = ib
          if (ips(1) .gt. 0) then
            is = ips(ib)
          endif
          r1 = range(is)
        endif
C --- Shorten all pos relative to ctr(ib) ---
C ... Make pos-ctr(ib)
        call dpcopy ( pos_rv , pos , 1 , 3 * nbasp , 1d0 )
        call dpcopy ( ctr_rv , ctr , 1 , 3 * nsite , 1d0 )
        do  12  i = 1, 3
          p0(i)  = ctr(i,ib)
          do  14  jb = 1, nbasp
            pos(i,jb) = pos(i,jb) - p0(i)
   14     continue
   12   continue
C ... Shorten pos-ctr(ib)
        call shorps(nbasp,plat,mode2,pos,pos)
C ... Undo shift -ctr(ib) to restore shortened pos to absolute pos
        do   jb = 1, nbasp
        do   i = 1, 3
          pos(i,jb) = pos(i,jb) + p0(i)
        enddo
        enddo
C --- Find all sites in range of ctr ---
        call ppair2 ( nbasp , iltab , pltabp , moder , alat , qlat , 
     .  pos , p0 , range , ips , rtab , ib , r1 , nlat , lat_rv , 
     .  pos_rv , mxntab , nttab , iax )
C --- Sort table by increasing length ---
        call ppair3 ( nttab - mtab , iax(1, mtab ), rtab(1 , mtab))! , wk1_rv , wk2_rv , wk3_rv )
C --- Cleanup for this ib ---
        mtab = nttab
        ntab(ib+1) = nttab-1
        mxcsiz = max(mxcsiz,ntab(ib+1)-ntab(ib))
   10 continue
      nttab = nttab-1
C --- Restore original pos,ctr ---
      call dpcopy ( pos_rv , pos , 1 , 3 * nbasp , 1d0 )

      call dpcopy ( ctr_rv , ctr , 1 , 3 * nsite , 1d0 )

      if (allocated(ctr_rv)) deallocate(ctr_rv)
      if (allocated(pos_rv)) deallocate(pos_rv)


C --- Fill out iax table ---
      call ppair1(isw,ib1,ib2,nbasp,ips,alat,plat,pos,range,
     .nttab,ntab,iax,mxcsiz)

      if (allocated(lat_rv)) deallocate(lat_rv)

      end
      subroutine ppair1(isw,ib1,ib2,nbasp,ips,alat,plat,pos,range,
     .nttab,ntab,iax,mxcsiz)
C- Fill out parts of the aix table
C ----------------------------------------------------------------
Ci  Inputs
Ci   isw   :1's digit fixes how range is calculated.
Ci           0: vector length must be < range(i)+range(j)
Ci           1: include all connecting vecs w/ r < range(i)
Ci         :10's digit sets what part of iax table is calculated
Ci           0: make iax(6),iax(7),iax(10)
Ci           1: make iax(7),iax(10)
Ci           2: no change to iax: printout only
Ci           4: just make iax(6)
Ci   ib1   :fill out iax table for pairs ntab(ib1)+1..ntab(ib2)
Ci   ib2   :fill out iax table for pairs ntab(ib1)+1..ntab(ib2)
Ci   nbasp :size of padded basis (not needed)
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   range :maximum range for connecting vector, in a.u..
Ci          This quantity may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ips(1);
Ci          see description of ips.  See 1s digit of isw for
Ci          how range is used.
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   mxcsiz:maximum cluster size (for printout only)
Co  Outputs
Co   iax(6)   :index to conjugate (jb,ib) pair matching (ib,jb)
Co             NB: only matching pairs within site list can be found.
Co   iax(7)   :permutation index ordering cluster by increasing
Co             effective site index; see ppair4.f
Co   iax(10)  :effective site index; see siteid.f
C ----------------------------------------------------------------
      implicit none
      integer isw,ib1,ib2,nbasp,nttab,niax,ips(nbasp),ntab(ib1:ib2)
      parameter (niax=10)
      integer iax(niax,1),mxcsiz
      double precision alat,plat(3,3),pos(3,19),range(1)
C Local variables
      integer:: ib , is , jb , js , ipr , i , j , moder , it , jt , 
     .iprint , lgunit , stdo , nsite , isw1
      real(8) ,allocatable :: pos_rv(:)
      integer ,allocatable :: iwk_iv(:)

      double precision r1,r2,rr,rcut,vlat(3),drr2,tol
      parameter (tol=1d-5)

      isw1 = mod(isw/10,10)
      ipr = iprint()
      moder = mod(isw,10)
      nsite = ib2-ib1+1
      stdo = lgunit(1)
      if (isw1 .eq. 2) goto 80
      if (isw1 .eq. 4) goto 71

C --- Set iax(7) to sort this cluster ---
      call ppair5(ib1,ib2,plat,pos,tol,ntab,iax)

C --- For each pair, find matching pair, store in iax(6) ---
   71 continue
      do  74  it = 1,  nttab
        iax(6,it) = 0
   74 continue
      if (mod(isw1,2) .eq. 0) then
        do 170  ib = ib1, ib2
        do  70  it = ntab(ib)+1, ntab(ib+1)
          if (iax(6,it) .ne. 0) goto 70
          jb = iax(2,it)
C   ... No matching pair for padded sites
          if (jb .lt. ib1 .or. jb .gt. ib2) goto 70
          do  72  jt = ntab(jb)+1, ntab(jb+1)
C         if (iax(2,jt) .eq. ib) then
C         endif
            if (iax(2,jt) .eq. ib .and.
     .        iax(3,it) .eq. -iax(3,jt) .and.
     .        iax(4,it) .eq. -iax(4,jt) .and.
     .        iax(5,it) .eq. -iax(5,jt))  then
              iax(6,it) = jt
              iax(6,jt) = it
              goto 73
            endif
   72     continue
          call fexit3(-1,1,' Exit -1 pairc: cannot find pair'//
     .      ' matching sites (%i,%i), pair %i',ib,jb,it-ntab(ib))
   73     continue
   70   continue
  170   continue
      endif
      if (isw1 .eq. 4) return

C ... Assign a unique id for every different site in the cluster table
      allocate(iwk_iv(nttab))
      allocate(pos_rv(3*nttab))
      call siteid ( iax , nsite , ntab , plat , pos , pos_rv , iwk_iv 
     ., i )
      if (allocated(pos_rv)) deallocate(pos_rv)
      if (allocated(iwk_iv)) deallocate(iwk_iv)


C --- Printout ---
   80 if (ipr .lt. 30) goto 91
      if (ipr .le. 40) write(stdo,'(1x)')
      if (ipr .gt. 40) write(stdo,332)
  332 format(/'  ib  jb',9x,'--- r(jb)-r(ib) ---',10x,
     .'d       -x-plat-  map ord  id')
      i = 0
      do  90  it = 1, nttab
        ib = iax(1,it)
        jb = iax(2,it)
        rr = dsqrt(drr2(plat,pos(1,ib),pos(1,jb),
     .  iax(3,it),iax(4,it),iax(5,it),vlat))
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
          if (alat .ne. 1)
     .    write(stdo,345) ib,ntab(ib+1)-ntab(ib),rcut/alat,rcut
          if (alat .eq. 1) write(stdo,345) ib,ntab(ib+1)-ntab(ib),rcut
  345     format(' pairc, ib=',i3,':',i4,' neighbors in range',2f7.3)
        endif
        i = ib
        if (ipr .gt. 40) write(stdo,334) iax(1,it),iax(2,it),
     .  (vlat(j),j=1,3), rr, (iax(j,it), j=3,7),iax(10,it)
  334   format(i4,i4,3f11.6,f9.4,3x,3i3,i5,2i4)
   90 continue
   91 if (ipr .ge. 20) write(stdo,
     .'('' pairc:'',i8,'' pairs total'',i10,'' is max cluster size'')')
     .nttab, mxcsiz
      end

      subroutine ppair2(nbas,iltab,pltabp,moder,alat,qlat,pos,ctr,range,
     .ips,rtab,ib,r1,nlat,lat,trupos,mxntab,nttab,iax)
C- Kernel of pairc to find all sites in range of ctr
      implicit none
      integer nbas,ib,iltab,ips(nbas),pltabp(nbas),niax,nlat,moder,
     .mxntab,nttab
      parameter (niax=10)
      integer iax(niax,1)
      double precision alat,ctr(3),pos(3,nbas),range(nbas),rtab(3,1)
      double precision qlat(3,3),trupos(3,nbas),lat(3,*),r1
C Local variables
      integer i,ilat,jb,js
      double precision r2,rr,rcut,vlat(3),xx,rcutba,dpos(3)

      do  20  jb = 1, nbas
        if (iltab .ge. 0) then
          if (abs(pltabp(jb)-pltabp(ib)) .gt. iltab) goto 20
        endif
        r2 = range(1)
        if (ips(1) .ge. 0) then
          js = jb
          if (ips(1) .gt. 0) then
            js = ips(jb)
          endif
          r2 = range(js)
        endif
        if (moder .eq. 0) rcut = r1+r2
        if (moder .eq. 1) rcut = r1
        rcutba = (rcut / alat)**2
        dpos(1) = pos(1,jb)-ctr(1)
        dpos(2) = pos(2,jb)-ctr(2)
        dpos(3) = pos(3,jb)-ctr(3)

C   --- For each (ib,jb,ilat), do ---
        do  30  ilat = 1, nlat

          if (nttab .gt. mxntab) call rxi(
     .    'pairc: table exceeds input maximum size,',mxntab)

C ...   Add to list if connecting vector within range
          rtab(1,nttab) = dpos(1) + lat(1,ilat)
          rtab(2,nttab) = dpos(2) + lat(2,ilat)
          rtab(3,nttab) = dpos(3) + lat(3,ilat)
          rr = rtab(1,nttab)**2+rtab(2,nttab)**2+rtab(3,nttab)**2

c*        call awrit5('try ib,jb,ilat= %i %i %i rr=%;4d: %l',' ',80,
c*     .    6,ib,jb,ilat,rr,rr.lt.rcut)

C   --- Add to iax table if this pair in range ---
          if (rr .lt. rcutba) then

C     ... vlat += shortening vector
            do  32  i = 1, 3
              rtab(i,nttab) = alat*rtab(i,nttab)
C           rtab(i,nttab) = alat*(rtab(i,nttab)+ctr(i)-pos(i,ib))
              vlat(i) = lat(i,ilat) + pos(i,jb) - trupos(i,jb)
   32       continue

C     ... iax table for this pair
            iax(1,nttab) = ib
            iax(2,nttab) = jb
            do  33  i = 1, 3
              xx = vlat(1)*qlat(1,i)+vlat(2)*qlat(2,i)+vlat(3)*qlat(3,i)
              iax(2+i,nttab) = nint(xx)
   33       continue
            nttab = nttab+1

          endif

   30   continue
   20 continue
      end

      subroutine ppair3(nttab,iax,rtab) !,iwk,iwk2,rwk)
C- Sort neighbor table by distance
      implicit none
      integer nttab,niax,iwk2(nttab),i,j,k
      parameter (niax=10)
      integer iax(niax,nttab),iwk(niax,nttab)
      double precision rtab(3,nttab),rwk(3,nttab)

      do  10  i = 1, nttab
        rwk(1,i) = rtab(1,i)
        rwk(2,i) = rtab(2,i)
        rwk(3,i) = rtab(3,i)
        do  12  k = 1, niax
          iwk(k,i) = iax(k,i)
   12   continue
   10 continue
      call dvshel(3,nttab,rtab,iwk2,11)
      do  20  i = 1, nttab
        j = iwk2(i)+1
        rtab(1,i) = rwk(1,j)
        rtab(2,i) = rwk(2,j)
        rtab(3,i) = rwk(3,j)
        do  22  k = 1, niax
          iax(k,i) = iwk(k,j)
   22   continue
   20 continue
      end
      
      subroutine ppair4(iclus,nclus,plat,pos,ctr,iwk,rtab,tol,iax)
C- Sort cluster by increasing (x,y,z) relative to its center
C ----------------------------------------------------------------
Ci Inputs
Ci   iclus,nclus: sort iax(iclus..nclus)
Ci   plat :primitive lattice vectors
Ci    pos :basis vectors
Ci    ctr :cluster origin:does not affect the ordering, but shifts rtab
Ci    iwk :integer work array of length nclus-iclus+1
Ci    tol :tolerance to which positions are considered coincident
Ci         tol<0 => sort iax by iax(1..5)
Co Outputs
Co   iax(7,iclus..nclus) orders the cluster by increasing (x,y,z)
Co         (or increasing iax(1..5) if tol < 0
Co   rtab  :connecting vectors rtab(1..3,ip) = pos(jb)-ctr
Co          for pair ip and jb=iax(2,ip)
Cr Remarks
Cr  Each cluster is sorted by increasing (x,y,z),
Cr  sorted by x first, then by y, then by z, thus guaranteeing that
Cr  all sites common to any pair of clusters are ordered the same.
C ----------------------------------------------------------------
C     implicit none
      integer iclus,nclus,niax,iwk(15)
      parameter (niax=10)
      integer iax(niax,1)
      double precision plat(3,3),pos(3,1),ctr(3),rtab(3,29),tol
      integer ic,jb,ic0,ix,ia2,i,j,k
C Local variables
      double precision dx
C     integer jx
C     double precision wk2(3,nclus*3)
      dx(ia2,i,j,k) = pos(ix,ia2) +
     .plat(ix,1)*i + plat(ix,2)*j + plat(ix,3)*k

      ic0 = 0
      do  12  ic = iclus, nclus
        jb = iax(2,ic)
        ic0 = ic0+1
        do  14  ix = 1,3
          rtab(ix,ic0) = dx(jb,iax(3,ic),iax(4,ic),iax(5,ic)) - ctr(ix)
   14   continue
   12 continue

      if (tol .lt. 0) then
        call ivheap(niax,nclus-iclus+1,iax(1,iclus),iwk,1)
      else
        call dvheap(3,nclus-iclus+1,rtab,iwk,tol,1)
      endif
      do  20  ic = iclus, nclus
        iax(7,ic) = iwk(ic-iclus+1)
   20 continue
      end
      
      subroutine ppair5(ib1,ib2,plat,pos,tol,ntab,iax)
C- Sort a range of clusters according to tol
C ----------------------------------------------------------------------
Ci Inputs
Ci  ib1,ib2:range of clusters to sort
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   tol   :tolerance; see ppair4
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Co Outputs
Co   iax   :iax(7) is set to order cluster; see ppair4.
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer niax,ib1,ib2
      parameter (niax=10)
      integer iax(niax,1),ntab(ib2)
      double precision plat(3,3),pos(3,*),tol
      integer:: ib , nttab
      integer ,allocatable :: iwk1_rv(:)
      real(8) ,allocatable :: wk2_rv(:)
C --- Set iax(7) to sort this cluster ---
      do  10  ib = ib1, ib2
        nttab = ntab(ib+1)-ntab(ib)
        allocate(iwk1_rv(nttab*2))
        allocate(wk2_rv(nttab*3))
        call ppair4 ( ntab ( ib ) + 1 , ntab ( ib + 1 ) , plat , pos 
     .  , pos ( 1 , ib ) , iwk1_rv , wk2_rv , tol , iax )
        deallocate(iwk1_rv,wk2_rv)
   10 continue
      end
      end module m_lmaux

