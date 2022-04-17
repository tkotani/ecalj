      module m_lmaux
      public:: lmaux
      contains
      
      subroutine lmaux()        !main part of lmchk
      use m_mksym,only: ctrl_nclass,rv_a_oclabl,iv_a_oics,iv_a_oipc
      use m_lmfinit,only: iv_a_oips,str_mxnbr,str_rmax,ctrl_nbas,ctrl_nspec,ctrl_nspin,
     &     ctrl_nl,ctrl_modep,ctrl_omax1,ctrl_omax2,ctrl_wsrmax,slabl,sspec=>v_sspec,
     &     lat_avw,lat_alat
      use m_lattic,only: lat_nkd
      use m_lattic,only: lat_nkq
      use m_struc_def  !Cgetarg
      use m_pairs
      use m_ovmin , only: ovmin
      use m_lattic,only:lat_plat,rv_a_opos
C- lmchk 
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam:name of main program
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nclass nl nspec nspin lncol modep lpgf lrs,7
Ci                omax1 omax2 wsrmax
Ci     Stored:    *
Ci     Passed to: asars supot asamad asvsph aioxtn
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: oeula neula qss
Ci     Stored:    *
Ci     Passed to: asvsph
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: opnu oqnu opp oqc oqt oaamom orhrmx oves ovrmax
Ci     Stored:    *
Ci     Passed to: asars supot asamad asvsph
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: avw alat plat plat2 nkd odlv nkq oqlv
Ci     Stored:    *
Ci     Passed to: asars supot asamad plana asvsph aioxtn
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: idmod lmxa z rham mxcst rmt a nr
Ci     Stored:    z rmt
Ci     Passed to: asars spackv asamad plana pp2enu asvsph shopol aioxtn
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    pos spec
Ci     Passed to: asars iosits spackv plana aioxtn
Ci   sstr  :struct for parameters for screened strux; see routine ustr
Ci     Elts read: mxnbr rmax
Ci     Duplicate: mxnbr
Ci     Stored:    *
Ci     Passed to: *
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: ohave oics oclabl onrc npadl npadr nclasp opgfsl
Ci                ormax opos oipc oips
Ci     Stored:    ormax
Ci     Passed to: asamad asvsph aioxtn
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: asars
Ci   slabl :vector of species labels
Ci   sstrn :struct for global strings
Ci   mode  :a compound of bits, which are independent of each other
Ci         :  2**0 Show neighbors
Ci         :  2**1 Plane analysis
Ci         :  2**2 Generate input to xbs program
Ci         :  2**3 Shift moments, pp's to new linearization energy
Ci         :  2**4 Interpolate core to another mesh
Ci         :  2**5 Display poles of potential functions
Ci         :  2**6 Import data from other formats
Ci         :  2**7 Find empty spheres
Cu Updates
Cu   12 Aug 08 (L. Ke) empty sphere finder
Cu   04 Nov 04 Upgrade of rsta editor
Cu   26 Jan 03 Call to angtab changed
Cu   17 May 02 Modified MT radii scaling to lower priority for E.S.
Cu   23 Apr 02 Added option (--getwsr) to find MT radii
Cu   01 Mar 02 Updated Import data mode
Cu   05 Oct 01 Adapted mode 2**3 to work with lm v6.11
Cu   24 Nov 97 changed ovmin to run quickly
C ----------------------------------------------------------------------
      implicit none
c      character prgnam*8
c      data prgnam /'LMCHK'/
C ... Passed parameters
      integer:: mode=1 !,wksize
c      character*(*)  prgnam*8 !sstrn,
c      character*8 slabl(*)
c      type(s_array)::sarray
c      type(s_ctrl)::sctrl
c      type(s_ham)::sham
c      type(s_pot)::spot
c      type(s_lat)::slat
c      type(s_spec)::sspec(*)
c      type(s_site)::ssite(*)
c      type(s_str)::sstr
c      type(s_bz)::sbz
c      type(s_site),allocatable::ssite2(:)

C ... Local parameters
      character*120 outs,fnam*8
      integer NULLI
      logical cmdopt,T,F,swtmp
      parameter (T=.true., F=.false., NULLI=-99999)
      integer getdig,i,ip,j,k,m,ifi,iprint,lpbc,lgunit,
c !lncol,
     .nbas,nclasp,nclass,nl,nlspc,nsp,modep(3),parg,nbasp,
     .nbaspp,nkd,nkq,nspec,neul,nc,mxcsiz,nttab,igets, !npadl,npadr,
     .iosits,cmplat,ngrp,ival,stdo,irs(5),nclspp,bitand,igetss,
     .ngmx,nsgrp
      integer:: oeold  , olmx , opold , owk2 , orham , oamsh 
     .                       , onrmsh , oalpha , onpr , os , ormx , oip , opgfsl , mxclas
c      real(8),allocatable:: eula(:,:)
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
C ... Heap

cki      integer:: ifx,w_dummy=1
      integer:: ifx,w_dummy(1)=1
ckino 2012.Feb.13 needless      integer,allocatable:: ics(:)
      integer,allocatable:: lmxa(:)
c      real(8),allocatable::idmod(:,:), z(:),amsh(:),rmax(:), nrmsh(:)
      real(8),allocatable:: z(:),rmax(:)

      print *,' lmaux:'
ccc!! Allocate permanent class arrays, maps and other initialization ---
c      call clsprm(1,sspec)!, v_slat) !,v_sbz) !,v_sstr)!,sstrn)v_sctrl,v_sham,v_spot,,v_sarry
cc!! Read available class parameters from file ---
c      call aiocls(.false.,0, sspec,1,0)
c      call aiocls(.false.,17,sspec,1,0)
c      if(lc>=20) call clsprp(1,sspec)
c      call defspc(sspec)

      nbas=ctrl_nbas
      nclass=ctrl_nclass
      nl=ctrl_nl
      nspec=ctrl_nspec
      nsp=ctrl_nspin

c      i_copy_size=size(ctrl_modep)
c      call i8icopy(i_copy_size,ctrl_modep,1,modep,1)
      modep = ctrl_modep
      lpbc = 0
      nclasp=ctrl_nclass !sarray%nclasp
      avw=lat_avw
      alat=lat_alat
c      i_copy_size=size(lat_plat)
c      call dcopy(i_copy_size,lat_plat,1,plat,1)
      plat=lat_plat
      nkd=lat_nkd
      nkq=lat_nkq
c      mxnbr=str_mxnbr !sstr%mxnbr
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
      nbaspp = 2*nbasp - nbas
      stdo = lgunit(1)
      j = 10
      if (cmdopt('--shorten',j-1,0,outs)) then
        call shorps ( nbasp , plat , modep , rv_a_opos , rv_a_opos )
      endif

c$$$C ... Read from restart file
c$$$      if (cmdopt('--rs=',5,0,outs)) then
c$$$        irs ( 1 ) = irs11 !iand(7,int(ctrl_lrs))
c$$$        if (irs(1) .gt. 0) then
c$$$          ifi = fopna('rsta',-1,0)
c$$$          call asars ( 0 , sctrl , ssite , sspec , slat , spot , sbz ,
c$$$     .      spot%rv_a_opnu , spot%rv_a_oqnu , .false. , ifi )
c$$$          call fclr('rsta',ifi)
c$$$        endif
c$$$      endif

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
C        xx = dglob('lrel',1d0,1)
C        xx = dglob('nsp',1d0,1)
C       Initial estimate for sphere radii: overlapping atom potentials
          allocate(zz_rv(nspec))
          allocate(rmt_rv(nspec))
c          i_copy_size=1;
          do i_spackv=1,nspec
c            call spackv_array_copy_r8_r8 ( 'u' , sspec ( i_spackv ) %z , 
c     .      i_copy_size , i_spackv + 1 - 1 , zz_rv )
            zz_rv (i_spackv) = sspec(i_spackv)%z
            rmt_rv(i_spackv) = sspec(i_spackv)%rmt
          enddo
c          i_copy_size=1;
c          do i_spackv=1,nspec
c            call spackv_array_copy_r8_r8 ( 'u' , sspec ( i_spackv ) %rmt 
c     .      , i_copy_size , i_spackv + 1 - 1 , rmt_rv )
c          enddo
          allocate(lock_iv(nspec))
          lock_iv(:)=0
          do  i = 1, nspec
            call ivset ( lock_iv , i , i , bitand ( int ( sspec ( i ) 
     .      %mxcst ) , 2 ) )
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
c          do i_spackv=1,nspec
c            sspec(i_spackv)%rmt = rmt_rv(i_spackv)
c          enddo
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

C --- Show neighbors by shell ---
        outs = ' '
        j = 8
c$$$        if (cmdopt('--shell',j-1,0,outs)) then
c$$$          allocate(eula(nbas,3))
c$$$          call shoshl ( outs ( j: ) , nbas , rv_a_opos , plat , str_mxnbr ,
c$$$     .         z , slabl , sarray%rv_a_oclabl , sarray%iv_a_oips , sarray%iv_a_oipc , spot%rv_a_oves , 
c$$$     .         eula  , nclass )
c$$$          deallocate(eula)
c$$$        endif
C --- Show angles between neighbors ---
        j = 9
        if (cmdopt('--angles',j-1,0,outs)) then
          call shoang ( outs ( j: ) , nbas , rv_a_opos , plat , str_mxnbr ,
     .      slabl , iv_a_oips )
        endif
c$$$C --- Check whether basis w(opos) differs from site file by translation ---
c$$$        j = 9
c$$$        if (cmdopt('--basis=',j-1,0,outs)) then
c$$$C       call upack('lat nsgrp osymgr',slat,ngrp,og,0,0,0)
c$$$          fnam = outs(j:)
c$$$          call info(20,1,0,' checking whether basis equivalent to file '
c$$$     .    //fnam//'...',0,0)
c$$$          j = iosits(8070,3d0,0,fnam,ifi,slabl,alat,plat,nbas,nspec,
c$$$     .    sspec,ssite)
c$$$          allocate(pos2_rv(3,nbas))
c$$$          allocate(ips2_iv(nbas))
c$$$c          do i_spackv=1,nbas
c$$$c          i_copy_size=size(ssite(1)%pos)
c$$$c            call spackv_array_copy_r8_r8 ( 'u' , ssite ( i_spackv ) %pos 
c$$$c     .      , i_copy_size , i_spackv + 1 - 1 , pos2_rv )
c$$$c          i_copy_size=1;
c$$$c            call spackv_array_copy_i8_i ( 'u' , ssite ( i_spackv ) %spec 
c$$$c     .      , i_copy_size , i_spackv + 1 - 1 , ips2_iv )
c$$$c          enddo
c$$$          do i_spackv=1,nbas
c$$$            pos2_rv(:,i_spackv) = ssite(i_spackv)%pos 
c$$$            ips2_iv(i_spackv)   = ssite(i_spackv)%spec 
c$$$          enddo
c$$$          allocate(rv_a_og(9*ngmx))
c$$$          call symlat( plat , ngrp , rv_a_og , j )
c$$$          j = cmplat( nbas , plat , plat , ngrp , rv_a_og , iv_a_oips
c$$$     .       , rv_a_opos , ips2_iv , pos2_rv )
c$$$          call fexit(j,1,' Exit %i lmchk --basis= ...',j)
c$$$        endif

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
     .     , rmax , rv_a_oclabl , iv_a_oipc , modep , z , iv_a_ontab , iv_a_oiax ,
     .     rv_a_opos , j )
C ... Write positions to file
c        if ( cmdopt ( '--wpos=' , 7 , 0 , outs ) ) call iopos ( t , 0
c     .       , outs ( 8: ) , nbasp , rv_a_opos )
c         if ( cmdopt ( '--wpos=' , 7 , 0 , outs ) ) then
c            allocate(bas(3,nbas))
c            call dcopy(nbas*3, rv_a_opos,1,bas,1)
c            do i=1,nbas
c               write(6,*)'ibas pos=',i,bas(:,i)
c            enddo
c         endif
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



      subroutine spec2c(nspec,nclass,ics,rmts,rmtc,z,zc,lock,lockc)
C- Copy species data to class data
C ----------------------------------------------------------------------
Ci Inputs
Ci   nspec
Ci   nclass:number of inequivalent classes
Ci   ics   :species table: class ic belongs to species ics(ic)
Ci   rmtc  :rmt by species
Ci   z     :z by species
Co Outputs
Co   rmtc  :rmt by class
Co   zc    :Z by class
co   lockc :lock by class
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   11 Aug 08
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nspec,nclass,ics(nclass),lock(nspec),lockc(nclass)
      double precision rmts(nspec),rmtc(nclass),z(nspec),zc(nclass)
C ... Local parameters
      integer j,k

      do  k = 1, nclass
        j = ics(k)
        rmtc(k) = rmts(j)
        zc(k) = z(j)
        lockc(k) = lock(j)
C       if (iprint() .gt. 60) write(*,310) k,rmtc(k)
      enddo

C 310 format(1x,'class ',I3,T15,'rmt = ',f10.7)

      end subroutine spec2c



      


C================================================================
c     take subs/spinav.F here
      subroutine spinav(mode,nclass,nl,nsp,pnu,qnu)
C     - Averages up+down spin moments + pp's for all classes
C     ----------------------------------------------------------------------
C     i Inputs
C     i   mode  :0 average spins
C     i         :1 do not average, but exchange spins
C     i   nclass:number of inequivalent classes
C     i   nl    :(global maximum l) + 1
C     i   nsp   :2 for spin-polarized case, otherwise 1
C     i   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
C     i          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
C     i   qnu   :energy-weighted moments of the sphere charges
C     o Outputs :moments are spin-averaged
C     i   pnu   :spin-averaged (mode=0) or spin-flipped (mode=1)
C     i   qnu   :spin-averaged (mode=0) or spin-flipped (mode=1)
C     o   nsp   :set to 1 on output (mode=0)
C     l Local variables
C     l         :
C     r Remarks
C     r   
C     u Updates
C     u   10 Jan 06 Added mode
C     ----------------------------------------------------------------------
      implicit none
      integer mode,nclass,nl,nsp,ic
      double precision pnu(nl,nsp,nclass),qnu(3,nl,nsp,nclass)

      if (nsp .eq. 1) return

      if (mode .gt. 1) then
         call rx('spinav: bad mode')
      elseif (mode .eq. 1) then
         do  ic = 1, nclass
            call dswap(nl,pnu(1,2,ic),1,pnu(1,1,ic),1)
            call dswap(3*nl,qnu(1,1,2,ic),1,qnu(1,1,1,ic),1)
         enddo
         return
      endif

      do  ic = 1, nclass
         call daxpy(nl,1d0,pnu(1,2,ic),1,pnu(1,1,ic),1)
         call dscal(nl,.5d0,pnu(1,1,ic),1)
         call daxpy(3*nl,1d0,qnu(1,1,2,ic),1,qnu(1,1,1,ic),1)
      enddo

      do  ic = 2, nclass
         call dcopy(nl,pnu(1,1,ic),1,pnu(1,ic,1),1)
         call dcopy(3*nl,qnu(1,1,1,ic),1,qnu(1,1,ic,1),1)
      enddo
      nsp = 1
      end
      
      subroutine aiocls(lio,mode,sspec,ic1,ic2) !sctrl,sham,,slat spot,sarray,
      use m_clsprm,only: iv_a_oipqx,
     &     rv_a_ogrrme ,  rv_a_opmpol , rv_a_opnu , rv_a_opp , rv_a_opprel,
     &     rv_a_oqnu , rv_a_oqpp , rv_a_oqt , rv_a_orhrmx ,rv_a_opnu,
     &     rv_a_osop , rv_a_ovdif , rv_a_oves , rv_a_ovintr , rv_a_ovrmax 
      use m_mksym,only: rv_a_oclabl,iv_a_oics,ctrl_nclass
      use m_struc_def           !Cgetarg
      use m_lmfinit,only:  ctrl_nl,ctrl_nspin,lrel_g=>lrel


      use m_ext,only:sname
      use m_lattic,only: rv_a_opos
C     - File I/O atomic data for classes ic1..ic2
C     ----------------------------------------------------------------------
C     i Inputs
C     i   lio    F for read, T for write
C     i   mode   ones digit
C     i          1 use data from first class corresponding to same species,
C     i            if data from own class is missing
C     i          2 like 1, but attempt to read data from disk anyway
C     i          3 make no attempt to read from disk but copy info from
C     i            first class corresponding to same species, if data missing
C     i          4 Add 4 if to use default P,Q when not otherwise supplied
C     i          10s digit
C     i          1 assemble background rho
C     i   sarray:structure containing offsets to various arrays
C     i     Elts read: nclasp oclabl ohave oics onrc
C     i     Stored:
C     i     Passed to:
C     i   sctrl :struct for program flow parameters; see routine uctrl
C     i     Elts read: nl nspin zbak nbas nclass
C     i     Stored:    zbak
C     i     Passed to: lgors
C     i   sham  :struct for parameters defining hamiltonian; see routine uham
C     i          (not used now)
C     i   spot  :struct for information about the potential; see routine upot
C     i     Elts read: opp osop ogrrme opnu oqnu orhrmx ovrmax ovintr opmpol
C     i                oves
C     i     Stored:
C     i     Passed to:
C     i   sspec :struct for species-specific information; see routine uspec
C     i     Elts read: lmxa idmod p q
C     i     Stored:
C     i     Passed to:
C     i   slat  :struct for lattice information; see routine ulat
C     i     Elts read: vol
C     i     Stored:
C     i     Passed to:
C     i   ic1,ic2: range of classes to read data
C     r Remarks
C     r   Right now, read always takes data from file if available
C     u Updates
C     u   09 Nov 07 Corrected sign of default moment (paioc2)
C     u   29 Sep 04 Reads/writes relativistic ppar's
C     u   26 Apr 03 Added MPI calls
C     u   07 Feb 03 adjusted for redimensioned sop
C     u   30 May 02 Assign better default P
C     u   28 Apr 98 I/O of radial matrix elements of grad
C     u   28 Sep 00 Added setting default P,Q
C     ----------------------------------------------------------------------
      implicit none
C     Passed parameters
      logical lio
      integer mode,ic1,ic2
c      type(s_array)::sarray
c      type(s_ctrl)::sctrl
      type(s_spec)::sspec(*)
c      type(s_lat)::slat
c      type(s_ham)::sham
c      type(s_pot)::spot
C     Local variables
      character*8 clabl,alabel, outs1*20, outs2*20, outs3*20, outs*80
      logical sw,aiomom,aiopar,aiopot,lpot,lcor,scat,lgors,
     .     aiogen,aiocor,lrell,lgen
      logical havepq,havepp,haveso,haveop,haveva,havemp,
     .     readpq,readpp,readso,readop,readva,readmp
      integer::  onrc , n0 , nglob !olpp , olppr , olsop , olgrme , olva , olmp , opot  , ocor ,
      parameter(n0=10)
      integer ic,is,lmx,k,nl,nsp,ifi,jfi,lmxx,nspx,nrx,nr,
     .     idmod(n0),isw,bitand,i2,nclasp,iclbsj,icmap,jc,nbas,nclass,
     .     nclspp,iprint,lgunit,mode0,mode00,lrel,i_copy_size,ifile_handle
c     integer mpipid,procid
      double precision rhrmx,vrmax(2),ves,z,rmxx,ax,qc,dq,vrmxx(2),
     .     sumec,sumtc,sumev,thrpv,ekin,utot,rhoeps,etot,a,rmax,dval,
     .     zbak(2),dgets,pdf(n0,2),qdf(n0,2),pnuloc(100),qnuloc(100)
      real(8):: olpp(100),olppr(3000),olsop(500),olgrme(500),olva(100),olmp(3000),opot(3000),ocor(3000)
      
      nl=ctrl_nl
      nsp=ctrl_nspin
      nclass=ctrl_nclass
      nclasp=ctrl_nclass !sarray%nclasp
      lrel = lrel_g !globalvariables%lrel
      mode0 = mod(mode,10)
      mode00 = mod(mode0,4)
      i2 = ic2
      if (i2 .eq. 0) i2 = nclasp
      do  10  ic = ic1, i2
         is=iv_a_oics(ic)
         icmap = iclbsj ( is , iv_a_oics , - nclasp , 1 ) 
         if (icmap .eq. ic .or. mode00 .eq. 0) icmap = 0
         lmx=sspec(is)%lmxa
         i_copy_size=size(sspec(is)%idmod) 
         call icopy(i_copy_size,sspec(is)%idmod,1,idmod,1) 
         i_copy_size=size(sspec(is)%p) 
         call dcopy(i_copy_size,sspec(is)%p,1,pdf,1) 
         i_copy_size=size(sspec(is)%q) 
         call dcopy(i_copy_size,sspec(is)%q,1,qdf,1) 
         call dpzero(pnuloc,100)
         call dpzero(qnuloc,100)

C     --- Open the atom file ---
         outs1 = ' '
         outs2 = ' '
         outs3 = ' '
c         call r8tos8 ( dval ( sarray%rv_a_oclabl , ic ) , clabl ) 
         clabl = rv_a_oclabl(ic)
         if (mode00 .ne. 3) then
            ifi=ifile_handle()
            open(ifi,file=trim(clabl)//'.'//trim(sname))
         endif

C     --- Copy what is passed through to holding arrays ---
         k = nl*nsp
         call pvaioc ( iv_a_oipqx , 1 , ic , icmap , havepq , jc ) 
         readpq = .not. havepq .or. mode00 .eq. 2 .and. jc .ne. ic
         readpq = readpq .and. mode00 .lt. 3
         if (havepq) then
            call dpscop ( rv_a_opnu , pnuloc , k , 1 + ( jc - 1 ) * k , 1  , 1d0 ) 
            call dpscop ( rv_a_oqnu , qnuloc , 3 * k , 1 + ( jc - 1 ) * 3 
     .           * k , 1 , 1d0 ) 
            call dpscop ( rv_a_orhrmx , rhrmx , 1 , jc , 1 , 1d0 ) 
            call dpscop ( rv_a_ovrmax , vrmax , 2 , 2 * jc - 1 , 1 , 1d0 ) 
            if (ic.ne.jc) call awrit0('%a pq,',outs1,len(outs1),0)
         endif
         k = 6*nl*nsp
         call pvaioc ( iv_a_oipqx , 2 , ic , icmap , havepp , jc ) 
         readpp = .not. havepp .or. mode00 .eq. 2 .and. jc .ne. ic
         readpp = readpp .and. mode00 .lt. 3
         if (havepp) then
            call dpscop ( rv_a_opp , olpp , k , 1 + ( jc - 1 ) * k  , 1 , 1d0 ) 
            if (lrel .eq. 2) then
               k = 4*nl*2*nl*2*2
               call dpscop ( rv_a_opprel , olppr , k , 1 + ( jc - 1 ) 
     .              * k , 1 , 1d0 ) 
            endif
            call dpscop ( rv_a_oves , ves , 1 , jc , 1 , 1d0 ) 

            if (ic.ne.jc) call awrit0('%a pp,',outs1,len(outs1),0)
         endif
         k = (nl*nsp)**2
         call pvaioc ( iv_a_oipqx , 8 , ic , icmap , haveva , jc ) 
         haveva = haveva .and. allocated(rv_a_ovintr)
         readva = .not. haveva .or. mode00 .eq. 2 .and. jc .ne. ic
         readva = readva .and. allocated(rv_a_ovintr) .and. mode00 .lt. 3
         if (haveva) then
            call dpscop ( rv_a_ovintr , olva , k , 1 + ( jc - 1 ) * k , 1 , 1d0 ) 
            if (ic.ne.jc) call awrit0('%a va,',outs1,len(outs1),0)
         endif
         k = nl**2*(2*nl-1)*3*nsp
         call pvaioc ( iv_a_oipqx , 16 , ic , icmap , havemp , jc ) 
         havemp = havemp .and. allocated(rv_a_opmpol)
         readmp = .not. havemp .or. mode00 .eq. 2 .and. jc .ne. ic
         readmp = readmp .and. allocated(rv_a_opmpol) .and. mode00 .lt. 3
         if (havemp) then
            call dpscop ( rv_a_opmpol , olmp , k , 1 + ( jc - 1 ) * k 
     .           , 1 , 1d0 ) 
            if (ic.ne.jc) call awrit0('%a mp,',outs1,len(outs1),0)
         endif
         k = nl*nsp*nsp*9
         call pvaioc ( iv_a_oipqx , 4 , ic , icmap , haveso , jc ) 
         haveso = haveso .and. allocated(rv_a_osop)
         readso = .not. haveso .or. mode00 .eq. 2 .and. jc .ne. ic
         readso = readso .and. allocated(rv_a_osop) .and. mode00 .lt. 3
         if (haveso) then
            call dpscop ( rv_a_osop , olsop , k , 1 + ( jc - 1 ) * 
     .           k , 1 , 1d0 ) 

            if (ic.ne.jc) call awrit0('%a sop,',outs1,len(outs1),0)
         endif
         k = 16*nl*nsp
         call pvaioc ( iv_a_oipqx , 32 , ic , icmap , haveop , jc ) 
         haveop = haveop .and. allocated(rv_a_ogrrme)
         readop = .not. haveop .or. mode00 .eq. 2 .and. jc .ne. ic
         readop = readop .and. allocated(rv_a_ogrrme) .and. mode00 .lt. 3
         if (haveop) then
            call dpscop ( rv_a_ogrrme , olgrme , k , 1 + ( jc - 1 )* k , 1 , 1d0 ) 
            if (ic.ne.jc) call awrit0('%a opp,',outs1,len(outs1),0)
         endif

C     --- File WRITE ---
         if (lio) then
            lgen = .false.
            lpot = .false.
            lcor = .false.
C     ... Pick up GEN and POT, if available, to save again
            if (scat(iabs(ifi),'GEN:',':',.true.)) then
               lgen = aiogen(alabel,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,dq,
     .              vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
            endif
            if (scat(iabs(ifi),'POT:',':',.true.)) then
               read(ifi,102) nr,nsp,a,rmax
 102           format(2i5,2f12.5)
               lpot = aiopot(nr,nsp,a,rmax,-99d0,opot,ifi)
            endif
            lcor = aiocor(nr,nsp,a,rmxx,ocor,sumec,sumtc,ifi)

            rewind ifi
            jfi = -ifi
            if (lgen) sw = aiogen(clabl,z,rmxx,lmxx,nspx,lrell,nrx,ax,qc,
     .           dq,vrmxx,sumec,sumev,thrpv,ekin,utot,rhoeps,etot,jfi)
            if (havepq) sw = aiomom(clabl,pnuloc,qnuloc,idmod,
     .           nl,lmx,nsp,rhrmx,vrmax,jfi)
            if (havepp) sw = aiopar(clabl,lrel,olpp,olppr,ves,nl,
     .           lmx,nsp,jfi)
            if (haveva) sw = aiova(clabl,olva,nl,lmx,nsp,jfi)
            if (havemp) sw = aiomp(clabl,olmp,nl,2*nl-2,nsp,jfi)
            if (haveso) sw = aiosop(clabl,olsop,nl,lmx,nsp,jfi)
            if (haveop) sw = aiorme(clabl,olgrme,nl,nsp,jfi)
            if (lpot)   sw = aiopot(nr,nsp,a,rmax,-99d0,opot,jfi)
            if (lcor) lcor = aiocor(nr,nsp,a,rmxx,ocor,sumec,sumtc,jfi)

C     --- File READ ---
         else

C     ... Copy whatever is available on disk to holding arrays
            if (readpq .or. mode0 .ge. 4) then
               if (readpq) rewind ifi
               if (readpq) readpq = aiomom(clabl,pnuloc,qnuloc,idmod,
     .              nl,lmx,nsp,rhrmx,vrmax,ifi)
C     Couldn't read from atom file ; take default values
               if (readpq) call awrit0('%a pq,',outs2,len(outs2),0)
               if (mode0 .ge. 4 .and. .not. (readpq .or. havepq)) then
C     call dmcpy(pdf,n0,1,pnuloc,nl,1,nl,nsp)
                  call paioc2(nsp,nl,n0,pdf,qdf,pnuloc,qnuloc)
                  call awrit0('%a pq,',outs3,len(outs2),0)
                  call dvset(vrmax,1,2,-.7d0)
                  rhrmx = .1d0
                  readpq = .true.
               endif
            endif
            if (readpp) then
               readpp = aiopar(clabl,lrel,olpp,olppr,ves,nl,lmx,nsp,
     .              ifi)
               if (readpp) call awrit0('%a pp,',outs2,len(outs2),0)
            endif
            if (readso .and. allocated(rv_a_osop) ) then
               readso = aiosop(clabl,olsop,nl,lmx,nsp,ifi)
               if (readso) call awrit0('%a so,',outs2,len(outs2),0)
            endif
            if (readop .and. allocated(rv_a_ogrrme)) then
               readop = aiorme(clabl,olgrme,nl,nsp,ifi)
               if (readop) call awrit0('%a op,',outs2,len(outs2),0)
            endif
            if (readva .and. allocated(rv_a_ovintr) ) then
               readva = aiova(clabl,olva,nl,lmx,nsp,ifi)
               if (readva) call awrit0('%a va,',outs2,len(outs2),0)
            endif
            if (readmp .and. allocated(rv_a_opmpol) ) then
               readmp = aiomp(clabl,olmp,nl,2*nl-2,nsp,ifi)
               if (readmp) call awrit0('%a mp,',outs2,len(outs2),0)
            endif

C     ... Update what parameters are available
            iv_a_oipqx(ic) = isw(haveop.or.readop)*32+
     .           isw(havemp.or.readmp)*16+
     .           isw(haveva.or.readva)*8 +
     .           isw(haveso.or.readso)*4 +
     .           isw(havepp.or.readpp)*2 +
     .           isw(havepq.or.readpq)*1 +
     .           iv_a_oipqx(ic) - bitand(iv_a_oipqx(ic),63)
            k = nl*nsp
            if (havepq .or. readpq) then
               call dpscop ( pnuloc , rv_a_opnu , k , 1 , 1 + ( ic - 1 ) * k, 1d0 ) 
               call dpscop ( qnuloc , rv_a_oqnu , 3 * k , 1 , 1 + ( ic - 1 ) 
     .              * 3 * k , 1d0 ) 
               call dpscop ( rhrmx , rv_a_orhrmx , 1 , 1 , ic , 1d0 ) 
               call dpscop ( vrmax , rv_a_ovrmax , 2 , 1 , 2 * ic - 1 , 1d0 ) 
            endif
            k = 6*nl*nsp
            if (havepp .or. readpp) then
               call dpscop ( olpp , rv_a_opp , k , 1 , 1 + ( ic - 1 ) * k , 1d0 ) 
               if (lrel .eq. 2) then
                  k = 4*nl*2*nl*2*2
                  call dpscop ( olppr , rv_a_opprel , k , 1 , 1 + ( ic - 
     .                 1 ) * k , 1d0 ) 
               endif
               call dpscop ( ves , rv_a_oves , 1 , 1 , ic , 1d0 ) 
            endif
            k = (nl*nsp)**2
            if ( ( haveva .or. readva ) .and. allocated(rv_a_ovintr ) ) call 
     .           dpscop ( olva  , rv_a_ovintr , k , 1 , 1 + ( ic - 1 ) * k , 1d0 ) 

            k = nl**2*(2*nl-1)*3*nsp
            if ( ( havemp .or. readmp ) .and. allocated(rv_a_opmpol ) ) call 
     .           dpscop ( olmp  , rv_a_opmpol , k , 1 , 1 + ( ic - 1 ) * k , 1d0 ) 
            k = nl*nsp*nsp*9
            if ( ( haveso .or. readso ) .and. allocated(rv_a_osop ) ) 
     .           call dpscop ( olsop  , rv_a_osop , k , 1 , 1 + ( ic - 1) * k , 1d0 ) 

            k = 16*nl*nsp
            if ( ( haveop .or. readop ) .and. allocated(rv_a_ogrrme ) ) 
     .           call dpscop ( olgrme  , rv_a_ogrrme , k , 1 , 1 + ( ic - 1 ) * k , 1d0 ) 
         endif

         if (mode00 .lt. 3) call fclr(clabl,ifi)
         if (iprint() .gt. 40) then
            outs = ' '
            if (outs1 .ne. ' ') then
               call awrit1('%x '//clabl//'%a: copied '//outs1//'%a%b from '//
     .              'class %i',outs,len(outs),0,jc)
            endif
            if (outs2 .ne. ' ') then
               call awrit0('%x '//clabl//'%a: read '//outs2//'%a%b from '//
     .              'disk',outs,len(outs),0)
            endif
            if (outs3 .ne. ' ') then
               if (outs2 .eq.  ' ') then
                  call awrit0('%x '//clabl//'%a: use defaults for: '//outs3//
     .                 '%a%b',outs,len(outs),0)
               else
                  call awrit0('%a; use defaults for: '//outs3//'%a%b',
     .                 outs,len(outs),0)
               endif
            endif
            if (outs .eq. ' ')
     .           call awrit0(' '//clabl//'%a: nothing read',outs,len(outs),0)
            call awrit0(' aiocls class'//outs,' ',-len(outs),lgunit(1))
         endif
 10   continue
      end subroutine aiocls 

      subroutine pvaioc(initc,mask,ic0,icmap,lhave,ic)
C     - Find whether data avaiable either in class or mapped class
      implicit none
      logical lhave
      integer initc(1),mask,ic0,icmap
      integer ic
      ic = ic0
      lhave = mod(initc(ic)/mask,2) .eq. 1
      if (.not. lhave .and. icmap .ne. 0) then
         lhave = mod(initc(icmap)/mask,2) .eq. 1
         ic = icmap
      endif
      end subroutine pvaioc 

      subroutine paioc2(nsp,nl,n0,pat,qat,pnu,qnu)
C     - Widget to copy pat,qat to pnu,qnu
      implicit none
      integer n0,nl,nsp
      double precision pat(n0,2),qat(n0,2),pnu(nl,nsp),qnu(3,nl,nsp)
      integer i,il
      do  i = 1, nsp
         do  il = 1, nl
C     pnu(il,i) = int(pat(il,i)) + .5d0
            pnu(il,i) = pat(il,i)
            qnu(1,il,i) = qat(il,1)/nsp
            if (nsp .eq. 2) then
               if (pat(il,i) .eq. 0) pnu(il,i) = pnu(il,1)
               qnu(1,il,i) = qat(il,1)/nsp + qat(il,2)/2*dble(3-2*i)
            endif
            qnu(2,il,i) = 0d0
            qnu(3,il,i) = 0d0
         enddo
      enddo
      end subroutine paioc2 


      subroutine shoang(sopts,nbas,bas,plat,mxnbr0,slabl,ips)

C- Print bond angles
C ----------------------------------------------------------------
Ci Inputs
Ci   sopts :a set of modifiers, with the syntax
Ci         :  [:r=#][:spec=spec-list][:spec=spec-list]
Ci         :  :r=# sets range for shells
Ci         :  :sites=site-list collects angles only for sites within list
Ci         :  :spec=list       prints angles for bonds connecting to
Ci                             species in list
Ci   bas   :basis vectors, in units of alat
Ci   plat  :primitive lattice vectors, in units of alat
Ci   mxnbr0
Ci   slabl :list of species labels
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   nbas  :size of basis
Cu Updates
Cu   13 Sep 01 Added options sopts.  Altered argument list.
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,mxnbr0
      double precision bas(3,nbas),plat(3,3)
      character(*):: slabl(*)
      integer ips(nbas)
      character sopts*(*)
C ... Local parameters
      double precision avwsr,avw,range,xx,zz
      integer:: npr(2) , mxnbr , ib , nshell , nmxshl , niax , j , 
     .j1 , j2 , nspec , mxint , nsites ,  nbonds ,  parg  !oslist ,ilst ,
     ., m , iv(10) , ilast
      integer ,allocatable :: iax_iv(:)
      real(8) ,allocatable :: wk_rv(:)
      integer ,allocatable :: num_iv(:)
      real(8) ,allocatable :: ang_rv(:)
      real(8) ,allocatable :: d_rv(:)
      integer ,allocatable :: blist_iv(:),slist(:)

      parameter (niax=10)
      character dc*1
c      external nghbor,nwordg,pvang1,pvang2,rxs

C ... Setup
      range = 2.5d0
      nspec = mxint(nbas,ips)
      nsites = 0
      nbonds = 0
c      call defi(oslist,nbas)
      allocate(blist_iv(nbas),slist(nbas))


C ... Switches
      dc = sopts(1:1)
      if (dc .ne. ' ') then
        j2 = 0
C   ... Return here to resume parsing for arguments
   10   continue
        j2 = j2+1
        if (sopts(j2:j2) .eq. dc) goto 10
        j1 = min(len(sopts),j2)
        call nwordg(sopts,0,dc//' ',1,j1,j2)
        if (j2 .ge. j1) then
          if (.false.) then

C         range
          elseif (sopts(j1:j1+1) .eq. 'r=')  then
            m = 0
            j = parg('r=',4,sopts(j1:),m,len(sopts(j1:)),
     .      dc//' ',1,1,iv,range)
            if (j .le. 0) goto 999

C         Site list
          elseif (sopts(j1:j1+4) .eq. 'sites') then
            if (sopts(j1+5:j1+5) .eq. '=') sopts(j1+5:j1+5) = dc
            call baslst(0,10,sopts(j1+5:),ilast,ips,nbas,slabl,zz,0,' ',
     .      xx,nsites,slist)

C         Bond list
          elseif (sopts(j1:j1+4) .eq. 'bonds') then
            if (sopts(j1+5:j1+5) .eq. '=') sopts(j1+5:j1+5) = dc
            call baslst ( 0 , 10 , sopts ( j1 + 5: ) , ilast , ips , nbas 
     .      , slabl , zz , 0 , ' ' , xx , nbonds , blist_iv )

          endif
          goto 10

        endif
      endif

      if (mxnbr0 .eq. 0) then
        mxnbr = 2*range**3
      else
        mxnbr = mxnbr0
      endif
      allocate(iax_iv(niax*mxnbr))

      allocate(wk_rv(mxnbr))

      avw = avwsr(plat,1d0,xx,nbas)

C --- For every site in list, generate tables of bond angles ---
c      ilst = 0
      do  ib = 1, nbas
        if (nsites .gt. 0) then
c          if (w(oslist+ilst) .ne. ib) goto 31
          if (slist(ib) .ne. ib) goto 31
        endif
c        ilst = ilst+1

C   ... Get neighbor lists

        call nghbor ( nbas , plat , bas , range * avw , range * avw , 
     .  ib , mxnbr , npr , iax_iv , wk_rv )


C   ... Get shell dimensions
        allocate(num_iv(abs(-nspec*npr(1))))
        if (-nspec*npr(1)<0) num_iv(:)=0

        call pvang1 ( npr ( 1 ) , nbas , plat , bas , iax_iv , ips 
     .  , num_iv , nshell , nmxshl )

        if (allocated(num_iv)) deallocate(num_iv)

        allocate(num_iv(abs(-nspec**2*nshell**2)))
        if (-nspec**2*nshell**2<0) num_iv(:)=0

        allocate(ang_rv(nspec**2*nshell**2*nmxshl**2))

        allocate(d_rv(nshell))


C   ... Print bond angles
        call pvang2 ( npr ( 1 ) , nbas , nspec , nshell , nmxshl , plat 
     .  , bas , iax_iv , ips , slabl , nbonds , blist_iv , num_iv 
     .  , ang_rv , d_rv )

        if (allocated(d_rv)) deallocate(d_rv)
        if (allocated(ang_rv)) deallocate(ang_rv)
        if (allocated(num_iv)) deallocate(num_iv)

   31   continue
      enddo
      if (allocated(wk_rv)) deallocate(wk_rv)
      if (allocated(iax_iv)) deallocate(iax_iv)
      deallocate(blist_iv)
c      call rlse(oslist)
      deallocate(slist)
      return

  999 call rxs('shoang: failed to parse ',sopts)
      end subroutine shoang


      subroutine pvang1(npr,nbas,plat,bas,iax,ips,num,nshell,nmxshl)

C- Help routine for shoang
C ----------------------------------------------------------------------
Ci Inputs
Ci   npr   :number of pairs in neighbor table
Ci   nbas  :size of basis
Ci   plat  :primitive lattice vectors, in units of alat
Ci   bas   :basis vectors, in units of alat
Ci   iax   :neighbor table containing pair information for one site.
Ci         :Table must be sorted by increasing distance from iax(1)
Ci   ips   :species table: site ib belongs to species ips(ib)
Co Outputs
Co   num   :(ishell,is) number of pairs in shell ishell of species is
Co   nshell:number of shells
Co   nmxshl:max value of num
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer npr,nbas,nshell,nmxshl,niax,ips(1),num(npr,1)
      parameter (niax=10)
      integer iax(niax,*)
      double precision plat(3,3),bas(3,1)
C ... Local parameters
      integer i,is
      double precision d,wk(0:3),tol

      tol = 1d-6
C --- Get no. of shells and max no. of atoms in 1 shell and 1 class ---
      nshell = 0
      nmxshl = 0
      d = 0d0
      do  i = 2, npr
        is = ips(iax(2,i))
        call dlmn(nbas,plat,bas,iax(1,i),wk)
C       distance changed by more than tol ... new shell
        if (dabs(wk(0)-d) .gt. tol) then
          nshell = nshell + 1
          d = wk(0)
        endif
        num(nshell,is) = num(nshell,is) + 1
        nmxshl = max0(nmxshl,num(nshell,is))
      enddo
      end subroutine pvang1


      subroutine pvang2(npr,nbas,nspec,nshell,nmxshl,plat,bas,iax,
     .ips,slabl,nbonds,blist,num,ang,d)

C- Kernel called by shoang
C ----------------------------------------------------------------------
Ci Inputs
Ci   npr   :number of neighbors connecting site ib=iax(1,1)
Ci   nbas  :size of basis
Ci   nspec :number of species
Ci   nshell:number of shells
Ci   nmxshl:dimensions ang
Ci   plat  :primitive lattice vectors, in units of alat
Ci   bas   :basis vectors, in units of alat
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   slabl :struct containing global strings
Ci   num
Co Outputs
Co   ang   :table of angles
Co   d     :table of distances for each shell
Co   Angles and distances are printed out
Cu Updates
Cu   13 Sep 01
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer npr,nbas,nspec,nshell,nmxshl,niax,ips(1),
     .num(nspec,nspec,nshell,nshell),nbonds,blist(nbonds)
      parameter (niax=10)
      integer iax(niax,*)
      character*8 slabl(*)
C ... Local parameters
      double precision plat(3,3),bas(3,1),d(nshell),
     .ang(nmxshl**2,nspec,nspec,nshell,nshell)
      integer i,j,n,is,js,nsh1,nsh2,nmx2,k
      double precision rdtodg,d1,d2,dp,ddot,wk1(0:3),wk2(0:3)
C ... External calls
      external dlmn,rxx

      nmx2 = nmxshl**2
      rdtodg = 45d0 / datan(1.d0)

C --- Accumulate bond angles by shell and class ---
      nsh1 = 0
      d1 = 0d0
      do  20  i = 2, npr
        is = ips(iax(2,i))
        call dlmn(nbas,plat,bas,iax(1,i),wk1)
        if (dabs(wk1(0)-d1) .gt. 1d-6) then
          nsh1 = nsh1 + 1
          d1 = wk1(0)
          d(nsh1) = d1
        endif
        nsh2 = nsh1
        d2 = d1
        if (nbonds .gt. 0) then
          k = 0
          call hunti(blist,nbonds,iax(2,i),0,k)
          if (k .ge. nbonds) goto 20
          if (blist(k+1) .ne. iax(2,i)) goto 20
        endif
        do  10  j = i+1, npr
          js = ips(iax(2,j))
          call dlmn(nbas,plat,bas,iax(1,j),wk2)
          if (dabs(wk2(0)-d2) .gt. 1d-6) then
            nsh2 = nsh2 + 1
            d2 = wk2(0)
          endif
          if (nbonds .gt. 0) then
            k = 0
            call hunti(blist,nbonds,iax(2,j),0,k)
            if (k .ge. nbonds) goto 10
            if (blist(k+1) .ne. iax(2,j)) goto 10
          endif
          dp = ddot(3,wk1(1),1,wk2(1),1)
          if (dp .gt.  1d0) dp =  1d0
          if (dp .lt. -1d0) dp = -1d0
          if (nsh1 .eq. nsh2 .and. js .lt. is) then
            num(js,is,nsh1,nsh2) = num(js,is,nsh1,nsh2) + 1
            n = num(js,is,nsh1,nsh2)
            ang(n,js,is,nsh1,nsh2) = rdtodg*dacos(dp)
          else
            num(is,js,nsh1,nsh2) = num(is,js,nsh1,nsh2) + 1
            n = num(is,js,nsh1,nsh2)
            ang(n,is,js,nsh1,nsh2) = rdtodg*dacos(dp)
          endif
          call rxx(n .gt. nmx2,'PVANG2: num gt nmx2')
   10   continue
   20 continue
      call rxx(nsh1 .ne. nshell,'PVANG2: nsh1 ne nshell')

C --- Printout ---
      print 400, iax(1,1), slabl(ips(iax(1,1)))
  400 format(/' Bond angles for site',i4,', species ',a/
     .' shl1    d1    shl2    d2     cl1      cl2       angle(s) ...')

      do  60  nsh1 = 1, nshell
        do  50  nsh2 = nsh1, nshell
          do  40  is = 1, nspec
            do  30  js = 1, nspec
              n = num(is,js,nsh1,nsh2)
              if (n .ne. 0) then
                print 401, nsh1, d(nsh1), nsh2, d(nsh2), slabl(is),
     .          slabl(js), (ang(i,is,js,nsh1,nsh2), i = 1, n)
              endif
   30       continue
   40     continue
   50   continue
   60 continue

  401 format(2(1x,i3,1x,f9.6),1x,2(1x,a8),20(4(1x,f7.2):/47x))

      end subroutine pvang2

      subroutine psho1(ltab,lnn,iax,nbas,nttab,nlst1,lst1,nlst2,lst2,
     .rtab,plat,pos,pos2,ipc,dclabl,ifi)

C- Kernel called by supcel to displace pairs radially
C ----------------------------------------------------------------------
Ci Inputs
Ci   ltab  :style which to print out neighbor table
Ci         :1s digit
Ci         :0 do nothing ;return.  Else print table as:
Ci         :1 (standard mode)
Ci         :    ib jb dpos(1,jb) dpos(2,jb) dpos(3,jb)
Ci         :2  (just the positions)
Ci         :   dpos(1,jb) dpos(2,jb) dpos(3,jb)
Ci         :3 (sparse matrix format)
Ci         :     1 jb dpos(1,jb)
Ci         :     2 jb dpos(2,jb)
Ci         :     3 jb dpos(3,jb)
Ci         :10s digit
Ci         :1 print out neighbor table for pos and
Ci         :  displaced pos2 as well
Ci         :
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   nbas  :size of basis
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   nlst1 :number of sites of "center" type
Ci   lst1  :list of sites of "center" type
Ci   nlst2 :number of sites of "neighbor" type
Ci   lst2  :list of sites of "neighbor" type
Ci   rtab  :site positions corresponding to entries in a neighbor table
Ci          relative to some center
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   pos2  :displaced basis vectors (ltab >= 10)
Ci   mark  :work array of dimension nbas
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   dclabl:class name, packed as a real number
Ci   ifi   :file handle
Co Inputs/Outputs
Cu Updates
Cu   08 Aug 07 case ltab=12: allow for numerical imprecision in vector
Cu   19 Apr 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ltab,nlst1,lst1(nlst1),nlst2,lst2(nlst2),ifi,nbas,lnn
      integer niax,nttab,ipc(*),mark(nbas)
      double precision pos(3,*),pos2(3,*),rtab(3,*),plat(9)
      character(8):: dclabl(*)
      parameter (niax=10)
      integer iax(niax,*)
C ... Local parameters
      integer iat,low1,low2,i1,i2,ipr,stdo,lgunit,i,oldi1,oldiat,
     .ic,jc
      character clabl*8, clabl2*8
      double precision d,ddot,dsqrt,dpos(3),qlat(9),tol
C     logical latvec
      logical::isanrg, l_dummy_isanrg

      tol = 1d-6

      if (ltab .eq. 0) return
      call mkqlat(plat,qlat,d)

      call info5(20,1,1,' ... shoshl: making neighbor list for'//
     .' %?#n#%-1j%i%j#%i# atom(s), style %i',nlst1,nbas,ltab,0,0)

      if (ltab .eq. 3)
     .call awrit1('%% rows 3 cols %i sparse',' ',80,ifi,nbas)

      call getpr(ipr)
      ipr = 100
      stdo = lgunit(1)
Ckino isanrg is logical function,       call isanrg(mod(ltab,10),1,3,' shoshl:','tab',.true.)
      l_dummy_isanrg=isanrg(mod(ltab,10),1,3,' shoshl:','tab',.true.)

      if (ltab .eq. 3)
     .call awrit1('%% rows 3 cols %i sparse',' ',80,ifi,nbas)
      if (ltab .eq. 13)
     .call awrit1('%% rows 6 cols %i sparse',' ',80,ifi,nbas)

      low1 = 0
      low2 = 0
      oldiat = 0
      oldi1 = 0
      do  iat = 1, nttab

        i1 = iax(1,iat)
        i2 = iax(2,iat)
C   ... If site i1 isn't in the supplied list, skip this pair
        if (nlst1 .ne. 0) then
          call hunti(lst1,nlst1,i1,0,low1)
          if (i1 .ne. lst1(low1+1)) goto 10
        endif

C   ... New shell
        if (i1 .ne. oldi1) then
          call iinit(mark,nbas)
          oldi1 = i1
          oldiat = iat
          ic = ipc(i1)
c          call r8tos8(dclabl(ic),clabl)
          clabl=dclabl(ic)
          call awrit1('# neighbor list for site %i, class '//clabl,
     .    ' ',80,ifi,i1)
C         Don't print out on-site entry
C         goto 10
        endif

C   ... If site i2 isn't in the supplied lst2, skip this pair
        if (nlst2 .ne. 0) then
          call hunti(lst2,nlst2,i2,0,low2)
          if (i2 .ne. lst2(low2+1)) goto 10
        endif

C   ... If i2 already marked, skip this pair
        if (mark(i2) .ne. 0) goto 10

        if (mod(ltab,10) .eq. 1) then
          jc = ipc(i2)
c          call r8tos8(dclabl(jc),clabl2)
          clabl2=dclabl(jc)
          d = dsqrt(ddot(3,rtab(1,iat),1,rtab(1,iat),1))
          if (ltab .lt. 10) then
            write(ifi,331) i1,i2,(rtab(i,iat),i=1,3),d,clabl,clabl2
          else
            do  i = 1, 3
              dpos(i) = pos2(i,i2)-pos(i,i2)
            enddo
            call shorbz(dpos,dpos,plat,qlat)
C           Only print out if displacement nonzero
            if (ddot(3,dpos,1,dpos,1) .ne. 0) then
              write(ifi,1331) i1,i2,(rtab(i,iat),i=1,3),
     .        (dpos(i),i=1,3),d,clabl,clabl2
              mark(i2) = lnn
            endif
          endif
  331     format(1x,2i4,3f12.7,2x,f12.7,2x,a,1x,a)
 1331     format(1x,2i4,3f12.7,2x,3f12.7,2x,f12.7,2x,a,1x,a)
        elseif (ltab .eq. 2) then
          write(ifi,332) (rtab(i,iat),i=1,3)
          mark(i2) = lnn
  332     format(1x,3f12.7)
        elseif (ltab .eq. 12) then
          do  i = 1, 3
            dpos(i) = pos2(i,i2)-pos(i,i2)
          enddo
C         Only print out if displacement nonzero
          call shorbz(dpos,dpos,plat,qlat)
          if (ddot(3,dpos,1,dpos,1) .gt. tol*tol) then
C         if (.not. latvec(1,tol,qlat,dpos)) then
            write(ifi,1332)
     .      (rtab(i,iat),i=1,3),(dpos(i),i=1,3)
            mark(i2) = lnn
          endif
 1332     format(1x,3f12.7,2x,3f12.7)
        elseif (ltab .eq. 3) then
          do  i = 1, 3
            write(ifi,333) i, i2, rtab(i,iat)
            mark(i2) = lnn
          enddo
  333     format(1x,i3,i5,3f12.7)
        endif

   10   continue
      enddo

      end subroutine psho1


      logical function aiomp(alabl,pmpol,nl,lmxv,nsp,ifi)
C- File I/O for ASA multipole moments integrals
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   alabl,nl,lmxv,nsp
Ci   ifi:    logical unit: positive for read, negative for write
Cio  pmpol: intra-atomic density-density response matrix
Cr Remarks
Cr Bugs
Cr   Input not checked when file nsp mismatches passed nsp
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*8 alabl
      integer nl,lmxv,nsp,ifi
      double precision pmpol(nl,nl,0:lmxv,3,nsp)
C Local parameters
      double precision x2(10),damax
      integer i,n,ip,k,a2vec,ix2(10),ipr,nl2,nsp2,l1,l2,lf,idamax
      logical scat,rdstrn,pars1v
      character*100 s

      aiomp = .false.
      call getpr(ipr)
C --- File READ ---
      if (ifi .gt. 0) then
C   ... return unless file has PMPOL: category
        if (.not. scat(ifi,'PMPOL:',':',.true.)) return
C   ... read nl and nsp ... abort if missing
        backspace ifi
        read(ifi,'(a72)') s
        if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
        if (nl .ne. nl2 .and. ipr .ge. 10)
     .  print *, 'aiomp (warning) mismatch in nl, class '//alabl
        if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
        if (nsp .ne. nsp2 .and. ipr .ge. 10)
     .  print *, 'aiomp (warning) mismatch in nsp, class '//alabl
C   ... For now, abort
        if (nsp .ne. nsp2 .or. nl .ne. nl2) goto 18
        call dpzero(pmpol,nl*nl*(lmxv+1)*3*nsp)
C   ... read pmpol
        n = min(nsp,nsp2)
        do  100  i = 1, n
        do  10  k = 1, 3

          do  120  l1 = 0, nl-1
          do  12  l2 = 0, nl-1

   16       if (.not. rdstrn(ifi,s,len(s),.false.)) goto 18
            if (s .eq. ' ') goto 16
            ip = 0
C    ...  Abort if fail to read lmxv+1 numbers
            ip = a2vec(s,len(s),ip,4,' ',1,-2,-(lmxv+3),ix2,x2)
            if (ip .ne. lmxv+3) call rxs('AIOMP: failed to parse ',s)
            call dcopy(lmxv+1,x2(3),1,pmpol(l1+1,l2+1,0,k,i),nl**2)
   12     continue
 120      continue
   10   continue
 100    continue
        aiomp = .true.
        return
   18   continue
        print *, 'aiomp: (input skipped) bad syntax, class '//alabl
        return

C --- File WRITE ---
      else
        write(-ifi,21) alabl, nl, nsp
        do  121  i = 1, nsp
        do  20  k = 1, 3

          do  122  l1 = 0, nl-1
          do  22  l2 = 0, nl-1
C           x2(1) = damax(lmxv+1,pmpol(l1+1,l2+1,0,k,i),nl*nl)
            ip = idamax(lmxv+1,pmpol(l1+1,l2+1,0,k,i),nl*nl)
            x2(1) = pmpol(l1+1,l2+1,ip-1,k,i)
            if (x2(1) .gt. 9999d0 .or. x2(1) .lt. -999d0) then
              call dcopy(lmxv+1,pmpol(l1+1,l2+1,0,k,i),nl*nl,x2,1)
C             call awrit4('%,4i%,4i %n:1g',' ',100,6,l1,l2,lmxv+1,x2)
              call awrit4('%,4i%,4i %n:1;11F',' ',100,-ifi,l1,l2,lmxv+1,
     .            x2)
            else
              write(-ifi,333) l1,l2,(pmpol(l1+1,l2+1,lf,k,i),lf=0,lmxv)
            endif
   22     continue
  122     continue
  333     format(2i4,10f12.6)
          if (i .ne. nsp .or. k .ne. 3) write(-ifi,333)

   20   continue
  121   continue
        aiomp = .true.
      endif
   21 format('PMPOL: ',a4,'  nl=',i1,'  nsp=',i1)

      end


      logical function aiorme(alabl,radme,nl,nsp,ifi)
C- File I/O for matrix elements of radial gradients of w.f.
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   alabl,nl,nsp
Ci   ifi:  logical unit: positive for read, negative for write
Cio  radme(4,2,nl,2,nsp):  <g2_l grad g1_l'> with l' = l +/- 1
Cio        i = 1 for g1=phi,    g2=phidot
Cio            2 for g1=phidot, g2=phi
Cio            3 for g1=phi,    g2=phidot
Cio            4 for g1=phidot, g2=phidot
Cio        j = 1 for <g2 | grad g1> - (l+1) < g2 | 1/r g1 >
Cio            2 for <g2 | grad g1> +    l  < g2 | 1/r g1 >
Cio       ll = 1 for l' = l-1  for l' = l+1
Cr Remarks
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*8 alabl
      integer nl,nsp,ifi
      double precision radme(4,2,nl,2,nsp)
C Local parameters
      integer is1,nl2,nsp2,ipr
      logical scat,pars1v,lfdmp
      character s*72

      aiorme = .false.
      call getpr(ipr)
      if (ifi .gt. 0) then
C   ... return unless file has RGRAD category
        if (.not. scat(ifi,'RGRAD:',':',.true.)) return
C   ... read nl and nsp ... abort if missing
        backspace ifi
        read(ifi,'(a72)') s
        if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
        if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
        if (nl .ne. nl2 .or. nsp .ne. nsp2)
     .  call rxs('AIORME: mismatch in nl or nsp, class ',alabl)
        call dpzero(radme,4*2*nl*2*nsp)
        do  10  is1 = 1, nsp
          aiorme = lfdmp(radme(1,1,1,1,is1),4*2*nl*2,ifi)
   10   continue
        return
   18   continue
        print *, 'aiorme: (input skipped) bad syntax, class '//alabl
      else
        write(-ifi,21) alabl, nl, nsp
        do  20  is1 = 1, nsp
          aiorme = lfdmp(radme(1,1,1,1,is1),4*2*nl*2,ifi)
   20   continue
      endif

   21 format('RGRAD: ',a4,'  nl=',i1,'  nsp=',i1)
      end
      
      logical function aiosop(alabl,sop,nl,lmax,nsp,ifi)
C- File I/O for spin-orbit coupling parameters.
C ----------------------------------------------------------------------
Ci Inputs
Ci   alabl :class label
Ci   nl    :(global maximum l) + 1
Ci   lmax  :maximum l for a given site
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio Inputs/Outputs
Cio  sop   :spin-orbit coupling parameters (made in soprm.f)
Cio        :sop(l+1,is1,is2,i) are matrix elts between spins is1 and is2
Cio        :for quantum number l. Six types of integrals are stored.
Cio        :(i=1) <phi so phi> (i=2) <phi so dot> (i=3) <dot so dot>
Cio        :(i=4) <phi || phi> (i=5) <phi || dot> (i=6) <dot || dot>
Cio        :The first three are the terms need to make SO perturbation
Cio        :the last three are used when external field is applied.
Cr Remarks
Cr
Cu Updates
Cu    4 Apr 04 Extended  to include matrix elements for XC field
Cu   07 Feb 03 Extended  to include matrix elements for applied field
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      character*8 alabl
      integer nl,lmax,nsp,ifi
      double precision sop(0:nl-1,nsp,nsp,9)
C Local parameters
      integer l,is1,is2,ks1,ks2,j,k,lmx,nl2,nsp2,ipr,lmxf
      logical scat,pars1v
      character s*72,nam*8

      aiosop = .false.
      call getpr(ipr)
      if (ifi .gt. 0) then
C   ... return unless file has SO category
        if (.not. scat(ifi,'SO:',':',.true.)) return
C   ... read nl and nsp ... abort if missing
        backspace ifi
        read(ifi,'(a72)') s
        call word(s,2,is1,is2)
        nam = s(is1:is2)
        if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
        if (min(nl,lmax+1) .ne. nl2 .and. ipr .ge. 10)
     .  print *, 'aiosop (warning) mismatch in nl, class '//alabl
        if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
        if (nsp .ne. nsp2 .and. ipr .ge. 10)
     .  print *, 'aiosop (warning) mismatch in nsp, class '//alabl
        read(ifi,'(a72)') s
        call dpzero(sop,nl*nsp*nsp*9)
        lmx = min(nl,nl2)-1
        lmxf = nl2-1
C   ... read SO parms
        do  10  l = 1, lmxf
          do  1101  is2 = 1, min(nsp,nsp2)
          do  11  is1 = 1, min(nsp,nsp2)
            if (l .le. lmx) then
              read(ifi,*) k,ks1,ks2,(sop(l,is1,is2,j), j=1,3)
            else
              read(ifi,*) k,ks1,ks2
            endif
            if (ks1 .ne. is1 .or. ks2 .ne. is2)
     .        call rx('aiosop: spin mismatch')
   11     continue
 1101     continue
   10   continue
        read(ifi,'(a72)') s
        do  100  l = 0, lmxf
          do  1111  is2 = 1, min(nsp,nsp2)
          do  111  is1 = 1, min(nsp,nsp2)
            if (l .le. lmx) then
              read(ifi,*) k,ks1,ks2,(sop(l,is1,is2,j), j=4,6)
            else
              read(ifi,*) k,ks1,ks2
            endif
            if (ks1 .ne. is1 .or. ks2 .ne. is2)
     .        call rxs('aiosop: spin mismatch, class ',nam)
  111     continue
 1111  continue
  100   continue
        read(ifi,'(a72)') s
        do  110  l = 0, lmxf
          do  1156  is2 = 1, min(nsp,nsp2)
          do  115  is1 = 1, min(nsp,nsp2)
            if (l .le. lmx) then
              read(ifi,*) k,ks1,ks2,(sop(l,is1,is2,j), j=7,9)
            else
              read(ifi,*) k,ks1,ks2
            endif
            if (ks1 .ne. is1 .or. ks2 .ne. is2)
     .        call rxs('aiosop: spin mismatch, class ',nam)
  115     continue
 1156     continue
  110   continue

        if (nsp2 .lt. nsp) then
          do   j = 1, 9
          do   is2 = 1, nsp
          do   is1 = 1, nsp
          do   l = 0, nl-1
             sop(l,is1,is2,j) = sop(l,1,1,j)
          enddo   
          enddo   
          enddo   
          enddo   
        endif

        aiosop = .true.
        return
   18   continue
        print *, 'aiosop: (input skipped) bad syntax, class '//alabl
      else
        write(-ifi,21) alabl, lmax+1, nsp
        do   l = 1, lmax
        do   is2 = 1, nsp
        do   is1 = 1, nsp
          write(-ifi,333) l,is1,is2,(sop(l,is1,is2,j), j=1,3)
  333     format(i4,2i3,3f15.10)
        enddo   
        enddo   
        enddo   
C       write(-ifi,'(''   matrix elements of phi and phidot'')')
        write(-ifi,121)
        do  l = 0, lmax
        do  is2 = 1, nsp
        do  is1 = 1, nsp
          write(-ifi,333) l,is1,is2,(sop(l,is1,is2,j), j=4,6)
        enddo   
        enddo   
        enddo   
C       write(-ifi,'(''   matrix elements of phi B phidot'')')
        write(-ifi,221)
        do   l = 0, lmax
        do   is2 = 1, nsp
        do   is1 = 1, nsp
          write(-ifi,333) l,is1,is2,(sop(l,is1,is2,j), j=7,9)
        enddo   
        enddo   
        enddo   

        aiosop = .true.
      endif

   21 format('SO: ',a4,'  nl=',i1,'  nsp=',i1/
     .'   l is js  < phi so phi>  < phi so dot>  < dot so dot>')
  121 format(
     .'   l is js  < phi || phi>  < phi || dot>  < dot || dot>')
  221 format(
     .'   l is js  <phi|Bxc|phi>  <phi|Bxc|dot>  <dot|Bxc|dot>')
      end

      logical function aiova(alabl,vintra,nl,lmax,nsp,ifi)
C- File I/O for ASA intra-atomic density-density response matrix
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   alabl,nl,lmax,nsp
Ci   ifi:    logical unit: positive for read, negative for write
Cio  vintra: intra-atomic density-density response matrix
Cr Remarks
Cr Bugs
Cr   Input not checked when file nsp mismatches passed nsp
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*8 alabl
      integer nl,lmax,nsp,ifi
      double precision vintra(nl,nl,nsp,nsp)
C Local parameters
      double precision x2(0:10)
      integer i,i2,l,n,ip,k,a2vec,ix2(7),ipr,nl2,nsp2,lmx
      logical scat,rdstrn,pars1v
      character*120 s

      aiova = .false.
      call getpr(ipr)
C --- File READ ---
      if (ifi .gt. 0) then
C   ... return unless file has VINTRA: category
        if (.not. scat(ifi,'VINTRA:',':',.true.)) return
C   ... read nl and nsp ... abort if missing
        backspace ifi
        read(ifi,'(a120)') s
        if (.not. pars1v(s,len(s),'nl=','=',2,nl2)) goto 18
        if (lmax+1 .ne. nl2 .and. ipr .ge. 10)
     .  print *, 'aiova (warning) mismatch in nl, class '//alabl
        if (.not. pars1v(s,len(s),'nsp=','=',2,nsp2)) goto 18
        if (nsp .ne. nsp2 .and. ipr .ge. 10)
     .  print *, 'aiova (warning) mismatch in nsp, class '//alabl
        call dpzero(vintra,nl*nl*nsp*nsp)
        lmx = min(nl,nl2)-1
C   ... read vintra
        n = min(nsp,nsp2)
        do   i = 1, n
        do   l = 0, lmx
          if (.not. rdstrn(ifi,s,len(s),.false.)) goto 18
          ip = 0
          k = a2vec(s,len(s),ip,4,' ',1,-2,-(lmx+1)*n-1,ix2,x2)
C    ...  Abort if failed to read lmx+2 numbers
          if (k .ne. (lmx+1)*n+1) call rx('AIOVA: failed to parse '//s)
          call dcopy(lmx+1,x2(1),1,vintra(l+1,1,i,1),nl)
          if (n.eq.2) call dcopy(lmx+1,x2(lmx+2),1,vintra(l+1,1,i,2),nl)
          if (nsp2 .lt. nsp) then
            call dcopy(lmx+1,x2(1),1,vintra(l+1,1,2,2),nl)
            call dcopy(lmx+1,x2(1),1,vintra(l+1,1,1,2),nl)
            call dcopy(lmx+1,x2(1),1,vintra(l+1,1,2,1),nl)
          endif
        enddo
        enddo
        aiova = .true.
        return
   18   continue
        print *, 'aiova: (input skipped) bad syntax, class '//alabl
        return

C --- File WRITE ---
      else
        write(-ifi,21) alabl, lmax+1, nsp
        do  i = 1, nsp
        do  l = 0, lmax
          write(-ifi,333) l,((vintra(l+1,k,i,i2), k=1,1+lmax), i2=1,nsp)
 333      format(i4,8f12.6)
        enddo
        enddo
        aiova = .true.
      endif
   21 format('VINTRA: ',a4,'  nl=',i1,'  nsp=',i1)
      end
      end module m_lmaux
