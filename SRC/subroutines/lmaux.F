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
csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      module m_pairs
      public pairs,pairc,ppair1
      private
      contains
      subroutine pairs(nbas,nbasp,alat,plat,rmax,baspp,ipsp,nd,iltab,
     .pltab,nttab,iv_a_ontab,iv_a_oiax,mxcsiz)
C- Allocate memory for and create neighbor table for a crystal
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis (input)
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
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
      nbaspp = 2*nbasp - nbas
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
*     call prmx('plat',plat,3,3,3)
*     call prmx('starting pos',pos,3,3,nbasp)
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
c        allocate(wk1_rv((nttab-mtab)*niax))
c        allocate(wk2_rv((nttab-mtab)))
c        allocate(wk3_rv((nttab-mtab)*3))
        call ppair3 ( nttab - mtab , iax(1, mtab ), rtab(1 , mtab))! , wk1_rv , wk2_rv , wk3_rv )
c        if (allocated(wk3_rv)) deallocate(wk3_rv)
c        if (allocated(wk2_rv)) deallocate(wk2_rv)
c        if (allocated(wk1_rv)) deallocate(wk1_rv)
C --- Cleanup for this ib ---
        mtab = nttab
        ntab(ib+1) = nttab-1
        mxcsiz = max(mxcsiz,ntab(ib+1)-ntab(ib))
   10 continue
      nttab = nttab-1
C     call awrit2('xx ntab %n:1i',' ',80,6,nbasp+1,ntab)

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

C     call tcx('ppair1')

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

*        call awrit5('try ib,jb,ilat= %i %i %i rr=%;4d: %l',' ',80,
*     .    6,ib,jb,ilat,rr,rr.lt.rcut)

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

C ... Debugging ...
C      call ivprm(niax,nclus-iclus+1,iax(1,iclus),wk2,iwk,.false.)
C      call yprm('iax',0,wk2,0,niax,niax,nclus-iclus+1)
C      call awrit2('iwk %n:1i',' ',180,6,nclus-iclus+1,iwk)
C      do  30  ic = iclus, nclus
C      ic0 = ic-iclus+1
C   30 iwk(ic0) = iwk(ic0)-1
C      call dvperm(3,nclus-iclus+1,rtab,wk2,iwk,.true.)
C      do  32  ic = iclus, nclus
C        ic0 = ic-iclus+1
C        print 346, ic,(rtab(jx,ic0)+pos(jx,iax(1,ic)), jx=1,3)
C  346   format(i4,3f11.6)
C   32 continue
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

C heap:

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
c$$$      subroutine pair3c(offi,nvec,iaxc,iaxh,ntabh,ipiaxh)
c$$$C- Find connecting vectors within a cluster linked by a pair table
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   offi :offset to starting point in iaxc for this cluster
c$$$Ci   nvec :number of sites in cluster for which to find connecting vectors
c$$$Ci   iaxc :cluster pair table.  Only iaxc(2) and iaxc(7) are used.
c$$$Ci         iaxc(2,*) = field site index for this table
c$$$Ci         iaxc(7,*) = list of permutations that orders sites
c$$$Ci                     in iaxc and iaxh in the same way; see Remarks.
c$$$Ci   iaxh :pair (eg hamiltonian) neighbor table (pairc.f).
c$$$Ci         Only iaxh(2) and iaxh(7) are used. iaxh must consist of a
c$$$Ci         series of clusters grouped by site index, as in the
c$$$Ci         standard formed described in pairc.  The site index
c$$$Ci         (kept in iaxh(1,*) in the standard form) and the site
c$$$Ci         index iaxc(2,*) must refer to the same basis atom.
c$$$Ci   ntabh:ntabh(ib)=number of pairs in iaxh preceding site ib.
c$$$Co Outputs
c$$$Co  ipiaxh:indices to connecting vectors within a cluster.
c$$$Co         ipaxph(ip,kp) points to entry in iaxh table connecting
c$$$Co         ip to kp. ipaxph(ip,kp)=0 => no such vector in iaxh exists.
c$$$Cr Remarks
c$$$Cr   pair3c generates in ipiaxh the information needed to connect
c$$$Cr   three centers in a crystal.  For each pair of connecting vectors
c$$$Cr   in iaxc, a third vector is sought in iaxh that connects the two.
c$$$Cr
c$$$Cr   A schematic depiction of association of variables with their sites:
c$$$Cr
c$$$Cr            . origin
c$$$Cr
c$$$Cr    .  jp=field point
c$$$Cr    \
c$$$Cr     \ <- connecting vector ipiaxh(jp,ip) = entry in iaxh
c$$$Cr      \
c$$$Cr       \.  ip, source point
c$$$Cr
c$$$Cr   For efficient implementation, this routine relies on sorting
c$$$Cr   information in both iaxc and iaxh.  iaxc(7,*) and iaxh(7,*) must
c$$$Cr   hold a list of permutations that sorts each cluster in a unique
c$$$Cr   manner, so that connecting vectors common to different clusters are
c$$$Cr   ordered in the same way in those two clusters; this is done in
c$$$Cr   in ppair4, called with tol<0.
c$$$C ----------------------------------------------------------------------
c$$$      implicit none
c$$$      integer offi,nvec,niax,ntabh(8),ipiaxh(nvec,nvec)
c$$$      parameter (niax=10)
c$$$      integer iaxc(niax,1),iaxh(niax,9)
c$$$      integer ii,iii,jj,jjj,ip,jp,ib,nc,ich,icp,ic0
c$$$      integer iiax(5)
c$$$      call iinit(ipiaxh,nvec*nvec)
c$$$C --- Loop over 3C pairs i,j in permuted order ---
c$$$      ii = 0
c$$$      do  30  iii = 1, nvec
c$$$   31   ii = ii+1
c$$$        ip = iaxc(7,ii+offi)
c$$$        if (ip .gt. nvec) goto 31
c$$$        ib = iaxc(2,ip+offi)
c$$$C   ... ic0 = offset corresponding to iaxh(ic) = 1st pair for ib
c$$$        ic0 = ntabh(ib)
c$$$C   ... ich = current offset corresponding to iaxh.  Because we loop
c$$$C       in sorted order, we need not reset ich for each new pair below.
c$$$        ich = ic0+1
c$$$        nc  = ntabh(ib+1)
c$$$        jj = 0
c$$$        do  40  jjj = 1, nvec
c$$$   41     jj = jj+1
c$$$          jp = iaxc(7,jj+offi)
c$$$          if (jp .gt. nvec) goto 41
c$$$C         Do upper triangle only
c$$$C         if (lx6 .and. jp .lt. ip) goto 40
c$$$          iiax(3) = iaxc(3,jp+offi) - iaxc(3,ip+offi)
c$$$          iiax(4) = iaxc(4,jp+offi) - iaxc(4,ip+offi)
c$$$          iiax(5) = iaxc(5,jp+offi) - iaxc(5,ip+offi)
c$$$
c$$$C     ... Increment ich until no exhaust all pairs in iaxh for ib
c$$$C          if (jp .eq. 74 .and. ip .eq. 7) then
c$$$C            print *, 'test'
c$$$C          endif
c$$$          ich = ich-1
c$$$   42     ich = ich+1
c$$$C     ... Skip if no more vectors connecting this (ib,jb)
c$$$          if (ich .gt. nc) goto 40
c$$$          icp = ic0+iaxh(7,ich)
c$$$
c$$$C     ... Increment jp until iaxh is at least as large as iaxc
c$$$C         and ich  until iaxc is at least as large as iaxh
c$$$          if (iaxc(2,offi+jp) .gt. iaxh(2,icp)) goto 42
c$$$          if (iaxc(2,offi+jp) .lt. iaxh(2,icp)) goto 40
c$$$          if (iiax(3) .gt. iaxh(3,icp)) goto 42
c$$$          if (iiax(3) .lt. iaxh(3,icp)) goto 40
c$$$          if (iiax(4) .gt. iaxh(4,icp)) goto 42
c$$$          if (iiax(4) .lt. iaxh(4,icp)) goto 40
c$$$          if (iiax(5) .gt. iaxh(5,icp)) goto 42
c$$$          if (iiax(5) .lt. iaxh(5,icp)) goto 40
c$$$
c$$$C     ... A match was found
c$$$          ipiaxh(jp,ip) = icp
c$$$   40   continue
c$$$   30 continue
c$$$
c$$$C --- Construct lower triangle of table ---
c$$$C      if (lx6) then
c$$$C        do  50  i = 1, nvec
c$$$C        do  50  j = i+1, nvec
c$$$C          is = ipiaxh(i,j)
c$$$C          if (is .gt. 0) ipiaxh(j,i) = iaxh(6,is)
c$$$C   50   continue
c$$$C      endif
c$$$
c$$$      end
      end module m_pairs

!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine makrm0(opt,nspec,nbas,alat,plat,pos,slabl,ips,modep,
     .lock,z,rmt)
      use m_lmfinit,only: nsp
      use m_pairs,only:pairc

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
      use m_pairs,only: pairc
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
C     implicit none
C Passed variables:
      integer opt,nbas,nbasp,nspec,ips(nbas),modep(3),lock(nspec)
      double precision alat,basp(3,nbasp),volfac,wsmax,
     .omax1(3),omax2(3),plat(3,3),wsr(nspec),z(nspec)
      character*8 slabl(*)
C Local variables:
      integer:: i , k , nrspec(nspec) , ib , is , opt1
      real(8) ,allocatable :: wk_rv(:)

      integer niax,mxcsiz,mxnbr,nttab,ipr,stdo,lgunit,llock(nspec)
      integer ,allocatable :: ntab_iv(:)
      integer ,allocatable :: iax_iv(:)

      logical les
      double precision dovl1(3),dovl2(3),gamma,range(nbas),wsrs(nspec)
      double precision facr,tiny,avw,vol,avwsr,volnew,volnes,volold,
     .volsph
      parameter(niax=10,mxcsiz=200,facr=2d0,tiny=1d-5)
C ... Heap

      integer:: istdo,fopna

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
      call maknrs(nbas,ips,ib,nrspec)
      if (ib .gt. nspec) call rx('sclwsr: wrong number of species')
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
            call awrit2(' SCLWSR:  initial sphere packing = %;1d%%'//
     .      ' scaled to %;1d%% (no empty spheres)',
     .      ' ',120,stdo,100*volold,100*volnes)
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
        call awrit3(' SCLWSR:  vol = %,5g a.u..%3f'//
     .  'Initial sphere packing = %;1d%%  '//
     .  'scaled to %;1d%%',
     .  ' ',120,stdo,vol,100*volold,100*volnew)
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

      istdo = fopna('rmt',-1,2)
      do  i = 1, nspec
        write(istdo,"(a,f12.6)") slabl(i), wsr(i)
      enddo
      call fclose(istdo)

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


