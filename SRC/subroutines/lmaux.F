      subroutine lmaux() !main part of lmchk
      use m_mksym,only: ctrl_nclass,rv_a_oclabl,iv_a_oics,iv_a_oipc
      use m_lmfinit,only: iv_a_oips,str_mxnbr,str_rmax,ctrl_nbas,ctrl_nspec,ctrl_nspin,
     &     ctrl_nl,ctrl_modep,ctrl_omax1,ctrl_omax2,ctrl_wsrmax,rv_a_opos,slabl,sspec=>v_sspec,
     &     lat_avw,lat_alat,lat_nkd,lat_nkq
      use m_struc_def  !Cgetarg
      use m_pairs
      use m_ovmin , only: ovmin
      use m_lattic,only:lat_plat
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
      real(8),allocatable:: eula(:,:)
      integer,allocatable :: iv_a_ontab(:)
      integer,allocatable :: iv_a_oiax(:)
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
      subroutine ioqpp(lio)     !,spot)!,sctrl
      use m_clsprm,only:rv_a_oqpp
      use m_struc_def           !Cgetarg
      use m_lmfinit,only:   ctrl_nl,ctrl_nspin,ctrl_nbas
      use m_ext,only:sname
C     - File i/o for phi-phi, phi-dot, dot-dot products
C     ----------------------------------------------------------------------
C     i Inputs:
C     i   lio: true for write, false for read
C     i          <0 write
C     u Updates
C     u   08 Nov 07 (J. Xu) qpp is complex
C     ----------------------------------------------------------------------
      implicit none
C     ... Passed parameters
      logical lio
C     getarg       double precision sctrl(1),spot(1)
c      type(s_ctrl)::sctrl
c      type(s_pot)::spot
C     ... Local parameters
      integer:: i1 , i2 , nbas , nl , nsp , igets , ifi , rdm 
     .     , ipr ,ifile_handle
      call getpr(ipr)
      nbas=ctrl_nbas
      nl=ctrl_nl
      nsp=ctrl_nspin
      i1 = nl**2*(nl**2+1)
      i2 = 4*nsp*nbas
      ifi= ifile_handle()
      open(ifi,file='qpp.'//trim(sname),form='unformatted')
      if (lio) then
         call ywrm ( 1 , 'lmasa' , 1 , ifi , ' ' , rv_a_oqpp , 1 , i1 
     .        , i1 , i2 ) 
         if (ipr .ge. 30) print *, 'IOQPP:  wrote qpp to disk'
      else
         if (ipr .ge. 30) print *, 'IOQPP:  reading qpp from disk ...'
         call pshpr(0)
         if ( rdm ( ifi , 2 , i1 * i2 , ' ' , rv_a_oqpp , i1 , i2 ) .lt. 
     .        0 ) then 
            if (ipr .ge. 0) print *,'IOQPP:  (warning) failed to qpp file'
            call dvset ( rv_a_oqpp , 1 , 1 , - 1d0 ) 
         endif
         call poppr
      endif
      close(ifi)
      end subroutine ioqpp 


      


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
c     #define ONE_ATOM_UNIT 1
      
      subroutine aiocls(lio,mode,sspec,ic1,ic2) !sctrl,sham,,slat spot,sarray,
      use m_clsprm,only: iv_a_oipqx,
     &     rv_a_obxc , rv_a_ogrrme ,  rv_a_opmpol , rv_a_opnu , rv_a_opp , rv_a_opprel,
     &     rv_a_oqnu , rv_a_oqpp , rv_a_oqt , rv_a_orhrmx ,rv_a_opnu,
     &     rv_a_osop , rv_a_ovdif , rv_a_oves , rv_a_ovintr , rv_a_ovrmax 
      use m_mksym,only: rv_a_oclabl,iv_a_oics,ctrl_nclass
      use m_struc_def           !Cgetarg
      use m_lmfinit,only: rv_a_opos, ctrl_nl,ctrl_nspin,lrel_g=>lrel
      use m_supot,only:rv_a_omad
      use m_ext,only:sname
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
      logical sw,aiomom,aiopar,aiopot,aiova,lpot,lcor,scat,lgors,
     .     aiogen,aiosop,aiorme,aiocor,aiomp,lrell,lgen
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


