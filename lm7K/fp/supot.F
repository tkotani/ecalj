Cchp1 cmdvar= complex(8),zv_p_,osmrho
Cchp1 cmdvar= complex(8),zv_p_,osmpot
Cchp1 cmdvar= real(8),iv_p_,oips0
Cchp1 cmdvar= real(8),rv_p_,oqlv
#define FP 1
      subroutine supot(mode,sctrl,slat,spot)

      use m_struc_def  !Cgetarg

C- Initialization for potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1 make Madelung matrix for monopoles (ASA)
Ci   sctrl :struct containing parameters governing program flow
Ci     Elts read: nbas nbasp nspin lcd,4
Ci   slat  :struct containing information about the lattice
Ci     Elts read: alat vol opos awald nkd nkq odlv oqlv plat nabc gmax
Ci                nsgrp osymgr oag
Ci     Stored:    ng ogv okv gmax nabc oips0 obgv
Ci   spot  :struct containing information about the potential
Ci     Stored:    omad osmrho osmpot oorhat
Co Outputs
Cr Remarks
Cr   Permanent arrays are created to handle the potential, depending
Cr   on the mode
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode
      type(s_ctrl)::sctrl
      type(s_lat)::slat
      type(s_pot)::spot

C ... Local parameters
Cchp1       integer nbas,nbasp,nbaspp,nsp,omad,opos,nkd,nkq,odlv,oqlv,igets,
Cchp1      .ngabc(3),n1,n2,n3,k1,k2,k3,ngmx,ogv,okv,ng,osmrho,osmpot,oorhat,
Cchp1      .oips0,obgv,ngrp,og,oag,iprint
Cchp1        integer:: nbas , nbasp , nbaspp , nsp , omad , opos , nkd , nkq 
Cchp1      . , odlv , igets , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx 
Cchp1      . , ogv , okv , ng , osmrho , osmpot , oorhat , oips0 , obgv , 
Cchp1      . ngrp , og , oag , iprint 
Cchp1         integer:: nbas , nbasp , nbaspp , nsp , omad , opos , nkd , nkq 
Cchp1      .  , odlv , igets , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx 
Cchp1      .  , ogv , okv , ng , osmrho , osmpot , oorhat , obgv , ngrp , og 
Cchp1      .  , oag , iprint 
Cchp1          integer:: nbas , nbasp , nbaspp , nsp , omad , opos , nkd , nkq 
Cchp1      .   , odlv , igets , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx 
Cchp1      .   , ogv , okv , ng , osmrho , oorhat , obgv , ngrp , og , oag , 
Cchp1      .   iprint 
          integer:: nbas , nbasp , nbaspp , nsp , omad , opos , nkd , nkq 
     .    , odlv , igets , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx 
     .    , ogv , okv , ng , oorhat , obgv , ngrp , og , oag , iprint 
          complex(8),pointer :: zv_p_osmrho(:)

         complex(8),pointer :: zv_p_osmpot(:)

        integer,pointer :: iv_p_oips0(:)

       real(8),pointer :: rv_p_oqlv(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision awald,alat,vol,plat(9),gmax,xx
C ... Heap
      integer w(1)
      common /w/ w
      call info0(50,1,0,' supot : potential setup ... ')


      nbas=sctrl%nbas
      nbasp=sctrl%nbasp
      nsp=sctrl%nspin


      alat=slat%alat
      vol=slat%vol
      opos=slat%opos


      awald=slat%awald
      nkd=slat%nkd
      nkq=slat%nkq
      odlv=slat%odlv
Cchp1       oqlv=slat%oqlv
       rv_p_oqlv => slat%rv_p_oqlv 



      nbaspp = 2*nbasp - nbas

C --- Monopole Madelung matrix ---
      if (mod(mode,10) .eq. 1) then
        call defdr(omad,nbaspp**2)
Cchp1         call madmat(nbaspp,w(opos),awald,alat,vol,w(odlv),nkd,
Cchp1      .  w(oqlv),nkq,w(omad))
         call madmat ( nbaspp , w ( opos ) , awald , alat , vol , w ( 
     .   odlv ) , nkd , rv_p_oqlv , nkq , w ( omad ) ) 


        spot%omad=omad

      endif

#if FP
C --- Setup for FT charge density, potential representation ---
      if ( iand(4,int(sctrl%lcd)) .ne. 0 ) then


        alat=slat%alat
        i_copy_size=size(slat%plat)
        call dcopy(i_copy_size,slat%plat,1,plat,1)
        i_copy_size=size(slat%nabc)
        call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
        gmax=slat%gmax


C   ... Generate energy cutoff gmax or n1..n3
c        call pshpr(iprint()-0)
        call mshsiz(alat,plat,0,gmax,ngabc,ngmx)
        call fftz30(n1,n2,n3,k1,k2,k3)

C   ... Make list of lattice vectors within cutoff
        call defrr(ogv, ngmx*3)
        call defi (okv, ngmx*3)
        call gvlist(alat,plat,w,n1,n2,n3,gmax,8,ngmx,ng,w(okv),
     .  w(ogv),xx,xx)
        if (ng .ne. ngmx) call rx('supot: bug in gvlist')


        slat%ng=ng
        slat%ogv=ogv
        slat%okv=okv
        slat%gmax=gmax
        i_copy_size=size(slat%nabc)
        call ii8copy(i_copy_size,ngabc,1,slat%nabc,1)

C       call prmx('gv',w(ogv),ngmx,ng,3)

Cchp1         call defcc(osmrho, -k1*k2*k3*nsp)
         allocate(zv_p_osmrho(abs(-k1*k2*k3*nsp)))
         if (-k1*k2*k3*nsp<0) zv_p_osmrho(:)=0.0d0

Cchp1         call defcc(osmpot, -k1*k2*k3*nsp)
         allocate(zv_p_osmpot(abs(-k1*k2*k3*nsp)))
         if (-k1*k2*k3*nsp<0) zv_p_osmpot(:)=0.0d0

        call defi(oorhat,3*nbaspp)

Cchp1         spot%osmrho=osmrho
         spot%zv_p_osmrho => zv_p_osmrho 

Cchp1         spot%osmpot=osmpot
         spot%zv_p_osmpot => zv_p_osmpot 

        spot%oorhat=oorhat


C   ... Setup for symmetrization
Cchp1         call defi (oips0,   ng)
         allocate(iv_p_oips0(abs(ng)))
         if (ng<0) iv_p_oips0(:)=0.0d0

        call defcc(obgv,    ng)

        ngrp=slat%nsgrp
        og=slat%osymgr
        oag=slat%oag

Cchp1         call sgvsym(ngrp,w(og),w(oag),ng,w(ogv),w(oips0),w(obgv))
         call sgvsym ( ngrp , w ( og ) , w ( oag ) , ng , w ( ogv ) , 
     .   iv_p_oips0 , w ( obgv ) ) 


Cchp1         slat%oips0=oips0
         slat%iv_p_oips0 => iv_p_oips0 

        slat%obgv=obgv


      endif
#endif
      end subroutine supot


