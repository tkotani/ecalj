Cdwdef1 have_common_w_in supot , #_of_w_access= 1
Cchp1 cmdvar= real(8),rv_p_,opos
Cchp1 cmdvar= real(8),rv_p_,og
Cchp1 cmdvar= real(8),rv_p_,osymgr
Cchp1 cmdvar= integer,iv_p_,okv
Cchp1 cmdvar= real(8),rv_p_,ogv
Cchp1 cmdvar= real(8),rv_p_,odlv
Cchp1 cmdvar= complex(8),zv_p_,obgv
Cchp1 cmdvar= real(8),rv_p_,oag
Cchp1 cmdvar= integer,iv_p_,oorhat
Cchp1 cmdvar= real(8),rv_p_,omad
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
Cchp1           integer:: nbas , nbasp , nbaspp , nsp , omad , opos , nkd , nkq 
Cchp1      .    , odlv , igets , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx 
Cchp1      .    , ogv , okv , ng , oorhat , obgv , ngrp , og , oag , iprint 
Cchp1            integer:: nbas , nbasp , nbaspp , nsp , opos , nkd , nkq , odlv 
Cchp1      .     , igets , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ogv 
Cchp1      .     , okv , ng , oorhat , obgv , ngrp , og , oag , iprint 
Cchp1             integer:: nbas , nbasp , nbaspp , nsp , opos , nkd , nkq , odlv 
Cchp1      .      , igets , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ogv 
Cchp1      .      , okv , ng , obgv , ngrp , og , oag , iprint 
Cchp1              integer:: nbas , nbasp , nbaspp , nsp , opos , nkd , nkq , odlv 
Cchp1      .       , igets , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ogv 
Cchp1      .       , okv , ng , obgv , ngrp , og , iprint 
Cchp1               integer:: nbas , nbasp , nbaspp , nsp , opos , nkd , nkq , odlv 
Cchp1      .        , igets , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ogv 
Cchp1      .        , okv , ng , ngrp , og , iprint 
Cchp1                integer:: nbas , nbasp , nbaspp , nsp , opos , nkd , nkq , igets 
Cchp1      .         , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ogv , okv , 
Cchp1      .         ng , ngrp , og , iprint 
Cchp1                 integer:: nbas , nbasp , nbaspp , nsp , opos , nkd , nkq , igets 
Cchp1      .          , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , okv , ng , 
Cchp1      .          ngrp , og , iprint 
Cchp1                  integer:: nbas , nbasp , nbaspp , nsp , opos , nkd , nkq , igets 
Cchp1      .           , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ng , ngrp , 
Cchp1      .           og , iprint 
Cchp1                   integer:: nbas , nbasp , nbaspp , nsp , opos , nkd , nkq , igets 
Cchp1      .            , ngabc(3) , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ng , ngrp , 
Cchp1      .            iprint 
                   integer:: nbas , nbasp , nbaspp , nsp , nkd , nkq , igets , ngabc(3) 
     .             , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ng , ngrp , iprint 
                   real(8),pointer :: rv_p_opos(:)

                  real(8),pointer :: rv_p_og(:)

                 integer,pointer :: iv_p_okv(:)

                real(8),pointer :: rv_p_ogv(:)

               real(8),pointer :: rv_p_odlv(:)

              complex(8),pointer :: zv_p_obgv(:)

             real(8),pointer :: rv_p_oag(:)

            type(s_rv1),pointer :: sv_p_oorhat(:)

           real(8),pointer :: rv_p_omad(:)

          complex(8),pointer :: zv_p_osmrho(:)

         complex(8),pointer :: zv_p_osmpot(:)

        integer,pointer :: iv_p_oips0(:)

       real(8),pointer :: rv_p_oqlv(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision awald,alat,vol,plat(9),gmax,xx
C ... Heap
Cdwdef1       integer w(1)
Cdwdef1       common /w/ w
       integer ::iwdummy 

      call info0(50,1,0,' supot : potential setup ... ')


      nbas=sctrl%nbas
      nbasp=sctrl%nbasp
      nsp=sctrl%nspin


      alat=slat%alat
      vol=slat%vol
Cchp1       opos=slat%opos
       rv_p_opos => slat%rv_p_opos 



      awald=slat%awald
      nkd=slat%nkd
      nkq=slat%nkq
Cchp1       odlv=slat%odlv
       rv_p_odlv => slat%rv_p_odlv 

Cchp1       oqlv=slat%oqlv
       rv_p_oqlv => slat%rv_p_oqlv 



      nbaspp = 2*nbasp - nbas

C --- Monopole Madelung matrix ---
      if (mod(mode,10) .eq. 1) then
Cchp1         call defdr(omad,nbaspp**2)
         allocate(rv_p_omad(nbaspp**2))
         rv_p_omad(:)=0.0d0

Cchp1         call madmat(nbaspp,w(opos),awald,alat,vol,w(odlv),nkd,
Cchp1      .  w(oqlv),nkq,w(omad))
Cchp1          call madmat ( nbaspp , w ( opos ) , awald , alat , vol , w ( 
Cchp1      .   odlv ) , nkd , rv_p_oqlv , nkq , w ( omad ) ) 
Cchp1           call madmat ( nbaspp , w ( opos ) , awald , alat , vol , w ( 
Cchp1      .    odlv ) , nkd , rv_p_oqlv , nkq , rv_p_omad ) 
Cchp1            call madmat ( nbaspp , w ( opos ) , awald , alat , vol , rv_p_odlv 
Cchp1      .     , nkd , rv_p_oqlv , nkq , rv_p_omad ) 
            call madmat ( nbaspp , rv_p_opos , awald , alat , vol , rv_p_odlv 
     .      , nkd , rv_p_oqlv , nkq , rv_p_omad ) 





Cchp1         spot%omad=omad
         spot%rv_p_omad => rv_p_omad 


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
Cchp1         call defrr(ogv, ngmx*3)
         allocate(rv_p_ogv(abs(ngmx*3)))
c         if (ngmx*3<0) rv_p_ogv(:)=0.0d0

Cchp1         call defi (okv, ngmx*3)
         allocate(iv_p_okv(abs(ngmx*3)))
c         if (ngmx*3<0) iv_p_okv(:)=0

Cchp1         call gvlist(alat,plat,w,n1,n2,n3,gmax,8,ngmx,ng,w(okv),
Cchp1      .  w(ogv),xx,xx)
Cchp1          call gvlist ( alat , plat , w , n1 , n2 , n3 , gmax , 8 , ngmx 
Cchp1      .   , ng , w ( okv ) , rv_p_ogv , xx , xx ) 
Cdwdef1           call gvlist ( alat , plat , w , n1 , n2 , n3 , gmax , 8 , ngmx 
Cdwdef1      .    , ng , iv_p_okv , rv_p_ogv , xx , xx ) 
           call gvlist ( alat , plat , iwdummy , n1 , n2 , n3 , gmax , 8 
     .     , ngmx , ng , iv_p_okv , rv_p_ogv , xx , xx ) 



        if (ng .ne. ngmx) call rx('supot: bug in gvlist')


        slat%ng=ng
Cchp1         slat%ogv=ogv
         slat%rv_p_ogv => rv_p_ogv 

Cchp1         slat%okv=okv
         slat%iv_p_okv => iv_p_okv 

        slat%gmax=gmax
        i_copy_size=size(slat%nabc)
        call ii8copy(i_copy_size,ngabc,1,slat%nabc,1)

C       call prmx('gv',w(ogv),ngmx,ng,3)

Cchp1         call defcc(osmrho, -k1*k2*k3*nsp)
         allocate(zv_p_osmrho(k1*k2*k3*nsp))
c         if (-k1*k2*k3*nsp<0) zv_p_osmrho(:)=0.0d0
         zv_p_osmrho(:)=0d0

Cchp1         call defcc(osmpot, -k1*k2*k3*nsp)
         allocate(zv_p_osmpot(k1*k2*k3*nsp))
         zv_p_osmpot(:)=0.0d0

Cchp1         call defi(oorhat,3*nbaspp)
         allocate(sv_p_oorhat(3*nbaspp))
cki         sv_p_oorhat(:)=0
         do i=1,3*nbaspp
            nullify(sv_p_oorhat(i)%v)
         enddo


Cchp1         spot%osmrho=osmrho
         spot%zv_p_osmrho => zv_p_osmrho 

Cchp1         spot%osmpot=osmpot
         spot%zv_p_osmpot => zv_p_osmpot 

Cchp1         spot%oorhat=oorhat
         spot%sv_p_oorhat => sv_p_oorhat 



C   ... Setup for symmetrization
Cchp1         call defi (oips0,   ng)
         allocate(iv_p_oips0(ng))
         iv_p_oips0(:)=0.0d0

Cchp1         call defcc(obgv,    ng)
         allocate(zv_p_obgv(ng))
         zv_p_obgv(:)=0.0d0


        ngrp=slat%nsgrp
Cchp1         og=slat%osymgr
         rv_p_og => slat%rv_p_osymgr 

Cchp1         oag=slat%oag
         rv_p_oag => slat%rv_p_oag 


Cchp1         call sgvsym(ngrp,w(og),w(oag),ng,w(ogv),w(oips0),w(obgv))
Cchp1          call sgvsym ( ngrp , w ( og ) , w ( oag ) , ng , w ( ogv ) , 
Cchp1      .   iv_p_oips0 , w ( obgv ) ) 
Cchp1           call sgvsym ( ngrp , w ( og ) , rv_p_oag , ng , w ( ogv ) , iv_p_oips0 
Cchp1      .    , w ( obgv ) ) 
Cchp1            call sgvsym ( ngrp , w ( og ) , rv_p_oag , ng , w ( ogv ) , iv_p_oips0 
Cchp1      .     , zv_p_obgv ) 
Cchp1             call sgvsym ( ngrp , w ( og ) , rv_p_oag , ng , rv_p_ogv , iv_p_oips0 
Cchp1      .      , zv_p_obgv ) 
             call sgvsym ( ngrp , rv_p_og , rv_p_oag , ng , rv_p_ogv , iv_p_oips0 
     .       , zv_p_obgv ) 






Cchp1         slat%oips0=oips0
         slat%iv_p_oips0 => iv_p_oips0 

Cchp1         slat%obgv=obgv
         slat%zv_p_obgv => zv_p_obgv 



      endif
#endif
      end subroutine supot


