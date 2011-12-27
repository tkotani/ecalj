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
                   integer:: nbas , nbasp , nbaspp , nsp , nkd , nkq , igets , ngabc(3) 
     .             , n1 , n2 , n3 , k1 , k2 , k3 , ngmx , ng , ngrp , iprint 
ckino Dec.8.2011:                     real(8),pointer :: rv_p_opos(:) =>NULL()

ckino Dec.13.2011:                    real(8),pointer :: rv_p_og(:) =>NULL()

ckino Dec.8.2011:                   integer,pointer :: iv_p_okv(:) =>NULL()

ckino Dec.8.2011:                  real(8),pointer :: rv_p_ogv(:) =>NULL()

ckino Dec.8.2011:                 real(8),pointer :: rv_p_odlv(:) =>NULL()

ckino Dec.8.2011:                complex(8),pointer :: zv_p_obgv(:) =>NULL()

ckino Dec.8.2011:               real(8),pointer :: rv_p_oag(:) =>NULL()

ckino Dec.8.2011:              type(s_rv1),pointer :: sv_p_oorhat(:) =>NULL()

ckino Dec.8.2011:             real(8),pointer :: rv_p_omad(:) =>NULL()

ckino Dec.8.2011:            complex(8),pointer :: zv_p_osmrho(:) =>NULL()

ckino Dec.8.2011:           complex(8),pointer :: zv_p_osmpot(:) =>NULL()

ckino Dec.8.2011:          integer,pointer :: iv_p_oips0(:) =>NULL()

ckino Dec.8.2011:         real(8),pointer :: rv_p_oqlv(:) =>NULL()

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision awald,alat,vol,plat(9),gmax,xx
C ... Heap
       integer ::iwdummy 

      call info0(50,1,0,' supot : potential setup ... ')


      nbas=sctrl%nbas
      nbasp=sctrl%nbasp
      nsp=sctrl%nspin


      alat=slat%alat
      vol=slat%vol
ckino Dec.8.2011:         rv_p_opos => slat%rv_p_opos 



      awald=slat%awald
      nkd=slat%nkd
      nkq=slat%nkq
ckino Dec.8.2011:         rv_p_odlv => slat%rv_p_odlv 

ckino Dec.8.2011:         rv_p_oqlv => slat%rv_p_oqlv 



      nbaspp = 2*nbasp - nbas

C --- Monopole Madelung matrix ---
      if (mod(mode,10) .eq. 1) then
ckino Dec.8.2011:           allocate(rv_p_omad(nbaspp**2))
         allocate(spot%rv_p_omad(nbaspp**2))
ckino Dec.8.2011:           rv_p_omad(:)=0.0d0
         spot%rv_p_omad(:)=0.0d0

ckino Dec.8.2011:              call madmat ( nbaspp , rv_p_opos , awald , alat , vol , rv_p_odlv 
ckino Dec.8.2011:              call madmat ( nbaspp , slat%rv_p_opos , awald , alat , vol , rv_p_odlv 
ckino 12月.27.2011:              call madmat ( nbaspp , slat%rv_p_opos , awald , alat , vol , slat%rv_p_odlv 
            call madmat ( nbaspp , slat%rv_p_opos , awald , alat , vol , slat%rv_a_odlv 
ckino Dec.8.2011:       .      , nkd , rv_p_oqlv , nkq , rv_p_omad ) 
ckino Dec.8.2011:       .      , nkd , slat%rv_p_oqlv , nkq , rv_p_omad ) 
     .      , nkd , slat%rv_p_oqlv , nkq , spot%rv_p_omad ) 





ckino Dec.8.2011:           spot%rv_p_omad => rv_p_omad 


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
ckino Dec.8.2011:           allocate(rv_p_ogv(abs(ngmx*3)))
ckino 12月.27.2011:           allocate(slat%rv_p_ogv(abs(ngmx*3)))
         allocate(slat%rv_a_ogv(abs(ngmx*3)))
c         if (ngmx*3<0) rv_p_ogv(:)=0.0d0

ckino Dec.8.2011:           allocate(iv_p_okv(abs(ngmx*3)))
         allocate(slat%iv_p_okv(abs(ngmx*3)))
c         if (ngmx*3<0) iv_p_okv(:)=0

           call gvlist ( alat , plat , iwdummy , n1 , n2 , n3 , gmax , 8 
ckino Dec.8.2011:       .     , ngmx , ng , iv_p_okv , rv_p_ogv , xx , xx ) 
ckino Dec.8.2011:       .     , ngmx , ng , iv_p_okv , slat%rv_p_ogv , xx , xx ) 
ckino 12月.27.2011:       .     , ngmx , ng , slat%iv_p_okv , slat%rv_p_ogv , xx , xx ) 
     .     , ngmx , ng , slat%iv_p_okv , slat%rv_a_ogv , xx , xx ) 



        if (ng .ne. ngmx) call rx('supot: bug in gvlist')


        slat%ng=ng
ckino Dec.8.2011:           slat%rv_p_ogv => rv_p_ogv 

ckino Dec.8.2011:           slat%iv_p_okv => iv_p_okv 

        slat%gmax=gmax
        i_copy_size=size(slat%nabc)
        call ii8copy(i_copy_size,ngabc,1,slat%nabc,1)

C       call prmx('gv',w(ogv),ngmx,ng,3)

ckino Dec.8.2011:           allocate(zv_p_osmrho(k1*k2*k3*nsp))
         allocate(spot%zv_p_osmrho(k1*k2*k3*nsp))
c         if (-k1*k2*k3*nsp<0) zv_p_osmrho(:)=0.0d0
ckino Dec.8.2011:           zv_p_osmrho(:)=0d0
         spot%zv_p_osmrho(:)=0d0

ckino Dec.8.2011:           allocate(zv_p_osmpot(k1*k2*k3*nsp))
         allocate(spot%zv_p_osmpot(k1*k2*k3*nsp))
ckino Dec.8.2011:           zv_p_osmpot(:)=0.0d0
         spot%zv_p_osmpot(:)=0.0d0

ckino Dec.8.2011:           allocate(sv_p_oorhat(3*nbaspp))
         allocate(spot%sv_p_oorhat(3*nbaspp))
cki         sv_p_oorhat(:)=0
         do i=1,3*nbaspp
ckino Dec.8.2011:              nullify(sv_p_oorhat(i)%v)
            nullify(spot%sv_p_oorhat(i)%v)
         enddo


ckino Dec.8.2011:           spot%zv_p_osmrho => zv_p_osmrho 

ckino Dec.8.2011:           spot%zv_p_osmpot => zv_p_osmpot 

ckino Dec.8.2011:           spot%sv_p_oorhat => sv_p_oorhat 



C   ... Setup for symmetrization
ckino Dec.8.2011:           allocate(iv_p_oips0(ng))
         allocate(slat%iv_p_oips0(ng))
ckino Dec.8.2011:           iv_p_oips0(:)=0.0d0
         slat%iv_p_oips0(:)=0.0d0

ckino Dec.8.2011:           allocate(zv_p_obgv(ng))
ckino 12月.27.2011:           allocate(slat%zv_p_obgv(ng))
         allocate(slat%zv_a_obgv(ng))
ckino Dec.8.2011:           zv_p_obgv(:)=0.0d0
ckino 12月.27.2011:           slat%zv_p_obgv(:)=0.0d0
         slat%zv_a_obgv(:)=0.0d0


        ngrp=slat%nsgrp
ckino Dec.13.2011:           rv_p_og => slat%rv_p_osymgr 

ckino Dec.8.2011:           rv_p_oag => slat%rv_p_oag 


ckino Dec.8.2011:               call sgvsym ( ngrp , rv_p_og , rv_p_oag , ng , rv_p_ogv , iv_p_oips0 
ckino Dec.8.2011:               call sgvsym ( ngrp , rv_p_og , slat%rv_p_oag , ng , rv_p_ogv , iv_p_oips0 
ckino Dec.8.2011:               call sgvsym ( ngrp , rv_p_og , slat%rv_p_oag , ng , slat%rv_p_ogv , iv_p_oips0 
ckino Dec.13.2011:               call sgvsym ( ngrp , rv_p_og , slat%rv_p_oag , ng , slat%rv_p_ogv , slat%iv_p_oips0 
ckino 12月.27.2011:               call sgvsym ( ngrp , slat%rv_p_osymgr , slat%rv_p_oag , ng , slat%rv_p_ogv , slat%iv_p_oips0 
ckino 12月.27.2011:               call sgvsym ( ngrp , slat%rv_p_osymgr , slat%rv_a_oag , ng , slat%rv_p_ogv , slat%iv_p_oips0 
             call sgvsym ( ngrp , slat%rv_p_osymgr , slat%rv_a_oag , ng , slat%rv_a_ogv , slat%iv_p_oips0 
ckino Dec.8.2011:       .       , zv_p_obgv ) 
ckino 12月.27.2011:       .       , slat%zv_p_obgv ) 
     .       , slat%zv_a_obgv ) 






ckino Dec.8.2011:           slat%iv_p_oips0 => iv_p_oips0 

ckino Dec.8.2011:           slat%zv_p_obgv => zv_p_obgv 



      endif
#endif
      end subroutine supot


