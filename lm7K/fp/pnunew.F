      subroutine pnunew(nbas,nsp,ssite,sspec,pmin,pmax,lfrzw,hab,sab,
     .qbyl,hbyl,eferm)

      use m_struc_def  !Cgetarg
!!= Makes new boundary conditions pnu for phi,phidot =
!!* takao sep2010. P is setted to that at efermi if PZ is for semicore.


C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pnu pz ov0
Ci     Stored:    pnu pz
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa rmt idmod mxcst z a nr
Ci     Stored:    *
Ci     Passed to: *
Ci   lfrzw :0, float pnu to band CG, provided IDMOD=0 for that channel
Ci         :1 freeze all pnu for all species.
Ci         :  NB: pnu are also frozen for specific species
Ci         :      that have nonzero 4's bit of species->mxcst.
Ci   pmin  :lower bound for fractional part of P
Ci   pmax  :upper bound for fractional part of P
Ci   hab   :<u,s | H | u,s> for each pair uu, us, su, ss; see Remarks
Ci   sab   :<u,s | 1 | u,s>
Ci         :NB: hab and sab are no longer used.
Ci   qbyl  :l-decomposed charge
Ci   hbyl  :l-decomposed eigenvalue sum
Cl Local variables
Cl   lfrzv  if T, freeze valence pnu
Cl   lfrzz  if T, freeze local orbital pnz
Co Outputs
Ci   ssite->pnu :are floated to their band CG
Cr Remarks
Cu Updates
Cu   28 Jun 06 Handles case idmod=3; New constraints pmin,pmax
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu    9 May 02 Added species-specific freezing of pnu
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
Cu   20 Sep 01 Some patches to work in pathological cases
Cu   17 Sep 01 When local orbital present, allow semiore pnz to float
Cu   28 Aug 01 Extended to local orbitals.  For now, freeze pnu
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,nsp,lfrzw,n0,nab
      parameter(n0=10,nab=9)
      real(8):: pmin(n0) , pmax(n0) , qbyl(n0,nsp,nbas) , hbyl(n0,nsp,nbas) 
     ., sab(nab,n0,nsp,nbas) , hab(nab,n0,nsp,nbas)
cki      type(s_site)::ssite(*)
      type(s_site)::ssite(nbas)
      type(s_spec)::sspec(*)

C ... Local parameters
      logical lpz,lfrzv,lfrzz
       integer:: idmod(n0) , stdo , lgunit , ipr , ib , is , lmxa , 
     . l , ipqn , m , isp , nr , nn , mxcst 
       real(8),pointer :: rv_p_ov0(:) =>NULL()

      real(8) ,allocatable :: g_rv(:)
      real(8) ,allocatable :: gp_rv(:)
      real(8) ,allocatable :: rofi_rv(:)
      real(8) ,allocatable :: v0i_rv(:)

      double precision pi,rmt,p1,ebar,a,d0l,pfree,pold,ptry,z,
     .val(5),slo(5),dl,phi,dphi,pnu(n0,2),pnz(n0,2)
      double precision ez,umegam,phip,dphip,dlphi,dlphip,cz
      double precision pznew,fi(0:10),gi(0:10),xx,dnz
      character spid*8

      integer ::iwdummy ,i_copy_size,nnz,nnv
      real(8):: eferm,eee
      logical:: lsemicorepz

      call info(30,1,0,' Make new boundary conditions for phi,phidot..',
     .0,0)
      stdo = lgunit(1)
      call getpr(ipr)
      pi = 4d0*datan(1d0)

C --- For each site, do ---
      do  ib = 1, nbas

        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pnu)
        call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
        i_copy_size=size(ssite(ib)%pz)
        call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)
         rv_p_ov0 => ssite(ib)%rv_p_ov0 



        lmxa=sspec(is)%lmxa
        rmt=sspec(is)%rmt
        i_copy_size=size(sspec(is)%idmod)
        call i8icopy(i_copy_size,sspec(is)%idmod,1,idmod,1)
        mxcst=sspec(is)%mxcst

        if (lmxa .eq. -1) goto 10
c        do i_spacks=is,is
c          call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
c        enddo
        call r8tos8(sspec(is)%name,spid)

        if (mod(mxcst/4,2) .ne. 0) call ivset(idmod,1,n0,1)

        if (ipr .ge. 20) write(stdo,320) ib,is,spid
  320   format(/' site',i5,'   species',i4,':',a)
        if (ipr .ge. 20) write(stdo,311)
        do  l = 0, lmxa
          do  isp = 1, nsp
            m = l+1
            p1 = 2d10
            pznew = 2d10
C         Initially set lfrzv,lfrzz to external constraints.
            lfrzv = (mod(idmod(m),10).ne.0 .and. mod(idmod(m),10).ne.3)
     .      .or. lfrzw .ne. 0
            lfrzz = lfrzv
            lpz = pnz(m,1) .ne. 0

            lsemicorepz = .false.
            if(lpz)then
              nnv = int(mod(pnu(m,1),10d0))-l-1
              nnz = int(mod(pnz(m,1),10d0))-l-1
              if(nnz<nnv) then !this means semicore for local orbital.
c                write(*,"(a,2i3,l,5f12.4)") 'zzzzz ',m,isp,lpz,pnu(m,isp),pnz(m,isp)
                lsemicorepz = .true.
              endif
            endif  

C         if (dabs(qbyl(m,isp,ib)) .gt. 1d-8 .and. .not. lpz) then
            if (dabs(qbyl(m,isp,ib)) .gt. 1d-8) then
              ebar = hbyl(m,isp,ib)/qbyl(m,isp,ib)

C       ... Log derivative by direct num integration
              is=ssite(ib)%spec
              z=sspec(is)%z
              a=sspec(is)%a
              nr=sspec(is)%nr
              allocate(g_rv(nr*2))
              allocate(rofi_rv(nr*2))
              call radmsh ( rmt , a , nr , rofi_rv )
              allocate(v0i_rv(nr))
               call dpscop ( rv_p_ov0 , v0i_rv , nr , 1 + nr * ( isp - 1 ) , 
     .         1 , 1d0 ) 

              if (mod(idmod(m),10) .eq. 3) then
                val(1) = rmt
                dl = dtan(pi*(0.5d0-mod(pnu(m,isp),10d0)))
                slo(1) = dl + 1d0
                nn = int(mod(pnu(m,1),10d0))-l-1
                allocate(gp_rv(8*nr))

C             xx = dval(w(ov0i),nr)
                call phidx ( 0 , z , l , v0i_rv , 0d0 , 0d0 , rofi_rv , 
     .          nr , 2 , 1d - 12 , ebar , val , slo , nn , g_rv , gp_rv 
     .          , phi , dphi , phip , dphip , xx , xx , xx , xx , xx )

C         ... cz = estimate for energy of orbital with b.c. connecting
C             to Hankel of energy 0
                dlphi = rmt*dphi/phi
                dlphip = rmt*dphip/phip
                umegam = -(phi/phip)*(-l-1-dlphi)/(-l-1-dlphip)
                cz = ebar + umegam
                ebar = cz
C         ... estimate for log der. of wf for constant pot of value C
C             dh_l/dr = l*h_l/r - h_l+1, h=g_l/r**(l+1)
C             when cz -> 0,  dl -> -l-1
C             call bessl(cz*rmt**2,m,fi,gi)
C             dl = (l*gi(l) - gi(l+1))/gi(l)
C             p1 = 0.5d0 - datan(dl)/pi
C             val(1) = rmt
C             slo(1) = dl + 1d0
              endif

!! takao Sep24 2010
              if(lsemicorepz)then
                print *,' pnunew: valence with semicore ebar=efermi'
                ebar=eferm
              endif

              call phidx ( 2 , z , l , v0i_rv , 0d0 , 0d0 , rofi_rv , nr , 
     .         0 , 1d - 12 , ebar , val , slo , nn , g_rv , iwdummy , phi , 
     .         dphi , iwdummy , iwdummy , iwdummy , iwdummy , iwdummy , iwdummy 
     .         , iwdummy ) 
              nnv=nn

C           dphip = (slo(2)-phip)/rmt
              if (nn .eq. int(pnu(m,1))-l-1) then
                dl = rmt*slo(1)/val(1) - 1
                p1 = 0.5d0 - datan(dl)/pi
c              elseif (ipr .ge. 10) then
              elseif (ipr .ge. 10.and.(.not.lsemicorepz)) then
                print *,' node =',nn
                call info2(10,0,0,' (warning) failed to find proper '
     .          //'node count for l=%i  ebar=%;4d: pnu not calc',l,
     .          ebar)
              endif

              if(lsemicorepz) then
                dl = rmt*slo(1)/val(1) - 1
                p1 = 0.5d0 - datan(dl)/pi
c                print *,'zzzzzz p1=',p1,dl
              endif  

C       ... Estimate new pnz for semicore state
              if (lpz .and. int(mod(pnz(m,1),10d0)).lt.int(pnu(m,1))) then
                val(1) = rmt
                dnz = dtan(pi*(0.5d0-mod(pnz(m,isp),10d0)))
                slo(1) = dnz + 1d0
                nn = int(mod(pnz(m,1),10d0))-l-1
                allocate(gp_rv(8*nr))

                call phidx ( 0 , z , l , v0i_rv , 0d0 , 0d0 , rofi_rv , 
     .          nr , 2 , 1d - 12 , ez , val , slo , nn , g_rv , gp_rv , 
     .          phi , dphi , phip , dphip , xx , xx , xx , xx , xx )

C             dphip = (slo(2)-phip)/rmt
                dlphi = rmt*dphi/phi
                dlphip = rmt*dphip/phip
                umegam = -(phi/phip)*(-l-1-dlphi)/(-l-1-dlphip)
C         ... cz = estimate for energy of orbital with b.c. connecting
C             to Hankel of energy 0 (=> solution for interstitial
C             is constant potential, value C
                cz = ez + umegam
C         ... estimate for log der. of wf for constant pot of value C
C             Maybe better to recalc. val,slo for this cz like above?
                call bessl(cz*rmt**2,m,fi,gi)
C             dh_l/dr = l*h_l/r - h_l+1, h=g_l/r**(l+1)
C             when cz -> 0,  dl -> -l-1
                dl = (l*gi(l) - gi(l+1))/gi(l)
                pznew = 0.5d0 - datan(dl)/pi
c                lfrzv = .true. !takao commentout Sep2010
              else
                lfrzz = .true.
                pznew = 0.5d0 - datan(dble(l))/pi
                ez = 0
              endif
              if (allocated(gp_rv)) deallocate(gp_rv)
              if (allocated(v0i_rv)) deallocate(v0i_rv)
              if (allocated(rofi_rv)) deallocate(rofi_rv)
              if (allocated(g_rv)) deallocate(g_rv)

C       ... First root of quad equation is the one we want
C           huu = hab(1,m,isp,ib)
C           hus = (hab(2,m,isp,ib) + hab(3,m,isp,ib))/2d0
C           hss = hab(4,m,isp,ib)
C           suu = sab(1,m,isp,ib)
C           sus = (sab(2,m,isp,ib) + sab(3,m,isp,ib))/2d0
C           sss = sab(4,m,isp,ib)
C           a = hss - ebar*sss
C           b = 2d0*(hus - ebar*sus)
C           c = huu - ebar*suu
C           ddd = b*b - 4*a*c
C           if (ddd .ge. 0) then
C             x1 = (-b-dsqrt(ddd))/(2*a)
CC            q1 = suu+2*x1*sus + x1*x1*sss
CC            h1 = huu+2*x1*hus + x1*x1*hss
C             p1 = 0.5d0-datan(rmt*x1)/pi
C           endif

            endif

C     ... Free-electron value for pnu
            ipqn = pnu(m,isp)
            d0l = l
            pfree = ipqn + 0.5d0 - datan(d0l)/pi

C     --- Set the new pnu ---
            pold = pnu(m,isp)
            ipqn = pold
            ptry = pold
c            print *,'zzzz: ', ipqn,  nnv + l+1, nnv
            if (dabs(p1) .lt. 1d10) ptry = nnv + l+1 +p1 !takao
c            if (dabs(p1) .lt. 1d10) ptry = ipqn+p1

c            print *,' lfrzv=', lfrzv
c            if (.not. lfrzv) then

            if (.not. lfrzv) then
              pnu(m,isp) = ptry
C       ... Permit pnu no lower than free electron value or pmin
              if (ptry .lt. pfree) pnu(m,isp) = pfree
              if (pmin(m) .gt. 0 .and. pmin(m) .lt. 1) then
                if (ptry .lt. ipqn+pmin(m)) pnu(m,isp) = ipqn+pmin(m)
              endif
C       ... Permit pnu no higher than pmax
              if (pmax(m) .gt. 0 .and. pmax(m) .lt. 1) then
                if (ptry .gt. ipqn+pmax(m)) pnu(m,isp) = ipqn+pmax(m)
              endif
            endif

            if (ipr .ge. 20 .and. isp .eq. 1) write(stdo,310)
     .      l,idmod(m),qbyl(m,isp,ib),ebar,pold,ptry,pfree,pnu(m,isp)
            if (ipr .ge. 20 .and. isp .eq. 2) write(stdo,410)
     .      idmod(m),qbyl(m,isp,ib),ebar,pold,ptry,pfree,pnu(m,isp)
  310       format(i2,i6,6f12.6,l)
  410       format(' spn 2',i2,6f12.6,l)
  311       format(' l  idmod     ql',9x,'ebar',7x,' pold',8x,
     .      'ptry',8x,'pfree',8x,'pnew',8x)

C     --- Set the new pnz ---
            if (lpz) then
              pold = mod(pnz(m,isp),10d0)
              ipqn = pold
              ptry = pold
              pfree = ipqn + 0.5d0 - datan(d0l)/pi
              if (dabs(pznew) .lt. 1d10) ptry = ipqn+pznew
              if (.not. lfrzz) then
                pnz(m,isp) = ptry + (pnz(m,isp)-mod(pnz(m,isp),10d0))
C         ... Permit pnu no lower than free electron value
                d0l = l
                if (ptry .lt. pfree)
     .          pnz(m,isp) = pfree + (pnz(m,isp)-mod(pnz(m,isp),10d0))
              endif

              if (ipr .ge. 20 .and. isp .eq. 1) write(stdo,520)
     .        l,idmod(m),ez,pold,ptry,pfree,pnz(m,isp)
              if (ipr .ge. 20 .and. isp .eq. 2) write(stdo,620)
     .        idmod(m),ez,pold,ptry,pfree,pnz(m,isp)
  520         format(i2,i6,'    sc      ',6f12.6)
  620         format(' spn 2',i2,'    sc      ',6f12.6)

            elseif (lpz) then
            endif

          enddo
        enddo


        i_copy_size=size(ssite(ib)%pnu)
        call dcopy(i_copy_size,pnu,1,ssite(ib)%pnu,1)
        i_copy_size=size(ssite(ib)%pz)
        call dcopy(i_copy_size,pnz,1,ssite(ib)%pz,1)

   10   continue
      enddo

      end subroutine pnunew


