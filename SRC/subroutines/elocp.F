      module m_elocp
      use m_lmfinit,only: nspec,n0
      real(8),allocatable,protected:: ehl(:,:),rsml(:,:)
      
      contains
      subroutine elocp() 
      use m_struc_def  
      use m_lmfinit,only: nkaph,stdo,nspec,nbas,nsp,ssite=>v_ssite,sspec=>v_sspec,n0,ex=>ehl,rx=>rsml
      use m_uspecb,only: uspecb,nkapii
C- Make envlope parameters for extended local orbitals
C ----------------------------------------------------------------------
Ci Inputs
Ci         : 0 do nothing ; just return
Ci         : 1 make core and augmentation matrices
Co Outputs
Co   ehl,rsml
Co         : smoothing radius and energy set for extended local orbitals
Cr Remarks
Cu Updates
Cu   06 Jul 05 first created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
c      integer n0
c      parameter (n0=10)hl
cki      type(s_site)::ssite(*)
c      type(s_site)::ssite(nbas)
c      type(s_spec)::sspec(*)

C ... Local parameters
      character spid*8
      logical eloc
      integer ib,ibs,iclbsj,ipr,iprint,is,k,l,lmxa,lmxb,nglob,nkap0,
     .nr,nrmx,nrspec
      parameter (nrmx=1501, nkap0=3)
      integer:: lh(nkap0),nkape,idamax,nkapex
      integer ,allocatable :: ips_iv(:)
      double precision z,a,rmt,rs3,eh3,vmtz,xx
      double precision rofi(nrmx),vseli(4,n0),vsel(4,n0,nbas),
     .pnu(n0,2),pnz(n0,2),eh(n0,nkap0),rsmh(n0,nkap0),
     .     pnui(n0,2),pnzi(n0,2)
      real(8):: ehls(n0*2),rsmls(n0*2)
C ... Heap
      integer ::iwdummy,i_copy_size,i_spackv
C --- Setup ---
      call tcn('elocp')
      if(allocated(ehl)) deallocate(ehl,rsml)
      allocate( ehl(n0,nspec),rsml(n0,nspec))
!!      
      if(allocated(ex)) deallocate(ex,rx)
      allocate( ex(n0,nspec),rx(n0,nspec))
      
      ipr = iprint()
      eloc = .false.
      vsel=0d0
C --- Find val, slo, K.E. for all sites ---
      do  ib = 1, nbas
        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pnu)
        call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
        i_copy_size=size(ssite(ib)%pz)
        call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)
        spid = sspec(is)%name
        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt
        z=sspec(is)%z
        lmxa=sspec(is)%lmxa
        lmxb=sspec(is)%lmxb
        if (lmxa .eq. -1) goto 10
        if (pnz(idamax(lmxb+1,pnz,1),1) .lt. 10) goto 10
        eloc = .true.
        call radmsh(rmt,a,nr,rofi)
        call loctsh ( 1101 , spid , z , a , nr , nr , nsp , lmxa , rofi
     .       , ssite(ib)%rv_a_ov0 , pnu , pnz , xx , xx , vmtz , vsel ( 1 , 1 , ib ) 
     .       , rsml , ehl )
   10   continue
      enddo
      if (.not. eloc) goto 999
      if (ipr .ge. 30) write(stdo,199)
  199 format(/' elocp:')
C --- Determine shape of smooth Hankel tails for local orbitals ---
      allocate(ips_iv(nbas))
      i_copy_size=1;
      do i_spackv=1,nbas
        ips_iv(i_spackv)=ssite(i_spackv)%spec
      enddo
C ... Loop over species containing extended local orbitals
      do  is = 1, nspec
c        do i_spacks=is,is
c          call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
c        enddo
        spid=sspec(is)%name
        z=sspec(is)%z
        lmxa=sspec(is)%lmxa
        lmxb=sspec(is)%lmxb
        if (lmxa .eq. -1) goto 20
        nrspec = iabs ( iclbsj ( is , ips_iv , - nbas , nbas ) )
        if (nrspec .eq. 0) goto 20
        ib = iclbsj ( is , ips_iv , nbas , 1 )
        i_copy_size=size(ssite(ib)%pnu)
        call dcopy(i_copy_size,ssite(ib)%pnu,1,pnui,1)
        i_copy_size=size(ssite(ib)%pz)
        call dcopy(i_copy_size,ssite(ib)%pz,1,pnzi,1)
        if (pnzi(idamax(lmxb+1,pnzi,1),1) .lt. 10) goto 20

C   ... Average over sites within this species
        vseli=0d0 !call dpzero(vseli,4*n0)
        do  ibs = 1, nrspec
          ib = iclbsj ( is , ips_iv , nbas , ibs )
          i_copy_size=size(ssite(ib)%pz)
          call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)
          if (pnzi(idamax(lmxb+1,pnzi,1),1) .lt. 10) goto 22
c          call dpadd(vseli,vsel(1,1,ib),1,4*n0,1/dble(nrspec))
          vseli = vseli + vsel(:,:,ib)/dble(nrspec)
   22     continue
        enddo

C   ... Printout of input for parameters
        if (ipr .ge. 90/1) then
          write(stdo,261) spid
  261     format(/'  l  site    Eval        Val         Slo         K.E.',
     .    5x,'species ',a)
          do  l = 0, lmxb
            if (pnz(l+1,1) .lt. 10) goto 24
            do  ibs = 1, nrspec
              ib = iclbsj ( is , ips_iv , nbas , ibs )

              write (stdo,260)
     .        l,ib,vsel(4,l+1,ib),(vsel(k,l+1,ib),k=2,4)
  260         format(i3,i4,4f12.6:a)
  262         format(i3,' avg',4f12.6)
            enddo
            write (stdo,262) l,vseli(4,l+1),(vseli(k,l+1),k=2,4)

   24       continue
          enddo
        endif

C   ... Make parameters for this species
        call dpzero(rsmh,n0*nkap0)
        call dpzero(eh,n0*nkap0)
        
        call uspecb(0,2,sspec,is,is,lh,rsmh,eh,nkapex)
        nkape = nkapii(is)
c       call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkape)
        
        rs3=sspec(is)%rs3
        eh3=sspec(is)%eh3
        vmtz=sspec(is)%vmtz
        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt
        rsmls=0d0 !call dpzero(rsmls,n0*2)
        ehls=0d0  !call dpzero(ehls,n0*2)
        call radmsh(rmt,a,nr,rofi)
        call loctsh ( 1102 , spid , xx , a , nr , nr , nsp , lmxa , rofi !1102 means spin-averaged ehl,rsml
     .       , iwdummy , pnui , pnzi , rs3 , eh3 , vmtz , vseli , rsmls , ehls )
!!        
        rsml(1:n0,is) = rsmls(1:n0)
        ehl(1:n0,is) = ehls(1:n0)
!! copy xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        rx(1:n0,is) = rsml(1:n0,is)
        ex(1:n0,is) = ehls(1:n0)
c!!
        call dcopy(n0,rsml,1,rsmh(1,nkaph),1)
        call dcopy(n0,ehl, 1,eh(1,nkaph),1)
        call uspecb(1,1,sspec,is,is,lh,rsmh,eh,nkaph)

   20   continue
      enddo

  999 continue
      call tcx('elocp')
      if (allocated(ips_iv)) deallocate(ips_iv)
      end subroutine elocp
      end module m_elocp
