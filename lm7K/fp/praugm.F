Cgetarg...info...           structure ['praugm', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['praugm', 'm_struc_def']
      subroutine praugm(sspec,is)

      use m_struc_def  !Cgetarg

C-  Print species information
C ----------------------------------------------------------------------
Ci Inputs
Ci   sspec :struct containing species-specific information
Ci   is    :species index (use 0 to print info for all species)
Co Outputs
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer is
Cgetarg       double precision sspec(1)
      type(s_spec)::sspec(*)

C ... Local parameters
      integer is1,is2,js,kmxt,lmxa,nglob,lgunit,stdo
      integer kmxv,lmxl,lfoca
      double precision rmt,rsma,rfoca,rg,rsmv
      character spid*8
C      parameter (n0=10)
C      integer idmod(n0),l
C      double precision pnu(n0*2)
C      character*1  s2*80, idcode(0:1)

      stdo = lgunit(1)
C      idcode(0) = ' '
C      idcode(1) = '*'
      is1 = is
      is2 = is
      if (is .le. 0) then
        is1 = 1
        is2 = nglob('nspec')
      endif

      write (stdo,501)
  501 format(/' species data:  augmentation',27x,'density'/
     .     ' spec       rmt   rsma lmxa kmxa',5x,
     .     ' lmxl     rg   rsmv  kmxv foca   rfoca')
      do  10  js = is1,is2
Cgetarg         call spacks(0,'spec name',sspec,spid,js,js)
        do i_spacks=js,js
          call spacks_copy('u',sspec(i_spacks)%name,js,js,spid,i_spacks)
        enddo

Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg...info...           real(8) :: rsma , real(8) :: sspec%rsma 1
Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: kmxt , integer(8) :: sspec%kmxt 1
Cgetarg         call upack('spec rmt rsma lmxa kmxt',sspec,js,rmt,rsma,lmxa,
Cgetarg      .    kmxt)

        rmt=sspec(js)%rmt
        rsma=sspec(js)%rsma
        lmxa=sspec(js)%lmxa
        kmxt=sspec(js)%kmxt

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           real(8) :: rg , real(8) :: sspec%rg 1
Cgetarg...info...           real(8) :: rsmv , real(8) :: sspec%rsmv 1
Cgetarg...info...           integer :: kmxv , integer(8) :: sspec%kmxv 1
Cgetarg         call upack('spec lmxl rg rsmv kmxv',sspec,js,lmxl,rg,rsmv,kmxv)

        lmxl=sspec(js)%lmxl
        rg=sspec(js)%rg
        rsmv=sspec(js)%rsmv
        kmxv=sspec(js)%kmxv

Cgetarg...info...           integer :: lfoca , integer(8) :: sspec%lfoca 1
Cgetarg...info...           real(8) :: rfoca , real(8) :: sspec%rfoca 1
Cgetarg         call upack('spec lfoca rfoca',sspec,js,lfoca,rfoca,0,0)

        lfoca=sspec(js)%lfoca
        rfoca=sspec(js)%rfoca


C        call upack('spec p idmod',sspec,js,pnu,idmod,0,0)
C        write (s2,102) (pnu(l+1),idcode(idmod(l+1)),l=0,lmxa)
C  102   format(10(f5.2,a1))
C        write (stdo,300) spid,rmt,rsma,lmxa,kmxt,s2(1:36)
C  300   format(1x,a,f6.3,f7.3,i6,i6,2x,a)
        write (stdo,500) spid,rmt,rsma,lmxa,kmxt,
     .                   lmxl,rg,rsmv,kmxv,lfoca,rfoca
  500   format(1x,a,f6.3,f7.3,2i5,6x,i4,2f7.3,i6,i5,f8.3)

   10 continue

Cgetarg       end
      end subroutine praugm



