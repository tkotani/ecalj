      subroutine mchan(lmdim,ssite,sspec,nsp,nsites,lsites,ib,ilm,io,
     .ichan,lchan)

      use m_struc_def  !Cgetarg

C- set or get Mulliken channel for site ib and ilm index
C ----------------------------------------------------------------------
Ci Inputs: ssite, sspec, nsp
Ci         lmdim,nsites,ib,ilm,io < 0 poke ichan into lchan
Ci                                > 0 get ichan from lchan
Ci                                = 0 print channel table
Co Outputs
Co   ichan :(io > 0) lchan(ilm,ib) for supplied (ilm,ib) pair
Co         :Otherwise, ichan is not set
Co   lchan :(io < 0) set lchan(ilm,ib) to ichan
Co         :Otherwise, lchan is not set
Cr Remarks
Cr    For the Mulliken decomposition it is convenient to keep a table
Cr    lchan(ilm,ib) which holds the DOS channel number associated with
Cr    site ib, and lm channel ilm. If the DOS is site projected then
Cr    ilm is 1; if l projected then ilm=1,lmax+1; if lm-projected then
Cr    ilm=1,(lmax+1)**2. The leading dimension lmdim hence depends on
Cr    the mode (see sumlst) in this way.
Cu Updates
Cu   20 Mar 01 Written by ATP
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer lmdim,nsp,nsites,lsites(nsites),lchan(lmdim,nsites),
     .ib,ilm,io,ichan
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

C Local Variables
      integer i,j,lgunit,js,igetss,jb
      character clabl*8

      if (io .lt. 0) then
        lchan(ilm,ib) = ichan
      elseif (io .gt. 0) then
        ichan = lchan(ilm,ib)
      else
        call awrit1('%N mchan: channel table, %i sites',
     .  ' ',256,lgunit(1),nsites)
        if (nsp .eq. 2)
     .  call awrit0(' (each channel splits into two: up, down spin)',
     .  ' ',256,lgunit(1))
        call awrit0(' site  label          channels',' ',128,lgunit(1))
        do  j = 1, nsites
          jb = lsites(j)
          js = int(ssite(jb)%spec)

c          do i_spacks=js,js
c            call spacks_copy('u',sspec(i_spacks)%name,js,js,clabl,i_spacks)
c          enddo
          clabl=sspec(js)%name
          write (lgunit(1),1) j, clabl, (lchan(i,j),i=1,lmdim)
        enddo
      endif
    1 format (1x,i3,5x,a8,1x,256i3)
      end subroutine mchan


