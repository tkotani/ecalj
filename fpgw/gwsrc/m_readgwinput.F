      module m_readgwinput
      implicit none
!! Read values from GWinput
      real(8),protected:: egauss,ecut,ecuts,ebmx
      integer,protected:: nbcut,nbcut2,nbmx
      integer,protected:: mtet(3),nmbas
      integer,protected,allocatable:: imbas(:)

      contains
      subroutine ReadGWinputKeys()
      use m_keyvalue,only: Getkeyvalue
      use m_genallcf_v3,only: natom
      use m_readhbe,only: nband
c     use NaNum,only: NaN
c      integer:: nband
      integer:: istat
      integer,allocatable:: imbasd(:)
      call Getkeyvalue("GWinput","ecut_p" ,ecut, default=1d10 )
      call Getkeyvalue("GWinput","ecuts_p",ecuts,default=1d10 )
      
      call getkeyvalue("GWinput","nbcutlow",  nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2
      
      call getkeyvalue("GWinput","nband_chi0",nbmx, default=nband)
      call getkeyvalue("GWinput","emax_chi0", ebmx, default=1d10  )
      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/)) 
      call getkeyvalue("GWinput","GaussianFilterX0", egauss, default=0d0 )
      
!! this is effective only for chipm mode
      nmbas=natom
      allocate(imbas(nmbas),imbasd(nmbas))
      imbasd = -9999
      istat=-999               ! istat=-9999 means noumber of readin arguments is returened in istat.
      call getkeyvalue("GWinput","MagAtom", imbas,nmbas,status=istat,default=imbasd)
      nmbas = istat
      if(nmbas>0) write(6,"('Readin MagAtom nmbas =',i3,' imbas= ',10i3)") nmbas,imbas(1:nmbas)
      end subroutine
      end module

