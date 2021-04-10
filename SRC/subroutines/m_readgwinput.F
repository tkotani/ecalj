      module m_readgwinput
      implicit none
!! Read values from GWinput
      real(8),protected:: egauss,ecut,ecuts,ebmx ,ebmx_sig,ua_
      integer,protected:: nbcut,nbcut2,nbmx,nbcutlow_sig,nbmx_sig,iSigmode
      integer,protected:: mtet(3),nmbas
      integer,protected,allocatable:: imbas(:)
      logical,protected:: keeppositivecou
      contains
      subroutine SetIsigmode(ism)
      integer ism
      isigmode=ism
      end subroutine
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
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=999999 )
      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2
      
      call getkeyvalue("GWinput","nband_chi0",nbmx, default=nband)
      call getkeyvalue("GWinput","emax_chi0", ebmx, default=1d10  )
      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/)) 
      call getkeyvalue("GWinput","GaussianFilterX0", egauss, default=0d0 )
      call getkeyvalue("GWinput","nbcutlow_sig",nbcutlow_sig, default=0 )

      call getkeyvalue("GWinput","iSigMode"  ,iSigMode, default=3 )
      call getkeyvalue("GWinput","nband_sigm",nbmx_sig, default=9999999)
      call getkeyvalue("GWinput","emax_sigm", ebmx_sig, default=1d10)
      
c      call GETKEYVALUE("GWinput","EXonly",wex,default=0d0)
      call getkeyvalue("GWinput","gauss_img",ua_,default=1d0)

!     !! If keeppositevecou=F, it can cause negative eigenvalue for Coulomb matrix ===> this probably means
      !! problem in strxq
      call getkeyvalue("GWinput","KeepPositiveCou",KeepPositiveCou,default=.true.)
      
!! this is effective only for chipm mode
      nmbas=natom
      allocate(imbas(nmbas),imbasd(nmbas))
      imbasd = -9999
      istat=-9999   ! istat=-9999 means number of readin arguments is returened in istat.
      call getkeyvalue("GWinput","MagAtom", imbas,nmbas,status=istat,default=imbasd)
      nmbas = istat
      if(nmbas>0) write(6,"('Readin MagAtom nmbas =',i3,' imbas= ',10i3)") nmbas,imbas(1:nmbas)
      end subroutine
      end module

