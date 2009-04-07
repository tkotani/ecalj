      subroutine iax2fd(mode,nbas,ntab,iax,ntabib,iaxib,mxcsiz)
C- Converts an iax table in standard format from/to a fixed basis format
C  and purges table entries for which iax(1) or iax(2) are zero
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 copy ntab,iax to ntabib,iaxib
Ci         :1 copy ntabib,iaxib to ntab,iax
Ci   nbas  :size of basis
Cio Inputs/Outputs
Cio  ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Cio        :Input for mode=0, output for mode=1
Cio  iax   :neighbor table containing pair information (pairc.f)
Cio        :Input for mode=0, output for mode=1
Cio  ntabib:cluster size connected to site ib
Cio        :Input for mode=1, output for mode=0
Cio  iaxib :contains samie info as table iax, in fixed-basis format
Cio        :Input for mode=1, output for mode=0
Cio  mxcsiz:Input for mode=0: Third dimension of iaxib. Specifies
Cio        :the maximum allowed cluster size iaxib (max dim. of iaxib)
Cio        :Output for mode=1: number pairs in all clusters
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   29 Aug 02  First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer niax,mode,nbas,mxcsiz,ntabib(nbas),ntab(nbas+1)
      parameter (niax=10)
      integer iax(niax,1),iaxib(niax,nbas,1)
C ... Local parameters
      integer ib,is0,is,shrink

      if (mode .eq. 1) ntab(1) = 0

      do  ib = 1, nbas
C   ... Copy from iax to iaxib
        if (mode .eq. 0) then
          ntabib(ib) = ntab(ib+1)-ntab(ib)
          if (ntabib(ib) .gt. mxcsiz) call rxi(
     .      'iax2fd: iaxib dimensioned too small, need mxcsiz=',
     .      ntabib(ib))
          is0 = ntab(ib)
          shrink = 0
          do  is = ntab(ib)+1, ntab(ib+1)
            if (iax(1,is) .eq. 0 .or. iax(2,is) .eq. 0) then
              shrink = shrink+1
              is0 = is0+1
            else
              call icopy(niax,iax(1,is),1,iaxib(1,ib,is-is0),1)
            endif
          enddo
          ntabib(ib) = ntabib(ib)-shrink
C   ... Copy from iaxib to iax
        else
          is0 = ntab(ib)
          ntab(ib+1) = ntab(ib) + ntabib(ib)
          shrink = 0
          do  is = is0+1, is0+ntabib(ib)
            if (iaxib(1,ib,is-is0) .eq. 0 .or.
     .          iaxib(2,ib,is-is0) .eq. 0) then
              shrink = shrink+1
            else
              call icopy(niax,iaxib(1,ib,is-is0),1,iax(1,is-shrink),1)
            endif
          enddo
          ntab(ib+1) = ntab(ib+1)-shrink
        endif
      enddo

      if (mode .eq. 1) mxcsiz = ntab(nbas+1)

      end

