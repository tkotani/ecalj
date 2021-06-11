      subroutine tbhsi(sspec,nspec,nermx,net,et,ipet,nrt,rt,iprt,ltop)
      use m_uspecb,only:uspecb
      use m_struc_def,only:s_spec
C- Table of orbital envelope energies and smoothing radii
C ----------------------------------------------------------------------
Ci Inputs
Ci   sspec :struct containing species-specific information
Ci     Passed to: uspecb
Ci   nspec  :number of species
Ci   nermx :maximum allowed size of et,rt
Co Outputs
Co   net   :size of table et
Co   et    :table of all inequivalent energies
Co   ipet  :index to which entry in et a given orbital belongs
Co   nrt   :size of table rt
Co   rt    :table of all inequivalent smoothing radii
Co   iprt  :index to which entry in rt a given orbital belongs
Co   ltop  :largest l at any site
Cr Remarks
Cr   An orbital (l,ik,is) has smoothing radius index ir=iprt(l,ik,is)
Cr   and energy index ie=ipet(l,ik,is).  Its smoothing radius and
Cr   energy are rt(ir) and et(ie), respectively.
Cu Updates
Cu   12 Aug 04 First implementation of extended local orbitals
Cu   10 Apr 02 Redimensionsed ipet,iprt to accomodate larger lmax
Cu   18 May 00 Adapted from nfp tb_hsi.f
C ----------------------------------------------------------------------
c      implicit none
C     ... Passed parameters
      type(s_spec):: sspec(*) 
      integer nspec,nermx,net,nrt,ltop,n0,nkap0
      parameter (n0=10,nkap0=3)
      integer ipet(n0,nkap0,nspec),iprt(n0,nkap0,nspec)
      double precision et(nermx),rt(nermx)
C ... Local parameters
      integer is,j,k
      integer lh(nkap0),nkape,ik,l
      double precision x1,x2,rsmh(n0,nkap0),eh(n0,nkap0)
      logical dcmpre
      dcmpre(x1,x2) = dabs(x1-x2) .lt. 1d-8

      net = 0
      nrt = 0
      ltop = -1

C --- Loop over orbitals (is,io) ---
      do  is = 1, nspec
        call uspecb(is,lh,rsmh,eh,nkape)

        do  ik = 1, nkape
          ltop = max0(ltop,lh(ik))
C       l2 = -1
          do  l = 0, lh(ik)
            if (rsmh(l+1,ik) .gt. 0) then
C         if (l .gt. l2) then
C         call gtbsl(l,lh(ik),eh(1,ik),rsmh(1,ik),l2)

C     ... Find eh, or add it to list
              j = 0
              do  k = 1, net
                if (dcmpre(eh(l+1,ik),et(k))) then
                  j = k
                  goto 31
                endif
              enddo
  31          continue
              if (j .gt. 0) then
                ipet(l+1,ik,is) = j
              else
                net = net+1
                if (net .gt. nermx) call rx('tbhsi: nermx exceeded for et')
                et(net) = eh(l+1,ik)
                ipet(l+1,ik,is) = net
              endif

C     ... Find rmsh, or add it to list
              j = 0
              do  k = 1, nrt
                if (dcmpre(rsmh(l+1,ik),rt(k))) then
                  j = k
                  goto 32
                endif
              enddo
  32          continue
              if (j .gt. 0) then
                iprt(l+1,ik,is) = j
              else
                nrt = nrt+1
                if (nrt .gt. nermx) call rx('tbhsi: nermx exceeded for rt')
                rt(nrt) = rsmh(l+1,ik)
                iprt(l+1,ik,is) = nrt
              endif
            endif
          enddo
        enddo
      enddo

C ... Debugging printout
C      write (*,400) net,nrt,ltop
C  400 format(/' tbhsi: ',i4,' energies,',i4,'  radii,  ltop=',i3)
C      write (*,401) (et(j),j=1,net)
C  401 format('  E=',10f8.3)
C      write (*,402) (rt(j),j=1,nrt)
C  402 format('  R=',10f8.3)

      end

