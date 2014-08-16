      real(8) function nocctot(ispin, ef,esmr,qbz,wbz,eval,
     d                         nband,nqbz)
c Count the total number of electrons under Ef.
c use readeval
c ispin   = 1, paramagnetic
c           2, ferromagnetic
c ef      = fermi level 
c nband   = no. states
c nqbz    = no. k-points
      implicit none
      integer(4) ::it,is,k,ispin,nqbz,nband
      real(8):: wbz(nqbz),qbz(3,nqbz),ekt(nband),esmr,ef,wgt,wiocc,
     &          eval(nband,nqbz,ispin)
      nocctot = 0d0
      wgt       = 0d0
      do is = 1,ispin
        do k   = 1,nqbz
          ekt(:) = eval(:,k,is)
          wiocc = 0d0
        do it = 1,nband
         if(    ekt(it)  + 0.5d0*esmr < ef  ) then
            wiocc = wiocc  + 1d0
          elseif(ekt(it) - 0.5d0*esmr < ef  ) then
            wiocc  = wiocc + (ef- (ekt(it)-0.5d0*esmr))/esmr
          endif
        enddo
          nocctot = nocctot + wbz(k)* wiocc
          wgt       = wgt       + wbz(k)
        enddo
      enddo
      if(ispin==1) nocctot = nocctot*2
      print *,' Ef=',ef
      print *,' wgt nocc=',wgt,nocctot
      end

c------------------------------------------------------------------
