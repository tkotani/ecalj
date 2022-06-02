subroutine mkewgt(lmet,wgt,qval,ndimh,evl,&! & ef0,def,!esmear-->BZ_W and BZ_N instead,
  nevec,ewgt,sumev,sumqv)
  use m_lmfinit,only: bz_w,bz_n,stdo
  !- State-dependent weights for sampling BZ integration.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   lmet  :0 assume insulator
  !i         :1 save eigenvectors to disk
  !i         :2 read weights from file, if they are needed
  !i         :3 always make two band passes; weights never needed a priori
  !i         :4 BZ integration with 3-point scheme
  !i   wgt   :symmetry weight for this qp
  !i   qval  :valence charge; number of occupied states
  !i   ndimh :dimensions evl
  !i   evl   :eigenvalues
  ! xxx   ef0   :trial Fermi level
  ! xxx   def   :uncertainty in Fermi level; quantities are accumulated for
  ! xxx         :ef0-def, ef0, ef0+def
  !i   esmear:Parameter that describes gaussian broadening.
  !i         :Integer part >0 for for generalized gaussian broadening
  !i         :and is the the Methfessel-Paxton integration order
  !i         :Fractional part is the broadening width.
  !i         :Integer part <0 => Fermi-Dirac broadening used
  !i         :Fractional part is the temperature
  !i         :(see delstp.f)
  !i         :integer part above 100's digit is stripped.
  !o Outputs
  !o   nevec :number of states with nonzero weights
  !o   ewgt  :weights
  !o   sumev :sumev(*) are sumev(1..2) at Fermi levels
  !o         :sumev(1) = sum of eigenvalues
  !o         :sumev(2) = entropy
  !o   sumqv :sum of charges
  !r Remarks
  !r   Weights are made for generalized Gaussian broadening, where a
  !r   delta function is expressed as a hermite polynomial*gaussian,
  !r   or for Fermi-Dirac broadening (see delstp.f)
  !u Updates
  !u   17 Jan 05 Extension of esmear to Fermi distribution
  !u   31 May 00 Adapted from nfp mkewgt
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: lmet,ndimh,nevec!,numq
  double precision :: qval,wgt !esmear,
  double precision :: evl(ndimh),ewgt(1),sumev(2),sumqv
  ! ... Local parameters
  integer :: i,i1,i2,ie,iprint,iq,k,nord
  double precision :: dn,eff,fevec,fn,sn,width,wtop,wx,x
  if (qval/2 > ndimh) call fexit2(-1,111,'%N Exit -1 : MKEWGT: ' &
       //'Basis with ndimh=%i is insufficient to carry q/2=%d states', &
       ndimh,qval/2)
  ! --- Nonmetal case ---
  if (lmet == 0) then
     !        if (numq .ne. 1) call rxi('mkewgt: nonmetal but numq=',numq)
     fevec = qval/2d0
     nevec = fevec + 0.500001d0
     wtop = fevec-(nevec-1)
     if (dabs(wtop-1d0) > 1d-6) write(stdo,"(' mkewgt (warning):'// &
          ' uneven occupation for nonmet, top wgt=',d15.7)")wtop
     do  i = 1, nevec
        ewgt(i) = 1d0
        if (i == nevec) ewgt(i) = wtop
        sumev(1) = sumev(1) + wgt*ewgt(i)*evl(i)
        sumev(2) = 0d0
        sumqv    = sumqv   + wgt*ewgt(i)
     enddo
     ! --- Metal case: require three fermi levels ---
     !      else
     !         if(numq/=3) call rxi('mkewgt: metal but numq=',numq)
     !$$$        if (numq .ne. 3) call rxi('mkewgt: metal but numq=',numq)
     !$$$c        width = dabs(esmear) - int(dabs(esmear))
     !$$$c        nord = dsign(1d0,esmear) * mod(int(dabs(esmear)),100)
     !$$$        width= abs(bz_w)
     !$$$        nord = bz_n
     !$$$        nevec = 0
     !$$$        do  ie = 1, ndimh
     !$$$          wx = 0d0
     !$$$          do  k = 1, 3
     !$$$            eff = ef0 + (k-2)*def
     !$$$            x = (evl(ie)-eff)/width
     !$$$            call delstp(nord,x,dn,fn,sn)
     !$$$            ewgt(k,ie) = fn
     !$$$            sumqv(k) = sumqv(k) + wgt*ewgt(k,ie)
     !$$$            sumev(1,k) = sumev(1,k) + wgt*ewgt(k,ie)*evl(ie)
     !$$$            sumev(2,k) = sumev(2,k) - wgt*width*sn
     !$$$            wx = dmax1(wx,dabs(ewgt(k,ie)))
     !$$$          enddo
     !$$$          if (wx .lt. 1d-8) goto 80
     !$$$C          write(stdo,268) ie,evl(ie),wx,(ewgt(k,ie),k=1,numq)
     !$$$C  268     format(i5,f10.6,d13.2,8f12.6)
     !$$$          nevec = ie
     !$$$        enddo
     !$$$  80    continue
  endif

  !! this print out is confusing because we make truncation of sigm now in fpgw/exec/hqpe_sc
  !$$$C --- Printout weights near frontier ---
  !$$$      if (iprint() .gt. 40) then
  !$$$C       write(stdo,200) lmet,numq,nevec
  !$$$C 200   format(' mkewgt: lmet=',i2,'   numq=',i2,'  nevec=',i4)
  !$$$        i2 = nevec
  !$$$        i1 = max0(nevec-5,1)
  !$$$        write(stdo,501) (i,i=i1,i2)
  !$$$        write(stdo,502) (evl(i),i=i1,i2)
  !$$$        do  iq = 1, numq
  !$$$          write(stdo,503) iq,(ewgt(iq,i),i=i1,i2)
  !$$$        enddo
  !$$$  501   format(' state:',7i11)
  !$$$  502   format(' evl:  ',7f11.6)
  !$$$  503   format(' w',i1,':   ',7f11.6)
  !        do  ie = 1, nevec
  !          write(stdo,210) ie,evl(ie),(ewgt(iq,ie),iq=1,numq)
  !  210     format('  state',i3,'   evl',f10.6,'   wgts:',3f12.6)
  !        enddo
  !$$$      endif
end subroutine mkewgt

