      subroutine mkewgt(lmet,wgt,qval,ndimh,evl,!ef0,def,!esmear-->BZ_W and BZ_N instead,
     .     nevec,ewgt,sumev,sumqv)
      use m_lmfinit,only: bz_w,bz_n,stdo
C- State-dependent weights for sampling BZ integration.
C ----------------------------------------------------------------------
Ci Inputs
Ci   lmet  :0 assume insulator
Ci         :1 save eigenvectors to disk
Ci         :2 read weights from file, if they are needed
Ci         :3 always make two band passes; weights never needed a priori
Ci         :4 BZ integration with 3-point scheme
Ci   wgt   :symmetry weight for this qp
Ci   qval  :valence charge; number of occupied states
Ci   ndimh :dimensions evl
Ci   evl   :eigenvalues
Cixxx   ef0   :trial Fermi level
Cixxx   def   :uncertainty in Fermi level; quantities are accumulated for
Cixxx         :ef0-def, ef0, ef0+def
Ci   esmear:Parameter that describes gaussian broadening. 
Ci         :Integer part >0 for for generalized gaussian broadening
Ci         :and is the the Methfessel-Paxton integration order
Ci         :Fractional part is the broadening width.
Ci         :Integer part <0 => Fermi-Dirac broadening used
Ci         :Fractional part is the temperature
Ci         :(see delstp.f)
Ci         :integer part above 100's digit is stripped.
Co Outputs
Co   nevec :number of states with nonzero weights
Co   ewgt  :weights
Co   sumev :sumev(*) are sumev(1..2) at Fermi levels
Co         :sumev(1) = sum of eigenvalues
Co         :sumev(2) = entropy
Co   sumqv :sum of charges 
Cr Remarks
Cr   Weights are made for generalized Gaussian broadening, where a
Cr   delta function is expressed as a hermite polynomial*gaussian,
Cr   or for Fermi-Dirac broadening (see delstp.f)
Cu Updates
Cu   17 Jan 05 Extension of esmear to Fermi distribution
Cu   31 May 00 Adapted from nfp mkewgt
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmet,ndimh,nevec!,numq
      double precision qval,wgt !esmear,
      double precision evl(ndimh),ewgt(1),sumev(2),sumqv
C ... Local parameters
      integer i,i1,i2,ie,iprint,iq,k,nord
      double precision dn,eff,fevec,fn,sn,width,wtop,wx,x
      if (qval/2 .gt. ndimh) call fexit2(-1,111,'%N Exit -1 : MKEWGT: '
     .//'Basis with ndimh=%i is insufficient to carry q/2=%d states',
     .ndimh,qval/2)
C --- Nonmetal case ---
      if (lmet .eq. 0) then
c        if (numq .ne. 1) call rxi('mkewgt: nonmetal but numq=',numq)
        fevec = qval/2d0
        nevec = fevec + 0.500001d0
        wtop = fevec-(nevec-1)
        if (dabs(wtop-1d0) .gt. 1d-6) write(stdo,"(' mkewgt (warning):'//
     .  ' uneven occupation for nonmet, top wgt=',d15.7)")wtop
        do  i = 1, nevec
          ewgt(i) = 1d0
          if (i .eq. nevec) ewgt(i) = wtop
          sumev(1) = sumev(1) + wgt*ewgt(i)*evl(i)
          sumev(2) = 0d0
          sumqv    = sumqv   + wgt*ewgt(i)
        enddo
C --- Metal case: require three fermi levels ---
c      else
c         if(numq/=3) call rxi('mkewgt: metal but numq=',numq)
c$$$        if (numq .ne. 3) call rxi('mkewgt: metal but numq=',numq)
c$$$c        width = dabs(esmear) - int(dabs(esmear))
c$$$c        nord = dsign(1d0,esmear) * mod(int(dabs(esmear)),100)
c$$$        width= abs(bz_w)
c$$$        nord = bz_n
c$$$        nevec = 0
c$$$        do  ie = 1, ndimh
c$$$          wx = 0d0
c$$$          do  k = 1, 3
c$$$            eff = ef0 + (k-2)*def
c$$$            x = (evl(ie)-eff)/width
c$$$            call delstp(nord,x,dn,fn,sn)
c$$$            ewgt(k,ie) = fn
c$$$            sumqv(k) = sumqv(k) + wgt*ewgt(k,ie)
c$$$            sumev(1,k) = sumev(1,k) + wgt*ewgt(k,ie)*evl(ie)
c$$$            sumev(2,k) = sumev(2,k) - wgt*width*sn
c$$$            wx = dmax1(wx,dabs(ewgt(k,ie)))
c$$$          enddo
c$$$          if (wx .lt. 1d-8) goto 80
c$$$C          write(stdo,268) ie,evl(ie),wx,(ewgt(k,ie),k=1,numq)
c$$$C  268     format(i5,f10.6,d13.2,8f12.6)
c$$$          nevec = ie
c$$$        enddo
c$$$  80    continue
      endif

!! this print out is confusing because we make truncation of sigm now in fpgw/exec/hqpe_sc
c$$$C --- Printout weights near frontier ---
c$$$      if (iprint() .gt. 40) then
c$$$C       write(stdo,200) lmet,numq,nevec
c$$$C 200   format(' mkewgt: lmet=',i2,'   numq=',i2,'  nevec=',i4)
c$$$        i2 = nevec
c$$$        i1 = max0(nevec-5,1)
c$$$        write(stdo,501) (i,i=i1,i2)
c$$$        write(stdo,502) (evl(i),i=i1,i2)
c$$$        do  iq = 1, numq
c$$$          write(stdo,503) iq,(ewgt(iq,i),i=i1,i2)
c$$$        enddo
c$$$  501   format(' state:',7i11)
c$$$  502   format(' evl:  ',7f11.6)
c$$$  503   format(' w',i1,':   ',7f11.6)
C        do  ie = 1, nevec
C          write(stdo,210) ie,evl(ie),(ewgt(iq,ie),iq=1,numq)
C  210     format('  state',i3,'   evl',f10.6,'   wgts:',3f12.6)
C        enddo
c$$$      endif
      end

