      module m_psigd
      contains
      subroutine psigd(mode,ndimh,sig,eval,sigp,n123,sigd, q,isp,iout,iprx)

c psigd is also called from hambls.F. Then I found present(iout) is alwasy T even when iout is missing.
c This will be a bug in gfortran.
c GNU Fortran (GCC) 4.1.2 20080704 (Red Hat 4.1.2-44)
c

ctakao's modified version if iout exists.
c if iout exist. constant part of self-energy is overridden by ESEAVR, which should be given by
c hqpe_sc.
c
C- Approximate sigma for low,higher energies with diagonal part,
C- and further add constraints for the higher energies
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 sig is a complex square matrix; poke sig(i,i)
Ci         :  sigd is not used
Ci         :1 sigd is a real diagonal matrix; poke sigd(i)
Ci         :  for those elements that satisfy constraints.
Ci         :  sig is not used
Ci   ndimh :dimension of hamiltonian, and number of energies
Ci   eval  :list of (LDA) eigenvalues
Ci   sigp  :parameters for approximating self-energy sigma.  sigma
Ci         :is approximated by its diagonal part sigii for energies
Ci         :below a low-energy cutoff (specified nmin or emin) and
Ci         :above a low-energy cutoff (specified nmax or emax).
Ci         : arg 1: specifies how to set diagonal part sigii
Ci         :        for states above the high-energy cutoff nmax or emax
Ci         :        0 constrain sigii to be > asig+bsig*e
Ci         :        1 constrain sigii to be = asig+bsig*e
Ci         :        2 constrain sigii to be > asig and < bsig
Ci         :        3 constraint same as case 1; for this routine, there
Ci         :          is no difference.  Elsewhere,
Ci         :          arg1=3 differs in that the least-squares fit to
Ci         :          sigii (for informational purposes only, to help
Ci         :          estimate asig and bsig) is done for states between
Ci         :          efit and nmax or emax
Ci         : arg 2: nmin : usage depends on mode above.
Ci         :               mode = 0: for states 1..nmin, off-diagonal
Ci         :               parts of sig(1:nmin,1:nmin) are zeroed out.
Ci         :               mode = 1: sigd(1..nmin) is filled with emin
Ci         : arg 3: emin : usage depends on mode above.
Ci         :               mode = 0: for states e_i<emin, off-diagonal
Ci         :               parts of sig(1:i,1:i) are zeroed out.
Ci         :               mode = 1: sigd(1..nmin) is filled with emin
Ci         : arg 4: nmax : sigma for levels i>nmax are approximated by
Ci         :               sigii AND constrained according to arg 1
Ci         : arg 5: emax : (used only if nmax<=0)
Ci         :             : sigma for levels e<emax are approximated by
Ci         :               sigii AND constrained according to arg 1
Ci         : arg 6: asig : constraint used to approximate
Ci         :               sigii = asig + E * bsig  or
Ci         :               asig < sigii < bsig
Ci         : arg 7: bsig : constraint used to approximate
Ci         :               sigii = asig + E * bsig  or
Ci         :               asig < sigii < bsig
Ci         : arg 8: efit : fit sigii between efit and emax
Cio Inputs/Outputs
Cio  sig   :sigma, in LDA representation
Cio        :On output:
Cio        : *high and low states are replaced by diagonal part
Cio        :  of sigma
Cio        : *diagonal part may be altered to satisfy constraints
Co Outputs
Co   n123  :blocks sigma into lower, middle, high parts
Co         :n123(1) = 0
Co         :n123(2) = index to highest orbital in 'low' block
Co         :n123(3) = index to highest orbital in 'middle' block
Co         :n123(4) = ndimh
Cl Local variables
Cl   llow  :T if this eigenvalue is below minimum cutoff
Cl   lhigh :T if this eigenvalue is below above max cutoff
Cr Remarks
Cr
Cu Updates
Cu   19 May 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ndimh,n123(4)
      double precision eval(ndimh),sigp(10),sigd(ndimh)
      double complex sig(ndimh,ndimh)
C ... Local parameters
      logical llow,lhigh
      integer i,stdo,lgunit,ipr,nmin,nmax,modsgp,PRTE
      double precision emin,emax,asig,bsig,tol,siglin
      double complex zer,sigii
      parameter (zer=(0d0,0d0),tol=1d-7,PRTE=45)
ctakao
      integer,optional,intent(in):: iout,isp
      integer:: fopna
      integer:: ifiogw,ifiese,ispx
      real(8),optional,intent(in):: q(3)
      logical:: nexist
      real(8):: sigiir
      integer,optional,intent(in):: iprx

      logical:: isanrg, l_dummy_isanrg,oncewrite
c
      print *,'psigd',present(iout)
      print *,' iout=',iout

      call getpr(ipr)
      stdo = lgunit(1)
      call ivset(n123,1,3,0)
      n123(4) = ndimh

      modsgp = nint(sigp(1))
      if (modsgp .eq. 3) modsgp = 1
      nmin   = nint(sigp(2))
      emin   = sigp(3)
      nmax   = nint(sigp(4))
      emax   = sigp(5)
      asig   = sigp(6)
      bsig   = sigp(7)

      if(present(iout).and.iout==1) then
        continue
      else
        if (mode .eq. 0) then
          call info5(PRTE,1,0,
     .  ' hambls: approximate sigma'//
     .  '%?#(n<0)# for energies E(lda)<%d; and%-2j#%-1j#'//
     .  '%?#(n>0)# for states %-1jn=%i and below; and##%j'//
     .  '%?#(n<=0)# for energies E(lda)>%d%-2j#%-1j#'//
     .  '%?#(n>0)# for states above %-1jn=%i##%j'//
     .  '%N  state    E(lda)%8fsig_ii%4fconstraint%6fuse',
     .  nmin,emin,nmax,emax,0)
        elseif (mode .eq. 1) then
          call info5(PRTE,1,0,
     .  ' hambls: new diagonal sigma for:'//
     .  ' %?#(n>0)#%-1j %i DEEP states (E=%d)#%jno DEEP states# and'//
     .  ' %?#(n>0)#%-1j %i#no# HIGH states'//
     .  ' ',nmin,emin,ndimh-nmax,0,0)
          if (nmin .le. 0 .and. nmax .ge. ndimh) return
          call info0(PRTE,0,0,'  state    E(lda)%8fsig_ii')
          stop 'for now'
        else
          call rxi('psigd: bad mode ',mode)
        endif
      endif

      imin=ndimh
      imax=1
      do  i = 1, ndimh

C       Require evals to be ordered
        if (eval(i) .lt. eval(max(i-1,1))-tol)
     .  call rxi('psigd: eval %i not ordered',i)

C       Decide whether this eval is in low, middle, or high block
        llow  = (nmin .lt. 0 .and. eval(i) .lt. emin) .or.
     .  (nmin .ge. 0 .and. i .le. nmin)
        lhigh = (nmax .le. 0 .and. eval(i) .gt. emax) .or.
     .  (nmax .gt. 0 .and. i .gt. nmax)
        if (mode .eq. 1) llow = i .le. nmin
        if (llow) n123(2) = i
        if (.not. lhigh) n123(3) = i

C       Calculate new diagonal sigma that satisfies constraints
        sigii = 0
        if (mode .eq. 0) then
          sigii = sig(i,i)
        elseif (mode .eq. 1 .and. llow) then
          sigii = emin
        endif
        siglin = asig + bsig*eval(i)
        if (lhigh) then
          if (modsgp .eq. 0) then
            if (dble(sigii) .lt. siglin) sigii = siglin
          elseif (modsgp .eq. 1) then
            sigii = siglin
          elseif (modsgp .eq. 2) then
            if (dble(sigii) .lt. asig) sigii = asig
            if (dble(sigii) .gt. bsig) sigii = bsig
          else
Ckino isanrg is logical function,             call isanrg(modsgp,0,2,'hambls:','sig fit mode',.true.)
            l_dummy_isanrg=isanrg(modsgp,0,2,'hambls:','sig fit mode',.true.)
          endif
        endif

        if(llow .or. lhigh) then
          if(present(iout).and.iout==1) then
            inquire(file='ESEAVR',number=ifiese)
            rewind ifiese
            do iii=1,3
              read(ifiese,*) sigiir,ispx
              sigii=sigiir
              siglin=sigiir
              if(ndimh==i.and.present(iprx).and.oncewrite(2) ) write(6,"(a,2d13.6,i3)")'ESEAVR: ',sigii,ispx
              if(isp==ispx) goto 898
            enddo
            stop 'psigd: can not find ESEAVR for given isp'
 898        continue
          endif
        endif  

!!      Printout
c        if(present(iout).and.iout==1) then
c          continue
c        else
        if (mode .eq. 1) then
            if (ipr .ge. PRTE .and. (llow.or.lhigh)) then
              write(stdo,331) i,eval(i),dble(sigii)
  331         format(i6,f12.6,2x,f12.6)
            endif
        else
            if (ipr .ge. PRTE .and. lhigh .and. modsgp .eq. 2) then
              write(stdo,332)
     .      i,eval(i),dble(sig(i,i)),asig,bsig,dble(sigii)
  332         format(i6,f12.6,2x,f12.6,f7.2,',',f5.2,f12.6)
            elseif (ipr .ge. PRTE .and. lhigh) then
              write(stdo,333) i,eval(i),dble(sig(i,i)),siglin,dble(sigii)
  333         format(i6,f12.6,2x,2f12.6,f13.6)
            elseif (ipr .ge. PRTE) then
              write(stdo,334) i,eval(i),dble(sig(i,i)),dble(sigii)
  334         format(i6,f12.6,2x,f12.6,12x,f13.6)
            endif
        endif

C       Overwrite full sigma with diagonal matrix, or write to sigd
        if (llow .or. lhigh) then
          if (mode .eq. 0) then
            sig(i,:)=0d0 !call zscal(ndimh,zer,sig(1,i),1)     
            sig(:,i)=0d0 !call zscal(ndimh,zer,sig(i,1),ndimh) 
            sig(i,i) = sigii
          else
            sigd(i) = sigii
          endif
        endif
        if(.not.llow) then
          if(i<imin) imin=i
        endif
        if(.not.lhigh) imax=i
cccccccccccccccccccccccccccccccccccccccc
ctakaox
c           write(6,"('rrr:',i3,d13.5,2x,l,2x,l,2d13.5)")
c     &     i,eval(i)+dreal(sigii),llow,lhigh,sigii
cccccccccccccccccccccccccccccccccccccccc
      enddo
cccccccccccccccccccccccccccccccccccccccc
c       if (llow .or. lhigh) then
c         print *,'takaoxxx2:mode i=',mode,i
c         sig(i,:) = 0d0
c         sig(:,i) = 0d0
c         sig(i,i) = sigii
c       endif
c$$$      if(present(iout).and.iout==1) then
c$$$        inquire(file='NBANDGW',number=ifiogw)
c$$$        !print *,'ifiogw=',ifiogw
c$$$        write(ifiogw,"(3d23.15,i3,3x,2i8)") q,isp,imin,imax
c$$$c         print *,'xxx uuu xxx',ifio,q,isp,imin,imax
c$$$      endif
      end subroutine psigd
      end module m_psigd
