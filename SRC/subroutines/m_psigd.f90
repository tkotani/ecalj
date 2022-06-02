module m_psigd
contains
  subroutine psigd(mode,ndimh,sig,eval,sigp,n123,sigd, q,isp,iout,iprx)
    use m_lgunit,only:stdo

    ! psigd is also called from hambls.F. Then I found present(iout) is alwasy T even when iout is missing.
    ! This will be a bug in gfortran.
    ! GNU Fortran (GCC) 4.1.2 20080704 (Red Hat 4.1.2-44)


    ! akao's modified version if iout exists.
    ! if iout exist. constant part of self-energy is overridden by ESEAVR, which should be given by
    ! hqpe_sc.

    !- Approximate sigma for low,higher energies with diagonal part,
    !- and further add constraints for the higher energies
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :0 sig is a complex square matrix; poke sig(i,i)
    !i         :  sigd is not used
    !i         :1 sigd is a real diagonal matrix; poke sigd(i)
    !i         :  for those elements that satisfy constraints.
    !i         :  sig is not used
    !i   ndimh :dimension of hamiltonian, and number of energies
    !i   eval  :list of (LDA) eigenvalues
    !i   sigp  :parameters for approximating self-energy sigma.  sigma
    !i         :is approximated by its diagonal part sigii for energies
    !i         :below a low-energy cutoff (specified nmin or emin) and
    !i         :above a low-energy cutoff (specified nmax or emax).
    !i         : arg 1: specifies how to set diagonal part sigii
    !i         :        for states above the high-energy cutoff nmax or emax
    !i         :        0 constrain sigii to be > asig+bsig*e
    !i         :        1 constrain sigii to be = asig+bsig*e
    !i         :        2 constrain sigii to be > asig and < bsig
    !i         :        3 constraint same as case 1; for this routine, there
    !i         :          is no difference.  Elsewhere,
    !i         :          arg1=3 differs in that the least-squares fit to
    !i         :          sigii (for informational purposes only, to help
    !i         :          estimate asig and bsig) is done for states between
    !i         :          efit and nmax or emax
    !i         : arg 2: nmin : usage depends on mode above.
    !i         :               mode = 0: for states 1..nmin, off-diagonal
    !i         :               parts of sig(1:nmin,1:nmin) are zeroed out.
    !i         :               mode = 1: sigd(1..nmin) is filled with emin
    !i         : arg 3: emin : usage depends on mode above.
    !i         :               mode = 0: for states e_i<emin, off-diagonal
    !i         :               parts of sig(1:i,1:i) are zeroed out.
    !i         :               mode = 1: sigd(1..nmin) is filled with emin
    !i         : arg 4: nmax : sigma for levels i>nmax are approximated by
    !i         :               sigii AND constrained according to arg 1
    !i         : arg 5: emax : (used only if nmax<=0)
    !i         :             : sigma for levels e<emax are approximated by
    !i         :               sigii AND constrained according to arg 1
    !i         : arg 6: asig : constraint used to approximate
    !i         :               sigii = asig + E * bsig  or
    !i         :               asig < sigii < bsig
    !i         : arg 7: bsig : constraint used to approximate
    !i         :               sigii = asig + E * bsig  or
    !i         :               asig < sigii < bsig
    !i         : arg 8: efit : fit sigii between efit and emax
    ! o Inputs/Outputs
    ! o  sig   :sigma, in LDA representation
    ! o        :On output:
    ! o        : *high and low states are replaced by diagonal part
    ! o        :  of sigma
    ! o        : *diagonal part may be altered to satisfy constraints
    !o Outputs
    !o   n123  :blocks sigma into lower, middle, high parts
    !o         :n123(1) = 0
    !o         :n123(2) = index to highest orbital in 'low' block
    !o         :n123(3) = index to highest orbital in 'middle' block
    !o         :n123(4) = ndimh
    !l Local variables
    !l   llow  :T if this eigenvalue is below minimum cutoff
    !l   lhigh :T if this eigenvalue is below above max cutoff
    !r Remarks
    !r
    !u Updates
    !u   19 May 03 First created
    ! ----------------------------------------------------------------------
    implicit none
    ! ... Passed parameters
    integer :: mode,ndimh,n123(4),iii,imax,imin
    double precision :: eval(ndimh),sigp(10),sigd(ndimh)
    double complex sig(ndimh,ndimh)
    ! ... Local parameters
    logical :: llow,lhigh
    integer :: i,ipr,nmin,nmax,modsgp,PRTE
    double precision :: emin,emax,asig,bsig,tol,siglin
    double complex zer,sigii
    parameter (zer=(0d0,0d0),tol=1d-7,PRTE=45)
    ! akao
    integer,optional,intent(in):: iout,isp
    integer:: ifiese,ispx !ifiogw,
    real(8),optional,intent(in):: q(3)
    !      logical:: nexist
    real(8):: sigiir
    integer,optional,intent(in):: iprx

    logical:: isanrg, l_dummy_isanrg,oncewrite

    !     print *,'psigd:',present(iout)
    !      print *,' iout=',iout

    call getpr(ipr)
    !      stdo = lgunit(1)
    call ivset(n123,1,3,0)
    n123(4) = ndimh

    modsgp = nint(sigp(1))
    if (modsgp == 3) modsgp = 1
    nmin   = nint(sigp(2))
    emin   = sigp(3)
    nmax   = nint(sigp(4))
    emax   = sigp(5)
    asig   = sigp(6)
    bsig   = sigp(7)

    if(present(iout) .AND. iout==1) then
       continue
    else
       if (mode == 0) then
          call info5(PRTE,1,0, &
               ' hambls: approximate sigma'// &
               '%?#(n<0)# for energies E(lda)<%d; and%-2j#%-1j#'// &
               '%?#(n>0)# for states %-1jn=%i and below; and##%j'// &
               '%?#(n<=0)# for energies E(lda)>%d%-2j#%-1j#'// &
               '%?#(n>0)# for states above %-1jn=%i##%j'// &
               '%N  state    E(lda)%8fsig_ii%4fconstraint%6fuse', &
               nmin,emin,nmax,emax,0)
       elseif (mode == 1) then
          call info5(PRTE,1,0, &
               ' hambls: new diagonal sigma for:'// &
               ' %?#(n>0)#%-1j %i DEEP states (E=%d)#%jno DEEP states# and'// &
               ' %?#(n>0)#%-1j %i#no# HIGH states'// &
               ' ',nmin,emin,ndimh-nmax,0,0)
          if (nmin <= 0 .AND. nmax >= ndimh) return
          call info0(PRTE,0,0,'  state    E(lda)%8fsig_ii')
          stop 'for now'
       else
          call rxi('psigd: bad mode ',mode)
       endif
    endif

    imin=ndimh
    imax=1
    do  i = 1, ndimh

       !       Require evals to be ordered
       if (eval(i) < eval(max(i-1,1))-tol) &
            call rxi('psigd: eval %i not ordered',i)

       !       Decide whether this eval is in low, middle, or high block
       llow  = (nmin .lt. 0 .and. eval(i) .lt. emin) .or. &
            (nmin .ge. 0 .and. i .le. nmin)
       lhigh = (nmax .le. 0 .and. eval(i) .gt. emax) .or. &
            (nmax .gt. 0 .and. i .gt. nmax)
       if (mode == 1) llow = i <= nmin
       if (llow) n123(2) = i
       if ( .NOT. lhigh) n123(3) = i

       !       Calculate new diagonal sigma that satisfies constraints
       sigii = 0
       if (mode == 0) then
          sigii = sig(i,i)
       elseif (mode == 1 .AND. llow) then
          sigii = emin
       endif
       siglin = asig + bsig*eval(i)
       if (lhigh) then
          if (modsgp == 0) then
             if (dble(sigii) < siglin) sigii = siglin
          elseif (modsgp == 1) then
             sigii = siglin
          elseif (modsgp == 2) then
             if (dble(sigii) < asig) sigii = asig
             if (dble(sigii) > bsig) sigii = bsig
          else
             ! ino isanrg is logical function,             call isanrg(modsgp,0,2,'hambls:','sig fit mode',.true.)
             l_dummy_isanrg=isanrg(modsgp,0,2,'hambls:','sig fit mode',.true.)
          endif
       endif

       if(llow .OR. lhigh) then
          if(present(iout) .AND. iout==1) then
             inquire(file='ESEAVR',number=ifiese)
             rewind ifiese
             do iii=1,2
                read(ifiese,*,err=898) sigiir,ispx !error means we use sigii for isp=1 for isp=2 !sep2012takao
                sigii=sigiir
                siglin=sigiir
                if(ndimh==i .AND. present(iprx) .AND. oncewrite(2) ) write(6,"(a,2d13.6,i3)")' ESEAVR: ',sigii,ispx
                if(isp==ispx) goto 898
             enddo
             call rx('psigd: No ESEAVR file (given by hqpe_sc)! psigd can not find ESEAVR for given isp')
898          continue
          endif
       endif

       !!      Printout
       !        if(present(iout).and.iout==1) then
       !          continue
       !        else
       if (mode == 1) then
          if (ipr >= PRTE .AND. (llow .OR. lhigh)) then
             write(stdo,331) i,eval(i),dble(sigii)
331          format(i6,f12.6,2x,f12.6)
          endif
       else
          if (ipr >= PRTE .AND. lhigh .AND. modsgp == 2) then
             write(stdo,332) &
                  i,eval(i),dble(sig(i,i)),asig,bsig,dble(sigii)
332          format(i6,f12.6,2x,f12.6,f7.2,',',f5.2,f12.6)
          elseif (ipr >= PRTE .AND. lhigh) then
             write(stdo,333) i,eval(i),dble(sig(i,i)),siglin,dble(sigii)
333          format(i6,f12.6,2x,2f12.6,f13.6)
          elseif (ipr >= PRTE) then
             write(stdo,334) i,eval(i),dble(sig(i,i)),dble(sigii)
334          format(i6,f12.6,2x,f12.6,12x,f13.6)
          endif
       endif

       !       Overwrite full sigma with diagonal matrix, or write to sigd
       if (llow .OR. lhigh) then
          if (mode == 0) then
             sig(i,:)=0d0 !call zscal(ndimh,zer,sig(1,i),1)
             sig(:,i)=0d0 !call zscal(ndimh,zer,sig(i,1),ndimh)
             sig(i,i) = sigii
          else
             sigd(i) = dble(sigii)
          endif
       endif
       if( .NOT. llow) then
          if(i<imin) imin=i
       endif
       if( .NOT. lhigh) imax=i
       ! cccccccccccccccccccccccccccccccccccccc
       ! akaox
       !           write(6,"('rrr:',i3,d13.5,2x,l,2x,l,2d13.5)")
       !     &     i,eval(i)+dreal(sigii),llow,lhigh,sigii
       ! cccccccccccccccccccccccccccccccccccccc
    enddo
    ! cccccccccccccccccccccccccccccccccccccc
    !       if (llow .or. lhigh) then
    !         print *,'takaoxxx2:mode i=',mode,i
    !         sig(i,:) = 0d0
    !         sig(:,i) = 0d0
    !         sig(i,i) = sigii
    !       endif
    !$$$      if(present(iout).and.iout==1) then
    !$$$        inquire(file='NBANDGW',number=ifiogw)
    !$$$        !print *,'ifiogw=',ifiogw
    !$$$        write(ifiogw,"(3d23.15,i3,3x,2i8)") q,isp,imin,imax
    !$$$c         print *,'xxx uuu xxx',ifio,q,isp,imin,imax
    !$$$      endif
  end subroutine psigd
end module m_psigd
