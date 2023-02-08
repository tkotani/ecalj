module m_relax
  use m_ftox
  use m_lgunit,only:stdo,stdl
  public relax
  private
contains
  subroutine relax(it,indrlx,natrlx,force, p,w,basin,bas,icom)
    use m_lmfinit,only: nbas,nitrlx,slabl,ifrlx,ispec
    use m_lmfinit,only: lrlxr,rdhessr,nkillr,xtolr,gtolr,stepr
    use m_ext,only:     sname
    use m_gradzr,only :Gradzr
    use m_struc_def
    !- Relax atomic positions and volume using variable metric algorithm
    ! ----------------------------------------------------------------------
    !i Inputs:
    !i   basin: original position.
    !i   it:        iteration number
    !i   indrlx(1,i) points to the ith relaxing component and
    !i   indrlx(2,i) points to the corresponding site (see rlxstp)
    !i   natrlx:    # of relaxing degrees of freedom for atoms
    !i   force:     forces
    !i   p:         for gradzr (dimensioned in rlxstp)
    !i   w:         the hessian
    !i basin: original atomic position
    !o Outputs
    !o bas: new atomic position
    !o   icom       1 if relaxation is complete
    !o             -1 if relaxation is not complete and current gradient
    !o                is larger than prior case.
    !o                p(*,nm) holds positions for prior case
    !o              0 otherwise
    !l Local variables
    !l   wkg:       work array used by gradzr.
    !l              It should not be modified between iterations
    !r Remarks:
    !r   At version 6.10 volume and shear relaxations have been removed;
    !r   hence natrlx=natrlx.
    !r   When ctrl->mdprm(1) > 100, bas are shear coords
    !u Updates
    !u   21 Mar 06 mdprm(1)>100 signifies shear relaxation
    !u   08 Mar 06 Some improvements to gradzr; repackaged MPI
    !u   15 Feb 02 (ATP) Added MPI parallelization
    ! ----------------------------------------------------------------------
    implicit none
    integer,parameter:: nm=3
    integer :: it,nit,natrlx,icom,procid,master,mpipid,nitr
    integer :: indrlx(2,natrlx)
    real(8):: force(3,*), w(natrlx,natrlx) , p(natrlx,6), bas(:,:),basin(:,:)
    integer :: i,j,ipr,ifi,ix,lgunit,ltb,natrlx2,natrlx3, &
         ir,iprint,isw,rdm,is,idamax,nd,ns,lrlx
    parameter (nd=4,ns=6)
    logical :: lpp,cmdopt,lshr ,readhess !,a2bin
    double precision :: xtolrl,grfac,wkg(28),ddot, xv(6)
!    equivalence (stepr,mdprm(5)), (xtolr,mdprm(3)), (gtolr,mdprm(4))
    character clablj*8,dumstr*6,strn*128
    save ir,wkg
    data wkg /28*0d0/
    character(256)::lll=''
    !      stdo = lgunit(1)
    !      stdl = lgunit(2)
    master = 0
    procid = mpipid(1)
    if (procid == master) call pshpr(iprint()+30)
    call tcn('relax')
    bas=basin
    ! --- Setup ---
    nit =  nitrlx
!    lrlxr   = nint(mdprm(1))
!    rdhessr = nint(mdprm(2)) .eq. 1
!    xtolr   = mdprm(3)
!    gtolr   = mdprm(4)
!    stepr   = mdprm(5)
!    nkillr  = nint(mdprm(6))
    lshr   = lrlxr .gt. 100
    lrlx   = mod(lrlxr,100)
    call getpr(ipr)
    ! --- Make vector of positions and gradients ---
    j = 1
    if ( .NOT. lshr) then
       do  10  i = 1, natrlx
          p(j,1) = bas(indrlx(1,i),indrlx(2,i))
          p(j,2) = -force(indrlx(1,i),indrlx(2,i))
          j = j + 1
10     enddo
    else
       call dcopy(natrlx,bas,1,p,1)
       call dcopy(natrlx,force,1,p(1,2),1)
    endif
    ! --- Initialization ---
    if (it == 1) then
       ir = 0
       call dcopy(natrlx*natrlx,0d0,0,w,1)
       call dcopy(natrlx,1d0,0,w,natrlx+1)
       !   ... Read Hessian from disc
       if (rdhessr) then
          if (procid == master) then
             open(newunit=ifi,file='hssn.'//trim(sname),form='unformatted')
             readhess=.false.
             read(ifi,end=9888,err=9888) natrlx2,natrlx3 !,11
             if(natrlx2/=natrlx) call rx('relax:natlx2/=natlx')
             read(ifi,end=9888,err=9888) w
             close(ifi)
             readhess=.true.
9888         continue
             ! eadhess = rdm(ifi,2,natrlx*natrlx,' ',w,natrlx,natrlx)<0
             if (readhess) then
                write(stdo,*)' RELAX: Hessian read from disc'
             else
                write(stdo,*)' RELAX: (warning) failed to read Hessian; set to unity'
                w=0d0
                call dcopy(natrlx,1d0,0,w,natrlx+1)
             endif
          endif
          call mpibc1(w,natrlx*natrlx,4,.false.,'relax','hssn')
       else
          write(stdo,*)' RELAX: no Hessian read from disc'
       endif
    endif

    ! --- Relax ---
    grfac = 1.2d0
    ! if (cmdopt('-grfac=',6,0,strn)) then
    !    j = 6
    !    if ( .NOT. a2bin(strn,grfac,4,0,' ',j,len(strn))) then
    !       print *, 'RELAX: Ignored command line value of grfac'
    !       grfac = 1.2d0
    !    endif
    ! endif
    xtolrl = abs(xtolr)
    ! if (cmdopt('-xtolrl=',6,0,strn)) then
    !    j = 6
    !    if ( .NOT. a2bin(strn,xtolrl,4,0,' ',j,len(strn))) then
    !       print *, 'RELAX: Ignored command line value of xtolrl'
    !       xtolrl = abs(xtolr)
    !    endif
    ! endif
    ! ... the user sets xtolr _and_ gtolr.
    if (lrlx == 4) isw = 00021 + 40
    if (lrlx == 5) isw = 00121 + 40
    if (lrlx == 6) isw = 00221 + 00
    if (lrlx == 4 .OR. lrlx == 5 .OR. lrlx == 6) then
       if (xtolr == 0 .AND. gtolr == 0) call rx('RELAX: both xtolr and gtolr are zero')
       if (gtolr == 0) isw = isw-10
       if (xtolr == 0) isw = isw-20
    endif
    call gradzr(natrlx,p,w,xtolrl,stepr,xtolr,gtolr,grfac,wkg,isw,ir)
    if (ir > 0) then
       call rx1('RELAX: gradzr returned ir=%i ... aborting',ir)
    endif
    if ((lrlx == 4 .OR. lrlx == 5) .AND. ipr >= 20) then
       if (procid == master) then
          if(ir==-0) lll='converged to tolerance'
          if(ir==-1) lll='new line minimization'
          if(ir==-2) lll='bracketed root this line'
          if(ir==-3) lll='bracketed root this line'
          if(ir==-4) lll='extrapolated along this line'
          if(ir==-5) lll='extrapolated along this line'
          if(ir==-6) lll='is in trouble'
          write(stdl,ftox)' fp rlx ln ',ftof(wkg(19)),trim(lll), &
               '  dxmx=', ftof(wkg(1)*p(idamax(natrlx,p(1,nd),1),nd)), &
               '|g|=',ftof(dsqrt(ddot(natrlx,p(1,2),1,p(1,2),1)))
          if (ipr >= 20) &
               write(stdo,ftox)' fp rlx ln ',ftof(wkg(19)),trim(lll), &
               '  dxmx=', ftof(wkg(1)*p(idamax(natrlx,p(1,nd),1),nd)), &
               '|g|=',ftof(dsqrt(ddot(natrlx,p(1,2),1,p(1,2),1)))
       endif
    elseif (lrlx == 6) then
       if (procid == master) then
          write(stdl,ftox)' fp rlx Br ',-ir,'dxmx',ftof(p(idamax(natrlx,p(1,ns),1),ns)), &
               '|g|=',ftof(dsqrt(ddot(natrlx,p(1,2),1,p(1,2),1)))
          if(ipr>=20) write(stdo,ftox)' fp rlx Br',-ir,'dxmx',ftof(p(idamax(natrlx,p(1,ns),1),ns)), &
               '|g|=',ftof(dsqrt(ddot(natrlx,p(1,2),1,p(1,2),1)))
       endif
    endif
    if (ir == 0) icom = 1
    ! --- Printout  ---
    if (ipr >= 40) then
       write(stdo,*)'        Gradients:'
       print 100, (p(i,2), i = 1, natrlx)
       if (lrlx /= 4 .AND. (it /= 1 .OR. xtolr /= 0)) then
          write(stdo,*)'      Diagonal inverse Hessian:'
          print 100, (w(i,i), i = 1,natrlx)
       endif
    endif
100 format(10f8.3)
    ! --- Update atom positions ---
    call prelx1(1,1,lshr,natrlx,indrlx,p,bas)
    ! --- Write Hessian to disc ---
    if (rdhessr .AND. (icom == 1 .OR. it == nit) .OR. .TRUE. ) then
       if (procid == master) then
          open(newunit=ifi,file='hssn.'//trim(sname),form='unformatted')
          write(ifi) natrlx,natrlx,11
          write(ifi) w
          close(ifi)
       endif
    endif
    if (procid == master) call poppr
    call getpr(ipr)
    ! --- Periodically remove hessian ---
    if (nkillr > 0) then
       if (mod(it,nkillr) == 0) then
          write(stdo,ftox)'   ...  resetting hessian : iter=',it,' nkillr=',nkillr
          open(newunit=ifi,file='hssn.'//trim(sname))
          close(ifi,status="delete")
          w=0d0
          call dcopy(natrlx,1d0,0,w,natrlx+1)
       endif
    endif
    ! --- Printout ---
    lpp = (icom .eq. 1 .or. it .eq. nit) .and. ipr .ge. 30 .or.ipr .ge. 31
    if (natrlx > 0 .AND. lpp .AND. lshr) then
       xv=0d0
       call grdep2(1,natrlx,indrlx,bas,xv)
       write(stdo,ftox)' Update shear STRAIN=',ftof(xv)
    elseif (natrlx > 0 .AND. lpp) then
       write(stdo,"(/' Updated atom positions:')")
       print *,' Site   Class                      Position(relaxed)'
       do  70  j = 1, nbas
          clablj=slabl(ispec(j)) 
          write (stdo,130) j,clablj,(bas(ix,j),ifrlx(ix,j).eq.1,ix=1,3)
70     enddo
    endif
    ! ... Write new atom positions to LOG file, for use in CTRL file
    if (ipr >= 20 .AND. .NOT. lshr) then
       dumstr = 'SITE'
       do  80  j = 1, nbas
          if (j == 2) dumstr = ' '
          is=ispec(j) 
          clablj=slabl(is) 
          write(stdl,ftox)dumstr//'ATOM='//clablj,' POS=',ftof(bas(1:3,j),2),'RELAX=',ifrlx(:,j)
80     enddo
    endif
    if (icom == 0 .AND. wkg(28) < 0) icom = -1
    flush(stdo)
    flush(stdl)
130 format(i4,6x,a4,3x,3(f14.8,'(',l1,')'))
! 2023-jan moved from lmfp (not yet tested)
!    if (it==nit) then !Correct poss if this is last stepr. Set minimum gradient positions 
!       write(stdo,"(a)")'relax: restore positions for minimum g (given by relax)'
!       do i = 1, natrlx
!          bas(indrlx(1,i),indrlx(2,i)) = p(i,nm) !nm=3. this is given by relax-gradzr
!       enddo
!    endif
    call tcx('relax')
  end subroutine relax
  subroutine prelx1(mode,nm,lshr,natrlx,indrlx,p,bas)
    !- Copy vector of variables to be minimized from/to bas
    ! ----------------------------------------------------------------------
    !i Inputs
    !i   mode  :>0 copy p(*,mode) into bas(*)
    !i         :<0 copy bas(*) into p(*,mode)
    !i   nm    :which index to p that should be copied
    !i   lshr  :F p correspond to basis positions
    !i         :T p correspond to simple vector
    !i   natrlx:number of variables to relax
    !i   indrlx:permutation table mappng p -> bas
    !i   bas   :basis vectors
    !i   p     :vector of variables undergoing minimization process
    !o Outputs
    !r Remarks
    !r
    !u Updates
    !u   09 Mar 06 Extracted from relax.f for use by other routines
    ! ----------------------------------------------------------------------
    implicit none
    integer :: mode,nm,natrlx,indrlx(2,natrlx)
    double precision :: p(natrlx,nm),bas(3,*)
    logical :: lshr
    integer :: i,j,iat,ix
    ! --- Update atom positions, or vice-versa ---
    if ( .NOT. lshr) then
       j = 1
       do  60  i = 1, natrlx
          ix = indrlx(1,i)
          iat = indrlx(2,i)
          if (mode > 0) then
             bas(ix,iat) = p(j,nm)
          else
             p(j,nm) = bas(ix,iat)
          endif
          j = j + 1
60     enddo
    else
       if (mode > 0) then
          call dcopy(natrlx,p(1,nm),1,bas,1)
       else
          call dcopy(natrlx,bas,1,p(1,nm),1)
       endif
    endif
  end subroutine prelx1
  subroutine grdep2(i1,i2,indrlx,dstprm,dist)
    !- Build up actual distortion from distortion parms
    !r Distortion is added into dist
    !     implicit none
    integer :: i1,i2,indrlx(i2)
    double precision :: dstprm(i2),dist(6)
    integer :: i,ip
    do  i = i1, i2
       ip = indrlx(i)
       !       1 + 2 + 3
       if (ip == 1) then
          dist(1) = dist(1) + dstprm(i)
          dist(2) = dist(2) + dstprm(i)
          dist(3) = dist(3) + dstprm(i)
          !       1 - 2
       else if (ip == 2) then
          dist(1) = dist(1) + dstprm(i)
          dist(2) = dist(2) - dstprm(i)
          !       1-3
       else if (ip == 3) then
          dist(1) = dist(1) + dstprm(i)
          dist(3) = dist(3) - dstprm(i)
          !       4, 5, 6
       else
          dist(ip) = dist(ip) + dstprm(i)
       endif
    enddo
  end subroutine grdep2
end module m_relax
