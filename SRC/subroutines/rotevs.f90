subroutine rotevs(mode,ndimz,ndims,lc,lh,evl,sigm,sigii,Z,h)
  !- Various transformations on sigma, eigenvectors
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1  Rotate sigm from LDA to orbital basis
  !i         :   Given Z^LDA (orbital basis), sigm (LDA basis),
  !i         :   sig(orbital basis) = (Z^LDA+)^-1 sigm (Z^LDA)^-1
  !i         :2  Unitary transformation U LDA->QSGW eigenfunctions:
  !i             U = Z^QP (Z^LDA)^-1
  !i         :2  Make QP hamiltonian from Z^LDA, evl, U
  !i
  !i   ndimz :dimension of LDA hamiltonian generating Z^LDA
  !i   ndims :dimension of sig(LDA)
  !i
  !i   lc    :absolute value is dimension of 'core' block missing in
  !i         :either old or new basis.
  !i         :lc > 0  old basis does not contain this block.
  !i                  For mode 1, sig is enlarged by a (1:lc,1:lc)
  !i                  block depicted by 'c' in the figure in Remarks.
  !i                  The diaogonal of this block is filled by elements
  !i                  (1:lc) of given sigii.
  !i         :lc < 0  New basis does not contain this block.
  !i                  For mode 1, the first 1:lc rows and columns
  !i                  are not used.
  !i
  !i   lh    :absolute value is dimension of high-lying block missing in
  !i         :either old or new basis.
  !i         :lh > 0  old basis does not contain this block.
  !i                  For mode 1, sig is enlarged by a (1:lh,1:lh)
  !i                  block depicted by 'h' in the figure in Remarks.
  !i                  The diaogonal of this block is filled by elements
  !i                  (1:lh) of given sigii.
  !i         :lh < 0  New basis does not contain this block.
  !i                  For mode 1, sig is reduced by striking out its
  !i                  first 1:lh rows and columns.
  !i
  !i   Note: Inputs must satisfy constraint  ndims+lc+lh = ndimz
  !i
  !i   evl   :(mode 1) not used
  !i         :(mode 2) not used
  !i         :(mode 3) QP evals
  ! o Inputs/Outputs
  ! o  sig   :(mode 1,input)  sigma
  ! o        :(mode 1,output) LDA eigenvectors
  ! o        :(mode 2,input)  Z^QP  = QSGW eigenvectors
  ! o        :(mode 3,input)  U = Z^QP (Z^LDA)^-1
  ! o        :(mode 3,output) Hqp = [(U z^LDA)^-1]+ E(qp) (U Z^LDA)^-1
  ! o  z     :(mode 1,input)  Z^LDA = LDA  eigenvectors
  ! o        :(mode 1,output) (Z^LDA+)^-1 sigm (Z^LDA)^-1
  ! o        :(mode 2,input)  Z^LDA
  ! o        :(mode 3,input)  Z^LDA = LDA  eigenvectors
  ! o        :(mode 3,output) Z is DESTROYED
  ! o  h     :(mode 1) work array
  ! o        :(mode 2,output) Z^QP (Z^LDA)^-1
  ! o        :(mode 3) work array
  !l Local variables
  !l         :
  !r Remarks
  !r
  !r   mode 1: The self-energy matrix sig(1:ndims,1:ndims) may be dimensioned
  !r   differently from the eigenvector matrix Z.
  !r
  !r   There can be a block (1:1+lc,1:1+lc) missing from sigm (OR from Z), that
  !r   physically correspond to states included explicitly in the valence in
  !r   one case, and separately as cores in the other.
  !r
  !r   There can be a block (ndims+1:ndims+lh,ndims+1:ndims+lh) missing from sigm
  !r   (OR from Z), that originate from a larger basis in one or another case.
  !r
  !r   This code handles all four possibilities.
  !r   When sigm is missing both core and high-lying orbitals,
  !r   the sig and Z matrices have this structure:
  !r
  !r                  sigm(LDA)                                    Z(lda)
  !r              ______________________                 1 ______________________
  !r              | c                   |                  |                     |
  !r              |  c                  |                  |                     |
  !r            1 |   ______________    |             lc+1 |                     |
  !r              |   |            |    |                  |                     |
  !r              |   |            |    |                  |                     |
  !r              |   |    sigm    |    |                  |                     |
  !r              |   |            |    |                  |                     |
  !r              |   |            |    |                  |                     |
  !r        ndims |   --------------    |       ndims + lc |                     |
  !r              |                 h   |                  |                     |
  !r              |                  h  |                  |                     |
  !r   ndimz - lc -----------------------            ndimz -----------------------
  !r
  !r  When either the high or the low block is present in sigm, but missing in Z,
  !r  the sigm matrix is merely truncated and the information discarded.
  !r
  !r  When either the high or the low block is present in Z, but missing in sigm,
  !r  the sigm is enlarged by the supplied diagonal sigii for those blocks.
  !u Updates
  !u   21 Jun 07 Redesigned mode 1.  Other modes don't work.
  !u   03 Dec 07 Added mode 2
  !u   26 May 07 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,ndimz,ndims,lc,lh
  double precision :: evl(ndimz),sigii(ndimz)
  double complex sigm(ndimz,ndimz),Z(ndimz,ndimz),h(ndimz,ndimz)
  ! ... Local parameters
  integer :: ierr,ipr,is1,iz1,is2,iz2,nrow,i,j
  !     double precision xx
  double complex zer,one
  complex(8),allocatable:: wk(:,:)
  parameter (zer=(0d0,0d0),one=(1d0,0d0))

  call getpr(ipr)

  if (ndims+lc+lh-ndimz /= 0) call &
       rx('ROTEVS: constraint ndims+lc+lh = ndimz not satisfied')

  !     Mode 1: Rotate sigm(assumed to be stored LDA basis) to orbital basis
  !     Build with these steps:
  !     Step   evec     sig        Notes
  !     start  Z        sig(lda)   sigm dimensioned ndims, Z dim ndimz
  !       1    Z^-1     sig(lda)   inverse of Zlda, stored in Z
  !       1    Z^-1     sig(lda)   sigm Z^-1 stored in wk
  !       2    Z^-1     sig(orb)   (Z^-1)+ sigm Z^-1 stored in h
  if (mode == 1) then

     !       Overwrite Zlda with its inverse (step 1)
     !       call zprm('Z_LDA',2,Z,ndimz,ndimz,ndimz)
     !         call zqinv('n',Z,ndimz,-ndimz**2,ndimz,h,ndimz,ierr)
     call matcinv(ndimz,Z)
     !       call zprm('Z_LDA^-1',2,Z,ndimz,ndimz,ndimz)
     call rxx(ierr.lt.0,'rotevs: failed to invert Z^LDA')

     !       Step 2a
     !       Multiply middle block of sigm and Zlda; see figure in Remarks.
     !       Middle block of sigm is sigm(is1:is2,is1:is2) and corresponding
     !       block of (Zlda)^-1 is Z^-1(iz1:iz2,1:ndimz) where :
     !                             is1    is2        iz1    iz2
     !         lc>=0 and lh>=0 :   1      ndims      1+lc   ndimz-lh
     !         lc<0  and lh>=0 :   1-lc   ndims      1      ndimz-lh
     !         lc>=0 and lh<0  :   1      ndims+lh   1+lc   ndimz
     !         lc<0  and lh>=0 :   1-lc   ndims+lh   1      ndimz
     !       Note all cases satisfy constraint  ndims+lc+lh = ndimz
     if (lc >= 0) then
        is1 = 1
        iz1 = 1+lc
     else
        is1 = 1-lc
        iz1 = 1
     endif
     if (lh >= 0) then
        is2 = ndims
        iz2 = ndimz-lh
     else
        is2 = ndims+lh
        iz2 = ndimz
     endif
     nrow = is2-is1
     !       Middle block is a square array.  wk = wk(nrow,ndimz)
     !       call zprm('sigm_LDA',2,sigm,ndims,ndims,ndims)
     allocate(wk(ndimz,ndimz))
     call dpzero(wk,ndimz*ndimz*2)
     call zgemm('N','N',nrow,ndimz,nrow,one,sigm(is1,is1),ndims, &
          Z,ndimz,zer,wk,ndimz)
     !       Multiply into wk sigii(1..lc) Z^-1 if lc > 0
     if (lc > 0) then
        do  i = 1, lc
           do  j = 1, ndimz
              wk(i,j) = wk(i,j) + sigii(i)*Z(i,j)
           enddo
        enddo
     endif
     !       Multiply into wk sigii(ndimz-lh..lh) Z^-1 if lh > 0
     if (lh > 0) then
        do  i = ndimz-lh+1, ndimz
           do  j = 1, ndimz
              wk(i,j) = wk(i,j) + sigii(i)*Z(i,j)
           enddo
        enddo
     endif

     !       Step 3. h <- (Z^-1)+ wk Z^-1
     call zgemm('C','N',ndimz,ndimz,ndimz,one,wk, &
          ndimz,Z,ndimz,zer,h,ndimz)

     !       call zprm('sigm_orb',2,h,ndimz,ndimz,ndimz)

     deallocate(wk)

     !     Mode 2:  Make U = Z^QP (Z^LDA)^-1
     !      elseif (mode .eq. 2) then
     !C      call zprm('evec, QSGW',2,sigm,ndimz,ndimz,ndimz)
     !C      call zprm('evec, LDA',2,Z,ndimz,ndimz,ndimz)
     !        call zqinv('n',Z,ndimz,-ndimz**2,ndimz,h,ndimz,ierr)
     !        call rxx(ierr.lt.0,'rotevs: failed to invert Z^LDA')
     !        call zgemm('N','N',ndimz,ndimz,ndimz,one,sigm,ndimz,Z,ndimz,
     !     .    zer,h,ndimz)
     !C       call zprm('Z^QP (Z^LDA)^-1',2,h,ndimz,ndimz,ndimz)

     !C     Mode 3:  Make Hqp(orbital basis)
     !      elseif (mode .eq. 3) then
     !C       call prmx('evl, QP',evl,ndimz,ndimz,1)
     !C       call zprm('evec, LDA',2,Z,ndimz,ndimz,ndimz)
     !C       call zprm('U',2,sigm,ndimz,ndimz,ndimz)
     !C       Step 1: h <- U (Z^LDA)
     !        call zgemm('N','N',ndimz,ndimz,ndimz,one,sigm,ndimz,Z,ndimz,
     !     .    zer,h,ndimz)
     !C       call zprm('U Z^LDA',2,h,ndimz,ndimz,ndimz)
     !C       Step 2: sig <- r+ e r, where r = [U Z^LDA]^-1.  Z=work arrray
     !        call phmbls(32+64,ndimz,evl,xx,Z,xx,h,h,sigm)
     !C       call zprm('Hqp',12,sigm,ndimz,ndimz,ndimz)

     !C     Mode 4:  Make sigma(orbital basis)
     !C     Define U = Z^QP (Z^LDA)^-1
     !C     Note: Hqp(LDA bas) =
     !      elseif (mode .eq. 4) then

  else
     call rxi('rotevs: unknown mode %i',mode)
  endif

end subroutine rotevs
subroutine evcflg(dc,strn,lwevec,shftqp)
  !- Reads switches and data for evec IO
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   dc:   separator for switches
  !i   strn: string containing switches
  !o Outputs
  !o   lwevec:0 do nothing special
  !o         :1 Write evals, evecs of LDA hamiltonian to file 'evec'
  !o         :2 Write evals, evecs of full hamiltonian to file 'evec'
  !o         :10s digit => read qp from file
  !o   shftqp:shift qp by this amount (not compatible with 10s digit lwevec)
  !l Local variables
  !l         :
  !r Remarks
  !r
  !u Updates
  !u   26 May 07  First created
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  character dc*1
  character strn*(*)
  integer :: lwevec
  double precision :: shftqp(3)
  ! ... Local parameters
  integer :: parg
  integer :: i,j,j1,j2,iv(10),lfile

  logical:: l_dummy_isanrg, isanrg

  lfile = 0
  lwevec = 1
  call dpzero(shftqp,3)

  if (dc /= ' ') then
     !   ... Return here to resume parsing for arguments
     j2 = 0
50   continue
     j2 = j2+1
     if (strn(j2:j2) == dc) goto 50
     j1 = min(len(strn),j2)
     call nwordg(strn,0,dc//' ',1,j1,j2)
     if (j2 >= j1) then
        if ( .FALSE. ) then
        elseif (strn(j1:j1+4) == 'file')  then
           lfile = 1
        elseif (strn(j1:j1+4) == 'mode=')  then
           j = 0
           if (parg('mode=',2,strn(j1:),j,len(strn(j1:)), &
                dc//' ',1,1,i,lwevec) < 0) goto 999
           ! ino isanrg is logical function,             call isanrg(lwevec,0,2,'EVCFLG:','mode',.true.)
           l_dummy_isanrg=isanrg(lwevec,0,2,'EVCFLG:','mode',.true.)
        elseif (strn(j1:j1+5) == 'shftq=')  then
           j = 0
           if (parg('shftq=',4,strn(j1:),j,len(strn(j1:)), &
                ', '//dc,2,3,iv,shftqp) < 0) goto 999
        elseif (strn(j1:j1+4) == 'shftq')  then
           shftqp(1) = 1d-7
           shftqp(2) = 2d-7
           shftqp(3) = 3d-7
        else
           goto 59
        endif
        goto 50
59      continue
999     call rxs('EVCFLG: failed to parse option  ', strn(j1:))
     endif
  endif

  if (lfile /= 0) lwevec = lwevec + 10

end subroutine evcflg

