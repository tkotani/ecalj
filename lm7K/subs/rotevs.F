      subroutine rotevs(mode,ndimz,ndims,lc,lh,evl,sigm,sigii,Z,h)
C- Various transformations on sigma, eigenvectors
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1  Rotate sigm from LDA to orbital basis
Ci         :   Given Z^LDA (orbital basis), sigm (LDA basis),
Ci         :   sig(orbital basis) = (Z^LDA+)^-1 sigm (Z^LDA)^-1
Ci         :2  Unitary transformation U LDA->QSGW eigenfunctions:
Ci             U = Z^QP (Z^LDA)^-1
Ci         :2  Make QP hamiltonian from Z^LDA, evl, U
Ci
Ci   ndimz :dimension of LDA hamiltonian generating Z^LDA
Ci   ndims :dimension of sig(LDA)
Ci
Ci   lc    :absolute value is dimension of 'core' block missing in
Ci         :either old or new basis.
Ci         :lc > 0  old basis does not contain this block.
Ci                  For mode 1, sig is enlarged by a (1:lc,1:lc)
Ci                  block depicted by 'c' in the figure in Remarks.
Ci                  The diaogonal of this block is filled by elements
Ci                  (1:lc) of given sigii.
Ci         :lc < 0  New basis does not contain this block.
Ci                  For mode 1, the first 1:lc rows and columns
Ci                  are not used.
Ci
Ci   lh    :absolute value is dimension of high-lying block missing in
Ci         :either old or new basis.
Ci         :lh > 0  old basis does not contain this block.
Ci                  For mode 1, sig is enlarged by a (1:lh,1:lh)
Ci                  block depicted by 'h' in the figure in Remarks.
Ci                  The diaogonal of this block is filled by elements
Ci                  (1:lh) of given sigii.
Ci         :lh < 0  New basis does not contain this block.
Ci                  For mode 1, sig is reduced by striking out its
Ci                  first 1:lh rows and columns.
Ci
Ci   Note: Inputs must satisfy constraint  ndims+lc+lh = ndimz
Ci
Ci   evl   :(mode 1) not used
Ci         :(mode 2) not used
Ci         :(mode 3) QP evals
Cio Inputs/Outputs
Cio  sig   :(mode 1,input)  sigma
Cio        :(mode 1,output) LDA eigenvectors
Cio        :(mode 2,input)  Z^QP  = QSGW eigenvectors
Cio        :(mode 3,input)  U = Z^QP (Z^LDA)^-1
Cio        :(mode 3,output) Hqp = [(U z^LDA)^-1]+ E(qp) (U Z^LDA)^-1
Cio  z     :(mode 1,input)  Z^LDA = LDA  eigenvectors
Cio        :(mode 1,output) (Z^LDA+)^-1 sigm (Z^LDA)^-1
Cio        :(mode 2,input)  Z^LDA
Cio        :(mode 3,input)  Z^LDA = LDA  eigenvectors
Cio        :(mode 3,output) Z is DESTROYED
Cio  h     :(mode 1) work array
Cio        :(mode 2,output) Z^QP (Z^LDA)^-1
Cio        :(mode 3) work array
Cl Local variables
Cl         :
Cr Remarks
Cr
Cr   mode 1: The self-energy matrix sig(1:ndims,1:ndims) may be dimensioned
Cr   differently from the eigenvector matrix Z.
Cr
Cr   There can be a block (1:1+lc,1:1+lc) missing from sigm (OR from Z), that
Cr   physically correspond to states included explicitly in the valence in
Cr   one case, and separately as cores in the other.
Cr
Cr   There can be a block (ndims+1:ndims+lh,ndims+1:ndims+lh) missing from sigm
Cr   (OR from Z), that originate from a larger basis in one or another case.
Cr
Cr   This code handles all four possibilities.
Cr   When sigm is missing both core and high-lying orbitals,
Cr   the sig and Z matrices have this structure:
Cr
Cr                  sigm(LDA)                                    Z(lda)
Cr              ______________________                 1 ______________________
Cr              | c                   |                  |                     |
Cr              |  c                  |                  |                     |
Cr            1 |   ______________    |             lc+1 |                     |
Cr              |   |            |    |                  |                     |
Cr              |   |            |    |                  |                     |
Cr              |   |    sigm    |    |                  |                     |
Cr              |   |            |    |                  |                     |
Cr              |   |            |    |                  |                     |
Cr        ndims |   --------------    |       ndims + lc |                     |
Cr              |                 h   |                  |                     |
Cr              |                  h  |                  |                     |
Cr   ndimz - lc -----------------------            ndimz -----------------------
Cr
Cr  When either the high or the low block is present in sigm, but missing in Z,
Cr  the sigm matrix is merely truncated and the information discarded.
Cr
Cr  When either the high or the low block is present in Z, but missing in sigm,
Cr  the sigm is enlarged by the supplied diagonal sigii for those blocks.
Cu Updates
Cu   21 Jun 07 Redesigned mode 1.  Other modes don't work.
Cu   03 Dec 07 Added mode 2
Cu   26 May 07 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ndimz,ndims,lc,lh
      double precision evl(ndimz),sigii(ndimz)
      double complex sigm(ndimz,ndimz),Z(ndimz,ndimz),h(ndimz,ndimz)
C ... Local parameters
      integer ierr,ipr,is1,iz1,is2,iz2,nrow,i,j
C     double precision xx
      double complex zer,one
      complex(8),allocatable:: wk(:,:)
      parameter (zer=(0d0,0d0),one=(1d0,0d0))

      call getpr(ipr)

      if (ndims+lc+lh-ndimz .ne. 0) call
     .  rx('ROTEVS: constraint ndims+lc+lh = ndimz not satisfied')

C     Mode 1: Rotate sigm(assumed to be stored LDA basis) to orbital basis
C     Build with these steps:
C     Step   evec     sig        Notes
C     start  Z        sig(lda)   sigm dimensioned ndims, Z dim ndimz
C       1    Z^-1     sig(lda)   inverse of Zlda, stored in Z
C       1    Z^-1     sig(lda)   sigm Z^-1 stored in wk
C       2    Z^-1     sig(orb)   (Z^-1)+ sigm Z^-1 stored in h
      if (mode .eq. 1) then

C       Overwrite Zlda with its inverse (step 1)
C       call zprm('Z_LDA',2,Z,ndimz,ndimz,ndimz)
        call zqinv('n',Z,ndimz,-ndimz**2,ndimz,h,ndimz,ierr)
C       call zprm('Z_LDA^-1',2,Z,ndimz,ndimz,ndimz)
        call rxx(ierr.lt.0,'rotevs: failed to invert Z^LDA')
C
C       Step 2a
C       Multiply middle block of sigm and Zlda; see figure in Remarks.
C       Middle block of sigm is sigm(is1:is2,is1:is2) and corresponding
C       block of (Zlda)^-1 is Z^-1(iz1:iz2,1:ndimz) where :
C                             is1    is2        iz1    iz2
C         lc>=0 and lh>=0 :   1      ndims      1+lc   ndimz-lh
C         lc<0  and lh>=0 :   1-lc   ndims      1      ndimz-lh
C         lc>=0 and lh<0  :   1      ndims+lh   1+lc   ndimz
C         lc<0  and lh>=0 :   1-lc   ndims+lh   1      ndimz
C       Note all cases satisfy constraint  ndims+lc+lh = ndimz
        if (lc .ge. 0) then
          is1 = 1
          iz1 = 1+lc
        else
          is1 = 1-lc
          iz1 = 1
        endif
        if (lh .ge. 0) then
          is2 = ndims
          iz2 = ndimz-lh
        else
          is2 = ndims+lh
          iz2 = ndimz
        endif
        nrow = is2-is1
C       Middle block is a square array.  wk = wk(nrow,ndimz)
C       call zprm('sigm_LDA',2,sigm,ndims,ndims,ndims)
        allocate(wk(ndimz,ndimz))
        call dpzero(wk,ndimz*ndimz*2)
        call zgemm('N','N',nrow,ndimz,nrow,one,sigm(is1,is1),ndims,
     .    Z,ndimz,zer,wk,ndimz)
C       Multiply into wk sigii(1..lc) Z^-1 if lc > 0
        if (lc .gt. 0) then
          do  i = 1, lc
            do  j = 1, ndimz
              wk(i,j) = wk(i,j) + sigii(i)*Z(i,j)
            enddo
          enddo
        endif
C       Multiply into wk sigii(ndimz-lh..lh) Z^-1 if lh > 0
        if (lh .gt. 0) then
          do  i = ndimz-lh+1, ndimz
            do  j = 1, ndimz
              wk(i,j) = wk(i,j) + sigii(i)*Z(i,j)
            enddo
          enddo
        endif

C       Step 3. h <- (Z^-1)+ wk Z^-1
        call zgemm('C','N',ndimz,ndimz,ndimz,one,wk,
     .    ndimz,Z,ndimz,zer,h,ndimz)

C       call zprm('sigm_orb',2,h,ndimz,ndimz,ndimz)

        deallocate(wk)

C     Mode 2:  Make U = Z^QP (Z^LDA)^-1
C      elseif (mode .eq. 2) then
CC      call zprm('evec, QSGW',2,sigm,ndimz,ndimz,ndimz)
CC      call zprm('evec, LDA',2,Z,ndimz,ndimz,ndimz)
C        call zqinv('n',Z,ndimz,-ndimz**2,ndimz,h,ndimz,ierr)
C        call rxx(ierr.lt.0,'rotevs: failed to invert Z^LDA')
C        call zgemm('N','N',ndimz,ndimz,ndimz,one,sigm,ndimz,Z,ndimz,
C     .    zer,h,ndimz)
CC       call zprm('Z^QP (Z^LDA)^-1',2,h,ndimz,ndimz,ndimz)
C
CC     Mode 3:  Make Hqp(orbital basis)
C      elseif (mode .eq. 3) then
CC       call prmx('evl, QP',evl,ndimz,ndimz,1)
CC       call zprm('evec, LDA',2,Z,ndimz,ndimz,ndimz)
CC       call zprm('U',2,sigm,ndimz,ndimz,ndimz)
CC       Step 1: h <- U (Z^LDA)
C        call zgemm('N','N',ndimz,ndimz,ndimz,one,sigm,ndimz,Z,ndimz,
C     .    zer,h,ndimz)
CC       call zprm('U Z^LDA',2,h,ndimz,ndimz,ndimz)
CC       Step 2: sig <- r+ e r, where r = [U Z^LDA]^-1.  Z=work arrray
C        call phmbls(32+64,ndimz,evl,xx,Z,xx,h,h,sigm)
CC       call zprm('Hqp',12,sigm,ndimz,ndimz,ndimz)
C
CC     Mode 4:  Make sigma(orbital basis)
CC     Define U = Z^QP (Z^LDA)^-1
CC     Note: Hqp(LDA bas) =
C      elseif (mode .eq. 4) then

      else
        call rxi('rotevs: unknown mode %i',mode)
      endif

      end
      subroutine evcflg(dc,strn,lwevec,shftqp)
C- Reads switches and data for evec IO
C ----------------------------------------------------------------------
Ci Inputs
Ci   dc:   separator for switches
Ci   strn: string containing switches
Co Outputs
Co   lwevec:0 do nothing special
Co         :1 Write evals, evecs of LDA hamiltonian to file 'evec'
Co         :2 Write evals, evecs of full hamiltonian to file 'evec'
Co         :10s digit => read qp from file
Co   shftqp:shift qp by this amount (not compatible with 10s digit lwevec)
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   26 May 07  First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character dc*1
      character strn*(*)
      integer lwevec
      double precision shftqp(3)
C ... Local parameters
      integer parg
      integer i,j,j1,j2,iv(10),lfile

      lfile = 0
      lwevec = 1
      call dpzero(shftqp,3)

      if (dc .ne. ' ') then
C   ... Return here to resume parsing for arguments
        j2 = 0
   50   continue
        j2 = j2+1
        if (strn(j2:j2) .eq. dc) goto 50
        j1 = min(len(strn),j2)
        call nwordg(strn,0,dc//' ',1,j1,j2)
        if (j2 .ge. j1) then
          if (.false.) then
          elseif (strn(j1:j1+4) .eq. 'file')  then
            lfile = 1
          elseif (strn(j1:j1+4) .eq. 'mode=')  then
            j = 0
            if (parg('mode=',2,strn(j1:),j,len(strn(j1:)),
     .        dc//' ',1,1,i,lwevec) .lt. 0) goto 999
            call isanrg(lwevec,0,2,'EVCFLG:','mode',.true.)
          elseif (strn(j1:j1+5) .eq. 'shftq=')  then
            j = 0
            if (parg('shftq=',4,strn(j1:),j,len(strn(j1:)),
     .        ', '//dc,2,3,iv,shftqp) .lt. 0) goto 999
          elseif (strn(j1:j1+4) .eq. 'shftq')  then
            shftqp(1) = 1d-7
            shftqp(2) = 2d-7
            shftqp(3) = 3d-7
          else
            goto 59
          endif
          goto 50
   59     continue
  999     call rxs('EVCFLG: failed to parse option  ', strn(j1:))
        endif
      endif

      if (lfile .ne. 0) lwevec = lwevec + 10

      end

