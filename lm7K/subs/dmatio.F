      subroutine dmatio(ldmat,ldim,nl2,nbas,nsp,nspc,sdmat,xdmat,ifi)
C- I/O printout of density-matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   ldmat :1 site-diagonal density-matrix sdmat
Ci         :2 full crystal density-matrix  xdmat
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   nl2   :nl**2
Ci   nbas  :number of atoms in the basis (input)
Ci   nsp   :2 for spin-polarized case, otherwise 1 (input)
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Cio Inputs/Outputs
Cio  sdmat :site density-matrix
Cio  xdmat :full density-matrix
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ldmat,ldim,nl2,nbas,nsp,nspc,ifi
      double precision sdmat(nl2,nl2,nbas,nsp,nspc),
     .                 xdmat(ldim,nspc,ldim,nsp)
C ... Local parameters
      integer isp,jsp,ibas,lgunit,iprint,ldx
      character*80 outs

      if (ifi .gt. 0) then
        call rx('dmatio not ready for read')
      else
        rewind (-ifi)
        do  10  isp = 1, nsp
        do  10  jsp = 1, nspc
          if (ldmat .eq. 1) then
            do  12  ibas = 1, nbas
              call awrit5('%xdmat site %i'//
     .          '%?#n==2# spin %i#%j#%?#n==2# %i#%j#',
     .          outs,len(outs),0,ibas,nsp,isp,nspc,jsp)
              call ywrm(0,outs,1,-ifi,'(9f14.9)',
     .          sdmat(1,1,ibas,isp,jsp),nl2*nl2,nl2,nl2,nl2)
   12       continue
          elseif (ldmat .eq. 2 .and. nspc .eq. 1) then
            call awrit2('%xxtal dmat%?#n==2# spin %i',outs,len(outs),
     .        0,nsp,isp)
            call ywrm(0,outs,1,-ifi,'(9f14.9)',xdmat(1,1,1,isp),
     .        ldim*ldim,ldim,ldim,ldim)
          elseif (ldmat .eq. 2 .and. nspc .eq. 2) then
          else
            call rxi('dmatio: bad ldmat',ldmat)
          endif
   10   continue

C       Noncollinear density matrix is complex; write as one large array
        if (ldmat .eq. 2 .and. nspc .eq. 2) then
          ldx = ldim*2
          call ywrm(0,'xtal dmat',2,-ifi,'(9f14.9)',xdmat,ldx*ldx,ldx,
     .      ldx,ldx)
        endif

        if (iprint() .ge. 30) call awrit1('%N dmatio: '//
     .    'wrote %?;n==1;site ;;density-matrix',' ',80,lgunit(1),ldmat)
      endif
      end

