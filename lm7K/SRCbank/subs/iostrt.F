      subroutine iostrt(nbas,nf,iter,bas,vel,eps,zeta,zacc,zsqacc,veps,
     .                  time,ifi,ierr)
C- I/O positions and velocities
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas; nf, number of degrees of freedom (set to 3*nbas-3 in initv);
Ci   iter, current MD iteration number; bas, vel, current positions and
Ci   velocities; eps, current volume strain; zeta, Hoover viscosity;
Ci   zacc, ditto accumulated over the whole MD run; zsqacc, zeta^2
Ci   accumulated over the whole MD run; veps, barostat viscosity;
Ci   time, current MD time; ifi, file handle
Co Outputs:
Co   ierr, output error code: 0 on successful read; 1 if read failed
Cr Remarks
C ----------------------------------------------------------------------
Cu Updates
Cu   01 Jun 08 Rewritten by ATP
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,nf,iter,ifi,ierr
      double precision bas(3,nbas),vel(3,nbas),eps,zeta,zacc,zsqacc,
     .                 veps,time
C Local parameters
      integer nbas0,iprint,i1mach
      double precision autime
      data autime/0.048377d0/

      ierr = 0
      if (ifi .gt. 0) then
        iter = 0
        rewind ifi
        read(ifi,end=1,err=1) nbas0,nf,iter
        if (nbas .ne. nbas0) then
          goto 1
        endif
        read(ifi,err=1) eps,zeta,zacc,zsqacc,veps,time
        read(ifi,err=1) bas
        read(ifi,err=1) vel
        if (iprint() .ge. 20) then
          call awrit0('IOSTRT: read from STRT file',' ',120,i1mach(2))
          call awrit0('        ... continuing MD run',' ',120,i1mach(2))
          call awrit5('        eps=%d, zeta=%d, zacc=%d, veps=%d,'//
     .                'time=%d',' ',120,i1mach(2),eps,zeta,zacc,veps,
     .                time*autime)
        endif
      else
        write(-ifi) nbas,nf,iter
        write(-ifi) eps,zeta,zacc,zsqacc,veps,time
        write(-ifi) bas
        write(-ifi) vel
      endif
      return
    1 continue
      if (iprint() .ge. 10) then
        call awrit0(' IOSTRT: Error reading STRT file',' ',60,i1mach(2))
      endif
      ierr = 1
      iter = 0
      end

