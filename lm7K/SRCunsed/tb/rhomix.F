      subroutine rhomix(lov,nl,nspu,nbas,dclabl,ipc,it,itmax,cnvg,mmix,
     .                  neltst,beta,tm,tj,rhol0,rhol,rho0,rho,rhlits,
     .                  rhoits,a,rms)
C- Mixing rhoc or rhon for TB-L and TB+U
C ----------------------------------------------------------------------
Ci Inputs:
Ci   a: the usual Anderson work array, leading dimension neltst
Ci   rho: rhoc or rhon from tbfrce
Ci   rho0: The "input" rho, rho_in used at first iteration
Ci   rhol: rho from tbfrce, Mulliken charges unless lov=T these are
Ci         identical to the traces of rho and do not need to be
Ci         mixed.
Ci   rhlits, rhoits are work arrays holding previous iterations
Cc Confusion
Cc   Here we use "rho" to mean either rhoc (TB-L) or rhon (TB+U)
Cc   hence the usual rho (Mulliken charges) is here denoted "rhol"
Co Outputs:
Co   See remarks
Cr Remarks
Cr   For detailed remarks on mixing, see qmix.f
Cr   On input rho holds current rhoc or rhon from tbfrce.
Cr   On output it is overwritten with the mixed rho.
Cr   If lov=T, rhol are also mixed and returned.
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      logical lov
      integer nspu,nbas,nl,it,itmax,mmix,neltst,ipc(1)
      double precision rhol(nl,2,nbas),rhol0(nl,2,nbas),
     .              rhlits(nl,nspu,nbas,0:mmix+1,2),
     .                 rho(nl**2,nl**2,nbas,nspu),
     .                rho0(nl**2,nl**2,nbas,nspu),
     .              rhoits(nl**2,nl**2,nbas,nspu,0:mmix+1,2),
     .                   a(neltst,0:mmix+1,2),rms,dclabl(1)
      double precision cnvg,beta,tm,tj(1)
C Local Variables
      integer i,j,ic,ib,ispu,nelts,nelts0,neltso,npmix,ido,nmix,amix,
     .        ipr,ipr1,ipr2
      integer iprint,i1mach
      integer okpvt,onorm
      double precision b,dsum
      character*8 clabl
   
C Heap
      integer w(1)
      common /w/ w

      nelts  = nl**4*nbas*nspu
      if (lov) then
        neltso = nl*nspu*nbas
      else
        neltso = 0
      endif
      nelts0 = nelts + neltso
      call rxx(nelts0.ne.neltst,' RHOMIX: bug, nelts0.ne.neltst')
      call rxx(nspu .ne. 1,' RHOMIX not tested for spin pol')
      npmix = min(it-1,mmix)

C --- Roll back previous iterations ---
      do  i = mmix, 0, -1
        call dcopy(nelts,rhoits(1,1,1,1,i,1),1,rhoits(1,1,1,1,i+1,1),1)
        call dcopy(nelts,rhoits(1,1,1,1,i,2),1,rhoits(1,1,1,1,i+1,2),1)
        if (lov) then
          call dcopy(neltso,rhlits(1,1,1,i,1),1,rhlits(1,1,1,i+1,1),1)
          call dcopy(neltso,rhlits(1,1,1,i,2),1,rhlits(1,1,1,i+1,2),1)
        endif
      enddo

C --- Take new rho (and Mulliken charges if lov) for this iteration ---
      if (it .eq. 1) then
        call dcopy(nelts,rho0,1,rhoits(1,1,1,1,1,1),1)
        if (lov) then
          do  ib = 1, nbas
            do  ispu = 1, nspu
              call dcopy(nl,rhol0(1,ispu,ib),1,rhlits(1,ispu,ib,1,1),1)
            enddo
          enddo
        endif
      endif
      call dcopy(nelts,rho,1,rhoits(1,1,1,1,0,2),1)
      if (lov) then
        do  ib = 1, nbas
          do  ispu = 1, nspu
            call dcopy(nl,rhol(1,ispu,ib),1,rhlits(1,ispu,ib,0,2),1)
          enddo
        enddo
      endif

C --- Build work array for amix ---
      do  i = 0, npmix
        call dcopy(nelts,rhoits(1,1,1,1,i+1,1),1,a(1,i,2),1)
        call dcopy(nelts,rhoits(1,1,1,1,i,2),1,a(1,i,1),1)
        if (lov) then
          call dcopy(neltso,rhlits(1,1,1,i+1,1),1,a(1+nelts,i,2),1)
          call dcopy(neltso,rhlits(1,1,1,i,2),1,a(1+nelts,i,1),1)
        endif
        if (i .ne. 0) then
          call daxpy(neltst,-1d0,a(1,i,2),1,a(1,i,1),1)
        endif
      enddo

C --- Mix; don't chatter about it ---
      b = beta
C      if (it .lt. 3) b = b/10
      call pshprt(0)
      ipr = iprint()
      ido = 0
      call defdr(onorm,mmix**2)
      call defi(okpvt,mmix)
      nmix = amix(neltst,npmix,mmix,ido,b,ipr,tm,w(onorm),w(okpvt),
     .            a,tj,rms)
      call popprt
      call rlse(onorm)

C --- Get new rho (and Mulliken charges if lov) from work array ---
      call dcopy(nelts,a(1,0,2),1,rhoits(1,1,1,1,0,1),1)
      call dcopy(nelts,a(1,0,2),1,rho,1)
      if (lov) then
        call dcopy(neltso,a(1+nelts,0,2),1,rhlits(1,1,1,0,1),1)
        do  ib = 1, nbas
          do  ispu = 1, nspu
            call dcopy(nl,rhlits(1,ispu,ib,0,1),1,rhol(1,ispu,ib),1)
          enddo
        enddo
      endif

C --- If rhol not mixed then construct new Mulliken charges ---
      if (.not. lov) then
        do  ib = 1, nbas
          do ispu = 1, nspu
            rhol(1,ispu,ib) = rho(1,1,ib,ispu)
            if (nl .gt. 1) then
              rhol(2,ispu,ib) = dsum(3,rho(2,2,ib,ispu),nl**2+1)
            endif
            if (nl .gt. 2) then
              rhol(3,ispu,ib) = dsum(5,rho(5,5,ib,ispu),nl**2+1)
            endif
          enddo
        enddo
      endif

C --- Copy delta's into delL for next iteration of MD or static ---
C      if (it .eq. nitmax .or. rms .lt. cnvg)
C     .  call dcopy(nelts,delta(1,1,1,1,0,1),1,delL,1)

C --- Printout ---
      ipr1 = 30
      ipr2 = 30
      if (iprint() .lt. ipr1) return
      print 100
      call awrit7(
     .  ' Iteration %i. %i elements; mixed %i of %i, beta=%d, tol=%g, '/
     .  /'rms diff: %g',' ',90,i1mach(2),it,neltst,nmix,npmix,b,cnvg,rms
     .  )
      if (nmix .gt. 0) write (*,110) (tj(i),i=1,nmix)
      if (lov .or. iprint() .gt. ipr2) then
        print *, '   Input rho ... x_i'
        do  ib = 1, nbas
          ic = ipc(ib)
          call r8tos8(dclabl(ic),clabl)
          call awrit1('Atom %i '//clabl,' ',60,i1mach(2),ib)
          if (iprint() .gt. ipr2) then
            do  i = 1, nl**2
              write (*,300) (rhoits(i,j,ib,1,1,1),j=1,nl**2)
            enddo
          endif
          if (lov) then
            call awrit2('  Mulliken charges: %n:1d',' ',
     .                  120,i1mach(2),nl,rhlits(1,1,ib,1,1))
          endif
        enddo
      endif
      if (lov .or. iprint() .gt. ipr2) then
        print *, '   Output rho ... f(x_i)'
        do  ib = 1, nbas
          ic = ipc(ib)
          call r8tos8(dclabl(ic),clabl)
          call awrit1('Atom %i '//clabl,' ',60,i1mach(2),ib)
          if (iprint() .gt. ipr2) then
            do  i = 1, nl**2
              write (*,300) (rhoits(i,j,ib,1,0,2),j=1,nl**2)
            enddo
          endif
          if (lov) then
            call awrit2('  Mulliken charges: %n:1d',' ',
     .                  120,i1mach(2),nl,rhlits(1,1,ib,0,2))
          endif
        enddo
      endif
      print *, '   Mixed rho --> rho_i+1'
      do  ib = 1, nbas
        ic = ipc(ib)
        call r8tos8(dclabl(ic),clabl)
        call awrit1('Atom %i '//clabl,' ',60,i1mach(2),ib)
        if (iprint() .gt. ipr2) then
          do  i = 1, nl**2
            write (*,300) (rhoits(i,j,ib,1,0,1),j=1,nl**2)
          enddo
        endif
        call awrit2('  Mulliken charges: %n:1d',' ',
     .               120,i1mach(2),nl,rhol(1,1,ib))
      enddo
  100 format(/' RHOMIX mixing density matrix elements:')
  110 format(' t_j :',10f8.4)
  300 format (5x,9f10.6)
      end

