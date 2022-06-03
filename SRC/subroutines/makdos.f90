subroutine makdos(nqp,nband,nbmx,nsp,wgts,evl,n,w,tol,emin,emax, &
     ndos,dos)
  !- Make density of states from bands
  !-----------------------------------------------------------------------
  !i  Input
  !i   nqp   :number of q-points
  !i   nband :number of bands
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   wgts  :band weights
  !i   evl   :band eigenvalues
  !i   n,w   :Methfessel-Paxton order and broadening parameters
  !i   tol   :(tol>0) allowed error in DOS due to truncating the gaussian,
  !i         :        to a finite energy range (number of bins)
  !i         :(tol<0) dimensionless energy window specifying truncation
  !i         :        of gaussian.  Energy window for which gaussian is
  !i         :        taken to be nonzero is set to -tol*w
  !i   emin, emax, ndos: energy range and number of energy mesh points
  !i   nbmx  :leading dimension of evl
  !o  Ouput
  !o    dos: density of states
  !-----------------------------------------------------------------------
  !     implicit none
  integer :: nqp,nband,nbmx,nsp,n,ndos
  double precision :: wgts(nqp),evl(nbmx,nsp,nqp),dos(0:ndos-1,nsp), &
       w,emin,emax,tol,wt,emesh
  integer :: i,isp,iband,iq,meshpt,mesh1,mesh2,mrange=999999,iprint
  double precision :: e,x,range,test,step,d,s,xx
  external delstp

  call dpzero(dos,nsp*ndos)
  step = (emax - emin) / (ndos - 1)
  if ( tol > 0d0 ) then
     do  2  i = 0, ndos-1
        x = i * step / w
        call delstp(0,x,test,s,xx)
        if ( test < tol ) then
           mrange = i + 1
           goto 3
        endif
2    enddo
     if (iprint() > 30) print *,'makdos (warning) : tol too small'
3    continue
     range = 2 * mrange * step
     test = tol
  else
     range = -tol * w
     mrange = range / ( 2 * step )
     call delstp(0,-tol/2,test,s,xx)
  endif
  if (iprint() > 30) write (*,100) range/w,2*mrange,test
  do  7  iq = 1, nqp
     wt = abs(wgts(iq)) / nsp
     do  61  iband = 1, nband
        do  6  isp = 1, nsp
           e = evl(iband,isp,iq)
           meshpt = (e - emin) / step
           mesh1 = meshpt - mrange
           mesh2 = meshpt + mrange
           if (mesh2 >= ndos) mesh2 = ndos-1
           if (mesh1 < 0) mesh1 = 0
           do  4  meshpt = mesh1, mesh2
              emesh = emin + meshpt * step
              x = (emesh - e) / w
              call delstp(n,x,d,s,xx)
              dos(meshpt,isp) = dos(meshpt,isp) + wt * d / w
4          enddo
6       enddo
61   enddo
7 enddo
100 format(/1x,'MAKDOS :  range of gaussians is ',f5.2, &
       'W (',i4,' bins).' &
       /11x,'Error estimate in DOS : ',1pe9.2,' per state.')
end subroutine makdos

