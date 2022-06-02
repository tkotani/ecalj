subroutine maknos(nqp,nband,nbmx,nsp,wgts,evl,n,w,tol,emin,emax, &
     ndos,dos)
  use m_ftox
  use m_lmfinit,only:stdo
  !- Make density of states from bands
  !-----------------------------------------------------------------------
  !i  Input
  !i    nqp : number of q-points; nband : number of bands;
  !i    nsp : 2 for spin polarised bands, 1 otherwise;
  !i    wgts, evl : weights and bands (eigenvalues);
  !i    nbmx : first dimension of evl ;
  !i    n   : n>0 Methfessel-Paxton polynomial order
  !i        : n<0 sampling done with Fermi-Dirac statistics
  !i    w   : n>0 gaussian width in Methfessel-Paxton integration (Ry)
  !i        : n<0 Temperature for Fermi distribution (Ry)
  !i    tol : allowed error in DOS due to truncating the gaussian,
  !i          if negative on entry, range is set to -tol*W
  !i    emin, emax, ndos; energy range and number of energy mesh points
  !o  Ouput
  !o    dos: integrated DOS
  !u Updates
  !u   2 Nov 1995 (JEK) returns spin-polarized integrated dos
  !-----------------------------------------------------------------------
  !     implicit none
  integer :: nqp,nband,nbmx,nsp,n,ndos
  double precision :: wgts(nqp),evl(nbmx,nsp,nqp),dos(0:ndos-1,nsp), &
       w,emin,emax,tol,wt,emesh
  integer :: i,isp,iband,iq,meshpt,mesh1,mesh2,mrange,iprint,i1mach
  double precision :: e,x,range,test,step,d,s,xx
  external delstp

  call dpzero(dos,nsp*ndos)
  step = (emax - emin) / (ndos - 1)
  if ( tol > 0d0 ) then
     do  2  i = 0, ndos-1
        x = i * step / w
        call delstp(n,x,test,s,xx)
        if ( test < tol ) then
           mrange = i + 1
           goto 3
        endif
2    enddo
     if (iprint() > 30) print *,'maknos (warning) : tol too small'
3    continue
     range = 2 * mrange * step
     test = tol
  else
     range = -tol * w
     mrange = range / ( 2 * step )
     call delstp(n,-tol/2,test,s,xx)
  endif
  if (iprint() >= 40) write(stdo,ftox) ' MAKNOS: range=',ftof(range/w), &
       ' (',2*mrange,'bins) DOS error estimate=',ftof(test),'per state'
  do  7  iq = 1, nqp
     wt = abs(wgts(iq)) / nsp
     do  61  iband = 1, nband
        do  6  isp = 1, nsp
           e = evl(iband,isp,iq)
           meshpt = (e - emin) / step
           mesh1 = meshpt - mrange
           mesh2 = meshpt + mrange
           if (mesh2 >= ndos) mesh2 = ndos-1
           call rxx(mesh1 .lt. 0,'MAKNOS: emin too large')
           do  4  meshpt = mesh1, mesh2
              emesh = emin + meshpt * step
              x = (emesh - e) / w
              call delstp(n,x,d,s,xx)
              dos(meshpt,isp) = dos(meshpt,isp) + wt * (1d0 - s)
4          enddo
           do  5  meshpt = mesh2+1, ndos-1
              dos(meshpt,isp) = dos(meshpt,isp) + wt
5          enddo
6       enddo
61   enddo
7 enddo
  !     call prmx('nos',dos,ndos,ndos,nsp)
  !  100 format(1x,'MAKNOS:  range of gaussians is ',f5.2,'W (',i4,
  !     .  ' bins).'/10x,'Error estimate in DOS : ',1pe8.2,' per state.')
end subroutine maknos

