subroutine dostet(nbmx,nsp,nspx,nevmx,nchan,n1,n2,n3,ntet,idtet,eband,doswt,npts,emin,emax,lidos,wk,zos)! Density of states to third order by tetrahedron method
  !i Inputs:
  !i   nbmx, first dim. of eband;
  !i   nsp=1 spin degenerate, =2 non-deg;
  !i   nspx: 1 for spin up and down coupled, otherwise nsp
  !i   nevmx, no. of bands;
  !i   nchan, no. of DOS channels; n1,n2,n3;
  !i   ntet, idtet, o/p from tetirr
  !i   eband, work array for bands; doswt, work array for weights;
  !i   npts, no of points in energy range: [emin, emax];
  !i   lidos :F zos = dos
  !i         :T zos = energy integral of dos
  !o Outputs: zos : DOS (or integrated DOS for lidos=T) for each spin and nchan
  implicit none
  integer :: nchan,nsp,nspx,nbmx,npts,ntet,idtet(0:4,*),n1,n2,n3,nevmx,isp,ib,i,itet,ichan,iq1234(4),nspc,jsp,ksp
  real(8) :: eband(nbmx,nspx,1),emin,emax,wk(npts),zos(npts,nsp,nchan),doswt(nchan,nbmx,nsp,1),bin,eigen(4),v,wt,ebot,dmin1
  logical :: lidos
  if (npts <= 1 .OR. npts <= 2 .AND. .NOT. lidos) call rx1('dostet: npts(=%i) too small for DOS : require npts>2',npts)
  nspc = nsp / nspx
  zos=0d0
  bin = (emax - emin) / (npts - 1)
  v = ( 3d0  -  nsp ) / ( n1 * n2 * n3 * 6d0 )
  tetrahdronloop: do  5  itet = 1, ntet
     iq1234 = idtet(1:4,itet)
     do  4  isp = 1, nspx !Loop over spins and sum over bands ---
        do  3  ib = 1, nevmx
           eigen = eband(ib,isp,iq1234)
           if (minval(eigen) > emax) cycle
           do  jsp = 1, nspc
              ksp = max(jsp,isp) !ksp is isp for uncoupled spins, and jsp for coupled spins
              do  ichan = 1, nchan !       ... Accumulate no. states assuming constant wt from this tet
                 wt = sum( doswt(ichan,ib,ksp,iq1234) ) * idtet(0,itet) * v / 4d0
                 call slinz(wt,eigen,emin,emax,zos(1,ksp,ichan),npts)
              enddo
           enddo
3       enddo
4    enddo
5 enddo tetrahdronloop
  if (lidos) return
  do 111  isp  = 1, nsp !DOS from finite difference of NOS ---
     do ichan = 1, nchan
        zos(2:,  isp,ichan) = [((zos(i+1,isp,ichan) - zos(i-1,isp,ichan)) /bin/2d0, i=2,npts-1)]
        zos(1,   isp,ichan) = zos(2,isp,ichan)
        zos(npts,isp,ichan) = zos(npts-1,isp,ichan)
     enddo
111 enddo
end subroutine dostet
