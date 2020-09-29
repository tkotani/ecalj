      subroutine writedossawada()
!! DOS generator
! step 0.
!     prepare ctrl.cubic
!     ctrl.cubic should contain plat= 1 0 0 0 1 0 0 0 1.
!     Set nk1=10 nk2=10 nk3=10 or someghing in ctrl.cubic
!     Anything fine, because we just like to generate tetraf.dat qlistf.dat. Eg. ctrl.cubic I have.
! step 1.
!    > job_pdos  cubic -np 4 --tetraf 
!    gives tetraf.dat and qlistf.dat (-np 1 is fine)
! step 2.
!     Caluluate eigenvalues for given qlisf.dat
!       Write eigenf.dat. One line for eigenvalues for q in qplistf.dat
!       First line is for 'dimension, minimum energy, maximum energy, division'
!     See how to read eigenf.dat below.
! step 3.
!     Now we have tetraf.dat and eigenf.dat  .
!     Run this binary generated from program writedossawada.
!     >lmf-MPIK --wdsawada (see lmv7.F)
!     Then you get dosf.dat file
!! == readin dos input, print out dos
      implicit none
      integer:: ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,ndos,nkp
     & ,ibas,jsp,ifi,init,iend,ipts,j,isp,itet,ksp,i,ib,ifip
      integer,allocatable::idtete(:,:)
      real(8),allocatable:: evlall(:,:),pdosp(:,:),pdosalla(:,:),dwgtall(:,:,:)
      real(8)::eminp,emaxp,ef0,eee,eminp_,emaxp_
      character*100::strn
c      character*100::filenm(2)
      real(8):: rydberg=13.6058d0, eigen(4),wt,tot,bin,bin2
      logical:: mlog
      integer, dimension(:),allocatable :: kpproc
      complex(8),allocatable:: ham(:,:,:)
      integer::numprocs,procid,ierr,itete,iteti,mpipid,ifile_handle,ikp
#if MPIK
      include "mpif.h"
#endif
      open(newunit=ifip,form='unformatted',file='tetraf.dat')
      read(ifip) ndhamx,nkp,ntete
      allocate(idtete(0:4,6*nkp))
      read(ifip) idtete 
      close(ifip)
!!
      ndos=1000
      open(newunit=ifip,form='formatted',file='eigenf.dat')
      read(ifip,*) nevmin
      allocate(evlall(nevmin,nkp),dwgtall(nbas,nevmin,nkp))
      nbas=nevmin
      do ikp=1,nkp
         read(ifip,*) evlall(1:nevmin,ikp)
c         read(ifip,*) dwgtall(ibas,1:nevmin,ikp)
c         enddo   
c         write(6,"(i5,1000f10.3)") ikp,evlall(1:nevmin,ikp)
      enddo
      close(ifip)
!!      
      open(newunit=ifip,form='formatted',file='hamiltonian.dat')
      read(ifip,*)
      allocate(ham(nevmin,nevmin,nkp))
      read(ifip,*)ham(1:nevmin,1:nevmin,1:nkp)
      close(ifip)

      eminp =  minval(evlall)-0.5
      emaxp  = maxval(evlall)+0.5
      write(6,*) 'read eigenvalue data from eigenf.dat nkp=',nkp
      ef0=0d0
#if MPIK
      mlog=.false.
c      mlog = cmdopt('--mlog',6,0,strn) !--mlog here is taken by getarg.
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      allocate (kpproc(0:numprocs), stat=ierr)
      call dstrbp(ntete,numprocs,1,kpproc(0))
      iteti = kpproc(procid)
      itete = kpproc(procid+1)-1
c      print *,'ppppp',numprocs,procid,iteti,itete
#else
      iteti = 1
      itete = ntete
#endif
!! dostet ---
c      call dostet(ndhamx,nsp,nspx,nevmin,nchanp*nbas,nkk1,nkk2,nkk3,ntete,idtete,evlall,
c     &  dwgtall, ndos, eminp+ef0, emaxp+ef0,.false.,wkd,pdosall)
      bin = (emaxp - eminp) / (ndos - 1)
c      vvv = ( 3d0  -  nsp ) / ( nkk1 * nkk2 * nkk3 * 6d0 )/4d0
      allocate(pdosalla(ndos,0:nbas))
!! --- Loop over tetrahedra ---
      do itet = iteti, itete
      do ib = 1, nevmin
        eigen(1:4) = evlall(ib,idtete(1:4,itet))
        if( minval(eigen) > emaxp+ef0 ) cycle
        if( maxval(eigen) < eminp+ef0 ) cycle
        do ibas = 0,nbas
           if(ibas==0) wt = idtete(0,itet)
           if(ibas/=0) wt = sum(dwgtall(ibas,ib,idtete(1:4,itet))) * idtete(0,itet)
           call slinz(wt,eigen,eminp+ef0,emaxp+ef0,pdosalla(:,ibas),ndos)
        enddo
      enddo
      enddo
#if MPIK
      call mpibc2( pdosalla, ndos, 4 , mlog , 'writedossawada' , 'pdosalla' )
#endif
      if(procid==0) then        !master only
         allocate(pdosp (ndos,0:nbas))
         bin2 = 2d0 * bin
         ifi= ifile_handle()
         open(ifi,file='dosf.dat')
         do ibas=0,nbas
         do i = 2, ndos - 1
            pdosp(i,ibas)=(pdosalla(i+1,ibas) - pdosalla(i-1,ibas))/bin2
         enddo
         enddo
         pdosp(1,:)    = pdosp(2,:)
         pdosp(ndos,:) = pdosp(ndos-1,:)
         tot= sum(pdosp(1:ndos,:))*bin
         do ipts=1,ndos
            eee = eminp+ (ipts-1d0)*(emaxp-eminp)/(ndos-1d0)
            write(ifi,"(10000(f13.5,x))")eee,(pdosp(ipts,ibas)/tot,ibas=0,nbas)
         enddo
         close(ifi)
      endif
      end
