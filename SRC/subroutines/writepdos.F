      subroutine writepdos(ext)
!! == readin pdosinput, and print out pdos files. ==
      implicit none
      integer:: ifip,ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,ndos,nkp
     & ,ibas,jsp,ifi,init,iend,ipts,fopna,fopnn,j,ndos_,ichan,isp,itet,ksp,i,ib
      integer,allocatable::idtete(:,:),ipqe(:,:,:)
      real(8),allocatable:: evlall(:,:,:),dwgtall(:,:,:,:,:),pdosp(:,:),pdosalla(:,:,:,:)
      real(8)::eminp,emaxp,ef0,eee,eminp_,emaxp_
      character*3::charnum3
c      character*100::filenm(2)
      real(8):: rydberg=13.6058d0, bin,eigen(4),vvv,wt,bin2
      character strn*120
      character*(*)::ext
      logical::a2bin,cmdopt,mlog

      integer, dimension(:),allocatable :: kpproc
      integer::numprocs,procid,ierr,itete,iteti,mpipid,ifile_handle
#if MPIK
      include "mpif.h"
#endif
      print *,' pdosdata file=','pdosdata.'//trim(ext)
      open(newunit=ifip,form='unformatted',file='pdosdata.'//trim(ext))
      read(ifip) ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,nkp !ndos,nkp mar2015
      allocate(idtete(0:4,6*nkp),ipqe(nkk1,nkk2,nkk3))
      allocate(evlall(ndhamx,nspx,nkp))             
      allocate(dwgtall(nchanp,nbas,ndhamx,nspx,nkp)) 
      read(ifip) idtete !,ipqe 
      read(ifip) evlall         
      read(ifip) dwgtall        
      read(ifip) ef0 !may2015
      close(ifip)
#if MPIK
      mlog = cmdopt('--mlog',6,0,strn) !--mlog here is taken by getarg.
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
!! Default and external option.
      eminp =-25.0/rydberg
      emaxp = 30.0/rydberg
      ndos  = 5500
      if (cmdopt('-emin=',6,0,strn)) then
         read(strn(7:),*) eminp_
         eminp = eminp_/rydberg
      endif
      if (cmdopt('-emax=',6,0,strn)) then
        read(strn(7:),*) emaxp_
        emaxp = emaxp_/rydberg
      endif
      if(cmdopt('-ndos=',6,0,strn))  then
        read(strn(7:),*) ndos_
        ndos = ndos_
      endif

c$$$      if (cmdopt('-emin=',6,0,strn)) then
c$$$         j=6
c$$$        if (a2bin(strn,eminp_,4,0,' ',j,len(strn))) eminp=eminp_/rydberg
c$$$      endif
c$$$      if (cmdopt('-emax=',6,0,strn)) then
c$$$        j =6
c$$$        if (a2bin(strn,emaxp_,4,0,' ',j,len(strn))) emaxp=emaxp_/rydberg
c$$$      endif
c$$$      if(cmdopt('-ndos=',6,0,strn))  then
c$$$        j = 6
c$$$        if (a2bin(strn,ndos_,2,0,' ',j,len(strn))) ndos=ndos_
c$$$      endif
      
!! dostet ---
c      call dostet(ndhamx,nsp,nspx,nevmin,nchanp*nbas,nkk1,nkk2,nkk3,ntete,idtete,evlall,
c     &  dwgtall, ndos, eminp+ef0, emaxp+ef0,.false.,wkd,pdosall)
      bin = (emaxp - eminp) / (ndos - 1)
      vvv = ( 3d0  -  nsp ) / ( nkk1 * nkk2 * nkk3 * 6d0 )/4d0
      allocate(pdosalla(ndos,nsp,nchanp,nbas))
!! --- Loop over tetrahedra ---
      do itet = iteti, itete
      do isp = 1, nspx
      do ib = 1, nevmin
        eigen(1:4) = evlall(ib,isp,idtete(1:4,itet))
        if( minval(eigen) > emaxp+ef0 ) cycle
        do ibas = 1,nbas
        do ichan = 1, nchanp
          wt = sum(dwgtall(ichan,ibas,ib,isp,idtete(1:4,itet))) * idtete(0,itet) * vvv
          call slinz(wt,eigen,eminp+ef0,emaxp+ef0,pdosalla(1,isp,ichan,ibas),ndos)
        enddo
        enddo
      enddo
      enddo
      enddo
#if MPIK
      call mpibc2_real( pdosalla, ndos*nsp*nchanp*nbas, 'writepdos_pdosalla' )
#endif
      if(procid==0) then !master only
      allocate(pdosp (ndos,nchanp))
      bin2 = 2d0 * bin
      do isp =1,nspx
      do ibas=1,nbas
        open(newunit=ifi,file=trim('dos.isp'//char(48+isp)//'.site'//charnum3(ibas))//'.'//trim(ext))
        write(ifi,"('#lm ordering. See the end of lmdos. relative to efermi')")
!!      DOS from finite difference of NOS 
        do ichan=1,nchanp
          do i = 2, ndos - 1
            pdosp(i,ichan)=(pdosalla(i+1,isp,ichan,ibas) - pdosalla(i-1,isp,ichan,ibas))/bin2
          enddo
          pdosp(1,ichan)    = pdosp(2,ichan)
          pdosp(ndos,ichan) = pdosp(ndos-1,ichan)
        enddo
        do ipts=1,ndos
          eee = eminp+ (ipts-1d0)*(emaxp-eminp)/(ndos-1d0)
          write(ifi,"(255(f13.5,x))")eee,pdosp(ipts,1:nchanp)
        enddo
        call fclose(ifi)
      enddo
      enddo   
      deallocate(pdosp,pdosalla)
      endif
      end
!!
