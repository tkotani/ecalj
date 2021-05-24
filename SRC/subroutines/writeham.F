      subroutine writeham_init(ififft)
      use m_mkqp,only: bz_nabc
      use m_MPItk,only: strprocid
      use m_lattic,only: qlat=>lat_qlat,pos=>rv_a_opos,plat=>lat_plat
      use m_lmfinit,only: nbas,
     &     alat=>lat_alat
      integer:: nkp,iq,ik1,ik2,ik3,npairmx,ikpd,ni,ififft,ib1,ib2,ikp,ip,nkk1,nkk2,nkk3
      integer,allocatable:: npair(:,:),nlat(:,:,:,:),nqwgt(:,:,:)
      real(8),allocatable:: qplist(:,:)
      real(8)::rrrr,posp(3),qp(3),pi=4d0*datan(1d0)
      complex(8):: aaaa,img=(0d0,1d0)
      character*4:: cccx
      logical:: ok
      nkk1=bz_nabc(1)
      nkk2=bz_nabc(2)
      nkk3=bz_nabc(3)
c        if(master_mpi) then ! pair-table and weight generator.  nlat,qgt are not used yet...
          nkp = nkk1*nkk2*nkk3
          allocate(qplist(3,nkp))
          iq=0
          do ik3=1,nkk3
          do ik2=1,nkk2
          do ik1=1,nkk1
             iq=iq+1
             qplist(:,iq) =  qlat(:,1)*dble(ik1-1)/nkk1
     &        +   qlat(:,2)*dble(ik2-1)/nkk2
     &        +   qlat(:,3)*dble(ik3-1)/nkk3
          enddo
          enddo
          enddo
          do ip=1,nkp
            write(6,"(' qplist:',i9,3f10.4)")ip, qplist(:,ip)
          enddo
          allocate(npair(nbas,nbas))!,pos(3,nbas))
c          pos = rv_a_opos
          npairmx=nkk1*nkk2*nkk3*2 !initial size of npairmx
          do   !this loop is just to get reasonable npairmx satisfying npairmx >npair
            npairmx = npairmx + (nkk1*nkk2*nkk3+1)*.5 !+1 added Feb2014 for 1x1x1
            allocate( nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas) )
            call gennlat(pos,nbas,plat,nkk1,nkk2,nkk3,npairmx,ok,npair,nlat,nqwgt)
            if(ok) exit
            deallocate( nlat, nqwgt )
          enddo
          print *,'We got nlat and nqwgt'
          do ib1=1,nbas
          do ib2=1,nbas
              write(6,"(3i8,' !ib1 ib2 npair -------------')") ib1,ib2,npair(ib1,ib2)
c             if(abs(sum(qwgt(1:npair(ib1,ib2),ib1,ib2))-1d0)>1d-8) call rx('buggg:qwgt sum is not unity')
              do ni = 1,npair(ib1,ib2) 
                posp =  pos(:,ib1)-pos(:,ib2) + matmul(plat,nlat(:,ni,ib1,ib2)) ! R_i+T - R_j
                rrrr = sqrt(sum(posp**2))
                write(6,"(i6,3x,3i4,i3,x,f8.3)") ni, nlat(1:3,ni,ib1,ib2),nqwgt(ni,ib1,ib2),rrrr
              enddo
          enddo
          enddo
          open(newunit=ififft,file='HamiltonianPMTInfo',form='unformatted')
          write(ififft) plat,nkk1,nkk2,nkk3,nbas,qlat
          write(ififft) pos,alat
          write(ififft) qplist
          write(ififft) npair,npairmx
          write(ififft) nlat,nqwgt
!! delta fun check: k --> T --> k 
!!     \delta_{kk'} = \sum_{T \in T(i,j)} W_T exp( i (k-k') T)
          ikpd=1
          do ikp=1,nkp
            qp = qplist(:,ikp) - qplist(:,ikpd)
          do ib1=1,nbas
          do ib2=1,nbas
            aaaa=0d0
            do ni = 1,npair(ib1,ib2) 
                aaaa =  aaaa + 1d0/(nkp*nqwgt(ni,ib1,ib2))*exp(img*2d0*pi* sum(qp*matmul(plat,nlat(:,ni,ib1,ib2))))
            enddo
            cccx=''
            if(ikp==ikpd) cccx=' <--'
            write(6,"('\delta-fun test',3f10.4,2i3,2f23.15,a)") qplist(:,ikp),ib1,ib2,aaaa,cccx
          enddo
          enddo
          enddo
          deallocate(npair)
c        endif
      end subroutine
