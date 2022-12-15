module m_writeham
  use m_MPItk,only: master_mpi
  integer,private:: ififft,ifspec
  logical,private:: initset1=.true.,writeham=.false.
contains
  subroutine m_writeham_init()
    use m_MPItk,only: strprocid
    use m_lattic,only: qlat=>lat_qlat,pos=>rv_a_opos,plat=>lat_plat
    use m_lmfinit,only: nbas,  alat=>lat_alat
    use m_gennlat,only: nlat,npair,nqwgt,npairmx
    use m_mkqp,only: bz_nabc
    integer:: nkp,iq,ik1,ik2,ik3,ikpd,ni,ib1,ib2,ikp,ip,nkk1,nkk2,nkk3
    real(8),allocatable:: qplist(:,:) !this should be taken from m_qplist in future
    real(8)::rrrr,posp(3),qp(3),pi=4d0*datan(1d0)
    complex(8):: aaaa,img=(0d0,1d0)
    character*4:: cccx
    writeham=.true.
    nkk1=bz_nabc(1)
    nkk2=bz_nabc(2)
    nkk3=bz_nabc(3)
    nkp = nkk1*nkk2*nkk3
    allocate(qplist(3,nkp))
    iq=0
    do ik3=1,nkk3
       do ik2=1,nkk2
          do ik1=1,nkk1
             iq=iq+1
             qplist(:,iq)=matmul(qlat,[dble(ik1-1)/nkk1,dble(ik2-1)/nkk2,dble(ik3-1)/nkk3])
          enddo
       enddo
    enddo
    do ip=1,nkp
       write(6,"(' qplist:',i9,3f10.4)")ip, qplist(:,ip)
    enddo
    do ib1=1,nbas
       do ib2=1,nbas
          write(6,"(3i8,' !ib1 ib2 npair -------------')") ib1,ib2,npair(ib1,ib2)
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
    !     ! delta fun check: k --> T --> k
    !     !     \delta_{kk'} = \sum_{T \in T(i,j)} W_T exp( i (k-k') T)
    ikpd=1
    do ikp=1,nkp
       qp = qplist(:,ikp) - qplist(:,ikpd)
       do ib1=1,nbas
          do ib2=1,nbas
             aaaa=0d0
             do ni = 1,npair(ib1,ib2)
                aaaa=aaaa+1d0/(nkp*nqwgt(ni,ib1,ib2))&
                     *exp(img*2d0*pi*sum(qp*matmul(plat,nlat(:,ni,ib1,ib2))))
             enddo
             cccx=''
             if(ikp==ikpd) cccx=' <--'
             write(6,"('\delta-fun test',3f10.4,2i3,2f23.15,a)") qplist(:,ikp),ib1,ib2,aaaa,cccx
          enddo
       enddo
    enddo
  end subroutine m_writeham_init
  !--------------------------------------------
  subroutine m_writeham_write()
    use m_lmfinit,only: sspec=>v_sspec,lso,nsp,ispec
    use m_lmfinit,only: nlmto,stdo,slabl
    use m_hamindex, only: ngrp,symops,norbmto,ibastab,ltab,ktab,offl,ib_table,k_table,l_table
    integer:: ldim,iorb,ib,is,i,jobgw
    character spid*8
    !! --- get index for hamiltonian for m_hamindex takao june2009
    !! these are used in sigm mode(QSGW).
    !! memo:
    !!  ib = atom index
    !!  ltab= L (angular momentum index)
    !!  ktab=  =1 for EH, =2 for EH2, =3 for lo
    write(stdo,*)'mmmm m_writeham_init'
    if( .NOT. initset1  ) return
    initset1 = .false.
    ldim=nlmto  !     write(6,*) ' ib l  k offl(iorb)+1  offl(iorb)+2*l+1  trim(spec)'
    write(stdo,*)'mmmm m_writeham_init',ldim,norbmto,ibastab
    !      print *,'offl',offl,'ltab',ltab,'ktab',ktab
    print *,' norbmto mmmmmmmmmm',norbmto
    write(6,*) ' i  ib l  k trim(spec)'
    if(writeham) open(newunit=ifspec,file="atmspc")
    if(writeham) write(ifspec,*) ldim,nsp
    do i= 1, ldim
       ib   = ib_table(i)
       is   = ispec(ib) 
       spid = slabl(is) 
       write(6,"(i3,x,3i3,x,a)")i, ib_table(i),l_table(i),k_table(i),trim(spid)
       write(ifspec,"(i3,x,a,x,3i4)") i,trim(spid),ib_table(i),l_table(i),k_table(i) !sakakibara
    enddo
    if(writeham) write(ififft)ldim,lso,nsp
    if(writeham) write(ififft)ib_table,l_table,k_table
    if(writeham) close(ififft)
    if(writeham) close(ifspec)
!    deallocate(ib_table,l_table,k_table)
  end subroutine m_writeham_write
end module m_writeham
