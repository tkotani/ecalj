program initmlo
  use m_hamMTO,only:ReadHamMTOInfo,qplist,ib_table,l_table,k_table,m_table,nsp,nkp,NMTO,orbch,ipivot, &
       nqs,qprs,itor,niqs,nlo,nlo2,ilo,mkpivot,startprint,geninitmloInfo
  implicit none
  logical(4) :: debug=.false.
  integer(4) :: ideb,ifiin,ifiout,ifind,ifrot,ifoutmat,ifinput,ifovl,ifham,iflo
  integer(4) :: i,j,n,m,s,t,l,k,access,nq,LL,NS,jold,LL2
  integer(4) :: iq,NPMT,jsp,lso,lleg,rleg
  real(8)::qp(3), epsovl,rotrl,rotimg
  complex(8),dimension(:,:),allocatable:: ovlP,hamP,ovlR,hamR
  complex(8),allocatable,dimension(:,:):: dumm,uni,rmat1,rmat2,unir,unip,lmat1,lmat2
  character(200) inname,outname!,indexname

  call startprint(2)

  close(ifinput)
  open(newunit=ideb,file="initmlo.deb")

  call ReadHamMTOinfo(0) !HamiltonianPMTInfo

  call mkpivot
  call geninitmloInfo

  if(debug)then
     write(6,*) "j,atom,EH,l,magne"
     do j=1,NMTO
        i=ipivot(j)
        write(6,*) j,i,ib_table(i),k_table(i),l_table(i),m_table(i)
        write(6,*) j,ib_table(j),k_table(j),l_table(j),m_table(j)
     end do
  end if

  LL=NMTO
  LL2=LL*LL
  allocate(rmat1(1:LL,1:LL),rmat2(1:LL,1:LL))
  allocate(lmat1(1:LL,1:LL),lmat2(1:LL,1:LL))

  open(newunit=ifiin,file="HamiltonianPMT",form="unformatted")
  open(newunit=ifiout,file="HamI",form="unformatted")
  open(newunit=ifoutmat,file="Umat",form="unformatted")

  nq=0
  write(ideb,*) "nq,qp(vector),spin"
  do
     read(ifiin,end=2019) qp,NPMT,lso,epsovl,jsp
     write(ifiout)        qp,NPMT,lso,epsovl,jsp
     if(jsp==1) nq = nq +1
     if(jsp==2) write(ideb,*) " "
     write(ideb,"(i5,3(f5.2,x),i2)") nq,qp,jsp
     allocate(ovlP(1:NPMT,1:NPMT),hamP(1:NPMT,1:NPMT))

     read(ifiin,end=2020) ovlP(1:NPMT,1:NPMT) !ovl
     read(ifiin,end=2021) hamP(1:NPMT,1:NPMT) !ham

     allocate(uni(1:NPMT,1:NPMT),unip(1:NPMT,1:NPMT),unir(1:NPMT,1:NPMT))

     unip=0d0
     do i=1,NPMT
        s=i
        if(i <= LL) s = ipivot(i)
        do j=1,NPMT
           if(s==j) unip(i,j)=1d0
        end do
     end do

     uni=transpose(unip)
     write(ifoutmat) nq,NPMT
     write(ifoutmat) uni(1:NPMT,1:NPMT)

     allocate(ovlR(1:NPMT,1:NPMT),hamR(1:NPMT,1:NPMT))
     ovlR=matmul(transpose(uni),matmul(ovlP,uni))
     hamR=matmul(transpose(uni),matmul(hamP,uni))

     write(ifiout) ovlR(1:NPMT,1:NPMT)
     write(ifiout) hamR(1:NPMT,1:NPMT)

     deallocate(ovlP,hamP,ovlR,hamR,uni,unip,unir)
  end do
2019 continue
  if(nq /= nkp)then
     write(*,*) "nq,nkp=",nq,nkp
     stop "inconsistent nq!"
  end if

  stop "return 0"
2020 stop "error! faild during reading Overwrap files!"
2021 stop "error! faild during reading Hamiltonian files!"
end program initmlo