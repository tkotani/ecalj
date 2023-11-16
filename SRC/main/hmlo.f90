program hmlo ! Get bands for given k points for MLOGenerates MTO-only Hamiltonian
  use m_HamPMT, only: plat, npair,nlat,nqwgt, ldim, nkp,qplist, ib_table,alat, ReadHamPMTInfo !, HamPMTtoHamRsMTO
  ! note.  HamPMTtoHamRsMTO do not change variables. Only generate HamRsMTO file.
  use m_zhev,only:zhev_tk4
  use m_HamRsMTO,  only: hammr,ovlmr,ndimMTO,  ReadHamRsMTO,ix,npairmx,nspx
  use m_readqplist,only: eferm,qplistsy,ndat,xdat, Readqplistsy
  use m_mpi,only: MPI__hx0fp0_rankdivider2Q, MPI__Qtask, &
       MPI__Initialize, MPI__Finalize,MPI__root, &
       MPI__Broadcast, MPI__rank,MPI__size, MPI__consoleout,MPI__barrier
  implicit none
  integer:: i,j,ikp,ib1,ib2,it,nmx,nev,jsp
  complex(8)::img=(0d0,1d0),phase
  complex(8),allocatable:: hamm(:,:),ovlm(:,:),t_zv(:,:)
  real(8),allocatable:: evl(:)
  real(8)::qp(3),pi=4d0*atan(1d0),epsovl=0d0 !1d-8
  logical:: lprint=.true.,savez=.false.,getz=.false. !dummy
  integer:: ndatx,ifsy1,ifsy2,ifsy,iix(36)
  logical:: symlcase=.true.
  call MPI__Initialize()
  if(symlcase) then ! When symlcase=T, read qplist.dat (q points list, see bndfp.F). 
     call readqplistsy()
     open(newunit=ifsy1,file='band_lmfham_spin1.dat')
     if(nspx==2) open(newunit=ifsy2,file='band_lmfham_spin2.dat')
     write(6,*)  'ndat =',ndat
  endif
  call ReadHamPMTInfo()! Read infomation for Hamiltonian (lattice structures and index of basis).
  !call HamPMTtoHamRsMTO() ! HamRsMTO (real-space Hamiltonian hammr,ovlmr) is generated.
  call ReadHamRsMTO() !Read HamRsMTO. ! If HamRSMTO exist, you can skip 'call HamPMTtoHamRsMTO()'.
  
  GetEigenvaluesForSYML: block! Get Hamitonian at k points from hammr,ovlmr(realspace), then diagnalize.
    ! bands by original ecalj (by job_band), and TB hamiltonian read by ReadHamiltonianPMTInfo.
    write(6,*)  'eferm ndimMTO=',eferm,ndimMTO 
    allocate(ovlm(1:ndimMTO,1:ndimMTO),hamm(1:ndimMTO,1:ndimMTO))
    allocate(t_zv(ndimMTO,ndimMTO),evl(ndimMTO))
    nmx = ndimMTO
    ndatx = nkp
    if(symlcase) ndatx=ndat
    do ikp=1,ndatx
       if(symlcase) then
          qp= qplistsy(:,ikp)
       else
          qp = qplist(:,ikp) 
       endif
       write(6,"(' ikp along qplist, q=',i5,3f9.4)")ikp,qp! true q(i)= 2pi/alat * qplist(i,ikp)
       do jsp=1,nspx !nsp is the number of spin.  When lso=1(Lz.Sz), nspx=1
          ovlm = 0d0
          hamm = 0d0
          do i=1,ndimMTO
             do j=1,ndimMTO
                ib1 = ib_table(ix(i)) !atomic-site index in the primitive cell
                ib2 = ib_table(ix(j))
                do it =1,npair(ib1,ib2)
                   phase=1d0/nqwgt(it,ib1,ib2)*exp(-img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
                   hamm(i,j)= hamm(i,j)+ hammr(i,j,it,jsp)*phase
                   ovlm(i,j)= ovlm(i,j)+ ovlmr(i,j,it,jsp)*phase
                enddo
             enddo
          enddo
          call zhev_tk4(ndimMTO,hamm,ovlm,nmx,nev, evl,t_zv, epsovl)!Diangonale (hamm- evl ovlm) z=0
          if(symlcase) then
             if(jsp==1) ifsy=ifsy1
             if(jsp==2) ifsy=ifsy2
             do i=1,nev
                write(ifsy,"(f15.5,f15.5,2i4)") xdat(ikp),evl(i),jsp,i
             enddo
          endif
       enddo
    enddo
    close(ifsy1)
    if(nspx==2) close(ifsy2)
    write(6,"(a)")'!!! OK! band_lmfham_spin*.dat has generated! MTO only band plot'
    write(6,"(a)")'See README_MATERIALS.org for how to make a plot for band_lmfham_spin*.dat'
    write(6,"(a)")'  For example, gnuplot scrpt can be'
    write(6,"(a)")' plot \ '
    write(6,"(a)")' "bnd001.spin1" u ($2):($3) lt 1 pt 1 w lp,\ '
    write(6,"(a)")' "bnd002.spin1" u ($2):($3) lt 1 pt 1 w lp,\ '
    write(6,"(a)")' "bnd003.spin1" u ($2):($3) lt 1 pt 1 w lp,\ '
    write(6,"(a)")' "bnd004.spin1" u ($2):($3) lt 1 pt 1 w lp,\ '
    write(6,"(a)")' "bnd005.spin1" u ($2):($3) lt 1 pt 1 w lp,\ '
    write(6,"(a)")' "bnd006.spin1" u ($2):($3) lt 1 pt 1 w lp,\ '
    write(6,"(a)")' "band_lmfham_spin1.dat" u ($1):(13.605*($2-ef)) pt 2'
    write(6,"(a)")' pause -1'
    deallocate(ovlm,hamm)
  end block GetEigenvaluesForSYML
endprogram hmlo
