!! getq mode.  Current version is not for spin dependent nor many restrictions.! spin symmetic (or nspin=1, not 2 channell binded and so on...
subroutine getqmode()  !no output. getq mode just output. Not return variables.
  use m_lmfinit, only: nspec,ispec,nbas,stdo,lmxax
  use m_ext,only:sname   ! read veswavatm.* and qbyl.*
  implicit none
  logical:: debug,cmdopt0
  integer:: lmxa,is,ifiwv,il,isp,ifiqb,ib,ibas,idummy,npri
  real(8):: qsetsum,qlx(0:100),evll,qrmtx
  real(8),allocatable::qrmt(:,:),ql(:,:),qset(:,:),qatot(:)
  debug    = cmdopt0('--debugbndfp')
  write(stdo,"(a)") 'getqmode(): Q from ql given by lmf'
  write(stdo,"(a)") 'WARN current version is only for spin symmetric; and not more than 2(2l+1) occupancy'
  lmxa = lmxax 
  allocate(qrmt(0:lmxa,nspec),qset(0:lmxa,nspec),qatot(nspec))
  qrmt= -1d-10
  open(newunit=ifiwv,file='veswavatm.'//trim(sname))
  do is=1,nspec
     read(ifiwv,*) qatot(is)
     print *,'qatom tot=',qatot(is)
     read(ifiwv,*)
     do
        read(ifiwv,*)isp,il,npri,evll,qrmtx
        if(isp==9) exit     !find terminator line for each atom in veswavatm file
        if( isp > 1) cycle
        if(evll<0d0) then
           if (qrmt(il,is) >1d-2) then
              call rx('not yet implemented when two per l (eg. 3d and 4d) are bounded')
           else
              qrmt(il,is)=qrmtx
           endif
        endif
     enddo
     do il=0,lmxa
        if(qrmt(il,is)>0d0) print *, is,il,qrmt(il,is)
     enddo
  enddo
  close(ifiwv)
  open(newunit=ifiqb,file='qbyl.'//trim(sname))! WARN. not averaged yet for given specs...
  read(ifiqb,*)
  allocate(ql(0:lmxa,nspec))
  do ib=1,nbas
     read(ifiqb,*)ibas,isp,idummy,qlx(0:lmxa)
     is = ispec(ib)
     ql(0:lmxa,is)=qlx(0:lmxa)
     write(6,"(' ibas ql=',i3,20f12.6)") ib,ql(0:lmxa,is)
  enddo
  do is=1,nspec
     qsetsum=0d0
     do il=0,lmxa
        if(qrmt(il,is)<0d0) then
           qset(il,is)=0d0
        else
           qset(il,is)= ql(il,is)/qrmt(il,is)
           qsetsum = qsetsum+qset(il,is)
           if(debug) print *,'ddddddddd ', is,il ,ql(il,is),qrmt(il,is),qset(il,is)
        endif
     enddo
     write(*,"('not renormalized ',i3,' Q=',12f12.4)")is, qset(0:lmxa,is)
     qset = qatot(is)* qset/qsetsum
     write(*,"('    renormalized     Q=',12f12.4)") qset(0:lmxa,is)
  enddo
  close(ifiqb)
end subroutine getqmode
