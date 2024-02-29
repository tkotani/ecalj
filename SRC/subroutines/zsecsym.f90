module m_zsecsym
   use m_nvfortran,only:findloc
   private
   contains
   subroutine zsecsym(zsec,ntq,nq,nband,nbandmx,nspinmx,nspin, eibzsym,ngrp,tiii,q,is)!Symmetrize zsec for eibz4sig mode. Not recently checked!
   use m_readeigen,only: READEVAL
   implicit none
   complex(8),intent(inout)::zsec(ntq,ntq,nq)
  integer,intent(in)::ntq,nq,nspinmx,nband,nbandmx(nq,nspinmx),is,nspin
  integer,intent(in):: ngrp,eibzsym(ngrp,-1:1,nq)
  logical,intent(in):: tiii !time reversal switch
  real(8),intent(in):: q(3,nq)
  complex(8),allocatable::zsect(:,:)
  integer:: iqq
  integer:: procid,nrankv,ifvxc_,ifevec_,ifiproc,iqqxx, &
       isp,ixx,ixxx,nqixx,nspxx,ispxx,iqbz,i,igrp,iq
  character*256:: extn,ext
  character*256,allocatable:: extp(:)
  integer:: nsym,nhdim,it,nblk,iband,napw,ldim,ierr,ispx,nbsize,nbsizemx &
       ,iblk1,iblk2,ii1,ii2,ie1,ie2,ne1,ne2,iqxx, ndimhx, nspx
  integer,allocatable::iblki(:),iblke(:)
  complex(8),allocatable:: evec(:,:),evec_inv(:,:),evecrot(:,:),rmatjj(:,:,:)
  real(8),allocatable::evaliq(:)
  real(8)::tolry=1d-4,qqqx(3),qtarget(3),tolq=1d-8
  complex(8),allocatable:: ovl(:,:)
  integer::nev,j,nxx
  character(8):: xt
  write(6,*)'zsecsym:'
  allocate( zsect(ntq,ntq), evaliq(nband),iblki(nband),iblke(nband))
  iqq=0                     !iqq is to read multiple vxc.* evec.*
  do 3020 iq=1,nq           !nq means iq for which we will calculate sigma
     iqq=iqq+1
     do 3030 ispx=1,nspin !nspinmx  bugfix!!! at 2021-11-01 !ispx loop is to find isx=is
        if(ispx==is) then
           open(newunit=ifevec_,  file='evec'//trim(xt(iq))//trim(xt(is)),form='unformatted')
           read(ifevec_) nhdim, ldim
           allocate( evec(nhdim,nhdim),evecrot(nhdim,nhdim))
           read(ifevec_) qqqx(1:3), evec(1:nhdim,1:nhdim),nev !nev number of true bands nov2015
           zsect = 0d0
        else                  !skip isx/=is. Need to get access sequential files evec and v_xc.
           cycle
        endif
        do i=1,nq           !nq     !qqqx from evec v_xc.
           if(sum(abs(qqqx-q(:,i)))<tolq) then
              iqxx=i
              goto 3011
           endif
        enddo
        deallocate(evec,evecrot)
        call rx( 'hsfp0_sc: bug:qqqx can not find ...')
3011    continue
        if(tiii) call rx( 'timereversal is not yet implemented')
        !! evec_inv(ib1,iww)= \sum_ib2  ovlinv(ib1,ib2)*dconjg(evec(iww,ib2))  nov2015, we introduce nev. iww is for PMT basis. ib for band index.
        !! This is for converting rotated evec (=evecrot(ib)) in the representation of original evec(ib).
        allocate(ovl(nev,nev))
        do i=1,nev
           do j=1,nev
              ovl(i,j)=sum(dconjg(evec(:,i))*evec(:,j))
           enddo
        enddo
        call matcinv(nev,ovl) !ovl --> ovlinv
        allocate(evec_inv(nev,nhdim))
        evec_inv = matmul(ovl(1:nev,1:nev),dconjg(transpose(evec(:,1:nev)))) !note ovl means ovlinv
        deallocate(ovl)
        evaliq = READEVAL(q(:,iqxx), is)
        nsym = sum(eibzsym(:,:,iqxx))
        itloop:do it=1,1             !no-time reversal yet !it=1,-1,-2 !c.f. x0kf_v4h
           igrploop: do igrp=1,ngrp      !A-rotator
              if( eibzsym(igrp,it,iqxx)==0) cycle
              nblk=0
              iblki=0
              iblke=0
              iblki(1)=1
              !!  degeneracy divider for evaliq. See How to apply EIBZ to
              !! Is this procedure really make speed up so much?
              tolry= 0.2d0      !Degeneracy tol. if tolry is large,
              !! larger tolry is safer, although a little inefficient.
              !! If tolry is too small to divide degenerated values to different blocks --> then we have wrong results.
              ! NOTE that Hamiltonian can be not so symmetric in some reasons)
              nbsizemx=0
              do iband=2,nbandmx(iqxx,is)
                 ! nbandmx is the number of bands for which we calculate self-energy.
                 ! We assume nbandmx(iqxx,is) is well separated for degeneracy.
                 if(evaliq(iband)>evaliq(iband-1)+tolry .OR. iband==nbandmx(iqxx,is)) then
                    nblk=nblk+1
                    if(nblk>=2) iblki(nblk)=iblke(nblk-1)+1
                    if(iband==nbandmx(iqxx,is)) then
                       iblke(nblk)=iband
                    else
                       iblke(nblk)=iband-1
                    endif
                    nbsize = iblke(nblk)- iblki(nblk)+1
                    if( nbsize>nbsizemx ) nbsizemx = nbsize
                 endif
              enddo             ! iband
              !! rotation of evec. Generate evecrot. (Within degenerated block, evec are mapped).e
              allocate(rmatjj(nbsizemx,nbsizemx,nblk))
              napw=nhdim-ldim
              do iblk1=1,nblk
                 ii1=iblki(iblk1)
                 ie1=iblke(iblk1)
                 ne1=ie1-ii1+1
                 call rotwvigg2(igrp,q(:,iqxx),q(:,iqxx),nhdim,napw,ne1,evec(:,ii1:ie1),evecrot(:,ii1:ie1),ierr )
                 rmatjj(1:ne1,1:ne1,iblk1) = matmul(evec_inv(ii1:ie1,:),evecrot(:,ii1:ie1))
              enddo
              do iblk1=1,nblk
                 do iblk2=1,nblk
                    ii1=iblki(iblk1)
                    ie1=iblke(iblk1)
                    ne1=ie1-ii1+1
                    ii2=iblki(iblk2)
                    ie2=iblke(iblk2)
                    ne2=ie2-ii2+1
                    zsect(ii1:ie1,ii2:ie2)= zsect(ii1:ie1,ii2:ie2) + matmul(  dconjg(transpose(rmatjj(1:ne1,1:ne1,iblk1))), &
                         matmul(zsec(ii1:ie1,ii2:ie2,iqxx),rmatjj(1:ne2,1:ne2,iblk2))  )
                 enddo    
              enddo       
              deallocate(rmatjj)
           enddo igrploop
        enddo itloop
        deallocate(evec, evec_inv, evecrot)
        zsec(:,:,iqxx) = zsect(:,:)/dble(nsym)
        close(ifevec_)
3030 enddo
3020 enddo
deallocate(iblki,iblke,evaliq)
deallocate(zsect)
end subroutine zsecsym 
subroutine rotwvigg2(igg,q,qtarget,ndimh,napw_in,nband,evec,evecout,ierr) !wave funciton rotation. This is originally rotwvigg. 2023-9-13 copied from rotwv.f90 so as to modify origianl rotwvigg for another purpose. no shared I/O
   use m_hamindex,only: symops,invgx,miat,tiat,shtvg,qlat,plat,dlmm,ngrp,norbmto,ibastab,ltab,ktab,offl,offlrev
   use m_hamindex,only: igv2,igv2rev,napwk,nbas,pwmode
   use m_ftox
   implicit none
   intent(in)::        igg,q,qtarget,ndimh,napw_in,nband,evec
   intent(out)::                                              evecout,ierr
   !! ==  wave function rotator by space group operation.
   !! OUTPUT evecout, ierr
   !! NOTE:
   !! rotation of coefficients on PMT basis.
   !!  phi(r) = \sum_i evec(i,iband) |F_i> ==> Rotated[phi](r)=\sum_i evecout(i,iband) |F_i>  by sym(:,:,ig).
   !!  Rotated[phi](r)= phi[sym^-1(r)], where   sym(r)=r'= symops*r + shftvg.
   integer::ig,ndimh,napw_in,nband,ibaso,iorb,nnn(3),igx,init1,init2,iend1,iend2,nlmto,ierr &
        ,igg,ikt2,ikt,l,ibas,ig2,k
   real(8)::q(3),gout(3),delta(3),ddd(3),qpg(3),platt(3,3),qtarget(3),qx(3),det,qpgr(3),ddd2(3)
   complex(8):: evec(ndimh,nband),evecout(ndimh,nband),phase(nbas)
   real(8),parameter:: tolq=1d-4
   complex(8),parameter:: img=(0d0,1d0), img2pi=2*4d0*datan(1d0)*img
   platt = transpose(plat) !this is inverse of qlat
   ierr=1
   !! check q is really rotated to qtarget by symops(:,:,igg)
   call rangedq( matmul(platt,(qtarget-matmul(symops(:,:,igg),q)) ), qx)
   if(sum(abs(qx))>tolq) then
      write(6,"(a,3f7.3,2x,3f7.3)")'  rotwvigg: qtarget is not a star of q',q,qtarget
      call rx( 'rotwvigg: qtarget is not symops(:,:,ig)*q')
   endif
   evecout = 0d0
   nlmto = ndimh-napw_in
   !! mto part
   if(nlmto/=0) then
      phase = [(exp(-img2pi*sum(qtarget*tiat(:,ibas,igg))),ibas=1,nbas)]
      do iorb=1,norbmto !orbital-blocks are specified by ibas, l, and k.
         ibas = ibastab(iorb)
         l   = ltab(iorb)
         k   = ktab(iorb)
         init1 = offl(iorb)+1
         iend1 = offl(iorb)+2*l+1
         init2 = offlrev(miat(ibas,igg),l,k)+1
         iend2 = offlrev(miat(ibas,igg),l,k)+2*l+1
         evecout(init2:iend2,:)= matmul(dlmm(-l:l,-l:l,l,igg),evec(init1:iend1,:))*phase(ibas)
      enddo
   endif
   !! apw part
   if(napw_in/=0) then
      ikt  = getikt(q)       !index for q
      ikt2 = getikt(qtarget) !index for qtarget
      if(napw_in /= napwk(ikt) ) then
         call rx('rotwv: napw_in /= napw(ikt)')
      endif
      do ig = 1,napw_in
         qpg  = q + matmul( qlat(:,:),igv2(:,ig,ikt))  !q+G
         qpgr = matmul(symops(:,:,igg),qpg)            !rotated q+G
         nnn= nint(matmul(platt,qpgr-qtarget)) !integer representation of G= qpgr - qtarget
         ig2 = igv2rev(nnn(1),nnn(2),nnn(3),ikt2) !get index of G
         if(ig2>=999999) then
            block
              integer:: i1
            do i1=1,napwk(ikt)
               write(6,ftox)'yyy0 igv2', ftof(q,3),     ikt,i1, ' ',igv2(:,i1,ikt)
            enddo
            do i1=1,napwk(ikt2)
               write(6,ftox)'yyy1 igv2', ftof(qtarget,3),ikt2,i1,' ',igv2(:,i1,ikt2)
            enddo
            endblock
            write(6,ftox)'rotwvigg: q=',ftof( q,3),'qtarget=', ftof(qtarget,3)
            write(6,ftox)'rotwvigg  qr=',ftof(matmul(symops(:,:,igg),q),3)
            write(6,ftox)'rotwvigg: qpg=',ftof(qpg,3),'qpgr=', ftof(qpgr,3)
            write(6,ftox)'rorwvigg: igv2rev ikt2=',nnn(1),nnn(2),nnn(3)
            call rx('rotwvigg can not find index of mapped G vector ig2')
         endif
         evecout(nlmto+ig2,:)= evec(nlmto+ig,:) * exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
      enddo
   endif
   ierr=0
 end subroutine rotwvigg2
 integer function getikt(qin) !return !> get index ikt such that for qin(:)=qq(:,ikt)
   use m_hamindex,only: qtt,nqtt
   intent(in)::          qin
   integer::i
   real(8):: qin(3)
   getikt  = findloc([(sum(abs(qin-qtt(:,i)))<1d-8,i=1,nqtt)],value=.true.,dim=1)  !=index for q
   if(getikt<=0) call rx('getikt zsecsym can not find ikt for given q')
endfunction getikt
endmodule m_zsecsym
