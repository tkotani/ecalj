module m_rdsigm2
  use m_nvfortran,only:findloc
  public:: m_rdsigm2_init, Getsenex, Dsene, senex !getsenex returns the self-energy term.
  public:: ndimsig, sene
  integer,protected,public :: nk1,nk2,nk3
  complex(8),pointer,public:: hrr(:,:,:,:,:,:) ! self-energy on real-space mesh points.
  private
  complex(8),allocatable,target :: sfz(:,:,:,:,:,:)
  complex(8),allocatable,protected:: sene(:,:),senex(:,:)
  complex(8), allocatable,protected :: hhrs (:,:,:,:)
  integer, allocatable,protected :: iaxs(:,:)
  integer, allocatable,protected :: ntabs(:)
  integer,protected:: ndimsig,nspsigm,ndhrs,ham_nqsig
  real(8),protected,allocatable ::  rv_p_oqsig (:)
  logical,private:: debugmode=.true.
contains
  subroutine getsenex(qp,isp,ndimh,ovlm)! Return self-energy senex at qp,isp
    implicit none
    integer:: isp,ndimh,ispsigm
    real(8):: qp(3)
    complex(8),allocatable:: ovlmtoi(:,:),ovliovl(:,:)
    complex(8):: ovlm(ndimh,ndimh)
    call tcn('getsenex')
    allocate(sene(ndimsig,ndimsig))
    ispsigm=isp
    if(isp>nspsigm) ispsigm = nspsigm
    call bloch2(qp,ispsigm,sene) !return self-energy for given qp,ispsigm
    allocate( ovlmtoi(ndimsig,ndimsig),ovliovl(ndimsig,ndimh))
    ovlmtoi = ovlm(1:ndimsig,1:ndimsig)
    call matcinv(ndimsig,ovlmtoi)
    ovliovl = matmul(ovlmtoi,ovlm(1:ndimsig,1:ndimh))
    deallocate(ovlmtoi)
    allocate(senex(ndimh,ndimh))
    senex = matmul(transpose(dconjg(ovliovl)), matmul(sene,ovliovl))
    deallocate(ovliovl)
    call tcx('getsenex')
  end subroutine getsenex
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine dsene()
    deallocate(senex,sene)
  end subroutine dsene
  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine m_rdsigm2_init() !get hrr (fourier transformation of self-energy in full BZ)
    use m_lmfinit,only : nbas,pwmode=>ham_pwmode,ldim=>nlmto
    use m_lgunit,only:stdo
    use m_mksym,only: laf=>AFmode !symops_af!,napwmx
    use m_MPItk,only: procid,master
    use m_ext,only:sname
    use m_ftox
    integer:: ierr,ifi,ndimh_dummy,ifis2,ik1,ik2,ik3,is,iset,nqp
    logical:: mtosigmaonly,cmdopt0
    character strn*120
    real(8),allocatable:: qsmesh2(:,:,:,:)
    call tcn('m_rdsigm2_init')
    ndimsig= ldim      
    if(procid==master) then
       open(newunit=ifi,file='sigm.'//trim(sname),form='unformatted')
       read(ifi,err=9995,end=9995) nspsigm,ndimh_dummy,nk1,nk2,nk3,nqp
!       write(stdo,ftox) 'rrrrreading sigm nspsigm,nk1,nk2,nk3,nqp=', nspsigm,nk1,nk2,nk3,nqp
       write(stdo,"(' sigm file has ',i5,' irreducible QP: nk =',3i5)") nqp,nk1,nk2,nk3
!       laf=allocated(symops_af) !jun2015takao !reserved for future
       if(laf) nspsigm=2      ! we need recover laf mode
       allocate( qsmesh2(3,nk1,nk2,nk3) )
       if(mtosigmaonly()) then
          allocate(sfz(nk1,nk2,nk3,ndimsig,ndimsig,nspsigm))
       else
          call rx('mtosigmaonly()=T is needed in the current version sep2012')
       endif
       rewind ifi
       call rdsigm2(nspsigm,ifi,nk1,nk2,nk3,ldim,qsmesh2,mtosigmaonly(),ndimsig) ! Get self-energy sfz in the full BZ.
       close(ifi)
       WRITEsig_fbz: if(cmdopt0('--wsig_fbz')) then
          open(newunit=ifis2,file='sigm_fbz.'//trim(sname),form='unformatted')
          write(stdo,"(a)")' Writing sigm_fbz.* for SYMGRP e --wsig_fbz'
          write(ifis2) nspsigm,ndimsig,nk1,nk2,nk3,nk1*nk2*nk3,0,0,0
          do is=1,nspsigm
             do ik1=1,nk1
                do ik2=1,nk2
                   do ik3=1,nk3
                      write(ifis2) qsmesh2(1:3,ik1,ik2,ik3),is
                      write(ifis2) sfz(ik1,ik2,ik3,1:ndimsig,1:ndimsig,is)
                   enddo
                enddo
             enddo
          enddo
          close(ifis2) !call fclose(ifis2)
       endif WRITEsig_fbz
       call fftz3(sfz,nk1,nk2,nk3,nk1,nk2,nk3,ndimsig**2*nspsigm,iset,+1) !+1 backward ! FT hrr is on regular mesh points. For bloch2
       hrr => sfz ! rename sfz as hrr
    endif                  !procid==master
    if(cmdopt0('--wsig_fbz')) call rx0('end of --wsig_fbz mode')
    call mpibc1_int(nspsigm,1,'senebr:nspsigm')
    call mpibc1_int(nk1,1,'senebr:nk1')
    call mpibc1_int(nk2,1,'senebr:nk2')
    call mpibc1_int(nk3,1,'senebr:nk3')
    if(procid/=master) allocate(sfz(nk1,nk2,nk3,ndimsig,ndimsig,nspsigm))
    call mpibc1_complex(sfz, size(sfz),'senebroadcase:sfz')
    if(procid/=master) hrr => sfz    ! rename sfz as hrr
    if(procid==master) write(stdo,*)' --- end of sigmainit: rdsigm2 initialization section ---'
    call tcx('m_rdsigm2_init')
    return
9995 continue
    call rx('sigmainit: readin error of sigm file')
  end subroutine m_rdsigm2_init
  subroutine rdsigm2(nspsigm,ifis, nk1,nk2,nk3,ldim,qsmesh, mtosigmaonly,ndimsig) ! Expand self-energy (read by ifis) to all the q point on mesh.
    use m_mksym,only: symops,ngrp,ngrpAF,laf=>AFmode
    use m_hamindex,only: napwk
    use m_lgunit,only:stdo
    use m_lattic,only: plat=>lat_plat
    use m_ftox
    !! nbas is in this structure
    !! input
    !!    ifis:  file hundle for self-energy file sigm. only at irreducible q point.
    !! output
    !!   complex(8)::sfz(nk1,nk2,nk3,ndimsig,ndimsig,nsp):  self-energy (\Sigma-Vxc) for  all the q points on mesh.
    !!   real(8):: qsmesh(3,nk1,nk2,nk3)
    !!  Self-energy (\Sigma-Vxc) is read from ifis file.
    !!  It is stored into sigm_zv(ndimsig_r,ndimsig_r), which is rotated to be sfz in the full BZ by hamfb3k.
    !!  ndimsig<=ndimsig_r
    !!   * ndimsig<=ndimsig_r is because I expect compatibility with current hqpe_sc where ndimsig_r= nlmto+max(napw)
    !!     This should be corrected near future (written in 20sep2012).
    !!
    !!  We have to clean up this routine. The purpose of this routine is "read sigm file and expand it in full BZ".
    !!  Not do more than that. (in future, we do scaling of simga in bndfp.F.
    !!  Many un-used local variables are contained.
    !!  Especially qsmesh (regular q mesh for self-energy.) is very problematic. It should be given at a place, and then
    !!  it should be used somewhere else.
    implicit none
    integer:: ifis,ndimsig_r,lwsig,i,j,ifis2,ifiz,isp,nspsigm,nglob,lrsig, nkxyz(3),nk1,nk2,nk3,nsgrp,nsgrps,mxkp,nqp,nqps, &
         j1,iq1,nqsig,lrot,iprint,lssym,ndims,ndimz,iq,n123(4),lcore,lhigh,ohrss,osigm2,odelt,oistb2
    logical :: llshft(3),lphase,lsplts,lnwmsh, latvec,lfbzin,lfbzout
    character outs*80,out2*80,dc*1,rots*120
    integer,parameter::niax=10
    integer ,allocatable,target :: ipq(:,:,:)
    real(8) ,allocatable :: qp_rv(:,:), wgt_rv(:),evls(:),evlz(:),sigii(:)
    complex(8),allocatable :: wk_zv(:), sigm_zv(:,:), siglda(:,:),z(:,:),sigo(:,:)
    integer :: is(3),lshft(3),ifac(3), jj1,jj2,jj3,k,iwdummy,iaf
    real(8):: rb(3,3),qb(3,3), qp(3),tolq=1d-4,rsstol,rotm(3,3), qsmesh(3,nk1,nk2,nk3) 
    integer:: i1,i2,i3,ikt,ldim,napw_in,debugmode, ndimsig,ix
    logical:: isanrg, debug=.false.,mtosigmaonly
    real(8):: qir(3),diffq(3),platt(3,3)
    integer:: ii1,ii2,ii4,ispr,ig,nsp_,ndimh_,nk1_,nk2_,nk3_,nqp_
    character(300)::aaa
    character(8):: xt
    integer,allocatable,target:: ipqaf(:,:,:)
    integer,pointer:: ipq_pointer(:,:,:)
    call tcn('rdsigm2')
    write(stdo,*)'rdsigm2:'
    sfz=1d99
    lshft=0
    rewind ifis !Read sigma(orbital basis) from a file 
    read(ifis) nsp_,ndimsig_r,nk1_,nk2_,nk3_,nqp_
    nsgrp  = ngrp !lat_nsgrp
    nsgrps = nsgrp
!    print *,' lat_nsgrp=',lat_nsgrp
    llshft = .false. 
!    call bzmsh0(plat,llshft,nk1,nk2,nk3,is,ifac,rb,qb)!Get is,ifac,qb,qlat,qoff  Get list of irreducible k-points, and ipq and gstar arrays
    mxkp = nk1*nk2*nk3
    allocate(rv_p_oqsig(abs(3*mxkp)),qp_rv(3,mxkp),ipq(nk1,nk2,nk3))
    allocate(wgt_rv(mxkp))
    write(stdo,"(a)")'  points in Full BZ where sigma should be given ...'
    call bzmesh(plat,qb,ifac,nk1,nk2,nk3,llshft,iwdummy,0, ipq,rv_p_oqsig, wgt_rv, nqsig, mxkp)
    ham_nqsig=nqsig
    write(stdo,"(a)")'  Irreducible qp for which sigma is calculated ...'
    call bzmesh(plat,qb,ifac,nk1,nk2,nk3,llshft,symops, nsgrps,ipq, qp_rv,wgt_rv,nqps,mxkp)
    if(nqps/=nqp_) call rx('nqps/=nqp_ from sigm '//xt(nqps)//' '//xt(nqp))
    platt=transpose(plat)
    do concurrent(i1=1:nk1,i2=1:nk2,i3=1:nk3)
       qsmesh(:,i1,i2,i3) = matmul(qb(:,:), [(i1*ifac(1)-1), (i2*ifac(2)-1), (i3*ifac(3)-1)])
    enddo
    if(laf) then
       if(iprint()>10) write(stdo,*)'rdsimg2: AF mode, mapping from irr points to regular mesh point'
!       write(6,*)'ngrp ngrpAF=',ngrp,ngrpAF
!       do ig= ngrp+1,ngrp+ngrpAF !only AF symmetry (equivalent with symops_af)
!          write(stdo,ftox)ig,ftof(reshape(symops(:,:,ig),[9]),3)
!       enddo
       allocate(ipqaf(nk1,nk2,nk3))
       ipqaf=0
       do i1=1,nk1
          do i2=1,nk2
             do i3=1,nk3
                do 1111 iq1=1,nqps
                   qir = qp_rv(:,iq1)
                   do ig= ngrp+1,ngrp+ngrpAF !only AF symmetry (equivalent with symops_af)
                      call rangedq( matmul(platt,(qsmesh(:,i1,i2,i3) - matmul(symops(:,:,ig),qir))), diffq)
                      if(sum(abs(diffq))<tolq) then
                         ipqaf(i1,i2,i3) = iq1    !iq1 is pointer to the irreducible q point = qp_rv(:,iq1)
                         goto 1122
                      endif
                   enddo
1111            enddo
                write(aaa,"(3i5,3f13.5)") i1,i2,i3,qsmesh(:,i1,i2,i3)
                call rx('rdsigm2: 1111 loop can not find ipqaf'//trim(aaa))
1122            continue
             enddo
          enddo
       enddo
    endif
    ispsigmloop: do 2001 isp = 1, nspsigm !Generate hrs = sigma(T) from file sigma(k) ---
       if (isp==2 .AND. nsp_==1 ) then !If sigm file not spin polarized, use sigm from spin 1
          rewind ifis
          read(ifis)
       endif
       allocate(wk_zv(ndimsig_r**2),sigm_zv(ndimsig_r,ndimsig_r))
       do iq1 = 1, nqps ! look for a tag qp in sigm, where qp=qp_rv(:,iq1)  for given iq1
          do ix=0,1
             do
                read(ifis,end=468) qp, ispr ! ispr is added dec2013
                read(ifis) sigm_zv
                if(laf) then
                   if(ispr==2) cycle
                else
                   if(ispr/=isp) cycle
                endif
                if(sum(abs(matmul(transpose(plat),qp-qp_rv(:,iq1))))<tolq) goto 460
             enddo
468          continue
             rewind ifis
             read(ifis)
          enddo
          write(aaa,"(i5,3f13.5)") iq1, qp_rv(:,iq1)
          call rx(' rdsigm2: read error. In sigm, we did not find iq= '//trim(aaa))
460       continue
          if(mtosigmaonly .OR. ldim==ndimsig_r) then
             napw_in=0
             ikt=-9999
          else
             ikt = getikt(qp)
             napw_in= napwk(ikt)
          endif
          iaf=0
          ipq_pointer => ipq
          if(laf .AND. isp==1) then
             iaf=1
          elseif(laf .AND. isp==2) then
             iaf=2
             ipq_pointer => ipqaf
          endif
          write(stdo,"(a,2i5,' ',13f13.5)")' rdsigm2:Goto hamfb3k  xxx input isp,iaf,qp=', isp,iaf,qp
          if(iprint()>60) write(stdo,"(a,13f13.5)")' rdsigm2:Goto hamfb3k  xxx input qp=', qp
          call hamfb3k(qp, iq1, nk1,nk2,nk3,ipq_pointer, napw_in,ndimsig,ndimsig,ndimsig,qb,ldim,&
               ifac,sigm_zv(1:ndimsig,1:ndimsig),iaf, sfz(1,1,1,1,1,isp))
          if(debugmode>0) write(stdo,"(a,3f13.5)")'end of hamfbk3'
       enddo
       deallocate(sigm_zv,wk_zv)
2001 enddo ispsigmloop
    call tcx('rdsigm2')
    if (allocated(wgt_rv)) deallocate(wgt_rv)
    if (allocated(ipq)) deallocate(ipq)
    if (allocated(qp_rv)) deallocate(qp_rv)
!    if (allocated(gstar_iv)) deallocate(gstar_iv)
  end subroutine rdsigm2
  ! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine hamfb3k(qin,iq,nk1,nk2,nk3,ipq,napw_in,ndimh,ldima,ldimb,qb,ldim,ifac, hq,iaf, gfbz)! Generate gfbz (full BZ) from hq (1st BZ).
    use m_lgunit,only:stdo
    !i q iq: q and index for q
    !i   nk1,nk2,nk3:  no. divisions for the qp in 3 recip. latt. vecs
    !i   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
    !i          mesh point (i1,i2,i3) is mapped (bzmesh.f)
    !i   ndimh: dim of self energy hq at qin
    !i   ldima :dimensions gfbz; also the number of rows in gfbz to fill.
    !i         :usually dimension of lower (or l+i) block for crystal
    !i   ldimb :dimensions gfbz; also the number of columns in gfbz to fill
    !i         :usually dimension of lower (or l+i) block for crystal
    !i   qb    :vectors of a microcell in the Brillouin zone
    !i   hq    : Self-energy (or non-local potential for this iq
    !i read date in m_hamindex .
    !o Outputs
    !o   gfbz  : For those qp in star iq, hq stored
    ! ----------------------------------------------------------------------
    implicit none
    integer:: nk1,nk2,nk3,ipq(*),ndimh,ldima,ldimb,napw_in,debugmode
    real(8)::    qin(3),qb(3,3) !,plat(3,3),qlat(3,3)
    complex(8):: hq(ndimh,ndimh),gfbz(nk1,nk2,nk3,ldima,ldimb)
    integer:: i,i1,i2,i3,ig,iq,iq1,is,j,jj1,jj2,jj3,js,k,ierr,ifac(3),j1,j2,ik1,ik2,ik3,isp,ldim,iaf
    real(8):: q1(3),qk
    character(200)::aaa
    ! Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
    qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) + (jj2*ifac(2)-1)*qb(k,2) + (jj3*ifac(3)-1)*qb(k,3)
    !qk(k,jj1,jj2,jj3) = sum(qb(k,:)*[(jj1*ifac(1)-1),(jj2*ifac(2)-1), (jj3*ifac(3)-1)]) <=Not working in ifort
    call tcn('hamfb3k')
    if(debugmode>0) print *, 'hamfb3k: start...'
    iq1 = 0
    do  i3 = 1, nk3
       do  i2 = 1, nk2
          do  i1 = 1, nk1
             iq1 = iq1+1
             if(debugmode>0) print *,'iq iq1 ipq(iq1)',iq,iq1,ipq(iq1)
             !! ipq(iq1) ist gives a table to point irreducible point.
             !!          q1,iq1 is target on regular mesh <--- qin,iq is irreducible points; this mapping is by rotsig.
             !!          iq=ipq(iq1) shows iq for given iq1.
             if (ipq(iq1) /= iq) cycle !this must make things efficient
             q1(1) = qk(1,i1,i2,i3)
             q1(2) = qk(2,i1,i2,i3)
             q1(3) = qk(3,i1,i2,i3)
             if(debugmode>0) print *,q1
             if(debugmode>0) write(stdo,"(a,3f13.5)")' input          qin = ', qin
             if(debugmode>0) write(stdo,"(a,3f13.5)")' target a        q1 = ', q1 ! qin = g(:,:,ig)^{-1}*q1
             call rotsig(qin,q1,ndimh,napw_in,ldim,hq,gfbz(i1,i2,i3,:,:),ierr,iaf)
             if(ierr/=0) write(aaa,"(' qin=',3f13.6,' q1=',3f13.6)") qin,q1
             if(ierr/=0) call rx('hamfb3: rotsig do not map qin to q1;'//trim(aaa))
          enddo
       enddo
    enddo
    if(debugmode>0) print *, 'hamfb3k: end...'
    call tcx('hamfb3k')
  end subroutine hamfb3k
  
  subroutine rotsig(qin,qout,ndimh,napw_in,ldim,sig, sigout,ierr,iaf) !Sigm rotator. sig at qin to sig at qout. ===
    use m_lmfinit,only: ibastab,ltab,ktab,offl,offlrev,norbmto
    use m_mksym,only:   ngrp,ngrpAF,symops,invgx,miat,tiat,shtvg,dlmm
    use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
    use m_lgunit,only: stdo
    use m_hamindex,only: igv2,igv2rev,napwk
    implicit none
    intent(out)::                                    sigout,ierr
    !! obtain sigout for qout.
    !! a little confusing since qin=symops(qout), and basis rotation. Need to clean up.
    !! Both of q and qtarget should be in qq table(in m_hamindex) which is given by gen_hamindex
    !! (read from QGpsi).
    !! Used idea is   <base(qin)|sigma(qin)|base(qin)> = <g(base)|sigma |g(base)>.
    !! where qin=g(qout).  qtarget=qin= g(q=qout)
    integer:: ig,ndimh,napw_in,ibaso,iorb,nnn(3),igx,init1,init2,iend1,iend2,nlmto,ierr,igg,ikt2,ikt,l,ibas,ig2,k,ix
    real(8):: qin(3),qout(3)
    real(8)   :: q(3),delta(3),ddd(3),qpg(3),platt(3,3),qtarget(3),qx(3),det,qpgr(3),ddd2(3) !plat(3,3),qlat(3,3)
    complex(8):: phase,img=(0d0,1d0),img2pi, sig(ndimh,ndimh),sigout(ndimh,ndimh)
    complex(8),allocatable:: sigx(:,:)
    integer :: ldim,debugmode,iaf,ngini,ngend
    character(300)::aaa
    real(8),parameter:: tolq=1d-4
    call tcn('rotsig')
    img2pi=2*4d0*datan(1d0)*img
    ierr=1
    platt=transpose(plat) !this is inverse of qlat
    if(debugmode>0) write(stdo,"('rotsig: qin qout=',3f9.4,x,3f9.4)") qin,qout
    qtarget= qin
    q      = qout  
    AntiferroMechanism: if(iaf==2) then !AF isp=2. These are antiferro space-group operations
       ngini = ngrp+ 1
       ngend = ngrp+ ngrpAF
    else
       ngini = 1
       ngend = ngrp
    endif AntiferroMechanism
    GetiggOFsymops: do igx=ngini,ngend ! we find qtrget = symops(igx) * q    (this means qin = symops(igc) qout).
       if(debugmode>0) print *, 'ddd=',matmul(platt,(qtarget-matmul(symops(:,:,igx),q)))
       call rangedq(   matmul(platt,(qtarget-matmul(symops(:,:,igx),q))), qx)
       if(sum(abs(qx))<tolq) then
          igg=igx
          if(debugmode>0) then
             print *,'ttt: q      =',q
             print *,'ttt: qtarget=',qtarget
             print *,'ttt: matmul q =',matmul(symops(:,:,igx),q)
             print *,'ttt: rotsig: OK! igg=',igg
             print *
          endif
          goto 1012
       endif
    enddo Getiggofsymops
    write(aaa,"(a,3f7.3,2x,3f7.3)")' rotsig: qtarget is not a star of q',q,qtarget
    call rx(trim(aaa))
    print *
    return
1012 continue
    allocate(sigx(ndimh,ndimh))
    sigx=0d0
    nlmto=ldim
    if(debugmode>0) then
       print *,' tttt: invgx =',invgx(igg),shtvg(:,igg)
       print *,' tttt: ntorb napwin',norbmto,ndimh,napw_in,nlmto
    endif
    MTOpart: if(nlmto/=0 )then
       ibaso=-999
       OrbitalBlock: do iorb=1,norbmto !orbital-blocks are specified by ibas, l, and k.
          ! ndex of Hamiltonian is devided into these blocks.
          ibas = ibastab(iorb)
          if(ibas/=ibaso) phase = exp( -img2pi*sum(qtarget*tiat(:,ibas,igg)) )
          ibaso=ibas
          l   = ltab(iorb)
          k   = ktab(iorb)
          init1 = offl(iorb)+ 1
          iend1 = offl(iorb)+ 2*l+1
          init2 = offlrev(miat(ibas,igg),l,k)+ 1
          iend2 = offlrev(miat(ibas,igg),l,k)+ 2*l+1
          do ix=1,ndimh
             sigx(ix,init1:iend1)= matmul(sig(ix,init2:iend2),dlmm(-l:l,-l:l,l,igg))*phase
          enddo
       enddo OrbitalBlock
    endif MTOpart
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  This is just a hint for future.
    APWpart: if(napw_in/=0) then !unused block. We need bloch periodicity for interpolation.
       write(stdo,*) ' Probably OK-->Remove this stop to make this branch effective.' &
            //' But need to confirm two apw sections in this routines.(phase factors) ' &
            //' Idea of this routine: <i|\sigma|j>_qout= <g(i)|\sigma|g(i)>_qin, where qin=g(qout)'
       stop 'abort'
       ikt  = getikt(q)    !index for q
       ikt2 = getikt(qtarget) !index for qtarget
       !write(stdo,*)' rotsig ikt ikt2=',ikt,ikt2
       if(napw_in /= napwk(ikt) ) call rx('rotsig: napw_in /= napw(ikt)')
       do ig = 1,napw_in
          qpg = q + matmul( qlat(:,:),igv2(:,ig,ikt))      !q+G
          qpgr = matmul(symops(:,:,igg),qpg)               !rotated q+G
          nnn=nint(matmul(platt,qpgr-qtarget))
          !print *,ig,'nnn  ikt2=',nnn,ikt2
          ig2 = igv2rev(nnn(1),nnn(2),nnn(3),ikt2)
          phase= exp( -img2pi*sum(qpgr*shtvg(:,igg)) )
          do ix=1,ndimh
             sigx(ix,nlmto+ig) = sig(ix,nlmto+ig2) * phase
          enddo
       enddo
       if(debugmode>0) print *,' apw part end 111: ikt ikt2=',ikt,ikt2
    endif APWpart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    MTOpart2: if(nlmto/=0) then
       ibaso=-999
       do iorb=1,norbmto !orbital-blocks are specified by ibas, l, and k.
          ! ndex of Hamiltonian is devided into these blocks.
          ibas = ibastab(iorb)
          if(ibas/=ibaso) phase = exp( img2pi*sum(qtarget*tiat(:,ibas,igg)) )
          ibaso=ibas
          l   = ltab(iorb)
          k   = ktab(iorb)
          init1 = offl(iorb)+ 1
          iend1 = offl(iorb)+ 2*l+1
          init2 = offlrev(miat(ibas,igg),l,k)+ 1
          iend2 = offlrev(miat(ibas,igg),l,k)+ 2*l+1
          do ix=1,ndimh
             sigout(init1:iend1,ix)= phase * matmul(transpose(dlmm(-l:l,-l:l,l,igg)),sigx(init2:iend2,ix))
          enddo
       enddo
       if(debugmode>0) print *,' end of 2nd mto part q=',q
    endif MTOpart2
    
    APWpart2: if(napw_in/=0) then
       ikt  = getikt(q)    !index for q
       ikt2 = getikt(qtarget) !index for qtarget
       if(debugmode>0) print *,' rotsig 111 ikt ikt2=',ikt,ikt2
       if(napw_in /= napwk(ikt) ) then
          call rx('rotsig: napw_in /= napw(ikt)')
       endif
       do ig = 1,napw_in
          qpg = q + matmul( qlat(:,:),igv2(:,ig,ikt))      !q+G
          qpgr = matmul(symops(:,:,igg),qpg)               !rotated q+G
          nnn=nint( matmul(platt,qpgr-qtarget))
          if(debugmode>0) print *,ig,'nnn  ikt2=',nnn,ikt2
          ig2 = igv2rev(nnn(1),nnn(2),nnn(3),ikt2)
          phase=exp(img2pi*sum(qpgr*shtvg(:,igg)))
          do ix=1,ndimh
             sigout(nlmto+ig,ix) =   sigx(nlmto+ig2,ix) * phase
          enddo
       enddo
       if(debugmode>0) print *,' apw part end 222: ikt ikt2=',ikt,ikt2
    endif APWpart2
    deallocate(sigx)
    call tcx('rotsig')
    ierr=0
  end subroutine rotsig
  
  subroutine bloch2(qp,isp, sll)! Bloch transform of real-space matrix  sll(j,i) =\sum_T hrr(T,i,j,isp,T)*phase, where phase is the avarage for Tbar (shortest equivalent to T).  exp(ikT)
    use m_lmfinit, only: ib_table
    use m_gennlat_sig,only: nlatS,nlatE,nlat
    use m_lattic,only: plat=>lat_plat
    implicit none
    integer:: ik1,ik2,ik3,i,j,ib1,ib2,ix,nS,nE,nnn,isp,iset
    real(8):: qp(3)
    real(8),parameter:: twopi = 8*datan(1d0)
    complex(8):: sll(ndimsig,ndimsig), phase,img=(0d0,1d0)
    call tcn('bloch2')
    sll=0d0
    nnn=nk1*nk2*nk3
    do ik1=1,nk1
       do ik2=1,nk2
          do ik3=1,nk3
             do i=1,ndimsig
                do j=1,ndimsig
                   ib1 = ib_table(i) !atomic-site index in the primitive cell
                   ib2 = ib_table(j)
                   nS=nlatS(ik1-1,ik2-1,ik3-1,ib1,ib2)
                   nE=nlatE(ik1-1,ik2-1,ik3-1,ib1,ib2)
                   phase=0d0
                   do ix=nS,nE !T=(ik1,ik2,ik3) is mapped to nlat(:,ix,ib1,ib2). ix=nS,nE are all shortest.
                      phase=phase+exp(-img*twopi*sum(qp*matmul(plat,nlat(:,ix,ib1,ib2)))) !backward
                   enddo
                   phase= phase/(nE-nS+1) !phase averaged
                   sll(i,j) = sll(i,j) + hrr(ik1,ik2,ik3,i,j,isp)/nnn*phase
                enddo
             enddo
          enddo
       enddo
    enddo
    call tcx('bloch2')
  end subroutine bloch2
  integer function getikt(qin) !return !> get index ikt such that for qin(:)=qq(:,ikt)
    use m_hamindex,only: qtt,nqtt
    intent(in)::          qin
    integer::i
    real(8):: qin(3)
    getikt  = findloc([(sum(abs(qin-qtt(:,i)))<1d-8,i=1,nqtt)],value=.true.,dim=1)  !=index for q
    if(getikt<=0) call rx('getikt rdsigm2 can not find ikt for given q')
  endfunction getikt
end module m_rdsigm2
