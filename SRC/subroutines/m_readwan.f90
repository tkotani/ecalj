module m_readwan
  use m_keyvalue,only: getkeyvalue
  use m_iqindx_wan,only: iqindx2_wan
  implicit none

  public:: Write_qdata, Wan_readeigen, Wan_readeval, Wan_readeval2, Readscr, &
       Checkorb,Checkorb2, Diagwan, Diagwan_tr, Wan_imat, Writehmat, Writeddmat, Read_wandata
  integer,protected,public:: nwf,nsp_w,nqtt_w !! read by read_wandata

  private
  logical:: init=.true.
  logical:: readwan=.false.,lreadhrotr=.false.
  real(8):: alat,plat(3,3),ef
  integer:: natom,nrws,n1,n2,n3
  integer,allocatable:: irws(:),ibaswf(:)
  real(8),allocatable:: pos(:,:),rws(:,:),drws(:)
  complex(8),allocatable:: hrotr(:,:,:,:), evecc(:,:)
!!!!! wannier eigenvalue and eigenvector
  complex(8),allocatable:: evecc_w(:,:,:,:)
  real(8),allocatable:: eval_w(:,:,:)
!!!d-orbital check
  integer,allocatable:: idorb(:), nbasclass_mlwf(:)
  integer:: nclass_mlwf
  logical:: debug=.false.,sw1=.true.

contains

  !---------------------------------------------
  subroutine readhrotr()
    !!--- Required data set --------------------- (copied from htbplot.F)
    !! q(:,1:nq)
    !! ef: fermi energy or VBM
    !! alat: unit for primitive vector, atomic positions (in a.u.)
    !! plat: primitive vector
    !! rcut: real-space cutoff for tb parameters.
    !! nwf:  # of Wannier function
    !! nrws: # of R-R' pair
    !! rws(3,i): R-R' vector, i=1,nrws
    !! irws(i) : degeneracy, i=1,nrws
    !! drws(i) : distance, i=1,nrws
    !! hrotr(nwf,nwf,nrws,is) = Hmn(R) = <0m|H |Rn> (m=1,nwf; n=1,nwf, R=1,nrws)
    !! is: spin index
    !! natom: number of atoms in the primitive cell.
    !! ibaswf(nwf): atomic position for the Wannier orbital.
    !! pos(3,natom): atomic poisition in the cell.
    integer(4):: is, ifh, n1, n2, n3, ifile_handle,nq
    if ( .NOT. readwan) call rx("read_wandata have not yet done.") !! use nsp_w, nwf
    do 6001 is=1,nsp_w
       write (6,"('--- Reading HrotRS  :: isp=',I6)") is
       if(is==1) open(newunit=ifh,file='HrotRS.up',form='unformatted')
       if(is==2) open(newunit=ifh,file='HrotRS.dn',form='unformatted')
       read(ifh)alat,plat,natom
       if (is==1) allocate(pos(3,natom))
       read(ifh)pos
       read(ifh)ef
       read(ifh)nwf,nrws,n1,n2,n3

       write (6,"('nwf,nrws,n1,n2,n3=',I3,I6,3I4)") nwf,nrws,n1,n2,n3
       if (is==1) then
          allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8), &
               drws(n1*n2*n3*8),ibaswf(nwf), hrotr(nwf,nwf,nrws,nsp_w))
          ! eal space Hamiltonian in Wannier funciton basis
       endif
       read(ifh) irws,rws,drws,hrotr(:,:,:,is),ibaswf
       close(ifh)
6001 enddo
    lreadhrotr=.true.          !initialize
  end subroutine readhrotr

  !#############################################################3
  subroutine wan_readeigen(qbz,nqbz,isp)
    intent(in)::             qbz,nqbz,isp
    !! generate eigenvalue list
    integer::nqbz,isp
    real(8)::qbz(:,:)
    complex(8),allocatable:: hrotk(:,:,:),hrotkp(:,:),evec(:,:)
    real(8),allocatable:: eval(:)
    integer::iq
!!! fat band plot test
    integer::iffb,iband
    real(8):: rydberg
    logical::sw=.true.

    !! get HrotRS by file handling
    if ( .NOT. lreadhrotr) call readhrotr()

    if (init) then
       print *,"=== wan_readeigen::"
       allocate(hrotkp(nwf,nwf),evecc(nwf,nwf),eval(nwf))
       if (isp==1) allocate(eval_w(nwf,nqbz,nsp_w), &
            evecc_w(nwf,nwf,nqbz,nsp_w))
       do iq = 1,nqbz
          !     write(6,*)' got get_hrotkp_ws iq =',iq
          call get_hrotkp_ws2(qbz(:,iq),hrotkp,isp)
          call diag_hm2(hrotkp,nwf,eval,evecc)
          eval_w(1:nwf,iq,isp)=eval
          evecc_w(1:nwf,1:nwf,iq,isp)=evecc

          ! ccccccccc eigenvalue check
          !$$$            if(iq==1) open(iwf,file="waneval_check.data")
          !$$$!     if (kx==1 .and. iq > 100 ) open(iwf,file="wan_eval_check.data")
          !$$$            if (qbz(2,iq)==1.0) then
          !$$$               if (qbz(3,iq)==1.0) then
          !$$$!     if (iq > 100) then
          !$$$
          !$$$                  write(iwf,"('q(1) ev_w1(:,kx)',6f9.4)") qbz(1,iq),eval_w(:,iq,1)
          !$$$               endif
          !$$$            endif
          !$$$            if (iq==nqbz) close(iwf)
          !$$$            if (iq==nqbz) call rx("check: end readwan...")

          ! cccccccccc   sumcheck for eigenvector  cccccccccccccccc
          !$$$            if (iq==nqbz .and. is==2)  then
          !$$$!     if (iq==1 .and. is==1) then
          !$$$!     if (iq>50 .and. iq<70) then
          !$$$               print *," iq:",iq
          !$$$               do iwf=1,nwf
          !$$$                  write (6,"('sum evecc(up) evecc(dn) =',2f9.4)")
          !$$$     &                 sum(abs(evecc_w(:,iwf,iq,1)**2)),
          !$$$     &                 sum(abs(evecc_w(:,iwf,iq,2)**2))
          !$$$               enddo
          !$$$            endif
          ! cccccccccccccc  fat band plot (test: Gamma to H) ccccccccccccccccccc
          !$$$            if (qbz(2,iq)/=0.0 .and. qbz(2,iq)/=1.0) cycle
          !$$$            if (qbz(3,iq)/=0.0 .and. qbz(3,iq)/=1.0) cycle
          !$$$            if (qbz(1,iq)>1.05) cycle
          !$$$            write (6,"('SYML G-H  qbz',3f9.4)") qbz(:,iq)
          !$$$
          !$$$            !!! initial operation
          !$$$            if (sw) then
          !$$$               sw=.false.
          !$$$               iffb=ifile_handle()
          !$$$               if(is==1) open(iffb,file="fband.up.tmp")
          !$$$               if(is==2) open(iffb,file="fband.dn.tmp")
          !$$$               write(iffb,"('# ef nwf',f9.4,i4)") ef,nwf
          !$$$               write(iffb,*)
          !$$$            endif
          !$$$            !!! write eigenvalue and eigenvector for wannier
          !$$$            do iband=1,nwf
          !$$$               write(iffb,"(i5,f13.5,' ',f13.6,i5,' ')",ADVANCE='NO')
          !$$$     &          iq,qbz(1,iq),(eval_w(iband,iq,is)-ef)*rydberg(),iband
          !$$$               do iwf=1,nwf
          !$$$                  write(iffb,"(f13.6)",ADVANCE='NO')
          !$$$     &                 (abs(evecc_w(iwf,iband,iq,is)))**2
          !$$$               enddo
          !$$$               write(iffb,*)
          !$$$            enddo
          !$$$
          !$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       enddo
       sw=.true.
       close(iffb)
    endif
    deallocate(hrotkp,evecc,eval)
    if (isp==2) then
       init=.false.           !!!eigenvalue
       print *,"=== wan_readeigen end"
       ! eallocate(pos,rws,drws,hrotr)
    endif
  end subroutine wan_readeigen

  !--------------------------------
  subroutine wan_readeval(q,isp, ev_w, evc_w)
    intent(in)::            q,isp
    intent(out)::                  ev_w, evc_w
    integer :: isp
    real(8) :: q(3), ev_w(nwf)
    complex(8) :: evc_w(nwf,nwf)
    integer(4):: iq,iqindx,i
    real(8):: qu(3), shift_ev
    if(init) call rx('wan_readeval: wan_readeigen should be called')
    call iqindx2_wan(q, iq, qu) !qu is used q. q-qu is a G vector.
    !! shift: ev_w --> ev_w + shift_ev (change exchnage splitting by hand for test)
    call getkeyvalue("GWinput","shift_majority",shift_ev,default=0d0)
    if (shift_ev /= 0d0) write(6,*) "shift_ev [eV]: ",shift_ev
!!! return wannier eigenvalue : ev_w
    ev_w(1:nwf)  = eval_w(1:nwf,iq,isp) + shift_ev/13.605693009
    evc_w(1:nwf,1:nwf) = evecc_w(1:nwf,1:nwf,iq,isp)
  end subroutine wan_readeval

  !--------------------------------
  subroutine wan_readeval2(q,isp, ev_w, evc_w)
    intent(in)::             q,isp
    intent(out)::                   ev_w, evc_w
    integer :: isp
    real(8) :: q(3)
    real(8) :: ev_w(nwf)      !eigenvalue
    complex(8) :: evc_w(nwf,nwf) !eigenfunction
    real(8):: qu(3), shift_ev
    complex(8),allocatable:: hrotk(:,:,:),hrotkp(:,:),evec(:,:)
    real(8),allocatable:: eval(:)
!!! fat band plot test
    real(8):: rydberg
    !! need matrix
    if ( .NOT. lreadhrotr) call readhrotr()
    !! solve Hrotkp for given q (or q+k) --> eigenvalues and eigenfunction
    if ( .NOT. allocated(hrotkp)) then
       allocate(hrotkp(nwf,nwf),evecc(nwf,nwf),eval(nwf))
    endif
    call get_hrotkp_ws2(q ,hrotkp,isp)
    call diag_hm2(hrotkp,nwf,eval,evecc)
    !! shift: ev_w --> ev_w + shift_ev (change exchnage splitting by hand for test)
    call getkeyvalue("GWinput","shift_majority",shift_ev,default=0d0)
    if (shift_ev /= 0d0) write(6,*) "shift_ev [eV]: ",shift_ev
    ev_w(1:nwf) = eval(1:nwf) + shift_ev/13.605693009
    evc_w(1:nwf,1:nwf) = evecc(1:nwf,1:nwf)
    deallocate(hrotkp,evecc,eval)
    ! rite (6,"('qqq:',3E13.5)") q
    ! rite (6,"('eval_w=',9f9.4)") ev_w(1:nwf)
  end subroutine wan_readeval2

  !--------------------------------almpost same get_hrotkp_ws in maxloc3
!!!! hrotr (q) => hrotkp
  subroutine get_hrotkp_ws2(q,hrotkp,is)
    !! see Ref.[2] eq.26
    complex(8) :: hrotkp(nwf,nwf), ci,cikr,ceikr,ctmp
    real(8) :: q(3),pi,rk
    integer::im,in,ir,is
    pi = 4d0* atan(1d0)
    ci = (0d0,1d0)
    hrotkp = (0d0,0d0)
    do ir = 1,nrws
       rk = sum(rws(:,ir)*q(:))
       cikr = ci * 2d0 * pi * rk
       ceikr = exp(cikr) / dble(irws(ir))
       do im = 1,nwf
          do in = 1,nwf
             hrotkp(im,in) = hrotkp(im,in) + &
                  ceikr * hrotr(im,in,ir,is)
          enddo
       enddo
    enddo
  end subroutine get_hrotkp_ws2

  !--------------------------------
  subroutine diag_hm2(zmat,ndim,eval,evecc_)
    intent(in)::        zmat,ndim
    intent(out)::                 eval,evecc_
    integer:: ndim,i,nev,nmx
    complex(8),allocatable :: zmat2(:,:),ovlpc(:,:)
    complex(8):: zmat(ndim,ndim),evecc_(ndim,ndim)
    real(8):: eval(ndim),wk(ndim,11)
    integer :: iwk(ndim)
    allocate(zmat2(ndim,ndim),ovlpc(ndim,ndim))
    nev  = ndim
    nmx  = ndim
    zmat2 = zmat
    ovlpc = (0d0,0d0)
    do i=1,ndim
       ovlpc(i,i) = (1d0,0d0)
    enddo
    evecc_ = (0d0,0d0)
    eval = 0d0
    !      call diagno(ndim,zmat2,ovlpc,wk,iwk,evecc,eval)
    call diagcv(ovlpc,zmat2, evecc_, ndim, eval, nmx, 1d99, nev)
    deallocate(zmat2,ovlpc)
  end subroutine diag_hm2

  !---------------------------------------
  subroutine write_qdata(ginv,nqtt_in,qtt_in)
    intent(in)             ginv,nqtt_in,qtt_in
    integer::ifwqb
    integer::nqtt_in
    real(8)::qtt_in(3,nqtt_in)
    real(8)::ginv(3,3)
    open(newunit=ifwqb,file="wanqbz",form='unformatted')
    write(ifwqb) ginv
    write(ifwqb) nqtt_in
    write(ifwqb) qtt_in
    close(ifwqb)
  end subroutine write_qdata
  !---------------------------------------
  subroutine read_wandata() !nwf,nsp_w,nqtt_w
    integer::ifhamw
    open(newunit=ifhamw,file="wan4chi.d",form='unformatted')
    read(ifhamw) nwf,nsp_w,nqtt_w
    close(ifhamw)
    write(6,*) "nwf,nsp_w,nqtt",nwf,nsp_w,nqtt_w
    readwan=.true. !! initialize
  end subroutine read_wandata
  !---------------------------------------
  subroutine readscr(nwf,scrw_)
    intent(in)::       nwf
    intent(out)::          scrw_
    integer::nwf
    integer::ifscrwv,ifscrv!,ifd,ife,ifa
    integer::ir1,irws1
    character(len=9)::charadummy !dummy
    real(8)::rws1(3),freq,freq2 !dummy
    integer::is,iwf1,iwf2,iwf3,iwf4 !dummy
    integer::iwf,jwf,kwf,lwf,ijwf,klwf
    !      real(8),allocatable::rw_w(:,:,:,:), cw_w(:,:,:,:)
    complex(8),allocatable::scrw4(:,:,:,:),scrv4(:,:,:,:)
    complex(8),allocatable::scrw_(:,:)
    integer::idummy
    logical(8)::ijklmag
!!! set idorb: iwf ---> lorb
    call checkorb(1,nwf,idummy)
    !     allocate(rw_w(nwf,nwf,nwf,nwf),cw_w(nwf,nwf,nwf,nwf))
    allocate(scrw4(nwf,nwf,nwf,nwf), scrv4(nwf,nwf,nwf,nwf))
    allocate(scrw_(nwf*nwf,nwf*nwf));scrw_=0d0
    open(newunit=ifscrwv,file="Screening_W-v.UP",form="formatted") !only up
    open(newunit=ifscrv,file="Coulomb_v.UP",form="formatted") !only up
    !$$$      !!! write direct index (ijkl)
    !$$$      ifd=ifile_handle()
    !$$$      open(ifd,file="ijkl_direct.d",form="formatted")
    !$$$      !!! write exchange index (ijkl)
    !$$$      ife=ifile_handle()
    !$$$      open(ife,file="ijkl_exchange.d",form="formatted")
    !$$$      !!! write all index (i,j,k,l --> ijwf,klwf)
    !$$$      ifa=ifile_handle()
    !$$$      open(ifa,file="ijkl_all.d",form="formatted")
    !$$$      write(ifa,*) "# iwf jwf kwf lwf ijwf klwf"
    !     read Screening W, V
    write (6,*) "readscr: wan_ijkl index is wrriten ijkl_*.d"
    ijwf=0
    do 4001 iwf=1,nwf
       do 4002 jwf=1,nwf
          ijwf=ijwf+1
          klwf=0
          do 4003 kwf=1,nwf
             do 4004 lwf=1,nwf
                klwf=klwf+1
!!! Vare Coulomb (v)
                read(ifscrv,"(A,2i5, 3f12.6, 5i5,2f12.6)") &
                     charadummy,ir1, irws1, rws1 ,is,iwf1,iwf2,iwf3,iwf4 &
                     ,scrv4(iwf1,iwf2,iwf3,iwf4)
!!! Screened Coulomb (W-v)
                read(ifscrwv,"(A,2i5, 3f12.6, 5i5,4f12.6)") &
                     charadummy,ir1, irws1, rws1  &
                     ,is,iwf1,iwf2,iwf3,iwf4,freq, freq2 &
                     ,scrw4(iwf1,iwf2,iwf3,iwf4)
                !$$$         print *,"readscr     :",iwf1,iwf2,iwf3,iwf4
                !$$$         print *,"readscr  W-v:",scrw4(iwf1,iwf2,iwf3,iwf4)
                !$$$         print *,"readscr    v:",scrv4(iwf1,iwf2,iwf3,iwf4)
                !$$$         print *,"readscr    W:",scrw4(iwf1,iwf2,iwf3,iwf4)
                !$$$     &                          + scrv4(iwf1,iwf2,iwf3,iwf4)

                !         if (sw1) call rx("'checkorb should be done in advance")
                call checkorb2(iwf,jwf,kwf,lwf,ijklmag)
!!! ijklmag = F ---> W = 0 ; ijklmag = T ---> W = Wd
                !$$$         ijklmag=.true.
                if (ijklmag) then      !!! check iwf is derived from same atom
                   if (idorb(iwf)==2 .AND. idorb(jwf)==2  & !!only d-orbital
                      .and. idorb(kwf)==2 .and. idorb(lwf)==2) then
                      !$$$            if (.true.) then
                      !$$$            if (idorb(iwf)==idorb(jwf) .and. idorb(kwf)==idorb(lwf)) then !!!include sp orbitals
                      scrw_(ijwf,klwf)=scrw4(iwf1,iwf2,iwf3,iwf4) &
                           + scrv4(iwf1,iwf2,iwf3,iwf4)

                   else
                      scrw_(ijwf,klwf)=0d0
                   endif
                else
                   scrw_(ijwf,klwf)=0d0
                endif
                !         write(6,"('ijklmag:',4I3,2f9.4)") iwf,jwf,kwf,lwf,abs(scrw_(ijwf,klwf))


                !$$$         !!! write index i,j,k,l ---> ijwf,klwf
                !$$$         if (iwf==jwf .and. kwf==lwf) then
                !$$$            write(ifd,"('ijwf,klwf=',2i6)") ijwf,klwf
                !$$$         elseif (iwf==kwf .and. jwf==lwf) then
                !$$$            write(ife,"('ijwf,klwf=',2i6)") ijwf,klwf
                !$$$         endif
                !$$$         write(ifa,"(6i5)") iwf,jwf,kwf,lwf,ijwf,klwf

4004         enddo
4003      enddo
4002   enddo
4001 enddo
    !      close(ifd)
    !      close(ife)
    !      close(ifa)
    close(ifscrv)
    close(ifscrwv)
    call writescrw(scrw_) !! display matrix element of Wijkl
  end subroutine readscr

  !---------------------------------------
!!! identify if iwf is d-orbital or not
!!! checkorb iwf ---> lorb(1:s, 2:p, 3:d, 4:f, 5:g)
  subroutine checkorb(iwf_in,nwf_in,lorb_out)
    implicit none
    integer,intent(in)  ::iwf_in,nwf_in
    integer,intent(out) ::lorb_out   ! 2=d-orb
    integer::ifdorb,iiwf,ief,iwf

    if (sw1) then !!initialize

       print *,"checkorb nwf",nwf_in
       if ( .NOT. allocated(idorb)) allocate(idorb(nwf_in))
       open(newunit=ifdorb,file="Worb2lorb.d",form="unformatted")
       read(ifdorb) idorb(1:nwf_in)
       read(ifdorb) nclass_mlwf
       if ( .NOT. allocated(nbasclass_mlwf)) allocate(nbasclass_mlwf(nclass_mlwf+1))
       nbasclass_mlwf(0)=0
       read(ifdorb) nbasclass_mlwf(1:nclass_mlwf)
       close(ifdorb)

       !$$$         print *,"tttt",idorb
       !$$$         print *,"tttt",nclass_mlwf
       !$$$         print *,"tttt",nclass_mlwf,nbasclass_mlwf

       if (debug) then
          do iwf=1,nwf_in
             print *,"idorb check, iwf,idorb",sw1,idorb(iwf)
          enddo
       endif

       sw1=.false.
    endif
    lorb_out=idorb(iwf_in)

  end subroutine checkorb

  !---------------------------------------
!!!identify if iwf,jwf,kwf,lwf is derived from same MagAtom
!!! if not so, W=0 (or K=0) : return iwfmag=.False.
  subroutine checkorb2(iwf_in,jwf_in,kwf_in,lwf_in,iwfmag)
    intent(in) ::iwf_in,jwf_in,kwf_in,lwf_in
    intent(out)::                                    iwfmag
    integer ::   iwf_in,jwf_in,kwf_in,lwf_in,iwf
    logical(8) :: iwfmag
    integer::iclass,natom,r_mlwfs,r_mlwff
    logical::skip,sw_orb2
    if (sw1) call rx("checkorb2: call checkorb in advance")
    natom=0
    skip=.False.
    iwfmag=.False.
    sw_orb2=.False.
    do iclass=1,nclass_mlwf
!!! iwf
       if (sw_orb2) then
          continue
       elseif (nbasclass_mlwf(iclass-1) < iwf_in .AND. &
            iwf_in <= sum(nbasclass_mlwf(1:iclass)) ) then
          natom=iclass
          sw_orb2=.True.
       endif
    enddo
!!!   iwf_in comes from (r_mlwfs < iwf <= r_mlwff)
!!!   jwf (r_mlwfs < jwf_in <= r_mlwff)
    r_mlwfs=sum(nbasclass_mlwf(1:natom-1))
    r_mlwff=sum(nbasclass_mlwf(1:natom))
    if (r_mlwfs < jwf_in .AND. jwf_in <= r_mlwff ) then
       continue
    else
       skip=.True.
    endif
!!! kwf
    if ( .NOT. skip .AND. r_mlwfs < kwf_in .AND. kwf_in <= r_mlwff ) then
       continue
    else
       skip=.True.
    endif
!!! lwf
    if ( .NOT. skip .AND. r_mlwfs < lwf_in .AND. lwf_in <= r_mlwff ) then
       continue
    else
       skip=.True.
    endif
    if ( .NOT. skip) iwfmag = .TRUE. 
  end subroutine checkorb2
  !---------------------------------------
!!!identify iwf is from the intended atom(iddmat)
!!! if not so, return intended_iwf = .False.
!!! Now we can choose only 1 atom. It should be updated. Aug28,2019 Okumura
!!! it is useful to combine checkorb for identifying the intended orbital and atom
  subroutine checkorb3(iwf_in,iddmat,intended_iwf)
    implicit none
    integer,intent(in)  ::iwf_in
    integer(4),intent(in) ::iddmat
    logical(8),intent(out) ::intended_iwf
    integer::iclass,r_mlwfs,r_mlwff,iwf
    if (sw1) call rx("checkorb3: call checkorb in advance")

    natom=0
    intended_iwf=.False.

    ! c Example: Cu2MnAl
    ! c Cu1: 1 2 3 4 5 6 7 8 9; nbasclass_(1)=9
    ! c Cu2: 1 2 3 4 5 6 7 8 9; nbasclass_(2)=9
    ! c Mn : 1 2 3 4 5 6 7 8 9; nbasclass_(3)=9
    ! c Al : 1 2 3 4          ; nbasclass_(4)=4
    ! c if we want to extract Mn3d
    ! c r_mlwfs = (9+9), r_mlwff = (9+9+9) --> 18 < iwf <= 27 is a requirement
    ! c

!!!   iwf_in comes from (r_mlwfs < iwf <= r_mlwff)
    r_mlwfs=sum(nbasclass_mlwf(1:iddmat-1))
    r_mlwff=sum(nbasclass_mlwf(1:iddmat))
    if (r_mlwfs < iwf_in .AND. iwf_in<= r_mlwff ) then
       intended_iwf = .True.
    else
       intended_iwf = .False.
    endif
  end subroutine checkorb3


  subroutine writescrw(scrw_)
    implicit none
    integer::iwf,jwf,kwf,lwf,ijwf,klwf
    !      real(8),allocatable::rw_w(:,:,:,:), cw_w(:,:,:,:)
    integer:: iexc
    ! nteger::iwf,jwf,inwf,kwf,lwf,ijwf,klwf,nnwf,ijwf_j
    integer::ijwf_j
    complex(8),intent(in)::scrw_(:,:)

    do iexc=1,2
       do iwf=1,nwf
          ijwf=(1+nwf)*iwf-nwf

          if (iexc==1) then
             if (iwf==1) write (6,*) "--- W  (direct) -------------------"
             do jwf=1,nwf
                if (iwf==1) write (6,"(A6)",advance="no") "     "
                if (iwf==1) write (6,"(i3)",advance="no") jwf
             enddo
             if (iwf==1) write(6,*) ""

             write (6,"(i3)",advance="no") iwf
             do jwf=1,nwf
                ijwf_j=(1+nwf)*jwf-nwf
                write (6,"(f9.4)",advance="no") abs(scrw_(ijwf,ijwf_j))
             enddo
          else
             if (iwf==1) write (6,*) "--- W' (exchange) -----------------"
             do jwf=1,nwf
                if (iwf==1) write (6,"(A4)",advance="no") "    "
                if (iwf==1) write (6,"(i4)",advance="no") jwf
             enddo
             if (iwf==1) write(6,*) ""

             write (6,"(i3)",advance="no") iwf
             do jwf=1,nwf
                ijwf_j=(jwf-1)*nwf+iwf
                write (6,"(f9.4)",advance="no") abs(scrw_(ijwf_j,ijwf_j))
             enddo
          endif
          write(6,*) ""
       enddo
    enddo
  end subroutine writescrw
  !-------------------------------------------------

!!! create unit matrix
!!! imat(iwf,jwf)= 1 if iwf=jwf,kwf=lwf
!!!              = 0 otherwise
  subroutine wan_imat(nwf,imat_o)
    implicit none
    integer,intent(in)  :: nwf
    complex(8), intent(out) :: imat_o(nwf*nwf,nwf*nwf)
    integer:: iwf,kwf,ijwf,klwf

    imat_o=0d0
!!! diagonalize for iwf=jwf and kwf=lwf
    do iwf=1,nwf
       ijwf=(1+nwf)*iwf-nwf   !!! iwf=jwf
       do kwf=1,nwf
          klwf=(1+nwf)*kwf-nwf !!! kwf=lwf
          imat_o(ijwf,klwf)=(1d0,0d0)
       enddo
    enddo
    !! threshold
  end subroutine wan_imat

  !---------------------------------------
!!! extract zmat(ijwf,klwf) ---> eval_o(nnwf)
!!! sum of eval_o is trmat
  subroutine diagwan(zmat,eval_o)
    implicit none
    complex(8),intent(in) :: zmat(nwf*nwf,nwf*nwf)
    complex(8), intent(out) :: eval_o(nwf*nwf)
    integer:: iwf,kwf,ijwf,klwf

    eval_o=0d0
!!! diagonalize for iwf=jwf and kwf=lwf
    do iwf=1,nwf
       ijwf=(1+nwf)*iwf-nwf   !!! iwf=jwf
       do kwf=1,nwf
          klwf=(1+nwf)*kwf-nwf !!! kwf=lwf
          eval_o(nwf*(iwf-1)+kwf)=zmat(ijwf,klwf)
       enddo
    enddo
    !!
    call bubble_im(nwf*nwf,eval_o)
    !! threshold
  end subroutine diagwan
  !---------------------------------------
!!! diag for wannier matrix Im[K]
!!! sum of element in iwf=jwf and kwf=lwf
  subroutine diagwan_tr(zmat,trmat_o)
    implicit none
    complex(8),intent(in) :: zmat(nwf*nwf,nwf*nwf)
    complex(8), intent(out) :: trmat_o
    !      real(8):: pi,znorm
    integer:: iwf,kwf,ijwf,klwf

    !      pi = 4d0* atan(1d0)
    !      znorm=-1d0*pi
    trmat_o=0d0
!!! diagonalize for iwf=jwf and kwf=lwf
    do iwf=1,nwf
       ijwf=(1+nwf)*iwf-nwf   !!! iwf=jwf
       do kwf=1,nwf
          klwf=(1+nwf)*kwf-nwf !!! kwf=lwf
          trmat_o=trmat_o+zmat(ijwf,klwf)
       enddo
    enddo
    !! threshold
    !      if (abs(aimag(trmat_o)) < 1d-16) trmat_o=cmplx(dble(trmat_o),0d0,kind(0d0))
  end subroutine diagwan_tr
  !---------------------------------------
!!! write matrix element(i=j;k=l)
!!! For Hermite check
  subroutine writehmat(zmat,nwf,filename)
    implicit none
    complex(8),intent(in) :: zmat(nwf*nwf,nwf*nwf)
    integer, intent(in):: nwf
    character(*):: filename
    !      real(8):: pi,znorm
    integer:: iwf,jwf,ijwf,klwf,iffile
    open(newunit=iffile,file=filename(:len_trim(filename)))
    do iwf=1,nwf
       ijwf=(1+nwf)*iwf-nwf   !!! iwf=jwf
       do jwf=1,nwf
          klwf=(1+nwf)*jwf-nwf !!! kwf=lwf

          write(iffile,"(4i5,6E13.4)") iwf,jwf,ijwf,klwf, &
               zmat(ijwf,klwf),dconjg(zmat(klwf,ijwf)), &
               zmat(ijwf,klwf)-dconjg(zmat(klwf,ijwf))
       enddo
    enddo
    close(iffile)
  end subroutine writehmat

  !---------------------------------------
!!! write matrix element(i=j;k=l)
!!! For Hermite check
  subroutine writeddmat(zmat,nwf,filename,diag,zmat_o)
    intent(in)::          zmat,nwf,filename,diag
    intent(out)::                                zmat_o
    complex(8) ::  zmat(nwf*nwf,nwf*nwf)
    complex(8) ::  zmat_o(nwf*nwf,nwf*nwf)
    integer ::    nwf
    character(*):: filename
    logical:: diag
    logical(8)::intended_iwf
    real(8):: iddmat_in
    integer:: iwf,jwf,kwf,lwf, &
         ijwf,klwf,iffile, lorb
    integer:: iddmat
    open(newunit=iffile,file=filename(:len_trim(filename)))
!!! diagonal
    call getkeyvalue("GWinput","output_ddmat_atom",iddmat_in,default=1d0)
    iddmat=int(iddmat_in)
    write(6,"('output_ddmat_atom =',i5)") iddmat
    ijwf=0; klwf=0; zmat_o=(0d0,0d0)
    if (diag) then
       do 5001 iwf=1,nwf
          do 5002 jwf=1,nwf
             ijwf = ijwf + 1
             do 5003 kwf=1,nwf
                do 5004 lwf=1,nwf
                   klwf = klwf + 1

!!! checkorb3: identify whether iwf is from indicated atom (iddmat)
!!! checkorb:  identify whehter iwf is from 3d (lorb = .True.)
                   call checkorb3(iwf,iddmat,intended_iwf)
                   if ( .NOT. intended_iwf) then
                      zmat_o(ijwf,klwf) = (0d0,0d0)
                      cycle
                   endif
                   call checkorb3(jwf,iddmat,intended_iwf)
                   if ( .NOT. intended_iwf) then
                      zmat_o(ijwf,klwf) = (0d0,0d0)
                      cycle
                   endif
                   call checkorb3(kwf,iddmat,intended_iwf)
                   if ( .NOT. intended_iwf) then
                      zmat_o(ijwf,klwf) = (0d0,0d0)
                      cycle
                   endif
                   call checkorb3(lwf,iddmat,intended_iwf)
                   if ( .NOT. intended_iwf) then
                      zmat_o(ijwf,klwf) = (0d0,0d0)
                      cycle
                   endif

                   call checkorb(iwf,nwf,lorb)
                   if ( .NOT. lorb==2) then
                      zmat_o(ijwf,klwf) = (0d0,0d0)
                      cycle
                   endif
                   call checkorb(jwf,nwf,lorb)
                   if ( .NOT. lorb==2) then
                      zmat_o(ijwf,klwf) = (0d0,0d0)
                      cycle
                   endif
                   call checkorb(kwf,nwf,lorb)
                   if ( .NOT. lorb==2) then
                      zmat_o(ijwf,klwf) = (0d0,0d0)
                      cycle
                   endif
                   call checkorb(lwf,nwf,lorb)
                   if ( .NOT. lorb==2) then
                      zmat_o(ijwf,klwf) = (0d0,0d0)
                      cycle
                   endif
                   zmat_o(ijwf,klwf)=zmat(ijwf,klwf)

5004            enddo
5003         enddo
5002      enddo
5001   enddo

!!! off-diagonal
    else
       do 6001 iwf=1,nwf
          do 6002 jwf=1,nwf
             do 6003 kwf=1,nwf
                do 6004 lwf=1,nwf
!!! checkorb3: identify whether iwf is from indicated atom (iddmat)
!!! checkorb:  identify whehter iwf is from 3d (lorb = .True.)
                   call checkorb3(iwf,iddmat,intended_iwf)

                   call checkorb(iwf,nwf,lorb)
                   ijwf=(1+nwf)*iwf-nwf !!! iwf=jwf

6004            enddo
6003         enddo
6002      enddo
6001   enddo
    endif

    close(iffile)
  end subroutine writeddmat

  !---------------------------------------
  integer function ifile_handle() !! find open file handle
    implicit none
    integer:: i
    logical:: nexist
    do i=5001,9999
       inquire(unit=i,opened=nexist)
       if( .NOT. nexist) then
          ifile_handle=i
          return
       endif
    enddo
    stop 'ifile_handle: we did not find open file hundle'
  end function ifile_handle
  !---------------------------------------
  real(8) function det33(am)
    implicit none
    real(8),intent(in) :: am(3,3)
    det33= am(1,1)*am(2,2)*am(3,3) &
         -am(1,1)*am(3,2)*am(2,3) &
         -am(2,1)*am(1,2)*am(3,3) &
         +am(2,1)*am(3,2)*am(1,3) &
         +am(3,1)*am(1,2)*am(2,3) &
         -am(3,1)*am(2,2)*am(1,3)
  end function det33
  !---------------------------------------------------------------------
  subroutine bubble_im(n,array)
    ! ikinote, 2016/08/08
    implicit none
    integer,intent(in)::N
    complex(8),intent(inout)::array(1:N)
    integer::i,j
    complex(8)::t

    do i=1,N-1
       do j=i+1,N
          if(aimag(array(i)) > aimag(array(j)))then
             t=array(i)
             array(i)=array(j)
             array(j)=t
          end if
       end do
    end do

    return
  end subroutine bubble_im
end module m_readwan


