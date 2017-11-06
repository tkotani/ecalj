      module m_readwan

      use m_keyvalue,only: getkeyvalue
      use m_iqindx_wan,only: iqindx2_wan
      use m_read_bzdata,only: read_bzdata, ngrp_w=>ngrp
      use m_genallcf_v3,only: genallcf_v3, symgg_=>symgg

!     use m_hamindex, only: iqimap !!tmp use; need to modify
      
      implicit none
      logical,private:: init=.true.
      real(8),private::alat,plat(3,3),ef
      integer,private::nwf,natom,nrws,n1,n2,n3
      
      integer,private,allocatable:: irws(:),ibaswf(:)
      real(8),private,allocatable:: pos(:,:),rws(:,:),drws(:)
      complex(8),private,allocatable:: hrotr(:,:,:), evecc(:,:)

      !!!!! wannier eigenvalue and eigenvector
c$$$      complex(8),allocatable:: evecc_w(:,:,:,:)
c$$$      real(8),allocatable:: eval_w(:,:,:)
      complex(8),private,allocatable:: evecc_w(:,:,:,:)
      real(8),private,allocatable:: eval_w(:,:,:)
      integer::iwf

      !!!d-orbital check
      integer,private,allocatable:: idorb(:), nbasclass_mlwf(:)
      integer,private:: nclass_mlwf
      logical::debug=.false.,sw1=.true.

!!! for getsyml
      logical:: init2=.true.
      contains

c---------------------------------------------
      subroutine readhrotr(is)
      implicit none
      integer(4):: is, ifh, n1, n2, n3, ifile_handle,nq
            
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
!! hrotr(nwf,nwf,nrws) = Hmn(R) = <0m|H |Rn> (m=1,nwf; n=1,nwf, R=1,nrws)
!! 
!! natom: number of atoms in the primitive cell.
!! ibaswf(nwf): atomic position for the Wannier orbital.
!! pos(3,natom): atomic poisition in the cell.
      write (6,"('--- Reading HrotRS  :: isp=',I6)") is
      ifh=ifile_handle()
      if(is==1) open(ifh,file='HrotRS.up',form='unformatted')
      if(is==2) open(ifh,file='HrotRS.dn',form='unformatted')
      read(ifh)alat,plat,natom
      if (is==1) allocate(pos(3,natom))
      read(ifh)pos
      read(ifh)ef 
      read(ifh)nwf,nrws,n1,n2,n3

      write (6,"('nwf,nrws,n1,n2,n3=',I3,I6,3I4)") nwf,nrws,n1,n2,n3
      if (is==1) then
         allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8),
     &        drws(n1*n2*n3*8),ibaswf(nwf), hrotr(nwf,nwf,nrws))
          !real space Hamiltonian in Wannier funciton basis
      endif
      read(ifh) irws,rws,drws,hrotr,ibaswf
      close(ifh)
      
      end subroutine readhrotr
c--------------------------------
      !!! generate eigenvalue list
      subroutine wan_readeigen(qbz,nqbz,is,nspin)
      implicit none

      integer,intent(in)::nqbz,nspin
      real(8),intent(in)::qbz(:,:)
      complex(8),allocatable:: hrotk(:,:,:),hrotkp(:,:),evec(:,:)
      real(8),allocatable:: eval(:)
      integer::iq,is

      !!! fat band plot test
      integer::iffb,ifile_handle,iband
      real(8):: rydberg
      logical::sw=.true.
      
      !! get HrotRS by file handling
      call readhrotr(is)

      if (init) then
         print *,"=== wan_readeigen::"
         allocate(hrotkp(nwf,nwf),evecc(nwf,nwf),eval(nwf))
         if (is==1) allocate(eval_w(nwf,nqbz,nspin),
     &        evecc_w(nwf,nwf,nqbz,nspin))
         do iq = 1,nqbz
!     write(6,*)' got get_hrotkp_ws iq =',iq
            call get_hrotkp_ws2(qbz(:,iq),hrotkp)
            call diag_hm2(hrotkp,nwf,eval,evecc)
            eval_w(1:nwf,iq,is)=eval
            evecc_w(1:nwf,1:nwf,iq,is)=evecc

ccccccccccc eigenvalue check
c$$$            if(iq==1) open(iwf,file="waneval_check.data")
c$$$!     if (kx==1 .and. iq > 100 ) open(iwf,file="wan_eval_check.data")
c$$$            if (qbz(2,iq)==1.0) then
c$$$               if (qbz(3,iq)==1.0) then
c$$$!     if (iq > 100) then
c$$$
c$$$                  write(iwf,"('q(1) ev_w1(:,kx)',6f9.4)") qbz(1,iq),eval_w(:,iq,1)
c$$$               endif
c$$$            endif
c$$$            if (iq==nqbz) close(iwf)
c$$$            if (iq==nqbz) call rx("check: end readwan...")

cccccccccccc   sumcheck for eigenvector  cccccccccccccccc
c$$$            if (iq==nqbz .and. is==2)  then
c$$$!     if (iq==1 .and. is==1) then
c$$$!     if (iq>50 .and. iq<70) then
c$$$               print *," iq:",iq
c$$$               do iwf=1,nwf 
c$$$                  write (6,"('sum evecc(up) evecc(dn) =',2f9.4)")
c$$$     &                 sum(abs(evecc_w(:,iwf,iq,1)**2)),
c$$$     &                 sum(abs(evecc_w(:,iwf,iq,2)**2))
c$$$               enddo
c$$$            endif
cccccccccccccccc  fat band plot (test: Gamma to H) ccccccccccccccccccc
c$$$            if (qbz(2,iq)/=0.0 .and. qbz(2,iq)/=1.0) cycle
c$$$            if (qbz(3,iq)/=0.0 .and. qbz(3,iq)/=1.0) cycle
c$$$            if (qbz(1,iq)>1.05) cycle
c$$$            write (6,"('SYML G-H  qbz',3f9.4)") qbz(:,iq)
c$$$
c$$$            !!! initial operation
c$$$            if (sw) then
c$$$               sw=.false.
c$$$               iffb=ifile_handle()
c$$$               if(is==1) open(iffb,file="fband.up.tmp")
c$$$               if(is==2) open(iffb,file="fband.dn.tmp")
c$$$               write(iffb,"('# ef nwf',f9.4,i4)") ef,nwf
c$$$               write(iffb,*)
c$$$            endif
c$$$            !!! write eigenvalue and eigenvector for wannier
c$$$            do iband=1,nwf
c$$$               write(iffb,"(i5,f13.5,' ',f13.6,i5,' ')",ADVANCE='NO')
c$$$     &          iq,qbz(1,iq),(eval_w(iband,iq,is)-ef)*rydberg(),iband
c$$$               do iwf=1,nwf
c$$$                  write(iffb,"(f13.6)",ADVANCE='NO')
c$$$     &                 (abs(evecc_w(iwf,iband,iq,is)))**2 
c$$$               enddo
c$$$               write(iffb,*)            
c$$$            enddo
c$$$
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
         enddo
         sw=.true.
         close(iffb)
      endif
      deallocate(hrotkp,evecc,eval)
      if (is==2) then
         init=.false.           !!!eigenvalue 
         print *,"=== wan_readeigen end"
         deallocate(pos,rws,drws,hrotr)
      endif
      if (allocated(hrotkp)) deallocate(hrotkp)
      end subroutine wan_readeigen
c--------------------------------
      subroutine wan_readeval(q,isp, ev_w, evc_w)
!      use m_hamindex, only: iqimap
      implicit none
      integer(4), intent(in) :: isp
      real(8), intent(in) :: q(3)
      real(8), intent(out) :: ev_w(nwf)
      complex(8), intent(out) :: evc_w(nwf,nwf)

      integer(4):: iq,iqindx,i
      real(8):: qu(3)
      
      if(init) call rx('wan_readeval: wan_readeigen should be called')
      call iqindx2_ wan(q, iq, qu) !qu is used q. q-qu is a G vector.
c$$$      write (6,"('wan_readeval: iq=',I7,'  ---> iqimap, qu',I7,3f9.4)")
c$$$     &     iq,iqimap(iq),qu

c$$$      write (6,"('eval_w=',5f9.4)") eval_w(1:nwf,1,1)
c$$$      write (6,"('eval_w=',5f9.4)") eval_w(1:nwf,1,2)

      !!! return wannier eigenvalue : ev_w
      ev_w(1:nwf) = eval_w(1:nwf,iq,isp)
      !ev_w(1:nwf) = eval_w(1:nwf,iqimap(iq),isp)
      !!! return wannier eigenvalue : evc_w(orbital,band)
      evc_w(1:nwf,1:nwf) = evecc_w(1:nwf,1:nwf,iq,isp)
      end subroutine wan_readeval

c--------------------------------almpost same get_hrotkp_ws in maxloc3
      !!!! hrotr (q) => hrotkp
      subroutine get_hrotkp_ws2(q,hrotkp)

c see Ref.[2] eq.26
      implicit real*8(a-h,o-z)
      complex(8) :: hrotkp(nwf,nwf),
     &     ci,cikr,ceikr,ctmp
      real(8) :: q(3)
      integer::im,in,ir

      pi = 4d0* atan(1d0)
      ci = (0d0,1d0)

      hrotkp = (0d0,0d0)

      do ir = 1,nrws
         rk = sum(rws(:,ir)*q(:))
         cikr = ci * 2d0 * pi * rk
         ceikr = exp(cikr) / dble(irws(ir))
         do im = 1,nwf
            do in = 1,nwf
               hrotkp(im,in) = hrotkp(im,in) + 
     &              ceikr * hrotr(im,in,ir)
            enddo
         enddo
      enddo
      end subroutine get_hrotkp_ws2
c--------------------------------
      subroutine diag_hm2(zmat,ndim,eval,evecc_)
      implicit real*8(a-h,o-z)
      integer:: ndim,i,nev,nmx
      complex(8),allocatable :: zmat2(:,:),ovlpc(:,:)
      complex(8):: zmat(ndim,ndim),evecc_(ndim,ndim)
      real(8):: eval(ndim),wk(ndim,11)
      integer iwk(ndim)

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

c---------------------------------------
      subroutine write_qdata(ginv,nqtt_in,qtt_in)
      implicit none
      integer::ifwqb,ifile_handle
      integer::nqtt_in
      real(8)::qtt_in(3,nqtt_in)
      real(8)::ginv(3,3)
      
      open(ifwqb,file="wanqbz",form='unformatted')
      write(ifwqb) ginv
      write(ifwqb) nqtt_in
      write(ifwqb) qtt_in
      close(ifwqb)
      
      end subroutine write_qdata      
c---------------------------------------
      subroutine readscr(nwf,scrw_)
      implicit none
      integer,intent(in)::nwf
      integer::ifscrwv,ifscrv,ifd,ife,ifa,ifile_handle
      integer::ir1,irws1
      character(len=9)::charadummy !dummy
      real(8)::rws1(3),freq,freq2 !dummy
      integer::is,iwf1,iwf2,iwf3,iwf4 !dummy
      integer::iwf,jwf,kwf,lwf,ijwf,klwf
!      real(8),allocatable::rw_w(:,:,:,:), cw_w(:,:,:,:)
      complex(8),allocatable::scrw4(:,:,:,:),scrv4(:,:,:,:)
      complex(8),allocatable,intent(out)::scrw_(:,:)
      integer::idummy
      logical(8)::ijklmag
      
      !!! set idorb: iwf ---> lorb
      call checkorb(1,nwf,idummy)
      
!     allocate(rw_w(nwf,nwf,nwf,nwf),cw_w(nwf,nwf,nwf,nwf))
      allocate(scrw4(nwf,nwf,nwf,nwf), scrv4(nwf,nwf,nwf,nwf))
      allocate(scrw_(nwf*nwf,nwf*nwf));scrw_=0d0
      ifscrwv=ifile_handle()
      open(ifscrwv,file="Screening_W-v.UP",form="formatted") !only up
      ifscrv=ifile_handle()
      open(ifscrv,file="Coulomb_v.UP",form="formatted") !only up

      !!! write direct index (ijkl)
      ifd=ifile_handle()
      open(ifd,file="ijkl_direct.d",form="formatted")
      !!! write exchange index (ijkl)
      ife=ifile_handle()
      open(ife,file="ijkl_exchange.d",form="formatted")
      !!! write all index (i,j,k,l --> ijwf,klwf)
      ifa=ifile_handle()
      open(ifa,file="ijkl_all.d",form="formatted")
      write(ifa,*) "# iwf jwf kwf lwf ijwf klwf"
      
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
         read(ifscrv,"(A,2i5, 3f12.6, 5i5,2f12.6)")
     &        charadummy,ir1, irws1, rws1 !dummy
     &        ,is,iwf1,iwf2,iwf3,iwf4 !dummy
     &        ,scrv4(iwf1,iwf2,iwf3,iwf4)

         !!! Screened Coulomb (W-v)
         read(ifscrwv,"(A,2i5, 3f12.6, 5i5,4f12.6)")
     &        charadummy,ir1, irws1, rws1 !dummy
     &        ,is,iwf1,iwf2,iwf3,iwf4,freq, freq2 !dummy
     &        ,scrw4(iwf1,iwf2,iwf3,iwf4)


c$$$         print *,"readscr     :",iwf1,iwf2,iwf3,iwf4
c$$$         print *,"readscr  W-v:",scrw4(iwf1,iwf2,iwf3,iwf4)
c$$$         print *,"readscr    v:",scrv4(iwf1,iwf2,iwf3,iwf4)
c$$$         print *,"readscr    W:",scrw4(iwf1,iwf2,iwf3,iwf4)
c$$$     &                          + scrv4(iwf1,iwf2,iwf3,iwf4)

!         if (sw1) call rx("'checkorb' should be done in advance")
         call checkorb2(iwf,jwf,kwf,lwf,ijklmag)
         !!! ijklmag = F ---> W = 0 ; ijklmag = T ---> W = Wd
         if (ijklmag) then      !!! check iwf is derived from same atom
            if (idorb(iwf)==2 .and. idorb(jwf)==2 !!!only d-orbital
     &           .and. idorb(kwf)==2 .and. idorb(lwf)==2) then
               scrw_(ijwf,klwf)=scrw4(iwf1,iwf2,iwf3,iwf4)
     &              + scrv4(iwf1,iwf2,iwf3,iwf4)
            else
               scrw_(ijwf,klwf)=0d0
            endif
         else
            scrw_(ijwf,klwf)=0d0
         endif
c         write(6,"('ijklmag:',4I3,2f9.4)") iwf,jwf,kwf,lwf,abs(scrw_(ijwf,klwf))
            

         !!! write index i,j,k,l ---> ijwf,klwf
         if (iwf==jwf .and. kwf==lwf) then
            write(ifd,"('ijwf,klwf=',2i6)") ijwf,klwf
         elseif (iwf==kwf .and. jwf==lwf) then
            write(ife,"('ijwf,klwf=',2i6)") ijwf,klwf
         endif
         write(ifa,"(6i5)") iwf,jwf,kwf,lwf,ijwf,klwf
         
 4004 continue
 4003 continue
 4002 continue
 4001 continue
      close(ifd)
      close(ife)
      close(ifscrv)
      close(ifscrwv)

      end subroutine readscr

c---------------------------------------
      !!! identify if iwf is d-orbital or not
      !!! checkorb iwf ---> lorb(1:s, 2:p, 3:d, 4:f, 5:g) 
      subroutine checkorb(iwf_in,nwf_in,lorb_out) 
      implicit none
      integer,intent(in)  ::iwf_in,nwf_in
      integer,intent(out) ::lorb_out   ! 2=d-orb
      integer::ifdorb,iiwf

      if (sw1) then !!initialize

         print *,"checkorb nwf",nwf_in
         if (.not. allocated(idorb)) allocate(idorb(nwf_in))
         ifdorb=ifile_handle()
         open(ifdorb,file="Worb2lorb.d",form="unformatted")
         read(ifdorb) idorb(1:nwf_in)
         read(ifdorb) nclass_mlwf
         if (.not. allocated(nbasclass_mlwf)) allocate(nbasclass_mlwf(nclass_mlwf+1))
         nbasclass_mlwf(0)=0
         read(ifdorb) nbasclass_mlwf(1:nclass_mlwf)
         close(ifdorb)

c$$$         print *,"tttt",idorb
c$$$         print *,"tttt",nclass_mlwf
c$$$         print *,"tttt",nclass_mlwf,nbasclass_mlwf
         
         if (debug) then
            do iwf=1,nwf_in
               print *,"idorb check, iwf,idorb",sw1,idorb(iwf)
            enddo
         endif

         sw1=.false.
      endif      
      lorb_out=idorb(iwf_in)   
      
      end subroutine checkorb

c---------------------------------------
      !!!identify if iwf,jwf,kwf,lwf is derived from same MagAtom
      !!! if not so, W=0 (or K=0) : return logical iwfmag=F
      subroutine checkorb2(iwf_in,jwf_in,kwf_in,lwf_in,iwfmag) 
      implicit none
      integer,intent(in)  ::iwf_in,jwf_in,kwf_in,lwf_in
      logical(8),intent(out) ::iwfmag
      integer::iclass,natom,r_mlwf
      logical::skip

      
      if (sw1) call rx("checkorb2: call checkorb in advance")
      
      natom=0
      skip=.False.
      iwfmag=.False.
      do iclass=1,nclass_mlwf
         !!! iwf
         if (nbasclass_mlwf(iclass-1) < iwf_in .and.
     &        iwf_in <= sum(nbasclass_mlwf(1:iclass)) ) then
            natom=iclass
         endif
      enddo
!!!   jwf
      r_mlwf=sum(nbasclass_mlwf(1:natom))
!!!      print *,r_mlwf
      if (nbasclass_mlwf(natom-1) < jwf_in .and.
     &     jwf_in <= r_mlwf ) then
      else
         skip=.True.
      endif
!!! kwf
      if (.not. skip .and. nbasclass_mlwf(natom-1) < kwf_in .and.
     &     kwf_in <= r_mlwf ) then
      else
         skip=.True.
      endif
!!! lwf
      if (.not. skip .and. nbasclass_mlwf(natom-1) < lwf_in .and.
     &     lwf_in <= r_mlwf) then
         iwfmag=.True.
      endif
      end subroutine checkorb2

c---------------------------------------
      subroutine getqsym(qlat_in, ngrp_w, iq, symgg_, qin, qshort, llsym_o)
      real(8),intent(in) :: qlat_in(3,3), qin(3),symgg_(3,3,48)
      integer, intent(in) :: ngrp_w, iq
      real(8),intent(out):: qshort(3)
      logical(8), intent(out) :: llsym_o
      
      real(8)::qrot(3),cr=2d-5
      logical(8)::gskip
      logical(8)::d100=.false.,d110=.false.,d111=.false.,d1xx=.false.
     &     ,dhpb=.false.,dwlf=.false.,dxwf=.false.,dhnb=.false.
     &     ,dnpb=.false.,dxms=.false.
      real(8)::rlatp_(3,3),xmx2_(3),qqin(3),qshort2
      integer:: nlatout(3,48),nout,iout,ig
      character*5:: syml

      !!! get syml line
      call getkeyvalue("GWinput","qsym",syml,default="Gt100")          
      llsym_o=.true.      
      if (syml=="100") then
         d100=.true.           !0 0 0 --> 1 0 0
      else if (syml=="110") then
         d110=.true.           !0 0 0 --> 1 1 0
      else if (syml=="111") then
         d111=.true.           !0 0 0 --> 1 1 1
      else if (syml=="fccS") then
         d1xx=.true.           !S (fcc)
      else if (syml=="fccWL") then
         dwlf=.true.            ! W to L (fcc) !.5 1 0 --> .5 .5 .5
      else if (syml=="fccXW") then
         dxwf=.true.            ! X to W (fcc) !1 0 0 --> 1 .5 0
      else if (syml=="bccHN") then
         dhnb=.true.           ! G (bcc): H to N !1 0 0 --> .5 .5 0
      else if (syml=="bccNP") then
         dnpb=.true.           ! D (bcc): N to P !.5 .5 0 --> .5 .5 .5
      else if (syml=="bccHP") then
         dhpb=.true.           !F (bcc): H to P !1 0 0 --> .5 .5 .5
      else if (syml=="scXM") then
         dxms=.true.           !F (sc): X to M !.5 0 0 --> .5 .5 0
      else
         write(6,*) "not supported symmetry line"
      endif
         
      if (iq==1) then
         if(d1xx .or. dxwf .or. dhnb .or. dnpb .or.dhpb) gskip=.true.
      else
         do ig=1,ngrp_w
            if (.not. llsym_o) cycle
            qrot(1:3)=matmul(symgg_(:,:,ig),qin(:))
            if(init2) then
               call shortn3_initialize(qlat_in,rlatp_,xmx2_)
               init2=.false.
            endif
            qqin = matmul(transpose(plat),qrot)
            call shortn3(rlatp_,xmx2_, qqin,ngrp_w, nout,nlatout)
            iout=1
            qshort(1:3)=  matmul(qlat_in(:,:), qqin+nlatout(:,iout))
            
c$$$            write(6,"('xxxx iq q ig qshort'i5,3f8.4,i3,3f8.4)")
c$$$     &           iq,qin(:),ig,qshort(:)

            qshort2=sqrt(dot_product(qshort,qshort))
            qshort=int(1d6*qshort)/1d6

cccccccccccccccccsymline check
            if (d100) then      !(100)
               if (abs((qshort2)-abs(qshort(1))) > cr
     &              .and. abs(abs(qshort2)-abs(qshort(2))) > cr
     &              .and. abs(abs(qshort2)-abs(qshort(3))) > cr) then
                  llsym_o=.false.
               endif

            elseif (d110) then  !(110)
               if (abs(qshort(1)) > cr
     &              .and. abs(qshort(2)) > cr
     &              .and. abs(qshort(3)) > cr) then
                  llsym_o=.false.
                  cycle
               endif

!     ! exclude (100) direction
               if (abs(abs(qshort(1))-abs(qshort(2))) < cr
     &              .and. abs(qshort(1)) < cr) llsym_o=.false.
               if (abs(abs(qshort(2))-abs(qshort(3))) < cr
     &              .and. abs(qshort(2)) < cr) llsym_o=.false.
               if (abs(abs(qshort(3))-abs(qshort(1))) < cr
     &              .and. abs(qshort(3)) < cr) llsym_o=.false.
               
               if (abs(abs(qshort(1))-abs(qshort(2))) > cr
     &              .and. abs(abs(qshort(1))-abs(qshort(3))) > cr
     &              .and. abs(abs(qshort(2))-abs(qshort(3))) > cr) then
                  llsym_o=.false.
               endif

            elseif (d111) then  !(111)
               if (abs(abs(qshort(1))-abs(qshort(2))) > cr
     &              .or. abs(abs(qshort(1))-abs(qshort(3))) > cr) then
                  llsym_o=.false.
               endif

            elseif (d1xx) then  !(1xx) on S (fcc)
               if (abs(abs(qshort(1))-1d0) > cr
     &              .and. abs(abs(qshort(2))-1d0) > cr
     &              .and. abs(abs(qshort(3))-1d0) > cr) then
                  llsym_o=.false.
                  cycle
               endif
               if (abs(abs(qshort(1))-abs(qshort(2))) > cr
     &              .and. abs(abs(qshort(1))-abs(qshort(3))) > cr
     &              .and. abs(abs(qshort(2))-abs(qshort(3))) > cr) then
                  llsym_o=.false.
               endif

            elseif (dwlf) then  !  W to L (fcc): 48 lines (pm1/2, pm(1-pmx), pmx)
c$$$               if     (abs(abs(qshort(1))-5d-1) < cr ) then
c$$$               elseif (abs(abs(qshort(2))-5d-1) < cr ) then
c$$$               elseif (abs(abs(qshort(3))-5d-1) < cr ) then

               if     (abs(abs(qshort(1))-0.5d0) < cr .and.
     &              abs(abs(qshort(2))+abs(qshort(3))-1d0) < cr ) then
               elseif (abs(abs(qshort(2))-0.5d0) < cr .and.
     &              abs(abs(qshort(3))+abs(qshort(1))-1d0) < cr ) then
               elseif (abs(abs(qshort(3))-0.5d0) < cr .and.
     &              abs(abs(qshort(1))+abs(qshort(2))-1d0) < cr ) then
               else 
                  llsym_o=.false.
               endif
            
            elseif (dxwf) then     !(1xx) on S (fcc)
               if (abs(abs(qshort(1))-1d0) > cr
     &              .and. abs(abs(qshort(2))-1d0) > cr
     &              .and. abs(abs(qshort(3))-1d0) > cr) then
                  llsym_o=.false.
                  cycle
               endif
               if (abs(qshort(1)) > cr
     &              .and. abs(qshort(2)) > cr
     &              .and. abs(qshort(3)) > cr) then
                  llsym_o=.false.
               endif

            elseif (dhnb) then  ! G: H to N (bcc)
               if (abs(qshort(1)) > cr
     &              .and. abs(qshort(2)) > cr
     &              .and. abs(qshort(3)) > cr) then
                  llsym_o=.false.
                  cycle
               endif
               if (abs(abs(qshort(1))+abs(qshort(2))-1d0) < cr) then
               elseif (abs(abs(qshort(2))+abs(qshort(3))-1d0) < cr) then
               elseif (abs(abs(qshort(3))+abs(qshort(1))-1d0) < cr) then
               else 
                  llsym_o=.false.
               endif

            elseif (dnpb) then  ! D: N to P (bcc)
               if (abs(abs(qshort(1))-5d-1) < cr
     &              .and. abs(abs(qshort(2))-5d-1) < cr) then
               elseif (abs(abs(qshort(2))-5d-1) < cr
     &                 .and. abs(abs(qshort(3))-5d-1) < cr) then
               elseif (abs(abs(qshort(3))-5d-1) < cr
     &                 .and. abs(abs(qshort(1))-5d-1) < cr) then
               else
                  llsym_o=.false.
               endif
               
            elseif (dhpb) then  ! F: H to P (bcc)
               if (abs(abs(qshort(1))-abs(qshort(2))) > cr
     &              .and. abs(abs(qshort(1))-abs(qshort(3))) > cr
     &              .and. abs(abs(qshort(2))-abs(qshort(3))) > cr) then
                  llsym_o=.false.
                  cycle
               endif
               if (int(1d5*abs(qshort(1)))/=int(1d5*(1d0-abs(qshort(2))))
     &              .and. int(1d5*abs(qshort(2)))/=int(1d5*(1d0-abs(qshort(3))))
     &              .and. int(1d5*abs(qshort(3)))/=int(1d5*(1d0-abs(qshort(1))))) then
                  llsym_o=.false.
               endif           

            elseif (dxms) then  ! F: H to P (bcc)
               if (abs(qshort(1)) > cr
     &              .and. abs(qshort(3)) > cr
     &              .and. abs(qshort(2)) > cr) then
                  llsym_o=.false.
                  cycle
               endif

               if (abs(abs(qshort(1))-5d-1) > cr
     &              .and. abs(abs(qshort(2))-5d-1) > cr
     &              .and. abs(abs(qshort(3))-5d-1) > cr ) then
                  llsym_o=.false.
               endif
            endif
         enddo
      endif
      end subroutine getqsym

c---------------------------------------
      integer function ifile_handle() !! find open file handle
      implicit none
      integer:: i
      logical:: nexist
      do i=5001,9999
         inquire(unit=i,opened=nexist)
         if(.not.nexist) then
            ifile_handle=i
            return
         endif
      enddo
      stop 'ifile_handle: we did not find open file hundle'
      end function
c---------------------------------------            
      real(8) function det33(am)
      implicit none
      real(8),intent(in) :: am(3,3)
      det33= am(1,1)*am(2,2)*am(3,3)
     &       -am(1,1)*am(3,2)*am(2,3)
     &       -am(2,1)*am(1,2)*am(3,3)
     &       +am(2,1)*am(3,2)*am(1,3)
     &       +am(3,1)*am(1,2)*am(2,3)
     &       -am(3,1)*am(2,2)*am(1,3)
      end function

c---------------------------------------      
      end module m_readwan


