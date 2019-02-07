!! hx0fp0.m.F in order to develop by Okumura (2017/03/23)
!!  Calculate x0, \epsilon, spin susceptibility. 
!!
!! eps_lmf_cphipm mode is now commented out; you may need to recover this if necessary
!! (only epsPP_magnon_chipm mode works).
      program hmagnon
      use m_readwan,only: write_qdata, wan_readeigen, wan_readeval,
     &     readscr, checkorb, getqsym, checkorb2, diagwan, diagwan_tr, wan_imat,
     &     writehmat
      use m_ReadEfermi,only: readefermi,ef
      use m_readqg,only: readqg,readngmx
      use m_readeigen,only: readeval,init_readeigen,init_readeigen2
      use m_hamindex,only:qtt,nqtt
      use m_read_bzdata,only: read_bzdata,
     & ngrp2=>ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,qbas,ginv,
     & dq_,qbz,wbz,qibz,wibz,qbzw,
     & idtetf,ib1bz,idteti,
     & nstar,irk,nstbz
      use m_genallcf_v3,only: genallcf_v3,
     & nclass,natom,nspin,nl,nn, ngrp,
     & nlmto,nlnmx, nctot,niw, !nw_input=>nw,
     & alat, delta,deltaw,esmr,symgrp,clabl,iclass, !diw,dw,
     & invg, il,in,im,nlnm, 
     & plat, pos,ecore, symgg 
      use m_keyvalue,only: getkeyvalue
      use m_mpi,only: MPI__task,MPI__Initialize,MPI__Finalize,MPI__root,
     & MPI__hx0fp0_rankdivider2,
     & MPI__hmagnon_rankdivider, MPI__Initialize_magnon, MPI__consoleout_magnon,
     & MPI__Broadcast,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank=>mpi__rankMG,MPI__size=>mpi__sizeMG,
     & MPI__ranktab,MPI__consoleout,MPI__barrier,MPI__MEq,MPI__REAL8send,MPI__REAL8recv,MPI__AllreduceSum
!! frequency
      use m_freq,only: getfreq, !NOTE: call getfreq generate following data.
     & frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,wiw !, frhis0,nwhis0 !output of getfreq
!! tetwt
      use m_tetwt,only: tetdeallocate,gettetwt, !followings are output of 'L871:call gettetwt')
     &  whw,ihw,nhw,jhw,ibjb,nbnbx,nhwtot,n1b,n2b,nbnb 
!! w0 and w0i (head part at Gamma point)
      use m_w0w0i,only: w0w0i,
     & w0,w0i

      implicit none
!! ------------------------------------------------
!! We calculate chi0 by the follwoing three steps.
!!  gettetwt: tetrahedron weights
!!  x0kf_v4h: Accumlate Im part of the Lindhard function. Im(chi0) or Im(chi0^+-)
!!  dpsion5: calculate real part by the Hilbert transformation from the Im part
!!  eibz means extented irreducible brillowin zone scheme by C.Friedlich. (not so efficient in cases).
!!-------------------------------------------------

cccccc this may be wrong or correct cccccccccc
Cr Be careful for the indexing...
Cr      A routine idxlnmc(nindxv,nindxc,...  in index.f
Cr      specifies the order of the  (Core wave)+(Argumentation wave) in each MT.
Cr      The total number of the wave are mnl(ic)= mnlc(ic) + mnlv(ic).
Cr      The indexing starts with core first and then valence on top of core
Cr      So n-index in "in" for valence electron is different from "inv".
cccccccccccccccccccccccccccccccccccccccccccccccccc
      real(8):: q(3),  qgbin(3),qx(3)
      real(8):: ua=1d0 ! this is a dummy.
      integer:: ifrb(2),ifcb(2),ifrhb(2),ifchb(2) !,ifev(2)
      integer:: ndble=8
      integer:: nword
      real(8),allocatable:: vxcfp(:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:)
      complex(8),allocatable:: zxq(:,:,:),zxqi(:,:,:),zxq2(:,:,:)

!  tetrahedron method
      logical :: tetra=.true.
      ! so as to reduce the memory usage.
      complex(8) :: fff,img=(0d0,1d0)
      logical :: debug
      logical :: realomega=.true., imagomega=.true.

      real(8) :: omg2max, wemax

      integer::iopen,maxocc2,iclose,
     & ixc,iqxini,iqxend,iqxendx,
     &   nq0i,i,mxx,ini,ix,is
     &   ,iw,noccxv,noccx,iq,ngb
     &   ,nprecx,nblochpmx,mrecl,ifwd,nspinmx,ibas
     &   ,kx,isf,job
     &   ,ihis,ik,ibib,ib1,ib2,j

      integer:: incwfin,  verbose
      logical:: omitqbz=.false. !, noq0p

      real(8) :: erpaqw, trpvqw, trlogqw,rydberg,hartree
     &   ,pi,efz,qfermi,alpha,rs,voltot,ecelgas,efx,valn
      integer:: iqbz,iqindx,iflegas,nmx
     &   ,ifcor,nqitot,isx,ntot,ieclog,iww,iqq,ieceig,ecorr_on=-1 
      real(8) :: eclda_bh,eclda_pz,wk4ec,faca

      integer:: ifv,lxxilmx,ilm_r,nx_r,lb,nb,mb
      integer,allocatable:: nxx_r(:)
      real(8),allocatable:: svec(:,:),spinvec(:,:),consvec(:,:),cvec(:,:)
      character*3:: charnum3
      character*4:: charnum4
      character*5:: charnum5

      real(8)::schi=1d0
      integer:: new,nmxx,ii,iy,ipl1,ixx

      real(8),allocatable::SS(:),rwork(:),ss0(:)
      complex(8),allocatable:: UU(:,:),VT(:,:),work(:),ddd(:,:)
     & ,vtt(:,:),zzz(:,:),sqsvec(:),ooo(:,:),ppo(:,:) !,sqovlp(:,:),sqovlpi(:,:)
      integer::lwork,info,imin,ifzxq
      complex(8),allocatable:: UU0(:,:),VT0(:,:)

      logical ::  chipm=.false.,nolfco=.false. !sergeyv only ngczero=.false.,
     &  ,epsmode=.false.,normalm=.false., crpa=.false. 
      real(8):: qs,qt,ww,muu, ddq(3)

! Feb2006 time-reversal=off case
      logical :: timereversal, testtimer,onceww
      integer:: jpm,ncc

      integer:: istat,imb,imb1,imb2,nmbas_in
      integer,allocatable:: imbas(:), imbas_s(:),iibas(:)
      real(8):: fourpi,sqfourpi,tpioa,absq,vcou1,vcou1sq

      integer,allocatable:: nwgt(:,:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)
      integer:: ificlass,ifile_handle,k
      real(8)::ebmx
      integer:: nbmx,mtet(3)!,ifq0p

!     okumura
      integer::nwf,nsp_w,nqtt_w,iwf,jwf,inwf,iexc,kwf,lwf,ijwf,klwf,nnwf,ijwf_j
      integer::ifhamw,ifgas,ifchipmz_wan,ifchipmr_wan,ifchipmrk_wan
      real(8),allocatable::ev_w1(:,:),ev_w2(:,:)
      complex(8),allocatable::evc_w1(:,:,:),evc_w2(:,:,:)
      logical::wan=.true.,lhm,lsvd,nms,eibzmode
      logical(8):: ijklmag !for checkorb2
      complex(8),allocatable::kmat(:,:,:,:),wanmat(:,:)
     &     ,wkmat(:,:),wkmat2(:,:),rmat(:,:,:),rmat2(:,:)
!     &     ,swkwmat(:,:),swkwmat2(:,:)
!     &     ,sqemat(:,:),plmat(:,:),plpmat(:,:),plpmat_inv(:,:)
     &     ,wmat_check(:,:)
      complex(8)::trmat,trmatt,trmat1,trmat2
      complex(8),allocatable::imat(:,:) !unit matrix for 1-WK
      real(8):: wan_ecore(1), imagweight, www
      logical::npmtwo,diag=.false.,t2g
      complex(8)::wan_i,wan_j,wan_k,wan_l,wanijkl
      integer::igv
      real(8)::qlat(3,3),qsh1(3),qsh2(3),znorm,rnqbz, eta
      !For screening W
      complex(8),allocatable:: scrw(:,:) ,cmat2(:,:)
      complex(8),allocatable:: eval_wk(:),eval_k(:)
      integer:: it,itp,isdummy,lorb, size_lim
!!! q on symline
      real(8)::qrot(3),cr=2d-5
      logical(8)::init2=.true.,llsym,gskip
      logical(8)::d100,d110,d111,d1xx,dhpb,dxwf,dhnb,dnpb
      real(8)::rlatp(3,3),xmx2(3),qqin(3),qshort(3),qshort2
      integer:: nlatout(3,48),nout,iout,nqsym
      integer,allocatable:: iqlist(:)
!!! WK eigenvalue check
      integer::ifwkeigen,ifwkeigen2,ifwkeigen3
!!! temp variable
      logical(8)::lsmo,threshold=.true.
      integer::ifwanmat
!!! for MAX(Im[R])
      integer::imaximr=0
      real(8)::maximr,w_maximr
      complex(8)::sumrpa_maximr(1),summf_maximr(1)
      real(8),allocatable:: rpa_maximr(:),mf_maximr(:)

!! -------------------------------------------------------------------
      call getkeyvalue("GWinput","mpi_size_lim",size_lim,default=999)
      call MPI__Initialize_magnon()
      call MPI__consoleout_magnon('hmagnon',size_lim) ! size_lim for saving memory (avoid swapping)

      call cputid (0)
      hartree  = 2d0*rydberg()
      pi       = 4d0*datan(1d0)
      fourpi   = 4d0*pi
      sqfourpi = sqrt(fourpi)

      call cputid(0)

      write(6,*) " OK ixc=223    chipm_wannier sergey's "
      imagomega =.false.
      omitqbz=.false.
      epsmode = .true.
      chipm=.true.
      nolfco=.true.
      wan=.true.

!! We fix newaniso2=T now.

!! Prof.Naraga said " write(6,*)'Timereversal=',Timereversal()" 
!! here caused a stop in ifort ver.1x.x. Why? May be a compilar bug, and fixed now.

!! Readin BZDATA. See m_read_bzdata in gwsrc/rwbzdata.f
!! Read Bzdata; See use m_read_bzdata,only:... at the beginning of this routine.
      call read_BZDATA()
      
!! Read electron gas mode or not.
      call getkeyvalue("GWinput","lHermite",lhm,default=.false.)
      call getkeyvalue("GWinput","lsvd",lsvd,default=.false.)
      call getkeyvalue("GWinput","nms",nms,default=.false.)  !!! For NiMnSb

      write(6,*) "lsvd, lhm",lsvd,lhm
      write(6,*) "reduce mpi_size for saving memory (default =999) ",size_lim
      
c$$$!! check write
      if(MPI__root) then
         do i=1,nqbz
            if(i<10.or.i>nqbz-10) write(6,"('i qbz=',i8,3f8.4)") i,qbz(:,i)
            if(i==10.and.nqbz>18) write(6,"('... ')")
         enddo
         write(6,*)' !!nqbz nqibz =',nqbz,nqibz
      endif

!! EFERMI
cc need for gettetwt; automatically read file(okumura) OK
      call readefermi()

!! === Readin by genallcf. Set basic data for crystal
!! See "use m_genallcf_v3" at the begining of this routine
!!
      incwfin=0  !use ForX0 for core in GWIN
      call genallcf_v3(incwfin) !in module m_genallcf_v3
!     ! incwfin is for core handling
c     write(6,*)' ncore nw_input delta=', ncore,nw_input,delta
      write(6,"(' nqbz nqibz ngrp=',3i5)") nqbz,nqibz,ngrp
      if(ngrp/=ngrp2       ) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
      if(chipm.and.nspin==1) call rx( 'chipm mode is for nspin=2')

      tpioa=2d0*pi/alat

!     ! --- okumura Read dimensions of hamiltonian_wannier, spin, nqtt
c      if(ixc==223) then
      ifhamw=ifile_handle()
      open(ifhamw,file="wan4chi.d",form='unformatted')
      read(ifhamw) nwf,nsp_w,nqtt_w
      close(ifhamw)
      write(6,*) "nwf,nsp_w,nqtt",nwf,nsp_w,nqtt_w
c      endif

c     Weight for irreducible q-point (qibz)
c      do iq=1,nqibz
c         write(6,"('wibz',4f9.4)") wibz(iq),qibz(:,iq)
c      enddo
c
cccccccccccccccccc Screening W for magnon 
      nnwf=nwf*nwf
      allocate(scrw(nnwf,nnwf))
      call readscr(nwf,scrw)
cccc display W
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
                  write (6,"(f9.4)",advance="no") abs(scrw(ijwf,ijwf_j))
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
                  write (6,"(f9.4)",advance="no") abs(scrw(ijwf_j,ijwf_j))
               enddo
            endif
            write(6,*) ""
         enddo
      enddo
      scrw(:,:)=scrw(:,:)/hartree
      iqxend = nqibz !+ nq0i

!! Initialization of readEigen !readin m_hamindex
      ginv = transpose(plat)
      call minv33tp(plat,qlat)
      if(verbose()>50) print *,'eeee exit of init_readeigen2'

      do iq=1,nqibz
         iqbz = iqindx(qibz(:,iq),ginv,qbz,nqbz)
c     nstibz(iq) = nstbz(iqbz)
         write(6,"(' iq qibz nstibz=',2i5,3f9.4,i5)")iq,iqbz,qibz(:,iq) !,nstibz(iq)
      enddo


ccccccccccccccccccccccccccccccccccccc
!! We get frhis,freq_r,freq_i, nwhis,nw,npm,wiw  by getfreq
      wemax=    5d0 !max value for plot
      omg2max = wemax*.5d0+.2d0 ! (in Hartree) covers all relevant omega, +.2 for margin
c      if(MPI__root) write(6,"(' emin emax omega2max=',3f13.5)") emin, emax, omg2max
c
!! getfreq returun date given at " use m_freq,only:". 
c      if(.not.epsmode) call getwemax(lqall,wemax) !wemax is to determine nw !real axis divisions

!! NOTE: npmtwo=T sets npm=2
!! optional npmtwo is added aug2017
      call getfreq(epsmode,realomega,imagomega,tetra,omg2max,wemax,niw,ua,MPI__root, npmtwo=.true.)
      if(MPI__root) write(6,"(' nw npm=',2i5)") nw,npm
      if(.not.imagomega) niw=1

!! Tetra hedron initialization
      !if (wan) then
      noccxv = nwf
      print *,"noccxv ",noccxv
ccccccccccc
      noccx  = noccxv + nctot
cccccccccc
      nprecx = ndble  !We use double precision arrays only.
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()
      nspinmx = nspin
      iqxini=1
      mtet=(/1,1,1/) !dummy
c$$$      call getkeyvalue("GWinput","multitet",mtet,3,default=(/1,1,1/)) 
!okumura
      
!! === Use of symmetry. EIBZ procedure PRB81,125102 ===
 !! EIBZ mode memo for nolfco (right?) 
 !! If eibzmode=T, it is efficient but can slightly break crystal symmetry.(how much?)
 !! This is because band connectivity is judged by just from band ordering in tetrahedron weitht tetwt5.
 !!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
c     okumura
      eibzmode=.false.          !simple symmetry
c$$$c     end okumura
      iqxendx=iqxend
      allocate( nwgt(1,iqxini:iqxendx),igx(1,1,iqxini:iqxendx)
     &     ,igxt(1,1,iqxini:iqxendx), eibzsym(1,1,iqxini:iqxendx)) !dummy
      nwgt=1
            
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do is=1,nspinmx
         call wan_readeigen(qbz(:,:),nqbz,is,nsp_w)
         call write_qdata(ginv,nqbz,qbz(:,:))
      enddo      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      !!! For MPI test (intended symmetric iq)
      nqsym=0
      do iq=1,nqibz
         q = qibz(:,iq)         !qibze ! you can spefify any q, which you like.
         call getqsym(qlat,ngrp,iq,symgg,q,qshort,llsym)
         if (llsym) nqsym=nqsym+1
      enddo
      write(6,*) "nqsym:",nqsym
      call MPI__hmagnon_rankdivider(nqsym)

      write(6,*) "nqsym mpi_MEq(:)",nqsym,mpi__MEq
      allocate(rpa_maximr(mpi__MEq),mf_maximr(mpi__MEq)) !list of w(MAX(Im[R])): magnon peak
      rpa_maximr=0d0; mf_maximr=0d0

      !iqlist for syml
      !qibz(:,iqlist(iq)): intended q (iq=1,nqsym)
      allocate(iqlist(1:nqsym))
      nqsym=0
      do iq=1,nqibz
         q = qibz(:,iq)         !qibze ! you can spefify any q, which you like.
         call getqsym(qlat,ngrp,iq,symgg,q,qshort,llsym)
         if (llsym) then
            nqsym=nqsym+1
            iqlist(nqsym)=iq
            write(6,"('qshort (on syml) = ',I3,6f9.5)") iq,qshort,q
         endif
      enddo
c$$$      do iq=1,nqsym
c$$$         write(6,*) "iqlist",iqlist(iq),qibz(:,iqlist(iq))
c$$$      enddo

      !!! for MPI (finally nqsym ---> nqibz)
      do iq=1,nqsym
         if (MPI__task(iq)) write(6,'("iq,MPI_rank, mpi__ranktab",3I8)') iq,MPI__rank,mpi__ranktab(iq)
!         if (MPI__task(iq)) write(6,'("iq,qshort,I8,3f13.5")') iq,qibz(:,iqlist(iq))
      enddo
      call MPI__barrier()

!      call rx("aa")
c
!! ======== Loop over iq ================================
      do 1001 iq = 1,nqsym      ! NOTE: q=(0,0,0) is omitted when iqxini=2
         if(MPI__rank > size_lim) cycle !reduce mpi-size for test (skip 21-32)
         if(.not. MPI__task(iq) .and. iq .ne. 1) cycle
         q=qibz(:,iqlist(iq))
         if (iq == 1) then
            qshort=0d0; llsym=.true.
         else
            call getqsym(qlat,ngrp,iq,symgg,q,qshort,llsym) !!! exclude Gamma
            imaximr=imaximr+1 !!! imaximr=1,2,...,mpi__MEq
         endif
         if (llsym) write(6,"('qshort (on syml) = ',I3,3f9.5)") iq,qshort

!! Caution : confusing point
!!  ngc by QGcou is shown at the bottom of lqg4gw.
!!  ngc read from PPOVL are given by rdata4gw---> ngc(iq>nqibz )=ngc for q=0
         if( iq==1 ) then       ! *sanity check
            if(sum(q**2)>1d-10) then
               call rx( ' hx0fp0: sanity check. |q(iqx)| /= 0')
            endif
         endif

!!! write qshort instead of q
         write(6,"('===== do 1001: iq q(qshort)=',i7,3f9.4,' ========')")iqlist(iq),qshort !qq
         
         if (lhm) cycle

!     ! zxq and zxqi are the main output after Hilbert transformation
         allocate( zxq (nnwf,nnwf,nw_i:nw), zxqi (nnwf,nnwf,niw))
         zxq=0d0;  zxqi=0d0
         if(debug) write(6,*)' niw nw=',niw,nw
         
!     ! ==== spin chi_charge or chi_+- ====
         is=1
         isf=2
!! Tetrahedron weight.
!! output
!!     nbnbx
!!     ihw(ibjb,kx): omega index, to specify the section of the histogram.
!!     nhw(ibjb,kx): the number of histogram sections
!!     jhw(ibjb,kx): pointer to whw
!!     whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
!!     : histogram weights for given ib,jb,kx for histogram sections
!!     from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
c            write(6,*) ' --- goto x0kf_v4hz ---- newaniso= ',newaniso2

ccccccccccccccc simple band for check
         if(.not. allocated(ev_w1)) allocate(ev_w1(nwf,nqbz))
         if(.not. allocated(ev_w2)) allocate(ev_w2(nwf,nqbz))
         if(.not. allocated(evc_w1)) allocate(evc_w1(nwf,nwf,nqbz))
         if(.not. allocated(evc_w2)) allocate(evc_w2(nwf,nwf,nqbz))
          
         print *,"npm::",npm
         write(6,*) "nctot,ecore",nctot,ecore
         
         iwf=ifile_handle()
         print *,"wan_readeval start"
         
         do kx=1,nqbz           !!! ev_w1, ev_w2 unit: [Ry]
            call wan_readeval(  qbz(:,kx), is,
     o           ev_w1(1:nwf,kx), evc_w1(1:nwf,1:nwf,kx))
            call wan_readeval(q+qbz(:,kx), isf,
     o           ev_w2(1:nwf,kx), evc_w2(1:nwf,1:nwf,kx))            
         enddo
         
         write(6,"('= start wan_gettetwt =',2i6,3f9.4)") nwf,iq,q

ccccccccccctetrahedron part             
         call gettetwt(q,iq,isdummy,isdummy,nwgt(:,iq),frhis,nwhis,npm,
     i        qbas,ginv, ef, nqibz, nwf,ev_w1,ev_w2, nctot,wan_ecore,
     i        nqbz,qbz,nqbzw,qbzw,  ntetf,idtetf,ib1bz,
     i        nwf,ebmx,mtet,eibzmode,wan) !nov2016
         write(6,*) "=== end gettetwt"
         deallocate(ev_w1,ev_w2)
cccccccccnormalization of Im[K]:
         rnqbz=1/dble(nqbz)
         znorm=-1d0*pi

ccccccccc generate wanmat         
           if (.not. allocated(kmat)) allocate(kmat(1:nnwf,1:nnwf,1:nwhis,1:npm),
     &          wanmat(1:nnwf,1:nnwf)) !! (nwf^2*nwf^2)
!           kmat=(0d0,1d-15)
           kmat=(0d0,0d0)       !orbtest
           do 2011 kx=1,nqbz    !! discrete k-point loop
           do 2012 jpm=1,npm !! jpm=2: negative frequency
           do 2013 ibib=1,nbnb(kx,jpm) !! n,n' band loop
!!!     n1b(ibib,k,jpm) = n :band index for k (occupied)
!!!     n2b(ibib,k,jpm) = n':band index for q+k (unoccupied)
              if (debug) then
                 write(6,*) 'jpm,ibib',jpm,ibib
                 write(6,*) 'kx,nbnb',kx,nbnb(kx,jpm)
                 if (jpm==1) then
                    write(6,*) 'n1b(occ),n2b(unocc):',n1b(ibib,kx,jpm),n2b(ibib,kx,jpm)
                 else
                    write(6,*) 'n1b(unocc),n2b(occ):',n1b(ibib,kx,jpm),n2b(ibib,kx,jpm)
                 endif
              endif
              
c==========================caliculate M(:,:) : wanmat(nwf,nwf)
ccc   matrix elements from wannier eigenvector
c$$$  !!! evc_w1: upspin ; evc_w2: downspin
              it=n1b(ibib,kx,jpm) !index for n  for q
              itp=n2b(ibib,kx,jpm) !index for n' for q+k
              wanmat=(0d0,0d0)
              ijwf=0
c$$$ O(nwf^4)
              do 3002 iwf=1,nwf
              do 3003 jwf=1,nwf
                 ijwf=ijwf+1
                 klwf=0
              do 3004 kwf=1,nwf
              do 3005 lwf=1,nwf
                 klwf=klwf+1               
                 if ( ijwf <= klwf ) then !!! Herimite
c                   call checkorb2(iwf,jwf,kwf,lwf,ijklmag)
c                   ijklmag=.true. !test for FeFe
c                   print *,"checkorb2:",iwf,jwf,kwf,lwf,ijklmag !! ijklmag=T: {ijkl} in same atom
!                   if (ijklmag) then !!! same atom
!                 elseif (.true.) then !!! same atom
                   !!! calculate numerator of Kmatrix

                    wan_j=dconjg(evc_w1(jwf,it,kx)) !a_{Rj   beta}^{kn}*
                    wan_i=evc_w2(iwf,itp,kx) !a_{Ri  alpha}^{(k+q)n'}
                    wan_l=evc_w1(lwf,it,kx) !a_{R'l  beta}^{kn}
                    wan_k=dconjg(evc_w2(kwf,itp,kx)) !a_{R'k alpha}^{(k+q)n'}*
                    wanijkl=wan_j*wan_i*wan_k*wan_l

                    wanmat(ijwf,klwf)=wanijkl
                    wanmat(klwf,ijwf)=dconjg(wanijkl) !!! Suppose Hermite Kmat

c$$$                      write(ifwanmat,"('i j k l ijkltag abs(wan) Re(wan) Im(wan)',4I3,L4,5E13.5)") iwf,jwf,
c$$$     &                     kwf,lwf,ijklmag,abs(wanijkl),wanijkl
                   !!! wanmat (dimension:nnwf)
                 else
                    wanijkl=(0d0,0d0)
                 endif
c$$$                 if (abs(aimag(wanijkl))<1d-16) then
c$$$                    wanijkl=cmplx(dble(wanijkl),0d0,kind(0d0))
c$$$                 endif
!                 wanmat(ijwf,klwf)=wanijkl+(0d0,1d-9) ! +delta (~ 1d-9?)
                 ! delta makes a width of Im[R] peaks
                   
 3005  continue     !lwf
 3004  continue     !kwf
 3003  continue     !jwf
 3002  continue     !iwf
       debug=.False.
c========================== end caliculate M(:,:) : wanmat(nwf,nwf)
ccc   print *,"ihw,nhw",ihw(ibib,kx,jpm),nhw(ibib,kx,jpm)
       do iw=ihw(ibib,kx,jpm),ihw(ibib,kx,jpm)+nhw(ibib,kx,jpm)-1
          imagweight=whw(jhw(ibib,kx,jpm)+iw-ihw(ibib,kx,jpm))
!!!   accumulate Im[K]
          kmat(:,:,iw,jpm)=kmat(:,:,iw,jpm)+imagweight*wanmat(:,:)
       enddo
       
 2013 continue                  !ibib-loop
 2012 continue                  !jpm-loop
 2011 continue                  !k-loop

      deallocate(evc_w1,evc_w2,wanmat)
      call tetdeallocate()      !--> deallocate(ihw,nhw,jhw, whw,ibjb,n1b,n2b)

!! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
!     ! npm is set by npmtwo=T
      schi=1                    ! flip over maj/min spins.
      call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i     kmat, npm,nw_i, nnwf,nnwf,
     o     zxq, zxqi,
     i     .false., schi,1,1d99,1d99) !ecut(iecut),ecuts(iecut))
      deallocate(kmat)

      deallocate(zxqi) !!! not used in hmagnon (2018/06/20)

!     ! write matrix K(:,:,iq=200)
c      call writehmat(zxq(:,:,200) ,nwf, "zxqmat_check.dat")

      !!!Hermite for Kmat (zxq)
      if (lhm) then
         allocate(zxq2(nnwf,nnwf,nw_i:nw))
         zxq2=zxq
         zxq=0d0
         ijwf=0
         do 3006 iwf=1,nwf
         do 3007 jwf=1,nwf
            ijwf=ijwf+1
            klwf=0
            do 3008 kwf=1,nwf
            do 3009 lwf=1,nwf
               klwf=klwf+1               
               if ( ijwf < klwf ) then !!! Herimite
                  continue
                  
               elseif (ijwf == klwf) then
                  zxq(ijwf,ijwf,:)= zxq2(ijwf,ijwf,:)
               else             !!! ijwf > klwf 
                  zxq(ijwf,klwf,:)=( zxq2(ijwf,klwf,:) + dconjg(zxq2(klwf,ijwf,:)) )/2.0
                  zxq(klwf,ijwf,:)=dconjg( zxq(ijwf,klwf,:) )
               endif
 3009       enddo
 3008       enddo
 3007    enddo
 3006    enddo
         ijwf=0; klwf=0
         deallocate(zxq2)

      endif

      !!! threshold for Im[K] (zxq)
      threshold=.True.
      if (threshold) then
         do iw=nw_i,nw
            do iwf=1,nnwf
               do jwf=1,nnwf
                  if (log10(abs(aimag(zxq(iwf,jwf,iw)))) < -15) then
                     zxq(iwf,jwf,iw)=cmplx(dble(zxq(iwf,jwf,iw)),0d0,kind(0d0))
c$$$  if ( aimag(zxq(iwf,jwf,iw)) > 0d0 ) then
c$$$                     zxq(iwf,jwf,iw)=cmplx(dble(zxq(iwf,jwf,iw)),-1d0*aimag(zxq(iwf,jwf,iw)),kind(0d0))
                  endif
               enddo
            enddo
         enddo
      endif
      threshold=.False.
      allocate(wkmat(1:nnwf,1:nnwf),wkmat2(1:nnwf,1:nnwf)) !WKmatrix, WKmatrix_inv
      allocate(rmat(1:nnwf,1:nnwf,nw_i:nw)) !Rmatrix
      wkmat=0d0;wkmat2=0d0;rmat=0d0

!!! get eta for (1-eta*WK)
      if (iq==1) then

         if (.not. allocated(eval_wk)) allocate(eval_wk(nnwf))
         if (.not. allocated(eval_k)) allocate(eval_k(nnwf))

         wkmat(1:nnwf,1:nnwf)
     &        =matmul(scrw(1:nnwf,1:nnwf),zxq(1:nnwf,1:nnwf,0)) !omega=0

         call writehmat(scrw,nwf,"wmat_check.dat")
         
!!!   eval_wk is complex array because of Non-Hermite WK
         call diagcvuh3(wkmat(:,:),nnwf,eval_wk)
         !call diagwan(wkmat(:,:),eval_wk)

         if (debug) then
            write(6,"('eigenvalue of WK real&imag =',2E15.5)"
     &           ,advance='NO') eval_wk
            write(6,*)
         endif

!         eta= 1d0/minval(dble(eval_wk)) !moderate peak
         eta=-1d0/maxval(abs(eval_wk))
!         eta=  -1d0/abs(eval_wk(1)) !sharp peak
         write(6,*) "now eigenvalue abs(WK)",abs(eval_wk(1)),"is inversed"
         write(6,*) "check eigenvalue Re(WK)",real(eval_wk(1))
         write(6,*) "check eigenvalue Im(WK)",aimag(eval_wk(1))
         write(6,*) "wkmat calculated eta:", eta !negative value
      endif                     !iq==1

cccc open file each iq
      if (MPI__task(iq)) then
         ifchipmz_wan=ifile_handle()
         open(ifchipmz_wan,file="wan_ChiPMz.mat"//charnum4(iq))
         ifchipmr_wan=ifile_handle()
         open(ifchipmr_wan,file="wan_ChiPMr.mat"//charnum4(iq))          
         print *,'ifchipm=',ifchipmz_wan,ifchipmr_wan
c$$$         ifchipmrk_wan=ifile_handle()
c$$$         open(ifchipmrk_wan,file="wan_ChiPMr-k.mat"//charnum4(iq))          
      endif

!!! unit matrix (dimension: nwf*nwf), need for 1-WK matrix
      if (.not. allocated(imat)) then
         allocate(imat(1:nnwf,1:nnwf))
         imat=0d0
         if (nms) then
            do iwf=1,nnwf
               imat(iwf,iwf)=(1d0,-1d-6)
            enddo
            write(6,*) "imat(iwf,iwf)=",imat(1,1)
         else
            do iwf=1,nnwf
               imat(iwf,iwf)=(1d0,0d0)
            enddo
         endif
c$$$         imat=0d0
c$$$         call wan_imat(nwf,imat) !

      endif
ccccccc
      if (MPI__task(iq)) then
c$$$         ifwkeigen2=ifile_handle()
c$$$         open(ifwkeigen2,file="wk_eval_list.dat"//charnum4(iq))
         ifwkeigen3=ifile_handle()
         open(ifwkeigen3,file="1wk_eval_list.dat"//charnum4(iq))
      endif
cccccccdiagonalization for K
      print *,"nw_i,nw",nw_i,nw
      maximr=0d0; w_maximr=0d0 !!! search for Im[R] peak
      do 2050 iw = nw_i,nw
         trmat=(0d0,0d0)

         !!! diag 
         call diagcvuh3(zxq(:,:,iw),nnwf,eval_wk)
         !!! Hermite matrix diagonization
         !call diagcvh2(zxq(:,:,iw),nnwf,eval_wk) 
         !!!

         do iwf=1,nnwf
            if (abs(aimag(eval_wk(iwf))) < 1d-16) eval_wk(iwf)=cmplx(dble(eval_wk(iwf)),0d0,kind(0d0))
         enddo
c         call diagwan_tr(zxq(:,:,iw),trmat)
         !!!
         www=freq_r(iw)
         if(iw<0) www=-freq_r(-iw)
         if(.not. gskip .and. MPI__task(iq)) write(ifchipmz_wan,"(3f9.4,i6,E13.4,2x,2E13.4)") qshort,
!     &        iw,www*2d0,trmat/hartree
     &        iw,www*2d0,sum(eval_wk)/(hartree*znorm)

!!! K/(1-WK) = K(1-WK)^(-1)
!!! WK matrix
         wkmat(1:nnwf,1:nnwf)
     &        =eta*matmul(scrw(1:nnwf,1:nnwf),zxq(1:nnwf,1:nnwf,iw))

c??     &        =eta*matmul(scrw(1:nnwf,1:nnwf),abs(zxq(1:nnwf,1:nnwf,iw)))
cccc Hermite check (WKmat)


cccc
         if (iq==4 .and. iw==205) then
            eval_wk=0d0
            write(6,*) "wkmat check"
            call writehmat(wkmat,nwf,"wkmat_check.dat")
            if (lsvd) then
               call zgesvdnn(nnwf,zxq(:,:,iw),SS,UU,VT)
               do iwf=1,nnwf
                  write(6,'("wkmat eigenvalue:",f9.4)') SS(iwf)
               enddo
            endif
         endif
cccc
         if (debug) then
            call diagcvuh3(wkmat(:,:),nnwf,eval_wk)
            if (iw==45) then
               write(6,"('check sign eval WK:omega=45 case, iq=',i4)") iq
               write(6,"('check_sign eval WK:',2E13.5)"
     &              ,advance='NO') eval_wk(1:10)
            endif
            call diagwan(wkmat(:,:),eval_wk)
            if (iw==45) then
               write(6,"('check sign eval WK:omega=45 case, iq=',i4)") iq
               write(6,"('check_sign eval WK:',2E13.5)"
     &              ,advance='NO') eval_wk(1:10)
            endif
         endif

!!! 1-eWK ===> 1-WK
         wkmat2(1:nnwf,1:nnwf)=imat(1:nnwf,1:nnwf)-wkmat(1:nnwf,1:nnwf)
ccc Hermite check for 1-eWK
         if (iq==4 .and. iw==205) then
            call writehmat(wkmat2,nwf,"1wkmat_check.dat")
         endif

!!! eigenvalue of 1-WK
         if (debug) then
            if (iq==1 .and. iw==0) then
               call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
c               call diagwan(wkmat2(:,:),eval_wk)
               write(6,*) "q=0 and omega=0 case:"
               write(6,"('eigenvalue of (1-eta*WK) real&imag=',2E13.5)"
     &              ,advance='NO') eval_wk
               write(6,*)
            endif
         endif
cccc
         if (debug) then
            call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
            if (iw==45) then
               write(6,"('check sign eval 1-WK:omega=45 case, iq=',i4)") iq
               write(6,"('check_sign eval 1-WK:',2E13.5)"
     &              ,advance='NO') eval_wk(1:10)
            endif
            call diagwan(wkmat2(:,:),eval_wk)
            if (iw==45) then
               write(6,"('check sign eval 1-WK:omega=45 case, iq=',i4)") iq
               write(6,"('check_sign eval 1-WK:',2E13.5)"
     &              ,advance='NO') eval_wk(1:10)
            endif
         endif
cccc

ccc   WK eigenvalue writing (wk_eval_list.dat)
c$$$         call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)         
c$$$         if (MPI__task(iq)) write(ifwkeigen2,"(3f9.5,I6,10E13.5)") qshort,iw,www*2d0,
c$$$     &        eval_wk(1),abs(eval_wk(1))!,eval_wk(2),abs(eval_wk(2)),eval_wk(3),abs(eval_wk(3))
ccc   
         call matcinv(nnwf,wkmat2(1:nnwf,1:nnwf)) ! inv(1-WK)
         rmat(1:nnwf,1:nnwf,iw)=matmul(zxq(1:nnwf,1:nnwf,iw),
     &        wkmat2(1:nnwf,1:nnwf))

ccccc Hermite for Rmat
c$$$         allocate(rmat2(1:nnwf,1:nnwf))
c$$$         rmat2(1:nnwf,1:nnwf)=rmat(1:nnwf,1:nnwf,iw)
c$$$         ijwf=0
c$$$         do 4006 iwf=1,nwf
c$$$         do 4007 jwf=1,nwf
c$$$            ijwf=ijwf+1
c$$$            klwf=0
c$$$            do 4008 kwf=1,nwf
c$$$            do 4009 lwf=1,nwf
c$$$               klwf=klwf+1               
c$$$               if ( ijwf < klwf ) then !!! Herimite
c$$$                  continue
c$$$                  
c$$$               elseif (ijwf == klwf) then
c$$$                  rmat(ijwf,ijwf,:)= rmat2(ijwf,ijwf)
c$$$               else             !!! ijwf > klwf 
c$$$                  rmat(ijwf,klwf,:)=( rmat2(ijwf,klwf) + dconjg(rmat2(klwf,ijwf)) )/2.0
c$$$                  rmat(klwf,ijwf,:)=dconjg( rmat(ijwf,klwf,:) )
c$$$               endif
c$$$ 4009    enddo
c$$$ 4008    enddo
c$$$ 4007    enddo
c$$$ 4006    enddo
c$$$         deallocate(rmat2)
c$$$
ccccc


ccc   Taylor siries of (1-WK)inv (if WK << 1, (1-WK)inv = 1+WK )
c$$$         do ijwf=1,nnwf
c$$$            do klwf=1,nnwf
c$$$               if (abs(eta*aimag(wkmat(ijwf,klwf))) < 1d-16) then
c$$$                  wkmat2(ijwf,klwf)=1d0+eta*wkmat(ijwf,klwf)
c$$$               endif
c$$$            enddo
c$$$         enddo
c$$$         deallocate(wkmat)
cccq
!!! eigenvalue of (1-WK)inv
         if (debug) then
            if (iq==1 .and. iw==0) then
              call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
ccc
c               call diagwan(wkmat2(:,:),eval_wk)
ccc
               write(6,*) "q=0 and omega=0 case:"
               write(6,"('eigenvalue of (1-eta*WK)inv real&imag=',2E13.5)"
     &              ,advance='NO') eval_wk
               write(6,*)
            endif
ccc   find eigenvalue of (1-WK) at all omega (2018/04/24)
            if (iw==0) then
c               call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
               write(6,"('WKWK: eig 1-WK inv Re',3f9.4,f15.8)") qshort, dble(eval_wk(1))
               write(6,"('WKWK: eig 1-WK inv Im',3f9.4,f15.8)") qshort,aimag(eval_wk(1))
            endif
         endif

ccc 1-WK eigenvalue writing (1wk_eval_list.dat)
c          if (iq==1 .and. iw .ge. 0) then
c          if (iw .ge. 0) then ! omega >= 0

         call diagcvuh3(wkmat2(:,:),nnwf,eval_wk)
c         call diagwan(wkmat2(:,:),eval_wk)
         
!!! eigenvalue of rmat
         if (debug) then
            if (iq==1 .and. iw==0) then
c               call diagcvuh3(rmat(:,:,iw),nnwf,eval_wk)
               write(6,*) "q=0 and omega=0 case:"
               write(6,"('eigenvalue of rmat',i5,2f15.8)") iw,eval_wk(1)
               write(6,*)
            endif
         endif

cccc
         if (debug) then
            if (iw==45) then
               write(6,"('check sign eval Rmat: omega=45 case, iq=',i4)") iq
               write(6,"('check_sign eval Rmat:',2E13.5)"
     &              ,advance='NO') eval_wk(1:10)
               call diagwan(wkmat2(:,:),eval_wk)
               write(6,"('check sign eval Rmat: omega=45 case, iq=',i4)") iq
               write(6,"('check_sign eval Rmat:',2E13.5)"
     &              ,advance='NO') eval_wk(1:10)
            endif
         endif
cccc
         if (MPI__task(iq)) write(ifwkeigen3,"(3f9.5,I6,6E13.5)") qshort,iw,www*2d0,
     &        eval_wk(1) !,abs(eval_wk(1)),eval_wk(2),abs(eval_wk(2)),eval_wk(3),abs(eval_wk(3)) !iwf=1,2,3


!!! diagonalization for R
c         trmat=0d0
c         trmatt=0d0
c$$$         do iwf=1,nwf
c$$$            ijwf=(1+nwf)*iwf-nwf !!! iwf=jwf 
c$$$            do jwf=1,nwf
c$$$               klwf=(1+nwf)*jwf-nwf !!! kwf=lwf
c$$$c$$$               call checkorb2(iwf,iwf,jwf,jwf,ijklmag) ; ijklmag=.true.
c$$$c$$$               if(ijklmag) then
c$$$               trmat=trmat+rmat(ijwf,klwf,iw)
c$$$c               trmatt=trmatt+(rmat(ijwf,klwf,iw)-zxq(ijwf,klwf,iw))
c$$$               !endif
c$$$            enddo
c$$$         enddo

c         write(6,'("trmat_check3",81E12.4)') aimag(eval_wk)
         call diagwan(rmat(:,:,iw),eval_wk)
         trmat=sum(eval_wk(1:nnwf))
         call diagcvuh3(rmat(:,:,iw),nnwf,eval_wk)
c$$$         write(6,'("trmat_check diagwan and diagcvuh",3f9.4,i6,3E12.4)') qshort,iw,
c$$$     &        aimag(trmat),aimag(sum(eval_wk(1:nnwf))),aimag(sum(eval_wk(:)))
cccc
         if (MPI__task(iq)) then
            if(.not. gskip) write(ifchipmr_wan,"(3f9.4,i6,E12.4,x,12E12.4)")qshort,iw,
     &           www*2d0,trmat/hartree,
     &           www*2d0,sum(eval_wk(1:nnwf))/hartree
c            if(.not. gskip) write(ifchipmrk_wan,"(3f9.4,i6,E12.4,x,12E12.4)")qshort,iw,
c     &           www*2d0,trmatt/hartree

cccc search for MAX(Im[R]) 20180706 (0 - 1500 meV)
c            if (0d0 <www*hartree .and. www*hartree < 1.5) then
            if (0d0 < www) then
               if (maximr < -1d0*aimag(trmat)) then
                  maximr=-1d0*aimag(trmat)
                  w_maximr=www
               endif
            endif
cccc
         endif
ccc END K/(1-WK)

         gskip=.false.          !!! Skip Gamma point

cccc Hermite check (zxq and rmat)
cccc  
         if (iq==4 .and. iw==205) then
            call writehmat(zxq(:,:,iw),nwf,"kmat_check.dat")
            call writehmat(rmat(:,:,iw),nwf,"rmat_check.dat")
         endif
cccc  
     
 2050 continue
      deallocate(zxq, rmat, wkmat, wkmat2)

      if (.not. iq==1) then 
         !Tc(MF)
         ! BZweight*omega[eV] for sum(E(q))/N
         mf_maximr(imaximr)=wibz(iq)*w_maximr*hartree
         !Tc(RPA)
         ! BZweight*omega[eV] for sum(1/E(q))/N
         rpa_maximr(imaximr)=wibz(iq)/(w_maximr*hartree)
      endif

      if (MPI__task(iq)) then  !magnon peak
         write(6,"('AAAA',I4,3f9.4)") iq,qshort
         write(6,"('AAAA, w(MAX(im[R])), Im[R]',f13.5,E12.4)") w_maximr*2d0*13600,maximr/hartree
      endif

      if (MPI__task(iq)) then
c$$$         write(ifwkeigen2,*) 
         write(ifwkeigen3,*) 
         write(ifchipmz_wan,*)
         write(ifchipmr_wan,*)
c         write(ifchipmrk_wan,*)
      endif
      
      
c$$$  !! ImagOmega end =================
      
      debug=.False.
      
      if (MPI__task(iq)) then
c$$$         close(ifwkeigen2)
         close(ifwkeigen3)
         close(ifchipmz_wan)
         close(ifchipmr_wan)
c         close(ifchipmrk_wan)
      endif      
 1001 continue                  !q point loop
      write(6,*) "maximr RPA",rpa_maximr
      write(6,*) "maximr MFA",mf_maximr
      sumrpa_maximr=cmplx(sum(rpa_maximr(:)),0d0,kind(0d0))
      summf_maximr =cmplx(sum( mf_maximr(:)),0d0,kind(0d0))
c$$$!! =================== end of loop 1001 for q point ========================

      call MPI__barrier()

!! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===      
      write(6,*) "MPIcheck",MPI__size,MPI__rank,sumrpa_maximr
      if(MPI__size/=1) then
         call MPI__AllreduceSum(sumrpa_maximr(1),1)
         call MPI__AllreduceSum( summf_maximr(1),1)
      endif
      
      write(6,"('sum(E(q)) for MFA:',f9.4)") real(summf_maximr(1))
      write(6,"('[sum(1/E(q))]inv for RPA',f9.5)") 1d0/real(sumrpa_maximr(1))

      call cputid(0)
      call MPI__Finalize
      write(6,"('eta for 1-eta*WK:',f13.8)") eta
      call rx0( ' OK! hmagnon mode')

      end

c ----------------------------------------------------
      function matcinvf(a) result(b)
!!== Test routine for Inversion ==
      implicit none
      integer :: info,n,n2(2)
      integer,allocatable :: ipiv(:)
      complex(8):: a(:,:), b(1)
      complex(8),allocatable:: work(:)
      n2= SHAPE(a)
      n=n2(1)
      call zcopy(n,b,1,a,1)
      call zgetrf(n,n,a,n,ipiv,info)
      if(info/=0) then
        write(6,*)' matcinv: zegtrf info=',info
        call rx( ' matcinv: zegtrf ')
      endif
      allocate(work(n*n))
      call zgetri(n,a,n,ipiv,work,n*n,info)
      deallocate(work)
      if(info/=0) then
        write(6,*)'matcinv: zegtri info=',info
        call rx( 'matcinv: zegtri ')
      endif
      end
