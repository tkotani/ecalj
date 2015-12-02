      program hqpe_qsgw

c Jul,2000 takao mod from hqpe by ferdi. No routines required.
c calculates quasiparticle energies
c E(k,t) = e(k,t) + Z [SEx(k,t) + SEc(k,t) - xcLDA(k,t)]
c e(k,t) = LDA eigenvalue
c Z      = [1 - dSEc(e(k,t))/dw]^(-1)
c SEx(k,t)   = <psi(k,t)| SEx |psi(k,t)>
c SEc(k,t)   = <psi(k,t)| SEc |psi(k,t)>, SEc = GWc
c xcLDA(k,t) = <psi(k,t)| vxc |psi(k,t)>
c SEx and xcLDA are in file SEX
c SEc is in file SEC
      use keyvalue
      implicit real*8 (a-h,o-z)

c local data
      logical laf
      dimension ifsex(2),ifsexcore(2),ifxc(2),ifsec(2),ifqpe(2)
     & ,iftote(2),iftote2(2)
      integer(4),allocatable :: itxc(:),itc(:),itx(:)
      real(8),allocatable :: qxc(:,:,:),eldaxc(:,:),vxc(:,:),
     &    qc(:,:,:),eldac(:,:),sex(:,:),sexcore(:,:),
     &    qx(:,:,:),eldax(:,:),rsec(:,:,:),csec(:,:,:),zfac(:,:)

      integer:: ret,iix
c      logical:: readgwinput
      logical :: nozmode=.false.

! swj added
      integer(4):: n, iwcomg, nw, ninter, iinter, iw_int, iwini_int, iwend_int
      character(len=1024) :: cip, ci, sss
      character(len=1024) :: fname_sec
      dimension ifsecomg(2)

      integer(4),allocatable :: itcomg(:)
      real(8) :: rdenomi, cdenomi, xxxx!, polinta
      real(8),allocatable:: spect_ftn(:,:), spect_ftn_tot(:),
     &                   spect_ftn_int(:,:), spect_ftn_tot_int(:)
      real(8),allocatable :: qcomg(:,:),eldacomg(:,:), wibz(:),
     &        omega(:,:), rsecomg(:,:,:), csecomg(:,:,:), eh(:,:),
     &        omega_int(:,:), rsecomg_int(:,:,:), csecomg_int(:,:,:)
      complex(8),allocatable:: green_ftn(:,:,:),  sigma(:,:,:), constant_shift(:,:),
     &                    green_ftn_int(:,:,:),  sigma_int(:,:,:)
      real(8),allocatable:: rgreen_ftn2(:,:,:), cgreen_ftn2(:,:,:)
      integer :: ifoutgreen, ifoutsigma, ifoutspect1, ifoutspect2, ifoutsec_int
      integer :: ifoutgreen_int, ifoutsigma_int, ifoutspect1_int, ifoutspect2_int


! swj added



c shift quasiparticle energies (eV)
      write (*,*)' q+band index for zero?'
      read (*,*)jin
      if(jin>=1000) then
         jin=jin-1000
         nozmode=.true.
      endif
      call headver('hqpe',jin)

c open files
c> input files
c      ifqpnt     = iopenxx ('QPNT')
c      if(readgwinput()) then
      call getkeyvalue("GWinput","ninter_omega",ninter,default=1)
      call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
      !ninter = 1
c      else
c        ifqpnt     = iopenxx ('QPNT')
c      endif
      laf        = .false.

      call readx   (ifqpnt,10)
      read (ifqpnt,*) iqall,iaf
      if (iaf .eq. 1) laf = .true.


! swj added begin

      call readxxx(ifqpnt)     !skip to ***** for q point for spectrum function.
      omegamaxin = 1d70
      read (ifqpnt,*,err=2038,end=2038) dwplot, omegamaxin
      close(ifqpnt)
 2038 continue
      !print *, "dwplot", dwplot, "omegamaxin", omegamaxin

!        if(omegamaxin > omegamax) then
!          write(6,*)' --omegamaxin > omeagmax <- This case not implemented yet. -SWJ'
!          This part can cause the problem.
!          stop
!        endif
      write(6,*)' --Use readin dwplot and omegamaxin from <QPNT>'
      omegamax = omegamaxin
      if( omegamax <0) call rx( 'hsfp0 :strange omegamax <0 ')
      iwini =  -int( omegamax / dwplot )
      iwend =   int( omegamax/  dwplot )
      write(6,*)' iwini:iwend omegamax(Ry)=',iwini,iwend,omegamax


! swj added end


      ifsex(1)   = iopenxx ('SEXU')
      ifsexcore(1) = iopenxx ('SEXcoreU')
      ifsec(1)   = iopenxx ('SECU')
      ifxc(1)    = iopenxx ('XCU')


      call readx   (ifsex(1),50)
      read (ifsex(1),*) nspin,nq,ntq

! swj added
      read (ifsex(1),*)
      read (ifsex(1),*) deltaw
      read (ifsex(1),*) alat
      read (ifsex(1),*) ef

      allocate(omega(ntq,iwini:iwend))
! swj added


      if (nspin .eq. 2 .and. .not. laf) then
        ifsex(2)   = iopenxx ('SEXD')
        ifsexcore(2)   = iopenxx ('SEXcoreD')
        ifsec(2)   = iopenxx ('SECD')
        ifxc(2)    = iopenxx ('XCD')
      endif
      rewind (ifsex(1))

c> output file
      ifqpe(1)   = iopenxx ('QPU')
      iftote(1)  = iopenxx ('TOTE.UP')
      iftote2(1) = iopenxx ('TOTE2.UP')
      if (nspin == 2) then
        ifqpe(2)   = iopenxx ('QPD')
        iftote(2)  = iopenxx ('TOTE.DN')
        iftote2(2) = iopenxx ('TOTE2.DN')
      endif


c loop over spin
      do      is = 1,nspin
        write(6,*) ' --- is=',is
c read dimensions
        call readx   (ifsex(is),50)
        read (ifsex(is),*) nspinx,nqx,ntqx
c      call readx   (ifsex(is),50)
        read (ifsex(is),*)
        read (ifsex(is),*) deltaw
        read (ifsex(is),*) alat
        read (ifsex(is),*) ef
c
        call readx(ifsec(is),50)
        read (ifsec(is),*) nspinc,nqc,ntqc
c
        call readx   (ifxc(is),50)
        read (ifxc(is),*) nspinxc,nqxc,ntqxc
Cstop2rx 2013.08.09 kino        if (nspin .ne. nspinx)  stop 'hqpe: wrong nspin SEx'
        if (nspin .ne. nspinx)  call rx( 'hqpe: wrong nspin SEx')
Cstop2rx 2013.08.09 kino        if (nspin .ne. nspinc)  stop 'hqpe: wrong nspin SEc'
        if (nspin .ne. nspinc)  call rx( 'hqpe: wrong nspin SEc')
Cstop2rx 2013.08.09 kino        if (nspin .ne. nspinxc) stop 'hqpe: wrong nspin vxc'
        if (nspin .ne. nspinxc) call rx( 'hqpe: wrong nspin vxc')
Cstop2rx 2013.08.09 kino        if (nq .ne. nqx)        stop 'hqpe: wrong nq SEx'
        if (nq .ne. nqx)        call rx( 'hqpe: wrong nq SEx')
Cstop2rx 2013.08.09 kino        if (nq .ne. nqc)        stop 'hqpe: wrong nq SEc'
        if (nq .ne. nqc)        call rx( 'hqpe: wrong nq SEc')
Cstop2rx 2013.08.09 kino        if (nq .ne. nqxc)       stop 'hqpe: wrong nq vxc'
        if (nq .ne. nqxc)       call rx( 'hqpe: wrong nq vxc')
Cstop2rx 2013.08.09 kino        if (ntq .ne. ntqx)      stop 'hqpe: wrong ntq SEx'
        if (ntq .ne. ntqx)      call rx( 'hqpe: wrong ntq SEx')
Cstop2rx 2013.08.09 kino        if (ntq .ne. ntqc)      stop 'hqpe: wrong ntq SEc'
        if (ntq .ne. ntqc)      call rx( 'hqpe: wrong ntq SEc')
Cstop2rx 2013.08.09 kino        if (ntq .ne. ntqxc)     stop 'hqpe: wrong ntq vxc'
        if (ntq .ne. ntqxc)     call rx( 'hqpe: wrong ntq vxc')
c
        if(is==1) write(6,*)' ###  readin XCU'
        if(is==2) write(6,*)' ###  readin XCD'
        allocate( itxc(ntq),qxc(3,ntq,nq),eldaxc(ntq,nq),vxc(ntq,nq) )
        call readx (ifxc(is),50)
        read(ifxc(is),*)
        do ip = 1,nq
          do i  = 1,ntq
            read(ifxc(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &      itxc(i),ipxx,isxxx, qxc(1:3,i,ip), eldaxc(i,ip), 
     &      vxc(i,ip)
          enddo
        enddo
c
        if(is==1) write(6,*)' ###  readin SEXU'
        if(is==2) write(6,*)' ###  readin SEXD'
        allocate( itx(ntq), qx (3,ntq,nq),eldax (ntq,nq),sex(ntq,nq) )
        call readx   (ifsex(is),50)
        read(ifsex(is),*)
        do ip = 1,nq
          do i  = 1,ntq
            read(ifsex(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &      itx(i),ipxx,isxxx, qx(1:3,i,ip), eldax(i,ip), 
     &      sex(i,ip)
          enddo
        enddo
c
        if(is==1) write(6,*)' ###  readin SEXcoreU'
        if(is==2) write(6,*)' ###  readin SEXcoreD'
        allocate( sexcore(ntq,nq) )
        call readx   (ifsexcore(is),50)
        call readx   (ifsexcore(is),50)
        read(ifsexcore(is),*)
        do ip = 1,nq
          do i  = 1,ntq
            read(ifsexcore(is),"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &      ixx1,ixx2,ixx3, qxxx1,qxxx2,qxxx3, exxx, sexcore(i,ip)
          enddo
        enddo
c
        if(is==1) write(6,*)' ###  readin SECU'
        if(is==2) write(6,*)' ###  readin SECD'
        allocate( itc(ntq), qc (3,ntq,nq),eldac (ntq,nq)
     &                  ,rsec(3,ntq,nq),csec(3,ntq,nq),zfac(ntq,nq))
        call readx   (ifsec(is),50)
        read(ifsec(is),*)
        rsec=0d0
        csec=0d0
        do ip = 1,nq
          do i  = 1,ntq
          if(nozmode) then
            read(ifsec(is),*)
     &     itc(i),ipxxx,isxxx, qc(1:3,i,ip), eldac(i,ip), 
     &     rsec(2,i,ip),csec(2,i,ip)
           zfac(i,ip)=1.00
           write(*,*) i,ip,csec(2,i,ip)
          else
            read(ifsec(is),*)
c  "(3i5,3d24.16,3x,d24.16,3x,3d24.16,  3x,3d24.16,3x,d24.16)")
     &     itc(i),ipxxx,isxxx, qc(1:3,i,ip), eldac(i,ip), 
     &     rsec(1:3,i,ip),csec(1:3,i,ip),zfac(i,ip)
          endif  
c
c      write(6,*)" itc=",itc(i)
          enddo
        enddo

! swj added
        if(is==1) write(6,*)' ###  readin SEComg.UP'
        if(is==2) write(6,*)' ###  readin SEComg.DN'
        allocate( itcomg(ntq),qcomg(3,nq),eldacomg(ntq,nq),wibz(ip),
     &            rsecomg(iwini:iwend,ntq,nq), csecomg(iwini:iwend,ntq,nq) )
        do ip = 1,nq
          do i  = 1,ntq
          write (cip, "(I4)") ip
          write (ci, "(I4)") itc(i)
          fname_sec= "SEComg_n"//trim(adjustl(ci))//"k"//trim(adjustl(cip))//".UP"
          !print "SEComg_n"//trim(adjustl(ci))//"k"//trim(adjustl(cip))//".UP"
          ifsecomg(is) = iopenxx(fname_sec)
          !print *, fname_sec
            do iw = iwini,iwend
                 read(ifsecomg(is),*) 
     &              iwcomg, itcomg(i), ipcomg, iscomg, qcomg(1:3,ip), wibz(ip),
     &          eldacomg(i,ip), omega(i,iw),  rsecomg(iw,i,ip), csecomg(iw,i,ip)
            !print *, rsecomg(iw,i,ip), csecomg(iw,i,ip)
            enddo
          enddo
          close(ifsecomg(is))
        enddo

        ! starting for interpolation
        if (ninter > 1) then
        !ninter is the number of (interpolated points-1) within each energy bin.

          dwplot_int = dwplot*rydberg()/ninter
  
          iwini_int =  -int( omegamax*rydberg() / dwplot_int )
          iwend_int =   int( omegamax*rydberg()/  dwplot_int )
  
          !print *, iwini_int, iwend_int
          nw = iwend - iwini +1
          !print *, iwini_int, iwend_int
          ifoutsec_int=9300
  
  
          allocate(omega_int(ntq,iwini_int:iwend_int), rsecomg_int(iwini_int:iwend_int,ntq,nq),
     &           csecomg_int(iwini_int:iwend_int,ntq,nq) )

       
          if(is==1) sss='.UP'
          if(is==2) sss='.DN'
          open(ifoutsec_int,file='SEComg_int'//sss)
          do ip = 1, nq
            do i  = 1,ntq
              iw_int = iwini_int 
              do iw = iwini,iwend
                if (iw /= iwend) then
                  do iinter = 0,ninter-1
                    omega_int(i,iw_int) = omega(i,iw) + iinter*dwplot_int
                    !if (iinter == 0) goto 313
                    rsecomg_int(iw_int,i,ip) = polinta(omega_int(i,iw_int),
     &                      omega(i,:),rsecomg(:,i,ip),nw)
                    csecomg_int(iw_int,i,ip) = polinta(omega_int(i,iw_int),
     &                      omega(i,:),csecomg(:,i,ip),nw)
                    write(ifoutsec_int,"(4i5,3f10.6,3x,f10.6,2x,2f16.8,x,3f16.8)")
     &              iw_int,itcomg(i),ip,is, qcomg(1:3,ip), wibz(ip), eldacomg(i,ip), 
     &              omega_int(i,iw_int), rsecomg_int(iw_int,i,ip), csecomg_int(iw_int,i,ip)  !,sumimg

!                  rsecomg_int(iw_int,i,ip) = polinta(omega_int(i,-799),
!     &                        omega(i,:),rsecomg(:,i,ip),nw)
!                  csecomg_int(iw_int,i,ip) = polinta(omega_int(i,-799),
!     &                        omega(i,iw),csecomg_int(iw_int,i,ip),nw)
  313               continue
                    iw_int = iw_int + 1
                  end do
                else 
                    omega_int(i,iw_int) = omega(i,iw) 
                    rsecomg_int(iw_int,i,ip) = rsecomg(iw,i,ip)
                    csecomg_int(iw_int,i,ip) = csecomg(iw,i,ip)
                    write(ifoutsec_int,"(4i5,3f10.6,3x,f10.6,2x,2f16.8,x,3f16.8)")
     &                    iw_int,itcomg(i),ip,is, qcomg(1:3,ip), wibz(ip), eldacomg(i,ip), 
     &                    omega_int(i,iw_int), rsecomg_int(iw_int,i,ip), csecomg_int(iw_int,i,ip)  !,sumimg
                end if
              end do
            end do
          end do
          close(ifoutsec_int)

        end if
!        do iw_int = iwini_int, iwend_int
!           print *, omega_int(1,iw_int) 
!        enddo
        !stop


        if(is==1) sss='.UP'
        if(is==2) sss='.DN'

	  ifoutgreen = 9301
	  ifoutsigma = 9302
        ifoutspect1 = 9303
        ifoutspect2 = 9304

	  open(ifoutsigma,file='Self_energy'//sss)
	  open(ifoutgreen,file='Green_ftn'//sss)
        open(ifoutspect1,file='Spectral_ftn'//sss)
        open(ifoutspect2,file='Spectral_ftn_tot'//sss)



        allocate( eh(ntq,nq), constant_shift(ntq,nq),
     &           green_ftn(iwini:iwend,ntq,nq), sigma(iwini:iwend,ntq,nq),
     &           spect_ftn(iwini:iwend,nq), spect_ftn_tot(iwini:iwend))


        write(ifoutsigma,"(4a6,a25,3x,2a10,2x,a16,x,8a16)")
     &            "# iw","ib","ip","is", "q(1:3,ip)", "Weight",  "Ehk",
     &            "w-Ef",  "SEX", "SEXcore", "Re[SEComg]", "Im[SEComg]",
     &            "Re[Sig(n,k,w)]",  "Im[Sig(n,k,w)]" !,sumimg
        write(ifoutgreen,"(4a6,a25,3x,3a12,2x,a16,x,8a16)")
     &              "# iw","ib","ip","is", "q(1:3,ip)", "Weight",  "Ehk", "Vxc",
     &              "w-Ef", "SEX", "SEXcore", "Re[SEComg]", "Im[SEComg]", 
     &              "Re[G(n,k,w)]",  "Im[G(n,k,w)]"

        if (ninter > 1) then
	    ifoutgreen_int = 9305
	    ifoutsigma_int = 9306
          ifoutspect1_int = 9307
          ifoutspect2_int = 9308

    	    open(ifoutsigma_int,file='Self_energy_int'//sss)
  	    open(ifoutgreen_int,file='Green_ftn_int'//sss)
          open(ifoutspect1_int,file='Spectral_ftn_int'//sss)
          open(ifoutspect2_int,file='Spectral_ftn_tot_int'//sss)
          allocate(green_ftn_int(iwini_int:iwend_int,ntq,nq),
     &           sigma_int(iwini_int:iwend_int,ntq,nq),
     &           spect_ftn_int(iwini_int:iwend_int,nq),
     &           spect_ftn_tot_int(iwini_int:iwend_int))

          write(ifoutsigma_int,"(4a6,a25,3x,2a10,2x,a16,x,8a16)")
     &            "# iw","ib","ip","is", "q(1:3,ip)", "Weight",  "Ehk",
     &            "w-Ef",  "SEX", "SEXcore", "Re[SEComg]", "Im[SEComg]",
     &            "Re[Sig(n,k,w)]",  "Im[Sig(n,k,w)]" !,sumimg
          write(ifoutgreen_int,"(4a6,a25,3x,3a12,2x,a16,x,8a16)")
     &              "# iw","ib","ip","is", "q(1:3,ip)", "Weight",  "Ehk", "Vxc",
     &              "w-Ef", "SEX", "SEXcore", "Re[SEComg]", "Im[SEComg]", 
     &              "Re[G(n,k,w)]",  "Im[G(n,k,w)]"

        endif
  

        do ip = 1,nq
          do i  = 1,ntq
            do iw = iwini,iwend
              !wibz(ip) = 0.4

              eh(i,ip) = eldacomg(i,ip) - vxc(i,ip)
              sigma(iw,i,ip) = sex(i,ip) + sexcore(i,ip) +
     &                                 rsecomg(iw,i,ip) + (0.0,1.0)*csecomg(iw,i,ip)
              write(ifoutsigma,"(4i5,3f10.6,3x,f10.6,3x,f10.6,2x,f16.8,x,8f16.8)")
     &                    iw,itcomg(i),ip,is, qcomg(1:3,ip),  wibz(ip), eldacomg(i,ip),
     &                    omega(i,iw),  sex(i,ip), sexcore(i,ip),
     &                    rsecomg(iw,i,ip),  csecomg(iw,i,ip), sigma(iw,i,ip) !,sumimg

              constant_shift(i,ip) = sex(i,ip) + sexcore(i,ip) + 
     &                             rsec(2,i,ip) + (0.0, 1.0)*csec(2,i,ip) -vxc(i,ip)
              green_ftn(iw,i,ip) = 1./(omega(i,iw) - (eh(i,ip) + sigma(iw,i,ip)
     &                                     - constant_shift(i,ip)))

!              rdenomi = omega(i,iw) - (eh(i,ip)
!     &                    + sex(i,ip) + sexcore(i,ip) + rsecomg(iw,i,ip))
!              cdenomi = -csecomg(iw,i,ip)

              !print *, omega(i,iw), eh(i,ip), sex(i,ip), sexcore(i,ip), rsecomg(iw,i,ip)
              !print *, rdenomi, rdenomi**2
              !print *, " "

              !rgreen_ftn2(iw,i,ip) = rdenomi/(rdenomi**2+cdenomi**2)
              !cgreen_ftn2(iw,i,ip) = -cdenomi/(rdenomi**2+cdenomi**2)

              write(ifoutgreen,"(4i5,3f10.6,3x,f10.6,3x,2f12.6,2x,f16.8,x,10f16.8)")
     &                    iw,itcomg(i),ip,is, qcomg(1:3,ip),  wibz(ip),  eldacomg(i,ip),
     &                    vxc(i,ip), omega(i,iw),  sex(i,ip), sexcore(i,ip),
     &                    rsecomg(iw,i,ip),  csecomg(iw,i,ip),
     &                    green_ftn(iw,i,ip) !,sumimg
!     &                    rgreen_ftn2(iw,i,ip),  cgreen_ftn2(iw,i,ip)


            end do
            write(ifoutgreen,*)
            write(ifoutgreen,*)
            write(ifoutsigma,*)
            write(ifoutsigma,*)
            if (ninter >1) then
              do iw_int = iwini_int,iwend_int
                !wibz(ip) = 0.4
  
                sigma_int(iw_int,i,ip) = sex(i,ip) + sexcore(i,ip) +
     &                                 rsecomg_int(iw_int,i,ip) + (0.0,1.0)*csecomg_int(iw_int,i,ip)
                write(ifoutsigma_int,"(4i5,3f10.6,3x,f10.6,3x,f10.6,2x,f16.8,x,8f16.8)")
     &                    iw_int,itcomg(i),ip,is, qcomg(1:3,ip),  wibz(ip), eldacomg(i,ip),
     &                    omega_int(i,iw_int),  sex(i,ip), sexcore(i,ip),
     &                    rsecomg_int(iw_int,i,ip),  csecomg_int(iw_int,i,ip),
     &                    sigma_int(iw_int,i,ip) !,sumimg


                green_ftn_int(iw_int,i,ip) = 1./(omega_int(i,iw_int) - (eh(i,ip) + sigma_int(iw_int,i,ip)
     &                                     - constant_shift(i,ip)))

!              rdenomi = omega(i,iw) - (eh(i,ip)
!     &                    + sex(i,ip) + sexcore(i,ip) + rsecomg(iw,i,ip))
!              cdenomi = -csecomg(iw,i,ip)

              !print *, omega(i,iw), eh(i,ip), sex(i,ip), sexcore(i,ip), rsecomg(iw,i,ip)
              !print *, rdenomi, rdenomi**2
              !print *, " "

              !rgreen_ftn2(iw,i,ip) = rdenomi/(rdenomi**2+cdenomi**2)
              !cgreen_ftn2(iw,i,ip) = -cdenomi/(rdenomi**2+cdenomi**2)

                write(ifoutgreen_int,"(4i5,3f10.6,3x,f10.6,3x,2f12.6,2x,f16.8,x,10f16.8)")
     &                    iw_int,itcomg(i),ip,is, qcomg(1:3,ip),  wibz(ip),  eldacomg(i,ip),
     &                    vxc(i,ip), omega_int(i,iw_int),  sex(i,ip), sexcore(i,ip),
     &                    rsecomg_int(iw_int,i,ip),  csecomg_int(iw_int,i,ip),
     &                    green_ftn_int(iw_int,i,ip) !,sumimg
!     &                    rgreen_ftn2(iw,i,ip),  cgreen_ftn2(iw,i,ip)

              end do
            end if
          end do
        end do

        write(ifoutspect1,"(3a6,a25,a10,3x,a16,x,2a16)")
     &      "# iw", "ip", "is", "q(1:3,ip)","Weight","w-Ef",  "A(k,w)"
        write(ifoutspect2,"(2a6,3x,a16,x,2a16)")
     &      "# iw", "is", "w-Ef",  "A(w)"

!        do iw = iwini,iwend
!          do i  = 1,ntq+1
!             print *, "iw", iw, "iband", i, "omega", omega(i,iw)
!          end do
!        end do

        spect_ftn = 0d0
        spect_ftn_tot = 0d0

        if (ninter >1) then
          spect_ftn_int = 0d0
          spect_ftn_tot_int = 0d0
        end if

        do iw = iwini,iwend
          do ip = 1,nq
            do i  = 1,ntq
                spect_ftn(iw,ip) = spect_ftn(iw,ip) + 1./3.141592*abs(dimag(green_ftn(iw,i,ip)))
            end do
              spect_ftn_tot(iw) = spect_ftn_tot(iw) + spect_ftn(iw,ip)*wibz(ip)
          end do
          write(ifoutspect2,"(2i5,3x,f16.8,x,2f16.8)")
     &        iw, is, omega(1,iw),  spect_ftn_tot(iw)
        end do

        if (ninter >1) then
          do iw_int = iwini_int,iwend_int
            do ip = 1,nq
              do i  = 1,ntq
                  spect_ftn_int(iw_int,ip) = spect_ftn_int(iw_int,ip) + 1./3.141592*abs(dimag(green_ftn_int(iw_int,i,ip)))
              end do
                spect_ftn_tot_int(iw_int) = spect_ftn_tot_int(iw_int) + spect_ftn_int(iw_int,ip)*wibz(ip)
            end do
            write(ifoutspect2_int,"(2i5,3x,f16.8,x,2f16.8)")
     &        iw_int, is, omega_int(1,iw_int),  spect_ftn_tot_int(iw_int)
          end do
        end if

        do ip = 1,nq
          do iw = iwini,iwend
            write(ifoutspect1,"(3i5,3f10.6,3x,f10.6,f16.8,x,2f16.8)")
     &          iw,  ip, is, qcomg(1:3,ip), wibz(ip), omega(1,iw),  spect_ftn(iw,ip)
          end do

          if (ninter >1) then
          do iw_int = iwini_int,iwend_int
            write(ifoutspect1_int,"(3i5,3f10.6,3x,f10.6,f16.8,x,2f16.8)")
     &          iw_int,  ip, is, qcomg(1:3,ip), wibz(ip), omega_int(1,iw_int),  spect_ftn_int(iw_int,ip)
          end do
          end if
          write(ifoutspect1,*)
          write(ifoutspect1,*)
        end do

        close(ifoutgreen)
        close(ifoutsigma)
        close(ifoutspect1)
        close(ifoutspect2)

        if (ninter >1) then
          close(ifoutgreen_int)
          close(ifoutsigma_int)
          close(ifoutspect1_int)
          close(ifoutspect2_int)
        end if
! swj added

c check that LDA eigenvalues are the same
cmod 93.06.29
c     chev       = compvv (w(ieldax),w(ieldac),ntq*nq,tol,'eLDA x,c')
c     chev       = compvv (w(ieldax),w(ieldaxc),ntq*nq,tol,'eLDA x,xc')
c     chq        = compvv (w(iqx),w(iqc),3*nq,tol,'q-vec x,c')
c     chq        = compvv (w(iqx),w(iqxc),3*nq,tol,'q-vec x,xc')
c-----
c      itxc       = icompvv (w(iitx),w(iitc),ntq,'state-label x,c')
c      itxc       = icompvv (w(iitx),w(iitxc),ntq,'state-label x,xc')
        itxcxxx     = icompvv2 (itx, itc,ntq,' state-label x c')
        itxcxxx     = icompvv2 (itx,itxc,ntq,'state-label x xc')

c calculate quasiparticle energies and write them on file QPE
c      call qpe     (ifqpe(is),iftote,w(iitc),w(iqc),
c     i              w(ieldac),w(ivxc),w(isex),
c     i              w(irsec),w(icsec),w(izfac),eshift,
c     d              ntq,nq,is,
c     o              w(ieqp),w(iwqp))
        call qpe1_qsgw     (ifqpe(is),iftote(is),iftote2(is),itc,qc,
     i              eldac,vxc,sex,sexcore,
     i              rsec,csec,zfac,jin,deltaw,alat,ef,
     d              ntq,nq,is,
     o              eshift0,eshift02,eshtlda)
c     o              w(ieqp),w(iwqp))

        deallocate( itxc,qxc,eldaxc,vxc ,itc, qc ,eldac,
     &                 sexcore ,rsec,csec,zfac,
     &       itx, qx ,eldax,sex, 
     &    itcomg, qcomg, eldacomg, rsecomg, csecomg)

        if (laf) exit
        if (jin .gt. 0) jin = 999999
      end do
      deallocate( omega)

Cstop2rx 2013.08.09 kino      stop ' OK! hqpe '
      call rx0( ' OK! hqpe ')

      end

      subroutine readxxx(ifil)
      character*72 rchar
      integer(4):: n=1000
      do 10 i = 1,n
        read(ifil,'(a72)')rchar
        rchar=trim(adjustl(rchar))
        if(rchar(1:5) .eq. '*****') return
        if(rchar(1:5) .eq. '#####') return
   10 continue
Cstop2rx 2013.08.09 kino      stop 'readx: cannot find the string (gwsrc/wse.f/readxx)'
      call rx( 'readx: cannot find the string (gwsrc/wse.f/readxx)')
      end

