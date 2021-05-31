      module m_qplist
      integer,parameter,private:: nsymlmax=100

      real(8),allocatable,protected:: qplist(:,:)
      integer,protected:: nkp,nsyml=0
      integer,protected:: nqp_syml(nsymlmax),nqps_syml(nsymlmax),nqpe_syml(nsymlmax)
      integer,protected:: nqp2n_syml(nsymlmax)
      real(8),allocatable,protected:: xdatt(:)
      real(8),protected:: dqsyml(nsymlmax),etolv,etolc
      character*20,protected ::labeli(nsymlmax),labele(nsymlmax)
      
      contains
      subroutine m_qplist_init(plbnd) !,readeferm,eferm)
      use m_lmfinit,only: nspc,stdo
      use m_mkqp,only: rv_p_oqp,rv_a_owtkp
      use m_MPItk,only: master_mpi,mlog
      use m_mkqp,only: nkabc=> bz_nabc,bz_nkp
      use m_lattic,only: qlat=>lat_qlat
      use m_ext,only: sname
      intent(in)::             plbnd
c      intent(out)::                  readeferm,eferm
      integer:: nqp2_syml(nsymlmax),nqp2s_syml(nsymlmax),nqp2e_syml(nsymlmax)
      logical:: masslineon(nsymlmax)
      integer:: plbnd,ifqplist,nkk1,nkk2,nkk3,ik1,ik2,ik3,iq,ifi,ifisyml,isyml,i,ierr
      logical:: cmdopt0,fullmesh,PROCARon,readeferm,fsmode
      integer:: ikp,ifbnd,nsymln,onesp
      real(8):: qps_syml(3,nsymlmax), qpe_syml(3,nsymlmax),rq,eferm
      real(8)::  totxdatt
      integer:: ifile_handle, fxsts,ii
      character(512)::schar
      character(3):: charnum3
      character(50)::infoq
      nkk1=nkabc(1)
      nkk2=nkabc(2)
      nkk3=nkabc(3)
      fullmesh = cmdopt0('--fullmesh').or. cmdopt0('--fermisurface') !full mesh stop just after do 2010 iq loop.
      fsmode   = cmdopt0('--fermisurface') !full mesh stop just after do 2010 iq loop.
      PROCARon = cmdopt0('--mkprocar')
      readeferm=.false.
      if(allocated(qplist)) deallocate(qplist)
      if(master_mpi .and. plbnd==0) then
         nkp=bz_nkp
         allocate(qplist(3,nkp))
         qplist= rv_p_oqp !call dcopy(3*nkp, rv_p_oqp,1,qplist,1)
         open(newunit=ifqplist,file='QPLIST.IBZ')
         do iq=1,nkp
            write(ifqplist,"(i5,3f23.15,3x,3f23.15)")iq, qplist(:,iq), rv_a_owtkp(iq)/2d0
         enddo 
         close(ifqplist)
      elseif(master_mpi) then
         if(cmdopt0('--onesp') .and. nspc==1) onesp = 1
         if(fullmesh) then
            nkp = nkk1*nkk2*nkk3
            allocate(qplist(3,nkp))
            if(fsmode) then !!! Fermi surface version for xcrysden
               iq=0
               do ik1=1,nkk1    !ordering is differnt from procaron case...
                 do ik2=1,nkk2
                   do ik3=1,nkk3
                     iq=iq+1
                     qplist(:,iq)   =  qlat(:,1)*dble(ik1-1)/(nkk1-1)
     &                 +   qlat(:,2)*dble(ik2-1)/(nkk2-1)
     &                 +   qlat(:,3)*dble(ik3-1)/(nkk3-1)
                   enddo
                 enddo
               enddo
            endif   
!! pdos (--mkprocar and --fullmesh)
            if(procaron) then !to fit to the tetirr.F
               iq=0
               do ik3=1,nkk3
                 do ik2=1,nkk2
                   do ik1=1,nkk1
                     iq=iq+1
                     qplist(:,iq) =  qlat(:,1)*dble(ik1-1)/nkk1
     &                 +   qlat(:,2)*dble(ik2-1)/nkk2
     &                 +   qlat(:,3)*dble(ik3-1)/nkk3
                   enddo
                 enddo
               enddo
            endif
!! syml direct read for plbnd mode. See "call writeband" below. feb2015 
         else
c            readeferm=.false.
c            if(master_mpi) write(stdo,*)' --- Readin efermi.lmf --- '
c            open(newunit=ifi,file='efermi.lmf',status='old',err=1012)
c            read(ifi,*,err=1012) eferm
c            readeferm=.true.
c            close(ifi) 
c            goto 1013
c 1012       continue
c            call rx('No efermi.lmf!: Copy it, or run lmf-MPIK (sc mode) to get efermi.lmf.')
c 1013       continue

!! --- example of syml file ---
! ndiv qleft(1:3) qright(1:3) llabel rlabel  ndiv2 ninit2 nend2 etolv(Ry) etolc(Ry) 
! 5  0 0 0   .5 .5  .5        Gamma  L       1025  1  16     0.1      0.01 
! 5  0 0 0    1.  0  0        Gamma  X    
! 5  0 0 0   .75 .75 0        Gamma  K       1025  8  38     0.1      0.01 
!! As this shows you can add or not add line after ndiv2. These are for highly resolved calculations.
            write(6,*)' --- Readin syml file --- '
            open(newunit=ifisyml,status='old',file='syml.'//trim(sname))
            nsyml=0
            nsymln=0
            nqp2_syml=0
            do 
              if(nsyml+1>nsymlmax) call rx('bndfp: Enlarge nsymlmax')
              read(ifisyml,"(a)",end=1015) schar
              if(len(trim(schar))==0.or.schar(1:1)=='#' .or. schar(1:1)=='!' .or. schar(1:1)=='%') cycle !comment line
              read(schar,*,err=1014,end=1014) 
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         labeli(nsyml+1),labele(nsyml+1), nqp2_syml(nsyml+1),nqp2s_syml(nsyml+1),nqp2e_syml(nsyml+1),etolv,etolc
              masslineon(nsyml+1)=.true.
              nqp2n_syml(nsyml+1)= nqp2e_syml(nsyml+1)-nqp2s_syml(nsyml+1)+1
              nqps_syml(nsyml+1)=1
              nqpe_syml(nsyml+1)=nqp_syml(nsyml+1)
              write(6,"(' ',i4,3f9.4,' ',3f9.4,' ',a,' ',a,'  Massl:div,init,end=',3i5)")
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         trim(labeli(nsyml+1)),trim(labele(nsyml+1)),
     &         nqp2_syml(nsyml+1),nqp2s_syml(nsyml+1),nqp2e_syml(nsyml+1)
              goto 1025

 1014         continue
              read(schar,*,err=1015,end=1015) 
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         labeli(nsyml+1),labele(nsyml+1) 
              masslineon(nsyml+1)=.false.
              nqps_syml(nsyml+1)=1
              nqpe_syml(nsyml+1)=nqp_syml(nsyml+1)
              write(6,"(' ',i4,3f9.4,' ',3f9.4,' ',a,' ',a)")
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         trim(labeli(nsyml+1)),trim(labele(nsyml+1))
              nsymln = nsymln+1
              nqp2n_syml(nsyml+1)= 0
 1025         continue
              if(nqp_syml(nsyml+1)==0) exit
              nsyml = nsyml + 1
            enddo
 1015       continue
            close(ifisyml)
            nkp = sum(nqp_syml(1:nsyml)+nqp2n_syml(1:nsyml))
            allocate(qplist(3,nkp))
            if(allocated(xdatt)) deallocate(xdatt)
            allocate(xdatt(nkp))
            totxdatt=0d0
            ikp=0
            do isyml=1,nsyml
              dqsyml(isyml) = dsqrt(sum((qpe_syml(1:3,isyml) -qps_syml(1:3,isyml))**2))
              do i=1,nqp_syml(isyml)+nqp2n_syml(isyml)
                ikp= ikp+1
                if(i<=nqp_syml(isyml)) then
                  rq = dble(i-1)/(nqp_syml(isyml)-1)
                else  
                  ii= i-nqp_syml(isyml)-1
                  rq = dble(nqp2s_syml(isyml)-1+ii)/(nqp2_syml(isyml)-1)
                endif
                qplist(:,ikp)= (1d0-rq)*qps_syml(1:3,isyml) +rq*qpe_syml(1:3,isyml)
                xdatt(ikp) = totxdatt + dqsyml(isyml)*rq
              enddo  
              totxdatt = totxdatt + dqsyml(isyml)
            enddo  
            write(6,"('nsyml nkp=',3i5)") nsyml,nkp
         endif
         if (nkp <= 0) call rx('bndfp: nkp<=0') ! quit if nkp==0
      endif
      
!! broadcase nkp and qplist
      call mpibc1( nkp,1,2,mlog,'bndfp','nkp')
      call mpibc1( onesp,1,2,mlog,'bndfp','onesp')
      if(.not.master_mpi) allocate(qplist(3,nkp))
      call mpibc1 ( qplist, 3 * nkp , 4 , .false. , 'bndfp' , 'qp'  )
      if (master_mpi) then
         open(newunit=ifqplist,file='QPLIST')
         print *,'-------- qplist --------',nsyml
         iq=0
         do isyml=1,nsyml
         do i=1,nqp_syml(isyml) + nqp2n_syml(isyml)
            iq=iq+1
            infoq=''
            if(i==1) infoq=' <-- isyml= '//charnum3(isyml)
            if(i==nqp_syml(isyml)+1) infoq=' <-- isyml Mass= '//charnum3(isyml)
            write(6,"(i5,3f8.3,' ',a)")iq,qplist(:,iq),trim(infoq)
            write(ifqplist,"(i5,3f23.15,x,f12.6,' ',a)")iq,qplist(:,iq),xdatt(iq),trim(infoq)
         enddo
           write(ifqplist,*)
         enddo
         close(ifqplist)
      endif 
      if(procaron.and.nsyml==0.and.master_mpi) then !xdatt is dummy
        if(allocated(xdatt)) deallocate(xdatt)
        allocate(xdatt(nkp))
        xdatt=0d0
      endif
      if(PROCARon) then    
        if(.not.master_mpi) allocate(xdatt(nkp))
        call mpibc1( xdatt, nkp , 4 , .false. , 'bndfp' , 'xdatt'  )
      endif  
      end subroutine
      
      end module
