C-  Prepare spectrum function A(omega) and so.
C Need 
C     SEComg.UP (DN)
C
C     TOTE.UP  (DN) : only for oneshot=1 case
C
C     QIBZ for dos plot mode.
C --------------------------------------------------
      implicit none
      character(3)  :: charnum3,sss
      real(8),allocatable:: omg(:),wibz(:),qibz(:,:),eig(:,:),ez1(:,:)
     &  ,dos(:),omgx(:),dose(:,:),summ(:)
      complex(8):: sigm(:,:,:)
      character(30) fname
      character(130) fff
      character(2) fin
      integer(4):: ifiqibz,iqibz,nwmin,nwmax,ifoutsec,iw_r,itq_r,ip,is_r
     & ,iw,itq ,is,ntq ,ib_r,iq_r,nband,iq,ib,nwin,ndat,nqibz
     & ,iwx,nwx,im,ifx,oneshot,nbandx,nqibzx,iftote,ib_rx
      real(8):: q_r(3),q_rx(3)=-9999,eig_r,omg_r,egw,sse,ein,polinta,
     & s1,s2,sumchk,ef,sumchk2,eee,egw1,eeee
      complex(8):: sigm_r,ss,dd
      logical :: readqibz=.false.

c --- set parameter 
      integer(4),parameter:: nbandmx =200, nwplt=5000, nx=100,
     &     multi= 20
      real(8),parameter:: eps=0.001d0
      nqibz=100
C----------------------------------------------------------------
      write(6,*) "QSGW:0  or one-shot:1 ?"
      read(5,*) oneshot

C... Readin QIBZ
      ifiqibz = 6661
      open (ifiqibz, file='QIBZ',status='old',err=1011) !write q-points in IBZ.
      read(ifiqibz,*,err=1011,end=1011) nqibz
      readqibz=.true.
      allocate(qibz(1:3,nqibz),wibz(nqibz))
      do iqibz = 1,nqibz
        read(ifiqibz,"(3d24.16,3x,d24.16)") qibz(1:3,iqibz),wibz(iqibz)
c       write(6,"(3d24.16,3x,d24.16)") qibz(1:3,iqibz),wibz(iqibz)
      enddo
 1011 continue
      print *,' readqibz=',readqibz,nqibz
      ifoutsec=1015
      nband = nbandmx
      nwmin=0
      nwmax=0
      allocate( sigm(-nwplt:nwplt,nband,nqibz), omg(-nwplt:nwplt),
     & eig(nbandmx,nqibz), ez1(nbandmx,nqibz))
      eig=1d10

C-------------------------
      do 2100 is=1,2
        if(is==1) sss='.UP'
        if(is==2) sss='.DN'
        fname='SEComg'//sss
        open(ifoutsec, file=fname,status='old',err=2100)
        do  
          read(ifoutsec,*,end=2010)
     &    iw_r, itq_r, ip, is_r, q_r(1:3),  eig_r,  omg_r,  s1,s2
!          iw   band    q   is    q_r      
c          write(6,"(3d13.5)") q_r
c          write(6,"(3d13.5)") qibz(:,ip)
          if(readqibz) then
            if(sum(abs(q_r - qibz(:,ip))) >1d-6) then
c            write(6,"(3d13.5)") q_r
c            write(6,"(3d13.5)") qibz(:,ip)
              readqibz=.false.
c             stop 'q_r /= qibz(:,ip)'
            endif
          endif
          if(iw_r<nwmin) nwmin=iw_r
          if(iw_r>nwmax) nwmax=iw_r
          omg (iw_r) = omg_r
          eig(itq_r,ip)= eig_r
c        write(6,*) 'uuu:', iw_r,itq_r,ip,nwplt,s1,s2,nband,nqibz
          sigm(iw_r,itq_r,ip) =dcmplx(s1,s2)

C... If oneshot, get eQP(Z=1) from TOTE -----------
          if(oneshot==1) then
            if(sum(abs(q_rx-q_r)<1d-6).and.ib_rx==itq_r) then
              goto 1020
            endif
            iftote=1016
            fname='TOTE'//sss
            open(iftote, file=fname)
            read(iftote,*) nqibzx,nbandx
c            write(6,*) 'zzzzzzzzzzz',nqibzx,nbandx
            do iq=1,nqibzx
            do ib=1,nbandx
              read(iftote,*) q_rx, ib_rx, iq_r, eee,egw,egw1
c              write(6,"('xxx1:',3d13.6,i5)") q_rx, ib_rx
c              write(6,"('xxx2:',3d13.6,i5)")  q_r, itq_r
              if(sum(abs(q_rx-q_r)<1d-6).and.ib_rx==itq_r) then
c               write(6,"('xxxxxxxxxxx2:',3d13.6,i5)") 
                ez1(itq_r,ip) =egw1
                goto 1019
              endif
            enddo
            enddo
c            write(6,"('vvv2:',3d13.6,i5)") q_r, itq_r
            stop 'can not find eQP in TOTE.*'
 1019       continue
            close(iftote)
 1020       continue
          endif
C -----------------------------------------------------
        enddo
 2010   continue
        close(ifoutsec)  

C dos calculation  !   write (6,*) nwmin,nwmax
        ndat=nwmax-nwmin+1
        allocate(dos( ndat*multi),dose( ndat*multi,2),omgx( ndat*multi),
     &  summ( ndat*multi) )
        dos = 0d0
        dose = 0d0
        do ip =1, nqibz
        do itq=1, nband   
          if( eig(itq,ip)==1d10 ) cycle
          fff = 'SEComg'//'_ib'//charnum3(itq)
     &          //'_ip'//charnum3(ip)//trim(sss)

C for guuplot
          open(16,file= trim(fff)//'.gp')
          write(16,"(a)")"#!/usr/bin/gnuplot -persist"
          write(16,"(a)")"set title  '--- test spectrum ---'"
          write(16,"(a)")"set ylabel '---plot------'"
          write(16,"(a)")"#set xrange [-25:25]"
          write(16,"(a)")"#set yrange [-10: 10]"
          write(16,"(a)") 'plot '//'"'//trim(fff)//'"'
     &     //' using ($1):($2) with l, \\'
          write(16,"(a)")
     &    ' '//'"'//trim(fff)//'"'//' using ($1):($3) with l, \\'
          write(16,"(a)") 
     &    ' '//'"'//trim(fff)//'"'//' using ($1):(abs($4)) with l, \\'
          write(16,"(a)") '  '//'"'//trim(fff)//'"'
     &      //' using ($1):($4)/(($2-$3)**2+($4)**2) with l'
          write(16,"(a)") 'set term postscript enhanced color'
          write(16,"(a)") 'set output "pl.eps"'
          write(16,"(a)") 'replot'
          close(16)

          ifoutsec = 9300
          open(ifoutsec,file=trim(fff))
          write(6,*)
          write(6,*) trim(fff)
! Z=1 case : 
!  E_qp = eig +  (sigm_x +sigm(eig) - vxc)
!  Therefore
!     Denominator = \omega - (E_qp + sigm(omega) -sigma(E_qp) )
          sse= polinta (eig(itq,ip), omg(nwmin:nwmax), 
     &           dreal(sigm(nwmin:nwmax,itq,ip)) , nwmax-nwmin+1)
          iwx =0
          sumchk =0d0
          sumchk2=0d0
          do iw= nwmin,nwmax-1
          do im= 0,multi-1
            iwx= iwx+1
            omgx(iwx)= omg(iw) + (omg(iw+1)-omg(iw))*im/dreal(multi) 
            if(oneshot==1) then
              ein = omgx(iwx) + eig(itq,ip)- ez1(itq,ip) 
            else
              ein = omgx(iwx)  
            endif

! polinta is is interpolation routine to get value at ein.
c      f(ein) = polinta(ein,  x(1..n), f(1...n), n)
            ss = dcmplx( 
     &          polinta(ein, omg(nwmin:nwmax), 
     &               dreal(sigm(nwmin:nwmax,itq,ip)), nwmax-nwmin+1),
     &          polinta(ein, omg(nwmin:nwmax), 
     &               dimag(sigm(nwmin:nwmax,itq,ip)), nwmax-nwmin+1))

            if(oneshot==1) then
              eeee = ez1(itq,ip)
            else
              eeee = eig(itq,ip)
            endif
            dd = dcmplx(dreal( omgx(iwx) - (eeee + ss - sse)),
     &              abs(dimag(      ss)) + eps  )! eps is smearing
            if(readqibz) 
     &      dos(iwx) = dos(iwx) + wibz(ip)/2d0 * abs(dimag(1d0/dd))
            sumchk = sumchk +  abs(dimag(1d0/dd))
            summ(iwx) = sumchk

            write(ifoutsec,"(' ',15d15.7)") 
     &      omgx(iwx),  omgx(iwx)- eeee, (ss - sse), eeee, summ(iwx)

c           dd = dcmplx(dreal(omgx(iwx) - ez1(itq,ip)), eps  )  !smearing
c           dose(iwx,1) = dose(iwx,1) + wibz(ip)/2d0 * abs(dimag(1d0/dd))
            dd = dcmplx(dreal(omgx(iwx) - eeee), eps  )  !smearing
            if(readqibz) 
     &      dose(iwx,2) = dose(iwx,2) + wibz(ip)/2d0 *abs(dimag(1d0/dd))
            sumchk2 = sumchk2 +  abs(dimag(1d0/dd))
          enddo
          enddo
          write(6,"(' sss sumchk=',2i5,2d14.6)") itq,ip,sumchk,sumchk2
          close(ifoutsec)
          nwx=iwx
        enddo
        enddo

        if(readqibz) then
         ifx=106
         if(is==1) open(ifx,file='dos.dat.up')
         if(is==2) open(ifx,file='dos.dat.dn')
         do iwx=1,nwx
         write(ifx,"(12d23.15)")omgx(iwx),dos(iwx),dose(iwx,2)
         enddo  
         close(ifx)
        endif
        deallocate( dos,dose,omgx, summ)

 2100 continue
      stop ' end of dossig'
      end



      double precision function polinta(x,xa,ya,n)
c----------------------------------------------------------------------
c     Given arrays xa and ya, each of length n and given value x,
c     this function returns a value polint. If p(x) is the polynominal
c     of degree ndg such that p(xa(i))=ya(i), i=ns,..,ns+ndg then
c     the returned value polint=p(x). ns is obtained by hunting.
c     See Numerical Recipes
c     coded by H.Akai
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (ndgmx=4, nmx=ndgmx+1)
      dimension xa(n),ya(n),c(nmx),d(nmx)
      logical ascnd
      save jlo
      data jlo/0/ , small/1d-30/
      ndg=min(ndgmx,n-1)
      ndt=ndg+1
      ascnd=xa(n) .gt. xa(1)
      if(jlo .le. 0 .or. jlo .gt. n) then
      jlo=0
      jhi=n+1
      go to 30
      endif
      inc=1
      if(x .gt. xa(jlo) .eqv. ascnd) then
   10 jhi=jlo+inc
      if(jhi .gt. n) then
      jhi=n+1
      else if(x. gt. xa(jhi) .eqv. ascnd) then
      jlo=jhi
      inc=inc+inc
      go to 10
      endif
      else
      jhi=jlo
   20 jlo=jhi-inc
      if(jlo .lt. 1) then
      jlo=0
      else if(x .lt. xa(jlo) .eqv. ascnd) then
      jhi=jlo
      inc=inc+inc
      go to 20
      endif
      endif
   30 if(jhi-jlo .ne. 1) then
      jm=(jhi+jlo)/2
      if(x .gt. xa(jm) .eqv. ascnd) then
      jlo=jm
      else
      jhi=jm
      endif
      go to 30
      endif
      nlo=max(1,jlo-ndg/2)
      nhi=min(n,nlo+ndg)
      nlo=nhi-ndg
      if(jlo .eq. 0) then
      ns=1
      else if(jlo .eq. n) then
      ns=ndt
      else if(abs(x-xa(jlo)) .lt. abs(x-xa(jhi))) then
      ns=jlo-nlo+1
      else
      ns=jhi-nlo+1
      endif
      do 40 i=1,ndt
      ii=nlo+i-1
      c(i)=ya(ii)
   40 d(i)=ya(ii)
      polint=ya(nlo+ns-1)
      ns=ns-1
      do 60 m=1,ndg
      do 50 i=1,ndt-m
      ii=nlo+i-1
      ho=xa(ii)-x
      hp=xa(ii+m)-x
      w=c(i+1)-d(i)
      den=ho-hp
c
c     an error can occur if two xa's are identical
      if(abs(den) .lt. small) then
      write(6,1000)
 1000 format('   ***wrn in polint...data error')
      stop
      endif
c
      den=w/den
      d(i)=hp*den
   50 c(i)=ho*den
      if(2*ns .lt. ndt-m) then
      dy=c(ns+1)
      else
      dy=d(ns)
      ns=ns-1
      endif
   60 polint=polint+dy
c takao
      polinta=polint
      return
      end
      character(3) function charnum3(num)
      integer(4) ::num
      charnum3 = 
     &           char(48+mod(num/100,10))//
     &           char(48+mod(num/10,10))//
     &           char(48+mod(num,10))
      end
