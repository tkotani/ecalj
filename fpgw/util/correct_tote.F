      implicit none
      integer(4):: i,irow,icol,j,iq,idat,ndat,nq,idd
      character(300) ::iline
      character(10) :: label(20),flag='nul'
      real(8),allocatable:: dmat(:,:)

      character(3):: charnum3
      character(300):: ddd,fff
      integer(4),parameter:: ifin=10

      real(8):: q(3),SIG(4),e1,e2,eigen,e0,qq(3,1000),eee(3)
      integer(4):: is,isx(1000),iqx(1000),nqi,iqi

      real(8):: G(3)=(/0d0,0d0,0d0/)
      real(8):: X(3)=(/0d0,0d0,1d0/)
      real(8):: L(3)=(/-0.5d0,0.5d0,0.5d0/)
      character(200):: labelq(100),lqx,title,ttt
      integer(4):: iqpz,i1,ie,ix,nnn0
      integer(4),parameter:: nnnmx=10000
      
      real(8):: dse0(nnnmx)

      real(4):: qqq
C--- Get title
      do i=1,1000
        read(5,"(a)") iline
        if(iline(1:1) =='!'.or.iline(1:1) =='#'.or.iline(1:1)=='%')  then
          write(6,"(a)") trim(iline)
          cycle
        endif
        read(iline,"(a)") ttt
        exit
      enddo

C--- get qq(:,i) and label(i) i=1,nq
      do i=1,1000
        read(5,"(a)") iline
        if(iline(1:1) =='!'.or.iline(1:1) =='#'.or.iline(1:1)=='%')  then
          write(6,"(a)") trim(iline)
          cycle
        endif
        if(iline(1:1) =='$') exit
        read(iline,*) qq(:,i), labelq(i)
        write(6,"('% ',3f13.5,' ',a)") qq(:,i), trim(labelq(i))
      enddo
      nq=i-1

C--- Get iqx isx to specify eigen
      do i=1,1000
        read(5,"(a)") iline
        if(iline(1:1) =='!'.or.iline(1:1) =='#'.or.iline(1:1)=='%')  cycle
        if(iline(1:1) =='$') exit
        read(iline,*) lqx, isx(i)
        do iq=1,nq
          if(trim(lqx) ==trim(labelq(iq))) iqx(i) =iq
        enddo
        if(i==1) then
          write(6,"('% ',3f9.5,i5,' <--zerolevel')") qq(:,iqx(i)), isx(i)
        else
          write(6,"('% ',3f9.5,i5)") qq(:,iqx(i)), isx(i)
        endif
      enddo
      nqi=i-1

C--- label
      write(6,"(a18,' &',$)") trim(adjustl(ttt))
      do i=2,nqi
         iq=iqx(i)
        if(i==nqi) then
        write(6,"(a7,i2,' \\\\')") trim(labelq(iq)), isx(i)
        else
        write(6,"(a7,i2,'&',$)") trim(labelq(iq)), isx(i)
        endif
      enddo

c--- readin BASE Get dse0
        read(5,"(a)") ddd
        fff = trim(ddd) !ddd(1:index(ddd//' '," ")-1)
c        print *,fff
        write(6,"('%  BASE: TOTE.UP=',a)")
     &    trim(fff)
        open(ifin,file=trim(fff))
        read(ifin,*)ix
        ix=0
        do 
         read(ifin,*,end=1199) q, is,iq, eee !SIG,e1,e2, eigen
d          read(ifin,*,end=1199) qqq !, is,iq, eee !SIG,e1,e2, eigen
c          write(6,*) q
          ix=ix+1
          dse0(ix) = eee(3)-eee(1)
        enddo
 1199   continue
        close(ifin)
        nnn0=ix
c        write(6,*) 'zzzzzzzzzzzzz dat=',nnn0
        if(nnn0>nnnmx) stop 'nnn0>nnnmx'
c---------------------------------------
c        print *,' end of read in base------'
      do 
C--- Filename
        read(5,"(a)",end=9999) ddd
c        write(6,*) ' xxxxxx=',ddd(1:1)
        if(ddd(1:1) =='$') exit
        if(trim(ddd)=='') exit
        i1= index(ddd," ")-1
        title = ddd(1:i1)
        ddd=adjustl(trim(ddd(i1+1:))//' $')
c        write(6,*)' ddd=',trim(ddd)
        iqpz=1
        i = index(ddd," QP=")
c        write(6,*) 'i=',i
        if(i/=0) read(ddd(i+4:),*) iqpz 
        fff = ddd(1:index(ddd," ")-1)
        write(6,"('% ',a,' TOTE.UP=',a,' QP= ',i1)")
     &    trim(title),trim(fff),iqpz

c---get zerolevel
        open(ifin,file=trim(fff))
        read(ifin,*)
        ix=0
        do 
          ix=ix+1 
          read(ifin,*,err=199,end=199) q, is,iq, eee !SIG,e1,e2, eigen
c          write(6,*) q !,is,iq,eee
          eigen = eee(3) - dse0(ix)
          call pickzero( qq(:,iqx(1)), isx(1), e0) !,  q,is,eigen,e0) 
        enddo
  199   continue
c        write(6,*)'e0=',e0
cc          read(ifin,"(a)",end=999) iline
cc          write(6,"(3f8.3,i4,f10.3)") q,is,eigen
c----
        write(6,"(a18,' & ',$)")title

        do iqi = 2,nqi
          rewind(ifin)
          read(ifin,*)
          ix=0
          do 
            ix=ix+1
            read(ifin,"(a)") iline 
            iline=trim(iline)//" 0d0 0d0"
c            write(6,"(a)") trim(iline)
            read(iline,*,end=999) q, is,iq, eee !SIG,e1,e2, eigen
            eigen = eee(3) - dse0(ix)
            if(iqi==nqi) then
              call pickstate(qq(:,iqx(iqi)), isx(iqi), ' \\\\ ',ie) 
            else
              call pickstate(qq(:,iqx(iqi)), isx(iqi), ' & ',ie) 
            endif
            if(ie==1) exit
          enddo
c          if(ix/=nnn0) stop 'ix/=nnn0'
  999   continue
        enddo
        close(ifin)
        write(6,*)
      enddo
 9999 continue
      contains !-------------------------------------------------
      subroutine pickstate(qt,it,  comm,ie) !q,is,eigen,comm)
      character*(*) ::comm 
      integer(4):: it,ie
      real(8):: qt(3)
      ie=0
      if(it/=is)  return
      if(sum(abs(q-qt))>1d-6) return
c      write(6,"(f9.5,' ',a,$)")eigen-e0,comm
      write(6,"(f7.3,' ',a,$)")eigen-e0,comm
c      write(6,"(f6.2,' ',a,$)")eigen-e0,comm
      ie=1
      end subroutine

      subroutine pickzero(qt,it,  e0) !q,is,eigen,e0)
      integer(4):: it,ie
      real(8):: qt(3),e0
      if(it/=is)  return
      if(sum(abs(q-qt))>1d-10) return
      e0=eigen
c      write(6,"(f10.3,' ',a)") eigen,comm
      end subroutine

      end
c             q               state  SEx   SExcore SEc    vxc   ---  dSEnoZ  eLDA   ------  ------  eHF  Z=1  FWHM=2Z*Simg ReS(elda)
c   0.00000  0.00000  0.00000  1  -38.26 -18.46   9.65 -43.60   0.00  -3.46 -18.20   0.00   0.00 -32.31 1.00   0.00000    -47.05946
      
