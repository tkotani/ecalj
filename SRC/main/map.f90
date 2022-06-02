  allocate(vec(neh2,2))
  allocate(llen(njorb))
  do i=1,njorb
     read(ifjorb,*) llen(i)
     write(*,*)  "llen(i)=",i,llen(i)
  end do

  !         if(llen==1)
  !         read(ifjorb,*) jlist(i,1),jlist(i,2),jlist(i,3),jlist(i,4),jlist(i,5) !(5,4)=local orbital(1) or not(0)
  !        if(jlist(i,5).gt.2.or.jlist(i,5).le.0) stop "error in jlist(i,5)"
  !      end do
  !      close(ifjorb)

  write(*,*) "[map.F]njorb,NMTO,NM,neh2,nlo,",njorb,NMTO,NM,neh2,nlo

  open(newunit=ifmax,file="emax")

  lcut=200
  allocate(elist(0:lcut),dlist(0:lcut),thlist(0:lcut))
  do n=1,njorb
     if(allocated(jlist)) deallocate(jlist)
     allocate(jlist(llen(n),2))
     do j=1,llen(n)
        read(ifjorb,*)  jlist(j,1), jlist(j,2)
        write(*,*)"n,j,jlist(1),jlist(2=lo)",n,j,jlist(j,1), jlist(j,2)
     end do

     jp =jlist(1,1) !first orbital component
     ilo=jlist(1,2)

     j1=ir(jp,1)
     j2=ir(jp,2)
     j3=ir(jp,3)
     write(*,*) "jpj1j2j3",jp,j1,j2,j3

     emax=-1d0
     emin=1000d0
     do loop=0,lcut
        open(newunit=iffix,file="totvec1")
        !            write(*,*) "totvec1"
        do i=1,neh2
           read(iffix,*) vec(i,1),vec(i,2)
           !              write(*,*) i,vec(i,1),vec(i,2)
        end do
        close(iffix)

        !            theta=(dble(loop)/dble(lcut)*2d0-1d0)*pi*0.5 ![-pi/2,pi/2]
        theta=dble(loop)/dble(lcut)*pi ![0,pi]
        !     theta=dble(loop)/dble(lcut)*pi*0.5-pi*0.5d0
        do i=1,llen(n)
           vec(jlist(i,1),jlist(i,2))=theta
           !               write(*,"(f6.2)",advance="no") vec(jlist(i,1),jlist(i,2))
        end do

        do i=1,neh2
           write(*,"(f7.2)",advance="no") dsin(vec(i,1))
        enddo
        write(*,*)
        do i=1,neh2
           write(*,"(f7.2)",advance="no") dsin(vec(i,2))
        enddo
        write(*,*)
        etr=0d0
        dth=0d0

        do rq=1,nqs         !begin rq loop
           LD=llist(rq)
           jsp=jspl(rq)
           call mkdiff(jp,-ilo,neh2,nlo,ir,vec,NM,LD,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),et,dt)
           etr=etr+et/dble(nqs)
           dth=dth+dt(1,1)/dble(nqs)
        end do !end rq loop
        !     write(*,"(2f6.2)") vec(jorb,1),vec(jorb,2)
        if(etr > emax)then
           emax=etr
           thmax=theta
        end if
        if(etr <= emin) emin=etr
        elist(loop)=etr
        dlist(loop)=dth
        thlist(loop)=theta
     end do
     write(ifmax,*) "emax,theta-max,sin(theta-max)",emax,thmax,dsin(thmax)

     do loop=0,lcut
        etr=(elist(loop)-emin)/(emax-emin)-5d-1
        write(8000+n,"(4f16.10)") dsin(thlist(loop)),etr,dlist(loop),thlist(loop) !,n,jlist(n,:)
        write(7000+n,"(4f16.10)") dcos(thlist(loop)),etr,dlist(loop),thlist(loop) !,n,jlist(n,:)
        write(9000+n,"(4f16.10)") thlist(loop),etr,dlist(loop),thlist(loop) !,n,jlist(n,:)
     end do
  end do


  stop "checkprogram stoped. For return 0, change drawmap .false."

