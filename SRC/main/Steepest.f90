  if(jorb <= jskip) goto 2021

  lold=l_table(j1)
  ibold=ib_table(j1)
  if(lold==0)write(*,*) "steepest descent for s orbital of atm",ibold
  if(lold==1)write(*,*) "steepest descent for p orbital of atm",ibold
  if(lold==2)write(*,*) "steepest descent for d orbital of atm",ibold
  if(lold==3)write(*,*) "steepest descent for f orbital of atm",ibold

  lenlo=1
  do jdum=jp+1,jp+7
     if(l_table(ir(jdum,1))==lold .AND. ib_table(ir(jdum,1))==ibold)then
        lenlo=lenlo+1
        !            write(*,*) "jdum,lenlo",jdum,lenlo
     else
        jskip=jp-1+lenlo
        exit
     end if
  end do
  if(lenlo /= 1 .AND. lenlo /= 3 .AND. lenlo /= 5 .AND. lenlo /= 7) stop "error! lo in Steepest.F"
  write(*,*) "jskip",jskip

  if( .NOT. allocated(diff))allocate(diff(lenlo,2))

  vM(:,:)=totvec(:,:,jsp)
  do loop=1,400
     diff=0d0
     !         write(*,*) "jp,loop=",jp,loop,jskip
     do rq=1,nqs
        call mkdiff(jp,jskip,neh2,nlo,ir,vM,NM,LD,ovlM(1:LD,1:LD,rq),damM(1:LD,1:LD,rq),et,dt)
        diff(1:lenlo,1:2)=diff(1:lenlo,1:2)+dble(dt(1:lenlo,1:2))
     end do
     diff=diff/dble(nqs)
     !         dum=sum(diff(:,1))
     !         diff(:,1)=dum/dble(lenlo)
     !         dum=sum(diff(:,2))
     !         diff(:,2)=dum/dble(lenlo)
     vM(jp:jskip,1:2)= vM(jp:jskip,1:2)+diff(1:lenlo,1:2)*0.2d0
     if(maxval(abs(diff)) <= 0.001d0)then
        exit
     end if
  end do
  write(*,*) "LO loop=",loop

  dum=sum(vM(jp:jskip,1))
  vM(jp:jskip,1)=dum/dble(lenlo)
  dum=sum(VM(jp:jskip,2))
  vM(jp:jskip,2)=dum/dble(lenlo)

  do j=1,2
     if(j==1)write(*,*) "---------------normal orbital-------------------"
     if(j==2)write(*,*) "---------------local orbital-------------------"
     do i=1,lenlo
        write(*,"(i4,2f15.6,es12.2)") i,dsin(vM(i,j)),dsin(vM(i,j)-totvec(i,j,jsp)),diff(i,j)
     end do
  end do
  write(*,*) "---------------------------------------------"

  !      if(jsp.ne.1) stop "error in jsp![Steepest.F]"
  totvec(jp:jskip,1:2,1)=vM(jp:jskip,1:2)
  write(*,*) "exit steepest descent (goto 2021, not 2020)"
  goto 2021

