      program mkinput
      implicit none
      logical :: flg,newline
      character(2) :: atmspec, atmold,atmname(1000)
      character(200) :: cdum(8)
      integer :: access,i,j,k,ifspc,l,ifinp,ib,il,ik,lold,ind,m,NMTO,nsp,anum
     . ,cnum,onum,d,ii,kold,klist(1000),iblist(1000),ibold,ifim,ifix
      integer,allocatable :: num(:),cl(:)

      if(access("atmspc","").ne.0) stop "error! file atmspc is needed as the input!!"
      open(newunit=ifspc,file="atmspc") !main input. check the existence!!
      open(newunit=ifinp,file="mtoctrl")
      open(newunit=ifim,file="meta_mto")
      open(newunit=ifix,file="index_f")

      do i=1,8
         read(ifinp,"(a200)") cdum(i)
      end do
      close(ifinp)
      open(newunit=ifinp,file="mtoctrl",status="replace")
      do i=1,8
         write(ifinp,"(a)") cdum(i)
      end do

      write(ifinp,*) " "
      write(ifinp,"(a)") "<mkhamloc>         !DO NOT leave any spaces on the left of atomic symbol!"
c      write(ifinp,"(es9.2,2f6.1,2i4, 2i6,  es9.2,2i6,a,a)")
c     .     1d-5, -50d0,20d0, 1,0, 30,30,1d-7,100,100
c     .     ," !convergence, en-min, en-max, orb_sym(1=on), spin_sym(0=off), max-iter, "
c     .     ,"max iter(bisection), conv. of bisection,lcut, rcut"
      write(ifinp,"(2es9.2,a)") 1d-5,1d-7, " !convergence criterion of main and bisection part"
      write(ifinp,"(2f6.1,a)") -100d0,20d0, "       !E-min and E-max, for projection"
      write(ifinp,"(i0,2x,i0,a)") 1,0,"               !orbital and spin symmetrization(1=on)"
      write(ifinp,"(i0,2x,i0,a)") 30,30, "             !maximum iteration of main and bisection part"
      write(ifinp,"(i0,2x,i0,a)")
     . 100,100, "           !lcut and rcut(see error message if mkhamloc outputs)"

      
      atmold="PP"
      lold=0
      anum=0
      cnum=0
      onum=0
      kold=0      
      read(ifspc,*) NMTO,nsp
         
      allocate(num(NMTO),cl(NMTO))
      
      do j=1,NMTO
         read(ifspc,*,end=999) ind,atmspec,ib,l,k !main input!
         if(l==3)then
            call rareearth(atmspec,flg)
            if(.NOT.flg) write(ifix,"(i4)") j
            flg=.false.
         end if

         if((atmspec.ne.atmold.or.kold.ne.k).and.onum.ne.0)then
            anum=anum+1
            cl(anum)=onum
            atmname(anum)=trim(adjustl(atmold))
            klist(anum)=kold
            iblist(anum)=ibold
         end if
         onum=onum+1            

         if(l==lold.and.kold==k) m=m+1
         if(l/=lold.or. kold/=k) m=0
         
         if(l==0)then          !s orbital
            num(j)=1
         else if(l==1)then     !p orbital
            num(j)=m+2
         else if(l==2)then     !d orbital
            num(j)=m+5
         else if(l==3)then     !f orbital
            num(j)=m+10
         end if
         write(ifim,"(4i4)",advance="yes") j,ib,k,num(j)
         lold=l
         kold=k
         ibold=ib
         atmold=atmspec
      end do

 999  continue
      anum=anum+1
      cl(anum)=onum
      atmname(anum)=trim(adjustl(atmspec))
      klist(anum)=k
      iblist(anum)=ib

      do i=1,anum
         newline=.false.
         if(i==1)   j=1
         if(i.ne.1) j=cl(i-1)+1
         k=cl(i)

         if(klist(i)==1.or.nsp==2)then
            write(ifinp,"(a2,a3,i4,a3,i3,a3)",advance="no") atmname(i)," : ",iblist(i)," : ",klist(i)," : "
            newline=.true.            
         else if(klist(i)==3.or.nsp==2)then
            write(ifinp,"(a2,a3,i4,a3,i3,a3)",advance="no") atmname(i)," : ",iblist(i)," : ",klist(i)," : "
            newline=.true.
         end if         
         
         do ii=j,k
            if(klist(i)==3)  write(ifinp,"(i4)",advance="no") num(ii) !to be changed, but compromize in current version
            
            if(num(ii).le.4)then
               if(klist(i)==1)write(ifinp,"(i4)",advance="no") num(ii)
            else if(num(ii).le.9)then
               call typicalanion(atmname(i),flg)
               if(flg)cycle
               if(nsp==1.and.klist(i)==1) write(ifinp,"(i4)",advance="no") num(ii)
               if(nsp==2)                 write(ifinp,"(i4)",advance="no") num(ii)
            else !f-orbital: num=10--16
               call rareearth(atmname(i),flg)
               if(flg.and.nsp==1.and.klist(i)==1) write(ifinp,"(i4)",advance="no") num(ii)
               if(flg.and.nsp==2)                 write(ifinp,"(i4)",advance="no") num(ii)
            end if
         end do
         if(newline)write(ifinp,*) " "
      end do

      write(ifinp,*) " "
      write(ifinp,"(a)") "<lmfham>"      
      write(ifinp,"(f6.2,a)") 4.0d0," !projection exponent. (>2.0) tight fitting (~1.0) entangled cases."

      write(ifinp,*) " "
      write(ifinp,"(a)") "<plotmodel>"
      do i=1,nsp
         write(ifinp,"(a,i4,a)") "HamM",i," !Needed only once. Delete HamM after you obtain MTO bands."
         write(ifinp,"(a,i4,a)") "HamS",i," !The number indicates the spin (1:majority, 2:minority)"
         write(ifinp,"(a,i4)") "HamC",i
      end do

      stop "return 0(mkinput)"
      end program mkinput

      subroutine rareearth(name,flg)
      intent(in)  :: name
      intent(out) :: flg
      character(2) name
      logical(4) ::flg

      flg=.false.
c      write(*,*) "4f orbital check for [",name,"]"
      if(name=="La".or.name=="Ce".or.name=="Pr".or.name=="Nd")then
         flg=.true.
      else if(name=="Pm".or.name=="Sm".or.name=="Eu".or.name=="Gd")then
         flg=.true.
      else if(name=="Tb".or.name=="Dy".or.name=="Ho".or.name=="Er")then
         flg=.true.
      else if(name=="Tm".or.name=="Yb".or.name=="Lu")then
         flg=.true.
      end if
      if(flg)then
         write(*,*) "[Lathanoid] 4f-orbitals are considered explicitly for [",name,"]"
         return
      end if
         
c      write(*,*) "5f orbital check for [",name,"]"
      if(name=="Ac".or.name=="Th".or.name=="Pa".or.name=="U")then
         flg=.true.
      else if(name=="Np".or.name=="Pu".or.name=="Am".or.name=="Cm")then
         flg=.true.
      else if(name=="Bk".or.name=="Cf".or.name=="Es".or.name=="Fm")then
         flg=.true.
      else if(name=="Md".or.name=="No".or.name=="Lr")then
         flg=.true.
      end if
      if(flg)then
         write(*,*) "[Actinoid] 5f-orbitals are considered explicitly for [",name,"]"
      else
c         write(*,*) "[Non-f atms] NO f-orbitals are considered explicitly for [",name,"]"        
      end if      
      return
      end subroutine rareearth

      subroutine typicalanion(name,flg)
      intent(in)  :: name
      intent(out) :: flg
      character(2) name
      logical(4) ::flg

      flg=.false.
      if(name=="N".or.name=="O".or.name=="F".or.name=="Cl".or.name=="S")then
         flg=.true.
      end if
      return
      end subroutine typicalanion

