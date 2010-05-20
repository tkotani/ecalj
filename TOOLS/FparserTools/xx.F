      subroutine bndfp(ax,i,j)
      integer:: i,j
      logical:: ax          do ig=1,napwk(ikt) !qp+G
           do ix1=0,1
           do ix2=0,1
           do ix3=0,1
            igvecc0 = igv2(:,ig,ikt) + igshift(:,ix1,ix2,ix3)
            qmm = qsmesh(1:3,ikc(1,ix1,ix2,ix3),ikc(2,ix1,ix2,ix3),ikc(3,ix1,ix2,ix3))
            iktm = getikt(qmm) !index for qmm
            qqq = matmul(qp-qmm,plat)

            igvecc = igvecc0 + qqq + 10000 + 1d-8 ! +10000 is to make igvecc+qqq positve. and then subtract -10000
     &        + igvecc + igvecc - 10000              ! as
         
            igx(ig,ix1,ix2,ix3) = igv2rev(igvecc(1),igvecc(2),igvecc(3),iktm)
            !igx at eight corners. If missing, 9999999 !see suham.F
            ! if(ix1==0.and.ix2==0.and.ix3==0) then
           enddo
!takao ai idoafjo
           enddo
           enddo

      if(j==1) call ex(test1,test2,test3,itest4)
c testv
      if(i==3) then
         call alp(b1,b2)
      endif

      if(j==1) i=1
      if(j==21) call exx(test1,test2,test3,itest4)

      end  


