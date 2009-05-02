      integer(4) function checksymlon(qp,qq1,qq2,nline,qlat,plat,ifi)
C- check qp+G is on the line between qq1 to qq2 -------------
      implicit none
      logical chksyon,check
      integer(4):: nline,iline,ifi
      real(8)   :: qp(3),qq1(3,nline),qq2(3,nline)
     & ,qlat(3,3),plat(3,3),gvec(3)
      checksymlon=0
      do iline=1,nline
        check
     &  = chksyon(qp-qq1(1:3,iline),qp-qq2(1:3,iline),qlat,plat,gvec)
        if(check) then
          checksymlon = iline
          if(ifi/=0) then
            write(ifi,*)
            write(ifi,"('line init ',3f10.5,' iline=',i3)") 
     &      qq1(:,iline),iline
            write(ifi,"('  q-G   = ',3f10.5,' gvec='3f10.5)") 
     &      qp-gvec,gvec
            write(ifi,"('line end  ',3f10.5,3x,3f10.5,3x,3f10.5)") 
     &      qq2(:,iline)
          endif
          exit
        endif
      enddo
      end
c------------
      logical function chksyon(qq1,qq2,qlat,plat,gvec)
c Check G is between qq1 and qq2 or not.
      implicit none
      integer(4)::i,ix,iy,iz,nnmax(3),nnmin(3)
      real(8)  :: qlat(3,3),plat(3,3),qq1(3),qq2(3)
     & ,ee(3),vecn(3),aa1(3),aa2(3),dnorm,proj,gvec(3)
      chksyon=.false.
      do i=1,3
        aa1(i)  = sum(plat(:,i)*qq1)
        aa2(i)  = sum(plat(:,i)*qq2)
        nnmax(i)= nint(max(aa1(i),aa2(i))-0.5d0+1d-6)
        nnmin(i)= nint(min(aa1(i),aa2(i))+0.5d0-1d-6)
      enddo
      ee = (aa2-aa1)/sqrt(sum((aa1-aa2)**2))
      do ix = nnmin(1),nnmax(1)
        do iy = nnmin(2),nnmax(2)
          do iz = nnmin(3),nnmax(3)
            vecn = (/ix,iy,iz/) - aa1
            dnorm = sqrt(sum(vecn**2))
            proj  = abs(sum(vecn*ee))
            if(abs(proj-dnorm)<1d-6) then
              gvec= qlat(:,1)*ix+qlat(:,2)*iy+qlat(:,3)*iz
cccccccccccc
c          write(335,*)
c          write(335,"('aa1 aa2 ',3f10.5,3x,3f10.5)") aa1,aa2
c          write(335,"('ix iy iz=',3x,3i5)") ix,iy, iz
c          write(335,"('dnorm proj=',3x,3f10.5)") dnorm, proj
cccccccccccc
              chksyon=.true.
              goto 999
            endif
          enddo
        enddo
      enddo
 999  continue
      end

