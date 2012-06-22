      subroutine iqindx2sym(qin, ginv,sym,ngrp,qq,nq, iqindx,qu)
C- Find index as q=qq(:,iq) with modulo of premitive vector. symops version.
C- Find index as q=qq(:,iq) with modulo of premitive vector.
Ci ginv is the inverse of plat (premitive translation vector).
      implicit none
      integer(4):: nq,i_out, iq,iqx,iqindx,ngrp,ig !,saveiq
      real(8):: q(3),qq(3,nq),ginv(3,3),qx(3),qu(3),qin(3)
      integer(4),save :: isave=0  !accelation to find iqindx
      real(8)::sym(3,3,ngrp)
c      print *, ' iqindx2 q=',q
c      if(isave>nq.or.saveiq()==0) isave=0
      do ig=1,ngrp
        q=matmul(sym(:,:,ig),qin)
        if(isave>nq) isave=0
        do iqx = isave+1,isave+nq
          if(iqx > nq) then
            iq = iqx -nq
          else
            iq = iqx
          endif
          call rangedq(matmul(ginv,q-qq(:,iq)), qx)
          if(sum(abs(qx))< 1d-6) then
            iqindx=iq
            isave =iq
            qu=qq(:,iq)
c          print *,"isave=",isave
            return
          endif
        enddo
      enddo
! error exit
      write(6,"(' q  =  ',3f13.5,' ginv*q=',3f13.5)")q, matmul(ginv,q)
      write(6,"(' iq    qq    ginv*qq     qq-ginv*qq     err  ')")
      do iq = 1,nq
        call rangedq(matmul(ginv,q-qq(:,iq)), qx)
        write(6,"(i3,3f13.5,' | ',3f13.5,' | ',3f13.5,' diff= ',d13.6)")
     &       iq, qq(1:3,iq), matmul(ginv,qq(:,iq)), qx, sum(abs(qx))
      enddo
      print *,'iqindx2sym: ERROR! we can not find proper iq ###'
      stop    'iqindx2sym: ERROR! we can not find proper iq ###'
      end

      subroutine iqindx2(q, ginv,qq,nq, iqindx,qu)
C- Find index as q=qq(:,iq) with modulo of premitive vector.
Ci ginv is the inverse of plat (premitive translation vector).
      implicit none
      integer(4):: nq,i_out, iq,iqx,iqindx !,saveiq
      real(8):: q(3),qq(3,nq),ginv(3,3),qx(3),qu(3)
      integer(4),save :: isave=0  !accelation to find iqindx
      logical::debug=.false.
      if(debug) write(*,"(' iqindx2: q=',3f20.15)") q
c      if(isave>nq.or.saveiq()==0) isave=0
      if(isave>nq) isave=0
      do iqx = isave+1,isave+nq
        if(iqx > nq) then
          iq = iqx -nq
        else
          iq = iqx
        endif
        call rangedq(matmul(ginv,q-qq(:,iq)), qx)
c        print *,' iq,q-qq(:,iq)',iq,q-qq(:,iq)
        if(sum(abs(qx))< 1d-6) then
          iqindx=iq
          isave =iq
          qu=qq(:,iq)
c          print *,"iqindx2:isave iqindx=",isave,iqindx
          return
        endif
      enddo
! error exit
      write(6,"(' q  =  ',3f13.5,' ginv*q=',3f13.5)")q, matmul(ginv,q)
      write(6,"(' iq    qq    ginv*qq     qq-ginv*qq     err  ')")
      do iq = 1,nq
        call rangedq(matmul(ginv,q-qq(:,iq)), qx)
        write(6,"(i3,3f13.5,' | ',3f13.5,' | ',3f13.5,' diff= ',d13.6)")
     &       iq, qq(1:3,iq), matmul(ginv,qq(:,iq)), qx, sum(abs(qx))
      enddo
      print *,'iqindx2: ERROR! we can not find proper iq ###'
      stop    'iqindx2: ERROR! we can not find proper iq ###'
      end
c----------------------------------------------------------
      integer(4) function iqindx(q, ginv,qq,nq)
      implicit none
      integer(4):: nq,i_out, iq,iqx,iqindx0
      real(8):: q(3),qq(3,nq),ginv(3,3),qx(3)
      integer(4),save :: isave=0  !accelation to find iqindx
      if(isave>nq) isave=0
      iqindx= iqindx0(q, ginv,qq,nq,isave)
      end
      integer(4) function iqindx0(q, ginv,qq,nq,isave)
C- Find index as q=qq(:,iq) with modulo of premitive vector.
Ci ginv is the inverse of plat (premitive translation vector).
      implicit none
      integer(4):: nq,i_out, iq,iqx,isave
      real(8):: q(3),qq(3,nq),ginv(3,3),qx(3)
      do iqx = isave+1,isave+nq
        if(iqx > nq) then
          iq = iqx -nq
        else
          iq = iqx
        endif
        call rangedq(matmul(ginv,q-qq(:,iq)), qx)
        if(sum(abs(qx))< 1d-6) then
          iqindx0=iq
          isave =iq
c          print *,"isave=",isave
          return
        endif
      enddo
! error exit
      write(6,"(' q  =  ',3f13.5,' ginv*q=',3f13.5)")q, matmul(ginv,q)
      write(6,"(' iq    qq    ginv*qq     qq-ginv*qq     err  ')")
      do iq = 1,nq
        call rangedq(matmul(ginv,q-qq(:,iq)), qx)
        write(6,"(i3,3f13.5,' | ',3f13.5,' | ',3f13.5,' diff= ',d13.6)")
     &       iq, qq(1:3,iq), matmul(ginv,qq(:,iq)), qx, sum(abs(qx))
      enddo
      print *,'iqindx: ERROR! we can not find proper iq ###'
      stop    'iqindx: ERROR! we can not find proper iq ###'
      end
