      subroutine amagnc(nbas,nl,ipc,rhos,nrhos,qnu,eula,neul,mode,
     .  amag,aamom)
C- Printout magnetic moments in unit cell
C     implicit none
C Passed variables
      integer nbas,neul,nl,nrhos,ipc(*),mode
      double precision eula(nbas,neul,3),rhos(2,0:2,nrhos,2,2,1),
     .  qnu(3,nl,2,1),amag(3),aamom(nbas)
C Local variables

      call rx('amagnc require noncollinear package')

      end
      subroutine amagn2(mode,nl,nlo,nclass,nbas,ipc,eula,neul,bxc,qnu,
     .  nrhos,rhos)

C     implicit none
C ... Passed parameters
      integer mode,nl,nlo,nclass,nbas,neul,nrhos,ipc(nbas)
      double precision eula(nbas,neul,3),rhos(2,0:2,nrhos,2,2,nclass)
      double precision qnu(3,nl,2,1),bxc(3,*)

      call rx('amagnc require noncollinear package')
      end

