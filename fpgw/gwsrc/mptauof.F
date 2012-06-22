      subroutine eibzsym(nqibz,symgg,ngrp,qibze,iqxini,iqxend,qbz,nqbz,timereversal,ginv,
     o qeibz,neibz,nwgt,aik,aiktimereversal,ngrpt,igx)
!! === Obtain info for eibz symmetrization. See PRB81,125102-9 ===
!! For GaAs 4x4x4 (with timereversal), we have IBZxBZ=10x640 is reduced to be IBZxEBZ=286.
      implicit none
      integer:: nqibz,ngrp,iqxini,iqxend,iq,ig,igxx,neibzx,ieibz,ibz,nwgtsum,itimer,ntimer
      integer:: ngrpt(iqxini:iqxend),neibz(iqxini:iqxend),nwgt(nqbz,iqxini:iqxend),ik,nqbz,
     &  aiktimereversal(2*ngrp,iqxini:iqxend)
      logical:: timereversal
      real(8):: symgg(3,3,ngrp),qibze(3,iqxini:iqxend),qbz(3,nqbz),qeibz(3,nqbz,iqxini:iqxend),
     &  aik(3,3,2*ngrp,iqxini:iqxend),q(3),qx(3),qlat(3,3),qdiff(3),ginv(3,3),qxx(3),ddd,
     &  qout(3,nqbz,iqxini:iqxend,nqbz),sss(3,3),sumcheck1
      integer:: igx(nqbz,iqxini:iqxend,nqbz) ! ieibz,iq,iw.

!!    check if  symgg(:,:,1) 
      sss(1,:)=(/1d0,0d0,0d0/)
      sss(2,:)=(/0d0,1d0,0d0/)
      sss(3,:)=(/0d0,0d0,1d0/)
      sumcheck1= sum(abs(symgg(:,:,1)-sss))
      if(sumcheck1>1d-8) stop 'eibzsym: symgg(:,:,1) is not E'

      write(6,"('eibzsym: TimeReversal ngrp= ',l,i3)") timereversal,ngrp
      nwgt=0
      do iq=iqxini,iqxend
        q=qibze(:,iq) !q means k in eq.50 of PRB81,125102-9
        igxx=0
        ntimer=1
        if(timereversal) ntimer=-1
        do itimer=1,ntimer,-2
          do ig=1,ngrp
            if(sum(abs(dble(itimer)*matmul(symgg(:,:,ig),q)-q))<1d-8) then
              igxx = igxx+1
              aik(:,:,igxx,iq)=symgg(:,:,ig)
              aiktimereversal(igxx,iq)=itimer
            endif
          enddo
        enddo
        ngrpt(iq)=igxx

        neibzx=0
        do ibz=1,nqbz
          qx = qbz(:,ibz) !q means a candidate of q in EIBZ
          do ieibz=1,neibzx
            do igxx=1,ngrpt(iq)
              qdiff= - dble(aiktimereversal(igxx,iq))*matmul(aik(:,:,igxx,iq),qx) - qeibz(:,ieibz,iq)
              call rangedq(matmul(ginv,qdiff), qxx)
              if(sum(abs(qxx))<1d-8) then
                nwgt(ieibz,iq)=nwgt(ieibz,iq)+1
                qout(:,ieibz,iq, nwgt(ieibz,iq))=qx
                igx(ieibz,iq, nwgt(ieibz,iq))=igxx !symmetrizer
                goto 1012
              endif
            enddo
          enddo
          neibzx = neibzx+1
          nwgt(neibzx,iq)=1
          qeibz(:,neibzx,iq)=qx
          igx(neibzx,iq,nwgt(ieibz,iq))=1 ! At the beginig of this routine, sumcheck1 shows 1 means E.
 1012     continue
        enddo
        neibz(iq)=neibzx
        write(6,"('eibz: iq q ngrpt neibz nqbz= ',i3,3f10.4,3i3)") iq,q,ngrpt(iq),neibz(iq),nqbz
         do ik=1,neibz(iq)
           write(6,"('      ',i3,3f10.4,i4,2x,25i3)") ik, qeibz(:,ik,iq),nwgt(ik,iq)
     &       ,igx(ieibz,iq, 1:nwgt(ik,iq))
         enddo
      enddo
      nwgtsum = 0
      do iq = iqxini,iqxend
        if(iq>nqibz) cycle
        nwgtsum = nwgtsum + sum(nwgt(1:neibz(iq),iq))
        print *, iq, neibz(iq),sum(nwgt(1:neibz(iq),iq))
      enddo
      print *,' sum(neibz(iq)) sum(ngwt)=sum(nqbz)= ',
     &   sum(neibz(iqxini:iqxend-1)), nwgtsum, nqbz*(iqxend-iqxini+1-1)
      end

      subroutine mptauof(symops,ng,plat,nbas,bas,
     &             iclass,miat,tiat,invg,delta)
C- Mapping of each atom by points group operations------------------c
C Modded by okuda 1994 March.
C Simplified by kotani 1994 7/31.
Ci  Input
Ci     symops(1,ng),ng,plat,nbas,bas(3,nbas)
Ci     iclass(nbas); denote class for each atom
Co  Output
Co    miat(ibas  ,ig); ibas-th atom is mapped to miat-th atom, by the ig-th
Co    points group operation.  Origin is (0,0,0).
Co    tiat(k,ibas,ig);
Co    delta : shifting vector for non-symmorphic group.
Co            r' = matmul (am, r) + delta
Cr  Remarks
Cr
Cr (1) The ibas-th atom (position at bas(k,ibas) ) is mapped to
Cr
Cr    bas( k,miat(ibas,ig) )+ tiat(k,ibas,ig), k=1~3.
Cr
Cr (2) tiat= unit translation
Cr
C--------------------------------------------------------------------
      implicit none
      integer ng,nbas, miat(nbas,ng),iclass(nbas),invg(ng),
     &        nbmx, nsymx, ig,igd,i,j,ibas,mi,i1,i2,i3
      double precision SYMOPS(9,ng),plat(3,3),
     &                 tiat(3,nbas,ng),am(3,3),b1,b2,b3,bas(3,nbas),
     &                 tr1,tr2,tr3,ep, dd1,dd2,dd3,t1,t2,t3
      integer  iprint
      external iprint
      integer ires(3, nbas, ng)
      integer(4):: ib1,ib2
c
      real(8) ::tran(3),delta(3,ng)
      data ep/1.0d-3/
c      data ep/1.0d-7/
c
      write(6,*)' MPTAUO2: search miat tiat for wave function rotation'

      do 10 ig=1,ng
        do igd=1,ng
c seach for inverse  ig->igd
          if( abs( symops(1,ig)-symops(1,igd) ).le.ep.and.
     &        abs( symops(2,ig)-symops(4,igd) ).le.ep.and.
     &        abs( symops(3,ig)-symops(7,igd) ).le.ep.and.
     &        abs( symops(4,ig)-symops(2,igd) ).le.ep.and.
     &        abs( symops(5,ig)-symops(5,igd) ).le.ep.and.
     &        abs( symops(6,ig)-symops(8,igd) ).le.ep.and.
     &        abs( symops(7,ig)-symops(3,igd) ).le.ep.and.
     &        abs( symops(8,ig)-symops(6,igd) ).le.ep.and.
     &        abs( symops(9,ig)-symops(9,igd) ).le.ep  ) then
            invg(ig)=igd
            goto 16
          endif
        end do
 16     continue
c
        if(iprint() .ge.40) then
          print *,' '
          print *,' '
          print *,' **** group ops no. ig (igd)= ', ig, invg(ig)
          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
 1731     format (' ',3f9.4)
        endif

        do i=1,3
          do j=1,3
            am(i,j)=symops(i+3*(j-1),ig)
          end do
        end do
c
c trial shift vector tran
        do 20 ib1=1,nbas
        do 20 ib2=1,nbas
          tran =  bas(:,ib2)  - matmul(am,bas(:,ib1))
c
          do 30 ibas=1,nbas
            b1=am(1,1)*bas(1,ibas)
     &        +am(1,2)*bas(2,ibas)+am(1,3)*bas(3,ibas)
     &        +tran(1)
c     .        +( tr1*plat(1,1)+tr2*plat(1,2)+tr3*plat(1,3) )
            b2=am(2,1)*bas(1,ibas)
     &        +am(2,2)*bas(2,ibas)+am(2,3)*bas(3,ibas)
     &        +tran(2)
c     .        +( tr1*plat(2,1)+tr2*plat(2,2)+tr3*plat(2,3) )
            b3=am(3,1)*bas(1,ibas)
     &        +am(3,2)*bas(2,ibas)+am(3,3)*bas(3,ibas)
     &        +tran(3)
c     .        +( tr1*plat(3,1)+tr2*plat(3,2)+tr3*plat(3,3) )
c
            do 40 mi=1,nbas
              if( iclass(mi).ne.iclass(ibas) ) go to 40

              do 50 i1=-3,3
              do 50 i2=-3,3
              do 50 i3=-3,3
                dd1 = ( i1 *plat(1,1)+i2 *plat(1,2)+i3 *plat(1,3) )
                dd2 = ( i1 *plat(2,1)+i2 *plat(2,2)+i3 *plat(2,3) )
                dd3 = ( i1 *plat(3,1)+i2 *plat(3,2)+i3 *plat(3,3) )

                t1 = b1 - (bas(1,mi)+dd1)
                t2 = b2 - (bas(2,mi)+dd2)
                t3 = b3 - (bas(3,mi)+dd3)
                if(abs(t1).le.ep.and.abs(t2).le.ep.and.
     &             abs(t3).le.ep) go to 60
   50         continue
   40       continue
c seach failed, Not found mi and dd1. Try next (tr).
            goto 20

   60       continue
            miat(ibas,ig)  = mi
            tiat(1,ibas,ig)= dd1
            tiat(2,ibas,ig)= dd2
            tiat(3,ibas,ig)= dd3
            ires(1,ibas,ig)= i1
            ires(2,ibas,ig)= i2
            ires(3,ibas,ig)= i3
c
   30     continue
c When the do-30 loop has been completed, we get out of do-20 loop
          goto 21
   20   continue
        stop 'mptauo2: Can not find miat and tiat'
c
   21   continue
        delta(:,ig) = tran          ! r' = am(3,3) r +  delta  !Jun 2000

c- have gotten the translation-> check write --------------------
        if(iprint().ge.40) then
c          write(6,4658)tr1,tr2,tr3
          write(6,4658)tran
 4658     format('  Obtained translation operation=',3d12.4)
          do 123  ibas=1,nbas
            write(6,150) ibas, miat(ibas,ig), tiat(1,ibas,ig),
     &    tiat(2,ibas,ig), tiat(3,ibas,ig),
     &    ires(1,ibas,ig),ires(2,ibas,ig),ires(3,ibas,ig)
  150       format(' ibas=',i3,' miat=',i3,' tiat=',3f11.4,' i1i2i3=',3i3)
  123     continue
        endif
c---------------------------------------------------
   10 continue
      end
