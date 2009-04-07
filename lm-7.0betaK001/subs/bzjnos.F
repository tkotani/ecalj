      subroutine bzjnos(n1,n2,n3,ep,nq,nbmx,nfilm,nempm,nsp,nfilo,nfiup,
     .   nemlo,nemup,emin,emax,esciss,jdos,optmt,dos,nr,ef,ntet,idtet)
C- BZ integration of joint NOS by linear method
C ----------------------------------------------------------------------
Ci Inputs:
Ci  ep, energy bands;
Ci  ntet, No. of different tetrahedra
Ci  idtet(1-4,i), Identifies the i'th tetrahedron in terms of the four
Ci  irreducible k-points:
Ci  idtet(0,i), no. of tetrahedra of the i'th kind
Ci  nq, no. of irreducible k-points; nsp, no spins
Ci  nbmx, dimensions ep
Ci  emin, emax, nr: energy window and nr of points (see remarks)
Ci  ef, Fermi energy.
Co Outputs:
Co  dos, Integrated Joint Density of States (idos)
Cr Remarks
Cr   Adapted from bzints to make joint density of states.  
Cr   All energy differences between states below ef and states
Cr   above ef+emin are summed, and integrated over the BZ
Cr   Treatment near the critical points (ef and ef+emin) handled crudely
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer n1,n2,n3,nq,nbmx,nsp,idtet(0:4,*),nr,ntet,nfilm,nempm,
     .  nfilo,nfiup,nemlo,nemup
      double precision ep(nbmx,nsp,nq),dos(nr,3,nsp),
     .  emin,emax,ef,esciss,optmt(3,nfilm,nempm,nsp,*)
      logical jdos
C Local parameters
      integer lgunit,ib,jb,iprint,iq1,iq2,iq3,iq4,isp,itet,npol,k,ibf,
     .  ibm,jblo
      double precision ec(4),ecj(4),ed(4),ebot,etop,
     .  volwgt,eboti,ebotj,etopi,etopj,wt,wt2,wt3,wtm

      if ( iprint() .ge. 10 ) call awrit4(
     .  ' BZJNOS:  ef=%1;6d  emin=%1;6d;  emax=%1;6d;  %i bins',
     .  ' ',80,lgunit(1),ef,emin,emax,nr-1)
      if (emin .lt. 0) print *,' BZJNOS: (warning) emin<0 for joint DOS'

      npol = 3
      if (jdos) npol = 1
      call dpzero(dos,3*nsp*nr)
      volwgt = dble(3-nsp)/(n1*n2*n3*6)
      do  40  isp = 1, nsp
C --- Loop over tetrahedra ---
        do  20  itet = 1, ntet
         iq1 = idtet(1,itet)
         iq2 = idtet(2,itet)
         iq3 = idtet(3,itet)
         iq4 = idtet(4,itet)
         ibf = 0
         do  20  ib = nfilo, nfiup
           ibf = ibf+1
C --- Set up energies at 4 corners of tetrahedron ---
           ec(1) = ep(ib,isp,iq1)
           ec(2) = ep(ib,isp,iq2)
           ec(3) = ep(ib,isp,iq3)
           ec(4) = ep(ib,isp,iq4)
           etopi = dmax1(ec(1),ec(2),ec(3),ec(4))
           eboti = dmin1(ec(1),ec(2),ec(3),ec(4))
           if (eboti .gt. ef) goto 20
C ...      wt2 cludge for handling near ef
           if (dabs(etopi-eboti) .gt. 1d-8) then
             wt2 = dmin1(1d0,(ef-eboti)/(etopi-eboti))
           else
C             print *,  ' ***WARNING***  etopi=eboti, setting wt2=1'
C             print 100,' isp,itet,ib,etopi,ef=',isp,itet,ib,etopi,ef
C  100        format(a,i3,i7,i5,2f12.6)
             wt2 = 1d0
           endif
           ibm = 0
           jblo = max0(ib+1,nemlo)
           do  30  jb = jblo, nemup
             ibm = ibm+1
C ...        Set up energies at 4 corners of tetrahedron for jb
             ecj(1) = ep(jb,isp,iq1) + esciss
             ecj(2) = ep(jb,isp,iq2) + esciss
             ecj(3) = ep(jb,isp,iq3) + esciss
             ecj(4) = ep(jb,isp,iq4) + esciss
             etopj = dmax1(ecj(1),ecj(2),ecj(3),ecj(4))
             ebotj = dmin1(ecj(1),ecj(2),ecj(3),ecj(4))
             if ( etopj .lt. ef+emin) goto 30
C             print 335, ib,ec
C             print 335, jb,ecj
C  335        format(i4,4f12.6)
             ed(1) = ep(jb,isp,iq1)-ep(ib,isp,iq1)
             ed(2) = ep(jb,isp,iq2)-ep(ib,isp,iq2)
             ed(3) = ep(jb,isp,iq3)-ep(ib,isp,iq3)
             ed(4) = ep(jb,isp,iq4)-ep(ib,isp,iq4)
             etop = dmax1(ed(1),ed(2),ed(3),ed(4))
             ebot = dmin1(ed(1),ed(2),ed(3),ed(4))
             if (ebot .gt. emax) goto 30
C ...        wt3 cludge for handling near ef+emin
             if (dabs(etopj-ebotj) .gt. 1d-8) then
               wt3 = dmin1(1d0,(etopj-ef-emin)/(etopj-ebotj))
             else
C               print *,  ' ***WARNING***  etopj=ebotj, setting wt3=1'
C               print 100,' isp,itet,jb,etopj,ef+emin=',
C     .           isp,itet,jb,etopj,ef+emin
               wt3 = 1d0
             endif
C             print 336, 'emax,ebot,etop=      ',emax,ebot,etop
C             print 336, 'eboti,etopi, ef=     ',eboti,etopi,ef
C             print 336, 'ebotj,etopj, ef+emin=',ebotj,etopj,ef+emin
C  336        format(a,3f12.6)
C             if ( ebot .lt. emax  .and. (wt2+wt3 .ne. 2)) then
C               print 336, 'wt,wt2,wt3=',wt,wt2,wt3
C             endif
             if (wt2 .gt. 1 .or. wt3 .gt. 1) stop 'bug in bzjnos'
             do  10  k = 1, npol
               wt = volwgt*idtet(0,itet)
               if (.not. jdos) then
                 wtm = optmt(k,ibf,ibm,isp,iq1)+optmt(k,ibf,ibm,isp,iq2)
     .               + optmt(k,ibf,ibm,isp,iq3)+optmt(k,ibf,ibm,isp,iq4)
                 wt = wt*wtm / 4d0
               endif
               call slinz(wt*wt2*wt3,ed,emin,emax,dos(1,k,isp),nr)
   10        continue
   30      continue
   20    continue
   40  continue

      end

