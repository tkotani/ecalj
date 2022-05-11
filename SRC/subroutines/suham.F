      module m_suham
      integer,protected:: ham_ndham, ham_ndhamx,ham_nspx  !only these integers returned. 
      contains
      subroutine m_suham_init()
      use m_struc_def  
      use m_lmfinit,only:lso, ctrl_nbas,ctrl_nspec, ctrl_nl,nspc,nsp,nlmto,iprmb,nbas,nl,nspec,
     &     pwmode=>ham_pwmode,pwemax,pwemin
     &     ,alat=>lat_alat, lat_tolft,nkaph, pot_nlma, pot_nlml ,stdo,nlmto
      use m_supot,only: lat_ng, rv_a_ogv
      use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
      use m_MPItk,only:master_mpi
!! Hamiltonian dimensiton: APW part
!!       ham_ndham,ham_ndhamx,ham_nspx
Cl   ndim  :total number of lmto orbitals = nl**2 * nbas
Cl   npwmin:lower limit to number of PWs, estimated by sweeping
Cl         :over a fine mesh of qp.
Cl   npwmax:upper limit to number of PWs, estimated by sweeping
Cl         :over a fine mesh of qp.
Cl   nqdiv: loop over mesh of q-points to estimate npwmax
Cl        : nqdiv is fineness of q-mesh.
Cl  npwpad: a 'safety' padding to npwmax in case npwmax
Cl        : underestimates actual upper limit !required???
Cu     ndham = estimate for upper dimension of
Cu             hamiltonian, including possible PW part
C ----------------------------------------------------------------------
      implicit none
      integer hord,i,i1,i2,iprint, 
     . lidim,lihdim,nclasp,ndim,neul,
     .nlspcp,nspx,nttab,partok,igets,nvi,nvl, 
     .ib,is,lmxa,lmxl,isw,nbf,lmxax, j1,j2,j3,m,npw,npwmin,npwmax
      integer:: ndham=0
      double precision q(3),Gmin,Gmax,xx
      integer,parameter:: nqdiv=12
      integer ::npwpad
      double precision ckbas,cksumf,kap2(10)
      integer:: oalph , oiax , ontab , os , ng
      real(8),allocatable :: rv_a_ogvx(:)
      integer ,allocatable :: kv_iv(:)
      integer ,allocatable :: igv2_iv(:)
      integer:: obas , oo , opp ,  opti , obs , omagf,i_copy_size
      call tcn('m_suham_init')
C --- Hamiltonian offsets, orbital permutation table ---
      if (mod(pwmode,10)==2.and.master_mpi) write(stdo,"(a)")' suham: no MTO basis'
C   ... PW setup : estimate upper bound to number of G vectors
C     to set up upper bound to hamiltonian dimension
      ham_ndham = nlmto
      ndham     = nlmto
      if (pwemax>0 .and. mod(pwmode,10)>0) then
        Gmin = dsqrt(pwemin)
        Gmax = dsqrt(pwemax)
        if (mod(pwmode/10,10) .eq. 1) then
          write(stdo,*)' Estimate max size of PW basis fromcombinations of recip. lattice vectors ...'
          npwmax = -1
          npwmin = 99999
          do  j1 = 0, nqdiv
            do  j2 = 0, nqdiv
              do  j3 = 0, nqdiv
                do   m = 1, 3
                  q(m) = (qlat(m,1)/nqdiv)*j1 +(qlat(m,2)/nqdiv)*j2 +(qlat(m,3)/nqdiv)*j3
                enddo
                call pshpr(iprint()-40)
                call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,xx,xx,xx)
                call poppr
                npwmin = min(npwmin,npw)
                npwmax = max(npwmax,npw)
              enddo
            enddo
          enddo
          npwpad = max(nint((npwmax-npwmin)*0.2d0),3) !we need this for safe?
        else
          q=0d0
          call pshpr(0)
          call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,xx,xx,xx)
          call poppr
          npwmin = npw
          npwmax = npw
          npwpad = 0
        endif
        ndham = npwmax + npwpad
        if (mod(pwmode,10) .ne. 2) ndham = nlmto + npwmax + npwpad
        ham_ndham = ndham
!! Dimension maximum of Hamiltonian is ndhamx (ispx=1,nspx)
!! Note spinoffdiag=T case: we use nspc,nsp,nspx,ndhamx (a little complicated, I think).
!     nspx*nspc=nsp
        write(stdo,"(a,2f9.3,2i5)")' suham: PW basis Emin Emax npw ndham=',pwemin,pwemax,npw,ndham
C    ...  Printout
        if (iprint() .ge. 40) then
          call pshpr(0)
          call dpzero(q,3)
          call gvlst2(alat,plat,q,0,0,0,Gmin,Gmax,0,0,0,npw,xx,
     .    xx,xx,xx)
          call poppr
          write(stdo,*)' G vectors at the Gamma point:'
          allocate(igv2_iv(3*npw))
          allocate(kv_iv(3*npw))
          allocate(rv_a_ogvx(abs(3*npw))); rv_a_ogvx(:)=0.0d0
c          call pshpr(iprint())
c          if (iprint() .ge. 50) call setpr(100)
          call gvlst2 ( alat , plat , q , 0 , 0 , 0 , gmin , gmax , 0 ,
     .     8 + 2 , npw , npw , kv_iv , rv_a_ogvx , xx , igv2_iv )
c          call poppr
          if (allocated(kv_iv)) deallocate(kv_iv)
          if (allocated(igv2_iv)) deallocate(igv2_iv)
          deallocate(rv_a_ogvx)
        endif
      endif
      if(nspc==2.and.nsp==2) then 
c         spinoffdiag=.true. 
         ham_nspx= 1               
         ham_ndhamx = ndham*2
      elseif(nspc==1.and.nsp==2) then
c         spinoffdiag=.false.    !spin no-offdiagonal
         ham_nspx= 2 
         ham_ndhamx= ndham
      elseif(nspc==1.and.nsp==1) then
c         spinoffdiag=.false.    !paramagnetic
         ham_nspx= 1            !nsp/nspc
         ham_ndhamx= ndham
      else
         call rx('suham: nspc==2 but nsp=1')
      endif
      call tcx('m_suham_init')
      end subroutine m_suham_init
      end module m_suham



      
csssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine mptauof(symops,ng,plat,nbas,bas,
     &iclass,miat,tiat,invg,delta)

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
     &nbmx, nsymx, ig,igd,i,j,ibas,mi,i1,i2,i3
      double precision SYMOPS(9,ng),plat(3,3),
     &tiat(3,nbas,ng),am(3,3),b1,b2,b3,bas(3,nbas),
     &tr1,tr2,tr3,ep, dd1,dd2,dd3,t1,t2,t3
      integer::  iprintx=0
      integer ires(3, nbas, ng)
ckino delete integer(4) def.      integer(4):: ib1,ib2
      integer:: ib1,ib2
c
      real(8) ::tran(3),delta(3,ng)
      data ep/1.0d-3/
c      data ep/1.0d-7/
c
      if(iprintx>=46) write(6,*)'MPTAUOf: search miat tiat for wave function rotation'

      do 10 ig=1,ng
        do igd=1,ng
c seach for inverse  ig->igd
          if( abs( symops(1,ig)-symops(1,igd) ).le.ep.and.
     &    abs( symops(2,ig)-symops(4,igd) ).le.ep.and.
     &    abs( symops(3,ig)-symops(7,igd) ).le.ep.and.
     &    abs( symops(4,ig)-symops(2,igd) ).le.ep.and.
     &    abs( symops(5,ig)-symops(5,igd) ).le.ep.and.
     &    abs( symops(6,ig)-symops(8,igd) ).le.ep.and.
     &    abs( symops(7,ig)-symops(3,igd) ).le.ep.and.
     &    abs( symops(8,ig)-symops(6,igd) ).le.ep.and.
     &    abs( symops(9,ig)-symops(9,igd) ).le.ep  ) then
            invg(ig)=igd
            goto 16
          endif
        end do
 16     continue
c
c$$$        if(iprintx .ge.46) then
c$$$          print *,' '
c$$$          print *,' '
c$$$          print *,' **** group ops no. ig (igd)= ', ig, invg(ig)
c$$$          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
c$$$          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
c$$$          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
c$$$ 1731     format (' ',3f9.4)
c$$$        endif

        do i=1,3
          do j=1,3
            am(i,j)=symops(i+3*(j-1),ig)
          end do
        end do
c
c trial shift vector tran
        do 120 ib1=1,nbas
        do 121 ib2=1,nbas
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
              if( iclass(mi).ne.iclass(ibas) ) cycle
              do  i1=-3,3
              do  i2=-3,3
              do  i3=-3,3
                dd1 = ( i1 *plat(1,1)+i2 *plat(1,2)+i3 *plat(1,3) )
                dd2 = ( i1 *plat(2,1)+i2 *plat(2,2)+i3 *plat(2,3) )
                dd3 = ( i1 *plat(3,1)+i2 *plat(3,2)+i3 *plat(3,3) )

                t1 = b1 - (bas(1,mi)+dd1)
                t2 = b2 - (bas(2,mi)+dd2)
                t3 = b3 - (bas(3,mi)+dd3)
                if(abs(t1).le.ep.and.abs(t2).le.ep.and.
     &               abs(t3).le.ep) go to 60
             enddo
             enddo
             enddo
   40       continue
c seach failed, Not found mi and dd1. Try next (tr).
            goto 121

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
  121   continue
  120   continue
        call rx('mptauof: Can not find miat and tiat')
c
   21   continue
        delta(:,ig) = tran          ! r' = am(3,3) r +  delta  !Jun 2000

c- have gotten the translation-> check write --------------------
        if(iprintx.ge.46) then
c          write(6,4658)tr1,tr2,tr3
          write(6,4658)tran
 4658     format('  Obtained translation operation=',3d12.4)
          do 123  ibas=1,nbas
            write(6,150) ibas, miat(ibas,ig), tiat(1,ibas,ig),
     &      tiat(2,ibas,ig), tiat(3,ibas,ig),
     &      ires(1,ibas,ig),ires(2,ibas,ig),ires(3,ibas,ig)
  150       format(' ibas=',i3,' miat=',i3,' tiat=',3f11.4,' i1i2i3=',3i3)
  123     continue
        endif
c---------------------------------------------------
   10 continue
      end subroutine mptauof


      subroutine rotdlmm(symops,ng,nl ,dlmm)
c- Generate rotation matrix D^l_{m,m'} for L-representaiton, corresponding
c  to points group operations.
Ci symops(9,ng),ng; point ops.
Ci nl; num.of l =lmax+1
Co dlmm(2*nl-1,2*nl-1,0:nl-1,ng,2); D^l_{m,m'}. Indexes are for Real harmonics.
Cr   dlmmc is used as work area about 200kbyte used for  s,p,d,f -> nl=4
c-----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer:: is,i,ig,ikap,j,l,m,m1,m2,m3,md,mx,ix
      integer ng,nl
      double precision SYMOPS(9,ng),
     .am(3,3) ,fac1,fac2
      double precision dlmm( -(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)

      double complex   dlmmc(-(nl-1):(nl-1),-(nl-1):(nl-1),0:nl-1,ng)
      double precision det,igann,osq2
      double complex   msc(0:1,2,2), mcs(0:1,2,2),Img
     &,dum(2)
      parameter(Img=(0d0,1d0))
      integer:: debugmode
      real(8):: ep=1d-3 !ep was 1d-8 before feb2013
c      print *; print *,' ROTDLMM:'
      do 10 ig =1,ng
        do  i=1,3
        do  j=1,3
           am(i,j) = symops(i+3*(j-1),ig)
        enddo
        enddo
c calculate determinant(signature)
        det= am(1,1)*am(2,2)*am(3,3)
     &  -am(1,1)*am(3,2)*am(2,3)
     &  -am(2,1)*am(1,2)*am(3,3)
     &  +am(2,1)*am(3,2)*am(1,3)
     &  +am(3,1)*am(1,2)*am(2,3)
     &  -am(3,1)*am(2,2)*am(1,3)
        if(abs(abs(det)-1d0).ge.1d-10) then
          print *,' rotdlmm: det/=1 ig and det=',ig,det
          stop
        endif
c seek Euler angle   print *,' goto cbeta',ig,det
        cbeta = am(3,3)/det
c added region correction so as to go beyond domain error for functions, dsqrt and acos.
        if(abs(cbeta-1d0).le.1d-6) cbeta= 1d0
        if(abs(cbeta+1d0).le.1d-6) cbeta=-1d0
        beta = dacos(cbeta)
        sbeta= sin(beta)
c beta= 0~pi
        if(sbeta.le.1.0d-6) then
          calpha= 1d0
          salpha= 0d0
          alpha = 0d0
          cgamma= am(2,2)/det
          sgamma= am(2,1)/det
        else
          salpha =  am(2,3)/sbeta/det
          calpha =  am(1,3)/sbeta/det
          sgamma =  am(3,2)/sbeta/det
          cgamma = -am(3,1)/sbeta/det
        endif
        co2 = dcos(beta/2d0)
        so2 = dsin(beta/2d0)
c         print *,' calpha=',calpha
        if(abs(calpha-1.0d0).le.1.0d-6) calpha= 1.0d0
        if(abs(calpha+1.0d0).le.1.0d-6) calpha=-1.0d0
        if(abs(cgamma-1.0d0).le.1.0d-6) cgamma= 1.0d0
        if(abs(cgamma+1.0d0).le.1.0d-6) cgamma=-1.0d0
        alpha=dacos(calpha)
        if(salpha.lt.0d0) alpha=-alpha
        gamma=dacos(cgamma)
        if(sgamma.lt.0d0) gamma=-gamma
c         print *,'alpha beta gamma det=',alpha,beta,gamma,det
        do l =  0, nl-1
        do md= -l, l
        do m = -l, l
c  from 'Ele theo. ang. mom. by M. E. Rose 5th 1967 Wisley and Sons.  p.52 (4.13)
          fac1 = dsqrt( igann(l+m)*igann(l-m)*igann(l+md)*igann(l-md) )
          fac2 = 0d0
          do ikap=0,2*l
            if(l-md-ikap.ge.0 .and. l+m-ikap.ge.0
     &          .and.ikap+md-m.ge.0) then
              add= dble((-1)**ikap)/( igann(l-md-ikap)*igann(l+m-ikap)
     &            *igann(ikap+md-m)*igann(ikap) )
              if(2*l+m-md-2*ikap.ne.0) add=add*co2**(2*l+m-md-2*ikap)
              if(md-m+2*ikap.ne.0)     add=add*(-so2)**(md-m+2*ikap)
              fac2 = fac2+add
           endif
          enddo
c l-th rep. is odd or even according to (det)**l
          dlmmc(md,m,l,ig) = fac1*fac2*det**l*
     &        cdexp( -Img*(alpha*md+gamma*m) )
        enddo
        enddo
        enddo
        am(1,1)= cos(beta)*cos(alpha)*cos(gamma)-sin(alpha)*sin(gamma)
        am(1,2)=-cos(beta)*cos(alpha)*sin(gamma)-sin(alpha)*cos(gamma)
        am(1,3)= sin(beta)*cos(alpha)
        am(2,1)= cos(beta)*sin(alpha)*cos(gamma)+cos(alpha)*sin(gamma)
        am(2,2)=-cos(beta)*sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)
        am(2,3)= sin(beta)*sin(alpha)
        am(3,1)=-sin(beta)*cos(gamma)
        am(3,2)= sin(beta)*sin(gamma)
        am(3,3)= cos(beta)

        if(abs(am(1,1)*det-symops(1,ig))>ep.or.
     &  abs(am(2,1)*det-symops(2,ig))>ep.or.
     &  abs(am(3,1)*det-symops(3,ig))>ep.or.
     &  abs(am(1,2)*det-symops(4,ig))>ep.or.
     &  abs(am(2,2)*det-symops(5,ig))>ep.or.
     &  abs(am(3,2)*det-symops(6,ig))>ep.or.
     &  abs(am(1,3)*det-symops(7,ig))>ep.or.
     &  abs(am(2,3)*det-symops(8,ig))>ep.or.
     &  abs(am(3,3)*det-symops(9,ig))>ep) then
          print *,' rotdlmm: not agree. symgrp and one by eular angle'
          stop
        endif
cccccccccccccccccccccccc
c        if(iprint().ge.140) then
        if(debugmode()>9) then
          print *;print *;print *,' **** group ops no. ig=', ig
          write(6,1731)symops(1,ig),symops(4,ig),symops(7,ig)
          write(6,1731)symops(2,ig),symops(5,ig),symops(8,ig)
          write(6,1731)symops(3,ig),symops(6,ig),symops(9,ig)
          print *,' by Eular angle '
          write(6,1731)am(1,1)*det,am(1,2)*det,am(1,3)*det
          write(6,1731)am(2,1)*det,am(2,2)*det,am(2,3)*det
          write(6,1731)am(3,1)*det,am(3,2)*det,am(3,3)*det
        endif
 1731   format (' ',3f9.4)
cccccccccccccccccccccccc
   10 continue
c conversion to cubic rep. Belows are from csconvs
c  msc mcs conversion matrix generation 2->m 1->-m for m>0
      osq2 = 1d0/sqrt(2d0)
      do m = 0,1
        Msc(m,1,1)= osq2*(-1)**m
        Msc(m,1,2)=-osq2*Img*(-1)**m
        Msc(m,2,1)= osq2
        Msc(m,2,2)= osq2*Img

        Mcs(m,1,1)= osq2*(-1)**m
        Mcs(m,1,2)= osq2
        Mcs(m,2,1)= osq2*Img*(-1)**m
        Mcs(m,2,2)=-osq2*Img
      enddo
c
      if(debugmode()>1) print * ,' goto do 123'
      do 123 is=1,ng
c$$$        if(.false.) then
c$$$c        if(iprint().ge.150) then
c$$$          print *; print *,' **** group ops no. ig=', is
c$$$          write(6,1731) symops(1,is),symops(4,is),symops(7,is)
c$$$          write(6,1731) symops(2,is),symops(5,is),symops(8,is)
c$$$          write(6,1731) symops(3,is),symops(6,is),symops(9,is)
c$$$        endif
c convert to cubic rep.
      do 23   l =0,nl-1
        do  m2=-l,l
        do  m1= 1,l
          dum(1)= dlmmc(m2, m1,l,is)
          dum(2)= dlmmc(m2,-m1,l,is)
          mx    = mod(m1,2)
          dlmmc(m2,  m1,l,is)=
     &        dum(1)*msc(mx,1,1)
     &        +dum(2)*msc(mx,2,1)
          dlmmc(m2, -m1,l,is)=
     &        dum(1)*msc(mx,1,2)
     &         +dum(2)*msc(mx,2,2)
        enddo
        enddo
        do m2=  1,l
        do m1= -l,l
          dum(1)=dlmmc( m2, m1,l,is)
          dum(2)=dlmmc(-m2, m1,l,is)
          mx=mod(m2,2)
          dlmmc( m2, m1,l,is)=
     &        mcs(mx,1,1)*dum(1)
     &        +mcs(mx,1,2)*dum(2)
          dlmmc(-m2, m1,l,is)=
     &        mcs(mx,2,1)*dum(1)
     &         +mcs(mx,2,2)*dum(2)
        enddo
        enddo
        do m2=-l,l
        do m1=-l,l
          dlmm(m2,m1,l,is)=dreal( dlmmc(m2,m1,l,is) )
          if( abs(dimag(dlmmc(m2,m1,l,is))).ge.1.0d-12 ) stop
     &        ' rotdlmm: abs(dimag(dlmmc(m2,m1,l,is))).ge.1.0d-12'
        enddo
        enddo
ccccccccccccccccccccc
        if(.false.) then
c        if(.true.) then
c        if(iprint().ge.41) then
          print *; print *,'  points ops  ig, l=', is,l,' cubic   '
          do m2=-l,l
            write(6,"(28f10.5)")( dreal(dlmmc (m2, m1,l,is) ), m1=-l,l)
c    &    , ( dimag(dlmmc (m2, m1,l,is) ), m1=-l,l),( dlmm(m2, m1,l,is), m1=-l,l)
          enddo
        endif
cccccccccccccccccccccc
   23 continue
 123  continue
      if(debugmode()>1) print *,' end of rotdlmm'
      end subroutine rotdlmm

c--------------------------------------------
      double precision function igann(i)
      integer:: i,ix
      igann  = 1d0
      do ix =1,i
        igann=igann*dble(ix)
      enddo
      end function igann

