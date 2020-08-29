      subroutine gennlat(pos,nbas,plat,nk1,nk2,nk3,npairmx,ok,npair,  nlat,nqwgt)
!! find
!!i   pos,nbas,plot,nk,npairmx,range.
!!o   npair(rewritten), nlat,qwgt.
!!    qwgt is the weight for integration. Total number of qwgt is nk =nk1*nk2*nk3.
c
      implicit none
      integer:: nbas,nk1,nk2,nk3,nparimx,npair(nbas,nbas),npairmx
      integer:: ib1,ib2,nnn,ik,nmax(3),ix,i,j,ik1,ik2,ik3,nout,ni
      integer:: nlat(3,npairmx,nbas,nbas),nd(3)
      integer,parameter:: noutmx=100
      integer:: nqwgt(npairmx,nbas,nbas),nlatout(3,noutmx),debug=0
      real(8):: pos(3,nbas),plat(3,3), pi,qlat(3,3),q(3),pin(3),rmax2 !,nrmax(3)
      real(8):: eps=1d-8
      real(8):: posp(3),rxlat(3,3),rxprod(3,3),dummy,rrr,rmax
      logical :: ok
      real(8):: rlatp(3,3),xmx2(3)
!!
      call tcn('gennlat')
      print *,'gennlat:'
      ok=.true.
      pi=4d0*datan(1d0)
      call dinv33(plat,1,qlat,dummy)
! periodicity for R mesh for nk1 nk2 nk3 division.
      rxlat(:,1)=plat(:,1)*nk1
      rxlat(:,2)=plat(:,2)*nk2
      rxlat(:,3)=plat(:,3)*nk3
      call shortn3_initialize(rxlat, rlatp,xmx2)
c      do i=1,3
c      do j=1,3
c        rxprod(i,j) = sum(rxlat(:,i)*rxlat(:,j))
c      enddo
c      enddo
c      print *, 'goto ellips',rxprod
c      call ellipsoidxmax(rxprod,xmx2)
c      print *, 'end of ellips xmx2=',xmx2
      !Get Maxmum x_i**2 for ellipsoid for 1d0 = sum_{i,j} x_i rxprod(i,j) x_j
      ! Maximum value for x_i for ellipsoid
      ! Ellipsoid is given as 1d0 = sum x_i nn(i,j) x_j.
!
      npair=0
      do ib1=1,nbas
        do ib2=1,nbas
          if(debug>0) print *
          if(debug>0) print *,' ---- ib1 ib2=',ib1,ib2
          posp = matmul( pos(:,ib1)-pos(:,ib2), qlat(:,:)) !posp is on plat-coodinate
          nnn=0
          do ik1=0,nk1-1
            do ik2=0,nk2-1
              do ik3=0,nk3-1
                !print *,' ik1,ik2,ik3=',ik1,ik2,ik3
                nd=[ik1,ik2,ik3]
                if( ik1 >nk1/2) nd= nd-[nk1, 0,  0]
                if( ik2 >nk2/2) nd= nd-[0, nk2,  0]
                if( ik3 >nk3/2) nd= nd-[0, 0,  nk3]
                ! above procedures are not necessary, but these may accelarate shortn3 a little.
                pin(1)= (posp(1)+nd(1))/nk1
                pin(2)= (posp(2)+nd(2))/nk2
                pin(3)= (posp(3)+nd(3))/nk3
                call shortn3(rlatp,xmx2,pin,noutmx, nout,nlatout)

c          rmax2 = sum(pin*matmul(rxprod,pin)) + eps  ! eps is to make degeneracy safe.
c          nrmax(:) =  sqrt(rmax2*xmx2(:))+abs(pin(:)) ! range of ix used in shortn2
c          nmax =  nrmax
cccccccccccccc
c          print *,'pin=',pin
c          print *,' rmax2xmx=',rmax2,nrmax,nmax
c          nmax=4
cccccccccccccc
c          call shortn2(pin,nmax,rxprod,noutmx, nout,nlatout)
                ! find smallest nlat for mod(posp+dx, {rxlat(:,1),rxlat(:,2),rxlat(:,3)})
c          write(6,"(3i4,' nmax=',3i4,' nout=',i3)")nd(1:3),nmax(1:3),nout
                if(nnn+nout>npairmx) then
                  ok=.false.
                  return
                endif
                do ix=1,nout
                  nlat(:,nnn+ix,ib1,ib2)=
     &            nd + [nk1*nlatout(1,ix), nk2*nlatout(2,ix), nk3*nlatout(3,ix)]
c                  qwgt(nnn+ix,ib1,ib2)= 1d0/nout
                  nqwgt(nnn+ix,ib1,ib2)= nout
                enddo
                nnn = nnn+nout
              enddo
            enddo
          enddo
          npair(ib1,ib2)= nnn

c check write
          if(debug>0) then
            write(6,"(a,2i4,2x,i6)") ' ib1 ib2 npair=',ib1,ib2,npair(ib1,ib2)
            do ni = 1,npair(ib1,ib2)
              posp =  pos(:,ib1)-pos(:,ib2) + matmul(plat,nlat(:,ni,ib1,ib2))
              rrr = sqrt(sum(posp**2))
              write(6,"(i6,3x,3i3,f8.3,i5)") ni,nlat(1:3,ni,ib1,ib2),rrr,nqwgt(ni,ib1,ib2)
            enddo
          endif  
        enddo
      enddo
c      qwgt=qwgt/(nk1*nk2*nk3)
      do ib1=1,nbas
      do ib2=1,nbas
        if(abs(sum(1d0/nqwgt(1:npair(ib1,ib2),ib1,ib2))-nk1*nk2*nk3)>1d-8) call rx('bug:nqwgt sum is not unity')
      enddo
      enddo
      call tcx('gennlat')
      end


      subroutine pairtakao(job,pos,nbas,plat,nk,npair, npairmx,range, nlat,qwgt)
c Obtain pair table of atomic sites.
c    number of pairs are larger than nk. range(ib1,ib2) is calculated.
c job=1
ci   pos,nbas,plat,nk
co   npairmx,range
c job=2
ci   pos,nbas,plot,nk,npairmx,range.
co   npair(rewritten), nlat,qwgt.
cr   qwgt is the weight for integration. Total number of qwgt is nk =nk1*nk2*nk3.
c
      implicit none
      integer::nbas,job,nk,npair(nbas,nbas),npairmx,ibas,ib1,ib2,ni,isig,nnn,iwend,nend,iwinit
      real(8):: pos(3,nbas),plat(3,3), xx,pi,qlat(3,3),q(3),rrr(3),rx
      integer :: nlat(3,npairmx,nbas,nbas)
      real(8)::     qwgt(npairmx,nbas,nbas)
      real(8),allocatable:: rr(:,:,:),rrx(:)
      integer,allocatable:: iord(:),nlatx(:,:)
      real(8):: range(nbas,nbas),rangex
      logical:: exiton
      pi=4d0*datan(1d0)
      call dinv33(plat,1,qlat,xx)
      print *,'job=',job,nk,npairmx
      if(job==1) then !determine range,and npairmx
        npair=0
        do ib1=1,nbas
          do ib2=1,nbas
            exiton=.false.
            rangex=0d0
            do
              rangex= rangex+.2d0 !We will improve this to give some better guess.
              call getgv2( 2d0*pi,qlat,plat,pos(:,ib1)-pos(:,ib2), rangex,job,
     o        npair(ib1,ib2), nlat(:,:,ib1,ib2) )
              if(npair(ib1,ib2)>=nk) then
                rangex= rangex+.1d0 !We will improve this to give some better guess.
                call getgv2( 2d0*pi,qlat,plat,pos(:,ib1)-pos(:,ib2), rangex,job,
     o          npair(ib1,ib2), nlat(:,:,ib1,ib2) )
                range(ib1,ib2)=rangex ! rangex + 0.01d0 is to make degeneracy safe.
                print *,' ib1 ib2 rangex nk npair',ib1,ib2,rangex,nk,npair(ib1,ib2)
                exit
              endif
            enddo
          enddo
        enddo
        npairmx = maxval(npair)
        npair=-99999 !not the output
        return
      elseif(job==2) then
        allocate(rr(npairmx,nbas,nbas))
        do ib1=1,nbas
          do ib2=1,nbas
            call getgv2( 2d0*pi,qlat,plat,pos(:,ib1)-pos(:,ib2), range(ib1,ib2),job,
     o      npair(ib1,ib2), nlat(:,:,ib1,ib2) )
            nnn = npair(ib1,ib2)
            do ni=1,nnn
              rr(ni,ib1,ib2)= sqrt(sum((pos(:,ib1)-pos(:,ib2)+ matmul(plat, nlat(:,ni,ib1,ib2)))**2))
            enddo
c sort
            allocate(iord(nnn),nlatx(3,nnn),rrx(nnn))
            call sortea(rr(:,ib1,ib2),iord,nnn,isig)
            nlatx = nlat(:,1:nnn,ib1,ib2)
            rrx   = rr(1:nnn,ib1,ib2)
            do ni = 1,nnn
              nlat(:,ni,ib1,ib2) = nlatx(:,iord(ni))
              rr(ni,ib1,ib2)     = rrx(iord(ni))
            enddo
            deallocate(iord,nlatx,rrx)
c weight
            rx = rr(nk,ib1,ib2)
            iwinit=-99
            qwgt(:,ib1,ib2)=0d0
            do ni=1,nnn
              if(abs(rr(ni,ib1,ib2)-rx)< 1d-6.and.iwinit==-99) iwinit=ni
              if(abs(rr(ni,ib1,ib2)-rx)< 1d-6) iwend=ni
            enddo
            qwgt(1:iwinit-1,ib1,ib2)=1d0
            qwgt(iwinit:iwend,ib1,ib2)=dble(nk-iwinit+1)/(iwend-iwinit+1)
            npair(ib1,ib2)=iwend
            qwgt(:,ib1,ib2)=qwgt(:,ib1,ib2)/nk
            if(abs(sum(qwgt(:,ib1,ib2))-1d0)>1d-8) then
              stop 'paritakao: abs(sum(qwgt(:,ib1,ib2))-nk)/=0'
            endif
c check write
            write(6,"(a,2i4,2x,i6)") ' ib1 ib2 npair=',ib1,ib2,npair(ib1,ib2)
            do ni = 1,npair(ib1,ib2)
              write(6,"(i6,3x,3i3,2f8.3)") ni,nlat(1:3,ni,ib1,ib2),rr(ni,ib1,ib2),qwgt(ni,ib1,ib2)
            enddo
c          do ni = npair(ib1,ib2)+1,nnn
c            write(6,"(i6,3x,3i3,2f8.3)") ni,nlat(1:3,ni,ib1,ib2),rr(ni,ib1,ib2)
c          enddo
          enddo
        enddo
        deallocate(rr)
      endif
      end
c--- takao's getgv2. This is in fpgw/gwsrc/.
      subroutine getgv2(alat,plat,qlat,q, QpGcut,job,
     o     ng, ngvec)
C-  Set up a list of recip vectors within cutoff |Q+G| < QpGcut a.u.
c takao
c job==1 -> return ng (number of G ) only
c job==2 -> return ng and ngvec
c True G is given by
c    G(1:3,1:ng) = 2*pi/alat * matmul(qlat * ngvec(1:3,1:ng))
c
      implicit none
ckino delete integer(4) def.      integer(4) :: ng
      integer :: ng
      real(8):: s_lat(1),q(3),plat(3,3),qlat(3,3),qpg(3),enor(3)
     &,pi,alat,tpiba,QpGmax,QpGmax2,QpGcut,Qenor
ckino delete integer(4) def.      integer(4) ::
      integer ::  
     &n1max,n1min,n2max,n2min,n3max,n3min,
     &i1,i2,i3,ig,job,imx
ckino delete integer(4) def.     integer(4),optional ::  ngvec(3,ng)
ckino and delete optional def.
      integer ::  ngvec(3,ng)
      pi=4d0*datan(1d0)
c      call u_lat_vecs(s_lat,alat,plat,qlat)
      tpiba=2*pi/alat
      QpGmax   = QpGcut/tpiba  ! QpGcut in a.u.= tpiba*Qcut
      QpGmax2  = QpGmax**2
c search range of n1 n2 n3
      call eprod(qlat(1:3,2),qlat(1:3,3),enor)
      Qenor = sum(qlat(1:3,1)*enor)
      n1max =  QpGmax/abs(Qenor) - sum(q*enor)/Qenor +1
      n1min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1

      call eprod(qlat(1:3,3),qlat(1:3,1),enor)
      Qenor = sum(qlat(1:3,2)*enor)
      n2max =  QpGmax/abs(Qenor) - sum(q*enor)/Qenor +1
      n2min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1

      call eprod(qlat(1:3,1),qlat(1:3,2),enor)
      Qenor = sum(qlat(1:3,3)*enor)
      n3max =  QpGmax/abs(Qenor) - sum(q*enor)/Qenor +1
      n3min = -QpGmax/abs(Qenor) - sum(q*enor)/Qenor -1

c get ngvec within the limit.
      ig=0
      imx=-9999
      do i1 = n1min, n1max
        do i2 = n2min, n2max
          do i3 = n3min, n3max
ccccccccc TEST
c      do i1 = n1min-3, n1max+3
c      do i2 = n2min-3, n2max+3
c      do i3 = n3min-3, n3max+3
            qpg(1:3)= q(1:3) + 
     &      qlat(1:3,1)*i1 +qlat(1:3,2)*i2 +qlat(1:3,3)*i3
            if( sum(qpg(1:3)**2) < QpGmax2) then
              ig = ig+1
              if(job==2) ngvec(1:3,ig) = (/i1,i2,i3/)
              if(job==1) imx=max(imx,abs(i1),abs(i2),abs(i3))
c          if(job==2) write(1116,'(f8.4,3i3)')
c     &      tpiba*sqrt(sum(qpg(1:3)**2)),ngvec(1:3,ig) ! check write
            endif
          enddo
        enddo
      enddo
      ng = ig
      if(job==1) ngvec(1,1)=imx !mar2012takao
      end
!!-------------------------------------------------------------     
      subroutine eprod(a,b, c)
c- c gives normalized normal vector for a and b.
      real(8) :: a(3),b(3),c(3),cnorm !anorm,bnorm,
      c(1)= a(2)*b(3)-a(3)*b(2)
      c(2)= a(3)*b(1)-a(1)*b(3)
      c(3)= a(1)*b(2)-a(2)*b(1)
c      anorm = sqrt(sum(a(1:3)**2))
c      bnorm = sqrt(sum(b(1:3)**2))
c      c = c/anorm/bnorm
      cnorm = sqrt(sum(c(1:3)**2))
      c = c/cnorm
      end
      subroutine sortea(ea,ieaord,n,isig)
c mini-sort routine.
      implicit real*8(a-h,o-z)
      real(8)::        ea(n)
ckino delete integer(4) def.      integer(4):: ieaord(n)
      integer:: ieaord(n),n,isig,itmp,i,ix
c sorting of ea
      isig = 1
      do i = 1,n
        ieaord(i) = i
      enddo
      do ix= 2,n
        do i=ix,2,-1
          if( ea(ieaord(i-1)) >ea(ieaord(i) ) ) then
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp
            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo
      end
!!-----------------------------------------

      
      module m_pairs
      contains
      subroutine pairs(nbas,nbasp,alat,plat,rmax,baspp,ipsp,nd,iltab,
ckino Dec.14.2011:       .pltab,nttab,iv_p_ontab,iv_p_oiax,mxcsiz)
ckino Dec.14.2011:       .pltab,nttab,iv_p_ontab,iv_a_oiax,mxcsiz)
     .pltab,nttab,iv_a_ontab,iv_a_oiax,mxcsiz)
C- Allocate memory for and create neighbor table for a crystal
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis (input)
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat (input)
Ci   rmax  :maximum range for connecting vector, in a.u.
Ci          All connecting vectors with length < rmax(i)+rmax(j)
Ci          are retained.  rmax may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ipsp(1);
Ci          see description of ipsp
Ci   baspp :basis vectors, doubly padded for planar geometry
Ci   ipsp  :index to which species each site belongs, for padded basis;
Ci          identifies which rmax is associated with each site. NB:
Ci          ipsp(1) = -1 => rmax is a global scalar, independent of site
Ci          ipsp(1) =  0 => rmax is site-, not species-dependent.
Ci          In either of these cases, ipsp is not used.
Ci   nd    :number of dimensions for which periodic boundary conditions
Ci          are used
Ci   iltab :iltab<0, has no effect.  Otherwise, see pltabp.
Ci   pltabp:include only pairs for which pltabp(jb)-pltabp(ib) <= iltab
Ci   mxcsiz:if nonzero, use in place of internal formula for mxnbr
Co Outputs
Co   nttab   :total number of pairs in neighbor table and iax
Co   w(ontab):ntab array; see pairc below where it is generated
Co   w(oiax) :iax array; see pairc, below where it is generated
Co   mxcsiz  :size of the largest cluster encountered
C ----------------------------------------------------------------------
      implicit none
      integer:: nbas , nbasp , ipsp(1) , pltab(1) , nttab , nd , iltab,i_data_size
      integer,allocatable :: iv_a_ontab(:)
      integer,allocatable :: iv_a_oiax(:)
      double precision alat,plat(9),rmax(1),baspp(3,1)
      double precision avw,avwsr,vol
      integer:: modep(3) , nbaspp , mxcsiz , mxnbr , i , niax , isw
      real(8) ,allocatable :: wk_rv(:)
      integer,allocatable:: iv_a_tmp(:)
      parameter (niax=10)
C ... Set up input for call to pairc
      nbaspp = 2*nbasp - nbas
C ... Estimate an upper bound to the size of the neighbor table
      avw = avwsr(plat,alat,vol,nbas)
      mxnbr = 3*(2*rmax(1)/avw)**3*nbasp
      if (mxcsiz .gt. 0) mxnbr = mxcsiz
      if (allocated(iv_a_ontab)) deallocate(iv_a_ontab)
      allocate(iv_a_ontab(abs(nbasp+1)))
      if (allocated(iv_a_oiax)) deallocate(iv_a_oiax)
      allocate(iv_a_oiax(niax*mxnbr))
      iv_a_oiax(:)=0
      allocate(wk_rv(3*mxnbr))
      do  10  i = 1, 3
        modep(i) = 2
        if (i .gt. nd) modep(i) = 0
   10 continue

C ... This makes the neighbor table
      nttab = mxnbr
      isw = 0
      call pairc ( 1 , nbasp , nbaspp , modep , isw , ipsp , alat ,
     .  plat , baspp , baspp , rmax , iltab , pltab , nttab , iv_a_ontab
     .  , iv_a_oiax , wk_rv , mxcsiz )
C ... Allocate iax to proper size
      i_data_size=size(iv_a_oiax); allocate(iv_a_tmp(i_data_size))
      iv_a_tmp=iv_a_oiax; deallocate(iv_a_oiax)
      i_data_size=min(i_data_size,niax*nttab); allocate(iv_a_oiax(niax*nttab))
      iv_a_oiax(:i_data_size)=iv_a_tmp(:i_data_size); deallocate(iv_a_tmp)
      if (allocated(wk_rv)) deallocate(wk_rv)
      end subroutine pairs
      end module m_pairs

      subroutine pairc(ib1,ib2,nbasp,mode,isw,ips,alat,plat,pos,ctr,
     .range,iltab,pltabp,nttab,ntab,iax,rtab,mxcsiz)
C- Make a neighbor table (crystal version)
C ----------------------------------------------------------------
Ci Inputs:
Ci  ib1,ib2:range of sites for which to generate tables
Ci   nbasp :the size of the basis, plus possible extensions.
Ci          Usually nbasp=nbas, but will differ in special
Ci          cases, such as having padding sites to extend
Ci          to a semi-infinite geometry.
Ci   mode:  vector of length 3 governing how pos shortened (see shorps)
Ci   isw:   1's digit fixes how range is calculated.
Ci           0: vector length must be < range(i)+range(j)
Ci           1: include all connecting vecs w/ r < range(i)
Ci         10's digit sets what part of iax table is not calculated
Ci           1: do not calculate iax(6)
Ci              (may be needed when ctr and pos are different)
Ci           2: calculate only iax(1..5)
Ci   ips   :index to which species each site belongs, for padded basis;
Ci          identifies which rmax is associated with each site. NB:
Ci          ips(1) = -1 => rmax is a global scalar, independent of site
Ci          ips(1) =  0 => rmax is site-, not species-dependent.
Ci          In either of these cases, ips is not used.
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat (input)
Ci   pos   :site positions (doubly padded for planar geometry)
Ci   ctr   :ctr(1..3,ib) is the effective center of the cluster
Ci          associated with site ib for purposes of computing distance
Ci          pos(jb)-ctr(ib).  May point to the same address space as pos
Ci   range :maximum range for connecting vector, in a.u..
Ci          This quantity may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ips(1);
Ci          see description of ips.  Precisely what meaning range has
Ci          depends on mode and isw.
Ci   iltab :iltab<0, has no effect.  Otherwise, see pltabp.
Ci   pltabp:include only pairs for which pltabp(jb)-pltabp(ib) <= iltab
Ci   nttab :maximum dimension of iax table; used to guard against
Ci          generated table size exceeding dimension of iax.
Co Outputs:
Co   nttab    :total number of pairs generated
Co   iax      :neighbor table containing information about each pair ip
Co             For each pair ip, information is contained in iax(*,ip).
Co             as described below.  iax is ordered grouped by the basis
Co             atoms, so that all pairs connected to site ib are grouped
Co             together.  For each pair ip, iax(*,ip) contains:
Co   iax(1)   :site index to basis atoms ib=source;
Co             all pairs with common ib are contiguous
Co   iax(2)   :site index to jb=field of each pair
Co   iax(3..5):multiples of plat added the difference in site positions
Co             that connect the pair.
Co   iax(6)   :index to conjugate (jb,ib) pair matching (ib,jb)
Co             NB: no matching pairs outside (ib1..ib2) can be found.
Co   iax(7)   :permutation index ordering cluster by increasing
Co             effective site index; see ppair4.f
Co   iax(8)   :left untouched by pairc
Co   iax(9)   :left untouched by pairc
Co   iax(10)  :effective site index; see siteid.f
Co   ntab     :ntab(ib)=number of pairs in iax table preceding ib
Co             ntab is created for ib1:ib2+1.
Co   rtab     :rtab(1..3,ip) = pos(jb)-ctr(ib) for pair ip
Co   mxcsiz   :the largest cluster encountered
Cr Remarks
Cr   For each site ib=ib1..ib2, pairc finds all connecting vectors
Cr   for a lattice of points with periodic boundary conditions in
Cr   1, 2, or 3 dimensions, within a specified range of site ib.
Cr   The range can be defined in various ways, depending on isw.
Cu Updates
Cu   23 Apr 02 added option to make only iax(1..5) (isw=20)
C ----------------------------------------------------------------
      implicit none
      integer ib1,ib2,nbasp,mode(3),isw,nttab,niax,ips(nbasp),
     .ntab(ib1:ib2+1),iltab,pltabp(nbasp),mxcsiz
      parameter (niax=10)
      integer iax(niax,1),n,itmp,ixo
      double precision alat,plat(3,3),pos(3,nbasp),ctr(3,ib2),
     .range(1),rtab(3,1)
C Local variables
      integer:: ib , is , jb , mtab , i , moder , mode2(3) , nlat , 
     .mxntab , nsite
      real(8) ,allocatable :: wk1_rv(:)
      real(8) ,allocatable :: wk2_rv(:)
      real(8) ,allocatable :: wk3_rv(:)
      real(8) ,allocatable :: pos_rv(:)
      real(8) ,allocatable :: lat_rv(:)
      real(8) ,allocatable :: ctr_rv(:)

      double precision r1,rr,qlat(3,3),p0(3)
      integer ::iwdummy


C     call tcn('ppair1')

C --- Setup ---
      nsite = ib2-ib1+1
      mxntab = nttab
      moder = mod(isw,10)
      do  3  i = 1, 3
        mode2(i) = mode(i)
        if (mode2(i) .eq. 1) mode2(i) = 0
    3 continue
C ... Make r1 = 2*maximum range
      r1 = range(1)
      if (ips(1) .ge. 0) then
        do  5  ib = 1, nbasp
          is = ib
          if (ips(1) .gt. 0) then
            is = ips(ib)
          endif
          r1 = max(r1,range(is))
    5   continue
      endif
      if (moder .eq. 0) r1 = 2*r1
      r1 = 2*r1
C ... List of lattice vectors to add to pos(ib)-pos(jb)
      call xlgen ( plat , r1 / alat , 0d0 , 0 , 20 , mode , i , iwdummy 
     . )

      allocate(lat_rv(3*i))

      call xlgen ( plat , r1 / alat , 0d0 , i , 0 , mode , nlat , lat_rv 
     .)

C ... qlat = (plat^-1)^T so that qlat^T . plat = 1
*     call prmx('plat',plat,3,3,3)
*     call prmx('starting pos',pos,3,3,nbasp)
      call mkqlat(plat,qlat,rr)
C ... Save true pos in opos
C     and ctr in octr in case same address space used for ctr
      allocate(pos_rv(3*nbasp))

      call dpcopy ( pos , pos_rv , 1 , 3 * nbasp , 1d0 )

      allocate(ctr_rv(3*nsite))

      call dpcopy ( ctr , ctr_rv , 1 , 3 * nsite , 1d0 )


C --- For each ib, find all pairs for which dr < range ---
      nttab = 1
      ntab(ib1) = 0
      mtab = 1
      mxcsiz = 0
      do  10  ib = ib1, ib2
        r1 = range(1)
        if (ips(1) .ge. 0) then
          is = ib
          if (ips(1) .gt. 0) then
            is = ips(ib)
          endif
          r1 = range(is)
        endif

C --- Shorten all pos relative to ctr(ib) ---
C ... Make pos-ctr(ib)
        call dpcopy ( pos_rv , pos , 1 , 3 * nbasp , 1d0 )

        call dpcopy ( ctr_rv , ctr , 1 , 3 * nsite , 1d0 )

        do  12  i = 1, 3
          p0(i)  = ctr(i,ib)
          do  14  jb = 1, nbasp
            pos(i,jb) = pos(i,jb) - p0(i)
   14     continue
   12   continue
C ... Shorten pos-ctr(ib)
        call shorps(nbasp,plat,mode2,pos,pos)
C ... Undo shift -ctr(ib) to restore shortened pos to absolute pos
        do  16  jb = 1, nbasp
        do  16  i = 1, 3
          pos(i,jb) = pos(i,jb) + p0(i)
   16   continue

C --- Find all sites in range of ctr ---
        call ppair2 ( nbasp , iltab , pltabp , moder , alat , qlat , 
     .  pos , p0 , range , ips , rtab , ib , r1 , nlat , lat_rv , 
     .  pos_rv , mxntab , nttab , iax )


C --- Sort table by increasing length ---
c        allocate(wk1_rv((nttab-mtab)*niax))
c        allocate(wk2_rv((nttab-mtab)))
c        allocate(wk3_rv((nttab-mtab)*3))
        call ppair3 ( nttab - mtab , iax(1, mtab ), rtab(1 , mtab))! , wk1_rv , wk2_rv , wk3_rv )
c        if (allocated(wk3_rv)) deallocate(wk3_rv)
c        if (allocated(wk2_rv)) deallocate(wk2_rv)
c        if (allocated(wk1_rv)) deallocate(wk1_rv)
C --- Cleanup for this ib ---
        mtab = nttab
        ntab(ib+1) = nttab-1
        mxcsiz = max(mxcsiz,ntab(ib+1)-ntab(ib))
   10 continue
      nttab = nttab-1
C     call awrit2('xx ntab %n:1i',' ',80,6,nbasp+1,ntab)

C --- Restore original pos,ctr ---
      call dpcopy ( pos_rv , pos , 1 , 3 * nbasp , 1d0 )

      call dpcopy ( ctr_rv , ctr , 1 , 3 * nsite , 1d0 )

      if (allocated(ctr_rv)) deallocate(ctr_rv)
      if (allocated(pos_rv)) deallocate(pos_rv)


C --- Fill out iax table ---
      call ppair1(isw,ib1,ib2,nbasp,ips,alat,plat,pos,range,
     .nttab,ntab,iax,mxcsiz)

      if (allocated(lat_rv)) deallocate(lat_rv)

      end
      subroutine ppair1(isw,ib1,ib2,nbasp,ips,alat,plat,pos,range,
     .nttab,ntab,iax,mxcsiz)
C- Fill out parts of the aix table
C ----------------------------------------------------------------
Ci  Inputs
Ci   isw   :1's digit fixes how range is calculated.
Ci           0: vector length must be < range(i)+range(j)
Ci           1: include all connecting vecs w/ r < range(i)
Ci         :10's digit sets what part of iax table is calculated
Ci           0: make iax(6),iax(7),iax(10)
Ci           1: make iax(7),iax(10)
Ci           2: no change to iax: printout only
Ci           4: just make iax(6)
Ci   ib1   :fill out iax table for pairs ntab(ib1)+1..ntab(ib2)
Ci   ib2   :fill out iax table for pairs ntab(ib1)+1..ntab(ib2)
Ci   nbasp :size of padded basis (not needed)
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   range :maximum range for connecting vector, in a.u..
Ci          This quantity may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ips(1);
Ci          see description of ips.  See 1s digit of isw for
Ci          how range is used.
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   mxcsiz:maximum cluster size (for printout only)
Co  Outputs
Co   iax(6)   :index to conjugate (jb,ib) pair matching (ib,jb)
Co             NB: only matching pairs within site list can be found.
Co   iax(7)   :permutation index ordering cluster by increasing
Co             effective site index; see ppair4.f
Co   iax(10)  :effective site index; see siteid.f
C ----------------------------------------------------------------
      implicit none
      integer isw,ib1,ib2,nbasp,nttab,niax,ips(nbasp),ntab(ib1:ib2)
      parameter (niax=10)
      integer iax(niax,1),mxcsiz
      double precision alat,plat(3,3),pos(3,19),range(1)
C Local variables
      integer:: ib , is , jb , js , ipr , i , j , moder , it , jt , 
     .iprint , lgunit , stdo , nsite , isw1
      real(8) ,allocatable :: pos_rv(:)
      integer ,allocatable :: iwk_iv(:)

      double precision r1,r2,rr,rcut,vlat(3),drr2,tol
      parameter (tol=1d-5)

      isw1 = mod(isw/10,10)
      ipr = iprint()
      moder = mod(isw,10)
      nsite = ib2-ib1+1
      stdo = lgunit(1)
      if (isw1 .eq. 2) goto 80
      if (isw1 .eq. 4) goto 71

C --- Set iax(7) to sort this cluster ---
      call ppair5(ib1,ib2,plat,pos,tol,ntab,iax)

C --- For each pair, find matching pair, store in iax(6) ---
   71 continue
      do  74  it = 1,  nttab
        iax(6,it) = 0
   74 continue
      if (mod(isw1,2) .eq. 0) then
        do  70  ib = ib1, ib2
        do  70  it = ntab(ib)+1, ntab(ib+1)
          if (iax(6,it) .ne. 0) goto 70
          jb = iax(2,it)
C   ... No matching pair for padded sites
          if (jb .lt. ib1 .or. jb .gt. ib2) goto 70
          do  72  jt = ntab(jb)+1, ntab(jb+1)
C         if (iax(2,jt) .eq. ib) then
C         endif
            if (iax(2,jt) .eq. ib .and.
     .        iax(3,it) .eq. -iax(3,jt) .and.
     .        iax(4,it) .eq. -iax(4,jt) .and.
     .        iax(5,it) .eq. -iax(5,jt))  then
              iax(6,it) = jt
              iax(6,jt) = it
              goto 73
            endif
   72     continue
          call fexit3(-1,1,' Exit -1 pairc: cannot find pair'//
     .      ' matching sites (%i,%i), pair %i',ib,jb,it-ntab(ib))
   73     continue
   70   continue
      endif
      if (isw1 .eq. 4) return

C ... Assign a unique id for every different site in the cluster table
      allocate(iwk_iv(nttab))
      allocate(pos_rv(3*nttab))
      call siteid ( iax , nsite , ntab , plat , pos , pos_rv , iwk_iv 
     ., i )
      if (allocated(pos_rv)) deallocate(pos_rv)
      if (allocated(iwk_iv)) deallocate(iwk_iv)


C --- Printout ---
   80 if (ipr .lt. 30) goto 91
      if (ipr .le. 40) write(stdo,'(1x)')
      if (ipr .gt. 40) write(stdo,332)
  332 format(/'  ib  jb',9x,'--- r(jb)-r(ib) ---',10x,
     .'d       -x-plat-  map ord  id')
      i = 0
      do  90  it = 1, nttab
        ib = iax(1,it)
        jb = iax(2,it)
        rr = dsqrt(drr2(plat,pos(1,ib),pos(1,jb),
     .  iax(3,it),iax(4,it),iax(5,it),vlat))
        r1 = range(1)
        r2 = range(1)
        if (ips(1) .ge. 0) then
          is = ib
          if (ips(1) .gt. 0) then
            is = ips(ib)
          endif
          r1 = range(is)
          js = jb
          if (ips(1) .gt. 0) then
            js = ips(jb)
          endif
          r2 = range(js)
        endif
        if (moder .eq. 0) rcut = r1+r2
        if (moder .eq. 1) rcut = r1
        if (ib .ne. i) then
          if (alat .ne. 1)
     .    write(stdo,345) ib,ntab(ib+1)-ntab(ib),rcut/alat,rcut
          if (alat .eq. 1) write(stdo,345) ib,ntab(ib+1)-ntab(ib),rcut
  345     format(' pairc, ib=',i3,':',i4,' neighbors in range',2f7.3)
        endif
        i = ib
        if (ipr .gt. 40) write(stdo,334) iax(1,it),iax(2,it),
     .  (vlat(j),j=1,3), rr, (iax(j,it), j=3,7),iax(10,it)
  334   format(i4,i4,3f11.6,f9.4,3x,3i3,i5,2i4)
   90 continue
   91 if (ipr .ge. 20) write(stdo,
     .'('' pairc:'',i8,'' pairs total'',i10,'' is max cluster size'')')
     .nttab, mxcsiz

C     call tcx('ppair1')

      end
      subroutine ppair2(nbas,iltab,pltabp,moder,alat,qlat,pos,ctr,range,
     .ips,rtab,ib,r1,nlat,lat,trupos,mxntab,nttab,iax)
C- Kernel of pairc to find all sites in range of ctr
      implicit none
      integer nbas,ib,iltab,ips(nbas),pltabp(nbas),niax,nlat,moder,
     .mxntab,nttab
      parameter (niax=10)
      integer iax(niax,1)
      double precision alat,ctr(3),pos(3,nbas),range(nbas),rtab(3,1)
      double precision qlat(3,3),trupos(3,nbas),lat(3,*),r1
C Local variables
      integer i,ilat,jb,js
      double precision r2,rr,rcut,vlat(3),xx,rcutba,dpos(3)

      do  20  jb = 1, nbas
        if (iltab .ge. 0) then
          if (abs(pltabp(jb)-pltabp(ib)) .gt. iltab) goto 20
        endif
        r2 = range(1)
        if (ips(1) .ge. 0) then
          js = jb
          if (ips(1) .gt. 0) then
            js = ips(jb)
          endif
          r2 = range(js)
        endif
        if (moder .eq. 0) rcut = r1+r2
        if (moder .eq. 1) rcut = r1
        rcutba = (rcut / alat)**2
        dpos(1) = pos(1,jb)-ctr(1)
        dpos(2) = pos(2,jb)-ctr(2)
        dpos(3) = pos(3,jb)-ctr(3)

C   --- For each (ib,jb,ilat), do ---
        do  30  ilat = 1, nlat

          if (nttab .gt. mxntab) call rxi(
     .    'pairc: table exceeds input maximum size,',mxntab)

C ...   Add to list if connecting vector within range
          rtab(1,nttab) = dpos(1) + lat(1,ilat)
          rtab(2,nttab) = dpos(2) + lat(2,ilat)
          rtab(3,nttab) = dpos(3) + lat(3,ilat)
          rr = rtab(1,nttab)**2+rtab(2,nttab)**2+rtab(3,nttab)**2

*        call awrit5('try ib,jb,ilat= %i %i %i rr=%;4d: %l',' ',80,
*     .    6,ib,jb,ilat,rr,rr.lt.rcut)

C   --- Add to iax table if this pair in range ---
          if (rr .lt. rcutba) then

C     ... vlat += shortening vector
            do  32  i = 1, 3
              rtab(i,nttab) = alat*rtab(i,nttab)
C           rtab(i,nttab) = alat*(rtab(i,nttab)+ctr(i)-pos(i,ib))
              vlat(i) = lat(i,ilat) + pos(i,jb) - trupos(i,jb)
   32       continue

C     ... iax table for this pair
            iax(1,nttab) = ib
            iax(2,nttab) = jb
            do  33  i = 1, 3
              xx = vlat(1)*qlat(1,i)+vlat(2)*qlat(2,i)+vlat(3)*qlat(3,i)
              iax(2+i,nttab) = nint(xx)
   33       continue
            nttab = nttab+1

          endif

   30   continue
   20 continue
      end
      subroutine ppair3(nttab,iax,rtab) !,iwk,iwk2,rwk)
C- Sort neighbor table by distance
      implicit none
      integer nttab,niax,iwk2(nttab),i,j,k
      parameter (niax=10)
      integer iax(niax,nttab),iwk(niax,nttab)
      double precision rtab(3,nttab),rwk(3,nttab)

      do  10  i = 1, nttab
        rwk(1,i) = rtab(1,i)
        rwk(2,i) = rtab(2,i)
        rwk(3,i) = rtab(3,i)
        do  12  k = 1, niax
          iwk(k,i) = iax(k,i)
   12   continue
   10 continue
      call dvshel(3,nttab,rtab,iwk2,11)
      do  20  i = 1, nttab
        j = iwk2(i)+1
        rtab(1,i) = rwk(1,j)
        rtab(2,i) = rwk(2,j)
        rtab(3,i) = rwk(3,j)
        do  22  k = 1, niax
          iax(k,i) = iwk(k,j)
   22   continue
   20 continue
      end
      subroutine ppair4(iclus,nclus,plat,pos,ctr,iwk,rtab,tol,iax)
C- Sort cluster by increasing (x,y,z) relative to its center
C ----------------------------------------------------------------
Ci Inputs
Ci   iclus,nclus: sort iax(iclus..nclus)
Ci   plat :primitive lattice vectors
Ci    pos :basis vectors
Ci    ctr :cluster origin:does not affect the ordering, but shifts rtab
Ci    iwk :integer work array of length nclus-iclus+1
Ci    tol :tolerance to which positions are considered coincident
Ci         tol<0 => sort iax by iax(1..5)
Co Outputs
Co   iax(7,iclus..nclus) orders the cluster by increasing (x,y,z)
Co         (or increasing iax(1..5) if tol < 0
Co   rtab  :connecting vectors rtab(1..3,ip) = pos(jb)-ctr
Co          for pair ip and jb=iax(2,ip)
Cr Remarks
Cr  Each cluster is sorted by increasing (x,y,z),
Cr  sorted by x first, then by y, then by z, thus guaranteeing that
Cr  all sites common to any pair of clusters are ordered the same.
C ----------------------------------------------------------------
C     implicit none
      integer iclus,nclus,niax,iwk(15)
      parameter (niax=10)
      integer iax(niax,1)
      double precision plat(3,3),pos(3,1),ctr(3),rtab(3,29),tol
      integer ic,jb,ic0,ix,ia2,i,j,k
C Local variables
      double precision dx
C     integer jx
C     double precision wk2(3,nclus*3)
      dx(ia2,i,j,k) = pos(ix,ia2) +
     .plat(ix,1)*i + plat(ix,2)*j + plat(ix,3)*k

      ic0 = 0
      do  12  ic = iclus, nclus
        jb = iax(2,ic)
        ic0 = ic0+1
        do  14  ix = 1,3
          rtab(ix,ic0) = dx(jb,iax(3,ic),iax(4,ic),iax(5,ic)) - ctr(ix)
   14   continue
   12 continue

      if (tol .lt. 0) then
        call ivheap(niax,nclus-iclus+1,iax(1,iclus),iwk,1)
      else
        call dvheap(3,nclus-iclus+1,rtab,iwk,tol,1)
      endif

      do  20  ic = iclus, nclus
        iax(7,ic) = iwk(ic-iclus+1)
   20 continue

C ... Debugging ...
C      call ivprm(niax,nclus-iclus+1,iax(1,iclus),wk2,iwk,.false.)
C      call yprm('iax',0,wk2,0,niax,niax,nclus-iclus+1)
C      call awrit2('iwk %n:1i',' ',180,6,nclus-iclus+1,iwk)
C      do  30  ic = iclus, nclus
C      ic0 = ic-iclus+1
C   30 iwk(ic0) = iwk(ic0)-1
C      call dvperm(3,nclus-iclus+1,rtab,wk2,iwk,.true.)
C      do  32  ic = iclus, nclus
C        ic0 = ic-iclus+1
C        print 346, ic,(rtab(jx,ic0)+pos(jx,iax(1,ic)), jx=1,3)
C  346   format(i4,3f11.6)
C   32 continue
      end
      subroutine ppair5(ib1,ib2,plat,pos,tol,ntab,iax)
C- Sort a range of clusters according to tol
C ----------------------------------------------------------------------
Ci Inputs
Ci  ib1,ib2:range of clusters to sort
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   tol   :tolerance; see ppair4
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Co Outputs
Co   iax   :iax(7) is set to order cluster; see ppair4.
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer niax,ib1,ib2
      parameter (niax=10)
      integer iax(niax,1),ntab(ib2)
      double precision plat(3,3),pos(3,*),tol
      integer:: ib , nttab
      real(8) ,allocatable :: wk1_rv(:)
      real(8) ,allocatable :: wk2_rv(:)

C heap:

C --- Set iax(7) to sort this cluster ---
      do  10  ib = ib1, ib2
        nttab = ntab(ib+1)-ntab(ib)
        allocate(wk1_rv(nttab))
        allocate(wk2_rv(nttab*3))
        call ppair4 ( ntab ( ib ) + 1 , ntab ( ib + 1 ) , plat , pos 
     .  , pos ( 1 , ib ) , wk1_rv , wk2_rv , tol , iax )
        if (allocated(wk2_rv)) deallocate(wk2_rv)
        if (allocated(wk1_rv)) deallocate(wk1_rv)

   10 continue
      end
      subroutine pair3c(offi,nvec,iaxc,iaxh,ntabh,ipiaxh)
C- Find connecting vectors within a cluster linked by a pair table
C ----------------------------------------------------------------------
Ci Inputs
Ci   offi :offset to starting point in iaxc for this cluster
Ci   nvec :number of sites in cluster for which to find connecting vectors
Ci   iaxc :cluster pair table.  Only iaxc(2) and iaxc(7) are used.
Ci         iaxc(2,*) = field site index for this table
Ci         iaxc(7,*) = list of permutations that orders sites
Ci                     in iaxc and iaxh in the same way; see Remarks.
Ci   iaxh :pair (eg hamiltonian) neighbor table (pairc.f).
Ci         Only iaxh(2) and iaxh(7) are used. iaxh must consist of a
Ci         series of clusters grouped by site index, as in the
Ci         standard formed described in pairc.  The site index
Ci         (kept in iaxh(1,*) in the standard form) and the site
Ci         index iaxc(2,*) must refer to the same basis atom.
Ci   ntabh:ntabh(ib)=number of pairs in iaxh preceding site ib.
Co Outputs
Co  ipiaxh:indices to connecting vectors within a cluster.
Co         ipaxph(ip,kp) points to entry in iaxh table connecting
Co         ip to kp. ipaxph(ip,kp)=0 => no such vector in iaxh exists.
Cr Remarks
Cr   pair3c generates in ipiaxh the information needed to connect
Cr   three centers in a crystal.  For each pair of connecting vectors
Cr   in iaxc, a third vector is sought in iaxh that connects the two.
Cr
Cr   A schematic depiction of association of variables with their sites:
Cr
Cr            . origin
Cr
Cr    .  jp=field point
Cr    \
Cr     \ <- connecting vector ipiaxh(jp,ip) = entry in iaxh
Cr      \
Cr       \.  ip, source point
Cr
Cr   For efficient implementation, this routine relies on sorting
Cr   information in both iaxc and iaxh.  iaxc(7,*) and iaxh(7,*) must
Cr   hold a list of permutations that sorts each cluster in a unique
Cr   manner, so that connecting vectors common to different clusters are
Cr   ordered in the same way in those two clusters; this is done in
Cr   in ppair4, called with tol<0.
C ----------------------------------------------------------------------
      implicit none
      integer offi,nvec,niax,ntabh(8),ipiaxh(nvec,nvec)
      parameter (niax=10)
      integer iaxc(niax,1),iaxh(niax,9)
      integer ii,iii,jj,jjj,ip,jp,ib,nc,ich,icp,ic0
      integer iiax(5)
      call iinit(ipiaxh,nvec*nvec)
C --- Loop over 3C pairs i,j in permuted order ---
      ii = 0
      do  30  iii = 1, nvec
   31   ii = ii+1
        ip = iaxc(7,ii+offi)
        if (ip .gt. nvec) goto 31
        ib = iaxc(2,ip+offi)
C   ... ic0 = offset corresponding to iaxh(ic) = 1st pair for ib
        ic0 = ntabh(ib)
C   ... ich = current offset corresponding to iaxh.  Because we loop
C       in sorted order, we need not reset ich for each new pair below.
        ich = ic0+1
        nc  = ntabh(ib+1)
        jj = 0
        do  40  jjj = 1, nvec
   41     jj = jj+1
          jp = iaxc(7,jj+offi)
          if (jp .gt. nvec) goto 41
C         Do upper triangle only
C         if (lx6 .and. jp .lt. ip) goto 40
          iiax(3) = iaxc(3,jp+offi) - iaxc(3,ip+offi)
          iiax(4) = iaxc(4,jp+offi) - iaxc(4,ip+offi)
          iiax(5) = iaxc(5,jp+offi) - iaxc(5,ip+offi)

C     ... Increment ich until no exhaust all pairs in iaxh for ib
C          if (jp .eq. 74 .and. ip .eq. 7) then
C            print *, 'test'
C          endif
          ich = ich-1
   42     ich = ich+1
C     ... Skip if no more vectors connecting this (ib,jb)
          if (ich .gt. nc) goto 40
          icp = ic0+iaxh(7,ich)

C     ... Increment jp until iaxh is at least as large as iaxc
C         and ich  until iaxc is at least as large as iaxh
          if (iaxc(2,offi+jp) .gt. iaxh(2,icp)) goto 42
          if (iaxc(2,offi+jp) .lt. iaxh(2,icp)) goto 40
          if (iiax(3) .gt. iaxh(3,icp)) goto 42
          if (iiax(3) .lt. iaxh(3,icp)) goto 40
          if (iiax(4) .gt. iaxh(4,icp)) goto 42
          if (iiax(4) .lt. iaxh(4,icp)) goto 40
          if (iiax(5) .gt. iaxh(5,icp)) goto 42
          if (iiax(5) .lt. iaxh(5,icp)) goto 40

C     ... A match was found
          ipiaxh(jp,ip) = icp
   40   continue
   30 continue

C --- Construct lower triangle of table ---
C      if (lx6) then
C        do  50  i = 1, nvec
C        do  50  j = i+1, nvec
C          is = ipiaxh(i,j)
C          if (is .gt. 0) ipiaxh(j,i) = iaxh(6,is)
C   50   continue
C      endif

      end

