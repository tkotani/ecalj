CSFCPP#define F90 1
      subroutine gvlst2(alat,plat,q,n1,n2,n3,gmin,gmax,mshlst,job,ngmx,
     .     ng,kv,gv,igv,igv2)
      use m_ftox
      use m_shortn3,only: shortn3_initialize,shortn3
      use m_lgunit,only:stdo
C- Set up a list of recip vectors within cutoff |q+G| < gmax
C ----------------------------------------------------------------------
Ci Inputs
Ci   alat     Lattice constant
Ci   plat     Real-space primitive lattice vectors
Ci   q        vectors |q+G|<gmax are included in list
Ci   mshlst   if first entry is nonzero, a list of allowed values
Ci            n1,n2,n3 may take.  First entry is the size of
Ci            the list; then follows the list itself.
Ci   job      1s digit
Ci              0 return ng only
Ci              1 return kv and igv
Ci              2 return kv and igv2
Ci              4 return kv and gv
Ci              8 return kv and gv, and sort list
Ci                any combination of 1,2,4 is allowed
Ci            100s digit
Ci              1 to return internally generated values for n1,n2,n3;
Ci                see description of n1,n2,n3 below.
Ci            1000s digit
Ci              2  not used ... OLD Do not change input gmax if it is nonzero.
Ci   ngmx     Leading dimension of kv,gv,igv
Ci   gmin     Lower bound for reciprocal lattice vectors, in a.u.
Cio Inputs/Outputs
Cio   gmax    On input, cutoff for reciprocal lattice vectors, in a.u.
Cio           Energy cutoff is gmax**2.
Cio           If input gmax is zero, gvlst2 will generate it from n1..n3
Cio           (It is an error for both gmax and n1..n3 to be zero.)
Cio           On output, gmax may be altered; see Remarks.
Cio  n1..3    On input, max # divisions along the three lattice vectors.
Cio           (It is an error for both gmax and n1..n3 to be zero.)
Cio           Otherwise, input n1..n3 additionally constrain which
Cio           vectors are added to list; see Remarks.
Cio           On output, any n1..n3 initially zero are found from gmax
Co Outputs
Co   ng       Number of lattice vectors
Co   gv       list of reciprocal lattice vectors G
Co   igv      list of reciprocal lattice vectors G, represented as
Co            three integers (the multiples of qlat)
Co            gv and igv are related by:
Co              gv(1:3,1:ng) = 2*pi/alat * (qlat * igv(1:ng))
Co   igv2     same as igv except first and second columns are permuted
Co   kv       indices for gather/scatter operations.
Co            kv(ig,i=1,2,3) for vector ig point to which entry
Co            (i1,i2,i3) vector ig belongs
Cr Remarks
Cr   Collects a list of q + reciprocal lattice vectors (G+q) that lie
Cr   within a cutoff gmax.  List is optionally sorted by length.
Cr   Vectors G are integer multiples of primitive r.l.v. qlat.
Cr
Cr   Cutoff gmax may be input (preferred mode), or generated from
Cr   input values n1..n3.
Cr
Cr   Additional constraints may be imposed, namely that the number of
Cr   multiples in each of axes 1..3 not exceed input values n1..n3.
Cr
Cr   For copying the list of points to/from a regular mesh, use:
Cr     call gvgetf(ng,1,kv,k1,k2,k3,c,c0)
Cr       to collect elements into list c0 from 3D array c
Cr     call gvputf(ng,1,kv,k1,k2,k3,c0,c)
Cr       to poke elements from list c0 into 3D array c
Cr     call gvaddf(ng,kv,k1,k2,k3,c0,c)
Cr       to add elements from list c0 into 3D array c
Cr
Cu Updates
Cu   11 Jul 08 New argument gmin
Cu   01 Jun 01 revised gvlst2 together with gvlist.  They form similar
Cu             operations but with different functions; see gvlist.f
Cu   26 Mar 01 Another bug fix for input n1..n3 and gmax nonzero
Cu   06 Mar 01 Bug fix for input n1..n3 and gmax nonzero
Cu   Routine was adapted from T. Kotani, routine getgv2
C ----------------------------------------------------------------------
      implicit none
      integer n1,n2,n3,job,ng,ngmx,kv(ngmx,3),igv(ngmx,3),igv2(3,*),mshlst(0:*)
      double precision alat,gmin,gmax,gv(ngmx,3),plat(3,3),q(3)
      integer ig,n1max,n1min,n2max,n2min,n3max,n3min,i1,i2,i3,nn,i
      integer n1l,n2l,n3l
      integer PRTG,PRTG2,iset(3),ipr,job0,job1,job2,!job4,
     .job8,k1,k2,k3
      double precision qlat(3,3),vol,pi,tpiba,qpg(3),q2
      double precision gmin0,gmax0,gmin2,gmax2,h1,h2,h3,ddot,tol
      character*256 outs
      parameter (PRTG=30,PRTG2=100,tol=1d-8)
      logical::  lgpq,lgv,lsort,ligv,ligv2
      real(8):: plat1(3,3),qlat1(3,3),gg,gs(3)
      integer:: j1,j2,j3,noutmx,m,nout,jj1,jj2,jj3,nn1,nn2,nn3
      integer,allocatable:: nlatout(:,:)
      real(8):: rlatp(3,3),xmx2(3)
c      stdo = lgunit(1)
      call getpr(ipr)
      call dinv33(plat,1,qlat,vol)
      pi = 4d0*datan(1d0)
      tpiba = 2*pi/alat
      gmin0  = gmin/tpiba
      gmax0  = gmax/tpiba
      if (gmin .lt. 0) call rx('gvlst2: input gmin <= 0')
      if (gmax .le. 0) call rx('gvlst2: input gmax <= 0')
      job0 = mod(job,100)
      job1 = mod(job0,2)
      job2 = mod(job0/2,2)
      ligv  = mod(job0,2)/=0
      ligv2 = mod(job0/2,2)/=0
      lgv   = mod(job0/4,4)/=0
      lsort = mod(job0/8,2)/=0
      lgpq  = mod(job/100,10)>4 
!! ... Basis vectors for real-space mesh and recip-space supercell
      nn1=n1
      nn2=n2
      nn3=n3
      n1min=0; n1max=nn1-1
      n2min=0; n2max=nn2-1
      n3min=0; n3max=nn3-1
      if(nn1==0) call gvlstn(qlat(1,1),qlat(1,2),qlat(1,3),q,mshlst,gmax0,nn1) !,n1min,n1max)
      if(nn2==0) call gvlstn(qlat(1,2),qlat(1,3),qlat(1,1),q,mshlst,gmax0,nn2) !,n2min,n2max)
      if(nn3==0) call gvlstn(qlat(1,3),qlat(1,1),qlat(1,2),q,mshlst,gmax0,nn3) !,n3min,n3max)
      do m = 1, 3
        plat1(m,1) = plat(m,1)/nn1
        plat1(m,2) = plat(m,2)/nn2
        plat1(m,3) = plat(m,3)/nn3
        qlat1(m,1) = qlat(m,1)*nn1
        qlat1(m,2) = qlat(m,2)*nn2
        qlat1(m,3) = qlat(m,3)*nn3
      enddo
      n1l=nn1
      n2l=nn2
      n3l=nn3
c      print *,' jjjjjjj qlat qlat1=',qlat,qlat1
!! --- Loop through g vectors, shorten, count and keep if within gmax ---
      ig = 0
      gmax2 = (gmax0-tol)**2
      gmin2 = gmin0**2
      noutmx=48
      call shortn3_initialize(qlat1) !initialization for shoten3
      allocate(nlatout(3,noutmx))
c      print *,'wwwwwwwwwwwwwwwwwwwwwwwww  n1l =',n1l,n1min,n1max
c      print *,'wwwwwwwwwwwwwwwwwwwwwwwww  n2l =',n2l,n2min,n2max
c      print *,'wwwwwwwwwwwwwwwwwwwwwwwww  n3l =',n3l,n3min,n3max
c      print *,'wwwwwwwwwwwwwwwwwwwwwwwww  gmax2 =',gmax2
      do  212  j1 = 0,nn1-1 !n1min, n1max
      do  211  j2 = 0,nn2-1 !n2min, n2max
      do  21  j3 = 0,nn3-1 !n3min, n3max
          qpg= (/j1, j2, j3/) + matmul(q,plat(:,:)) 
          qpg= (/ qpg(1)/dble(nn1), qpg(2)/dble(nn2), qpg(3)/dble(nn3) /)
             !qpg on the supercell corrdinate of qlat1.
          call shortn3(qpg, noutmx, nout,nlatout)
          gs= matmul(qlat1(:,:), (qpg+nlatout(:,1))) 
          gg = (gs(1)**2+gs(2)**2+gs(3)**2)
cccccccccccccccccccccccccccccccccccccccc
c          write(*,"(a,2d13.5,3i5,3d13.5)") 'gggggggggggg gg=',gg,gmax2,nlatout(:,1),sum(abs(qlat1)),sum(abs(qlat1)),sum(abs(qpg))
cccccccccccccccccccccccccccccccccccccccc
              ! qpg+nlatout(:,ix) are shortest on the corrdinate of
              ! qlat1 (number of vectors are nlatout).
              ! We use nlatout(1:3,1).
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$          print *,'jjjjj nlatout',j1,j2,j3,nout,nlatout(:,1)
c$$$            do  i = 1, 3
c$$$              qpg(i)= q(i) + qlat(i,1)*j1 + qlat(i,2)*j2 + qlat(i,3)*j3
c$$$            enddo
c$$$            gg = qpg(1)**2+qpg(2)**2+qpg(3)**2
cccccccccccccccccccccccccccccccccccccccccccccccc
          k1 = j1+1 !mod(j1+n1l,n1l)+1
          k2 = j2+1 !mod(j2+n2l,n2l)+1
          k3 = j3+1 !mod(j3+n3l,n3l)+1
c          write(*,"('wwwww1111 j1j2j3=',7i4,f13.5)")ig,j1,j2,j3,k1,k2,k3,gg
          if(gmin2<= gg .and. gg < gmax2) then
            ig = ig+1
c            print *,' zzz ig=',ig,ngmx
            if (job0 .ne. 0) then
              if (ig .gt. ngmx) then
                 print *,' ig ngmx=',ig,ngmx
                 call rx('gvlist: ng exceeds ngmx')
              endif   
              kv(ig,1) = k1!j1+1
              kv(ig,2) = k2!j2+1
              kv(ig,3) = k3!j3+1
              if (ligv .or. ligv2) then
                jj1= j1 + nn1*nlatout(1,1)
                jj2= j2 + nn2*nlatout(2,1)
                jj3= j3 + nn3*nlatout(3,1)
                if (ligv) then
                  igv(ig,1) = jj1
                  igv(ig,2) = jj2
                  igv(ig,3) = jj3
                endif
                if (ligv2) then
                  igv2(1,ig) = jj1
                  igv2(2,ig) = jj2
                  igv2(3,ig) = jj3
                endif
              endif
              if (lgv) then
                if (lgpq) then
                  gv(ig,1) = gs(1)
                  gv(ig,2) = gs(2)
                  gv(ig,3) = gs(3)
                else
                  gv(ig,1) = gs(1) - q(1)
                  gv(ig,2) = gs(2) - q(2)
                  gv(ig,3) = gs(3) - q(3)
                endif
              endif
            endif
          endif
 21   continue
 211  continue
 212  continue
      deallocate(nlatout)
      ng = ig
C --- Printout ---
      if (ipr .ge. PRTG .and. n1l*n2l*n3l .eq. 0) then
         write(stdo,ftox)' GVLST2: gmax=',ftof(gmax0*tpiba,3),'created',ng,' recip. lattice vectors'
      elseif (ipr .ge. PRTG) then
        h1 = alat*sqrt(ddot(3,plat(1,1),1,plat(1,1),1))/n1l
        h2 = alat*sqrt(ddot(3,plat(1,2),1,plat(1,2),1))/n2l
        h3 = alat*sqrt(ddot(3,plat(1,3),1,plat(1,3),1))/n3l
        write(stdo,ftox)' GVLST2: gmax = ',ftof(gmax,3),
     .  ' a.u. created ',ng,' vectors of ',n1l*n2l*n3l,'(',(ng*100)/(n1l*n2l*n3l),'%)'
c$$$        iset(1) = n1
c$$$        iset(2) = n2
c$$$        iset(3) = n3
c$$$        i = iset(1)*iset(2)*iset(3)
c$$$        call awrit7('%a%N%9f%?#n#(input) ##mesh has %i x %i x %i'//
c$$$     .  ' divisions; length %,3;3d, %,3;3d, %,3;3d',outs,
c$$$     .  len(outs),0,i,n1l,n2l,n3l,h1,h2,h3)
c$$$        if (i .eq. 0 .and. iset(1)**2+iset(2)**2+iset(3)**2 .ne. 0) then
c$$$          call awrit3('%a%N%9fgenerated from input mesh with ('//
c$$$     .    '%?#n#%-1j%i#*#,%?#n#%-1j%i#*#,'//
c$$$     .    '%?#n#%-1j%i#*#) divisions',
c$$$     .    outs,len(outs),0,iset(1),iset(2),iset(3))
c$$$        endif
c$$$        call awrit0('%a',outs,len(outs),-stdo)
      endif
!! --- Sort the list of vectors --
      if(lsort) then
        call gvlsts(ngmx,ng,gv,kv,igv,igv2,job1,job2)
      endif
      if (ipr.ge.PRTG2 .and. ng.gt.0 .and. job1+job2.ne.0) then
        write(stdo,333)
  333   format(' G vectors (multiples of reciprocal lattice vectors)'/
     .  '   ig    G1   G2   G3     E')
        do  ig = 1, ng
          if (job1 .ne. 0) then
            i1 = igv(ig,1)
            i2 = igv(ig,2)
            i3 = igv(ig,3)
          endif
          if (job2 .ne. 0) then
            i1 = igv2(1,ig)
            i2 = igv2(2,ig)
            i3 = igv2(3,ig)
          endif
          do  i = 1, 3
            qpg(i)= q(i) + qlat(i,1)*i1 + qlat(i,2)*i2 + qlat(i,3)*i3
          enddo
          q2 = (qpg(1)**2+qpg(2)**2+qpg(3)**2) *tpiba**2
          write(stdo,334)  ig,i1,i2,i3,q2
          write(stdo,"('q  qpg=',3d13.5,3x,3d13.5)")  q,qpg
  334     format(i5,1x,3i5,2x,f8.4)
        enddo
      endif
      end

      subroutine gvlstn(q0,q1,q2,qp,mshlst,gmax0,nn) !,nmin,nmax)
C- Multiples of r.l.v. that bound cutoff gmax0  (r.l.v. reciplocal lattice vector)
C ----------------------------------------------------------------------
Ci Inputs
Ci   q0    :r.l.v. for which bounds nmin and nmax are to be computed
Ci   q1    :first  r.l.v. different from q0
Ci   q2    :second r.l.v. different from q0
Ci   qp    :k-point added to G vectors: sphere is centered G=qp.
Ci   mshlst:An ordered list of integers used to restrict the assignment
Ci         :of nn to one of an allowed list of points, should a value
Ci         :of nn be assigned (see description for nn below).
Ci         :If first entry is nonzero, mshlst(1..) = list of allowed
Ci         :values nn may take.  mshlst(0) is the size of mshlst.
Ci         :Certain fourier transform routines have restrictions
Ci         :on the allowed mesh sizes; this constraint is designed
Ci         :to handle that restriction.
Ci Inputs/Outputs
Ci   nn    :On input, if nonzero, maximum number of mesh points
Ci         :allowed.  Thus upper limit of nmax is nn/2; lower limit of nmin is -nn/2.
Ci         :On output, if input is zero, nn = 2*max(|nmin|,nmax)+1
Cio  gmax0 :On input, cutoff G
Cio        :On output, gmax0 may be reduced because constraints
Cio        :on nmin,nmax cause lattice vectors to be neglected
Cio        :that are smaller than input gmax0 (if input nn nonzero).
Co Outputs
Co     nn :   meshsize
Cxxo   nmin  :search for lattice vectors limited to (nmin..nmax)*q0
Cxxo   nmax  :search for lattice vectors limited to (nmin..nmax)*q0
Cr Remarks
Cr   q0,q1,q2,qp and G are all dimensionless (units of 2*pi/a)
Cu Updates
Cu   15 Apr 05 Bug fix when nn > max mshlst
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nmin,nmax,nn,mshlst(0:*)
      double precision q0(3),q1(3),q2(3),qp(3),gmax0
C ... Local parameters
      integer indx
      double precision qperp(3),ddot,qqperp
C ... qperp = q1 x q2 / |q1 x q2| ; qqperp = q . qperp
      qperp(1)  = q1(2)*q2(3) - q1(3)*q2(2)
      qperp(2)  = q1(3)*q2(1) - q1(1)*q2(3)
      qperp(3)  = q1(1)*q2(2) - q1(2)*q2(1)
C     anorm = sqrt(ddot(3,q1,1,q1,1))
C     bnorm = sqrt(ddot(3,q2,1,q2,1))
C     call dscal(3,1/anorm/bnorm,qperp,1)
      call dscal(3,1/sqrt(ddot(3,qperp,1,qperp,1)),qperp,1)
      qqperp = ddot(3,q0,1,qperp,1)
C   10 continue
C     print *, gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp + 1
C     print *, -gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp - 1
      nmax =  gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp + 1
      nmin = -gmax0/abs(Qqperp) - ddot(3,qp,1,qperp,1)/Qqperp - 1
C ... Assign nn, if input value is zero
      if (nn .eq. 0) then
        nn = 2*max(iabs(nmin),nmax)+1
        if (mshlst(0) .ne. 0) then
          indx = 1
          call hunti(mshlst(1),mshlst,nn,0,indx)
          nn = mshlst(min(indx+1,mshlst(0)))
        endif
c      else !bugfix added oct2013
c         print *,' nnnnn1111 nmax nmin=',nmax,nmin
c         nmax= min(nmax,  nn/2)
c         nmin= nmax-nn+1
c         print *,' nnnnn2222 nmax nmin=',nmax,nmin
      endif
      end

      subroutine gvlsts(ngmx,ng,gv,kv,igv,igv2,job1,job2)
C- Kernel called by gvlst2 to sort gv and kv
C     implicit none
      integer ngmx,ng,kv(ngmx,3),igv(ngmx,3),igv2(3,ngmx),job1,job2
      double precision gv(ngmx,3)
C Local variables
      integer ig,m,jg
CSFCPP#if F90 | AUTO_ARRAY
      integer kk(ngmx),iprm(ngmx)
      double precision gg(ngmx)
CSFCPP#else
CSFCPP      integer ngmxx
CSFCPP      parameter (ngmxx=20000)
CSFCPP      integer kk(ngmxx),iprm(ngmxx)
CSFCPP      double precision gg(ngmxx)
CSFCPP      if (ng .gt. ngmxx) call rxi('gvlst2: increase ngmx, need',ng)
CSFCPP#endif
      do  ig = 1, ng
        gg(ig) = gv(ig,1)**2 + gv(ig,2)**2 + gv(ig,3)**2
ctakao test case1 Apr2009
        gg(ig)= gg(ig) *(1d0 + 1d-15*ig)
      enddo
      call dvshel(1,ng,gg,iprm,1)
C     call dvheap(1,ng,gg,iprm,0d0,11)
C ... Rearrange gv,kv
      do  20  m = 1, 3
        do  22  ig = 1, ng
          jg = iprm(ig)+1
          gg(ig) = gv(jg,m)
          kk(ig) = kv(jg,m)
   22   continue
        do  24  ig = 1, ng
          gv(ig,m) = gg(ig)
          kv(ig,m) = kk(ig)
   24   continue
   20 continue
C ... Rearrange igv
      if (job1 .ne. 0) then
        do  30  m = 1, 3
          do  32  ig = 1, ng
            jg = iprm(ig)+1
            kk(ig) = igv(jg,m)
   32     continue
          do  34  ig = 1, ng
            igv(ig,m) = kk(ig)
   34     continue
   30   continue
      endif
C ... Rearrange igv2
      if (job2 .ne. 0) then
        do  40  m = 1, 3
          do  42  ig = 1, ng
            jg = iprm(ig)+1
            kk(ig) = igv2(m,jg)
   42     continue
          do  44  ig = 1, ng
            igv2(m,ig) = kk(ig)
   44     continue
   40   continue
      endif
C      do   ig = 1, ng
C        gg(ig) = gv(ig,1)**2 + gv(ig,2)**2 + gv(ig,3)**2
C        print 550, ig,gv(ig,1),gv(ig,2),gv(ig,3),
C     .     kv(ig,1),kv(ig,2),kv(ig,3),sqrt(gg(ig))
C  550   format(i5,3f11.5,3i6,f11.5)
C      enddo
C      pause
      end

