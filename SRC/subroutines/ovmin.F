      module m_ovmin
      contains
      subroutine ovmin(sovmin,nbas,nbasp,alat,plat,rmax,rmt,dclabl,
     .ips,mode,z,iv_a_ontab,iv_a_oiax,pos,iprtbl)
C- Check volume and sphere overlaps, optionally minimizing them
C ----------------------------------------------------------------
Ci Inputs
Ci   sovmin: a set of modifiers, with the syntax
Ci          -mino[:dxmx=#][:xtol=#][:style=#]:site-list
Ci   nbas  :size of basis
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   rmax  :potential radius, in a.u.
Ci   rmt   :augmentation radius, in a.u.
Ci   dclabl:class name, packed as a real number
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   mode:  vector of length 3 governing how pos shortened (see shorps)
Ci   z     :nuclear charge
Ci   w(ontab):ntab(ib)=# pairs in iax table preceding ib (pairc.f)
Ci   w(oiax):neighbor table containing pair information (pairc.f)
Ci   pos   :basis vectors
Ci   iprtbl: nonzero if to call ovlchk and print table of overlaps
Co Outputs
Co   Sphere overlaps are printed out
Cr Remarks
Cr   rmt(1)  not used now
Cu Updates
Cu   22 Oct 02  weight ES-ES and atom-ES overlaps differently when
Cu              minimizing sphere overlap positions
Cu    9 Dec 98  replace call to frpmin with call to gradzr.
Cu    8 Sep 98  small patches in minimizing algorithm
Cu   24 Nov 97  changed ovmin to call fovlp for fast execution
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer:: nbas , nbasp , iprtbl
ckino Dec.14.2011:          integer,pointer :: iv_p_ontab(:)
      integer,allocatable :: iv_a_ontab(:)

ckino Dec.14.2011:         integer,pointer :: iv_p_oiax(:)
      integer,allocatable :: iv_a_oiax(:,:)

C     integer ntab(nbas+1),iax(niax,1)
      double precision plat(3,3),pos(3,nbasp),rmax(1),rmt(1),z(1),alat
      character(8):: dclabl(*) !double precision
      integer ips(1),mode(3)
      character sovmin*(*)
C static:
      double precision alato,plato(9),xx
      integer:: nbaso , nbaspo , mxlst , nlst , modeo(3) , novl
ckino Dec.14.2011 manual correction        integer,pointer :: iv_p_ontabo(:) =>NULL()

ckino Dec.14.2011 manual correction       integer,pointer :: iv_p_oiaxo(:) =>NULL()

      parameter (mxlst=256)
      integer ilst(mxlst)
c      common /ovstat/ plato,alato,nbaso,nbaspo,nlst,ilst,opos,oz,
c     .ormax,oips,modeo,ontabo,oiaxo

C heap:
      integer ::iwdummy

C Local parameters
      double precision fovl,xtol,gtol,dxmn,dxmx,fovmx
      double precision wk(0:27)
      integer i1mach,isw,ir,i,j,j1,j2,ls,m,lstyle, !op
     .iv,parg,nlstc,mxint,nclass,ib,ic,iclbsj,maxit,ipr,n
      character dc*1
      external mxint
      integer,allocatable:: olist(:)
      real(8),allocatable:: w_opos(:,:),w_oz(:),w_ormax(:),w_oips(:),w_op(:)

C --- Print out positions and overlaps ---
      call getpr(ipr)
      if (iprtbl .gt. 0) call ovlchk(nbas,nbasp,pos,alat,rmax,0d0,
     .dclabl,ips,mode,plat,fovmx,xx)
ckino Dec.14.2011:          call fovlp ( 1 , nbas , iv_p_ontab , iv_p_oiax , plat , pos ,
ckino Dec.14.2011:          call fovlp ( 1 , nbas , iv_p_ontab , iv_a_oiax , plat , pos ,
      call fovlp ( 1 , nbas , iv_a_ontab , iv_a_oiax , plat , pos , 
     .  ips , alat , rmax , z , 6d0 , 1d0 , .75d0 , .5d0 , fovmx , fovl 
     .  , novl )


      if (novl .eq. 0) novl = 1
      if (ipr .ge. 10 .or. iprtbl .gt. 0)
     .call awrit3('%N OVMIN:     fovl = %;6g   <ovlp> = %;1d%%'//
     .'   max ovlp = %;1d%%',' ',80,
     .i1mach(2),fovl/novl,(fovl/novl)**(1/6d0)*100,fovmx*100)

C --- Minimize overlaps wrt positions in list ---
      if (sovmin .ne. ' ') then
C   ... Default values for gradzr call
        xtol = 2d-4
        gtol = 1d-5
        dxmn = 1d-6
        dxmx = .10d0
        maxit = 20
        isw = 10051

        ls = len(sovmin)
        j1 = 1
        dc = sovmin(j1:j1)
        j1 = j1+1
        lstyle = 0

C   ... Return here to resume parsing for arguments
   40   continue
        call nwordg(sovmin,0,dc//' ',1,j1,j2)

C   ... Parse special arguments
        if (sovmin(j2+1:j2+1) .ne. ' ')  then
          m = j1-1
          i = parg('dxmx=',4,sovmin,m,ls,dc,1,1,iv,dxmx)
          m = j1-1
          i = parg('xtol=',4,sovmin,m,ls,dc,1,1,iv,xtol)
          m = j1-1
          i = parg('style=',2,sovmin,m,ls,dc,1,1,iv,lstyle)
          m = j1-1
          i = parg('maxit=',2,sovmin,m,ls,dc,1,1,iv,maxit)
          j1 = j2+2
          goto 40
        endif

C   ... List of all sites to move
        if (lstyle .gt. 0) then
          nclass = mxint(nbas,ips)
c          call defi(olist, nclass)
          allocate(olist(nclass))
          call clist(lstyle,sovmin(j1:j2+1),dclabl,z,nclass,nlstc,
     .    olist)
          nlst = 0
          do  12  i = 1, nlstc
            ic = olist(i) !w(olist+i-1)
            do  14  j = 1, nbas
              ib = iclbsj(ic,ips,-nbas,j)
              if (ib .lt. 0) goto 12
              nlst = nlst+1
              ilst(nlst) = ib
   14       continue
   12     continue
          deallocate(olist)
        elseif (sovmin(j1:j1+1) .eq. 'z ' .or.
     .  sovmin(j1:j1+1) .eq. 'Z ') then
          nlst = 0
          do  10  ib = 1, nbasp
            ic = ips(ib)
            if (z(ic) .eq. 0) then
              nlst = nlst+1
              ilst(nlst) = ib
            endif
   10     continue
        else
          call mkilst(sovmin(j1:),nlst,ilst)
        endif
        call awrit2(' min wrt:  %n:1i',' ',80,i1mach(2),nlst,ilst)
        call awrit3(' setup:     xtol = %,2g   dxmx = %,2g   maxit = %i'
     .  ,' ',80,i1mach(2),xtol,dxmx,maxit)
        if (nlst .le. 0) then
          print *, 'no sites in list ... no minimization'
          return
        endif

C  ...  set up static block for ovcall
        alato = alat
        nbaso = nbas
        nbaspo = nbasp
cki        ontabo = ontab
ckino Dec.14.2011 manual correction       iv_p_ontabo => iv_p_ontab
cki        oiaxo = oiax
ckino Dec.14.2011 manual correction        iv_p_oiaxo => iv_p_oiax

c        call defrr(opos,3*nbasp)
        allocate(w_opos(3,nbasp))
        call dpcopy(pos,w_opos,1,3*nbasp,1d0)
        nclass = mxint(nbas,ips)
c        call defrr(oz,nclass)
        allocate(w_oz(nclass))
        call dpcopy(z,w_oz,1,nclass,1d0)
c        call defrr(ormax,nbasp)
        allocate(w_ormax(nbasp))
        call dpcopy(rmax,w_ormax,1,nbasp,1d0)
c        call defi(oips,nbasp)
        allocate(w_oips(nbasp))
        call icopy(nbasp,ips,1,w_oips,1)
        call icopy(3,mode,1,modeo,1)
        call dpcopy(plat,plato,1,9,1d0)

C  ...  initialization for gradzr call
        n = 3*nlst
c        call defrr(op,-10*n)
c        call defrr(og,n)
        allocate(w_op(10*n))
        w_op=0d0
        ir = 0
        do  20  i = 1, nlst
          j = ilst(i)
          call dpscop(w_opos,w_op,3,3*j-2,3*i-2,1d0)
   20   continue
        xx = ovcall(n,0d0,w_op,ir,
     &  plato,alato,nbaso,nbaspo,nlst,ilst,w_opos,w_oz,
ckino Dec.14.2011 manual correction     .  w_ormax,w_oips,modeo,iv_p_ontabo,iv_p_oiaxo)
ckino Dec.14.2011:       .  w_ormax,w_oips,modeo,iv_p_ontab,iv_p_oiax)
ckino Dec.14.2011:       .  w_ormax,w_oips,modeo,iv_p_ontab,iv_a_oiax)
     .  w_ormax,w_oips,modeo,iv_a_ontab,iv_a_oiax)

        call pshpr(ipr-5)
 22     call gradzr ( n , w_op , iwdummy , dxmn , dxmx , xtol , gtol 
     ., 1.0d0 , wk , isw , ir )

        xx = ovcall(n,0d0,w_op,ir,
     &  plato,alato,nbaso,nbaspo,nlst,ilst,w_opos,w_oz,
ckino Dec.14.2011 manual correction     .  w_ormax,w_oips,modeo,iv_p_ontabo,iv_p_oiaxo)
ckino Dec.14.2011:       .  w_ormax,w_oips,modeo,iv_p_ontab,iv_p_oiax)
ckino Dec.14.2011:       .  w_ormax,w_oips,modeo,iv_p_ontab,iv_a_oiax)
     .  w_ormax,w_oips,modeo,iv_a_ontab,iv_a_oiax)

        if (ir .lt. 0) goto 22
        call poppr

C ...   Update positions
        do  30  i = 1, nlst
          j = ilst(i)
          call dpscop(w_op,pos,3,3*i-2,3*j-2,1d0)
   30   continue

C --- Print out updated positions and overlaps ---
        print '(/'' OVMIN:  updated site positions:'')'
        if (iprtbl .gt. 0) call ovlchk(nbas,nbasp,pos,alat,rmax,0d0,
     .  dclabl,ips,mode,plat,fovmx,xx)
ckino Dec.14.2011:            call fovlp ( 1 , nbas , iv_p_ontab , iv_p_oiax , plat , pos ,
ckino Dec.14.2011:            call fovlp ( 1 , nbas , iv_p_ontab , iv_a_oiax , plat , pos ,
        call fovlp ( 1 , nbas , iv_a_ontab , iv_a_oiax , plat , pos , 
     .    ips , alat , rmax , z , 6d0 , 1d0 , .75d0 , .5d0 , fovmx , fovl 
     .    , novl )


        if (novl .eq. 0) novl = 1
        if (ipr .ge. 10)
     .  call awrit3(' minimized: fovl = %;6g   <ovlp> = %;1d%%'//
     .  '   max ovlp = %;1d%%',' ',80,
     .  i1mach(2),fovl/novl,(fovl/novl)**(1/6d0)*100,fovmx*100)
      endif

      end subroutine ovmin


      double precision function ovcall(n,x,p,ir,
     &  plato,alato,nbaso,nbaspo,nlst,ilst,w_opos,w_oz,
ckino Dec.14.2011:       .  w_ormax,w_oips,modeo,iv_p_ontabo,iv_p_oiaxo)
ckino Dec.14.2011:       .  w_ormax,w_oips,modeo,iv_p_ontabo,iv_a_oiaxo)
     .  w_ormax,w_oips,modeo,iv_a_ontabo,iv_a_oiaxo)

C- Generic function call for projection grad fovl in a spec'd direction
Ci x,p,ir see gradzr
C     implicit none
      integer ir,n
      double precision x,p(3*n)
C static:
      double precision alato,plato(9)
      integer:: nbaso , nbaspo , mxlst , nlst , modeo(3) , oposb , 
     .  novl , novlp , novlm
ckino Dec.14.2011:          integer,pointer :: iv_p_ontabo(:)
      integer,allocatable :: iv_a_ontabo(:)

ckino Dec.14.2011:         integer,pointer :: iv_p_oiaxo(:)
      integer,allocatable :: iv_a_oiaxo(:,:)

cki
      real(8) :: w_ormax(:)

      parameter (mxlst=256)
      integer ilst(mxlst)
c      common /ovstat/ plato,alato,nbaso,nbaspo,nlst,ilst,opos,oz,
c     .ormax,oips,modeo,ontabo,oiaxo
C heap:
C Local
      logical cmdopt
C     character*8 clabl(10)
      integer j,i,ix,ipr,lgunit,novl0
      double precision fovl,ddot,dx,val,fovp,fovm,pos(3),xx,fov0
      character*120 outs
      parameter (dx=1d-4)
      real(8):: w_opos(*),w_oz(*),w_oips(*)
      real(8),allocatable:: w_oposb(:,:)
C ... Save pos, other initialization
      call getpr(ipr)
      allocate(w_oposb(3,nbaspo))
      call dpcopy(w_opos,w_oposb,1,3*nbaspo,1d0)
      call pshpr(0)

      do  12  i = 1, nlst
        j = ilst(i)
        call dpscop(p,w_opos,3,3*i-2,3*j-2,1d0)
 12   continue
      
      call ovlchk(nbaso,nbaspo,w_opos,alato,w_ormax,0d0,0d0,
     .w_oips,modeo,plato,fovl,xx)
      call fovlp ( 1 , nbaso , iv_a_ontabo , iv_a_oiaxo , plato , w_opos 
     .  , w_oips , alato , w_ormax , w_oz , 6d0 , 1d0 , .75d0 , .5d0 
     .  , xx , fovl , novl )
      if (fovl .eq. 0) then
        print *, 'ovmin: no spheres overlap:'
        call poppr
        call ovlchk(nbaso,nbaspo,w_opos,alato,w_ormax,0d0,0d0,
     .  w_oips,modeo,plato,fovp,xx)
      endif

C ... Gradient of fovl wrt pos
      do  201  i = 1, nlst
        j = ilst(i)
        call fovlp ( j , j , iv_a_ontabo , iv_a_oiaxo , plato , w_opos 
     .    , w_oips , alato , w_ormax , w_oz , 6d0 , 1d0 , .75d0 , .5d0 
     .    , xx , fov0 , novl0 )
      do  20  ix = 1, 3
        val = p(3*i-3+ix)
        call dvset(w_opos,3*j-3+ix,3*j-3+ix,val+dx)
        call fovlp ( j , j , iv_a_ontabo , iv_a_oiaxo , plato , w_opos 
     .      , w_oips , alato , w_ormax , w_oz , 6d0 , 1d0 , .75d0 , .5d0 
     .      , xx , fovp , novlp )
        call dvset(w_opos,3*j-3+ix,3*j-3+ix,val-dx)
        call fovlp ( j , j , iv_a_ontabo , iv_a_oiaxo , plato , w_opos 
     .      , w_oips , alato , w_ormax , w_oz , 6d0 , 1d0 , .75d0 , .5d0 
     .      , xx , fovm , novlm )
        call dvset(w_opos,3*j-3+ix,3*j-3+ix,val)
        fovp = fovl + 2*(fovp-fov0)
        fovm = fovl + 2*(fovm-fov0)
        p(n+3*i-3+ix) = dlog(fovp/fovm)/2/dx
   20 continue
 201  continue
      ovcall = ddot(n,p(n+1),1,p(2*n+1),1)
      if (ipr .ge. 50) then
        call awrit5('  ovcall: x=%d  f %;4g  lf %;4g  |glf| %;4g  '//
     .  'glf.x %;4g',' ',80,lgunit(1),x,fovl/novl,dlog(fovl/novl),
     .  dsqrt(ddot(n,p(n+1),1,p(n+1),1)),ddot(n,p(n+1),1,p(2*n+1),1))
        call awrit5('  ovcall: x=%d  f %;4g  lf %;4g  |glf| %;4g  '//
     .  'glf.x %;4g',' ',80,lgunit(2),x,fovl/novl,dlog(fovl/novl),
     .  dsqrt(ddot(n,p(n+1),1,p(n+1),1)),ddot(n,p(n+1),1,p(2*n+1),1))
        do  30  i = 1, nbaspo
          call dpscop(w_opos,pos,3,3*i-2,1,1d0)
          write(lgunit(2),140) pos
  140     format(3f12.6)
   30   continue
      endif
C ... restore pos
      call dpcopy(w_oposb,w_opos,1,3*nbaspo,1d0)
c      call rlse(oposb)
      call poppr
      end function ovcall
      end module m_ovmin


      subroutine clist(lstyle,slist,dclabl,z,nclass,nlist,list)
C- Generates a list of classes from a string specification
C ----------------------------------------------------------------
Ci Inputs
Ci   slist:  string specifying list of classes
Ci   lstyle: style of slist specification; see Remarks
Ci   nclass  number of classes.
Co Outputs
Co   nlist,list a list of classes satisfying specifications
Cr Remarks
Cr *Syntax of slist: depends on one of three styles (lstyle)
Cr
Cr *lstyle=1 : a list of integers; see mkilst.f for complete syntax.
Cr             Example: '1,4:6,11' generates a list of five numbers,
Cr             1,4,5,6,11.
Cr
Cr *lstyle=2 : the list is specified according to an expression.
Cr             The expression can involve the class index ic and
Cr             atomic number z.  Any class satisfying expression is
Cr             included in the list.  Example:  'ic<6&z==14'
Cr
Cr *lstyle=3 : is specifically for unix systems.  slist is a filename
Cr             with the usual unix wildcards, eg a[1-6].  'clist'
Cr             makes a system 'ls' call for that string; any class
Cr             which 'ls' finds is included in the list.
C ----------------------------------------------------------------
C     implicit none
      integer lstyle,nlist,nclass,list(1)
      character*(*) slist
      double precision z(nclass)
C Local variables
      integer iv0,ic,ival,is,i,j,ls,ls1
      logical a2bin,sw
      character strn*120,filnam*72,cnam*72
      character(8):: clabl,dclabl(nclass)
C Heap:
c      integer w(1)
c      common /w/ w
      integer,allocatable:: oilst(:)

      ls = len(slist)
      nlist = 0
c      goto (10,20,30) lstyle
      select case(lstyle*10)
C -- lstyle=1 ---
c   10 continue
      case(10)
      call mkils0(slist,nlist,i)
c      call defi(oilst, nlist)
      allocate(oilst(nlist))
      call mkilst(slist,nlist,oilst)
      if (nlist .eq. 0) return
      call ishell(nlist,oilst)
c      list(1) = w(oilst)
      list(1) = oilst(1)
      j = 1
      do  12  i = 2, nlist
c        if (w(oilst+i-1) .gt. list(j)
c     .      .and. w(oilst+i-1) .le. nclass) then
        if (oilst(i) .gt. list(j)
     .    .and. oilst(i) .le. nclass) then
          j = j+1
c          list(j) = w(oilst+i-1)
          list(j) = oilst(i)
        endif
   12 continue
      nlist = j
c       call rlse(oilst)
      return

C --- lstyle=2 ---
c   20 continue
      case(20)
      call numsyv(iv0)
      nlist = 0
      do  42  ic = 1, nclass
        call lodsyv('ic',1,dble(ic),ival)
        call lodsyv('z',1,z(ic),ival)
        is = 0
        if (a2bin(slist,sw,0,0,slist(ls:ls),is,ls)) then
          if (sw) then
            nlist = nlist+1
            list(nlist) = ic
          endif
C   ... Abort if a2bin can't parse expression
        else
          call rxs('clist: failed to parse',slist)
        endif
   42 continue
      call clrsyv(iv0)
      return

c$$$C --- lstyle=3 ---
c$$$c   30 continue
c$$$      case(30)
c$$$      nlist = 0
c$$$      call skpblb(slist,ls,ls1)
c$$$      call ffnam(slist(1:ls1+1),filnam)
c$$$      do  44  ic = 1, nclass
c$$$        call pshpr(0)
c$$$c        call r8tos8(dclabl(ic),clabl)
c$$$        clabl=dclabl(ic)
c$$$        call ffnam(clabl,cnam)
c$$$        call poppr
c$$$        call awrit0('%xls ' // filnam //'%a|grep -s '
c$$$     .    // cnam // '%a>/dev/null',strn,len(strn),0)
c$$$        call locase(strn)
c$$$        call fsystm(strn,j)
c$$$        if (j .eq. 0) then
c$$$          nlist = nlist+1
c$$$          list(nlist) = ic
c$$$        endif
c$$$   44 continue
c$$$
      case default
      call rxi('clist: bad style',lstyle)
      return

      end select

      end

      
