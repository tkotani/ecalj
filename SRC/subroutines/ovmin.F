      module m_ovmin
      contains
      subroutine ovmin(sovmin,nbas,nbasp,alat,plat,rmax,rmt,dclabl,
     .     ips,mode,z,iv_a_ontab,iv_a_oiax,pos,iprtbl)
      use m_gradzr,only:gradzr
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
      integer:: nbas , nbasp , iprtbl
      integer,allocatable :: iv_a_ontab(:)
      integer,allocatable :: iv_a_oiax(:,:)
      double precision plat(3,3),pos(3,nbasp),rmax(1),rmt(1),z(1),alat
      character(8):: dclabl(*) !double precision
      integer ips(1),mode(3)
      character sovmin*(*)
      double precision alato,plato(9),xx
      integer:: nbaso , nbaspo , mxlst , nlst , modeo(3) , novl
      parameter (mxlst=256)
      integer ilst(mxlst)
      real(8) ::wdummy(1,1)
      double precision fovl,xtol,gtol,dxmn,dxmx,fovmx
      double precision wk(0:27)
      integer i1mach,isw,ir,i,j,j1,j2,ls,m,lstyle, !op
     .iv,parg,nlstc,mxint,nclass,ib,ic,iclbsj,maxit,ipr,n,lgunit
      character dc*1
      external mxint
      integer,allocatable:: olist(:)
      real(8),allocatable:: w_opos(:,:),w_oz(:),w_ormax(:),w_oips(:),w_op(:)
C --- Print out positions and overlaps ---
      call getpr(ipr)
      if (iprtbl .gt. 0) call ovlchk(nbas,nbasp,pos,alat,rmax,0d0,
     .dclabl,ips,mode,plat,fovmx,xx)
      call fovlp ( 1 , nbas , iv_a_ontab , iv_a_oiax , plat , pos , 
     .  ips , alat , rmax , z , 6d0 , 1d0 , .75d0 , .5d0 , fovmx , fovl 
     .  , novl )
      if (novl .eq. 0) novl = 1
      if (ipr .ge. 10 .or. iprtbl .gt. 0)
     .write(lgunit(1),"(' OVMIN:  fovl= ',f5.1,' <ovlp> = ',f5.1,'%',
     .'   max ovlp = ',f5.1,'%')") fovl/novl,(fovl/novl)**(1/6d0)*100,fovmx*100
C     --- Minimize overlaps wrt positions in list ---
      if (sovmin .ne. ' ') then
            call rx('ovmin: need to recover if necessary.takao')
c$$$C   ... Default values for gradzr call
c$$$        xtol = 2d-4
c$$$        gtol = 1d-5
c$$$        dxmn = 1d-6
c$$$        dxmx = .10d0
c$$$        maxit = 20
c$$$        isw = 10051
c$$$
c$$$        ls = len(sovmin)
c$$$        j1 = 1
c$$$        dc = sovmin(j1:j1)
c$$$        j1 = j1+1
c$$$        lstyle = 0
c$$$
c$$$C   ... Return here to resume parsing for arguments
c$$$   40   continue
c$$$        call nwordg(sovmin,0,dc//' ',1,j1,j2)
c$$$
c$$$C   ... Parse special arguments
c$$$        if (sovmin(j2+1:j2+1) .ne. ' ')  then
c$$$          m = j1-1
c$$$          i = parg('dxmx=',4,sovmin,m,ls,dc,1,1,iv,dxmx)
c$$$          m = j1-1
c$$$          i = parg('xtol=',4,sovmin,m,ls,dc,1,1,iv,xtol)
c$$$          m = j1-1
c$$$          i = parg('style=',2,sovmin,m,ls,dc,1,1,iv,lstyle)
c$$$          m = j1-1
c$$$          i = parg('maxit=',2,sovmin,m,ls,dc,1,1,iv,maxit)
c$$$          j1 = j2+2
c$$$          goto 40
c$$$        endif
c$$$
c$$$C   ... List of all sites to move
c$$$        if (lstyle .gt. 0) then
c$$$          nclass = mxint(nbas,ips)
c$$$c          call defi(olist, nclass)
c$$$          allocate(olist(nclass))
c$$$          call clist(lstyle,sovmin(j1:j2+1),dclabl,z,nclass,nlstc,
c$$$     .    olist)
c$$$          nlst = 0
c$$$          do  12  i = 1, nlstc
c$$$            ic = olist(i) !w(olist+i-1)
c$$$            do  14  j = 1, nbas
c$$$              ib = iclbsj(ic,ips,-nbas,j)
c$$$              if (ib .lt. 0) goto 12
c$$$              nlst = nlst+1
c$$$              ilst(nlst) = ib
c$$$   14       continue
c$$$   12     continue
c$$$          deallocate(olist)
c$$$        elseif (sovmin(j1:j1+1) .eq. 'z ' .or.
c$$$     .  sovmin(j1:j1+1) .eq. 'Z ') then
c$$$          nlst = 0
c$$$          do  10  ib = 1, nbasp
c$$$            ic = ips(ib)
c$$$            if (z(ic) .eq. 0) then
c$$$              nlst = nlst+1
c$$$              ilst(nlst) = ib
c$$$            endif
c$$$   10     continue
c$$$        else
c$$$          call mkilst(sovmin(j1:),nlst,ilst)
c$$$        endif
c$$$        call awrit2(' min wrt:  %n:1i',' ',80,i1mach(2),nlst,ilst)
c$$$        call awrit3(' setup:     xtol = %,2g   dxmx = %,2g   maxit = %i'
c$$$     .  ,' ',80,i1mach(2),xtol,dxmx,maxit)
c$$$        if (nlst .le. 0) then
c$$$          print *, 'no sites in list ... no minimization'
c$$$          return
c$$$        endif
c$$$
c$$$C  ...  set up static block for ovcall
c$$$        alato = alat
c$$$        nbaso = nbas
c$$$        nbaspo = nbasp
c$$$cki        ontabo = ontab
c$$$ckino Dec.14.2011 manual correction       iv_p_ontabo => iv_p_ontab
c$$$cki        oiaxo = oiax
c$$$ckino Dec.14.2011 manual correction        iv_p_oiaxo => iv_p_oiax
c$$$
c$$$c        call defrr(opos,3*nbasp)
c$$$        allocate(w_opos(3,nbasp))
c$$$        call dpcopy(pos,w_opos,1,3*nbasp,1d0)
c$$$        nclass = mxint(nbas,ips)
c$$$c        call defrr(oz,nclass)
c$$$        allocate(w_oz(nclass))
c$$$        call dpcopy(z,w_oz,1,nclass,1d0)
c$$$c        call defrr(ormax,nbasp)
c$$$        allocate(w_ormax(nbasp))
c$$$        call dpcopy(rmax,w_ormax,1,nbasp,1d0)
c$$$c        call defi(oips,nbasp)
c$$$        allocate(w_oips(nbasp))
c$$$        call icopy(nbasp,ips,1,w_oips,1)
c$$$        call icopy(3,mode,1,modeo,1)
c$$$        call dpcopy(plat,plato,1,9,1d0)
c$$$
c$$$C  ...  initialization for gradzr call
c$$$        n = 3*nlst
c$$$c        call defrr(op,-10*n)
c$$$c        call defrr(og,n)
c$$$        allocate(w_op(10*n))
c$$$        w_op=0d0
c$$$        ir = 0
c$$$        do  20  i = 1, nlst
c$$$          j = ilst(i)
c$$$          call dpscop(w_opos,w_op,3,3*j-2,3*i-2,1d0)
c$$$   20   continue
c$$$        xx = ovcall(n,0d0,w_op,ir,
c$$$     &  plato,alato,nbaso,nbaspo,nlst,ilst,w_opos,w_oz,
c$$$ckino Dec.14.2011 manual correction     .  w_ormax,w_oips,modeo,iv_p_ontabo,iv_p_oiaxo)
c$$$ckino Dec.14.2011:       .  w_ormax,w_oips,modeo,iv_p_ontab,iv_p_oiax)
c$$$ckino Dec.14.2011:       .  w_ormax,w_oips,modeo,iv_p_ontab,iv_a_oiax)
c$$$     .  w_ormax,w_oips,modeo,iv_a_ontab,iv_a_oiax)
c$$$
c$$$        call pshpr(ipr-5)
c$$$ 22     call gradzr ( n , w_op , wdummy , dxmn , dxmx , xtol , gtol 
c$$$     ., 1.0d0 , wk , isw , ir )
c$$$
c$$$        xx = ovcall(n,0d0,w_op,ir,
c$$$     &  plato,alato,nbaso,nbaspo,nlst,ilst,w_opos,w_oz,
c$$$ckino Dec.14.2011 manual correction     .  w_ormax,w_oips,modeo,iv_p_ontabo,iv_p_oiaxo)
c$$$ckino Dec.14.2011:       .  w_ormax,w_oips,modeo,iv_p_ontab,iv_p_oiax)
c$$$ckino Dec.14.2011:       .  w_ormax,w_oips,modeo,iv_p_ontab,iv_a_oiax)
c$$$     .  w_ormax,w_oips,modeo,iv_a_ontab,iv_a_oiax)
c$$$
c$$$        if (ir .lt. 0) goto 22
c$$$        call poppr
c$$$
c$$$C ...   Update positions
c$$$        do  30  i = 1, nlst
c$$$          j = ilst(i)
c$$$          call dpscop(w_op,pos,3,3*i-2,3*j-2,1d0)
c$$$   30   continue
c$$$
c$$$C --- Print out updated positions and overlaps ---
c$$$        print '(/'' OVMIN:  updated site positions:'')'
c$$$        if (iprtbl .gt. 0) call ovlchk(nbas,nbasp,pos,alat,rmax,0d0,
c$$$     .  dclabl,ips,mode,plat,fovmx,xx)
c$$$ckino Dec.14.2011:            call fovlp ( 1 , nbas , iv_p_ontab , iv_p_oiax , plat , pos ,
c$$$ckino Dec.14.2011:            call fovlp ( 1 , nbas , iv_p_ontab , iv_a_oiax , plat , pos ,
c$$$        call fovlp ( 1 , nbas , iv_a_ontab , iv_a_oiax , plat , pos , 
c$$$     .    ips , alat , rmax , z , 6d0 , 1d0 , .75d0 , .5d0 , fovmx , fovl 
c$$$     .    , novl )
c$$$
c$$$
c$$$        if (novl .eq. 0) novl = 1
c$$$        if (ipr .ge. 10)
c$$$     .  call awrit3(' minimized: fovl = %;6g   <ovlp> = %;1d%%'//
c$$$     .  '   max ovlp = %;1d%%',' ',80,
c$$$     .  i1mach(2),fovl/novl,(fovl/novl)**(1/6d0)*100,fovmx*100)
c$$$      endif
      endif
      end subroutine ovmin

c$$$      double precision function ovcall(n,x,p,ir,
c$$$     &  plato,alato,nbaso,nbaspo,nlst,ilst,w_opos,w_oz,
c$$$     .  w_ormax,w_oips,modeo,iv_a_ontabo,iv_a_oiaxo)
c$$$C- Generic function call for projection grad fovl in a spec'd direction
c$$$Ci x,p,ir see gradzr
c$$$C     implicit none
c$$$      integer ir,n
c$$$      double precision x,p(3*n)
c$$$      double precision alato,plato(9)
c$$$      integer:: nbaso , nbaspo , mxlst , nlst , modeo(3) , oposb , 
c$$$     .  novl , novlp , novlm
c$$$      integer,allocatable :: iv_a_ontabo(:)
c$$$      integer,allocatable :: iv_a_oiaxo(:,:)
c$$$      real(8) :: w_ormax(:)
c$$$      parameter (mxlst=256)
c$$$      integer ilst(mxlst)
c$$$      logical cmdopt
c$$$      integer j,i,ix,ipr,lgunit,novl0
c$$$      double precision fovl,ddot,dx,val,fovp,fovm,pos(3),xx,fov0
c$$$      character*120 outs
c$$$      parameter (dx=1d-4)
c$$$      real(8):: w_opos(*),w_oz(*),w_oips(*)
c$$$      real(8),allocatable:: w_oposb(:,:)
c$$$C ... Save pos, other initialization
c$$$      call getpr(ipr)
c$$$      allocate(w_oposb(3,nbaspo))
c$$$      call dpcopy(w_opos,w_oposb,1,3*nbaspo,1d0)
c$$$      call pshpr(0)
c$$$      do  12  i = 1, nlst
c$$$        j = ilst(i)
c$$$        call dpscop(p,w_opos,3,3*i-2,3*j-2,1d0)
c$$$ 12   continue
c$$$      call ovlchk(nbaso,nbaspo,w_opos,alato,w_ormax,0d0,0d0,
c$$$     .w_oips,modeo,plato,fovl,xx)
c$$$      call fovlp ( 1 , nbaso , iv_a_ontabo , iv_a_oiaxo , plato , w_opos 
c$$$     .  , w_oips , alato , w_ormax , w_oz , 6d0 , 1d0 , .75d0 , .5d0 
c$$$     .  , xx , fovl , novl )
c$$$      if (fovl .eq. 0) then
c$$$        print *, 'ovmin: no spheres overlap:'
c$$$        call poppr
c$$$        call ovlchk(nbaso,nbaspo,w_opos,alato,w_ormax,0d0,0d0,
c$$$     .  w_oips,modeo,plato,fovp,xx)
c$$$      endif
c$$$
c$$$C ... Gradient of fovl wrt pos
c$$$      do  201  i = 1, nlst
c$$$        j = ilst(i)
c$$$        call fovlp ( j , j , iv_a_ontabo , iv_a_oiaxo , plato , w_opos 
c$$$     .    , w_oips , alato , w_ormax , w_oz , 6d0 , 1d0 , .75d0 , .5d0 
c$$$     .    , xx , fov0 , novl0 )
c$$$      do  20  ix = 1, 3
c$$$        val = p(3*i-3+ix)
c$$$        call dvset(w_opos,3*j-3+ix,3*j-3+ix,val+dx)
c$$$        call fovlp ( j , j , iv_a_ontabo , iv_a_oiaxo , plato , w_opos 
c$$$     .      , w_oips , alato , w_ormax , w_oz , 6d0 , 1d0 , .75d0 , .5d0 
c$$$     .      , xx , fovp , novlp )
c$$$        call dvset(w_opos,3*j-3+ix,3*j-3+ix,val-dx)
c$$$        call fovlp ( j , j , iv_a_ontabo , iv_a_oiaxo , plato , w_opos 
c$$$     .      , w_oips , alato , w_ormax , w_oz , 6d0 , 1d0 , .75d0 , .5d0 
c$$$     .      , xx , fovm , novlm )
c$$$        call dvset(w_opos,3*j-3+ix,3*j-3+ix,val)
c$$$        fovp = fovl + 2*(fovp-fov0)
c$$$        fovm = fovl + 2*(fovm-fov0)
c$$$        p(n+3*i-3+ix) = dlog(fovp/fovm)/2/dx
c$$$   20 continue
c$$$ 201  continue
c$$$      ovcall = ddot(n,p(n+1),1,p(2*n+1),1)
c$$$      if (ipr .ge. 50) then
c$$$        call awrit5('  ovcall: x=%d  f %;4g  lf %;4g  |glf| %;4g  '//
c$$$     .  'glf.x %;4g',' ',80,lgunit(1),x,fovl/novl,dlog(fovl/novl),
c$$$     .  dsqrt(ddot(n,p(n+1),1,p(n+1),1)),ddot(n,p(n+1),1,p(2*n+1),1))
c$$$        call awrit5('  ovcall: x=%d  f %;4g  lf %;4g  |glf| %;4g  '//
c$$$     .  'glf.x %;4g',' ',80,lgunit(2),x,fovl/novl,dlog(fovl/novl),
c$$$     .  dsqrt(ddot(n,p(n+1),1,p(n+1),1)),ddot(n,p(n+1),1,p(2*n+1),1))
c$$$        do  30  i = 1, nbaspo
c$$$          call dpscop(w_opos,pos,3,3*i-2,1,1d0)
c$$$          write(lgunit(2),140) pos
c$$$  140     format(3f12.6)
c$$$   30   continue
c$$$      endif
c$$$C ... restore pos
c$$$      call dpcopy(w_oposb,w_opos,1,3*nbaspo,1d0)
c$$$      call poppr
c$$$      end function ovcall
      end module m_ovmin


c$$$      subroutine clist(lstyle,slist,dclabl,z,nclass,nlist,list)
c$$$C- Generates a list of classes from a string specification
c$$$C ----------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   slist:  string specifying list of classes
c$$$Ci   lstyle: style of slist specification; see Remarks
c$$$Ci   nclass  number of classes.
c$$$Co Outputs
c$$$Co   nlist,list a list of classes satisfying specifications
c$$$Cr Remarks
c$$$Cr *Syntax of slist: depends on one of three styles (lstyle)
c$$$Cr
c$$$Cr *lstyle=1 : a list of integers; see mkilst.f for complete syntax.
c$$$Cr             Example: '1,4:6,11' generates a list of five numbers,
c$$$Cr             1,4,5,6,11.
c$$$Cr
c$$$Cr *lstyle=2 : the list is specified according to an expression.
c$$$Cr             The expression can involve the class index ic and
c$$$Cr             atomic number z.  Any class satisfying expression is
c$$$Cr             included in the list.  Example:  'ic<6&z==14'
c$$$Cr
c$$$Cr *lstyle=3 : is specifically for unix systems.  slist is a filename
c$$$Cr             with the usual unix wildcards, eg a[1-6].  'clist'
c$$$Cr             makes a system 'ls' call for that string; any class
c$$$Cr             which 'ls' finds is included in the list.
c$$$C ----------------------------------------------------------------
c$$$C     implicit none
c$$$      integer lstyle,nlist,nclass,list(1)
c$$$      character*(*) slist
c$$$      double precision z(nclass)
c$$$C Local variables
c$$$      integer iv0,ic,ival,is,i,j,ls,ls1
c$$$      logical a2bin,sw
c$$$      character strn*120,filnam*72,cnam*72
c$$$      character(8):: clabl,dclabl(nclass)
c$$$C Heap:
c$$$c      integer w(1)
c$$$c      common /w/ w
c$$$      integer,allocatable:: oilst(:)
c$$$
c$$$      ls = len(slist)
c$$$      nlist = 0
c$$$c      goto (10,20,30) lstyle
c$$$      select case(lstyle*10)
c$$$C -- lstyle=1 ---
c$$$c   10 continue
c$$$      case(10)
c$$$      call mkils0(slist,nlist,i)
c$$$c      call defi(oilst, nlist)
c$$$      allocate(oilst(nlist))
c$$$      call mkilst(slist,nlist,oilst)
c$$$      if (nlist .eq. 0) return
c$$$      call ishell(nlist,oilst)
c$$$c      list(1) = w(oilst)
c$$$      list(1) = oilst(1)
c$$$      j = 1
c$$$      do  12  i = 2, nlist
c$$$c        if (w(oilst+i-1) .gt. list(j)
c$$$c     .      .and. w(oilst+i-1) .le. nclass) then
c$$$        if (oilst(i) .gt. list(j)
c$$$     .    .and. oilst(i) .le. nclass) then
c$$$          j = j+1
c$$$c          list(j) = w(oilst+i-1)
c$$$          list(j) = oilst(i)
c$$$        endif
c$$$   12 continue
c$$$      nlist = j
c$$$c       call rlse(oilst)
c$$$      return
c$$$
c$$$C --- lstyle=2 ---
c$$$c   20 continue
c$$$      case(20)
c$$$      call numsyv(iv0)
c$$$      nlist = 0
c$$$      do  42  ic = 1, nclass
c$$$        call lodsyv('ic',1,dble(ic),ival)
c$$$        call lodsyv('z',1,z(ic),ival)
c$$$        is = 0
c$$$        if (a2bin(slist,sw,0,0,slist(ls:ls),is,ls)) then
c$$$          if (sw) then
c$$$            nlist = nlist+1
c$$$            list(nlist) = ic
c$$$          endif
c$$$C   ... Abort if a2bin can't parse expression
c$$$        else
c$$$          call rxs('clist: failed to parse',slist)
c$$$        endif
c$$$   42 continue
c$$$      call clrsyv(iv0)
c$$$      return
c$$$
c$$$c$$$C --- lstyle=3 ---
c$$$c$$$c   30 continue
c$$$c$$$      case(30)
c$$$c$$$      nlist = 0
c$$$c$$$      call skpblb(slist,ls,ls1)
c$$$c$$$      call ffnam(slist(1:ls1+1),filnam)
c$$$c$$$      do  44  ic = 1, nclass
c$$$c$$$        call pshpr(0)
c$$$c$$$c        call r8tos8(dclabl(ic),clabl)
c$$$c$$$        clabl=dclabl(ic)
c$$$c$$$        call ffnam(clabl,cnam)
c$$$c$$$        call poppr
c$$$c$$$        call awrit0('%xls ' // filnam //'%a|grep -s '
c$$$c$$$     .    // cnam // '%a>/dev/null',strn,len(strn),0)
c$$$c$$$        call locase(strn)
c$$$c$$$        call fsystm(strn,j)
c$$$c$$$        if (j .eq. 0) then
c$$$c$$$          nlist = nlist+1
c$$$c$$$          list(nlist) = ic
c$$$c$$$        endif
c$$$c$$$   44 continue
c$$$c$$$
c$$$      case default
c$$$      call rxi('clist: bad style',lstyle)
c$$$      return
c$$$
c$$$      end select
c$$$
c$$$      end
c$$$
c$$$      
