c$$$#define USEOPTIONAL
c$$$
c$$$      module m_ioden2
c$$$      contains
c$$$      subroutine ioden2(mode,sopts,nsp,plat,alat,nbas,pos,z,k1,k2,k3,
c$$$     .ngabc,smrho,wk,nwk,modrhg)
c$$$
c$$$C- Kernel called by ioden
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   mode  :0 just return these arguments parsed from sopts:
c$$$Ci         :  modrhg
c$$$Ci         :1 copy and save density in appropriate plane
c$$$Ci   sopts :string containing options; see ioden above.
c$$$Ci   plat  :primitive lattice vectors, in units of alat (for printout)
c$$$Ci   alat  :lattice constant
c$$$Ci   nbas  :number of atoms in unit cell
c$$$Ci   pos   :atomic positions in unit cell
c$$$Ci   z     : atomic number of atoms in cell
c$$$Ci   nsp   :number of spin channels
c$$$Ci   k1..3 :dimensions smrho
c$$$Ci   ngabc :number of divisions in smrho along plat(1..3)
c$$$Ci   smrho :smooth density on uniform mesh
c$$$Ci   wk    :2D work array
c$$$Ci   nwk   :dimensions wk
c$$$Co Outputs
c$$$Co   modrhg:controls what part of core density is added
c$$$Co   The projection of smrho onto a plane is written to disk
c$$$Co   or if option is g3d the full 3d grid is written to disk
c$$$Cl Local variables
c$$$Cl   lic   :lic(1..3,1..2) = step length in each vector axes 1..2
c$$$Cl         :lic(4,1..2)    = number of points along axes 1..2
c$$$Cl         :lic(5..7,1)    = starting element for each lattice vector
c$$$Cr Remarks
c$$$Cr   Choosing p1, p2, p3 sets the origin for that plat and also
c$$$Cr   creates a default set of vectors, choosing the other lattice
c$$$Cr   vectors in the system.  Thus:
c$$$Cr   p1=# => lic(1..7,1) = (0 1 0 nb 1 1 1) lic(1..7,2) = (0 0 1 nc **)
c$$$Cr   p2=# => lic(1..7,1) = (0 0 1 nc 1 1 1) lic(1..7,2) = (1 0 0 na **)
c$$$Cr   p3=# => lic(1..7,1) = (1 0 0 na 1 1 1) lic(1..7,2) = (0 1 0 nb **)
c$$$Cu Updates
c$$$Cu   09 Feb 01 Revised how a plane is defined.
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      character*(*) sopts
c$$$      integer mode,ngabc(3),k1,k2,k3,nwk,nsp,modrhg,nbas
c$$$#ifdef USEOPTIONAL
c$$$      double complex ,optional:: smrho(k1,k2,k3,nsp)
c$$$      double precision plat(3,3),alat,pos(3,nbas)
c$$$      real(8),optional:: wk(nwk,nwk,nsp)
c$$$#else
c$$$      double complex smrho(k1,k2,k3,nsp)
c$$$      double precision plat(3,3),wk(nwk,nwk,nsp),alat,pos(3,nbas)
c$$$#endif
c$$$      double precision z(nbas)
c$$$C ... Local parameters
c$$$      logical lok
c$$$      integer lgunit,stdo,j2,j1,fopn,ifi,i
c$$$      integer lic(7,2),j,k,iv(3),jv(3),kv(4),a2vec,isp
c$$$      double precision vecs(3,2),vv,ddot
c$$$      character*120 dc*1, fn, prfmt*40, sout
c$$$      logical lg3d
c$$$      integer i1,i2,i3
c$$$
c$$$      logical:: l_dummy_isanrg,isanrg
c$$$
c$$$ccccccccccccccccccccccccccccc
c$$$c      print *,' zzxxx=',z
c$$$c      stop 'zzzzzzzzzzzzzzzzzzzzzzzzzzzz'
c$$$ccccccccccccccccccccccccccccc
c$$$
c$$$      stdo = lgunit(1)
c$$$      dc = sopts(1:1)
c$$$      fn = 'smrho'
c$$$C     prfmt = '(8f10.6)'
c$$$      prfmt = '(8e14.6)'
c$$$      lg3d = .false.
c$$$C     Default: origin at (1,1,1), (na+1,nb+1) points along 1st,2nd axes
c$$$      call iinit(lic,14)
c$$$      call ivset(lic,5,7,1)
c$$$      lic(1,1) = 1
c$$$      lic(4,1) = ngabc(1)+1
c$$$      lic(2,2) = 1
c$$$      lic(4,2) = ngabc(2)+1
c$$$
c$$$      if (dc .ne. ' ') then
c$$$C ... Return here to resume parsing for arguments
c$$$        j2 = 0
c$$$   10   continue
c$$$        j2 = j2+1
c$$$        if (sopts(j2:j2) .eq. dc) goto 10
c$$$        j1 = min(len(sopts),j2)
c$$$        call nwordg(sopts,0,dc//' ',1,j1,j2)
c$$$        if (j2 .ge. j1) then
c$$$          if (.false.) then
c$$$          elseif (sopts(j1:j1+2) .eq. 'fn=')  then
c$$$            if (j1+3 .gt. j2) goto 99
c$$$            fn = sopts(j1+3:j2)
c$$$C   ... option g3d
c$$$          elseif (sopts(j1:j1+2).eq. 'g3d') then
c$$$            lg3d = .true.
c$$$C   ... option o=#,#,#
c$$$          elseif (sopts(j1:j1+1) .eq. 'o=')  then
c$$$            if (j1+2 .gt. j2) goto 99
c$$$            i = j1+1
c$$$            if (a2vec(sopts,j2,i,2,', '//dc,3,2,3,kv,iv) .ne. 3)
c$$$     .      goto 99
c$$$            lic(5,1) = mod(iv(1)+ngabc(1),ngabc(1))+1
c$$$            lic(6,1) = mod(iv(2)+ngabc(2),ngabc(2))+1
c$$$            lic(7,1) = mod(iv(3)+ngabc(3),ngabc(3))+1
c$$$C   ... option l[12]=#,#,#,#
c$$$          elseif (sopts(j1:j1+2) .eq. 'l1=' .or.
c$$$     .    sopts(j1:j1+2) .eq. 'l2=')  then
c$$$            if (j1+3 .gt. j2) goto 99
c$$$            i = 0
c$$$            call chrps2(sopts(j1+1:),'12',2,0,i,iv)
c$$$C         this check should never be necessary
c$$$Ckino isanrg is logical function, C         call isanrg(iv,1,2,' ','iv in ioden ... bug ..',.true.)
c$$$C         l_dummy_isanrg=isanrg(iv,1,2,' ','iv in ioden ... bug ..',.true.)
c$$$            i = j1+2
c$$$            j = iv(1)
c$$$            k = a2vec(sopts,j2,i,2,', '//dc,3,2,4,kv,lic(1,j))
c$$$            if (k .ne. 3 .and. k .ne. 4) goto 99
c$$$          elseif (sopts(j1:j1+4) .eq. 'core=') then
c$$$            if (j1+5 .gt. j2) goto 99
c$$$            i = j1+4
c$$$            if (a2vec(sopts,j2,i,2,' '//dc,2,1,1,kv,modrhg) .ne. 1)
c$$$     .      goto 99
c$$$          else
c$$$            call rxs('ioden: unrecognised option ... ',sopts(j1:j2))
c$$$          endif
c$$$          goto 10
c$$$        endif
c$$$      endif
c$$$      if (mode .eq. 0) return
c$$$
c$$$C --- Printout ---
c$$$      if (lg3d) then
c$$$        ifi = fopn(fn)
c$$$        rewind ifi
c$$$        call awrit3('      Writing smooth density to file'//fn//
c$$$     .  ' %a : full 3d grid (%i,%i,%i).',' ',80,stdo,k1,k2,k3)
c$$$        do isp = 1, nsp
c$$$C           call awrit4('%% 3D block data %i %i %i, spin=%i',' ',80,ifi,
c$$$C     .     k1,k2,k3,isp)
c$$$          write(ifi,'("CRYSTAL")')
c$$$          write(ifi,'("PRIMVEC")')
c$$$          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3)
c$$$     .    ,i2=1,3)
c$$$          write(ifi,'("PRIMCOORD")')
c$$$          write(ifi,'(2i5)') nbas,1
c$$$          do i = 1, nbas
c$$$            write(ifi,'(i4,2x,3f10.5)') int(z(i)),
c$$$     .      (pos(i2,i)*alat*0.529177208,i2=1,3)
c$$$          enddo
c$$$          write(ifi,'("BEGIN_BLOCK_DATAGRID_3D")')
c$$$          write(ifi,'("charge_density_spin_",i1)') isp
c$$$          write(ifi,'("BEGIN_DATAGRID_3D_isp_",i1)') isp
c$$$          write(ifi,'(3i4)') k1,k2,k3
c$$$          write(ifi,'(3f10.5)') 0.,0.,0.
c$$$          write(ifi,'(3f10.5)') ((plat(i1,i2)*alat*0.529177208,i1=1,3)
c$$$     .    ,i2=1,3)
c$$$          write(ifi,'(8e14.6)')  
c$$$     .    (((dble(smrho(i1,i2,i3,isp)),i1=1,k1),i2=1,k2),i3=1,k3)
c$$$          write(ifi,'("END_DATAGRID_3D_isp_",i1)') isp
c$$$          write(ifi,'("END_BLOCK_DATAGRID_3D")')
c$$$        enddo
c$$$        call fclose(ifi)
c$$$        return
c$$$      else
c$$$        call awrit3('         Writing smooth density to file '//fn//
c$$$     .  '%a : origin at (%i,%i,%i).',' ',80,stdo,lic(5,1)-1,
c$$$     .  lic(6,1)-1,lic(7,1)-1)
c$$$        call dpzero(vecs,3)
c$$$        do  i = 1, 3
c$$$          call daxpy(3,dble(lic(4+i,1)-1)/ngabc(i),plat(1,i),1,vecs,1)
c$$$        enddo
c$$$        call ivset(kv,1,3,2)
c$$$        call shorps(1,plat,kv,vecs,vecs)
c$$$        write(stdo,368) (vecs(i,1), i=1,3)
c$$$  368   format(9x,'In cartesian coordinates, origin is',3f10.6)
c$$$
c$$$        call dpzero(vecs,6)
c$$$        do  k = 1, 2
c$$$          do  i = 1, 3
c$$$            vv = dble((lic(4,k)-1)*lic(i,k))/dble(ngabc(i))
c$$$            call dpadd(vecs(1,k),plat(1,i),1,3,vv)
c$$$          enddo
c$$$          vv = dsqrt(ddot(3,vecs(1,k),1,vecs(1,k),1))
c$$$          sout = ' '
c$$$          call awrit8('%9fv%i: (%i pts) = (%i,%i,%i)p -> '//
c$$$     .    '(%,6;6d,%,6;6d,%,6;6d)',sout,80,0,k,lic(4,k),
c$$$     .    lic(1,k),lic(2,k),lic(3,k),vecs(1,k),vecs(2,k),vecs(3,k))
c$$$          call awrit1('%a l=%,6;6d',sout,-80,-stdo,vv)
c$$$        enddo
c$$$
c$$$C ... Sanity checks
c$$$        lok = .true.
c$$$        lok = lok .and. lic(4,1) .gt. 0
c$$$        lok = lok .and. lic(4,2) .gt. 0
c$$$        if (.not. lok) call fexit2(-1,1,' Exit -1 ioden: number of '//
c$$$     .  'points along axes (%i,%i) are not > 0',lic(4,1),lic(4,2))
c$$$        if (nwk.lt.lic(4,1) .or. nwk.lt.lic(4,2)) call rx('increase nwk')
c$$$
c$$$C --- Copy points to wk from plane ---
c$$$        call icopy(3,lic(5,1),1,kv,1)
c$$$        do  k = 1, lic(4,2)
c$$$          call icopy(3,kv,1,iv,1)
c$$$C         PRINT *, K, KV(1), KV(2), KV(3)
c$$$          do  i = 1, lic(4,1)
c$$$C         PRINT *, ' ',I, IV(1), IV(2), IV(3)
c$$$            do  j = 1, 3
c$$$              jv(j) = mod(iv(j)-1,ngabc(j))+1
c$$$   44         if (jv(j) .le. 0) then
c$$$                jv(j) = jv(j) + ngabc(j)
c$$$                goto 44
c$$$              endif
c$$$            enddo
c$$$            if (jv(1) .le. 0 .or. jv(2) .le. 0 .or. jv(3) .lt. 0) then
c$$$              call rx('bug in ioden')
c$$$            endif
c$$$            do  isp = 1, nsp
c$$$              wk(i,k,isp) = dble(smrho(jv(1),jv(2),jv(3),isp))
c$$$            enddo
c$$$C          if (i .eq. 1 .or. i .eq. 16) then
c$$$C          PRINT *, I,K, JV(1), JV(2), JV(3),WK(I,K)
c$$$C          endif
c$$$            do  j = 1, 3
c$$$              iv(j) = iv(j) + lic(j,1)
c$$$            enddo
c$$$          enddo
c$$$          do  j = 1, 3
c$$$            kv(j) = kv(j) + lic(j,2)
c$$$          enddo
c$$$        enddo
c$$$
c$$$C --- Write the file ---
c$$$        ifi = fopn(fn)
c$$$        rewind ifi
c$$$        do  isp = 1, nsp
c$$$          call ywrm(0,' spin 1',1,ifi,prfmt,wk(1,1,isp),0,nwk,lic(4,1),
c$$$     .    lic(4,2))
c$$$        enddo
c$$$        call fclose(ifi)
c$$$        return
c$$$      endif
c$$$
c$$$   99 continue
c$$$      call rxs('ioden: failed to parse option ... ',sopts(j1:j2))
c$$$
c$$$      end subroutine ioden2
c$$$      end module m_ioden2
c$$$
c$$$      subroutine ioden ( sopts , slat , ssite , sspec , sv_p_orhat 
c$$$     . , smrho )
c$$$
c$$$
c$$$      use m_struc_def  !Cgetarg
c$$$      use m_globalvariables
c$$$      use m_ioden2
c$$$C- File I/O charge density on a uniform mesh in a plane or full 3d mesh
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   sopts :character string specifying plane and other options.
c$$$Ci         :Specifications and options are separated by a delimiter,
c$$$Ci         :which is the first character in sopts.
c$$$Ci         :
c$$$Ci         :if option is g3d then the density is written to
c$$$Ci         :disk on the full 3d grid.   Otherwise:
c$$$Ci         :The density is written to disk for a uniform of mesh of
c$$$Ci         :points in one plane.  This information is specified by three
c$$$Ci         :groups of numbers: the origin, a first direction vector with
c$$$Ci         :its number of points, and a second direction vector with its
c$$$Ci         :number of points.
c$$$Ci         :
c$$$Ci         :At present, these points must be taken from the points on
c$$$Ci         :the smooth mesh, smrho.  In this case, all three groups of
c$$$Ci         :information are sets of integers.  For example,
c$$$Ci         :specify the origin by three numbers:
c$$$Ci         :    o=#1,#2,#3
c$$$Ci         :The point (#1,#2,#3) corresponds to the Cartesian coordinate
c$$$Ci         :   #1/n1 p1 + #2/n2 p2 + #3/n3 p3
c$$$Ci         :where (n1,n2,n3) are the number of divisions in the
c$$$Ci         :mesh along the three lattice vectors (p1,p2,p3).
c$$$Ci         :o=0,0,0 corresponds to the origin.
c$$$Ci         :
c$$$Ci         :Specify the direction vectors by
c$$$Ci         :    l1=#1,#2,#3[,#4]
c$$$Ci         :    l2=#1,#2,#3[,#4]
c$$$Ci         :
c$$$Ci         :l1 and l2 specify the first and second direction vectors,
c$$$Ci         :respectively.  #1,#2,#3 select the
c$$$Ci         :increments in mesh points along each of the three lattice
c$$$Ci         :vectors that define the direction vector.  Thus in Cartesian
c$$$Ci         :coordinates a direction vector is
c$$$Ci         :   #1/n1 p1 + #2/n2 p2 + #3/n3 p3
c$$$Ci         :where as before (n1,n2,n3) are the number of divisions in
c$$$Ci         :the mesh along the three lattice vectors (p1,p2,p3).
c$$$Ci         :The last number (#4) specifies how many points to take
c$$$Ci         :in that direction.
c$$$Ci         :
c$$$Ci         :Other options:
c$$$Ci         :  g3d       output of 3d grid plus headers in xsf format ready for xcrysden
c$$$Ci         :  fn=name   specifies the file name for file I/O
c$$$Ci         :
c$$$Ci         :   core=#   specifies how local rho is to be included
c$$$Ci         :            #=0 include core densities - nuclear charge
c$$$Ci         :            #=1 include core densities
c$$$Ci         :            #=2 (default) exclude core densities
c$$$Ci         :            #=-1 no local densities to be included
c$$$Ci         :            #=-2 true local density, no smoothed part
c$$$Ci         :            #=-3 istl-local sm densities, no true part
c$$$Ci         :
c$$$Ci         :Example: using '~' as the delimiter, and suppose
c$$$Ci         :n1=n2=48 and n3=120, the specification
c$$$Ci         :  ~fn=myrho~o=0,0,60~l1=1,1,0,49~l2=0,0,1,121
c$$$Ci         :writes 'myrho.ext' a mesh (49,121) points.
c$$$Ci         :The origin (first point) lies at (p3/2)
c$$$Ci         :The first vector points along (p1+p2), and has that length;
c$$$Ci         :the second vector points along p3, and has that length.
c$$$Ci   slat  :struct for lattice information; see routine ulat
c$$$Ci     Elts read: plat nabc ng okv vol
c$$$Ci     Stored:    *
c$$$Ci     Passed to: rhgcmp
c$$$Ci   ssite :struct for site-specific information; see routine usite
c$$$Ci     Elts read: *
c$$$Ci     Stored:    *
c$$$Ci     Passed to: rhgcmp
c$$$Ci   sspec :struct for species-specific information; see routine uspec
c$$$Ci     Elts read: *
c$$$Ci     Stored:    *
c$$$Ci     Passed to: rhgcmp
c$$$Ci   orhat :offsets to arrays containing local atom densities
c$$$Ci   smrho :smooth density on uniform mesh
c$$$Co Outputs
c$$$Co   The sum local gaussian densities and smrho is written to disk
c$$$Cl Local variables
c$$$Cl   modrhg:controls what part of core density is added
c$$$Cl         :-1 Add no local densities
c$$$Cl         : 0 include core density - sm nuclear charge
c$$$Cl         : 1 include core density, but no nuclear charge
c$$$Cl         : 2 exclude core density
c$$$Cr Remarks
c$$$Cr    sopts specifies which plane(s) are written to disk
c$$$Cr Bugs
c$$$Cr   Routines create smoothe approximation to density, not true density
c$$$Cu Updates
c$$$Cu   July 7 05 (Walter Lambrecht) New option g3d
c$$$Cu   25 Aug 04 New modes -2, -3
c$$$Cu   24 May 03 Corrections to errors in messages
c$$$Cu   23 Oct 01 Local densities are more accurately represented
c$$$Cu             in G_kL expansion:  k=0..kmax.  Added core= option.
c$$$Cu   25 Apr 01 Simplified the plane specification
c$$$Cu   02 Mar 01 Spin polarized
c$$$Cu   09 Feb 01 Added local gaussian densities to file I/O
c$$$Cu             and extended how a plane may be specified
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$C ... Passed parameters
c$$$      character*(*) sopts
c$$$      type(s_lat)::slat
c$$$      type(s_site)::ssite(*)
c$$$      type(s_spec)::sspec(*)
c$$$
c$$$      type(s_rv1) :: sv_p_orhat
c$$$
c$$$      double complex smrho(*)
c$$$C ... Local parameters
c$$$      integer ngabc(3),n1,n2,n3,k1,k2,k3,kmax
c$$$      integer kkk,lgunit,ng,nglob,nsp,nwk,stdo,modrhg,nbas,i
c$$$ckino Dec.9.2011:          real(8),pointer :: rv_p_opos(:) =>NULL()
c$$$
c$$$ckino Dec.9.2011:         integer,pointer :: iv_p_okv(:) =>NULL()
c$$$
c$$$      complex(8) ,allocatable :: cn_zv(:)
c$$$      complex(8) ,allocatable :: psrho_zv(:)
c$$$ckino not complex but real      complex(8) ,allocatable :: wk_zv(:)
c$$$      real(8) ,allocatable :: wk_rv(:)
c$$$
c$$$      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
c$$$      double precision vol,xx0,xx,dval,plat(3,3),alat
c$$$      integer ib,is,igetss
c$$$C ... Heap
c$$$      integer ::iwdummy
c$$$
c$$$      logical:: l_dummy_isanrg,isanrg
c$$$
c$$$C ... External calls
c$$$ckino del upack,upack2
c$$$      external awrit2,dcopy,fclose,fftz3,fftz30,gvgetf,
c$$$     .gvputf,icopy,ivset,mkils0,mkilst,nwordg,poppr,
c$$$     .pshpr,rhgcmp,rhomom,rx
c$$$
c$$$      real(8),allocatable:: zz(:)
c$$$C ... Unpack and setup
c$$$Changenglob      nsp = nglob('nsp')
c$$$      nsp = globalvariables%nsp
c$$$Changenglob      nbas = nglob('nbas')
c$$$      nbas = globalvariables%nbas
c$$$      stdo = lgunit(1)
c$$$
c$$$      i_copy_size=size(slat%plat)
c$$$      call dcopy(i_copy_size,slat%plat,1,plat,1)
c$$$      i_copy_size=size(slat%nabc)
c$$$      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
c$$$      ng=slat%ng
c$$$ckino Dec.9.2011:         iv_p_okv => slat%iv_p_okv
c$$$
c$$$      vol=slat%vol
c$$$
c$$$ccccccccccccccccccccccccccccc
c$$$c      allocate(zz(nbas))
c$$$c      zz=sspec(1:nbas)%z
c$$$c      print *,' zz1=',zz
c$$$c      print *,' zz2=',sspec(1:nbas)%z
c$$$c      stop 'zzzzzzzzzzzzzzzzzzzzzzzzzzzz'
c$$$ccccccccccccccccccccccccccccc
c$$$
c$$$
c$$$      alat=slat%alat
c$$$ckino Dec.9.2011:         rv_p_opos => slat%rv_p_opos
c$$$
c$$$
c$$$c      call defrr(oz,nbas)
c$$$c      do ib = 1,nbas
c$$$c        is = int(ssite(ib)%spec)
c$$$cc        w(oz-2+ib*2)=sspec(is)%z !but Kino's email on 5June2010.
c$$$c        call dcopy(1,sspec(is)%z,1,w(oz-2+ib*2),1)
c$$$c      enddo
c$$$      call fftz30(n1,n2,n3,k1,k2,k3)
c$$$      kkk = k1*k2*k3
c$$$      kmax = 3
c$$$      modrhg = 2
c$$$c      call ioden2(0,sopts,nsp,plat,alat,nbas,w(opos),w(oz),k1,k2,k3,
c$$$c     .ngabc,w,w,nwk,modrhg)
c$$$Ckino  mode=0 checks format
c$$$#ifdef USEOPTIONAL
c$$$      call ioden2 ( mode=0 , sopts=sopts ,nsp= nsp ,plat= plat ,
c$$$ckino Dec.9.2011:       .   alat=alat , nbas=nbas , pos=rv_p_opos
c$$$ckino Jan.04.2012:       .   alat=alat , nbas=nbas , pos=slat%rv_p_opos 
c$$$     .   alat=alat , nbas=nbas , pos=slat%rv_a_opos 
c$$$     .  , z=sspec ( ssite ( 1:nbas ) %spec ) %z , k1=k1 , k2=k2 , k3=k3
c$$$     .   , ngabc=ngabc 
c$$$     .  , nwk=nwk , modrhg=modrhg )
c$$$#else
c$$$ckino Dec.9.2011:          call ioden2 ( 0 , sopts , nsp , plat , alat , nbas , rv_p_opos
c$$$ckino Jan.04.2012:          call ioden2 ( 0 , sopts , nsp , plat , alat , nbas , slat%rv_p_opos
c$$$      call ioden2 ( 0 , sopts , nsp , plat , alat , nbas , slat%rv_a_opos 
c$$$     .  , sspec ( ssite ( 1:nbas ) %spec ) %z , k1 , k2 , k3 , ngabc 
c$$$     .  , iwdummy , iwdummy , nwk , modrhg )
c$$$#endif
c$$$
c$$$
c$$$Ckino isanrg is logical function,       call isanrg(modrhg,-3,2,'ioden','core option',.true.)
c$$$      l_dummy_isanrg=isanrg(modrhg,-3,2,'ioden','core option',.true.)
c$$$
c$$$C ... Overwrite smrho+, smrho- with smrho, smrho+ - smrho-
c$$$Ckino  smrho is complex(8), dsumdf() is for real(8)
c$$$Ckino  so kkk*2 is necessary to use dsumdf().
c$$$Ckino
c$$$Ckino  return values: smrho(:kkk)=spin(1)+spin(2)
c$$$Ckino  and            smrho(kkk+1:kkk*2)=spin(1)-spin(2)
c$$$      if (nsp .eq. 2) then
c$$$        call dsumdf(kkk*2,1d0,smrho,0,1,smrho(1+kkk),0,1)
c$$$      endif
c$$$
c$$$C ... Put n0(G) into psrho and cn
c$$$      allocate(psrho_zv(kkk*nsp))
c$$$
c$$$      allocate(cn_zv(ng*nsp))
c$$$
c$$$Ckino smrho is complex(8) , so kkk*2 is necessary.
c$$$      call dcopy ( kkk * 2 * nsp , smrho , 1 , psrho_zv , 1 )
c$$$
c$$$      call fftz3 ( psrho_zv , n1 , n2 , n3 , k1 , k2 , k3 , nsp 
c$$$     ., 0 , - 1 )
c$$$
c$$$ckino Dec.9.2011:         call gvgetf ( ng , nsp , iv_p_okv , k1 , k2 , k3 , psrho_zv ,
c$$$ckino Dec.28.2011:         call gvgetf ( ng , nsp , slat%iv_p_okv , k1 , k2 , k3 , psrho_zv ,
c$$$      call gvgetf ( ng , nsp , slat%iv_a_okv , k1 , k2 , k3 , psrho_zv , 
c$$$     . cn_zv )
c$$$Ckino now cn_zv(G,nsp)
c$$$
c$$$
c$$$Ckino xx0 = real and G=0 component
c$$$      xx0 = dval ( cn_zv , 1 )
c$$$
c$$$
c$$$C ... Add sum of local gaussian densities to mesh density
c$$$      print *,' modrhg=',modrhg
c$$$Ckino core=#   specifies how local rho is to be included
c$$$Ckino         :            #=0 include core densities - nuclear charge
c$$$Ckino         :            #=1 include core densities
c$$$Ckino         :            #=2 (default) exclude core densities
c$$$Ckino         :            #=-1 no local densities to be included
c$$$Ckino         :            #=-2 true local density, no smoothed part
c$$$Ckino         :            #=-3 istl-local sm densities, no true part
c$$$      if (modrhg .ge. 0 .or. modrhg .le. -2) then
c$$$        if (modrhg .eq.  0) i = 131
c$$$        if (modrhg .eq.  1) i =  31
c$$$        if (modrhg .eq.  2) i =   1
c$$$        if (modrhg .eq. -2) i =   2
c$$$        if (modrhg .eq. -3) i =   3
c$$$        if ( i .eq. 2 ) call dscal ( ng * nsp * 2 , 0d0 , cn_zv , 
c$$$     .  1 )
c$$$Ckino                      -> cn_zv=0 if i==2
c$$$
c$$$        if ( i .eq. 3 ) call dscal ( ng * nsp * 2 , - 1d0 , cn_zv 
c$$$     .  , 1 )
c$$$Ckino                      -> cn_zv=-cn_zv if i==3
c$$$
c$$$        call rhgcmp ( i , 1 , nbas , ssite , sspec , slat , sv_p_orhat 
c$$$     .   , kmax , ng , cn_zv )
c$$$
c$$$
c$$$        if ( i .eq. 3 ) call dscal ( ng * nsp * 2 , - 1d0 , cn_zv 
c$$$     .  , 1 )
c$$$
c$$$      endif
c$$$      xx = dval ( cn_zv , 1 )
c$$$
c$$$
c$$$C ... FFT (n0 + gaussians) (G) to real-space mesh
c$$$      if (modrhg .ge. 0) then
c$$$ckino Dec.9.2011:           call gvputf ( ng , nsp , iv_p_okv , k1 , k2 , k3 , cn_zv , psrho_zv
c$$$ckino Dec.28.2011:           call gvputf ( ng , nsp , slat%iv_p_okv , k1 , k2 , k3 , cn_zv , psrho_zv
c$$$        call gvputf ( ng , nsp , slat%iv_a_okv , k1 , k2 , k3 , cn_zv , psrho_zv 
c$$$     .   )
c$$$
c$$$
c$$$        call awrit2('%N ioden : local densities + envelope density,'
c$$$     .  //' Qloc=%,6;6d  Q=%,6;6d',' ',80,stdo,(xx-xx0)*vol,xx*vol)
c$$$        if (modrhg .lt. 2) then
c$$$          call awrit1('%9fLocal densities include core'//
c$$$     .    '%?#n==0#+nuclear## contributions.',' ',80,stdo,modrhg)
c$$$        endif
c$$$      elseif (modrhg .eq. -2) then
c$$$        call info2(0,0,0,'%N ioden : local'//
c$$$     .  ' densities (true-smooth terms),  Qloc=%,6;6d',xx*vol,0)
c$$$      elseif (modrhg .eq. -3) then
c$$$        call info2(0,0,0,'%N ioden : smooth (envelope - local)'//
c$$$     .  ' density,  Qs=%,6;6d  Qs-Qs(loc)=%,6;6d',
c$$$     .  xx0*vol,xx*vol)
c$$$      else
c$$$        call info2(0,0,0,'%N ioden : smooth density only'//
c$$$     .  ' (no local densities added) Qsm=%,6;6d',xx0*vol,0)
c$$$      endif
c$$$
c$$$      call fftz3 ( psrho_zv , n1 , n2 , n3 , k1 , k2 , k3 , nsp 
c$$$     ., 0 , 1 )
c$$$Ckino now psrho_zv is real space mesh.
c$$$
c$$$
c$$$      if (allocated(cn_zv)) deallocate(cn_zv)
c$$$
c$$$
c$$$C ... File I/O
c$$$      nwk = 12*max(k1,k2,k3)
c$$$
c$$$ckino not complex but real,      allocate(wk_zv(nwk**2*nsp))
c$$$      allocate(wk_rv(nwk**2*nsp))
c$$$
c$$$Ckino mode=1 write smrho to the file
c$$$ckino Dec.9.2011:         call ioden2 ( 1 , sopts , nsp , plat , alat , nbas , rv_p_opos
c$$$ckino Jan.04.2012:         call ioden2 ( 1 , sopts , nsp , plat , alat , nbas , slat%rv_p_opos
c$$$      call ioden2 ( 1 , sopts , nsp , plat , alat , nbas , slat%rv_a_opos 
c$$$     . , sspec ( ssite ( 1:nbas ) %spec ) %z , k1 , k2 , k3 , ngabc
c$$$ckino not complex but real     . , psrho_zv , wk_zv , nwk , modrhg ) 
c$$$     . , psrho_zv , wk_rv , nwk , modrhg )
c$$$
c$$$
c$$$ckino not complex but real      if (allocated(wk_zv)) deallocate(wk_zv)
c$$$      if (allocated(wk_rv)) deallocate(wk_rv)
c$$$      if (allocated(psrho_zv)) deallocate(psrho_zv)
c$$$
c$$$c      call rlse(oz)
c$$$
c$$$
c$$$C ... Restore smrho+, smrho-
c$$$      if (nsp .eq. 2) then
c$$$        call dsumdf(kkk*2,.5d0,smrho,0,1,smrho(1+kkk),0,1)
c$$$      endif
c$$$      end subroutine ioden
c$$$
c$$$
c$$$
c$$$
