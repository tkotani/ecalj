      subroutine suctrl(prgnam)
C- Creates ctrl file from init file
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam:name of calling program
Co Outputs
Co   an input file is created, default name of `ctrl'
Cr Remarks
Cr
Cu Updates
Cu   08 Jul 07 Reads in effective command-line arguments from cat CMD
Cu   13 Dec 03 Uses plat(3) when scaling translation part of symgrp
Cu   06 Nov 01 Initially created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character prgnam*8
C ... Local parameters
      logical cmdopt,lout
      integer ifin,ifout,fopna,fxst,swt(7),fgtcat,rdtok,catlen,ngmx
      integer mxrec,recl0,nspmx,nbmx,infonl,is,kount
      double precision xx,xv(10)
      parameter (mxrec=200,recl0=128,nspmx=200,nbmx=2000,ngmx=48)
      integer nrmx,nxi0,n0
      parameter (nrmx=1501, nxi0=10, n0=10)
      character*50 outs,ext,errmsg, cxx*10
      character*8 slabl(nspmx),alabl
      double precision z(nspmx),rmt(nspmx),pos(3,nbmx)
      character*(recl0) instr(mxrec),aa
      integer ips(nbmx),nrc(nspmx),lcart(nbmx),lock(nspmx),modep(3)
      integer helppr,infopr,lgunit,stdo,iprint,i,k,nspec,ib,nbas,iz,
     .  a2vec,ng,oistab,oiprm,owk,nggen

C ... For the lattice and space group
      integer havep,quswt
      logical aunits,linfo,lshow,a2bin
C     Following Stuttgart:
C     plat(1) = Lattice vectors actually used
C     plat(2) = Lattice vectors of the most compact unit cell.
C     plat(3) = Lattice vectors of the conventional unit cell.
C               i.e. that with maximum number of right angles.
C               The conventional cell may not be a primitive cell
C               But it is the set which the chemists use and
C               the set for which basis vectors are scaled
C               to convert to cartesian coordinates. 
C              *If plat are specified, plat(3) assigned to plat(1)
C     plat(4) = Lattice vectors of the standard unit cell; see stplat
C     plat(5) = Lattice vectors of the standard unit cell, possibly
C               rotated (roplat.f)
      double precision a,b,c,alat(2),plat(3,3,5),qlat(3,3),
     .  qlatcv(9),rbohr
      double precision dorig(3),alpha,beta,gam,fptol,dglob
      double precision gen(9,ngmx),g(9,ngmx),ag(3,ngmx)
      integer isym(7),ngen,ivec(10),it(10)
      character*30 csym(7), grpnam*10, gens*128, strn*128
C     For determining sphere MT radii 
      double precision omax1(3),omax2(3)
C     for cnvsop
C     character*72 symops(ngmx)
      parameter(rbohr=0.529177d0)

C ... Heap
      integer w(1)
      common /w/ w

C ... data initializations
      data csym/7*'undef'/

C ... Setup
      stdo = lgunit(1)
      gens = ' '
      call stswt('init io,read lunit',stdo,swt)
      if (cmdopt('--pr',4,0,outs) .or. cmdopt('-pr',3,0,outs)) then
        i = 4
        if (cmdopt('-pr',3,0,outs)) i = 3
        i = a2vec(outs,len(outs),i,2,', ',2,2,1,it,ivec)
        call pshpr(ivec)
      endif
      linfo = .false.
      infonl = -1
C     Turn on 'help' mode
      if (cmdopt('--input',7,0,outs)) then
        call stswt('io,help',0,swt)
        linfo = .true.
C     Print out as read in; infonl = number of 'info' newlines
      elseif (cmdopt('--show',6,0,outs)) then
        call stswt('io,rw',0,swt)
        infonl = 0
      endif
      lshow = quswt('io,w',swt) .eq. 1
      if (cmdopt('--iactiv',7,0,outs)) call initqu(.true.)
      infopr = 30
      helppr = 999
      if (linfo) then
        helppr = 0
        infopr = 999
      endif

C ... Open init file for reading and ctrl file for output
      call fextg(ext)
      if (fxst('INIT') .ne. 1) then
        errmsg = 'missing setup file, init'//ext
        goto 999
      endif
      ifin  = fopna('INIT',-1,1)

C ... get command-line arguments in CMD
      call stswt('cat,opt,mxrec',0,swt)
      instr = ' '
      k = fgtcat(ifin,'CMD ','command-line switches can be specified '
     .  //'in this category',swt,.true.,mxrec,recl0,instr)
      if (k .gt. 0) then
        call acmdop(instr(1)(4:),len(instr)*mxrec-4,0)
        if (lshow) call acmdop(instr(1)(4:),len(instr)*mxrec-4,1)
      endif

C --- Lattice vectors and symmetry group ---
      call info(infopr,1,infonl,' ... Reading lattice data',0,0)

      call stswt('cat,reqd,mxrec',0,swt)
      k = fgtcat(ifin,'LATTICE ',
     .  'Specifies the lattice vectors and symmetry group',
     .  swt,.true.,mxrec,recl0,instr)

      call info(helppr,1,0,'   The lattice vectors may be defined '
     .  //'indirectly through%N%5ftoken SPCGRP=, or explicitly using '//
     .  'PLAT=.',0,0)
      call info(helppr,0,0,'%3fIn the former case, the program'//
     .  ' completes the basis%N%5fto make it compatible with the'//
     .  ' space group.',0,0)
      call info(helppr,0,1,'%3fIn the latter case, the basis is '//
     .  'enlarged only if the%N%5fuser additionally '//
     .  'supplies group generators (token GENS=)',0,0)

C ... Look for space group
      call stswt('token,alt,fixlen',0,swt)
C     Read space group symbol
      havep = -1
      k = rdtok('SPCGRP=',instr,
     .  'The space group symbol or space group number (1..230)',
     .  ' ',', ',1,swt,1,1,xx,grpnam)
      if (k .eq. 1) havep = 0
C     Or primitive lattice vectors
      call stswt('token,reqd,fixlen',0,swt)
      if (k .eq. 0 .or. linfo) then
        k = rdtok('PLAT=',instr,
     .  ' %;7d,,Primitive lattice vectors (dimensionless)',
     .  ' ',' ,',4,swt,1,9,plat,cxx)

        havep = 1
        if (.not. lshow)
     .    call info(infopr,0,infonl,' ... read plat',0,0)

        call stswt('token,opt',0,swt)
        k = rdtok('GENS=',instr,
     .  'Generators of the space group',' ',', ',1,swt,1,1,xx,gens)
        if (k .eq. 0) then
          call info(infopr,0,0,' ... no group generators',0,0)
        elseif (.not. lshow) then
          call info(infopr,0,0,' ... read generators',0,0)
        endif
      endif

C ... Check whether units are in AA or a.u.
      call stswt('token,opt',0,swt)
      aunits = (1 .eq. rdtok('UNITS=A',instr,
     .  'Specifies that lattice dimensions are in Angstrom',' ',', ',1,
     .  swt,1,0,xx,cxx))

      call info(helppr,1,0,'   ... If the lattice vectors are '
     .  //'defined through token SPCGRP=,%N%7fone or more of the '//
     .  'following is required.',0,0)
      a = 0
      b = 0
      c = 0
      alpha = 0
      beta = 0
      gam = 0
      if (havep .ne. 1 .or. linfo) then
      call stswt('token,opt,fixlen',0,swt)
      k = rdtok('A=',instr,
     .  '%;7d,,`A'' lattice parameter, in a.u. (Angstrom if UNITS=A)',
     .  ' ',', ',4,swt,1,1,a,cxx)
      k = rdtok('B=',instr,'%;7d,,`B'' lattice parameter',
     .  ' ',', ',4,swt,1,1,b,cxx)
      k = rdtok('C=',instr,'%;7d,,`C'' lattice parameter',
     .  ' ',', ',4,swt,1,1,c,cxx)
      k = rdtok('ALPHA=',instr,'%;7d,,angle(b,c), in degrees',
     .  ' ',', ',4,swt,1,1,alpha,cxx)
      k = rdtok('BETA=',instr,'%;7d,,angle(c,a), in degrees',
     .  ' ',', ',4,swt,1,1,beta,cxx)
      k = rdtok('GAMMA=',instr,'%;7d,,angle(a,b), in degrees',
     .  ' ',', ',4,swt,1,1,gam,cxx)
      endif

      call info(helppr,1,0,'   ... If the lattice vectors are '
     .  //'defined through token PLAT=,%N%7fthe '//
     .  'following is required.',0,0)
      alat(1) = 0
      call stswt('token,o,fixlen',0,swt)
      if (havep .eq. 0) call stswt('token,i',0,swt)
      if (havep .eq. 1) call stswt('token,r',0,swt)
      k = rdtok('ALAT=',instr,
     .  '%;8g,,scaling of lattice vectors PLAT, in a.u.'//
     .  ' (Angstrom if UNITS=A)',' ',', ',4,swt,1,1,alat,cxx)

C ... Construct the lattice vectors from input data
C     Lattice vectors are supplied; construct alternate sets
      if (havep .eq. 1 .and. .not.linfo) then
        alat(2) = alat(1)
        isym(3) = 1
        call mkplat(alat,csym,isym,.true.,plat)
C       If plat specified, use them to define pos = plat * X
        call dcopy(9,plat(1,1,1),1,plat(1,1,3),1)

C     Find space group; plat from isym, a,b,c and alpha.,beta,gamma
      elseif (havep .eq. 0 .and. .not.linfo) then
        call info(infopr,1,0,
     .    ' ... Generating space group from SPCGRP='//grpnam,0,0)
        k = 1
        gens = ' '
        call gengrp(csym,dorig,gens,grpnam,-1,k,isym,ngen,plat,qlat)
        if (isym(7) .eq. 0) call rxs('illegal space group',grpnam)
        call roplat(a,alat,alpha,b,beta,c,gam,isym,
     .    plat(1,1,3),plat(1,1,5))
        call dcopy(9,plat(1,1,5),1,plat(1,1,1),1)
        call dcopy(9,plat(1,1,5),1,plat(1,1,4),1)
        call cpplat(plat(1,1,1),plat(1,1,2))
      endif
C ... Convert to atomic units, if not using them already
      if (aunits) alat(1) = alat(1)/rbohr

C --- SITE category ---
      call info(infopr,1,infonl,' ... Reading site data ...',0,0)

      call stswt('cat,reqd,mxrec',0,swt)
      k = fgtcat(ifin,'SITE ',
     .  'Specifies basis.  At least one site must be specified.',
     .  swt,.true.,mxrec,recl0,instr)

C ... Read in site positions
      nspec = 0
      ib = 1
      catlen = swt(2)
      swt(2) = 1
   40 continue
        swt(1) = swt(2)
        swt(2) = catlen
        call stswt('token,opt,fixlen,mult',0,swt)
        k = rdtok('ATOM=',instr,
     .    'Species label.  There is one token for each site.',
     .    ' ',', ',1,swt,1,1,xx,alabl)
        call stswt('token,nomult',0,swt)
C   ... No further labels found -> terminates infinite loop
        if (k .eq. 0) goto 42

C   ... Add new species, if this label is new
        if (.not.linfo) then
          call tokmat(alabl,slabl,nspec,len(alabl),' ',i,k,.false.)
          if (i .lt. 0) then
             i = nspec
             nspec = nspec+1
             slabl(nspec) = alabl
           endif
          ips(ib) = i+1
        endif
        call info(helppr,0,0,'%5fFollowing each label, '//
     .    'information about the particular site is supplied:',0,0)
        call stswt('token,alt,fixlen',0,swt)
C       Basis vectors as multiples of plat
        lcart(ib) = 0
        k = rdtok('X=',instr,'4;;Site positions, as (fractional) '//
     .    'multiples of the lattice vectors',' ',' ,',4,swt,1,3,
     .    pos(1,ib),cxx)
C       Basis vectors in cartesian coordinates
        call stswt('token,reqd,fixlen',0,swt)
        if (k .eq. 0 .or. linfo) then
          k = rdtok('POS=',instr,'4;;Site positions, cartesian '//
     .      'coordinates',' ',' ,',4,swt,1,3,pos(1,ib),cxx)
          lcart(ib) = 1
        endif

C       print *, ib,slabl(ib),z(ib),rmt(ib)
        ib = ib+1
        if (ib .gt. nbmx) call rxi(
     .    'too many sites. Increase nbmx in suctrl; now',nbmx)
        if (nspec .gt. nspmx) call rxi(
     .    'too many species. Increase nspmx in suctrl; now',nspmx)
        if (.not.linfo) goto 40
C ... End of loop
   42 continue
      nbas = ib-1

      if (.not.linfo) then
      if (nbas .eq. 0) then
        call rx('No sites found')
      endif
      if (iprint() .ge. infopr) call
     .  awrit3('  %i site%?#(n==1)##s#, %i species found',
     .  ' ',80,stdo,nbas,nbas,nspec)

      call info2(infopr+10,0,0,
     .  '%5fScaling and shortening basis vectors',0,0)
C     qlat = (plat+)^-1
      call dinv33(plat(1,1,3),1,qlat,xx)
      do  ib = 1, nbas
C       posc+ = plat posp+
        if (lcart(ib) .eq. 0) then
          call dcopy(3,pos(1,ib),1,xv,1)
          call dmpy(plat(1,1,3),3,1,xv,3,1,pos(1,ib),3,1,3,1,3)
        endif
        if (cmdopt('--noshorten',11,0,outs)) then
          continue
        else
          call shorbz(pos(1,ib),pos(1,ib),plat(1,1,3),qlat)
        endif
      enddo

C     Initialize values z,rmt to mark they haven't been set
      call dvset(z,1,nspec,-1d0)
      call dvset(rmt,1,nspec,0d0)
      endif

C --- SPEC category ---
      call info(infopr,1,infonl,' ... Reading species data ...',0,0)

      call stswt('cat,opt,mxrec',0,swt)
      k = fgtcat(ifin,'SPEC ',
     .  'Specifies additional species information '//
     .  '%N%3f(overriding default values or possible input from '//
     .  'SITE category)'//
     .  '%N%3fLabels should coincide with labels in SITE category%N',
     .  swt,.true.,mxrec,recl0,instr)

C ... Make a list of species
      catlen = swt(2)
      swt(2) = 1
      is = 0
      kount = 0
   30 continue
        swt(1) = swt(2)
        swt(2) = catlen
        call stswt('token,opt,fixlen,mult',0,swt)
        k = rdtok('ATOM=',instr,
     .    'Species label.  Enter one for each species.',
     .    ' ',', ',1,swt,1,1,xx,alabl)
        call stswt('token,nomult',0,swt)
C   ... No further labels found -> terminates infinite loop
        if (k .eq. 0) goto 32
        call info(helppr,0,0,'%5f... Following each label, the'//
     .    ' following may be supplied:',0,0)

C       If this species missing from list of species:
        if (.not.linfo) then
          call tokmat(alabl,slabl,nspec,len(alabl),' ',is,k,.false.)
          if (is .lt. 0) then
            call info(infopr-20,0,0,'%N%5f(warning) ... label `'
     .        //alabl//'%a'' does not match any sites ... '//
     .        'species ignored',0,0)
            goto 30
          else
            is = is+1
          endif
        endif
        kount = kount+1

        k = rdtok('Z=',instr,'%;7d,,4;;Atomic number.  If not '//
     .    'supplied, Z is inferred from the label',' ',', ',4,swt,
     .    is,1,z,cxx)
        k = rdtok('R=',instr,'%;7d,,4;;Sphere radius',
     .    ' ',', ',4,swt,is,1,rmt,cxx)
        if (.not.linfo) goto 30
C ... End of loop
   32 continue
      if (iprint() .ge. infopr) call awrit1
     .  (' species data read for %i species',' ',80,stdo,kount)

      call info(infopr+20,0,0,'     Find missing atomic numbers ...',
     .  0,0)
      if (.not.linfo) then
      do  i = 1, nspec
        if (z(i) .eq. -1) then
          call zslabl(1,slabl(i),iz)
          if (iz .eq. -1) call rxs(
     .    'species label does not correspond to any formula: ',slabl(i))
          z(i) = iz
        endif
      enddo
      endif

      if (linfo) return

C --- Complete basis for compatibility with supplied symmetry group ---
      call info(infopr,1,0,' ... Complete basis for supplied symmetry',
     .  0,0)
      fptol = 1d-5
      if (cmdopt('--fixpos',8,0,strn)) then
        j1 = 8+1
        if (strn(9:13) .eq. ':tol=') then
          j1 = 13
        endif
        if (strn(9:9) .ne. ':' .or.
     .    .not. a2bin(strn,fptol,4,0,' ',j1,len(strn))) fptol = 1d-5
      endif
      call defrr(oistab, ngmx*nbmx)
      if (gens .eq. ' ') gens = 'i*i'
      call gensym(slabl,gens,0,.true.,.false.,fptol,.false.,nbas,nspec,
     .  ngmx,plat,plat(1,1,3),0,xx,pos,ips,nrc,ng,g,ag,ngen,gen,aa,
     .  nggen,isym,w(oistab))
      call rlse(oistab)
C     At this point we can determine the space group name number
C     but to do so we need the following routines taken from
C     Stuttgart code:
C       cnvsop to make a set of generators (Stuttgart conventions)
C       pntgrp fills out isym(4..6), in particular isym(6)
C     Then we can call gengrp to get space group number
C     But this takes too much stuff for something not essential.
C     if (havep .eq. 1) then
C       call cnvsop(.true.,plat,symops,g,ag,ng,1)
C       call pntgrp(csym,isym,ng,symops)
C       call gengrp(csym,dorig,gens,grpnam,1,k,isym,ngen,plat,qlat)
C     endif

      if (.not. cmdopt('--noshorten',11,0,outs)) then
        call info2(infopr+10,1,0,' ... Shortening basis vectors',0,0)
        call dinv33(plat,1,qlat,xx)
        do  ib = 1, nbas
          call shorbz(pos(1,ib),pos(1,ib),plat,qlat)
        enddo
      endif

C --- Sort basis according to species index ---
      call defi(oiprm,nbas)
      call defi(owk,3*nbas)
      call ivheap(1,nbas,ips,w(oiprm),101)
      call ivprm(1,nbas,ips,w(owk),w(oiprm),.true.)
      call dvprm(3,nbas,pos,w(owk),w(oiprm),.true.)
      call rlse(oiprm)

C --- Printout final lattice information ---
      if (iprint() .ge. infopr) then
        call dinv33(plat,1,qlat,xx)
        call dinv33(plat(1,1,3),1,qlatcv,xx)
        call awrit2('%N Lattice vectors : '//
     .    'alat = %,6,6g a.u.  vol=%,6,6g a.u.',
     .    ' ',80,stdo,alat(1),xx*alat(1)**3)
        write (stdo,302) ((plat(i,k,1),i=1,3),(plat(i,k,3),i=1,3),k=1,3)
  302   format(
     .    /t18,'Plat',t47,'Conventional unit cell'/
     .    3(1x,3f11.7,5x,3f11.7/),/1x,
     .    'Basis vectors after sorting and shortening:')
        write(stdo,357)
  357   format(' site spec',8x,'pos (Cartesian coordinates)',9x,
     .  'pos (multiples of plat)')
        do  i = 1, nbas
C         posp = posc (plat+)^-1  and  posp+ = (plat)^-1 posc+
          call dgemm('T','N',3,1,3,1d0,qlat,3,pos(1,i),3,0d0,xv,3)
          print 345, i, slabl(ips(i)), (pos(k,i),k=1,3), (xv(k),k=1,3)
          call dgemm('T','N',3,1,3,1d0,qlatcv,3,pos(1,i),3,0d0,xv,3)
          if (iprint() .ge. infopr+10) print 346, (xv(k),k=1,3)
  345     format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
  346     format(13x,'in terms of conventional unit cell',3f11.6)
        enddo
      endif

      if (cmdopt('--wpos=',7,0,strn)) then
          call iopos(.true.,0,strn(8:),nbas,pos)
C         call fclr(strn(8:),-1)
        endif

      if (linfo) return

C ----------------------- Determine sphere radii ----------------
      call info(10,1,0,' ... Make sphere radii',0,0)
      xx = dglob('lrel',1d0,1)
      xx = dglob('nsp',1d0,1)
      xx = dglob('lxcf',2d0,1)
C     Determine initial guess for sphere radii
      call iinit(lock,nspec)
      do  i = 1, 3
        modep(i) = 2
      enddo
      call makrm0(101,nspec,nbas,alat,plat,pos,slabl,ips,modep,lock,z,
     .  rmt)
C     Scale sphere radii to touching
      call dpzero(omax1,3)
      call dpzero(omax2,3)
      call sclwsr(0,nbas,nbas,nspec,alat,plat,pos,ips,modep,slabl,z,
     .  lock,1d0,3.5d0*0,omax1,omax2,rmt)

C ----------------------- Create the ctrl file ----------------
      call info(10,1,0,' ... Create input file actrl'//ext,0,0)
      ifout = fopna('ACTRL',-1,2)
      rewind ifout
      write(ifout,'(''# Input file autogenerated from init file '',a)')
     .  'init'//ext
      lout = .true.
C     ifout = stdo
      call stswt('io,write lunit',ifout,swt)

C --- VERS ---
      write(ifout,'(''VERS    LMASA-6.0'')')
      write(ifout,
     .  '(''IO      SHOW=f HELP=f WKP=F IACTIV=f VERBOS=31 OUTPUT=*'')')

C --- SYMGRP ---
C     Not needed, but written for informational purposes
      call skpblb(aa,len(aa),k)
      k = k+1
      if (k .gt. 0) then
        i = 0
        call skipbl(aa,k,i)
        write(ifout,'(''#SYMGRP '',a)') aa(i+1:k)
      else
        write(ifout,'(''SYMGRP  find'')')
      endif

C --- STRUC (v6.xx format) ---
      call stswt('cat,reqd,mxrec',0,swt)
      k = fgtcat(ifin,'STRUC ',
     .  'Specifies the lattice vectors basis size',
     .  swt,.false.,mxrec,recl0,instr)

      call stswt('lunit',0,swt)
      k = rdtok('NSPEC=',instr,
     .  'number of species',' ',' ',2,swt,1,1,nspec,cxx)
      k = rdtok('NBAS=',instr,
     .  '%i,,number of species',' ',' ',2,swt,1,1,nbas,cxx)
      call stswt('lunit',ifout,swt)
      k = rdtok('NL=',instr,
     .  'default value of 1+l for basis and augmentation',' ',' ',
     .  2,swt,1,1,3,cxx)
      call stswt('lunit',0,swt)
      k = rdtok('ALAT=',instr,
     .  '%;8g,,scaling of lattice vectors, in a.u.',
     .  ' ',', ',4,swt,1,1,alat,cxx)
C      k = rdtok('DALAT=',instr,
C     .  '%;8g,,is added to alat after input is read',
C     .  ' ',', ',4,swt,1,1,0d0,cxx)
      call stswt('lunit',ifout,swt)
      k = rdtok('PLAT=',instr,
     .  ' %;7d,,Primitive lattice vectors (dimensionless)',
     .  ' ',' ,',4,swt,1,9,plat,cxx)

C --- SPEC category ---
      call stswt('cat,reqd,mxrec',0,swt)
      k = fgtcat(ifin,'SPEC ',
     .  'Specifies additional species information '//
     .  '%N%3f(overriding default values or possible input from '//
     .  'SITE category)'//
     .  '%N%3fLabels should coincide with labels in SITE category%N',
     .  swt,.true.,mxrec,recl0,instr)

      call stswt('token,opt,fixlen,mult',0,swt)
      catlen = swt(2)
      swt(2) = 1
      is = 0
      kount = 0
  130 continue
        if (lout) then
          if (kount .ge. nspec) goto 132
          is = kount
          alabl = slabl(is+1)
        else
          swt(1) = swt(2)
          swt(2) = catlen
        endif
        call stswt('lunit',0,swt)
        k = rdtok('ATOM=',instr,
     .    '%16p,,Species label.  Enter one for each species.',
     .    ' ',', ',1,swt,1,1,xx,alabl)

C       No further labels found -> terminates infinite loop
        if (k .eq. 0) goto 132
        call info(helppr,0,0,'%5f... Following each label, the'//
     .    ' following may be supplied:',0,0)

C       If this species missing from species list, do
        if (.not.linfo) then
          call tokmat(alabl,slabl,nspec,len(alabl),' ',is,k,.false.)
          if (is .lt. 0) then
            call info(infopr-20,0,0,'%N%5f(warning) ... label `'
     .        //alabl//'%a'' does not match any sites ... '//
     .        'species ignored',0,0)
            goto 130
          else
            is = is+1
          endif
        endif
        kount = kount+1

C       Tokens for this species
        k = rdtok('Z=',instr,'%;4,0D%b %b,,4;;Atomic number.  If not '//
     .    'supplied, Z is inferred from the label',' ',', ',4,swt,
     .    is,1,z,cxx)
        call stswt('lunit',ifout,swt)
        k = rdtok('R=',instr,'%;9,6D,,4;;Sphere radius',
     .    ' ',', ',4,swt,is,1,rmt,cxx)
        if (.not.linfo) goto 130
C ... End of loop
  132 continue

C --- SITE category ---
      call stswt('cat,reqd,mxrec',0,swt)
      k = fgtcat(ifin,'SITE ',
     .  'Specifies basis.  '//
     .  'At least one must be specified.%N%3f'//
     .  'Labels must coincide with specifications'//
     .  ' in SPEC category.%N',
     .  swt,.true.,mxrec,recl0,instr)

      ib = 1
      catlen = swt(2)
      swt(2) = 1
      is = 0
  140 continue
        if (lout) then
          if (ib .gt. nbas) goto 142
          is = ips(ib)
          alabl = slabl(is)
        else
          swt(1) = swt(2)
          swt(2) = catlen
        endif
        call stswt('lunit',0,swt)
        call stswt('token,opt,fixlen,mult',0,swt)
        k = rdtok('ATOM=',instr,'%16p,,Species label.',
     .    ' ',', ',1,swt,1,1,xx,alabl)
C   ... No further labels found -> terminates infinite loop
        if (k .eq. 0) goto 142
C   ... Add new species, if this label is new
        if (.not.linfo .and. .not.lout) then
          call tokmat(alabl,slabl,nspec,len(alabl),' ',i,k,.false.)
          if (i .lt. 0) then
             i = nspec
             nspec = nspec+1
             slabl(nspec) = alabl
           endif
          ips(ib) = i+1
        endif
        is = ips(ib)

        call info(helppr,0,0,'%5f... Following each label, the'//
     .    ' following must be supplied:',0,0)
        call stswt('token,reqd,fixlen',0,swt)
        call stswt('lunit',ifout,swt)
        k = rdtok('POS=',instr,'%,;11,7D,,4;;Site positions, Cartesian '
     .    //'coordinates',' ',' ,',4,swt,1,3,pos(1,ib),cxx)

        ib = ib+1
        if (ib .gt. nbmx) call rxi(
     .    'too many sites. Increase nbmx in suctrl; now',nbmx)
        if (nspec .gt. nspmx) call rxi(
     .    'too many species. Increase nspmx in suctrl; now',nspmx)
        if (.not.linfo) goto 140
C ... End of loop
  142 continue
      nbas = ib-1

      return
C --- Error exit ---
  999 continue
      outs = prgnam // '%a : ' // errmsg
      call fexit(1,9,outs,0)

      end

