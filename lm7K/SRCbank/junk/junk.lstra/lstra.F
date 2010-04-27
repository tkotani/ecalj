      subroutine lstra(sname,indxs,offe,noffe)
C- Returns a list of indices corresponding to names in struc
C ----------------------------------------------------------------
Ci Inputs
Ci   sname   struc type: must be one of those listed in styps, below,
Ci           followed by a ist of elements sought in struc;
Ci           see Remarks, below
Co Outputs
Co   indxs   index to which entry in styps sname corresponds
Co   offe    list of indices in struc corresponding to names
Co   noffe   number of entries in offe
Ci   ... entry lstrs
Ci   offe(1) index in struc for which name is sought
Co   sname   element name associated with index offe(1)
Co   noffe   number of elements in structure
Ci   ... entry lstrp
Ci   offi     offset in work array to this structure
Cr Remarks
Cr   lstra makes indices associated with names to entries in
Cr   structures it knows about.  For lstra to work, the indices
Cr   elements in a given structure are sorted alphabetically by
Cr   name.  Internal tables within each structure keep track
Cr   of the offset within the structure.
Cr   NB: the structures named here should be consistent with
Cr       those in spack and spackv, below.
Cr
Cr  *To create a new structure (e.g. 'ugw')
Cr     1.  Add name to 'styps'.  Must add new nlist,list,lists.
Cr         list must be alphabetically ordered,
Cr         so you must shift all indices list[k,...] to
Cr         list[k+1,...] in the data statements fo list[k]
Cr         and similarly parameter statements in nlist[k]
Cr         Update nlists; add entry in data statement for 'offs'
Cr     2.  Create a new subroutine, e.g. 'ugw'  The number of
Cr         elements defined there must match nlist[k]
Cr     3.  Create a new data statement for list[k]
Cr     4.  Add a new entry in spack.  Again, shift the
Cr         computed 'goto' statements for k+1...
Cr         If the structure has a species index, first argument
Cr         in 'call pack(...') should be 10+lpack, otherwise
Cr         just lpack.
Cr     5.  In routine spackv, make sure that the indices
Cr         corresponding to those structures with species
Cr         (eg uspec, usite) correspond with those in spack,
Cr         and if your new struc has species, add lines there.
Cu Updates
Cu   29 Jul 08 New ctrl->pfloat
Cu   05 Jul 08 New ham->(lmxax,ndham,npwmin,npwpad,pwemax,pwemin,pwmode,oveps)
Cu   27 Jun 08 New str->rmaxg,ivl spec->ehvl
Cu   07 Dec 07 Add sham->qpoff
Cu   17 Nov 07 LDA+U parameters sham->nlibu,lmaxu,udiag; remove sctrl->nlibu
Cu   08 Jun 07 New stb and sspec parameters for JEK's additions to tbe
Cu   01 May 07 New pot->opdel
Cu   27 Mar 07 New ctrl->nlibu
Cu   27 Jun 06 New sham->pmax
Cu   21 Mar 06 New ctrl->defm
Cu   01 Feb 06 New spec->rcfa
Cu   21 Dec 05 Increase spec->coreq to length 2
Cu    1 Aug 05 Parameters for ES finder
Cu   27 Apr 05 (WRL) parameters for LDA+U
Cu   12 Sep 04 Increment size of bz->semsh
Cu   25 Jun 04 Add uspec->eh3,rs3,vmtz
Cu   19 May 04 Add sstr->loka
Cu   13 May 04 Add spec->(colxbs,lxi,radxbs,rcut,rint)
Cu             Add sham->alfsi,dabc
Cu             Add ssite->mpole,dpole
Cu             add str->nbisi,nalf,ncupl,ndust,adec,wx
Cu    6 Apr 04 Add pot->nrhos and pot->obxc
Cu    5 Mar 04 Enlarge sham->eterms to size 20
Cu    2 Mar 04 Added rpad to LAT structure
Cu    6 Jan 04 Added rsstol to HAM structure
Cu   24 Jun 03 Added ecuts to GW structure
Cu   24 Jun 03 (Chantis) relativistic arrays in spot
Cu   21 May 03 Added sham->sigp
Cu   09 Mar 03 Enlarged sham->eterms
Cu   28 Feb 03 Potential function arrays for fully rel case (upot)
Cu   27 Jan 03 parms for magnetic field in (ham, site, spot) strux
Cu   22 Jul 02 More parameters to sham for r.s. sigma
Cu   22 Jul 02 Added parameters to sham and sgw for r.s. sigma
Cu   10 Apr 02 Redimensionsed sspec to accomodate larger lmax
Cu   12 Oct 00 new structure sgw.
Cu    2 Feb 00 igets returns nint(strux entries with double cast)
Cu             if specified range, lgors uses first entry
C ----------------------------------------------------------------
C     implicit none
      character*(*) sname
      integer indxs,offe(1),noffe
C Local variables
      integer nmax,ntyps,
     .  nlist1,nlist2,nlist3,nlist4,nlist5,nlist6,nlist7,nlist8,nlist9,
     .  nlista,nlistb,nlistc,nlistd,nliste,nlistf,nlistg
C     ntyps = number of structures; nmax = max length of structure
      parameter (ntyps=16,nmax=74)
      parameter (nlist1=32,nlist2=30,nlist3=56,nlist4=19,nlist5=56,
     .  nlist6=43,nlist7=23,nlist8=11,nlist9=13,nlista=11,nlistb=58,
     .  nlistc=27,nlistd=74,nliste=30,nlistf=9,nlistg=7)
      integer nlists(ntyps),offs(ntyps),i,j,offi
      character*6 lists(nmax,ntyps),styps(ntyps),
     .  list1(nlist1),list2(nlist2),list3(nlist3),list4(nlist4),
     .  list5(nlist5),list6(nlist6),list7(nlist7),list8(nlist8),
     .  list9(nlist9),lista(nlista),listb(nlistb),listc(nlistc),
     .  listd(nlistd),liste(nliste),listf(nlistf),listg(nlistg)
      equivalence
     .  (lists(1,1),list1),(lists(1,2),list2),(lists(1,3),list3),
     .  (lists(1,4),list4),(lists(1,5),list5),(lists(1,6),list6),
     .  (lists(1,7),list7),(lists(1,8),list8),(lists(1,9),list9),
     .  (lists(1,10),lista),(lists(1,11),listb),(lists(1,12),listc),
     .  (lists(1,13),listd),(lists(1,14),liste),(lists(1,15),listf),
     .  (lists(1,16),listg)
      save styps,nlists,offs
C     styps = structure names
      data styps /'array','bz','ctrl','gw','ham','lat','mix','move',
     .  'optic','ordn','pot','site','spec','str','strn','tb'/
      data offs  /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data list1 /'ncl','nclasp','nofgl','nofgr','npadl','npadr','npl!',
     .  'oclabl','ocllst','oclp','oclssl','ogroup','ogrp2','ohave',
     .  'oics','oifrlx','oipc','oipcp','oips','omxcst','onrc','onrcp',
     .  'opgfsl','opgfvl','opgord','opgplp','opos','ormax','osordn',
     .  'ostni','ozos','size'/
      data list2 /'def','dosw','ef','efmax','fsmom','lcond','lio',
     .  'lmet','lmull','lopt','lshft','n','ndos','nevmx','nkabc',
     .  'nkp','ntet','odos','oidtet','oipq','opdos','oqp','ostar',
     .  'owtkp','range','semsh','size','stnr','w','zval'/
      data list3 /'defm','elin','lasa','lbas','lcd','lcgf','ldos','lfp',
     .  'lfrce','lgen3','lham','lmet','lncol','loptc','lordn','lpgf',
     .  'lqp','lrel','lrs','lscr','lstonr','lstr','lsx','ltb',
     .  'lves','lxcf','maxit','mdprm','modep','nbas','nbasp','nclass',
     .  'nesabc','nitmv','nl','nmap','npl','nsite','nspec','nspin',
     .  'nvario','omax1','omax2','pfloat','quit','rmaxes','rmines',
     .  'sclwsr','sdmod','sdprm','sdxsi','size','smalit','tol',
     .  'wsrmax','zbak'/
      data list4 /'delre','deltaw','deltax','ecuts','gcutb','gcutx',
     .  'gsmear','lgw','lmxpb','mksig','nband','nime','nkabc','pb1',
     .  'pb2','pbtol','qoffp','rdsig','size'/
      data list5 /'alfsi','amgm','bandw','dabc','ehf','ehk','elind',
     .  'eterms','hord','kmto','lasa','ldham','lgen3','lham','lmaxu',
     .  'lmxax','lncol','lsig','ltb','lxcf','nbf','ndham','ndhrs',
     .  'ndofH','neula','nkaph','nlibu','nmto','npwmin','npwpad',
     .  'nqsig','obdots','oeula','ohrs','oiaxs','oindxo','olmxa',
     .  'omagf','ontabs','ooffH','oqsig','oveps','pmax','pmin','pwemax',
     .  'pwemin','pwmode','qpoff','qss','rsrnge','rsstol','seref',
     .  'sigp','size','thrpv','udiag'/
      data list6 /'alat','as','avw','awald','dist','gam','gmax','ldist',
     .  'nabc','ng','nkd','nkdmx','nkq','nkqmx','npgrp','nsgrp','oag',
     .  'obgv','ocg','ocy','odlv','ogv','oidxcg','oips0','oistab',
     .  'ojcg','okv','opos','oqlv','osymgr','plat','plat0','plat2',
     .  'plate','platl','platr','qlat','rpad','size','slat','tol',
     .  'tolft','vol'/
      data list7 /'b','bl','bv','elind','fn','kill','lxpot','mmix',
     .  'mode','model','n','nitu','nmix','nsave','r','rms1','rms2',
     .  'size','tj','tolu','umix','w','wc'/
      data list8 /'ct','gyro','kt','modt','nmodt','prmint','size',
     .  'tnow','ts','tsequ','tstot'/
      data list9 /'axes','cll','cln','cls','clsite','esciss','mode',
     .  'nchi2','ne','ocrng','size','unrng','window'/
      data lista /'efre','mode','ncl','ndofh','oclp','oclssl','oiaxg',
     .  'omapgv','ooffch','rmaxg','size'/
      data listb /'bfield','nlma','nlml','nrhos','oaamom','obxc',
     .  'odddpf','oddpf','oddpfr','odel','odpf','odpfr','ofes','ogma',
     .  'ogmar','ogrrme','ohab','oivso','omad','oorhat','opalp','opapg',
     .  'opdel','opf','opfnc','opfr','opmpol','opnu','opp','oppi',
     .  'oppn','opprel','opti','oqc','oqmom','oqnu','oqpp','oqt',
     .  'orhos','orhrmx','osab','osgw','osmpot','osmrho','osop',
     .  'osoptc','osrout','otau','ovab','ovdif','oves','ovintr',
     .  'ovrmax','ovshf','size','vconst','vmtz','vmtz0'/
      data listc /'bfield','clabel','class','cli','delta','dpole',
     .  'eula','force','mpole','ndelta','norb','offh','orho','ov0',
     .  'ov1','pl','plv','pnu','pos','pos0','pz','relax','sid','size',
     .  'spec','vel','vshft'/
      data listd /'a','alpha','beta','chfa','colxbs','coreh','coreq',
     .  'ctail','dv','eh3','ehvl','enu','eref','etail','etf','exi',
     .  'group','grp2','hcr','idmod','idu','idxdn','iq1','ivso','jh',
     .  'kmxh','kmxt','kmxv','lfoca','lmxa','lmxb','lmxf','lmxl',
     .  'lmxpb','lxi','mass','mxcst','name','naug','ngcut','norb',
     .  'norp','nr','ntorb','nxi','orbp','orhoc','p','pb1','pb2','pz',
     .  'q','qc','qpol','radxbs','rcfa','rcut','rfoca','rg','rham',
     .  'rint','rmt','rs3','rsma','rsmfa','rsmv','size','stc','stni',
     .  'stnm','uh','vmtz','vso','z'/
      data liste /'adec','amode','drwats','iinv','ivl','kaps','lequiv',
     .  'lmaxw','loka','lshow','mxnbr','nalf','nbisi','ncupl','ndust',
     .  'nitab','nkaps','nttab','oadot','oalp','oiax','ontab','os',
     .  'osdot','rfit','rmax','rmaxg','size','skmsh','wztcf'/
      data listf /'amix','gemb','gfopt','jobid','map','mix','mmham',
     .  'sxopt','symg'/
      data listg /'alam','alsc','ebfit','fmode','nbfit','rmfit','size'/
      data nlists /nlist1,nlist2,nlist3,nlist4,nlist5,nlist6,
     .             nlist7,nlist8,nlist9,nlista,nlistb,nlistc,
     .             nlistd,nliste,nlistf,nlistg/

C      print 333, listd
C  333 format(a)
C      stop

      call word(sname,1,i,j)
C     Find index indxs to structure with ascii repsn snam(i:j)
      call strxls(styps,ntyps,sname(i:j),indxs,noffe)
      if (noffe.lt.0) call rxs('lstra:  unknown structure: ',sname)

      if (j+1 .lt. len(sname)) then
        call strxls(lists(1,indxs),nlists(indxs),sname(j+2:),offe,noffe)
      else
        noffe = 0
      endif
      return

C --- Entry lstrs ---
      entry lstrs(sname,indxs,offe,noffe)
      sname =  lists(offe(1),indxs)
      noffe = nlists(indxs)
      return

C --- Entry lstrp ---
      entry lstrp(sname,offi)

      call strxls(styps,ntyps,sname,i,j)
      offs(i) = offi

      end
      double precision function dgets(sname,struc)
C- Returns one double precision entry in a structure, referred by name
C     implicit none
      character*(*) sname
      integer mask(5),range(2,5),nmax
      parameter (nmax=30)
      double precision struc(1),
     .  x1(nmax),x2(nmax),x3(nmax),x4(nmax),x5(nmax)
      character*40 trunam,alias,switch

      call salias(sname,trunam,alias,switch,range,mask)
      call spack(0,trunam,struc,x1,x2,x3,x4,x5)
      if (range(1,1) .eq. -1) then
        dgets = x1(1)
      else
        if (range(1,1) .gt. nmax) call rx('dgets: increase nmax')
        dgets = x1(range(1,1))
      endif

      end
      double precision function dgetss(sname,is,struc)
C- Returns one double precision entry in a structure, referred by name
C     implicit none
      integer is
      character*(*) sname
      double precision struc(1),x1(20),x2(20),x3(20),x4(20),x5(20)

      call spacki(0,sname,struc,is,x1,x2,x3,x4,x5)
      dgetss = x1(1)

      end
      subroutine scalss(sname,is1,is2,struc,val)
C- Scales one entry in a structure by val
C     implicit none
      double precision val
      integer is1,is2
      character*(*) sname
      double precision struc(1)
      integer off(5),cast(5),nelt(5),i,is

      do  10  is = is1, is2
        call spacki(2,sname,struc,is,off,cast,nelt,i,i)
        if (cast(1) .ne. 4) call rx('scalss not ready')
        off(1) = off(1) + nint(struc(1))*max((is-1),0)
        call dscal(nelt(1),val,struc(off(1)),1)
   10 continue

      end
      integer function igets(sname,struc)
C- Returns one integer entry in a structure, referred by name
C     implicit none
      character*(*) sname
      double precision struc(1)
      integer nmax
      parameter (nmax=30)
      double precision x1(nmax),x2(nmax),x3(nmax),x4(nmax),x5(nmax)
      integer mask(5),range(2,5),bitand,offi(5),cast(5),nelti(5)
      integer ix1,ival,k
      character*40 trunam,alias,switch

      call salias(sname,trunam,alias,switch,range,mask)
C     get the internal representation for the cast
      call spack(2,trunam,struc,offi,cast,nelti,0,0)
      call spack(0,trunam,struc,x1,x2,x3,x4,x5)
      if (range(1,1) .gt. nmax) call rx('igets: increase nmax')
      k = 1
      if (range(1,1) .gt. 1) k = range(1,1)
      if (cast(1) .eq. 2) then
        ix1 = ival(x1,k)
      else
        ix1 = nint(x1(k))
      endif
      igets = ix1
      if (mask(1) .gt. 0) igets = bitand(mask,ix1)

      end
      integer function igetss(sname,is,struc)
C- Returns one integer entry in one species of a structure, referred by name
C     implicit none
      character*(*) sname
      integer ix1(30),is
      double precision struc(1),x2(30),x3(30),x4(30),x5(30)

      call spacki(0,sname,struc,is,ix1,x2,x3,x4,x5)
      igetss = ix1(1)

      end
      logical function lgors(sname,struc)
C- Logical OR of bits in one (masked) entry in a structure
C     implicit none
      character*(*) sname
      integer mask(5),bitand,range(2,5),ival,off(5),cast(5),nelt(5),i,k
      double precision struc(1),x1(20),x2(20),x3(20),x4(20),x5(20)
      character*40 trunam,alias,switch

      call salias(sname,trunam,alias,switch,range,mask)
      call spack(0,trunam,struc,x1,x2,x3,x4,x5)
      if (mask(1) .le. 0) then
        lgors = ival(x1,1) .ne. 0
      else
        call spack(2,trunam,struc,off,cast,nelt,i,i)
        if (range(1,1) .le. 1) range(1,1) = 1
        k = range(1,1)-1
        if (cast(1) .eq. 2) then
          lgors = bitand(mask,ival(x1,1+k)) .ne. 0
        else
          lgors = bitand(mask,nint(x1(1+k))) .ne. 0
        endif
      endif
      end
      logical function lgands(sname,struc)
C- Logical AND of all bits corresponding to mask in a structure entry
C     implicit none
      character*(*) sname
      integer ix1(30),mask(5),bitand,range(2,5)
      double precision struc(1),x2(30),x3(30),x4(30),x5(30)
      character*40 trunam,alias,switch

      call salias(sname,trunam,alias,switch,range,mask)
      call spack(0,trunam,struc,ix1,x2,x3,x4,x5)
      lgands = bitand(mask,ix1) .eq. mask(1)
      end

      subroutine lsets(sname,struc,lval,mask)
C- Set a logical bits to a structure entry
C ----------------------------------------------------------------
Ci Inputs
Ci   lval is a logical T or F, of length of sname entry
Ci   mask should be an integer multiple of 2.
Ci        Only the lowest bit of mask is used.
Co Outputs
Co  struc element corresponding to label 'name' is modified.
Co        The mask bit of that entry is set to lval.
C ----------------------------------------------------------------
C     implicit none
      logical lval(1)
      integer mask
      character*(*) sname
      integer ix1(30),lmask,bitlow,bitor,bitand,
     .  off(5),cast(5),nelt(5),i,ires
      double precision struc(1),x2(30),x3(30),x4(30),x5(30)

      if (mask .le. 0) call rxi('lsets: mask must be >0 but given',mask)
      call spack(0,sname,struc,ix1,x2,x3,x4,x5)
      call spacki(2,sname,struc,0,off,cast,nelt,i,i)
      lmask = bitlow(mask)
      if (nelt(1) .ne. 1) call rx('lsets needs be checked here')
      do  10  i = 1, nelt(1)
        if (lval(i)) then
          struc(off(1)) = bitor(int(struc(off(1))),lmask)
        else
          ires = int(struc(off(1)))
          struc(off(1)) = ires - bitand(ires,lmask)
        endif
        off(1) = off(1)+1
   10 continue

      end

C      subroutine lsets(sname,struc,lval,mask)
CC- Set a logical bits to integer entry in a structure, referred by name
CC ----------------------------------------------------------------
CCi Inputs
CCi   lval is a logical T or F
CCi   mask should be an integer multiple of 2.
CCi        Only the lowest bit of mask is used.
CCo Outputs
CCo  struc element corresponding to label 'name' is modified.
CCo        The mask bit of that entry is set to lval.
CC ----------------------------------------------------------------
C      implicit none
C      logical lval
C      integer mask
C      character*(*) sname
C      integer ix1(30),maski,nmask,t2n
C      double precision struc(1),x2(30),x3(30),x4(30),x5(30)
C
C      if (mask .le. 0) call rxi('lsets: mask must be >0 but given',mask)
C      call spack(0,sname,struc,ix1,x2,x3,x4,x5)
C      maski = mask
C      nmask = ix1(1)
C      t2n = 1
CC ... Find lowest nonzero bit of mask, corresponding ix1
C   10 continue
C      if (mod(maski,2) .eq. 0) then
C        t2n = t2n*2
C        maski = maski/2
C        ix1(1) = ix1(1)/2
C        goto 10
C      endif
C      if (lval) then
C        nmask = nmask + t2n*(1-mod(ix1(1),2))
C      else
C        nmask = nmask + t2n*(0-mod(ix1(1),2))
C      endif
C      call spack(1,sname,struc,nmask,x2,x3,x4,x5)
C
C      end

      subroutine spacki(lpack,sname,struc,is,x1,x2,x3,x4,x5)
C- Pack/unpack integer or double precision elements in a struc by name
C  if is=0, calls spack(x1..x5); else calls spack(is,x1..x4)
C     implicit none
      integer lpack,is
      character*(*) sname
      double precision struc(1),x1(1),x2(1),x3(1),x4(1),x5(1)

C ... Pick up offsets in structure
      if (is .eq. 0) call spack(lpack,sname,struc,x1,x2,x3,x4,x5)
      if (is .gt. 0) call spack(lpack,sname,struc,is,x1,x2,x3,x4)
      end

      subroutine spack(lpack,sname,struc,x1,x2,x3,x4,x5)
C- Pack/unpack integer or double precision elements in a struc by name
C ----------------------------------------------------------------
Ci Inputs
Ci   lpack   0 for unpack, 1 for pack, 2 to return parameters
Ci   sname   'strucnam member member ...' : both strucnam and its
Ci           members must belong to those listed in lstra.f.
Ci   x1      (integer) species index, for strucs containing species
Cio  struc   structure containing elements
Cio  x1..x5  elements to pack/unpack (x2..x5 if x1 used for input is)
Cio          depending on value of lpack; see packi
Cr Remarks
Cr   Structure types must be consistent with those in lstra, above
C ----------------------------------------------------------------
C     implicit none
      integer lpack
      character*(*) sname
      double precision struc(1),x1(1),x2(1),x3(1),x4(1),x5(1)
C Local variables
      integer is1,is2,i1,i2,j1,j2,offe(20),noffe,indxs
      external ubz,uarray,uctrl,ugw,uham,ulat,umix,umove,uoptic,uordn,
     .  upot,usite,uspec,ustr,ustrn,utb

C ... Call packi with appropriate function name and list
      call lstra(sname,indxs,offe,noffe)

      goto (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16) indxs
      call rxi('spack: unrecognized index',indxs)

    1 call packi(lpack,uarray,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

    2 call packi(lpack,ubz,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

    3 call packi(lpack,uctrl,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

    4 call packi(lpack,ugw,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

    5 call packi(lpack,uham,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

    6 call packi(lpack,ulat,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

    7 call packi(lpack,umix,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

    8 call packi(lpack,umove,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

    9 call packi(lpack,uoptic,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

   10 call packi(lpack,uordn,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

   11 call packi(lpack,upot,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

   12 call packi(10+lpack,usite,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

   13 call packi(10+lpack,uspec,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

   14 call packi(lpack,ustr,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

   15 call packi(lpack,ustrn,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

   16 call packi(lpack,utb,struc,offe,noffe,x1,x2,x3,x4,x5)
      goto 20

   20 continue

      if (noffe .lt. 0) then
        call word(sname,1,i1,i2)
        call word(sname,-noffe+1,j1,j2)
        call rxs4('spack: unrecognized element "',sname(i1:i2),'->',
     .    sname(j1:j2),'"')
      endif
      return

      entry spackv(lpack,sname,struc,is1,is2,x1)
C- Pack/unpack a vector of one element type from a struc containing spec
C ----------------------------------------------------------------
Ci Inputs
Ci   lpack  1s digit
Ci           0 unpacks, 1 pack
Ci         10s digit
Ci           0 x1 is same for all species (packing only)
Ci           1 x1 is different for each species
Ci           2 x1 and its size is different for each species
Ci             (not implemented)
Ci is1,is2
Ci    x1    packed into is1..is2.
Cr Remarks
Cr   Structure types must be consistent with those in lstra, above
C ----------------------------------------------------------------

C ... Call packi with appropriate function name and list
      call lstra(sname,indxs,offe,noffe)
      if (noffe .lt. 0) then
        call word(sname,1,i1,i2)
        call word(sname,-noffe-1,j1,j2)
        call rxs4('spack: unrecognized element "',sname(i1:i2),
     .    '->',sname(j1:j2),'"')
      endif

      if (indxs .eq. 12) then
        call packv(lpack,usite,struc,offe,is1,is2,x1)
      elseif (indxs .eq. 13) then
        call packv(lpack,uspec,struc,offe,is1,is2,x1)
      endif

      end
      subroutine upack(sname,struc,x1,x2,x3,x4,x5)
C- Unpack integer or double precision elements in a struc by name
C     implicit none
      character*(*) sname
      double precision struc(1),x1,x2,x3,x4,x5
      call spack(0,sname,struc,x1,x2,x3,x4,x5)
      return

      entry upack2(sname,struc,x1,x2)
      call spack(0,sname,struc,x1,x2,x2,x2,x2)
      return

      entry upack3(sname,struc,x1,x2,x3)
      call spack(0,sname,struc,x1,x2,x3,x3,x3)
      return

      entry upack4(sname,struc,x1,x2,x3,x4)
      call spack(0,sname,struc,x1,x2,x3,x4,x4)
      return

      entry pack2(sname,struc,x1,x2)
      call spack(1,sname,struc,x1,x2,x2,x2,x2)
      return

      entry pack3(sname,struc,x1,x2,x3)
      call spack(1,sname,struc,x1,x2,x3,x3,x3)
      return

      entry pack4(sname,struc,x1,x2,x3,x4)
      call spack(1,sname,struc,x1,x2,x3,x4,x4)
      return

      entry pack5(sname,struc,x1,x2,x3,x4,x5)
      call spack(1,sname,struc,x1,x2,x3,x4,x5)
      return

C ... This returns lists of offset, cast and number of elts for list
      entry upackp(sname,struc,x1,x2,x3,x4)
      call spack(2,sname,struc,x1,x2,x3,x4,x4)
      return

      end
      subroutine upack1(sname,struc,x1)
C- Sets/retrives a subblock of one entry in a structure
C     implicit none
      character*(*) sname
      double precision struc(1),x1(1)
      integer mask(5),range(2,5),offi,casti,nelti
      character*40 trunam,alias,switch
      logical lupack

      lupack = .true.
      goto 10

      entry pack1(sname,struc,x1)
      lupack = .false.

   10 continue
      call salias(sname,trunam,alias,switch,range,mask)
      call spack(2,trunam,struc,offi,casti,nelti,0,0)
      if (range(1,1) .le. 1) range(1,1) = 1
      if (range(2,1) .lt. 0) range(2,1) = nelti
      offi = offi+range(1,1)-1
      nelti = min(nelti,range(2,1))-range(1,1)+1
c     print *, range(1,1),offi,nelti
      if (casti .eq. 2) then
        if (lupack) call discop(struc(offi),x1,nelti,1,1,0)
        if (.not. lupack) call idscop(x1,struc(offi),nelti,1,1)
      else
        if (lupack) call dcopy(nelti,struc(offi),1,x1,1)
        if (.not. lupack) call dcopy(nelti,x1,1,struc(offi),1)
      endif

      end

      subroutine upacks(sname,is1,is2)
C- Unpacks one entry in string structure
C     implicit none
      character*(*) sname
      integer is1,is2,i,offs,j,xx,leni,casti

      call lstra(sname,i,offs,j)
      if (j .eq. -1)
     .  call rxs('upacks: unrecognized element "',sname)
      call ustrn(xx,offs,1,is1,casti,leni)
      is2 = is1 + leni-1

      end
      subroutine spacks(lpack,sname,struc,name,is1,is2)
C- Pack/unpack string element in a struc, element referenced by name
C ----------------------------------------------------------------
Ci   lpack   0 for unpack, 1 for pack, 2 to return parameters
Ci   is1..is2 range of species is for which to pack
C             For structure without species, set both to 0
C ----------------------------------------------------------------
C     implicit none
      integer is,lpack
      character*(*) sname,name(1)
      double precision struc(1)
C Local variables
      integer offi(5),nelti(5),casti(5),is1,is2,i

      call spacki(2,sname,struc,is1,offi,casti,nelti,i,i)
      offi(1) = offi(1) + nint(struc(1))*max((is1-1),0)

      i = 1
      do  10  is = is1, is2
        if (lpack .eq. 1) then
          call s8tor8(name(i),struc(offi(1)))
        elseif (lpack .eq. 0) then
          call r8tos8(struc(offi(1)),name(i))
        endif
        i = i+1
        offi(1) = offi(1) + nint(struc(1))
   10 continue
      end
      subroutine shstru(sname,struc,is1,is2)
C- Debugging printout for struc, species is1..isw
C ----------------------------------------------------------------
Ci   sname    structure name, which lstra must know about.
Ci   struc    array containing structure information
Ci   is1,is2  range of species to print
Ci            is2 should be 0 if struc has no species
C ----------------------------------------------------------------
C     implicit none
      integer is1,is2
      character *(*) sname
      double precision struc(1)
C Local variables
      integer ifi,lgunit,is,js1,js2,ielt,indxs,nelt,
     .  offi,casti,nelti,ival,ks
      character s*120, eltnam*6, strn*8, name*8
      double precision xx(100)
      logical lspec
*     integer ix(100)
      integer NULLI
      parameter (NULLI=-99999)

      lspec =  (is2 .ne. 0)
      if (lspec) then
        js1 = is1
        js2 = is2
      else
        js1 = 1
        js2 = 1
      endif
      ifi = lgunit(1)
      print '(/'' shstru: parameters for structure "'',a,''"'')', sname

C ... Pull out structure type and number of elements in this structure
      call lstra(sname,indxs,ielt,nelt)
      call lstrs(eltnam,indxs,1,nelt)

      do  10  is =  js1, js2

C   ... Header info for this species, if structure has a species
        ks = 0
        if (is2 .gt. 0) then
          s = sname // ' name'
          call lstra(s,indxs,ielt,nelti)
          name = ' '
          if (nelti .gt. 0) call spacks(0,s,struc,name,is,is)
          call awrit2('%u%N ... species %i%?#n>0#, "'//name//'%a"#',s,
     .      len(s),ifi,is,nelti)
          ks = is
        endif

        do  20  ielt = 1, nelt
          call lstrs(eltnam,indxs,ielt,nelt)
          s = sname // ' ' // eltnam
          if (lspec) then
            call upackp(s,struc,is,offi,casti,nelti)
            call upack2(s,struc,is,xx)
          else
            call upackp(s,struc,offi,casti,nelti,xx)
            call upack1(s,struc,xx)
          endif
*         print *, eltnam, offi, casti, nelti

          if (nelti .gt. 0) then
            if (casti .eq. 2) then
*         print *, ival(xx,1)
              if (ival(xx,1) .ne. NULLI) call
     .          awrit2(' '//eltnam//'%u%8p%n:2i ',s,len(s),ifi,nelti,xx)
            elseif (casti .eq. 4) then
              if (xx(1) .ne. dble(NULLI)) call awrit2(' '//eltnam//
     .          '%u%8p%n:2,6;6d ',s,len(s),ifi,nelti,xx)
            elseif (casti .eq. 1) then
              call spacks(0,s,struc,strn,ks,ks)
              if (xx(1) .ne. dble(NULLI))
     .          call awrit0(' '//eltnam//'%10p'//strn,s,len(s),ifi)
            endif
          endif
   20   continue
   10 continue

      end

      subroutine ubz(sbz,offe,nli,off,cast,nelt)
C- Returns offsets, casts and number for Brillouin zone integration
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe  list of parameters to extract
Ci         offe(1) = -1 => initialize sbz
Ci   nli   number of entries in offe
Ci         offe(1) = -1 => use nli as initialization value
Co   off   for each parameter in offe, returns offset to sbz
Co   cast  for each parameter in offe, returns cast of entry
Co   nelt  for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   def      Uncertainty in Fermi level
Cr    3    2   dosw     Energy window over which DOS accumulated
Cr    5    3   ef       Fermi level
Cr    6    4   efmax    eigenvalues above efmax are not needed
Cr    7    5   fsmom    fixed-spin moment (fixed spin moment method)
Cr    8    6   lcond    conductivity, and directions
Cr   12    7   lio      1  read qp
Cr                      2  write qp
Cr                      4  write dos
Cr                      8  density matrix
Cr                      16 save evcs
Cr                      32 save dos weights
Cr                      NB: 4,32 out of date ; see ldos
Cr   13    8   lmet     0 insulator   1 3-point scheme  2 save evcs
Cr                      3 two passes  4 wt from prior iter
Cr                      (numbers, not powers of 2)
Cr   14    9   lmull    not used
Cr   15   10   lopt     2 BZMAP option
Cr   16   11   lshft    shift mesh in each of three directions
Cr   19   12   n        Polynomial order for Methfessel-Paxton sampling
Cr   20   13   ndos     No. DOS points (sampling integration, and lmdos)
Cr   21   14   nevmx    eigenvectors above nevmx are not needed
Cr   22   15   nkabc    number of divisions in qp mesh
Cr   25   16   nkp      Number of qp
Cr   26   17   ntet     number of tetrahedra
Cr   27   18   odos     offset to total dos
Cr   28   19   oidtet   offset to tetrahedron corners
Cr   29   20   oipq     offset to ipq array made in bzmesh
Cr   30   21   opdos    offset to partial dos
Cr   31   22   oqp      offset to qp array
Cr   32   23   ostar    offset to array containing info on star of k
Cr   33   24   owtkp    offset to array containing qp weights
Cr   34   25   range    number of FWHM for sampling integration
Cr   35   26   semsh    nz, modec, emin,emax, ecc,eps delta
Cr   45   27   size     size of this structure
Cr   46   28   stnr     parameters for Stoner model (not used)
Cr   47   29   w        Line broadening for sampling integration
Cr   48   30   zval     No. electrons to accumulate in BZ integration
C ----------------------------------------------------------------
C     implicit none
      integer nli,offe(nli),off(nli),cast(nli),nelt(nli)
C Local variables
      integer i,size,nlists,psize
      parameter (size=48, nlists=30, psize=45)
      double precision sbz(size)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the offe
      data ilists /2,3,5,6,7,8,12,13,14,15,16,19,20,21,22,25,26,27,28,
     .  29,30,31,32,33,34,35,45,46,47,48,49/
C ... Casts of each of the entries in offe
      data casts/4,4,4,4,4,4,18*2,4,4,-1,4,4,4/

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(nli),0,sbz,1)
        sbz(1) = size
        sbz(psize) = size
        return
      endif

      do  10  i = 1, nli
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue
      end
      subroutine uarray(sarray,lsti,nlsti,off,cast,nelt)
C- Returns offsets, casts and number for global array offsets and data
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci     is  site index
Ci   lsti  list of parameters to extract
Ci          lsti(1) = -1 => initialize sarray
Ci   nlsti  number of entries in lsti
Co   off    for each parameter in lsti, returns offset to sarray
Co   cast   for each parameter in lsti, returns cast of entry
Co   nelt   for each parameter in lsti, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr    off lsti  name    purpose
Cr    2    1   ncl     (O-N) number of clusters
Cr    3    2   nclasp  (pgf) number of classes including PL -1,npl
Cr    4    3   nofgl   (pgf) number of left  PL for off-diagonal GF
Cr    5    4   nofgr   (pgf) number of right PL for off-diagonal GF
Cr    6    5   npadl   (pgf) number of basis atoms in PL -2,-1,0
Cr    7    6   npadr   (pgf) number of basis atoms in PL npl-1,npl,npl+1
Cr    8    7   npl     (pgf) number of principal layers
Cr    9    8   oclabl  offset to class label array (r8tos8 format)
Cr   10    9   ocllst  (O-N) list of sites belonging to each cluster
Cr   11   10   oclp    (O-N) cluster-related dimensioning parameters
Cr   12   11   oclssl  (O-N) site ib belongs to cluster clssl(ib)
Cr   13   12   ogroup
Cr   14   13   ogrp2
Cr   15   14   ohave   Switches containing what info available, by class
Cr                     1 P,Q  2 pp  4 sop  8 vintra  16 pmpol  32 gradrm
Cr   16   15   oics    Index to which species each class belongs
Cr   17   16   oifrlx  (dynamics)
Cr   18   17   oipc    index to which class each site belongs
Cr   19   18   oipcp   index to which class each site belongs
Cr   20   19   oips    index to which species each site belongs
Cr   21   20   omxcst  mixing constraint
Cr   22   21   onrc    number of sites in each class
Cr   23   22   onrcp   number of sites in each class
Cr   24   23   opgfsl
Cr   25   24   opgfvl
Cr   26   25   opgord
Cr   27   26   opgplp
Cr   28   27   opos
Cr   29   28   ormax   sphere radius, by class
Cr   30   29   osordn  structure for order-N information
Cr   31   30   ostni   Stoner
Cr   32   31   ozos    Stoner DOS
Cr   33   32   size    size of this structure
C ----------------------------------------------------------------
C     implicit none
      integer nlsti,lsti(nlsti),off(nlsti),cast(nlsti),nelt(nlsti)
C Local variables
      integer i,size,nlists,psize
      parameter (size=33, nlists=32, psize=33)
      double precision sarray(size)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the lsti
      data ilists /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     .  22,23,24,25,26,27,28,29,30,31,32,33,34/
C ... Casts of each of the entries in lsti
      data casts/ 31*2,-1/

      if (lsti(1) .eq. -1) then
        call dcopy(size,dble(nlsti),0,sarray,1)
        sarray(1)  = size
        sarray(psize) = size
        return
      endif

      do  10  i = 1, nlsti
        il = lsti(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue

      end
      subroutine uctrl(sctrl,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for control parameters
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe  list of parameters to extract
Ci         offe(1) = -1 => initialize sctrl
Ci   noffe number of entries in offe
Ci         offe(1) = -1 => use noffe as initialization value
Co   off   for each parameter in offe, returns offset to sctrl
Co   cast  for each parameter in offe, returns cast of entry
Co   nelt  for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   defm    Lines for deformation when relaxing plat
Cr    8    2   elin    ccor linearization energy, 2C hamiltonian
Cr    9    3   lasa    1 Make V from P,Q  2 Make pp  4 ccor  8 free atm
Cr                     16 map  32 nonspherical mpol moms 64 MT corr
Cr                     128 interpretation of sphere Q2; see newrho.f
Cr                     256 how ASA Q1,Q2 are accumulated; see makwts.f
Cr                     512 (spin pol) alpha repsn = (gamma(1) + gamma(nsp))/2
Cr   10    4   lbas    1 Hamiltonian has no screening transformation
Cr                     2 Hamiltonian is nfp style
Cr                     16 freeze phi,phidot for all species
Cr   11    5   lcd     1 freeze core
Cr                     2 non-self-consistent Harris
Cr                     4 represent full potential density on a uniform mesh
Cr                     8 represent full potential density via TCF
Cr                    16 unused
Cr                    32 unused
Cr                    64 (molecules) XC potential by FFT
Cr   12    6   lcgf    Green's function
Cr                    >0 Green's function calculation
Cr   13    7   ldos    1 make dos
Cr                     2 generate weights for partial dos
Cr                     4 generate weights for m-decompos'n of pdos
Cr                     8 BZ map
Cr   14    8   lfp     switches for full-potential
Cr                     1 program uses full potential
Cr   15    9   lfrce   How forces are calculated
Cr                     0 do not calculate forces
Cr                     1 shift in FA density
Cr                     2 shift in core+nucleus
Cr                    10 added for shift to be screened
Cr   16   10   lgen3   switches for third generation LMTO
Cr                     1 3rd generation LMTO
Cr                     2 make structure constants by Ewald summation
Cr   17   11   lham    1 (ASA) 2-center
Cr                     1 (molecules) two-panel
Cr                     2 (ASA) 2-c + pert. corr
Cr                     4 (ASA) auto-down-fold
Cr                     8 (ASA) change rep interactively
Cr                    16 (ASA) suppress d hybridization
Cr                    32 (ASA) preserve ortho. evecs
Cr                    64 (ASA) save evecs to disk
Cr                   128 (ASA) gamma-rep
Cr                   256       use true spherical harmonics
Cr                   512       not used
Cr                  1024       not used
Cr   18   12   lmet    1 metal     2 tetrahedron
Cr                     4 (GF) V-shift1 is zero  8 (GF) V-shift2 is zero
Cr   19   13   lncol   1 noncollinear magnetism
Cr                     2 spin spirals
Cr                     4 spin-orbit coupling
Cr                     8 External magnetic field
Cr                    16 mag. forces
Cr                    32 spin-orbit coupling, LzSz only
Cr                    64 spin-orbit coupling, LzSz + (L.S-LzSz) pert
Cr   20   14   loptc   1 generate Im(eps) (optics package)
Cr                     2 generate Im(eps) w/ on-the-fly sampling integration
Cr                    -1 generate joint DOS
Cr                     Add 10 for SHG
Cr   21   15   lordn   1 Embedded GF  4 Vanderbilt  8 Embedded cluster
Cr   22   16   lpgf    layer GF.  First entry:
Cr                     1: Generate diagonal GF and output density
Cr                     2: Find Fermi left- and right-bulk Ef
Cr                     3: Trace out band structure for spec'd energy mesh
Cr                        In this case make Im(z)=0.
Cr                     4: Ditto for right bulk.
Cr                     5: Calculate current through structure
Cr                        Second entry:
Cr                     0  Make GF using standard approach
Cr                     1  Use LU decomposition
Cr   24   17   lqp     1 do not add inversion 2 inverse iteration
Cr   25   18   lrel    1 scalar relativistic
Cr                     2 fully relativistic
Cr   26   19   lrs     switches concerning restart mode.
Cr                     1 Read from restart file
Cr                     2 Read from restart file, ascii mode
Cr                     4 Read from restart file, invoke smshft
Cr                     8 Write new density to restart file
Cr                    16 Write new density to restart file, ascii format
Cr                    32 read site positions from input file
Cr                    64 read starting fermi level from input file
Cr                   128 read starting pnu level from input file
Cr                   256 rotate local density after reading
Cr   27   20   lscr    0 do nothing
Cr                     1 Make P0(0)
Cr                     2 Screen output q and ves
Cr                     3 Screen output ves only
Cr                     4 Use model response to screen output q
Cr                       Add 1 to combine mode 1 with another mode
Cr                       Add 10*k to compute intra-site contribution to
Cr                       vbare each kth iteration
Cr                       Add 100*k to compute response function only
Cr                       each kth iteration
Cr   28   21   lstonr  second digit for graphical output
Cr   31   22   lstr    Note: no longer used; see str->lshow,str->lequiv
Cr                     1 print strux   2 find equiv strux
Cr   32   23   lsx     1 Calculate screened exchange sigma
Cr                    10's digit nonzero make vintra to include in sigma
Cr   33   24   ltb     1 overlap        2 crystal-field     4 ovlp+CF
Cr                     8 add ebarLL    16 forces           32 fij
Cr                    64 not used     128 pressure        256 evdisc
Cr                   512 pair pot    1024 TrH & local E  2048 local rho
Cr                  2^12 Hubbard U   2^13 No Madelung    2^14 wgt avg U
Cr                  2^15 L>0 estat   2^16 disc read incr 2^17 gamma-pt
Cr   34   25   lves    1 take ves as input
Cr   35   26   lxcf    parameter defining XC functional
Cr                     1s digit:
Cr                     1 for Ceperly-Alder
Cr                     2 for Barth-Hedin (ASW fit)
Cr                     3 for PW91
Cr                     4 for PBE
Cr                     10s digit for GGAs
Cr                     0 for LSDA
Cr                     1 for LMH
Cr                     2 for PW91
Cr                     3 for PBE
Cr                     4 for PBE with Becke exchange
Cr   36   27   maxit   max. no.  iterations in self-consistency cycle
Cr   37   28   mdprm   arg 1: 1 new dynamics  2  restart dynamics
Cr                            4 relax with conjugate gradients
Cr                            5 relax with variable metric
Cr                            6 relax with Broyden
Cr                     arg 2: statics: switch
Cr                            1 read hessian matrix
Cr                            dynamics:
Cr                              number of iterations between printouts.
Cr                     arg 3: (stat) relaxation x-tolerance
Cr                            (dyn)  temperature
Cr                     arg 4: (stat) relaxation g-tolerance
Cr                            (dyn)  time step
Cr                     arg 5: (stat) step length
Cr                            (dyn)  relaxation time
Cr                     arg 6: (stat) Remove hessian after this many steps
Cr                            (dyn)  --
Cr   43   29   modep.. which dimensions are periodic
Cr   46   30   nbas    size of basis
Cr   47   31   nbasp   size of padded basis (layer geometry)
Cr   48   32   nclass  size of class
Cr   49   33   nesabc  mesh for determining empty spheres
Cr   52   34   nitmv   max number of mol-dynamics iterations
Cr   53   35   nl      1 + maximum lmxa
Cr   54   36   nmap    number of maps (ASA)
Cr   55   37   npl     number of principal layers (PGF)
Cr   56   38   nsite   number of sites
Cr   57   39   nspec   number of species
Cr   58   40   nspin   number of spins
Cr   59   41   nvario  number of variables to output
Cr   60   42   omax1   sphere overlap constraints, type 1
Cr   63   43   omax2   sphere overlap constraints, type 2
Cr   66   44   pfloat  stop after:  1 show  2 atom  4 bands
Cr   67   45   quit    stop after:  1 show  2 atom  4 bands
Cr   68   46   rmaxes  upper limit to ES radius when finding new empty spheres
Cr   69   47   rmines  lower limit to ES radius when finding new empty spheres
Cr   70   48   sclwsr  scale wsr until reaching this fractional vol
Cr                     10s digit used for assymetric treatment of ES:
Cr                       0 ES and other sites are treated symmetrically
Cr                       1 all sites with z>0 are resized first; then
Cr                         all sites are resized.
Cr                       2 all sites with z>0 are resized first; then
Cr                         the ES sites only are resized
Cr   71   49   sdmod   spin dynamics mode
Cr   72   50   sdprm   spin dynamics parameters
Cr   77   51   sdxsi   Bulgac and Kusnezov global deamons
Cr   81   52   size    size of this structure
Cr   82   53   smalit  parameters for small iterations
Cr   84   54   tol     1 q- tolerance for self-consistency
Cr                     2 e- tolerance for self-consistency
Cr   87   55   wsrmax  constraint on size of largest WS sphere
Cr   88   56   zbak    background charge and MT correction parameters
Cr                     zbak(1) = uniform background charge included
Cr                               in electrostatics but not in xc pot.
Cr                     zbak(2) = charge used to make MT correction
Cr                               (ASA only)
C ----------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,size,nlists,psize
      parameter (size=90, nlists=56, psize=81)
      double precision sctrl(size)
      integer ilists(nlists+1),casts(nlists),il
C      data ilists /2,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,
C     .  26,27,28,31,32,33,34,35,36,37,43,46,47,48,49,52,53,54,55,56,57,
C     .  58,59,60,63,66,67,68,69,70,71,76,80,81,83,86,87,89/
      data ilists /2,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,25,
     .  26,27,28,31,32,33,34,35,36,37,43,46,47,48,49,52,53,54,55,56,57,
     .  58,59,60,63,66,67,68,69,70,71,72,77,81,82,84,87,88,90/

C ... Pointers to locations in structure for each of the offe
C ... Casts of each of the entries in offe
      data casts/ 4,4,25*2,4,13*2,4,4,2,2,4,4,4,2,4,4,-1,2,4,4,4/

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(noffe),0,sctrl,1)
        sctrl(1)  = size
        sctrl(psize) = size
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue
      end
      subroutine ugw(sgw,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for GW and dielectric response
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe  list of parameters to extract
Ci         offe(1) = -1 => initialize sgw
Ci   noffe number of entries in offe
Ci         offe(1) = -1 => use noffe as initialization value
Co   off   for each parameter in offe, returns offset to sgw
Co   cast  for each parameter in offe, returns cast of entry
Co   nelt  for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   delre   spacing on real axis for energy contour
Cr                     delre(1) = dw(0)
Cr                     delre(2) = omegac (coff to parabolic term)
Cr    4    2   deltaw  width in finite diff for sigma energy derivative
Cr    5    3   deltax  Im. E (broadening delta) in computation of x0
Cr    7    4   ecuts   energy cutoff for calculating sigma_nn'
Cr    7    5   gcutb   G-vector cutoff for basis
Cr    8    6   gcutx   G-vector cutoff for response function
Cr    9    7   gsmear  smearing in pole of Green's function
Cr   10    8   lgw     specifies mode for GW computation
Cr                     1s digit
Cr                     0 standard mode
Cr                     1 Faleev's fast energy integration
Cr   11    9   lmxpb   default l-cutoff for product basis
Cr   12   10   mksig   Mode for generating self-energy sigma matrix
Cr                     0 do not create sigma matrix
Cr                     1 Make sigma_nn'(E_f)
Cr                     2 Make sigma_nn'(E_n/2+E_n'/2)
Cr                     3 Make sigma_nn'(E_n)/2 + sigma_nn'(E_n')/2
Cr                     lsig>0 implies 1s digit lgw>0
Cr   13   11   nband   number of bands in perturbation calculation
Cr   14   12   nime    # energies along Im axis for energy contour
Cr   15   13   nkabc   the number of k-points for GW calculation
Cr   18   14   pb1     default product basis 1
Cr   19   15   pb2     default product basis 2
Cr   20   16   pbtol   product-basis overlap cutoff
Cr   21   17   qoffp   qp offset parameter
Cr   22   18   rdsig   parameters for how to read,interpolate SE matrix
Cr   23   19   size    size of this structure
C -----------------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,size,nlists,psize
      parameter (size=23, nlists=19, psize=23)
      double precision sgw(size)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the offe
      data ilists /2,4,5,6,7,8,9,10,11,12,13,14,15,18,19,20,21,22,23,24/
C ... Casts of each of the entries in offe
      data casts /4,4,4,4,4,4,4,2,2,2,2,2,2,1,1,4,4,2,-1/

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(noffe),0,sgw,1)
        sgw(1)  = size
        sgw(psize) = size
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue
      end
      subroutine uham(sham,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for hamiltonian parameters
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe  list of parameters to extract
Ci         offe(1) = -1 => initialize sham
Ci   noffe number of entries in offe
Ci         offe(1) = -1 => use noffe as initialization value
Co   off   for each parameter in offe, returns offset to sham
Co   cast  for each parameter in offe, returns cast of entry
Co   nelt  for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   alfsi  stability factor alfsi*(S+)*S (molecules)
Cr    3    2   amgm   magnetization
Cr    4    3   bandw  Maximum size of off-diagonal, packed storage
Cr    5    4   dabc   Spacing for real-space mesh (molecules)
Cr    8    5   ehf    Harris-Foulkes energy
Cr    9    6   ehk    Hohnberg-Kohn energy
Cr   10    7   elind  Lindhard screening parameter
Cr   11    8   eterms terms making up the total energy.   For FP:
Cr                     1  =  ehf
Cr                     2  =  eks
Cr                     3  =  utot
Cr                     4  =  rhoves int (rhoval x Ves) (true-sm)
Cr                     5  =  cpnves int (core+nucleus x Ves) (true-sm)
Cr                     6  =  rhoeps int (rho * exc)
Cr                     7  =  rhomu  int (rho * vxc)
Cr                     8  =  sumec
Cr                     9  =  sumtc
Cr                     10 =  xcore  rhoc * total potential
Cr                     11 =  valvef int (rhov * vef) (true-sm)
Cr                     12 =  sumt0  sum of foca core energies
Cr                     13 =  sumev  band structure energy
Cr                     13 =  dq1
Cr                     14 =  dq2
Cr                     15 =  amom   system magnetic moment
Cr                     16 =  sumev  sum of single-particle evals
Cr                     17 =  rinvxt input density * external pot
Cr                     18 =  rouvxt output density * external pot
Co                     19 =  rhosig trace of self-energy over occ states
Cr   31    9   hord   order of polynomial approximation to ham.
Cr                    In 2nd gen LMTO, relevant only in GF context.
Cr   32   10   kmto   envelope kinetic energies making up chord LMTO
Cr   38   11   lasa   not used; see ctrl lasa
Cr   39   12   ldham  vector describing hamiltonian dimensions:
Cr                    1: ldim   = dimension of lmto basis
Cr                    2: lidim  = ldim + size of downfolding block
Cr                    3: lidhim = lidim + size of higher block
Cr                    4: nspc   = number of coupled spins
Cr                              = 1 unless noncollinear magnetism
Cr                    5: ldimc  = ldim * nspc
Cr                    6: lidimc = lidim * nspc
Cr                    7: lihdimc = lihdim * nspc
Cr                    8: nspx   = number of separate spin channels
Cr                              = nsp if nspc is not 2; else 1
Cr   55   13   lgen3   shouldn't be used; just a copy of ctrl->lgen3
Cr   56   14   lham    1 (ASA) 2-center
Cr                     1 (molecules) two-panel
Cr                     2 (ASA) 2-c + pert. corr
Cr                     4 (ASA) auto-down-fold
Cr                     8 (ASA) change rep interactively
Cr                    16 (ASA) suppress d hybridization
Cr                    32 (ASA) preserve ortho. evecs
Cr                    64 (ASA) save evecs to disk
Cr                   128 (ASA) gamma-rep
Cr                   256       use true spherical harmonics
Cr                   512       not used
Cr                  1024       not used
Cr                       NB: ham->lham and ctrl->lham should be duplicates
Cr   57   15   lmaxu   dimensioning parameter for LDA+U potential
Cr   58   16   lmxax   largest augmentation lmax in basis
Cr   59   17   lncol   1 noncollinear magnetism
Cr                     2 spin spirals
Cr                     4 spin-orbit coupling
Cr                     8 External magnetic field
Cr                    16 mag. forces
Cr                    32 spin-orbit coupling, LzSz only
Cr                    64 spin-orbit coupling, LzSz + (L.S-LzSz) pert
Cr                    NB: ham->lncol and ctrl->lncol should be duplicates
Cr   60   18   lsig   parameters concerning self-energy
Cr                     1 read sigma, Assume r.s. sigma is real
Cr                     2 read sigma
Cr   61   19   ltb    not used; see ctrl->ltb
Cr   62   20   lxcf   parameter defining local XC functional.
Cr                    Not used: see ctrl->lxcf
Cr   63   21   nbf    number of channels for magnetic field
Cr   64   22   ndham  Largest dimension of lmto+PW basis
Cr   65   23   ndhrs  dimension of site block for r.s. h or sigma
Cr   66   24   ndofH  leading dimension to ooffH
Cr   67   25   neula  number of channels for euler angles
Cr                    0 -> No euler angles defined
Cr   68   26   nkaph  number of repetitions of one l-quant.n. in basis
Cr   69   27   nlibu  Number of LDA+U blocks
Cr   70   28   nmto   Polynomial order to NMTO
Cr   71   29   npwmin Estimate for minimum number of PWs to be calculated
Cr                    (PW basis is q-dependent; max size not known a priori)
Cr   72   30   npwpad Padding to be added to estimated max basis dimension
Cr                    and subtracted from min basis dimension
Cr                    (PW basis is q-dependent; max size not known a priori)
Cr   73   31   nqsig  Number of k-points for which self-energy known
Cr                    (Used for interpolating between qpoints)
Cr   74   32   obdots mag. field in local coordinates, dotted with
Cr                    Paulie matrices (pointer in w array)
Cr   75   33   oeula  Euler angles (pointer in w array)
Cr   76   34   ohrs   r.s. h or sigma (pointer in w array)
Cr   77   35   oiaxs  neighbor table (pointer in w array)
Cr   78   36   oindxo orbital permutation table (pointer in w array)
Cr   79   37   olmxa  vector of augmentation l-cutoffs (pointer in w array)
Cr   80   38   omagf  magnetic field (pointer in w array)
Cr   81   39   ontabs table of number-of-pairs per site (pointer in w array)
Cr   82   40   ooffH  table of offsets to hamiltonian matrix (pointer in w array)
Cr                    (pointer in w array)
Cr                    (Used for interpolating between qpoints)
Cr   83   41   oqsig  list of qp at which sigma can be computed
Cr   84   42   oveps  When diagonalizing hamiltonian, discard part of hibert space
Cr                    corresponding to evals of overlap < oveps
Cr   85   43   pmax   global minimum allowed values for pnu
Cr   95   44   pmin   global minimum allowed values for pnu
Cr  105   45   pwemax High Energy cutoff for PW part of basis
Cr  106   46   pwemin Low Energy cutoff for PW part of basis
Cr  107   47   pwmode Controls PW part of basis
Cr                    0 => no PW part of basis
Cr                    1 => include PWs in basis
Cr                    2 => include only PWs in basis
Cr  108   48   qpoff  qp offset when generating transformation of self-energy
Cr  111   49   qss    spin spiral
Cr  115   50   rsrnge cutoff in connecting vector length for r.s. sigma
Cr  116   51   rsstol tolerance in Bloch sum error for r.s. sigma
Cr  117   52   seref  Sum of reference energies
Cr  118   53   sigp   parameters for approximating self-energy sigma
Cr                    arg 1: mode : specifies how to set its diagonal part
Cr                           for states above the high-energy cutoff
Cr                           0 constrain sigii to be > asig+bsig*e
Cr                           1 constrain sigii to be = asig+bsig*e
Cr                           2 constrain sigii to be > asig and < bsig
Cr                           3 constraint same as mode 1.
Cr                             Mode 3 differs in that the least-squares fit to
Cr                             sigii (for informational purposes only, to help
Cr                             estimate asig and bsig) is done for states between
Cr                             efit and nmax or emax
Cr                    arg 2: nmin : sigma for states 1..nmin are approximated by sigii
Cr                    arg 3: emin : (used only if nmin<0)
Cr                                : sigma for levels e<emin are approximated by sigii
Cr                    arg 4: nmax : sigma for levels i>nmax are approximated by
Cr                                  sigii AND constrained according to mode
Cr                    arg 5: emax : (used only if nmax<=0)
Cr                                : sigma for levels e<emax are approximated by
Cr                                  sigii AND constrained according to mode
Cr                    arg 6: asig : constraint used to approximate
Cr                                  sigii = asig + E * bsig  or
Cr                                  asig < sigii < bsig
Cr                    arg 7: bsig : constraint used to approximate
Cr                                  sigii = asig + E * bsig  or
Cr                                  asig < sigii < bsig
Cr                    arg 8: efit : (mode 3) energy minimium (not used?)
Cr                                  for fitting asig and bsig
Cr                    arg 9: ? lwrite:T, write sigii to file (not used?)
Cr                    10 not used now
Cr  128   54   size   size of this structure
Cr  129   55   thrpv  3 PV for cell
Cr  130   56   udiag  diagonal-only LDA+U
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C ... Local parameters
      integer i,size,nlists,psize
      parameter (size=130, nlists=56, psize=128)
      double precision sham(size)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the offe
C      data ilists /2,3,4,5,8,9,10,11,31,32,38,39,55,56,57,58,59,60,61,
C     .  62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,
C     .  83,84,94,104,105,106,107,110,114,115,116,117,127,128,129,130/
      data ilists /2,3,4,5,8,9,10,11,31,32,38,39,55,56,57,58,59,60,61,
     .  62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,
     .  83,84,85,95,105,106,107,108,111,115,116,117,118,128,129,130,131/

C ... Casts of each of the entries in offe
      data casts/4,4,2,4,4,4,4,4,2,4,31*2,5*4,2,6*4,-1,4,2/

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(noffe),0,sham,1)
        sham(1)  = size
        sham(psize) = size
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue
      end
      subroutine ulat(slat,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for lattice structure
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe   list of parameters to extract
Ci          offe(1) = -1 => initialize slat
Ci   noffe  number of entries in offe
Co   off    for each parameter in offe, returns offset to slat
Co   cast   for each parameter in offe, returns cast of entry
Co   nelt   for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1    alat    lattice parameter, in a.u.
Cr    3    2    as      dimensionless Ewald smoothing parameter
Cr    4    3    avw     average MT radius
Cr    5    4    awald   Ewald smoothing parameter
Cr    6    5    dist    deformation parameters (cf lattdf.f)
Cr   15    6    gam     lattice shear parms: gam, gx,gy,gz
Cr   19    7    gmax    cutoff gmax for Fourier transforms
Cr   20    8    ldist   switch specifying what kind of dist (cf lattdf.f)
Cr   21    9    nabc    no. divisions for F.T. mesh
Cr   24   10    ng      no. G vectors
Cr   25   11    nkd     no. direct latt. vecs. for Ewald sum
Cr   26   12    nkdmx   dimensioning for arrays holding latt. vecs
Cr   27   13    nkq     no. reciprocal latt. vecs. for Ewald sum
Cr   28   14    nkqmx   dimensioning for arrays holding latt. vecs
Cr   29   15    npgrp   Number of point symmetry group operations
Cr   30   16    nsgrp   Number of space symmetry group operations
Cr   31   17    oag     offset to symmetry group translations
Cr   32   18    obgv    phase factor sum for symmetrization of mesh rho
Cr   33   19    ocg     offset to Clebsch Gordan coeffs
Cr   34   20    ocy     offset to Ylm normalization constants
Cr   35   21    odlv    offset to direct lattice vector
Cr   36   22    ogv     offset to list of F.T. G vectors
Cr   37   23    oidxcg  offset to Clebsch Gordan indxcg
Cr   38   24    oips0   pointer to first vec in star (for symm mesh rho)
Cr   39   25    oistab  offset to site permutations table for group ops
Cr   40   26    ojcg    offset to Clebsch Gordan jcg
Cr   41   27    okv     offset to indices in list of F.T. G vectors
Cr   42   28    opos    offset to site positions
Cr   43   29    oqlv    offset to Ewald reciprocal lattice vectors
Cr   44   30    osymgr  offset to symmetry group rotation matrices
Cr   45   31    plat..  lattice vectors, units of alat
Cr   54   32    plat0.. lattice vectors before distortion
Cr   63   33    plat2.. secondary lattice vecs used in various contexts
Cr   72   34    plate.. order-N
Cr   81   35    platl.. pgf
Cr   90   36    platr.. pgf
Cr   99   37    qlat..  reciprocal lattice vectors, units 2pi/a
Cr  108   38    rpad..  truncate Ewald to rpad*rmax when lattice vector
Cr                      list has to be padded in order to include at
Cr                      least one lattice vector
Cr  109   39    size    size of this structure
Cr  110   40    slat    superlattice vectors
Cr  119   41    tol     Ewald tolerance
Cr  120   42    tolft   FT mesh tolerance
Cr  121   43    vol     cell volume
C ----------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,size,nlists,psize
      parameter (size=121, nlists=43, psize=109)
      double precision slat(size)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the offe
      data ilists /2,3,4,5,6,15,19,20,21,24,25,26,27,28,29,30,31,32,33,
     .  34,35,36,37,38,39,40,41,42,43,44,45,54,63,72,81,90,99,108,109,
     .  110,119,120,121,122/
C      data ilists /2,3,4,5,6,15,19,20,21,24,25,26,27,28,29,30,31,32,33,
C     .  34,35,36,37,38,39,40,41,42,43,44,45,54,63,72,81,90,99,108,109,
C     .  118,119,120,121/
C ... Casts of each of the entries in offe
      data casts/ 4,4,4,4,4,4,4,23*2,4,4,4,4,4,4,4,4,-1,4,4,4,4/

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(noffe),0,slat,1)
        slat(1) = size
        slat(psize) = size
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue
      end
      subroutine umix(smix,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for mixing parameters
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe  list of parameters to extract
Ci          offe(1) = -1 => initialize slat
Ci   noffe    number of entries in offe
Co   off    for each parameter in offe, returns offset to slat
Co   cast   for each parameter in offe, returns cast of entry
Co   nelt   for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   b       mixing beta
Cr    3    2   bl      previous mixing beta
Cr    4    3   bv      extra potential mixing
C     5    4   elind   Lindhard energy for model screening
Cr    6    5   fn      mixing file name
Cr    7    6   kill    kill the mixing file after k iterations
Cr    8    7   lxpot   decouple potential and the charge
Cr                     1: mix vin and v(qmix)
Cr                     2: mix vin and v(qout)
Cr    9    8   mmix    maximum number to mix
Cr   10    9   mode    1 Anderson 2 Broyden
Cr   11   10   model   previous mixing mode
Cr   12   11   n       Number of iterations for this species
Cr   13   12   nitu    max number of LDA+U itreations
Cr   14   13   nmix    actual number mixed
Cr   15   14   nsave   number of iterations to save on disk
Cr   16   15   r..     expression for rmscst
Cr   19   16   rms1    1st rms error
Cr   20   17   rms2    2nd rms error
Cr   21   18   size    total size of this struc
Cr   22   19   tj..    Anderson t's
Cr   32   20   tolu    tolerance for LDA+U
Cr   33   21   umix    mixing parameter for LDA+U
Cr   34   22   w..     Linear mixing weights
Cr   37   23   wc      Broyden weight
C ----------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,size,nlists
      parameter (size=37, nlists=23)
      double precision smix(size)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the offe
C      data ilists /2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19,20,21,31,34,
C     .  35/
      data ilists /2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,19,20,21,22,32,
     .  33,34,37,38/
C ... Casts of each of the entries in offe
      data casts /4,4,4,4,1,2,2,2,2,2,2,2,2,2,1,4,4,-1,4,4,4,4,4/

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(noffe),0,smix,1)
        smix(1) = size
        smix(21) = size
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue

      end
      subroutine umove(smove,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for dynamics
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe  list of parameters to extract
Ci         offe(1) = -1 => initialize smove
Ci   noffe number of entries in offe
Ci         offe(1) = -1 => use noffe as initialization value
Co   off   for each parameter in offe, returns offset to smove
Co   cast  for each parameter in offe, returns cast of entry
Co   nelt  for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   ct      coefficients for global deamons thermostat
Cr    5    2   gyro    gyromagnetic ratio (magnetic dynamics)
Cr    6    3   kt      temperature, atomic units
Cr    7    4   modt    mode for thermostat
Cr   10    5   nmodt   number of thermostat modes
Cr   11    6   prmint  Parameters for numerical integration
Cr                     For Bulirsch-Stoer integration:
Cr                     1:   mode: (1 for BS-integration)
Cr                     2:   rst:  1 for start new integration, 0 to continue
Cr                     3:   ts0:  minimum time step size
Cr                     4:   tol:  tolerance in integration errors
Cr                     5:   mx:   order of rational function extrapolation
Cr                     6:   mi:   number of midpoint rules
Cr                     7-17 nseq: sequence of no. midpoint divisions
Cr                     18:        offset to bs workspace
Cr   31    7   size
Cr   32    8   tnow    duration of simulation to this point, in units of ts
Cr   33    9   ts      suggested time step size
Cr   34   10   tsequ   initial equilibration time before stats accumulated
Cr   35   11   tstot   duration of total simulation, in units of ts
C ----------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,size,nlists,psize
      parameter (size=35, nlists=11, psize=31)
      double precision smove(size)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the offe
      data ilists /2,5,6,7,10,11,31,32,33,34,35,36/
C ... Casts of each of the entries in offe
      data casts/ 4,4,4,2,2,4,-1,4,4,4,4 /

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(noffe),0,smove,1)
        smove(1)  = size
        smove(psize) = size
        return
      endif


      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue
      end
      subroutine uoptic(soptic,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for optics
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe  list of parameters to extract
Ci         offe(1) = -1 => initialize soptic
Ci   noffe number of entries in offe
Ci         offe(1) = -1 => use noffe as initialization value
Co   off   for each parameter in offe, returns offset to soptic
Co   cast  for each parameter in offe, returns cast of entry
Co   nelt  for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   axes   (abc) axes in xi^(abc) for nonlinear optics
Cr   20    2   cll    core level l quantum number
Cr   21    3   cln    core level n quantum number
Cr   22    4   cls    1 or 2 for EELS or XANES
Cr   23    5   clsite site id for core level
Cr   24    6   esciss Energy shift (scissors operator)
Cr   25    7   mode   1 calculate band-to-band contribution
Cr   26    8   nchi2  number of (abc) for which to calculate xi^(abc)
Cr   27    9   ne     number of energy points
Cr   28   10   ocrng  range of occupied bands
Cr   30   11   size   size of this structure
Cr   31   12   unrng  range of unoccupied bands
Cr   33   13   window energy window for epsilon
Cr   35   ----------- points past end of structure
C ----------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,size,nlists,psize
      parameter (size=34, nlists=13, psize=30)
      double precision soptic(size)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the offe
      data ilists /2,20,21,22,23,24,25,26,27,28,30,31,33,35/
C ... Casts of each of the entries in offe
      data casts/ 2,2,2,2,2,4,2,2,2,2,-1,2,4/

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(noffe),0,soptic,1)
        soptic(1)  = size
        soptic(psize) = size
        return
      endif


      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
        if (il .eq. 1) then
C     ... Sanity check
          call rxx(soptic(25).gt.6,'uoptic: increase size of axes')
          nelt(i) = 3*max(nint(soptic(26)),0)
        endif
   10 continue
      end
      subroutine uordn(sordn,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for order-N parameters
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe  list of parameters to extract
Ci         offe(1) = -1 => initialize sordn
Ci   noffe number of entries in offe
Ci         offe(1) = -1 => use noffe as initialization value
Co   off   for each parameter in offe, returns offset to sordn
Co   cast  for each parameter in offe, returns cast of entry
Co   nelt  for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   efre    effective energy at which free electrons
Cr                     start propagating
Cr    3    2   mode    order-N mode
Cr    4    3   ncl     Number of clusters
Cr    5    4   ndofh   leading dimension to offcH
Cr    6    5   oclp
Cr    7    6   oclssl
Cr    8    7   oiaxg   neighbor table for Green's function
Cr    9    8   omapgv  maps GF pairs into iax table
Cr   10    9   ooffch  offset to table of Hamiltonian offsets
Cr   11   10   rmaxg   ranges of Green's function:
Cr                     (1) = range for direct inversion
Cr                     (2) = range for iterative inversion
Cr                     (3) = range for free-electron GF
Cr   14   11   size    size of this structure
C ----------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,size,nlists
      parameter (size=14, nlists=11)
      double precision sordn(size)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the offe
      data ilists /2,3,4,5,6,7,8,9,10,11,14,15/
C ... Casts of each of the entries in offe
      data casts/ 4,2,2,2,2,2,2,2,2,4,-1/

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(noffe),0,sordn,1)
        sordn(1) = size
        sordn(size) = size
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue
      end
      subroutine upot(spot,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for density, potential, matrix elt
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe  list of parameters to extract
Ci         offe(1) = -1 => initialize spot
Ci   noffe number of entries in offe
Ci         offe(1) = -1 => use noffe as initialization value
Co   off   for each parameter in offe, returns offset to spot
Co   cast  for each parameter in offe, returns cast of entry
Co   nelt  for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   bfield  global magnetic field direction and amplitude (4)
Cr    6    2   nlma    number of augmentation channels in system
Cr    7    3   nlml    (FP) total number of site density channels
Cr    8    4   nrhos   (FP) total number of site density channels
Cr    9    5   oaamom  offset to array containing local ASA mag.mom
Cr   10    6   obxc    offset to array containing local ASA Bxc dir.
Cr   11    7   odddpf  3rd derivative of potential functions (mkptfp)
Cr   12    8   oddpf   2nd derivative of potential functions (mkptfp)
Cr   13    9   oddpfr   2nd derivative of potential functions (mkptfp),
Cr                     fully relativistic case
Cr   14   10   odel    (tbe) potential shifts.
Cr   15   11   odpf    energy derivative of potential functions (mkptfp)
Cr   16   12   odpfr   energy derivative of potential functions (mkptfp),
Cr                     fully relativistic case
Cr   17   13   ofes    electrostatic contribution to forces
Cr   18   14   ogma    gamma-alpha in potential functions
Cr   19   15   ogmar   (gamma-alpha)*P^alpha/P^gamma fully relativistic
Cr   20   16   ogrrme  radial matrix elements of gradient.
Cr   21   17   ohab    not used now
Cr   22   18   oivso   (tbe) spin orbit
Cr   23   19   omad    Madelung matrix
Cr   24   20   oorhat  pointers to atom densities
Cr   25   21   opalp   ASA potential functions, alpha repsn (mkptfp)
Cr   26   22   opapg   p^alpha/p^gamma fully relativistic
Cr   27   23   opdel   ASA Delta, in hamiltonian order (for LDA+U)
Cr   28   24   opf     ASA potential functions (mkptfp)
Cr   29   25   opfnc   noncollinear ASA potential functions
Cr   30   26   opfr    ASA potential functions (mkptfp),
Cr                     fully relativistic case
Cr   31   27   opmpol  multipole integrals of phi,phidot
Cr   32   28   opnu    P-nu (Methfessel's log derivative function)
Cr   33   29   opp     potential parameters
Cr   34   30   oppi    not used now
Cr   35   31   oppn    NMTO generation potential parameters
Cr   36   32   opprel  Dirac potential parameters
Cr   37   33   opti    inverse potential functions
Cr   38   34   oqc     Sphere core charge
Cr   39   35   oqmom   offset to multipole moments
Cr   40   36   oqnu    Energy moments of the charge density
Cr   41   37   oqpp    Multipole moments of the nonspherical charge
Cr   42   38   oqt     Sphere total charge
Cr   43   39   orhos   spin-density matrix
Cr   44   40   orhrmx  sphere electron density at rmax
Cr   45   41   osab    not used now
Cr   46   42   osgw    structure containing GW parameters
Cr   47   43   osmpot  offset to smoothed potential
Cr   48   44   osmrho  offset to smoothed density
Cr   49   45   osop    spin-orbit parameters
Cr   50   46   osoptc  structure containing optical matrix elements
Cr   51   47   osrout  offset to smoothed output density
Cr   52   48   otau    not used now
Cr   53   49   ovab    not used now
Cr   54   50   ovdif   difference ves(q,rmax) - vactual(rmax)
Cr   55   51   oves    electrostatic potential at rmax
Cr   56   52   ovintr  intra-atomic W for screening, monopole approx
Cr   57   53   ovrmax  total potential at rmax
Cr   58   54   ovshf   ASA constant potential shifts
Cr   59   55   size    size of this structure
Cr   60   56   vconst  Constant estat potential shifts, used where
Cr                     the Fermi level is specified and the potential
Cr                     adjusts to it.
Cr                     vconst(1) = potential shift
Cr                     vconst(2) = potential shift of L end region (PGF)
Cr                     vconst(3) = potential shift of R end region (PGF)
Cr   63   57   vmtz    ASA muffin-tin zero
Cr   64   58   vmtz0   Input ASA muffin-tin zero
C -----------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,size,nlists,psize
      parameter (size=64, nlists=58, psize=59)
      double precision spot(size)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the offe
      data ilists /2,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
     .  24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,
     .  45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,63,64,65/
C      data ilists /2,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
C     .  24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,
C     .  45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,62,63,64/
C ... Casts of each of the entries in offe
      data casts/ 4,53*2,-1,3*4/

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(noffe),0,spot,1)
        spot(1) = size
        spot(psize) = size
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue
      end

      subroutine uspec(sspec,offe,noffe,is,off,cast,nelt)
C- Returns offsets, casts and number for species structure
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci      is  species index
Ci   offe  list of parameters to extract
Ci          offe(1) = -1 => initialize sspec
Ci   noffe number of entries in offe
Co   off    for each parameter in offe, returns offset to sspec
Co   cast   for each parameter in offe, returns cast of entry
Co   nelt   for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr   20     1  a      a for mesh
Cr  111..   2  alpha  screening alphas
Cr   24     3  beta   mixing beta
Cr  211..   4  chfa   coefficients to fit of free-atom density tails
Cr   46     5  colxbs color of atom (for graphics)
Cr   17     6  coreh  core hole channel (char: eg '1s')
Cr   18     7  coreq  coreq(1) = charge in core hole channel
Cr                    coreq(2) = moment in core hole channel (nsp=2 only)
Cr   26     8  ctail  coefficients to fit of free-atom core tail by unsm. Hankel
Cr   37     9  dv     constant shift in potential
Cr   53    10  eh3    sm Hankel energy for high local orbitals
Cr  341    11  ehvl   l-dependent sm Hankel energies for val-lap basis
Cr  291    12  enu    linearization energies
Cr   22    13  eref   reference energy
Cr   27    14  etail  energy to fit of free-atom core tail
Cr   23    15  etf    Thomas-Fermi screening length
Cr  171..  16  exi    Hankel energies for fit to c.d.
Cr                    For free atoms, fit to free-atom density tails.
Cr   39    17  group  group
Cr   40    18  grp2   grp2
Cr  161..  19  hcr    kinky wave hard sphere radii, atomic units
Cr  121..  20  idmod  idmol(l) controls how linearization energy is
Cr                    determined for an orbital of angular momentum l
Cr                    0 float to center of gravity of occupied part of band
Cr                    1 freeze to spec'd logarithmic derivative
Cr                    2 freeze to spec'd linearization energy
Cr                    3 float to natural center of band (C parameter)
Cr  321    21  idu    identifies l-channels with Hubbard U (LDA+U)
Cr  131..  22  idxdn  idxdn
Cr  333..  23  iq1    (tbfit) switches to freeze TB site energies
Cr  337..  24  ivso   (tbfit) switches to freeze TB spin-orbit params
Cr  329    25  jh     LDA+U J parameters for each l-channel
Cr   12    26  kmxh   k cutoffs for head augmentation expansion
Cr   14    27  kmxt   k cutoffs for tail augmentation expansion
Cr   56    28  kmxv   k-cutoff for 1-center projection of free-atom rho
Cr   25    29  lfoca  switch specifying treatment of core density
Cr   11    30  lmxa   l cutoff for augmentation expansion
Cr   31    31  lmxb   highest l in basis
Cr   32    32  lmxf   fitting l
Cr   13    33  lmxl   cutoff for local density
Cr   43    34  lmxpb  l-cutoff for product basis
Cr   49    35  lxi    l-cutoffs for interstital product basis
Cr    3    36  mass   atomic mass
Cr   38    37  mxcst  species-specific constraints of various types:
Cr                    bit function
Cr                    1   suppresses c.d. mixing for this species (ASA)
Cr                    2   Exclude this species when auto-resizing sphere radii
Cr                    4   Freeze augmentation w.f. for this species (FP)
Cr   28    38  name   species name
Cr    6    39  naug   number of k,L for augmentation
Cr  181    40  ngcut  orbital-dependent G cutoffs (for nfp basis)
Cr    7    41  norb   number of orbitals in basis
Cr   16    42  norp   number of parameters needed to define orbital
Cr   21    43  nr     nr for mesh
Cr    5    44  ntorb  number of orbital types in basis
Cr   41    45  nxi    Number of energies in fit of free-atom density tails
Cr  231..  46  orbp   (lmf) parameters defining orbital shape
Cr   29    47  orhoc  pointer to core density
Cr   71..  48  p      log derivative parameter
Cr   44    49  pb1    product basis for linear response
Cr   45    50  pb2    second product basis for linear response
Cr  301... 51  pz     log derivative parameter, semicore states
Cr   91..  52  q      starting q's (charges)
Cr   42    53  qc     core charge
Cr  181..  54  qpol   (tbe) 10 polarisability parameters
Cr   50    55  radxbs radius of atom (for graphics)
Cr   57    56  rcfa   renormalization radius of free atom density, and width
Cr   51    57  rcut   range over which true, smoothed Hankels differ---may
Cr                    be used in the calculation of XC potential
Cr   36    58  rfoca  smoothing radius for frozen core overlap approx
Cr   10    59  rg     rsm for gaussians to fix multipole moments
Cr   35    60  rham   range of hamiltonian
Cr   52    61  rint   range of interstitial product basis
Cr    4    62  rmt    augmentation radius
Cr   54    63  rs3    Lower bound to rsm for local orbital
Cr    9    64  rsma   rsm for augmentation expansion
Cr    8    65  rsmfa  rsm to fit free atom density
Cr   15    66  rsmv   rsm for one-center projection of free-atom density
Cr                    In TCF, smoothing radius to calculate exchange
Cr                    (artificially smoothed density)
Cr    1    67  size   leading dimension of sspec
Cr   30    68  stc    core kinetic energy
Cr   33    69  stni   (tbe) Stoner I
Cr   34    70  stnm   Stoner mmax
Cr  325    71  uh     Hubbard U
Cr   55    72  vmtz   Asymptotic potential for fitting functions at rmt
Cr  341    73  vso    (tbe) spin-orbit parameters
Cr    2    74  z      atomic number
C ----------------------------------------------------------------
C     implicit none
      integer is,noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,spsiz,nlists
      integer orbp,ochfa,op,oq,oidx,opz,ongcut,ocxbs,oidu,ojh,ouh
      integer ocoreq,oorcfa,oiq1,oivso,ovso
      parameter (spsiz=350, nlists=74)
C     Next 3 lines correspond to right column offsets
      parameter (ochfa=4,op=48,opz=51,oq=52,ocxbs=5)
      parameter (oidx=22,orbp=46,ongcut=40,ocoreq=7)
      parameter (oidu=21,ojh=25,ouh=71,oiq1=23,oivso=24,ovso=73)
C     Next line correspond to left column offsets
      parameter (oorcfa=57)
      double precision sspec(spsiz,1)
      integer ilists(nlists),casts(nlists),il,n0,nkap0
      parameter(n0=10,nkap0=3)
C ... Pointers to locations in structure for each of the offe
C      data ilists /20,111,24,211,46,17,18,26,37,53,291,22,27,23,171,39,
C     .  40,161,121,321,131,333,337,329,12,14,56,25,11,31,32,13,43,49,3,
C     .  38,28,6,181,7,16,21,5,41,231,29,71,44,45,301,91,42,181,50,57,51,
C     .  36,10,35,52,4,54,9,8,15,1,30,33,34,325,55,341,2/
      data ilists /20,111,24,211,46,17,18,26,37,53,341,291,22,27,23,171,
     .  39,40,161,121,321,131,333,337,329,12,14,56,25,11,31,32,13,43,49,
     .  3,38,28,6,181,7,16,21,5,41,231,29,71,44,45,301,91,42,181,50,57,
     .  51,36,10,35,52,4,54,9,8,15,1,30,33,34,325,55,341,2/

C ... Casts of each of the entries in offe
      data casts/ 5*4,1,10*4,2,2,4,5*2,4,10*2,4,2,1,7*2,4,2,
     .            4,1,1,16*4,-1,7*4/

      if (offe(1) .eq. -1) then
        call dcopy(spsiz*is,dble(noffe),0,sspec,1)
        sspec(1,1) = spsiz
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        nelt(i) = 1
        if (casts(il) .lt. 0) nelt(i) = 0
        if (ilists(il) .gt. 70) nelt(i) = n0
        if (ilists(il) .eq. oorcfa) nelt(i) = 2
        if (il.eq.ocxbs) nelt(i) = 3
        if (il.eq.op .or. il.eq.opz .or. il.eq.oq .or. il.eq.ochfa)
     .    nelt(i) = 2*n0
        if (il.eq.oidx .or. il.eq.ongcut) nelt(i) = nkap0*n0
        if (il.eq.orbp) nelt(i) = nkap0*n0*2
        if (il.eq.oidu .or. il.eq.ojh .or. il.eq.ouh .or. il.eq.oiq1
     .    .or. il.eq.oivso .or. il.eq.ovso) nelt(i) = 4
        if (il.eq.ocoreq) nelt(i) = 2
        cast(i) = casts(il)
        off(i) = ilists(il)
   10 continue
      end
      subroutine usite(ssite,offe,noffe,is,off,cast,nelt)
C- Returns offsets, casts and number for site structure
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci      is  site index
Ci   offe  list of parameters to extract
Ci          offe(1) = -1 => initialize sspec
Ci   noffe  number of entries in offe
Co   off    for each parameter in offe, returns offset to sspec
Co   cast   for each parameter in offe, returns cast of entry
Co   nelt   for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   bfield  site magnetic field (not used)
Cr    5    2   clabel  string labelling class
Cr    6    3   class   class index
Cr    7    4   cli     (RGF) cluster index
Cr    8    5   delta   (tbe) vector of on-site energy shifts
Cr   14    6   dpole   point charge dipole (molecules)
Cr   17    7   eula    Euler angles (noncollinear magnetism)
Cr   20    8   force   force (dynamics)
Cr   23    9   mpole   point charge monopole (molecules)
Cr   24   10   ndelta  size of delta
Cr   25   11   norb    dimension contributed by this site to hamilt.
Cr   26   12   offh    hamiltonian offset
Cr   27   13   orho    pointer to MT density
Cr   28   14   ov0     pointer to potential that defines wave functions
Cr   29   15   ov1     pointer to spherical part of MT potential
Cr   30   16   pl      (PGF) principal layer index
Cr   31   17   plv     (PGF) principal layer potential index
Cr   32   18   pnu     log derivative parameter
Cr   52   19   pos     true coordinates of atom
Cr   55   20   pos0    coordinates of atom before distortion
Cr   58   21   pz      log derivative parameter, semicore states
Cr   78   22   relax   (dynamics) flags which coordinates to relax
Cr   81   23   sid     (RGF) Site id, used to merge coincident clusters
Cr   82   24   size    size of this structure
Cr   83   25   spec    species index
Cr   84   26   vel     velocity
Cr   87   27   vshft   constant potential shift for this site
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer is,noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C ... Local parameters
      integer i,size,nlists,odelt,ndelti,psize
      parameter (size=87,nlists=27,ndelti=5,odelt=24,psize=82)
      double precision ssite(size,1)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the offe
C      data ilists /2,5,6,7,8,14,17,20,21,22,23,24,25,26,27,28,48,51,54,
C     .  74,77,78,79,80,83,84/
      data ilists /2,5,6,7,8,14,17,20,23,24,25,26,27,28,29,30,31,32,52,
     .  55,58,78,81,82,83,84,87,88/
C ... Casts of each of the entries in offe
      data casts/ 4,1,2,2,5*4,8*2,4*4,2,2,-1,2,4,4/

      if (offe(1) .eq. -1) then
        call dcopy(size*is,dble(noffe),0,ssite,1)
        ssite(1,1) = size
        ssite(psize,1) = size
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
C   ... Handle special delta
        if (il .eq. ndelti) nelt(i) = nint(ssite(odelt,is))
   10 continue
      end
      subroutine ustr(sstr,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for tb and TCF strux parameters
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe  list of parameters to extract
Ci         offe(1) = -1 => initialize sstr
Ci   noffe number of entries in offe
Ci         offe(1) = -1 => use noffe as initialization value
Co   off   for each parameter in offe, returns offset to sstr
Co   cast  for each parameter in offe, returns cast of entry
Co   nelt  for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   adec    log spacing in distance for TCF fit
Cr    3    2   amode   mode for calculating screening alpha's
Cr                    1s digit
Cr                     0: alpha=alp0, (which is input)
Cr                     1: alpha=gamma, calculated from pot pars.
Cr                        NO LONGER SUPPORTED
Cr                     2: alpha calculated for hard spheres.
Cr                        NO LONGER SUPPORTED; see mktral
Cr                    10s digit 0: adot=Andersen's localized
Cr                     0: adot=Andersen's localized
Cr                     1: adot so Kanpur notes Eq 3.87 is zero
Cr                     2: adot calculated for hard spheres.
Cr                        NO LONGER SUPPORTED; see mktral
Cr    4    3   drwats  padding beyond cluster for Watson sphere
Cr    5    4   iinv    Parameters for iterative inversion of tb-strx
Cr                     arg 1: number of sites for cutoff
Cr                     arg 2: number of iterations
Cr                     arg 3: allowed tolerance
Cr                     arg 4-5: not used now.
Cr   10    5   ivl     Type of val-lap functions
Cr   11    6   kaps    Kinetic energies at which strux evaluated (NMTO)
Cr   17    7   lequiv  look for equivalent connecting vectors
Cr   18    8   lmaxw   lmax for Watson sphere
Cr   19    9   loka    conventions for Hankel and Bessel functions
Cr                     0 Methfessel's conventions
Cr                     1 Andersen's 2nd gen LMTO conventions
Cr   20   10   lshow   show strux
Cr   21   11   mxnbr   maximum number of neighbors
Cr   22   12   nalf    polynomial order for d-dependence of TCF fit coffs
Cr   23   13   nbisi   mesh parameters for TCF
Cr   26   14   ncupl   reduced polynomial order in TCF
Cr   27   15   ndust   no. iterations to improve on TCF
Cr   28   16   nitab   ?? same as nttab
Cr   29   17   nkaps   number of kaps
Cr   30   18   nttab   number of pairs in this table
Cr   31   19   oadot   offset to screening alpha's, e. derivative
Cr   32   20   oalp    offset to screening alpha's
Cr   33   21   oiax    offset to neighbor table
Cr   34   22   ontab   offset to number-of-pairs table
Cr   35   23   os      offset to strux
Cr   36   24   osdot   offset to strux energy derivative
Cr   37   25   rfit    ?? (not used) default hard core radius
Cr   38   26   rmax    range of strux (tbe: range of hamiltonian)
Cr   39   27   rmaxg   range of val-lap strux
Cr   40   28   size    size of this structure
Cr   41   29   skmsh   parameters for 3rd gen LMTO energy interpolation
Cr                     not used?
Cr   47   30   wztcf   weighting of z-points in TCF
C ----------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,size,nlists,psize
      parameter (size=47, nlists=30, psize=40)
      double precision sstr(size)
      integer ilists(nlists+1),casts(nlists),il
C ... Pointers to locations in structure for each of the offe
C      data ilists /2,3,4,5,10,16,17,18,19,20,21,22,25,26,27,28,29,30,31,
C     .  32,33,34,35,36,37,38,39,45,46/
      data ilists /2,3,4,5,10,11,17,18,19,20,21,22,23,26,27,28,29,30,31,
     .  32,33,34,35,36,37,38,39,40,41,47,48/
C ... Casts of each of the entries in offe
      data casts/ 4,2,4,4,2,4,18*2,4,4,4,-1,4,4/

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(noffe),0,sstr,1)
        sstr(1)  = size
        sstr(psize) = size
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue
      end
      subroutine utb(stb,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for TB s
C ----------------------------------------------------------------
Ci  Inputs/Outputs
Ci   offe  list of parameters to extract
Ci         offe(1) = -1 => initialize stb
Ci   noffe number of entries in offe
Ci         offe(1) = -1 => use noffe as initialization value
Co   off   for each parameter in offe, returns offset to stb
Co   cast  for each parameter in offe, returns cast of entry
Co   nelt  for each parameter in offe, returns number of entries
Cr Remarks.
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr    2    1   alam   (tbfit) lambda variable, Levenburg-Marquardt (Numerical Recipes)
Cr    3    2   alsc   (tbfit) scale factor for Levenburg-Marquardt lambda
Cr    4    3   ebfit  (tbfit) energy window for fitting eigenvalues
Cr                    lower and upper limit of eigenvalues to be fit
Cr    6    4   fmode  (tbfit) mode for fitting TB parmeters
Cr                    =0 fit parameters used in tbham, =1 not implemented yet
Cr    7    5   nbfit  (tbfit) min/max band numbers to fit
Cr                    lower and upper band indices to be fit for all k-points
Cr    9    6   rmfit  (tbfit) real space range to fit TB parameters (fmode=1)
Cr   11    7   size   size of this structure
C ----------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,size,nlists
      parameter (size=11, nlists=7)
      double precision stb(size)
      integer ilists(nlists+1),casts(nlists),il
      integer psize
      parameter (psize=11)
C ... Pointers to locations in structure for each of the offe
      data ilists /2,3,4,6,7,9,11,12/
C ... Casts of each of the entries in offe
      data casts /4,4,4,2,2,4,-1/

      if (offe(1) .eq. -1) then
        call dcopy(size,dble(noffe),0,stb,1)
        stb(1)  = size
        stb(psize) = size
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = casts(il)
        off(i) = ilists(il)
        nelt(i) = ilists(il+1)-ilists(il)
   10 continue

      end
      subroutine ustrn(sstrn,offe,noffe,off,cast,nelt)
C- Returns offsets, casts and number for string structure
C ----------------------------------------------------------------
Ci Inputs (see Remarks)
Ci   sstrn  not used; present only for compatibility with other
Ci                    unpacking routines.
Ci   offe   list of indices offe(i) for i=1..offe (see Remarks)
Ci          offe(1) = 0 => do nothing
Ci          offe(i)>0 unpacks offset and length of string in off,nelt
Ci                    Here offe(i) are the list of i=1..noffe indices
Ci          offe(i)<0 assigns offset and length of string from off,nelt
Ci                    Here -offe(i) are the list of i=1..noffe indices
Ci   noffe  number of entries in offe
Ci   nelt   (for each offe, if offe<0) string length
Co  Outputs
Co   For each i, i=1..noffe, return the following. If offe(i) > 0:
Co     off    for each parameter in offe, returns offset to sstrn
Co     cast   for each parameter in offe, returns cast=1 of entry
Co     nelt   for each parameter in offe, returns length of string
Co   If offe(i) < 0:
Co     off    for each parameter in offe, returns offset to sstrn
Co     cast   for each parameter in offe, returns cast=1 of entry
Cr Remarks.
Cr   sstrn holds a set of starting and ending points in a large char
Cr   array (not passed here) corresponding to strings labeled below.
Cr   lstra associates string labels with index (eg 'mix' with '6')
Cr
Cr   The string structure is different from the others, because
Cr   of the screwy way fortran handles string arguments.
Cr
Cr   ustrn takes a list of noffe indices (passed in offe)
Cr   and, if offe>0
Cr       returns a list of starting positions (in off(1..noffe)) and
Cr       string lengths (in nelt(1..offe)) in the global character array
Cr       (not passed here) that consitute strings associated w/ indices
Cr   or, if offe<0
Cr       assigns offsets and lengths (in off(1..noffe) and nelt(1...noffe)
Cr       based for future storage.
Cr
Cr   Caller must actually copy the string into or out of the structure.
Cr
Cr   Example: copy string instr(i..strsiz) to 'mmham' slot
Cr   Here sstrng is name of char array holding all strings.
Cr       ... Get offs(1) = index to mmham
Cr       call lstra('strn mmham',tmp,offs,i)
Cr       ... Allocate string of length 1 to get offset into offi (i1
Cr       call ustrn(tmp,-offs,1,offi,tmp,1)
Cr       ... copy instr(i..strsiz) to sstrng, compressing spaces.
Cr       call cpstr(instr,strsiz,101,'}',i,j,sstrng(offi(1):))
Cr       ... j characters were copied; reallocate string of proper length
Cr       call ustrn(tmp,-offs,1,offi,i1,1)
Cr
Cr   The index to the elements are arranged alphabetically by name:
Cr   The following maps the internal location to the index:
Cr   off offe  name    purpose
Cr         1   amix    amix mode
Cr         2   gemb    file names for embedding GF
Cr         3   gfopt   GF options
Cr         4   jobid   job id
Cr         5   map     mapping
Cr         6   mix     mixing mode
Cr         7   mmham   micromagnetics hamiltonian
Cr         8   sxopt   SX options
Cr         9   symg    symmetry group
C ----------------------------------------------------------------
C     implicit none
      integer noffe,offe(noffe),off(noffe),cast(noffe),nelt(noffe)
C Local variables
      integer i,nlists
      parameter (nlists=9)
      double precision sstrn(1)
      integer ilists(nlists),ilen(nlists),il,offmx,k
      save ilists,ilen,offmx
      data offmx /0/ ilists /nlists*0/

      if (offe(1) .eq. 0) then
C       size = noffe
        return
      endif

      do  10  i = 1, noffe
        il = offe(i)
        cast(i) = 1
        if (il .lt. 0) then
C     ... See whether to re-allocate string.
          if (ilists(-il) .ne. 0) then
            do  12  k = 1, nlists
C       ... Branch if not the top one; make a new allocation
   12       if (ilists(-il) .lt. ilists(k)) goto 14
C       ... It is the top one; just remake ilen and offmx
            goto 16
          endif
C     ... Entry point for new string allocation
   14     ilists(-il) = offmx+1
C     ... Entry point for last string re-allocation
   16     ilen(-il)   = nelt(i)
          offmx  =  ilists(-il)+nelt(i)
          off(i) =  ilists(-il)
        else
          off(i) =  ilists(il)
          nelt(i) = ilen(il)
        endif
   10 continue

      end

