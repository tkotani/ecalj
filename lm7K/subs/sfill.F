Cgetarg...info...           structure ['sfill', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['sfill', 'array', 'sarray']
Cgetarg...info...           structure ['sfill', 'lat', 'slat']
Cgetarg...info...           structure ['sfill', 'site', 'ssite']
Cgetarg...info...           structure ['sfill', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['sfill', 'm_struc_def']
      subroutine sfill(sfargs,slabl,sctrl,sarray,slat,sspec,ssite)
       
       use m_struc_def  !Cgetarg

C- Adjust sphere sizes to fill to specified fraction of cell volume
C ----------------------------------------------------------------
Ci Inputs
Ci   sfargs string specifying how rmax to be scaled.
Ci          First character is the delimiter separating switches.
Ci          Last argument specifies which classes are to be included.
Ci          For now, there are no special switches.
Ci   slabl :vector of species labels
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nclass nspec modep sclwsr omax1 omax2 wsrmax
Ci     Stored:    *
Ci     Passed to: *
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: npadl npadr nclasp oipc oics oclabl ormax
Ci     Stored:    ormax
Ci     Passed to: *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: avw alat plat vol opos
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: z mxcst rmt
Ci     Stored:    rmt
Ci     Passed to: *
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: *
Co Outputs
Co   sspec->rmax and sarray->rmax
Cr Remarks
Cr   Using automatic scaling routine sclwsr, program tries to
Cr   rescale sphere radii so that sum-of-sphere volumes = specified
Cr   volume fraction of cell volume.
Cr   Constraints are imposed through the following (see routine sclwsr)
Cr     1 maximum overlaps may not exceed omax1 and omax2
Cr     2 specific radii maybe locked (2's digit of spec->mxcst)
Cr     3 If wsrmax is nonzero, no sphere may exceed wsrmax.
Cr
Cr   OLD method: user supplies class list and program scales those
Cr   classes to satisfy volume requirements.  This method is kept
Cr   for now, but will be discarded in future.
Cr   A subset of the nclass rmax are scaled to fulfill some criterion.
Cr   For now rmax are scaled fill the cell volume.
Cu Updates
Cu   17 May 02  Modified MT radii scaling to lower priority for E.S.
Cu   23 Apr 02  Added automatic rescaling using sclwsr.  New arg list.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
Cgetarg       double precision sctrl(1),sarray(1),slat(1),sspec(1),ssite(1)
       type(s_ctrl)::sctrl
       type(s_array)::sarray
       type(s_lat)::slat
       type(s_spec)::sspec(*)
       type(s_site)::ssite(*)

      character*(*) sfargs
      character*8 slabl(*)
C Heap:
      integer w(1)
      common /w/ w
C Local:
      character dc*1
      integer nbas,nbasp,nclass,nspec,npadl,npadr,nclasp
      integer ls,lstyle,i,iv,j1,j2,m,parg,nlist,lgunit,ipr,ic,
     .  ormax,ormx2,olist
      integer modep(3),bitand,igetss,isw
      integer oips,oipc,opos,oics,oclabl,oz,ormt,olock,ormax2,opt
      double precision plat(3,3),alat,avw,vol,rmax,wsrmax,omax1(3),
     .  omax2(3),fovl,vol0,vols,scale,dval,volfac
      character*8 clabl

Cgetarg...info...           integer :: nbas , integer(8) :: sctrl%nbas 1
Cgetarg...info...           integer :: nclass , integer(8) :: sctrl%nclass 1
Cgetarg...info...           integer :: nspec , integer(8) :: sctrl%nspec 1
Cgetarg...info...           integer :: modep(3) , integer(8) :: sctrl%modep 3
Cgetarg       call upack('ctrl nbas nclass nspec modep',sctrl,nbas,nclass,nspec,
Cgetarg      .  modep,0)
       
       nbas=sctrl%nbas
       nclass=sctrl%nclass
       nspec=sctrl%nspec
       i_copy_size=size(sctrl%modep) 
       call i8icopy(i_copy_size,sctrl%modep,1,modep,1) 

Cgetarg...info...           integer :: npadl , integer(8) :: sarray%npadl 1
Cgetarg...info...           integer :: npadr , integer(8) :: sarray%npadr 1
Cgetarg...info...           integer :: nclasp , integer(8) :: sarray%nclasp 1
Cgetarg       call upack('array npadl npadr nclasp',sarray,npadl,npadr,nclasp,0,
Cgetarg      .  0)
       
       npadl=sarray%npadl
       npadr=sarray%npadr
       nclasp=sarray%nclasp

Cgetarg...info...           integer :: oipc , integer(8) :: sarray%oipc 1
Cgetarg...info...           integer :: oics , integer(8) :: sarray%oics 1
Cgetarg...info...           integer :: oclabl , integer(8) :: sarray%oclabl 1
Cgetarg...info...           integer :: ormax , integer(8) :: sarray%ormax 1
Cgetarg       call upack('array oipc oics oclabl ormax',sarray,oipc,oics,oclabl,
Cgetarg      .  ormax,0)
       
       oipc=sarray%oipc
       oics=sarray%oics
       oclabl=sarray%oclabl
       ormax=sarray%ormax

Cgetarg...info...           real(8) :: avw , real(8) :: slat%avw 1
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: plat(3,3) , real(8) :: slat%plat 9
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg...info...           integer :: opos , integer(8) :: slat%opos 1
Cgetarg       call upack('lat avw alat plat vol opos',slat,avw,alat,plat,vol,
Cgetarg      .  opos)
       
       avw=slat%avw
       alat=slat%alat
       i_copy_size=size(slat%plat) 
       call dcopy(i_copy_size,slat%plat,1,plat,1) 
       vol=slat%vol
       opos=slat%opos

      call sp2cls('spec z',sspec,w(oics),1,1,nclasp,oz)
      nbasp = nbas + npadl + npadr
      call getpr(ipr)

C ... Do nothing if no specifications
      if (sfargs .eq. ' ') return
      ls = len(sfargs)
      if (sfargs(1:4) .ne. 'auto') goto 100

C ... Use sclwsr
C  10 continue
Cgetarg...info...           real(8) :: volfac , real(8) :: sctrl%sclwsr 1
Cgetarg...info...           real(8) :: omax1(3) , real(8) :: sctrl%omax1 3
Cgetarg...info...           real(8) :: omax2(3) , real(8) :: sctrl%omax2 3
Cgetarg...info...           real(8) :: wsrmax , real(8) :: sctrl%wsrmax 1
Cgetarg       call upack('ctrl sclwsr omax1 omax2 wsrmax',sctrl,volfac,omax1,
Cgetarg      .  omax2,wsrmax,0)
       
       volfac=sctrl%sclwsr
       i_copy_size=size(sctrl%omax1) 
       call dcopy(i_copy_size,sctrl%omax1,1,omax1,1) 
       i_copy_size=size(sctrl%omax2) 
       call dcopy(i_copy_size,sctrl%omax2,1,omax2,1) 
       wsrmax=sctrl%wsrmax

      opt = volfac/10
      opt = 10*opt
      volfac = volfac-opt
      if (volfac .eq. 0) return
      call defrr(opos,3*nbasp)
Cgetarg       call spackv(10,'site pos',ssite,1,nbasp,w(opos))
       i_copy_size=size(ssite(1)%pos)
       do i_spackv=1,nbasp 
       call spackv_array_copy_r8_r8('u',ssite(i_spackv)%pos,i_copy_size,i_spackv+1-1,w(opos))
       enddo

      call defrr(oips,nbasp)
Cgetarg       call spackv(10,'site spec',ssite,1,nbasp,w(oips))
       i_copy_size=1; 
       do i_spackv=1,nbasp 
       call spackv_array_copy_i8_i('u',ssite(i_spackv)%spec,i_copy_size,i_spackv+1-1,w(oips))
       enddo

      call defrr(oz,nspec)
Cgetarg       call spackv(10,'spec z',sspec,1,nspec,w(oz))
       i_copy_size=1; 
       do i_spackv=1,nspec 
       call spackv_array_copy_r8_r8('u',sspec(i_spackv)%z,i_copy_size,i_spackv+1-1,w(oz))
       enddo

      call defrr(ormt,nspec)
Cgetarg       call spackv(10,'spec rmt',sspec,1,nspec,w(ormt))
       i_copy_size=1; 
       do i_spackv=1,nspec 
       call spackv_array_copy_r8_r8('u',sspec(i_spackv)%rmt,i_copy_size,i_spackv+1-1,w(ormt))
       enddo

      call defi(olock,-nspec)
      do  i = 1, nspec
Cgetarg         call ivset(w(olock),i,i,bitand(igetss('spec mxcst',i,sspec),2))
         call ivset ( w ( olock ) , i , i , bitand ( int(sspec(i)%mxcst) 
     .   , 2 ) ) 

      enddo
      i = 3
      if (modep(3) .ne. 2) i = 2
      if (ipr .ge. 30) then
        call awrit2('%N SFILL:  automatic scaling of sphere radii.'//
     .    '%?#n#  Constrain rmax<=%d##',
     .    ' ',160,lgunit(1),isw(wsrmax.ne.0),wsrmax)
      endif
      call sclwsr(opt,nbas,nbasp,nspec,alat,plat,w(opos),w(oips),modep,
     .  slabl,w(oz),w(olock),volfac,wsrmax,omax1,omax2,w(ormt))
      call rlse(opos)
      if (i .eq. 2) call rx('sfill needs update for 2D')
C     Poke back into spec->rmt
Cgetarg       call spackv(11,'spec rmt',sspec,1,nspec,w(ormt))
       i_copy_size=1; 
       do i_spackv=1,nspec 
       call spackv_array_copy_r8_r8('p',sspec(i_spackv)%rmt,i_copy_size,i_spackv+1-1,w(ormt))
       enddo

C     nclspp = max(2*nclasp-nclass,nspec)
C     Poke back into array->rmax
      call sp2cls('spec rmt',sspec,w(oics),1,1,nclasp,ormax2)
      call dcopy(nclasp,w(ormax2),1,w(ormax),1)
Cgetarg...info...           integer :: ormax , integer(8) :: sarray%ormax 1
Cgetarg       call pack1('array ormax',sarray,ormax)
       
       sarray%ormax=ormax 

      call rlse(opos)
      return

C --- OLD style ---
  100 continue
      j1 = 0
      call skipbl(sfargs,ls,j1)
      j1 = j1+1
      dc = sfargs(j1:j1)
      j1 = j1+1
      lstyle = 1

C ... Return here to resume parsing for arguments
  101 continue
      call nwordg(sfargs,0,dc//' ',1,j1,j2)
      j2 = j2+1

C ... Parse special arguments
      if (sfargs(j2:j2) .ne. ' ')  then
        m = j1-1
        i = parg('style=',2,sfargs,m,ls,dc,1,1,iv,lstyle)
        j1 = j2+1
        goto 101
      endif

C --- Change subset of rmax to fill cell volume ---
C ... Make the list of classes to include in the scaling, backup rmax
      call defi(olist, nclass)
      call defrr(ormx2,nclass)
      call clist(lstyle,sfargs(j1:j2),w(oclabl),w(oz),nclass,nlist,
     .  w(olist))
      call dcopy(nclass,w(ormax),1,w(ormx2),1)

C ... Cell volume
C ... Sum of sphere volumes
      call pshpr(0)
      call ovlchk(nbas,nbasp,w(opos),alat,w(ormax),w(ormax),w(oclabl),
     .  w(oipc),modep,plat,fovl,vols)

C ... vol0 <- volume of all spheres not in list
      call psfil2(nclass,nlist,w(olist),1d0,w(ormx2),0d0,w(ormx2),
     .  w(ormax))
      call ovlchk(nbas,nbasp,w(opos),alat,w(ormax),w(ormax),w(oclabl),
     .  w(oipc),modep,plat,fovl,vol0)

C ... scale <- relative change in sphere sizes to make vols = vol
      scale = ((vol-vol0)/(vols-vol0))**(1d0/3d0)
      call psfil2(nclass,nlist,w(olist),1d0,w(ormx2),scale,w(ormx2),
     .  w(ormax))

C ... Poke back into sspec
      do  120  ic = 1, nclasp
        rmax = dval(w(ormax),ic)
Cgetarg...info...           real(8) :: rmax , real(8) :: sspec%rmt 1
Cgetarg         call pack2('spec rmt',sspec,ic,rmax)
         
         sspec(ic)%rmt=rmax 

  120 continue

      if (ipr .ge. 10) then
        call awrit6('%N SFILL:  Cell vol %;8g   Sum of sphere vol '//
     .    '%;8g (%1;5d)%N%9fRelative vol of %i classes to resize: %;4d'
     .    //'  scale=%;5d',' ',160,lgunit(1),vol,vols,vols/vol,
     .    nlist,(vols-vol0)/vols,scale)
C        call awrit2(' change rmax in classes:  %n:1i',' ',80,
C     .    lgunit(1),nlist,w(olist))
        if (ipr .ge. 30 .and. dabs(scale-1d0) .gt. 1d-5) then
          call awrit0('  Class%8fOld rmax    New rmax',' ',80,
     .      lgunit(1))
          do  130  i = 1, nlist
            ic = w(olist+i-1)
            call r8tos8(dval(w(oclabl),ic),clabl)
            print 351, ic,clabl,dval(w(ormx2),ic),dval(w(ormax),ic)
  351       format(i3,2x,a,f10.6,f12.6)
  130     continue
        endif
      endif

C ... Debugging to make sure it worked
      call ovlchk(nbas,nbasp,w(opos),alat,w(ormax),w(ormax),w(oclabl),
     .  w(oipc),modep,plat,fovl,vols)

      call poppr
Cgetarg       end
       end subroutine sfill 

      subroutine psfil2(nclass,nlist,list,s1,rmax1,s2,rmax2,rmax)

C- Copy s1*rmax1 to rmax, except for elements in list copy rmax2*scale
C     implicit none
      integer nclass,nlist,list(nlist)
      double precision s1,s2,rmax1(nclass),rmax2(nclass),rmax(nclass)
      integer i

      call dpcopy(rmax1,rmax,1,nclass,s1)
      do  10  i = 1, nlist
   10 rmax(list(i)) = rmax2(list(i))*s2

Cgetarg       end
       end subroutine psfil2 


