Cdelw1 var ['ocn', '0', 'ocn', 'ng', 'complex(8)', 'zv_w_', 'mixrho']
Cdelw1 var ['owk', '0', 'owk', 'kkk', 'complex(8)', 'zv_w_', 'mixrho']
Cdelw1 var ['okv', '-1', 'undef', '0', 'undef', 'x', 'mixrho']
Cdelw1 var ['ocg1', '0', 'ocg1', '-ng', 'complex(8)', 'zv_w_', 'mixrho']
Cdelw1 var ['ocg2', '0', 'ocg2', '-ng', 'complex(8)', 'zv_w_', 'mixrho']
Cdelw1 var ['ogv', '-1', 'undef', '0', 'undef', 'x', 'mixrho']
Cdelw1 var ['ofkl', '0', 'ofkl', '(k0+1)*nlm0*nbas', 'complex(8)', 'zv_w_', 'mixrho']
Cdelw1 var ['orofi', '1', 'orofi', 'nr', 'real(8)', 'rv_w_', 'mixrho']
Cdelw1 var ['orwgt', '1', 'orwgt', 'nr', 'real(8)', 'rv_w_', 'mixrho']
Cdelw1 var ['oips0', '0', 'oips0', 'ng', 'integer', 'iv_w_', 'mixrho']
Cdelw1 var ['oco', '0', 'oco', 'ng*2*nsp', 'real(8)', 'rv_w_', 'mixrho']
Cdelw1 var ['oa', '0', 'oa', '-nda*nsp*(mxsav+2)*2', 'real(8)', 'rv_w_', 'mixrho']
Cdelw1 var ['oqkl', '0', 'oqkl', '-(kmxr+1)*nlm0*nsp*4*nbas', 'complex(8)', 'zv_w_', 'mixrho']
Cdelw1 var ['oaa', '1', 'oaa', '-naa*(mxsav+2)*2', 'real(8)', 'rv_w_', 'mixrho']
Cdelw1 var ['owk1', '2', 'owk1', 'nr', 'real(8)', 'rv_w_', 'mixrho']
Cdelw1 var ['owk2', '2', 'owk2', 'nr', 'real(8)', 'rv_w_', 'mixrho']
Cdelw1 var ['ori', '1', 'ori', 'nr', 'real(8)', 'rv_w_', 'pvmix3']
Cdelw1 var ['orwgt', '1', 'orwgt', 'nr', 'real(8)', 'rv_w_', 'pvmix3']
Cdelw1 var ['orsm', '3', 'orsm(m)', 'nr*nlml*nsp', 'real(8)', 'rv_w_', 'pvmix3']
Cdelw1 var ['orofi', '2', 'orofi', 'nr', 'real(8)', 'rv_w_', 'pvmix5']
Cdelw1 var ['orofi', '1', 'orofi', 'nr', 'real(8)', 'rv_w_', 'pvmix7']
Cdelw1 var ['orwgt', '1', 'orwgt', 'nr', 'real(8)', 'rv_w_', 'pvmix7']
Cdelw1 var ['orsm', '3', 'orsm(m)', 'nr*nlml*nsp', 'real(8)', 'rv_w_', 'pvmix7']
Cdelw1 undel [['undef', 'wref'], 'okv', 'mixrho']
Cdelw1 undel [['undef', 'wref'], 'ogv', 'mixrho']
Cdelw1 undel [['call'], 'orhnew', 'mixrho']
Cdelw1 undel [['call'], 'orhold', 'mixrho']
Cdelw1 undel [['redef', 'wref'], 'ocn', 'mixrho']
Cdelw1 undel [['redef', 'wref'], 'owk', 'mixrho']
Cdelw1 undel [['wref'], 'oa', 'mixrho']
Cdelw1 undel [['wref'], 'oaa', 'mixrho']
Cdelw1 undel [['array', 'wref', 'call'], 'orsm', 'pvmix3']
Cdelw1 undel [['call'], 'orhold', 'pvmix3']
Cdelw1 undel [['call'], 'orhnew', 'pvmix3']
Cdelw1 undel [['call'], 'orhold', 'pvmix5']
Cdelw1 undel [['call'], 'orhnew', 'pvmix5']
Cdelw1 undel [['call'], 'orhold', 'pvmix7']
Cdelw1 undel [['array', 'wref', 'call'], 'orsm', 'pvmix7']
Cdelw1 undel [['wref'], 'oqkl', 'mixrho']
Cgetarg...info...           structure ['mixrho', 'lat', 'slat']
Cgetarg...info...           structure ['mixrho', 'site', 'ssite']
Cgetarg...info...           structure ['mixrho', 'spec', 'sspec']
Cgetarg...info...           structure ['pvmix3', 'site', 'ssite']
Cgetarg...info...           structure ['pvmix3', 'spec', 'sspec']
Cgetarg...info...           structure ['pvmix5', 'site', 'ssite']
Cgetarg...info...           structure ['pvmix5', 'spec', 'sspec']
Cgetarg...info...           structure ['pvmix7', 'site', 'ssite']
Cgetarg...info...           structure ['pvmix7', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['pvmix3', 'm_struc_def']
Cgetarg...info...           use_to_add ['pvmix7', 'm_struc_def']
Cgetarg...info...           use_to_add ['mixrho', 'm_struc_def']
Cgetarg...info...           use_to_add ['pvmix5', 'm_struc_def']
      subroutine mixrho(ssite,sspec,slat,nsp,iter,mixmod,qval,elind,
     .  orhnew,orhold,k1,k2,k3,dmxp,smrnew,smrho)
       
       use m_struc_def  !Cgetarg

C- Mix old and new charge densities
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: rhgcmp ftlxp pvmix5 pvmix3 pvmix7 rhogkl
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: a nr rmt lmxl rsmv kmxv rg
Ci     Stored:    *
Ci     Passed to: rhgcmp ftlxp pvmix5 pvmix3 pvmix7 rhogkl
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: nabc ng ogv okv vol alat
Ci     Stored:    *
Ci     Passed to: rhgcmp
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   iter  :iteration number
Ci   mixmod:string containing iteration-dependent mixing (parmxp)
Ci   qval  :total valence charge, used to estimate Lindhard parameter
Ci         :elind if it is not supplied
Ci   elind :default Lindhard screening parameter
Ci   k1,k2,k3 dimensions smrho
Cio Inputs/Outputs
Cio   dmxp  :Default and other values of mixing parameters.
Cio         :The input defaults may be overridden by the string mixmod.
Cio         :Inputs (defaults)
Cio         :(1)  broy  :mixing scheme: 0->Anderson 1->Broyden
Cio         :(2)  beta  :linear mixing beta: nmix = (1-beta)nin + beta*nout
Cio         :(3)  wc    :Broyden weights for prior iterations
Cio         :(4,5) weights for n^+ + n^- and n^+ - n^- in spin pol cases
Cio                Avg, diff weighted by wt(1) and wt(2) in the mixing.
Cio                As special cases,
Cio                wt(2)=0 freezes the moments and mixes the charge only
Cio                wt(1)=0 freezes the charge and mixes the moments only
Cio         :(6)  mxsav : maximum number of prior iter to save on disk
Cio         :(7)  mmix  : maximum number of prior iter to mix
Cio         :(8)  nkill : kill mixing file after multiples nkill
Cio         :(9)  betv  : special potential mixing parameter
Cio         :(10) rmscst: not used here
Cio         :Outputs
Cio         :(11) rmsdel : rms (output-input) density
Cio         :(12) rmsl   : rms (output-input) local density
Cio         :(13) nmix   : number of prior iter actually mixed
Cio         :(14) broy   : actual mixing scheme used
Cio         :(15) beta   : actual mixing beta used
Cio         :(16-24) tj  : Anderson mixing coefficients
Cio         :(25)        : (spin polarized case only)
Cio                      : 1  if weight for n^+ + n^- is nonzero
Cio                      : 10 if weight for n^+ - n^- is nonzero
Cio                      : 11 if both nonzero
Cio   orhnew:On input, local parts of the density that generated the ham.
Cio   orhold:local parts of the density that the hamiltonian generated
Cio         :On output, the mixed density
Cio   smrho :On input, smooth density that generated the hamiltonian H
Cio         :On output mixed smooth density
Cio   smrnew:smooth density generated by the H(rhold,smrho)
Cio         :DESTROYED on output
Cl Local variables
Cl   locmix :governs linear transformation of local densities for mixing
Cl          :Note see Bugs, below
Cl          :0  This was the only mode in versions 6.11 and earlier.
Cl          :   The smooth part of the density and the spherical part of
Cl          :   the local densities are mixed in an Anderson or Broyden scheme.
Cl          :   In the remaining options,
Cl          :   local densities are rotated into
Cl          :     rho1,rho2 -> w(r)*(rho1+rho2), rho1-rho2
Cl          :   where w(r) is defined in pvmix9, mode controlled by k9.
Cl          :   For now, k9=10.
Cl          :   This has no effect for linear mixing, but affects the
Cl          :   Anderson or Broyden weights in the spherical part.
Cl          :1  Same as locmix=0, but for the transformation
Cl              CAUTION: not tested in a long time
Cl          :2  Local densities are further mapped onto a G_kL expansion.
Cl          :   whose coefficients are included in the Anderson/Broyden
Cl          :   mixing schemes, including the nonspherical densities.
Cl          :   the residual (what is left of rho beyond the G_kL expansion)
Cl          :   is linearly mixed.
Cl              CAUTION: not tested in a long time
Cl          :3  Similar to locmix=2, but l=0 part is of the G_kL expansion
Cl          :   is projected out.
Cl          :
Cl   rms2   : <(rout-rin)**2>, where rin and rout are the input and
Cl          : screened output densities represented as a vector,
Cl          : including smoothed and local parts.  rmsdel changes with
Cl          : the choice of transformation of local densities;see locmix
Cl          : above.
Cl   rmsdel : (linked to dmxp(11)).
Cl          : On input, rmsdel=<(rout-rin)**2>, defined as rms2, above.
Cl          : On output, rmsdel is overwritten by rms2.
Cl   rmsuns : RMS difference in output-input sm. density, unscreened
Cl   rms    : RMS difference in output-input sm. density, screened
Cl   sumo   : integral of smrho
Cl   sums   : integral of screened smoothed output rho
Cl   kmxs   : k-cutoff for G_kL expansion of screening charge projected
Cl          : into local densities; see Remarks
Cl   kmxr   : k-cutoff for G_kL expansion used to smooth local densities
Cl          : for mixing
Cl   ng0    : condensed number of G vectors. That rho(G) is hermitian
Cl          : is exploited to reduce ng to ng0
Cl   qcell  : cell charge
Cb Bugs
Cb   Mixing scheme needs some revision.  For now, uses locmix=3 and k9=10
Cb   Problem: G-kL expansion doesn't represent rim density well; not
Cb          : good for rho1+rho2.
Cb   Better:  1. Use simpler trans: (rho1+rho2)/100 or so, and (rho1-rho2)
Cb            2. For nonspherical parts, mix coefficient to multipole in
Cb               the Anderson/Broyden scheme; linearly mix the rest.
Cr Remarks
Cr   To estimate the self-consistent density a model (Lindhard) dielectric
Cr   function is used.  A smooth "total density" is assembled on a mesh
Cr   from the smooth density + local contributions; the latter are
Cr   smoothed into a G_kL expansion of rho1-rho2 and added to the mesh
Cr   density.  This "mesh total density" is assembled for the input and
Cr   output densities smrho and smrnew, and the screening charge
Cr   (eps^-1 - 1) (smrnew(tot) - smrho(tot)) is obtained using the
Cr   Lindhard approximation for eps^-1.  The screening charge is then:
Cr
Cr     (1)  added to smrnew
Cr
Cr     (2)  added back onto local densities as a G_kL expansion
Cr          to both the true and smooth local mixed densities.
Cr
Cr   Finally, the (screened) output density is mixed with the input density
Cr   using an Anderson or Broyden mixing scheme.  Both mesh part and local
Cr   parts are mixed, the latter consists of the true, smooth,
Cr   and core parts.
Cr   How the mixing proceeds depends on locmix; see Local variables above.
Cr
Cr   The steps in this routine follow this path (see numbers in comments)
Cr   1.   smooth rho+,rho-  ->  smrho, smmag
Cr   2.   FT smrnew(r)-smrho(r) -> cn = smrho(G)
Cr   3.   rho1-rho2 -> smoothed, FT; added to smrho(G)
Cr   4.   cn = screening (eps^-1 - 1) smrho(G)
Cr   5.   Add screening to smrnew; smrnew(G) -> smrnew(r)
Cr   6.   smrnew = sm. screened nout, smrho,smrnew -> smooth rho+,rho-
Cr        At this stage, smrho = input smooth rho+,rho-
Cr                       smrnew = screened output rho+,rho-
Cr   7.   Add screening to rho1, rho2; scale if locmix > 0
Cr   8.   smrho,smrnew(r) -> smrho,smrnew(G)
Cr   9.   Pack in,out densities into large matrix a;
Cr        Read prior iterations into a
Cr   10.  Linear mixing of local densities
Cr   11.  Make modified a when wt(1) or wt(2) restricts q or spin
Cr   12.  Mix input, output densities
Cr   13.  Inverse of step 11.
Cr   14.  Poke contents of a into smrho,orhold
Cr   15.  Undo scaling of local rho (step 7)
Cr
Cm MPI
Cm   master process handles the mix files and broadcasts. All processes
Cm   then mix.
Cu Updates
Cu   18 Dec 08 properly treats constrained cases w=1,0 and w=0,1
Cu   30 Jun 06 Bug fix (memory leak)
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   21 Nov 01 Linear transformation of local densities for mixing.
Cu             New modes locmix.
Cu   23 Oct 01 Local densities used for screening charge are more
Cr             accurately represented in G_kL expansion, k=0..kmxs
Cu   19 Oct 01 correctly uses elind override from mixing string
Cu   22 Mar 01 Bug fix for case mixing string is empty
Cu   22 Jun 00 spin polarized
Cu   30 May 00 Adapted from nfp mix_rho
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
#if MPI | MPIK
      include "mpif.h"
      integer numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
      integer iunit
      character*120 strn
      logical mlog,cmdopt
#endif
Cdelw1 do not change  [['call'], 'orhold', 'mixrho']
Cdelw1 do not change  [['call'], 'orhnew', 'mixrho']
      integer nsp,k1,k2,k3,orhold(3,1),orhnew(3,1),iter,procid,master
      character*(*) mixmod
Cgetarg       double precision ssite(1),sspec(1),slat(1)
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

      double precision qval,elind,dmxp(25)
      double complex smrnew(k1,k2,k3,nsp),smrho(k1,k2,k3,nsp)
C ... Local parameters
      integer stdo,stdl,lgunit,nbas,i,i1,i2,i3,ib,ipl,ipr,is,k0,k9,kmxv,
     .  lmxl,n1,n2,n3,ng,nglob,nlml,nr,nmixr,nmix,nda,mxsav,ifi,nlm0,
     .  kkk,nnnew,nnmix,ngabc(3),igetss,fopna,broy,nx,
     .  nkill,isw,naa,kmxr,kmxs,locmix,offx,off2
Cdelw1 do not change  [['redef', 'wref'], 'ocn', 'mixrho']
Cdelw1 do not change  [['undef', 'wref'], 'ogv', 'mixrho']
Cdelw1 do not change  [['undef', 'wref'], 'okv', 'mixrho']
Cdelw1 do not change  [['wref'], 'oa', 'mixrho']
Cdelw1 do not change  [['wref'], 'oaa', 'mixrho']
Cdelw1 do not change  [['redef', 'wref'], 'owk', 'mixrho']
Cdelw1 do not change  [['wref'], 'oqkl', 'mixrho']
Cdelw1       integer oips0,ng0,oco,ocn,ocg1,ocg2,ofkl,ogv,okv,orofi,orwgt,oa,
Cdelw1      .  oaa,owk,oqkl,owk1,owk2
       integer:: ng0 , ocn , ogv , okv , oa , oaa , owk , oqkl 
       integer ,allocatable :: iv_w_oips0(:)
       real(8) ,allocatable :: rv_w_oco(:)
       complex(8) ,allocatable :: zv_w_ocg1(:)
       complex(8) ,allocatable :: zv_w_ocg2(:)
       complex(8) ,allocatable :: zv_w_ofkl(:)
       real(8) ,allocatable :: rv_w_orofi(:)
       real(8) ,allocatable :: rv_w_orwgt(:)
       real(8) ,allocatable :: rv_w_owk1(:)
       real(8) ,allocatable :: rv_w_owk2(:)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      parameter (nlm0=49)
      double precision a,beta0,beta,dif,difx,difxu,fac,rms,rmt,
     .  sumo,summ,sums,top,vol,alat,tpiba,pi,dquns,rmsuns,ddot,q1,
     .  qin(2),qout(2),qscr(2),qmix(2),qcell,rms2,rms2f,wc,rsmv,qmx,
     .  dgets,wt(3),rmsdel,elinl,srfpi,xx
      double complex xxc
      logical parmxp
      character sout*80,fnam*8
C ... Heap
      integer w(1)
      common /w/ w

      complex(8),allocatable:: smrnewbk(:,:,:,:)
C     print *, '!!' ; call pshpr(40)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [okv ogv orhnew orhold ocn owk oa oaa oqkl]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (mixmod .eq. 'none') return
      call tcn('mixrho')

ccccccccccccccccccccccccccccccccccccc
c      allocate(smrnewbk(k1,k2,k3,nsp))
c      smrnewbk= (smrnew+smrho)/2d0
cccccccccccccccccccccccccccccccccccccc
      nnnew=0
c      open(1898,file='test2vvv_smrho')
      do  i3 = 1, k3
      do  i2 = 1, k2
      do  i1 = 1, k1
        if (dreal(smrho(i1,i2,i3,1)) .lt. 0) nnnew = nnnew+1
c        if( (i1==1.and.i2==1).or.
c     &      (i3==1.and.i2==1).or.
c     &      (i1==1.and.i3==1) ) then
c        write(1898,"(3i5,2d13.5)")i1,i2,i3,smrho(i1,i2,i3,1)
c        endif
      enddo
      enddo
      enddo
c      close(1898)
c      print *,'vvv'
      print *,'mixrho: sum smrho  init =',sum(smrho(:,:,:,1)),sum(abs(smrho(:,:,:,1))),nnnew
      nnnew=0
c      open(1898,file='test2vvv_smrnew')
      do  i3 = 1, k3
      do  i2 = 1, k2
      do  i1 = 1, k1
        if (dreal(smrnew(i1,i2,i3,1)) .lt. 0) nnnew = nnnew+1
c        if( (i1==1.and.i2==1).or.
c     &      (i3==1.and.i2==1).or.
c     &      (i1==1.and.i3==1) ) then
c          write(1898,"(3i5,2d13.5)")i1,i2,i3,smrnew(i1,i2,i3,1)-smrho(i1,i2,i3,1)
c        endif
      enddo
      enddo
      enddo
c      close(1898)
      print *,'mixrho: sum smrnew new  =',sum(smrnew(:,:,:,1)),sum(abs(smrnew(:,:,:,1))),nnnew
ccccccccccccccccccccccccccccccccccccc


#if MPI | MPIK
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
#else
      procid = 0
      master = 0
#endif

C --- Setup ---
      nbas = nglob('nbas')
      stdo = nglob('stdo')
      stdl = nglob('stdl')
      call getpr(ipr)
      ipl = 1
      nx = 0
      kkk = k1*k2*k3
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
ccccccccccccccccccccccccccccccccc
      locmix = 3 !original lmv7 by mark.
c      locmix = 3 
c      print *,'vvv mixrho locmix=',locmix
ccccccccccccccccccccccccccccccccc
      qmix(2) = 0
Cgetarg...info...           integer :: ngabc(3) , integer(8) :: slat%nabc 3
Cgetarg...info...           integer :: ng , integer(8) :: slat%ng 1
Cgetarg...info...           integer :: ogv , integer(8) :: slat%ogv 1
Cgetarg...info...           integer :: okv , integer(8) :: slat%okv 1
Cgetarg...info...           real(8) :: vol , real(8) :: slat%vol 1
Cgetarg       call upack('lat nabc ng ogv okv vol',slat,ngabc,ng,ogv,okv,vol)
       
       i_copy_size=size(slat%nabc) 
       call i8icopy(i_copy_size,slat%nabc,1,ngabc,1) 
       ng=slat%ng
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 314
Cdelw1 do not change ogv because of ['undef', 'wref']
       ogv=slat%ogv
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 315
Cdelw1 do not change okv because of ['undef', 'wref']
       okv=slat%okv
       vol=slat%vol

Cgetarg       alat = dgets('lat alat',slat)
       alat = slat%alat 

      elinl = elind
      kmxs = 3
      kmxr = 8

C --- Iteration-dependent mixing parameters ---
      broy  = dmxp(1)
      beta  = dmxp(2)
      wc    = dmxp(3)
      wt(1) = dmxp(4)
      wt(2) = dmxp(5)
      wt(3) = dmxp(4)
C     Flags parmxp that there are no extra elements to mix
      if (nx .eq. 0) wt(3) = -9
      mxsav = nint(dmxp(6))
      nmix  = nint(dmxp(7))
      nkill = nint(dmxp(8))
      fnam  = 'mixm'
      rmsdel = dmxp(11)
      rms2 = 0
      if (ipr .ge. 20) write(stdo,*) ' '
      if (.not. parmxp(iter,mixmod,len(mixmod),broy,nmix,wt,beta,elinl,
     .  fnam,wc,nkill,dmxp(9),rms2)) call rx(
     .  'MIXRHO: parse in parmxp failed')

C     In case parmxp doesn't touch wt, unset flag
      if (wt(3) .eq. -9) wt(3) = 0
      if (nmix .eq. 0) broy = 0
C ... Initial printout
C      call awrit7(' Mix density with beta=%;4d'//
C     .  '%?;n;, Lindhard E=%;3d;;'//
C     .  '%?;(n==1);, Broyden scheme using n=%,1i;%j;'//
C     .  '%?;(n==0);, Anderson scheme using n=%,1i;%j;',
C     .  ' ',80,stdo,beta,isw(elinl.ne.0),elinl,broy,nmix,
C     .  broy,nmix)
C      call rx('done')

C ... Interactively reset mix-parms block
      i1 = 0
      i2 = 0
      call query('mix-parms-block',2,i1)
      call query('block-iter',2,i2)
      if (i1 .gt. 0 .or. i2 .gt. 0) then
        if (i2 .le. 0) i2 = 1
        call parmx0(i1,i2,0d0)
        if (.not. parmxp(iter,mixmod,len(mixmod),broy,nmix,wt,beta,
     .    elinl,fnam,wc,nkill,dmxp(9),rms2))call rx(
     .    'MIXRHO: parse in parmxp failed')
      endif

      if (nsp .eq. 1) wt(2) = 0
      if (wt(1)**2+wt(2)**2+wt(3)**2 .eq. 0)
     .  call fexit(-1,111,' Exit -1 MIXRHO: '//
     .  'bad mixing weights w =%3:1;6d',wt)
      call dscal(3,1/dsqrt(wt(1)**2+wt(2)**2+wt(3)**2),wt,1)

C      wt(1) = 1 ; wt(2) = 1; print * , '!! wt=', wt(1:2)

      if (nmix .lt. 0) nmix = mxsav
      if (elinl .lt. 0d0) elinl = -(3*pi**2*qval/vol)**.66666d0*elinl
      if (elinl .lt. .01d0) elinl = 0

C ... Initial charges
      call rhoqm(smrho,k1,k2,k3,n1,n2,n3,nsp,vol,qin)
      call rhoqm(smrnew,k1,k2,k3,n1,n2,n3,nsp,vol,qout)

C --- Screen sm output rho; distribute screening rho over local rho ---
C ... 1. Overwrite smrho+, smrho- with smrho, smrho+ - smrho-
      if (nsp .eq. 2) then
C       call zprm3('smrho, spin1',0,smrho,k1,k2,k3)
C       call zprm3('smrho, spin2',0,smrho(1,1,1,2),k1,k2,k3)
C       call zprm3('smrnw, spin1',0,smrnew,k1,k2,k3)
C       call zprm3('smrnw, spin2',0,smrnew(1,1,1,2),k1,k2,k3)
        call dsumdf(kkk*2,1d0,smrho,0,1,smrho(1,1,1,2),0,1)
        call dsumdf(kkk*2,1d0,smrnew,0,1,smrnew(1,1,1,2),0,1)
C       call zprm3('smrho, charge',0,smrho,k1,k2,k3)
C       call zprm3('smrho, moment',0,smrho(1,1,1,2),k1,k2,k3)
C       call zprm3('smrnw, charge',0,smrnew,k1,k2,k3)
C       call zprm3('smrnw, moment',0,smrnew(1,1,1,2),k1,k2,k3)
      endif

C ... <smrout-smrhoin> and <(smrout-smrhoin)**2>, unscreened smrout
      dquns = 0
      rmsuns = 0
      do  12  i3 = 1, n3
      do  12  i2 = 1, n2
      do  12  i1 = 1, n1
      dquns  = dquns + dble(smrnew(i1,i2,i3,1)-smrho(i1,i2,i3,1))
   12 rmsuns = rmsuns + dble(smrnew(i1,i2,i3,1)-smrho(i1,i2,i3,1))**2

C ... 2. cn <- (smrnew - smrho)(G)
Cdelw1 do not change ocn because of ['redef', 'wref']
      call defcc(ocn, ng)
Cdelw1 do not change owk because of ['redef', 'wref']
      call defcc(owk, kkk)
Cdelw1 do not change owk because of ['redef', 'wref']
      call dcopy(kkk*2,smrnew,1,w(owk),1)
Cdelw1 do not change owk because of ['redef', 'wref']
      call dpadd(w(owk),smrho,1,kkk*2,-1d0)
Cdelw1 do not change owk because of ['redef', 'wref']
      call fftz3(w(owk),n1,n2,n3,k1,k2,k3,1,0,-1)
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 do not change owk because of ['redef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 416
Cdelw1 do not change okv because of ['undef', 'wref']
      call gvgetf(ng,1,w(okv),k1,k2,k3,w(owk),w(ocn))
Cdelw1 rlse name= owk old_list= ocn owk 
Cdelw1 rlse name= owk new_list= ocn 
Cdelw1 not deallocate owk because of [redef wref]
Cdelw1       call rlse(owk)
       call rlse(owk)


C ... 3. Sum of local densities on Fourier mesh (up+down)
Cdelw1       call defcc(ocg1,  -ng)
       allocate(zv_w_ocg1(abs(-ng)))
       if (-ng<0) zv_w_ocg1(:)=0.0d0

Cdelw1       call defcc(ocg2,  -ng)
       allocate(zv_w_ocg2(abs(-ng)))
       if (-ng<0) zv_w_ocg2(:)=0.0d0

Cdelw1       call rhgcmp(131,1,nbas,ssite,sspec,slat,orhold,kmxs,ng,w(ocg1))
       call rhgcmp ( 131 , 1 , nbas , ssite , sspec , slat , orhold 
     . , kmxs , ng , zv_w_ocg1 ) 

Cdelw1       call rhgcmp(131,1,nbas,ssite,sspec,slat,orhnew,kmxs,ng,w(ocg2))
       call rhgcmp ( 131 , 1 , nbas , ssite , sspec , slat , orhnew 
     . , kmxs , ng , zv_w_ocg2 ) 

C     call gvputf(ng,1,w(okv),k1,k2,k3,w(ocg1),smrnew)
C     call fftz3(smrnew,n1,n2,n3,k1,k2,k3,1,0,1)
C     call zprm3('input local densities (r)',0,smrnew,k1,k2,k3)

C ... Add output-input difference in gaussian rho to make total rho
C     call dpzero(w(ocn),ng*2)
      if (.true.) then
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1       call daxpy(2*ng,1d0,w(ocg2),1,w(ocn),1)
       call daxpy ( 2 * ng , 1d0 , zv_w_ocg2 , 1 , w ( ocn ) , 1 ) 

Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1       call daxpy(2*ng,-1d0,w(ocg1),1,w(ocn),1)
       call daxpy ( 2 * ng , - 1d0 , zv_w_ocg1 , 1 , w ( ocn ) , 1 ) 

      endif
Cdelw1 rlse name= ocg1 old_list= ocn ocg1 ocg2 
Cdelw1 rlse name= ocg1 new_list= ocn 
Cdelw1       call rlse(ocg1)
       if (allocated(zv_w_ocg2)) deallocate(zv_w_ocg2)
       if (allocated(zv_w_ocg1)) deallocate(zv_w_ocg1)


C ... 4. Make cn = screening charge = (eps^-1 - 1) (n_out - n_in)
      tpiba = 2*pi/alat
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 438
Cdelw1 do not change ogv because of ['undef', 'wref']
      call lindsc(2,ng,w(ogv),tpiba,elinl,w(ocn))

C ... 5. Add cn = screening charge into smrnew
Cdelw1 do not change owk because of ['redef', 'wref']
      call defcc(owk,kkk)
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 do not change owk because of ['redef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 442
Cdelw1 do not change okv because of ['undef', 'wref']
      call gvputf(ng,1,w(okv),k1,k2,k3,w(ocn),w(owk))
Cdelw1 do not change owk because of ['redef', 'wref']
      call fftz3(w(owk),n1,n2,n3,k1,k2,k3,1,0,1)
C     call zprm3('unscreened smrnew',0,smrnew,k1,k2,k3)
Cdelw1 do not change owk because of ['redef', 'wref']
      call dpadd(smrnew,w(owk),1,kkk*2,1d0)
Cdelw1 rlse name= owk old_list= ocn owk 
Cdelw1 rlse name= owk new_list= ocn 
Cdelw1 not deallocate owk because of [redef wref]
Cdelw1       call rlse(owk)
       call rlse(owk)

C     call zprm3('screened smrnw, charge',0,smrnew,k1,k2,k3)
C     call zprm3('screened smrnw, moment',0,smrnew(1,1,1,2),k1,k2,k3)

C ... RMS differences screened and input sm density; count neg points
C     summ
      sumo = 0d0
      sums = 0d0
      summ = 0d0
      top = 0d0
      rms = 0d0
      nnnew = 0
      nnmix = 0
      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1
        sumo = sumo + dble(smrho(i1,i2,i3,1))
        sums = sums + dble(smrnew(i1,i2,i3,1))
        dif = dble(smrnew(i1,i2,i3,1)-smrho(i1,i2,i3,1))
        rms = rms + dif*dif
        top = dmax1(top,dabs(dif))
        xxc = beta*smrnew(i1,i2,i3,1) + (1d0-beta)*smrho(i1,i2,i3,1)
        if (dble(smrnew(i1,i2,i3,1)) .lt. 0) nnnew = nnnew+1
        if (dble(xxc) .lt. 0) nnmix = nnmix+1
c       smrho(i1,i2,i3) = xxc
        summ = summ + dble(xxc)
      enddo
      enddo
      enddo
      fac = vol/(n1*n2*n3)
      rms = dsqrt(fac*rms/vol)
      rmsuns = dsqrt(fac*rmsuns/vol)
C     Defer printout until after file read prior iterations; see '100' below

C ... 6. Restore smrho+, smrho-
      if (nsp .eq. 2) then
        call dsumdf(kkk*2,.5d0,smrho,0,1,smrho(1,1,1,2),0,1)
        call dsumdf(kkk*2,.5d0,smrnew,0,1,smrnew(1,1,1,2),0,1)
      endif
C     call zprm3('smrnew',0,smrnew,k1,k2,k3)

C --- 7. Project cn = screening density into local densities ---
      k0 = 20
Cdelw1       call defcc(ofkl,(k0+1)*nlm0*nbas)
       allocate(zv_w_ofkl((k0+1)*nlm0*nbas))
Ctakao_ZeroClear_NotRequiered        if ((k0+1)*nlm0*nbas<0) zv_w_ofkl(:)=0.0d0

Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 490
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 490
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1       call ftlxp(nbas,ssite,sspec,alat,ng,w(ogv),w(ocn),k0,nlm0,w(ofkl))
       call ftlxp ( nbas , ssite , sspec , alat , ng , w ( ogv ) , w 
     . ( ocn ) , k0 , nlm0 , zv_w_ofkl ) 

      do  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg         call upack('spec a nr rmt',sspec,is,a,nr,rmt,0)
         
         a=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           real(8) :: rsmv , real(8) :: sspec%rsmv 1
Cgetarg...info...           integer :: kmxv , integer(8) :: sspec%kmxv 1
Cgetarg         call upack('spec lmxl rsmv kmxv',sspec,is,lmxl,rsmv,kmxv,0)
         
         lmxl=sspec(is)%lmxl
         rsmv=sspec(is)%rsmv
         kmxv=sspec(is)%kmxv

        if (lmxl .lt. 0) cycle
        nlml = (lmxl+1)**2
        if (nlml.gt.nlm0) call rxi('mixrho: nlml.gt.nlm0, need',nlml)
Cdelw1         call defrr(orofi, nr)
         allocate(rv_w_orofi(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) rv_w_orofi(:)=0.0d0

Cdelw1         call defrr(orwgt, nr)
         allocate(rv_w_orwgt(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) rv_w_orwgt(:)=0.0d0

Cdelw1         call radmsh(rmt,a,nr,w(orofi))
         call radmsh ( rmt , a , nr , rv_w_orofi ) 

Cdelw1         call radwgt(rmt,a,nr,w(orwgt))
         call radwgt ( rmt , a , nr , rv_w_orwgt ) 

C       Overwrite rho+, rho- with rho, rho+ - rho-
Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1 3 is digit, dropped
        call splrho(0,nsp,nr,nlml,
     .    w(orhnew(1,ib)),w(orhnew(2,ib)),w(orhnew(3,ib)))
C   ... Add site-projected screening density to rhn1,rhn2
Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1         call pkl2ro(110,ib,rsmv,kmxv,nr,nlml,1,w(orofi),w(orwgt),
Cdelw1      .    k0,nlm0,w(ofkl),w,w(orhnew(1,ib)),w(orhnew(2,ib)),qmx)
         call pkl2ro ( 110 , ib , rsmv , kmxv , nr , nlml , 1 , rv_w_orofi 
     .   , rv_w_orwgt , k0 , nlm0 , zv_w_ofkl , w , w ( orhnew ( 1 , ib 
     .   ) ) , w ( orhnew ( 2 , ib ) ) , qmx ) 

C       Restore rho+, rho-
Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1 3 is digit, dropped
        call splrho(1,nsp,nr,nlml,
     .    w(orhnew(1,ib)),w(orhnew(2,ib)),w(orhnew(3,ib)))

C   ...  Always work with rho1+rho2, rho1-rho2
C        print *, 'ib=',ib
C        call prrmsh('rhold1',w(orofi),w(orhold(1,1)),nr,nr,nlml*nsp)
C        call prrmsh('rhold2',w(orofi),w(orhold(2,1)),nr,nr,nlml*nsp)

Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1         call pvmix9(1,0,nr,nlml*nsp,0,0d0,w(orofi),w(orhold(1,ib)),
Cdelw1      .    w(orhold(2,ib)))
         call pvmix9 ( 1 , 0 , nr , nlml * nsp , 0 , 0d0 , rv_w_orofi 
     .   , w ( orhold ( 1 , ib ) ) , w ( orhold ( 2 , ib ) ) ) 

Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1         call pvmix9(1,0,nr,nlml*nsp,0,0d0,w(orofi),w(orhnew(1,ib)),
Cdelw1      .    w(orhnew(2,ib)))
         call pvmix9 ( 1 , 0 , nr , nlml * nsp , 0 , 0d0 , rv_w_orofi 
     .   , w ( orhnew ( 1 , ib ) ) , w ( orhnew ( 2 , ib ) ) ) 


C       debugging
C       call prrmsh('1 post t',w(orofi),w(orhnew(1,1)),nr,nr,nlml*nsp)
C       call prrmsh('2 post t',w(orofi),w(orhnew(2,1)),nr,nr,nlml*nsp)
C       Check inverse op restores original
C       call pvmix9(k9,-1,nr,nlml*nsp,0,rf,w(orofi),w(orhold(1,ib)),
C    .    w(orhold(2,ib)))
C       call pvmix9(k9,-1,nr,nlml*nsp,0,rf,w(orofi),w(orhnew(1,ib)),
C    .    w(orhnew(2,ib)))
C       call prrmsh('1 rest',w(orofi),w(orhold(1,1)),nr,nr,nlml*nsp)
C       call prrmsh('2 rest',w(orofi),w(orhold(2,1)),nr,nr,nlml*nsp)
C       stop

Cdelw1 rlse name= orofi old_list= ocn ofkl orofi orwgt 
Cdelw1 rlse name= orofi new_list= ocn ofkl 
Cdelw1         call rlse(orofi)
         if (allocated(rv_w_orwgt)) deallocate(rv_w_orwgt)
         if (allocated(rv_w_orofi)) deallocate(rv_w_orofi)

      enddo                     ! Loop over sites
Cdelw1 rlse name= ocn old_list= ocn ofkl 
Cdelw1 rlse name= ocn new_list= (None)
Cdelw1 not deallocate ocn because of [redef wref]
Cdelw1       call rlse(ocn)
       if (allocated(zv_w_ofkl)) deallocate(zv_w_ofkl)
       call rlse(ocn)


C --- Fancy mixing of smoothed + some representation of local rho ---
C ... Count number of elts from local densities for fancy mixing scheme
      nda = 0
      do  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg         nr = igetss('spec nr',is,sspec)
         nr = int(sspec(is)%nr) 

Cgetarg         lmxl = igetss('spec lmxl',is,sspec)
         lmxl = int(sspec(is)%lmxl) 

        nlml = (lmxl+1)**2
        if (lmxl .eq. -1) goto 20
C       include spherical part of local densities only
        if (locmix .le. 1) then
          nda = nda + 2*nr
C       also include nonspherical part, G_kL expansion
        elseif (locmix .ge. 2) then
          nda = nda + 2*(kmxr+1)*nlml
C         and spherical part of rho1-rho2
          if (locmix .eq. 3) nda = nda + 2*nr
        else
          call rxi('bad locmix,',locmix)
        endif
   20   continue
      enddo
C ... Find and include number of smooth mesh points.
Cdelw1       call defi(oips0,ng)
       allocate(iv_w_oips0(ng))
Ctakao_ZeroClear_NotRequiered        if (ng<0) iv_w_oips0(:)=0

Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 586
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 586
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1       call lgstar(0,ng,1,w(ogv),ng0,w(oips0),w)
       call lgstar ( 0 , ng , 1 , w ( ogv ) , ng0 , iv_w_oips0 , w ) 

C     if (ng .ne. (ng0-1)*2+1) call rx('bug in mixrho')
      nda = nda + 2*ng0

C ... 8. FFT smooth densities -> rho(G) for mixing
Cdelw1 w_varlistall_append_uniq ERROR, define again with another type: ocn complex(8) and real(8)
Cdelw1 do not change ocn because of ['redef', 'wref']
      call defrr(ocn,ng*2*nsp)
Cdelw1       call defrr(oco,ng*2*nsp)
       allocate(rv_w_oco(ng*2*nsp))
Ctakao_ZeroClear_NotRequiered        if (ng*2*nsp<0) rv_w_oco(:)=0.0d0

Cdelw1 w_varlistall_append_uniq ERROR, define again with another type: owk complex(8) and real(8)
Cdelw1 do not change owk because of ['redef', 'wref']
      call defrr(owk,kkk*2*nsp)
Cdelw1 do not change owk because of ['redef', 'wref']
      call dpcopy(smrnew,w(owk),1,kkk*2*nsp,1d0)
Cdelw1 do not change owk because of ['redef', 'wref']
      call fftz3(w(owk),n1,n2,n3,k1,k2,k3,nsp,0,-1)
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 do not change owk because of ['redef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 596
Cdelw1 do not change okv because of ['undef', 'wref']
      call gvgetf(ng,nsp,w(okv),k1,k2,k3,w(owk),w(ocn))
C     call zprm('cn-u',2,w(ocn),ng,ng,nsp)
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 598
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 598
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1       call lgstar(1,ng,nsp,w(ogv),i,w(oips0),w(ocn))   !ocn new
       call lgstar ( 1 , ng , nsp , w ( ogv ) , i , iv_w_oips0 , w ( 
     . ocn ) ) 

C     call zprm('cn-c',2,w(ocn),ng,ng0,nsp)
C     call lgstar(2,ng,nsp,w(ogv),i,w(oips0),w(ocn))
C     call zprm('cn-u',2,w(ocn),ng,ng,nsp)
Cdelw1 do not change owk because of ['redef', 'wref']
      call dpcopy(smrho,w(owk),1,kkk*2*nsp,1d0)
Cdelw1 do not change owk because of ['redef', 'wref']
      call fftz3(w(owk),n1,n2,n3,k1,k2,k3,nsp,0,-1)
Cdelw1 do not change owk because of ['redef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 604
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1 do not change owk because of ['redef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 604
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1       call gvgetf(ng,nsp,w(okv),k1,k2,k3,w(owk),w(oco))
       call gvgetf ( ng , nsp , w ( okv ) , k1 , k2 , k3 , w ( owk ) 
     . , rv_w_oco ) 

C     call zprm('co-u',2,w(oco),ng,ng,nsp)
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 606
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 606
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1       call lgstar(1,ng,nsp,w(ogv),i,w(oips0),w(oco))   !oco init
       call lgstar ( 1 , ng , nsp , w ( ogv ) , i , iv_w_oips0 , rv_w_oco 
     . ) 

ccccccccccc
c      print *,'vvv lgstar i=',i
ccccccccccc
C     call zprm('co-c',2,w(oco),ng,ng0,nsp)
Cdelw1 rlse name= owk old_list= oips0 ocn oco owk 
Cdelw1 rlse name= owk new_list= oips0 ocn oco 
Cdelw1 not deallocate owk because of [redef wref]
Cdelw1       call rlse(owk)
       call rlse(owk)

Cdelw1 do not change oa because of ['wref']
      call defrr(oa,-nda*nsp*(mxsav+2)*2)

C ... Allocate memory for coffs to G_kL expansion of rho
Cdelw1 warning(1) , probably  oqkl  is not defined yet at linenumber= 615
Cdelw1 do not change oqkl because of ['wref']
      oqkl = 1
Cdelw1 do not change oqkl because of ['wref']
      if (locmix .ge. 2) call defcc(oqkl,-(kmxr+1)*nlm0*nsp*4*nbas)

C --- 9. Read prior iterations from disk; update with current iter ---
      if (procid .eq. master) then
        ifi = fopna(fnam,-1,4)
        rewind ifi
      endif
#if MPI | MPIK
      call MPI_BCAST(ifi,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit3(' mixrho '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' bcast ifi = %i',' ',
     .      256,lgunit(3),procid,numprocs,ifi)
      endif
#endif

C     NB: For now, nx must be zero
      k9 = 10
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 do not change oa because of ['wref']
Cdelw1 do not change oqkl because of ['wref']
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 do not change oa because of ['wref']
Cdelw1 do not change oqkl because of ['wref']
Cdelw1       call pvmix5(nmix,mxsav,fnam,ifi,.true.,rmsdel,locmix,k9,nbas,kmxr,
Cdelw1      .  nlm0,ssite,sspec,nsp,orhold,orhnew,w(oco),w(ocn),ng,ng0,0,nda,
Cdelw1      .  w(oa),w(oqkl),rms2,nmixr)
       call pvmix5 ( nmix , mxsav , fnam , ifi , .true. , rmsdel , locmix 
     . , k9 , nbas , kmxr , nlm0 , ssite , sspec , nsp , orhold , orhnew 
     . , rv_w_oco , w ( ocn ) , ng , ng0 , 0 , nda , w ( oa ) , w ( 
     . oqkl ) , rms2 , nmixr ) 


C     debugging check
C      ib = 1
C      is = igetss('site spec',ib,ssite)
C      call upack('spec a nr rmt',sspec,is,a,nr,rmt,0)
C      call upack('spec lmxl rsmv kmxv',sspec,is,lmxl,rsmv,kmxv,0)
C      nlml = (lmxl+1)**2
C      call defrr(orofi, nr)
C      call defrr(orwgt, nr)
C      call radmsh(rmt,a,nr,w(orofi))
C      call radwgt(rmt,a,nr,w(orwgt))
C      call upack2('spec rg',sspec,is,rg)
C      call pkl2ro(1001,ib,rg,kmxr,nr,nlml,nsp,w(orofi),w(orwgt),
C     .  kmxr,nlm0,w,w(oqkl),w(orhold(1,ib)),w,qmx)
C      call prrmsh('rhots1',w(orofi),w(orhold(1,ib)),nr,nr,nlml*nsp)
C      stop

      rmsdel = rms2
      nmix = min(nmix,nmixr)
C ... Write this and prior iterations onto disk
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 do not change oa because of ['wref']
Cdelw1 do not change oqkl because of ['wref']
Cdelw1 do not change ocn because of ['redef', 'wref']
Cdelw1 do not change oa because of ['wref']
Cdelw1 do not change oqkl because of ['wref']
Cdelw1       call pvmix5(nmix,mxsav,fnam,-ifi,.true.,rmsdel,locmix,k9,nbas,
Cdelw1      .  kmxr,nlm0,ssite,sspec,nsp,orhold,orhnew,w(oco),w(ocn),ng,ng0,0,
Cdelw1      .  nda,w(oa),w(oqkl),rms2,nmixr)
       call pvmix5 ( nmix , mxsav , fnam , - ifi , .true. , rmsdel , 
     . locmix , k9 , nbas , kmxr , nlm0 , ssite , sspec , nsp , orhold 
     . , orhnew , rv_w_oco , w ( ocn ) , ng , ng0 , 0 , nda , w ( oa 
     . ) , w ( oqkl ) , rms2 , nmixr ) 

      if (nkill .lt. 0) then
        if (ipr .ge. 20) then
          sout = ' ...  deleting file '//fnam
          call awrit1('%a : nkill=%i',sout,80,-lgunit(1),-nkill)
        endif
        call dfclos(ifi)
      else
        call fclose(ifi)
      endif

C --- Printout smooth density mixing data ---
      if (nnnew+nnmix.gt.0 .and. ipr.ge.20) call awrit2(
     .  ' mixrho: (warning) scr. and lin-mixed densities had'//
     .  ' %i and %i negative points',' ',80,stdo,nnnew,nnmix)
      call rhoqm(smrnew,k1,k2,k3,n1,n2,n3,nsp,vol,qscr)
C     if (ipr .gt. 30) write(stdo,100)
C    .  sumo*fac,sumu*fac,sums*fac,rms,summ*fac
      if (ipr .gt. 30) then
        write(stdo,100)
     .    qin(1),qout(1),qscr(1),rms,beta*qscr(1)+(1-beta)*qin(1)
        if (nsp .eq. 2) write(stdo,101)
     .    qin(2),qout(2),beta*qscr(2)+(1-beta)*qin(2)
      endif
  100 format(' charges:',7x,'old',11x,'new',9x,'screened',
     .  6x,'rms diff',7x,'lin mix'/' smooth ',5f14.6)
  101 format(' mmom   ',2f14.6,28x,f14.6)

C --- 10. Linear mixing of local densities  ---
Cdelw1 do not change oqkl because of ['wref']
      call pvmix3(ssite,sspec,nbas,nsp,beta,locmix,wt,kmxr,nlm0,
     .  k9,w(oqkl),orhold,orhnew,difx)
      difxu = difx

C --- 11. Spin polarized case: separate weighting for spin channels ---
      if (nsp .eq. 2 .or. nx .gt. 0) then
        if (wt(2).ne.0 .and. nsp.eq.1)
     .    call rx1('mixrho wt(2)=%d invalid when nsp=1',wt(2))
C       Map to ASA call:
C             dim. of a     extra    #-mix   #-elts/spin  dim of aa
C       ASA    nda           nx      naa      npq            naa
C       FP     nda*nsp       0       naa      nda            naa
        if (nx .gt. 0) call rx('mixrho not ready for nx=0')
        naa = 0
        if (wt(1) .ne. 0) naa = naa+nda
        if (wt(2) .ne. 0) naa = naa+nda
        if (wt(3) .ne. 0) naa = naa+nx
        offx = 0                !offset to extra elements (none now)
        off2 = (nsp-1)*nda      !offset to spin down part of a
Cdelw1 do not change oaa because of ['wref']
        call defdr(oaa,-naa*(mxsav+2)*2)
Cdelw1 do not change oa because of ['wref']
Cdelw1 do not change oaa because of ['wref']
        call pqsclf(0,nda*nsp,nda,offx,off2,nx,naa,mxsav,wt,
     .    w(oa),w(oaa),rms2f)
        if (wt(1)*wt(2) .ne. 0) then
C          rms2f = rms2f * 2
        else
          call info2(20,0,0,' Constrained spin mixing, wt =%2:1;5d.'//
     .      '  Constrained rms DQ=%1,3;3e',wt,rms2f)
          rmsdel = rms2f
        endif
      else
Cdelw1 do not change oa because of ['wref']
Cdelw1 do not change oaa because of ['wref']
        oaa = oa
        naa = nda
      endif

C --- 12. Mix the soup of densities ---
      beta0 = beta
Cdelw1 do not change oaa because of ['wref']
      call pvmix6(broy,nmix,nmixr,mxsav,beta,wc,naa,w(oaa))

C ... 13. Restore matrix a to rho+, rho-
      if (nsp .eq. 2 .or. nx .gt. 0) then
Cdelw1 do not change oa because of ['wref']
Cdelw1 do not change oaa because of ['wref']
        call pqsclb(nda*nsp,nda,offx,off2,nx,naa,mxsav,wt,w(oa),w(oaa))
Cdelw1 do not change oa because of ['wref']
        call pqsclc(nda*nsp,nda*nsp,nx,mxsav,w(oa))
Cdelw1 rlse name= oaa old_list= oips0 ocn oco oa oqkl oaa 
Cdelw1 rlse name= oaa new_list= oips0 ocn oco oa oqkl 
Cdelw1 not deallocate oaa because of [wref]
Cdelw1         call rlse(oaa)
         call rlse(oaa)

      endif

C      nlml = 16
C      call prrmsh('rho1, ante 7',w(orofi),w(orhold(1,1)),nr,nr,nlml*nsp)
C      call prrmsh('rho2, ante 7',w(orofi),w(orhold(2,1)),nr,nr,nlml*nsp)

C ... 14. Poke mixed smooth and local densities into smrho,rhoold
Cdelw1 do not change owk because of ['redef', 'wref']
      call defcc(owk,kkk)
Cdelw1 do not change owk because of ['redef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 741
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 741
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 do not change oa because of ['wref']
Cdelw1 do not change oqkl because of ['wref']
Cdelw1 do not change owk because of ['redef', 'wref']
Cdelw1 warning(1) , probably  okv  is not defined yet at linenumber= 741
Cdelw1 do not change okv because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ogv  is not defined yet at linenumber= 741
Cdelw1 do not change ogv because of ['undef', 'wref']
Cdelw1 do not change oa because of ['wref']
Cdelw1 do not change oqkl because of ['wref']
Cdelw1       call pvmix7(ssite,sspec,nbas,nsp,nda,w(oa),n1,n2,n3,k1,k2,k3,
Cdelw1      .  locmix,wt,k9,kmxr,nlm0,w(oqkl),ng,ng0,w(okv),w(oips0),w(ogv),
Cdelw1      .  w(oco),w(owk),orhold,smrho)
       call pvmix7 ( ssite , sspec , nbas , nsp , nda , w ( oa ) , n1 
     . , n2 , n3 , k1 , k2 , k3 , locmix , wt , k9 , kmxr , nlm0 , w 
     . ( oqkl ) , ng , ng0 , w ( okv ) , iv_w_oips0 , w ( ogv ) , rv_w_oco 
     . , w ( owk ) , orhold , smrho ) 


Cdelw1 rlse name= oa old_list= oips0 ocn oco oa oqkl owk 
Cdelw1 rlse name= oa new_list= oips0 ocn oco 
Cdelw1 not deallocate owk because of [redef wref]
Cdelw1 not deallocate oqkl because of [wref]
Cdelw1 not deallocate oa because of [wref]
Cdelw1       call rlse(oa)
       call rlse(oa)

Cdelw1 rlse name= ocn old_list= oips0 ocn oco 
Cdelw1 rlse name= ocn new_list= oips0 
Cdelw1 not deallocate ocn because of [redef wref]
Cdelw1       call rlse(ocn)
       if (allocated(rv_w_oco)) deallocate(rv_w_oco)
       call rlse(ocn)

Cdelw1 rlse name= oips0 old_list= oips0 
Cdelw1 rlse name= oips0 new_list= (None)
Cdelw1       call rlse(oips0)
       if (allocated(iv_w_oips0)) deallocate(iv_w_oips0)


C      nlml = 16
C      call prrmsh('rho1 post 7',w(orofi),w(orhold(1,1)),nr,nr,nlml*nsp)
C      call prrmsh('rho2 post 7',w(orofi),w(orhold(2,1)),nr,nr,nlml*nsp)
C      stop

C      Debugging
C      call dsumdf(kkk*2,1d0,smrho,0,1,smrho(1,1,1,2),0,1)
C      call zprm3('mixed sm rho, charge',0,smrho,k1,k2,k3)
C      call zprm3('mixed sm rho, spin',0,smrho(1,1,1,2),k1,k2,k3)
C      call dsumdf(kkk*2,.5d0,smrho,0,1,smrho(1,1,1,2),0,1)
C      call zprm3('mixed sm rho, spin1',0,smrho,k1,k2,k3)
C      call zprm3('mixed sm rho, spin2',0,smrho(1,1,1,2),k1,k2,k3)

C ... 15. Restore local densities: rho+ +/- rho-  -> rho+, rho-'
      call dpzero(qmix,2)
      do  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg         call upack('spec a nr rmt',sspec,is,a,nr,rmt,0)
         
         a=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt

Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           real(8) :: rsmv , real(8) :: sspec%rsmv 1
Cgetarg...info...           integer :: kmxv , integer(8) :: sspec%kmxv 1
Cgetarg         call upack('spec lmxl rsmv kmxv',sspec,is,lmxl,rsmv,kmxv,0)
         
         lmxl=sspec(is)%lmxl
         rsmv=sspec(is)%rsmv
         kmxv=sspec(is)%kmxv

        if (lmxl .lt. 0) cycle
        nlml = (lmxl+1)**2
        if (nlml.gt.nlm0) call rxi('mixrho: nlml.gt.nlm0, need',nlml)
Cdelw1         call defrr(orofi, nr)
         allocate(rv_w_orofi(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) rv_w_orofi(:)=0.0d0

Cdelw1         call defrr(orwgt, nr)
         allocate(rv_w_orwgt(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) rv_w_orwgt(:)=0.0d0

Cdelw1         call radmsh(rmt,a,nr,w(orofi))
         call radmsh ( rmt , a , nr , rv_w_orofi ) 

Cdelw1         call radwgt(rmt,a,nr,w(orwgt))
         call radwgt ( rmt , a , nr , rv_w_orwgt ) 


Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1         call pvmix9(1,-1,nr,nlml*nsp,0,0d0,w(orofi),w(orhold(1,ib)),
Cdelw1      .    w(orhold(2,ib)))
         call pvmix9 ( 1 , - 1 , nr , nlml * nsp , 0 , 0d0 , rv_w_orofi 
     .   , w ( orhold ( 1 , ib ) ) , w ( orhold ( 2 , ib ) ) ) 

Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1         call pvmix9(1,-1,nr,nlml*nsp,0,0d0,w(orofi),w(orhnew(1,ib)),
Cdelw1      .    w(orhnew(2,ib)))
         call pvmix9 ( 1 , - 1 , nr , nlml * nsp , 0 , 0d0 , rv_w_orofi 
     .   , w ( orhnew ( 1 , ib ) ) , w ( orhnew ( 2 , ib ) ) ) 


C       debugging
C        print *, 'ib=',ib
C        call prrmsh('1 final',w(orofi),w(orhold(1,1)),nr,nr,nlml*nsp)
C        call prrmsh('2 final',w(orofi),w(orhold(2,1)),nr,nr,nlml*nsp)

C   ... Add net local charge to qmix
        do  i = 1, nsp
          off2 = 1 + nr*nlml*(i-1)
Cdelw1           call defrr(owk1,nr)
           allocate(rv_w_owk1(nr))
Ctakao_ZeroClear_NotRequiered            if (nr<0) rv_w_owk1(:)=0.0d0

Cdelw1           call defrr(owk2,nr)
           allocate(rv_w_owk2(nr))
Ctakao_ZeroClear_NotRequiered            if (nr<0) rv_w_owk2(:)=0.0d0

Cdelw1 1 is digit, dropped
Cdelw1           call dpscop(w(orhold(1,ib)),w(owk1),nr,off2,1,1d0)
           call dpscop ( w ( orhold ( 1 , ib ) ) , rv_w_owk1 , nr , off2 
     .     , 1 , 1d0 ) 

Cdelw1 2 is digit, dropped
Cdelw1           call dpscop(w(orhold(2,ib)),w(owk2),nr,off2,1,1d0)
           call dpscop ( w ( orhold ( 2 , ib ) ) , rv_w_owk2 , nr , off2 
     .     , 1 , 1d0 ) 

Cdelw1           qmx = srfpi*(ddot(nr,w(owk1),1,w(orwgt),1) -
Cdelw1      .                 ddot(nr,w(owk2),1,w(orwgt),1))
           qmx = srfpi * ( ddot ( nr , rv_w_owk1 , 1 , rv_w_orwgt , 1 ) 
     .     - ddot ( nr , rv_w_owk2 , 1 , rv_w_orwgt , 1 ) ) 

C         print *, 'spin, qmx',i,qmx
          if (i .eq. 2) qmix(2) = qmix(2) + qmx-q1
          q1 = qmx
          qmix(1) = qmix(1) + qmx
Cdelw1 rlse name= owk1 old_list= orofi orwgt owk1 owk2 
Cdelw1 rlse name= owk1 new_list= orofi orwgt 
Cdelw1           call rlse(owk1)
           if (allocated(rv_w_owk2)) deallocate(rv_w_owk2)
           if (allocated(rv_w_owk1)) deallocate(rv_w_owk1)

        enddo

Cdelw1 rlse name= orofi old_list= orofi orwgt 
Cdelw1 rlse name= orofi new_list= (None)
Cdelw1         call rlse(orofi)
         if (allocated(rv_w_orwgt)) deallocate(rv_w_orwgt)
         if (allocated(rv_w_orofi)) deallocate(rv_w_orofi)

      enddo                     !Loop over sites

C ... Force density positive
C     call rhopos(smrho,k1,k2,k3,n1,n2,n3)

C ... Net interstitial charge and cell charge
      summ = 0
      do  60  i = 1, nsp
      do  60  i3 = 1, n3
      do  60  i2 = 1, n2
      do  60  i1 = 1, n1
C       smrho(i1,i2,i3,i) = smrho(i1,i2,i3,i) + .001
        summ = summ + dble(smrho(i1,i2,i3,i))
   60 continue
C     Net system charge
      qcell  = qval - summ*fac - qmix(1)

C     Add constant to smrho to force charge neutrality
      xx = qcell/(vol*nsp)
      summ = 0

cccccccccccccccccccccccccc
c       print *,'vvv no xx test'
c       xx=0d0
ccccccccccccccccccccccccc
      do  64  i  = 1, nsp
      do  64  i3 = 1, n3
      do  64  i2 = 1, n2
      do  64  i1 = 1, n1
        smrho(i1,i2,i3,i) = smrho(i1,i2,i3,i) + xx
C       summ = summ + dble(smrho(i1,i2,i3,i))
   64 continue
C     qcell  = qval - summ*fac - qmix(1)
C     print *, 'qcell after shift',qcell; stop

C --- Clean up ---
      if (nmix .lt. 0) dmxp(7) = -nmix
      if (beta0 .ne. beta) dmxp(2) = beta
      dmxp(11) = rmsdel
      dmxp(12) = difx
      dmxp(13) = iabs(nmix)
      dmxp(14) = broy
      dmxp(15) = beta
      if (broy .eq. 1) dmxp(15) = 1
      dmxp(25) = 0
      if (wt(1) .ne. 0 .or. wt(2) .ne. 0) dmxp(25) = 1
      if (wt(3) .ne. 0) dmxp(25) = dmxp(25) + 10

C ... Printout
      if (ipr .ge. 10 .and. abs(qcell) .gt. 1d-6) write(stdo,
     .  '('' add q='',f10.6,'' to preserve neutrality'')') qcell

      if (ipr .ge. 30) then
        write(stdo,450) ' unscreened ',rmsuns,difxu
        if (elinl .ne. 0) write(stdo,450) '   screened ',rms,difx,rmsdel
  450   format(a,'rms difference:  smooth',f10.6,'   local',f10.6:
     .    '   tot',f10.6)
      elseif (ipr .ge. 20) then
        call awrit3(' rms smooth dq=%;3g  max local dq=%;3g  dq=%;3g',
     .    ' ',80,stdo,rms,difx,rmsdel)
      endif
      if (ipl .gt. 0 .and. ipr .gt. 0) then
        call awrit6('%xfp %?;n;elind %;4g  ;%j;'//
     .    '%?#(n==1)#Broyden n %1i wc %;3g#%2j#%-3j'//
     .    '%?#(n==0)#Anderson n %i beta %j%;3g#%3j#',
     .    sout,80,0,isw(elinl.ne.0),elinl,broy,nmix,wc,beta)
        call awrit3('%a  sm-dq %;3g  mx loc %;3g  dq %;3g',sout,80,
     .    -stdl,rms,difx,rmsdel)
      endif

C      call zprm3('exit sm rho, spin1',0,smrho(1,1,1,1),k1,k2,k3)
C      call zprm3('exit sm rho, spin2',0,smrho(1,1,1,2),k1,k2,k3)

      call tcx('mixrho')
C     print *, '!!' ; call poppr
C     stop 'for now'
      nnnew = 0
      nnmix = 0
      do  i3 = 1, n3
      do  i2 = 1, n2
      do  i1 = 1, n1
        sumo = sumo + dble(smrho(i1,i2,i3,1))
        sums = sums + dble(smrnew(i1,i2,i3,1))
        dif = dble(smrnew(i1,i2,i3,1)-smrho(i1,i2,i3,1))
        rms = rms + dif*dif
        top = dmax1(top,dabs(dif))
        xxc = beta*smrnew(i1,i2,i3,1) + (1d0-beta)*smrho(i1,i2,i3,1)
        if (dble(smrnew(i1,i2,i3,1)) .lt. 0) nnnew = nnnew+1
        if (dble(xxc) .lt. 0) nnmix = nnmix+1
c       smrho(i1,i2,i3) = xxc
        summ = summ + dble(xxc)
      enddo
      enddo
      enddo
cccccccccccccccccccccccccccccccccccc
      nnnew=0
c      open(1898,file='test3vvv')
      do  i3 = 1, k3
      do  i2 = 1, k2
      do  i1 = 1, k1
c        write(1898,"(3i5,2d13.5)")i1,i2,i3,smrho(i1,i2,i3,1)
        if (dreal(smrho(i1,i2,i3,1)) .lt. 0) nnnew = nnnew+1
      enddo
      enddo
      enddo
c      close(1898)
      print *,'mixrho: sum smrho output=',sum(smrho(:,:,:,1)),sum(abs(smrho(:,:,:,1))),nnnew
c      nnnew=0
c      do  i3 = 1, k3
c      do  i2 = 1, k2
c      do  i1 = 1, k1
c        if(dreal(smrho(i1,i2,i3,1))<1d-10) smrho(i1,i2,i3,1)=1d-10
c        if (dreal(smrho(i1,i2,i3,1)) .lt. 0) nnnew = nnnew+1
c      enddo
c      enddo
c      enddo
c      print *,'vvv mixrho sum smrho output 222=',sum(smrho(:,:,:,1)),sum(abs(smrho(:,:,:,1))),nnnew
ccccccccccccccccccccccccccccccccc
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [okv ogv orhnew orhold ocn owk oa oaa oqkl]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine mixrho 


      subroutine pvmix3(ssite,sspec,nbas,nsp,beta,locmix,wt,kmxr,nlm0,
     .  k9,qkl,orhold,orhnew,difx)
       
       use m_struc_def  !Cgetarg

C- Linearly mix local densities, possibly subtracting G_kL expansion
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: a nr rmt rg lmxl
Ci     Stored:    *
Ci     Passed to: *
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   beta  :linear mixing parameter
Cl Local variables
Cl   locmix :governs linear transformation of local densities for mixing
Cl          :see description in Local variables, routine mixrho.
Ci   kmxr   :k-cutoff for G_kL expansion when locmix=2; also dim. qkl
Ci   nlm0   :second dimension of qkl
Ci   qkl    :(locmix>=2) G_kL expansion of local densities
Ci          : qkl(:,:,isp,1,ib) = rho1(rhold), G_kL expansion
Ci          : qkl(:,:,isp,2,ib) = rho2(rhold), G_kL expansion
Ci          : qkl(:,:,isp,3,ib) = rho1(rhnew), G_kL expansion
Ci          : qkl(:,:,isp,4,ib) = rho2(rhnew), G_kL expansion
Cio Inputs/Outputs
Cio  orhold :On input, local densities generating hamiltonian
Cio         :transformed into unscaled rho1+rho2 and rho1-rho2 (pvmix9)
Cio         :If a G_kL expansion has been generated (if locmix>=2)
Cio         :this expansion is first subtracted from w(orhold).
Cio         :On output, w(orhold) is overwritten by the linear
Cio         :combination (1-beta)*w(orhold) + beta*w(orhnew)
Cio  orhnew :On input, local densities gen. by ham. (maybe screened)
Cio         :transformed into unscaled rho1+rho2 and rho1-rho2 (pvmix9)
Cio         :If a G_kL expansion has been generated (if locmix>=2)
Cio         :this expansion is first subtracted from w(orhnew).
Co Outputs
Co   difx   :maximum rms difference in rhonew-rhoold
Cr Remarks
Cr
Cu Updates
Cu   21 Nov 01 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,nsp,kmxr,nlm0,locmix,k9
Cdelw1 do not change  [['call'], 'orhold', 'pvmix3']
Cdelw1 do not change  [['call'], 'orhnew', 'pvmix3']
      integer orhold(3,nbas),orhnew(3,nbas)
Cgetarg       double precision ssite(1),sspec(1),difx,beta,wt(2),
Cgetarg      .  qkl(0:kmxr,nlm0,nsp,4,nbas)
       real(8):: difx , beta , wt(2) , qkl(0:kmxr,nlm0,nsp,4,nbas) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)

C ... Local parameters
      integer ib,is,igetss,nr,nlml,m,lmxl,k9l
Cdelw1 do not change  [['array', 'wref', 'call'], 'orsm', 'pvmix3']
Cdelw1       integer ori,orwgt,orsm(4),i
       integer:: orsm(4) , i 
       real(8) ,allocatable :: rv_w_ori(:)
       real(8) ,allocatable :: rv_w_orwgt(:)

      double precision a,rmt,rg,difa,rf
C ... Heap
      integer w(1)
      common /w/ w

      difx  = 0
      do  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg...info...           real(8) :: rg , real(8) :: sspec%rg 1
Cgetarg         call upack('spec a nr rmt rg',sspec,is,a,nr,rmt,rg)
         
         a=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt
         rg=sspec(is)%rg

Cgetarg         lmxl = igetss('spec lmxl',is,sspec)
         lmxl = int(sspec(is)%lmxl) 

        if (lmxl .eq. -1) goto 10
        nlml = (lmxl+1)**2
Cdelw1         call defrr(ori, nr)
         allocate(rv_w_ori(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) rv_w_ori(:)=0.0d0

Cdelw1         call defrr(orwgt, nr)
         allocate(rv_w_orwgt(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) rv_w_orwgt(:)=0.0d0

Cdelw1         call radmsh(rmt,a,nr,w(ori))
         call radmsh ( rmt , a , nr , rv_w_ori ) 

Cdelw1         call radwgt(rmt,a,nr,w(orwgt))
         call radwgt ( rmt , a , nr , rv_w_orwgt ) 

        rf = rmt/3

        do  m = 1, 4
          if (locmix .ge. 2) then
C           Generate on a mesh the smooth density to be subtracted
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
            call defrr(orsm(m),nr*nlml*nsp)
            i = 1001
            if (wt(1) .eq. 0) i = 11001
            if (wt(2) .eq. 0) i = 21001
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
Cdelw1             call pkl2ro(i,1,rg,kmxr,nr,nlml,nsp,w(ori),w(orwgt),
Cdelw1      .        kmxr,nlm0,w,qkl(0,1,1,m,ib),w(orsm(m)),w,difa)
             call pkl2ro ( i , 1 , rg , kmxr , nr , nlml , nsp , rv_w_ori 
     .       , rv_w_orwgt , kmxr , nlm0 , w , qkl ( 0 , 1 , 1 , m , ib ) , 
     .       w ( orsm ( m ) ) , w , difa ) 


C           Undo scaling of rho1+rho2 for linear mix
            if (locmix .eq. 3 .and. (m .eq. 2 .or. m .eq. 4)) then

C             qlk(1,2) -> scaled rho1+rho2, rho1-rho2
C             qkl(3,4) -> scaled rhn1+rhn2, rhn1-rhn2
C             if (m .eq. 2) print *, m, 'old'
C             if (m .eq. 4) print *, m, 'new'
C             call prrmsh('scaled rhosm1+2 from qkl',w(ori),w(orsm(m-1)),
C     .         nr,nr,nlml*nsp)
C             call prrmsh('scaled rhosm1-2 from qkl',w(ori),w(orsm(m)),
C     .         nr,nr,nlml*nsp)
              k9l = k9 - mod(k9,10)
Cdelw1 1 is digit, dropped
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
Cdelw1               call pvmix9(k9l,-1,nr,nlml*nsp,0,rf,w(ori),
Cdelw1      .            w(orsm(m-1)),w(orsm(m)))
               call pvmix9 ( k9l , - 1 , nr , nlml * nsp , 0 , rf , rv_w_ori 
     .         , w ( orsm ( m - 1 ) ) , w ( orsm ( m ) ) ) 

C             call prrmsh('unscaled rhosm1+2',w(ori),w(orsm(m-1)),
C     .         nr,nr,nlml*nsp)
C             call prrmsh('unscaled rhosm1-2',w(ori),w(orsm(m)),
C     .         nr,nr,nlml*nsp)
            endif

          else
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
            orsm(m) = 1
          endif
        enddo

C   ... Constraints on rnew, smoothed rho
C        i = 10
C        if (wt(1) .eq. 0) i = 110
C        if (wt(2) .eq. 0) i = 210
C        call pvmix8(i,nr,nlml,nsp,w(orhnew(1,ib)),w(orhnew(2,ib)))
C        call pvmix8(i,nr,nlml,nsp,w(orsm(1)),w(orsm(2)))
C        call pvmix8(i,nr,nlml,nsp,w(orsm(3)),w(orsm(4)))

C   ... Linear mix, this site
Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1 3 is digit, dropped
Cdelw1 4 is digit, dropped
Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1 3 is digit, dropped
Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1 3 is digit, dropped
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
Cdelw1         call pvmix4(nr,nlml,nsp,w(ori),locmix,wt,beta,ib,
Cdelw1      .    w(orsm(1)),w(orsm(2)),w(orsm(3)),w(orsm(4)),
Cdelw1      .    w(orhold(1,ib)),w(orhold(2,ib)),w(orhold(3,ib)),
Cdelw1      .    w(orhnew(1,ib)),w(orhnew(2,ib)),w(orhnew(3,ib)),w(orwgt),difa)
         call pvmix4 ( nr , nlml , nsp , rv_w_ori , locmix , wt , beta 
     .   , ib , w ( orsm ( 1 ) ) , w ( orsm ( 2 ) ) , w ( orsm ( 3 ) ) 
     .   , w ( orsm ( 4 ) ) , w ( orhold ( 1 , ib ) ) , w ( orhold ( 2 
     .   , ib ) ) , w ( orhold ( 3 , ib ) ) , w ( orhnew ( 1 , ib ) ) 
     .   , w ( orhnew ( 2 , ib ) ) , w ( orhnew ( 3 , ib ) ) , rv_w_orwgt 
     .   , difa ) 

        difx  = dmax1(difx,difa)

Cdelw1 rlse name= ori old_list= ori orwgt orsm 
Cdelw1 rlse name= ori new_list= (None)
Cdelw1 not deallocate orsm because of [array wref call]
Cdelw1         call rlse(ori)
         if (allocated(rv_w_orwgt)) deallocate(rv_w_orwgt)
         if (allocated(rv_w_ori)) deallocate(rv_w_ori)
         call rlse(orsm)

   10   continue
      enddo
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [orsm orhold orhnew]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pvmix3 


      subroutine pvmix4(nr,nlml,nsp,ri,locmix,wt,beta,ib,
     .  rhos1,rhos2,rhns1,rhns2,
     .  rho1,rho2,rho3,rhn1,rhn2,rhn3,rwgt,dif)

C- Linear mixing of local densities for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for rho1,rho2
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ri    :radial mesh points
Ci   locmix:>2 => work with difference rho-rhos
Ci   beta  :Linear mixing beta
Ci   ib    :site index (printout only)
Ci   rhos1 :input local true density, G_kL expansion
Ci   rhos2 :input local smoothed density, G_kL expansion
Ci   rhns1 :output local true density, G_kL expansion
Ci   rhns2 :output local smoothed density, G_kL expansion
Ci   rho1  :input local true density (generates output rho)
Ci   rho2  :input local smoothed density (generates output rho)
Ci   rho3  :input core charge (generates output rho)
Ci   rhn1  :output local true density (generated by input rho)
Ci   rhn2  :output local smoothed density (generated by input rho)
Ci   rhn3  :output core charge (generated by input rho)
Ci   rwgt  :radial mesh weights
Co Outputs
Co   dif   :integral (rho_out - rho_in)**2
Co   rho1,rho2,rho3 are overwritten by the linearly mixed densities
Cr Remarks
Cu Updates
Cu   09 Dec 08 Suppress q,moment mixing when wt(1) or wt(2) is zero
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,nlml,nsp,ib,locmix
      double precision dif,beta,wt(3)
      double precision ri(nr),rwgt(nr),
     .  rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),rho3(nr,nsp),
     .  rhn1(nr,nlml,nsp),rhn2(nr,nlml,nsp),rhn3(nr,nsp),
     .  rhos1(nr,nlml,nsp),rhos2(nr,nlml,nsp),
     .  rhns1(nr,nlml,nsp),rhns2(nr,nlml,nsp)
C ... Local parameters
      integer i,ir,iprint,stdo,lgunit,np,npc
      double precision pi,srfpi,vsph,ddot,qnew,qold,qmix,rmt,
     .  amnew,amold,amix

C      print *, 'start of pvmix4'
C      call prrmsh('rho1',ri,rho1,nr,nr,nlml*nsp) ! rho1+rho2
C      call prrmsh('rho2',ri,rho2,nr,nr,nlml*nsp) ! rho1-rho2
C      call prrmsh('rhn1',ri,rhn1,nr,nr,nlml*nsp) ! rhscr1+rhscr2
C      call prrmsh('rhn2',ri,rhn2,nr,nr,nlml*nsp) ! rhscr1-rhscr2
C
C      call prrmsh('rhos1',ri,rhos1,nr,nr,nlml*nsp) ! rhos1+rhos2
C      call prrmsh('rhos2',ri,rhos2,nr,nr,nlml*nsp) ! rhos1-rhos2
C      call prrmsh('rhns1',ri,rhns1,nr,nr,nlml*nsp) ! rhscrs1+rhscrs2
C      call prrmsh('rhns2',ri,rhns2,nr,nr,nlml*nsp) ! rhscrs1-rhscrs2


C ... Setup
      stdo = lgunit(1)
      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)
      rmt = ri(nr)

C ... rho+,rho- -> rho,amag  and some integrated quantities
      call splrho(0,nsp,nr,nlml,rho1,rho2,rho3)
      call splrho(0,nsp,nr,nlml,rhn1,rhn2,rhn3)
C      qnew = srfpi*(ddot(nr,rhn1,1,rwgt,1) - ddot(nr,rhn2,1,rwgt,1))
C      qold = srfpi*(ddot(nr,rho1,1,rwgt,1) - ddot(nr,rho2,1,rwgt,1))
C      amnew = srfpi*(ddot(nr,rhn1(1,1,nsp),1,rwgt,1) -
C     .               ddot(nr,rhn2(1,1,nsp),1,rwgt,1))
C      amold = srfpi*(ddot(nr,rho1(1,1,nsp),1,rwgt,1) -
C     .               ddot(nr,rho2(1,1,nsp),1,rwgt,1))
      qnew = srfpi*ddot(nr,rhn2,1,rwgt,1)
      qold = srfpi*ddot(nr,rho2,1,rwgt,1)
      amnew = srfpi*ddot(nr,rhn2(1,1,nsp),1,rwgt,1)
      amold = srfpi*ddot(nr,rho2(1,1,nsp),1,rwgt,1)

      dif = 0
      do  i  = 1, nsp
        if (wt(1) .eq. 0 .and. i .eq. 1) cycle
        if (wt(2) .eq. 0 .and. i .eq. 2) cycle
        do  ir = 1, nr
          dif  = dif + rwgt(ir)*(rhn2(ir,1,i)-rho2(ir,1,i))**2
        enddo
      enddo
      vsph = (4d0*pi/3d0)*rmt**3
      dif = dsqrt(dif/vsph)

C --- Overwrite rho with (1-beta)*rho + beta*rhn ---
C     wt(1)=0 => only mix spin part (set i=2)
C     wt(2)=0 => only mix charge part (set i=1)
      np = nr*nlml*nsp
      i = 1
      if (wt(1) .eq. 0) then
        np = nr*nlml
        i = 2
      elseif (wt(2) .eq. 0) then
        np = nr*nlml
        i = 1
      endif
      call dscal(np,1-beta,rho1(1,1,i),1)
      call daxpy(np,beta,rhn1(1,1,i),1,rho1(1,1,i),1)
      call dscal(np,1-beta,rho2(1,1,i),1)
      call daxpy(np,beta,rhn2(1,1,i),1,rho2(1,1,i),1)
      npc = nr*nsp
      if (wt(1) .eq. 0 .or. wt(2) .eq. 0) npc = nr
      call dscal(npc,1-beta,rho3(1,i),1)
      call daxpy(npc,beta,rhn3(1,i),1,rho3(1,i),1)

C ... Charges and mag. moments of mixed density
      qmix = srfpi*(ddot(nr,rho2,1,rwgt,1))
      if (nsp .eq. 2) amix = srfpi*ddot(nr,rho2(1,1,nsp),1,rwgt,1)

C --- Subtract smoothed (1-beta)*rhos + beta*rhns ---
      if (locmix .ge. 2) then
        call splrho(20,nsp,nr,nlml,rhos1,rhos2,rho3)
        call splrho(20,nsp,nr,nlml,rhns1,rhns2,rhn3)
        np = nr*nlml*nsp
        call daxpy(np,-(1-beta),rhos1,1,rho1,1)
        call daxpy(np,-beta,rhns1,1,rho1,1)
        call daxpy(np,-(1-beta),rhos2,1,rho2,1)
        call daxpy(np,-beta,rhns2,1,rho2,1)
        call splrho(21,nsp,nr,nlml,rhos1,rhos2,rho3)
        call splrho(21,nsp,nr,nlml,rhns1,rhns2,rhn3)
      endif

C ... Restore rho,amag -> rho+,rho-
      call splrho(1,nsp,nr,nlml,rho1,rho2,rho3)
      call splrho(1,nsp,nr,nlml,rhn1,rhn2,rhn3)

      if (iprint() .gt. 30) then
        write(stdo,100) ib,qold,qnew,qnew,dif,qmix
  100   format(' site ',i4,f12.6,5f14.6)
C  100   format(' site ',i4,f12.6,f14.6,14x,2f14.6)
        if (nsp .eq. 2) write(stdo,101) amold,amnew,amix
  101   format(' mmom   ',2f14.6,28x,f14.6)
      endif

C       print *, 'end of pvmix4'
C       call prrmsh('rho1 after lin mix',ri,rho1,nr,nr,nlml*nsp)
C       call prrmsh('rho2 after lin mix',ri,rho2,nr,nr,nlml*nsp)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pvmix4 


      subroutine pvmix5(nmix,mxsav,fnam,ifi,lbin,rmsdel,locmix,k9,nbas,
     .  kmxr,nlm0,ssite,sspec,nsp,orhold,orhnew,co,cn,ng,ng0,cnst,nda,a,
     .  qkl,rms2,nmixr)
       
       use m_struc_def  !Cgetarg

C- Copy rho into holding array, read prior iterations from disk
C ----------------------------------------------------------------------
Ci Inputs
Ci   nmix  :number of prior iterations sought (for printout)
Ci   mxsav :max no. prior iteration to read, also dimensions a
Ci   fnam  :file name (for printout)
Ci   ifi   :file logical unit; sign delineates whether copying into
Ci         :or out of holding array
Ci   lbin  :T if binary read
Ci  rmsdel :Same as rms2 (see Outputs), from prior iteration.
Ci         :If no prior iteration, rsmdel=0.  For printout only.
Ci  locmix :switch governing linear transformation of local densities
Ci          for mixing; see Local variables in subroutine mixrho.
Ci  nbas   :size of basis
Ci  ssite  :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: rhogkl
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: a nr rmt lmxl
Ci     Stored:    *
Ci     Passed to: rhogkl
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   orhold:input local density this iteration,
Ci         :transformed into unscaled rho1+rho2 and rho1-rho2 (pvmix9)
Ci   orhnew:output local density this iteration,
Ci         :transformed into unscaled rho1+rho2 and rho1-rho2 (pvmix9)
Ci   co    :input smooth density this iteration, in FT form
Ci   cn    :output smooth density this iteration, in FT form
Ci   ng    :leading dimension of co,cn
Ci   ng0   :number of Fourier coefficients
Ci   cnst  :site-dependent vector marking whether to constrain
Ci         :selected mixing (not used here)
Co Inputs/Outputs
Cio  nmixr :(output for ifi>0) number of prior iter read from file
Cio        :(input  for ifi<0) number of prior iter written to file
Cio        :on output number of prior iterations read from file
Cio  a,nda :work array a and leading dimension (nda must be >= na)
Cio         a(*,i,1) holds f(xi) for iteration i-1
Cio         a(*,i,2) holds   xi  for iteration i-1
Co Outputs
Co  qkl    :(locmix>=2) G_kL expansion of local densities
Co         : qkl(:,:,isp,1,ib) = rho1(rhold), G_kL expansion
Co         : qkl(:,:,isp,2,ib) = rho2(rhold), G_kL expansion
Co         : qkl(:,:,isp,3,ib) = rho1(rhnew), G_kL expansion
Co         : qkl(:,:,isp,4,ib) = rho2(rhnew), G_kL expansion
Co  rms2   :rms difference in f-x, where f=output density and x the
Co          input density.  f and x include smoothed and local
Co          contributions, with the latter depending on choice of
Co          linear transformation; see locmix above.  Thus, there
Co          is some arbitrariness in rms2.
Cr Remarks
Cr   Content of a(*) = a(*,isp,:,:):
Cr   Case locmix = 0:
Cr     a(1:2*ng0) = Fourier coffs of sm density (hermitian => ng0<ng)
Cr     a(1+2*ng0:2*ng0+2*nr) = l=0 part of rho(1), rho(2)
Cr     ... nda = 2*ng0 + 2*nr
Cr   Case locmix = 1:
Cr     a(1:2*ng0) = Fourier coffs of sm density (hermitian => ng0<ng)
Cr     a(1+2*ng0:2*ng0+2*nr) = Weighted l=0 part of rho(1), rho(2)
Cr     ... nda = 2*ng0 + 2*nr
Cr   Case locmix = 2
Cr     a(1:2*ng0) = Fourier coffs of sm density (hermitian => ng0<ng)
Cr     a(1+2*ng0:2*ng0+(1+kmxr)*nlml) = P_kL expansion of rho1-rho2
Cr     ... nda = 2*ng0 + (1+kmxr)*nlml
Cr   Case locmix = 3
Cr     a(1:2*ng0) = Fourier coffs of sm density (hermitian => ng0<ng)
Cr     a(1+2*ng0:2*ng0+2*nr) = Weighted l=0 part of rho(1), rho(2)
Cr     a(1+k:2*ng0:2*ng0+2*nr+(1+kmxr)*nlml) = P_kL exp. of rho1-rho2
Cr                                             w/ l=0 part zeroed out
Cr     ... nda = 2*ng0 + 2*nr + (1+kmxr)*nlml
C -----------------------------------------------------------------
C     implicit none
C ... Passed parameters
#if MPI | MPIK
      include "mpif.h"
      integer procid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
      character*120 strn
      logical mlog,cmdopt
      integer iunit
      logical readerror
#endif
      logical lbin
      integer ng,ng0,nda,nmix,mxsav,ifi,nbas,nr,nsp,cnst(0:1),locmix,
     .  k9,kmxr,nlm0
Cdelw1 do not change  [['call'], 'orhold', 'pvmix5']
Cdelw1 do not change  [['call'], 'orhnew', 'pvmix5']
      integer orhold(3,1),orhnew(3,1)
      double precision a(nda,nsp,mxsav+2,2),rms2,rmsdel
Cgetarg       double precision ssite(1),sspec(1),co(2,ng,nsp),cn(2,ng,nsp),
Cgetarg      .  qkl(0:kmxr,nlm0,nsp,4,1)
       real(8):: co(2,ng,nsp) , cn(2,ng,nsp) , qkl(0:kmxr,nlm0,nsp,4,1) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)

      character fnam*8
C ... Local parameters
      integer stdo,lgunit,ib,na,i,j,k,m,np,iprint,nmixr,is,igetss,
     .  awrite,off,nlml,lmxl,k9l
Cdelw1       integer orofi
       real(8) ,allocatable :: rv_w_orofi(:)

      logical lcnst,lfdmp,lddump
      double precision ddot,rmt,aat,rf

      character outs*80
C ... Heap
      integer w(1)
      common /w/ w

#if MPI | MPIK
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
#endif

      stdo = lgunit(1)
      lcnst = cnst(0) .gt. 0

C --- Copy rho to a; calc. RMS error; read prior iter from ifi ---
      if (ifi .gt. 0) then

C   ... FT Mesh density
        na = 1
        do  i = 1, nsp
          call dpcopy(co(1,1,i),a(na+00,i,1,2),1,2*ng0,1d0)
          call dpcopy(cn(1,1,i),a(na+00,i,1,1),1,2*ng0,1d0)
        enddo
        na = 1 + 2*ng0

C   ... Site densities, depending on locmix
        do  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           real(8) :: aat , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg         call upack('spec a nr rmt lmxl',sspec,is,aat,nr,rmt,lmxl)
         
         aat=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt
         lmxl=sspec(is)%lmxl

        if (lmxl .eq. -1) cycle
        rf = rmt/3
        nlml = (lmxl+1)**2
Cdelw1         call defrr(orofi, nr)
         allocate(rv_w_orofi(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) rv_w_orofi(:)=0.0d0

Cdelw1         call radmsh(rmt,aat,nr,w(orofi))
         call radmsh ( rmt , aat , nr , rv_w_orofi ) 


C   ... If site is constrained, freeze density at rhold (not checked)
        if (lcnst) then
          if (cnst(ib) .ne. 0) then
            do  i = 1, nsp
            off = 1+nr*nlml*(i-1)
Cdelw1 1 is digit, dropped
Cdelw1 1 is digit, dropped
            call dpscop(w(orhold(1,ib)),w(orhnew(1,ib)),nr,off,off,1d0)
Cdelw1 2 is digit, dropped
Cdelw1 2 is digit, dropped
            call dpscop(w(orhold(2,ib)),w(orhnew(2,ib)),nr,off,off,1d0)
            enddo
          endif
        endif

C   --- Copy local densities to mix matrix a, depending on locmix ---

C   ... Overwrite rho1,rho2 with linear transf. according to k9
CXX     mch tmp/rho1 tmp/rho2 -+ out.fe -- -coll 2,18  -px
C        call prrmsh('rhold1',w(orofi),w(orhold(1,1)),nr,nr,nlml*nsp)
C        call prrmsh('rhold2',w(orofi),w(orhold(2,1)),nr,nr,nlml*nsp)

C   ... locmix=0,1,3 : copy spherical part of local rho to a
C       NB: Increment na after all actions takend w/ local densities
        if (locmix .le. 1 .or. locmix .eq. 3) then

C         Scaling rho1+rho2, old and new densities
          k9l = k9 - mod(k9,10)
Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1           call pvmix9(k9l,0,nr,nlml*nsp,0,rf,w(orofi),w(orhold(1,ib)),
Cdelw1      .      w(orhold(2,ib)))
           call pvmix9 ( k9l , 0 , nr , nlml * nsp , 0 , rf , rv_w_orofi 
     .     , w ( orhold ( 1 , ib ) ) , w ( orhold ( 2 , ib ) ) ) 

Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1           call pvmix9(k9l,0,nr,nlml*nsp,0,rf,w(orofi),w(orhnew(1,ib)),
Cdelw1      .      w(orhnew(2,ib)))
           call pvmix9 ( k9l , 0 , nr , nlml * nsp , 0 , rf , rv_w_orofi 
     .     , w ( orhnew ( 1 , ib ) ) , w ( orhnew ( 2 , ib ) ) ) 

C         call prrmsh('1 post',w(orofi),w(orhold(1,1)),nr,nr,nlml*nsp)
C         call prrmsh('2 post',w(orofi),w(orhold(2,1)),nr,nr,nlml*nsp)
C         Copy (scaled) spherical part of rho1+rho2, rho1-rho2 to a
          do  i = 1, nsp
            off = 1+nr*nlml*(i-1)
Cdelw1 1 is digit, dropped
            call dpscop(w(orhold(1,ib)),a(na+00,i,1,2),nr,off,1,1d0)
Cdelw1 2 is digit, dropped
            call dpscop(w(orhold(2,ib)),a(na+nr,i,1,2),nr,off,1,1d0)
Cdelw1 1 is digit, dropped
            call dpscop(w(orhnew(1,ib)),a(na+00,i,1,1),nr,off,1,1d0)
Cdelw1 2 is digit, dropped
            call dpscop(w(orhnew(2,ib)),a(na+nr,i,1,1),nr,off,1,1d0)
          enddo

C         Undo scaling transformation of rho1+rho2
C          call pvmix9(k9l,-1,nr,nlml*nsp,0,rf,w(orofi),w(orhold(1,ib)),
C     .      w(orhold(2,ib)))
C          call pvmix9(k9l,-1,nr,nlml*nsp,0,rf,w(orofi),w(orhnew(1,ib)),
C     .      w(orhnew(2,ib)))
C         call prrmsh('1 rest',w(orofi),w(orhold(1,1)),nr,nr,nlml*nsp)
C         call prrmsh('2 rest',w(orofi),w(orhold(2,1)),nr,nr,nlml*nsp)

        endif
C       print *, 'a(nr)'
C       print *, a(na+nr-1,:,1,:)

C   ... locmix>=2 : Make G_kL expansion of rho; copy coffs to a
        if (locmix .ge. 2) then
          do  i = 1, nsp
            call pshpr(iprint()-30)
            call rhogkl(ib,ib,i,2,ssite,sspec,orhold,kmxr,
     .        qkl(0,1,i,1,ib))
            call rhogkl(ib,ib,i,3,ssite,sspec,orhold,kmxr,
     .        qkl(0,1,i,2,ib))
            call rhogkl(ib,ib,i,2,ssite,sspec,orhnew,kmxr,
     .        qkl(0,1,i,3,ib))
            call rhogkl(ib,ib,i,3,ssite,sspec,orhnew,kmxr,
     .        qkl(0,1,i,4,ib))
            call poppr
C           Convert qkl for 2nd spin channel into 2nd channel only,
C           since rhogkl makes expansion combining channels 1+2
            if (i .eq. 2) then
              j = (1+kmxr)*nlm0
C             qkl(nsp=2) -> qkl(nsp=1) - qkl(nsp=2)
              do  m = 1, 4
               call dpsadd(qkl(0,1,1,m,ib),qkl(0,1,1,m,ib),j,1+j,1,-1d0)
              enddo
            endif

C           locmix=3 : zero out spherical part of G_kL
            if (locmix .eq. 3) then
C             print *, '!! skip zero qkl'
              do  m = 1, 4
                call dpzero(qkl(0,1,i,m,ib),kmxr+1)
              enddo
            endif

C           locmix=3 : increment index by 2*nr
            if (locmix .eq. 3) na = na + 2*nr

C           Copy qkl into a
            np = (1+kmxr)*nlml
            call dpscop(qkl(0,1,i,1,ib),a(na+00,i,1,2),np,1,1,1d0)
            call dpscop(qkl(0,1,i,2,ib),a(na+np,i,1,2),np,1,1,1d0)
            call dpscop(qkl(0,1,i,3,ib),a(na+00,i,1,1),np,1,1,1d0)
            call dpscop(qkl(0,1,i,4,ib),a(na+np,i,1,1),np,1,1,1d0)

C           locmix=3, : undo increment unless final spin
            if (locmix .eq. 3 .and. i .lt. nsp) na = na - 2*nr

          enddo                 ! end of spin loop
C         call prmx('qkl',qkl,1+kmxr,1+kmxr,nlm0*nsp)
        endif

C       Undo scaling transformation of rho1+rho2
        if (locmix .le. 1 .or. locmix .eq. 3) then
Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1           call pvmix9(k9l,-1,nr,nlml*nsp,0,rf,w(orofi),w(orhold(1,ib)),
Cdelw1      .      w(orhold(2,ib)))
           call pvmix9 ( k9l , - 1 , nr , nlml * nsp , 0 , rf , rv_w_orofi 
     .     , w ( orhold ( 1 , ib ) ) , w ( orhold ( 2 , ib ) ) ) 

Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1           call pvmix9(k9l,-1,nr,nlml*nsp,0,rf,w(orofi),w(orhnew(1,ib)),
Cdelw1      .      w(orhnew(2,ib)))
           call pvmix9 ( k9l , - 1 , nr , nlml * nsp , 0 , rf , rv_w_orofi 
     .     , w ( orhnew ( 1 , ib ) ) , w ( orhnew ( 2 , ib ) ) ) 

C         call prrmsh('1 rest',w(orofi),w(orhold(1,1)),nr,nr,nlml*nsp)
C         call prrmsh('2 rest',w(orofi),w(orhold(2,1)),nr,nr,nlml*nsp)
        endif

C   ... Increment index na to a
        if (locmix .le. 1) then
          na = na + 2*nr
        elseif (locmix .ge. 2) then
          na = na + 2*np
        endif

Cdelw1 rlse name= orofi old_list= orofi 
Cdelw1 rlse name= orofi new_list= (None)
Cdelw1         call rlse(orofi)
         if (allocated(rv_w_orofi)) deallocate(rv_w_orofi)

        enddo                     !Loop over sites
        na = na-1
        if (nda .ne. na) call rx('mixrho: bug in pvmix5')
        na = na*nsp
C       call prmx('a',a,nda,na,1)
        rms2 = dsqrt(dabs(ddot(na,a,1,a,1)     -
     .             2*ddot(na,a,1,a(1,1,1,2),1) +
     .               ddot(na,a(1,1,1,2),1,a(1,1,1,2),1))/(na-0))*nsp
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [orhold orhnew]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
        if (mxsav .eq. 0) return
C   ... File I/O of previous iterations; nmixr is number actually read
        nmixr = 0
#if ! (MPI | MPIK)
        if (lbin) read (ifi,err=31,end=31) nmixr, na
        if (.not. lbin) read (ifi,*,err=31,end=31) nmixr, na
        if (nda*nsp .ne. na) then
          call awrit2(' mixrho:  expecting %i elements but found %i ...'
     .      //' discarding file',' ',80,stdo,nda*nsp,na)
          nmixr = 0
          goto 31
        endif
        k = min(mxsav,nmixr)
        do  30  j = 1, k
          if (lbin) then
            if (.not. lddump(a(1,1,j+1,1),na,ifi)) goto 31
            if (.not. lddump(a(1,1,j+1,2),na,ifi)) goto 31
          else
            if (.not. lfdmp(a(1,1,j+1,1),na,ifi)) goto 31
            if (.not. lfdmp(a(1,1,j+1,2),na,ifi)) goto 31
          endif
          nmixr = j
   30   continue
   31   continue
#else
        readerror = .true.
        if (procid .eq. master) then
          readerror = .false.
          if (lbin) read (ifi,err=311,end=311) nmixr, na
          if (.not. lbin) read (ifi,*,err=311,end=311) nmixr, na
          goto 312
  311     continue
          readerror = .true.
  312     continue
        endif
        call MPI_BCAST(readerror,1,MPI_LOGICAL,
     .                 master,MPI_COMM_WORLD,ierr)
        if (mlog) then
          call gettime(datim)
          call awrit3(' mixrho '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' bcast readerror = %l',' ',
     .      256,lgunit(3),procid,numprocs,readerror)
        endif
        if (readerror) goto 31
        call MPI_BCAST(nmixr,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
        if (mlog) then
          call gettime(datim)
          call awrit3(' mixrho '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' bcast nmixr = %i',' ',
     .      256,lgunit(3),procid,numprocs,nmixr)
        endif
        call MPI_BCAST(na,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
        if (mlog) then
          call gettime(datim)
          call awrit3(' mixrho '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' bcast na = %i',' ',
     .      256,lgunit(3),procid,numprocs,na)
        endif
        if (nda*nsp .ne. na) then
          if (procid .eq. master) then
            call awrit2
     .        (' mixrho:  expecting %i elements but found %i ...'
     .        //' discarding file',' ',80,stdo,nda*nsp,na)
          endif
          nmixr = 0
          goto 31
        endif
        k = min(mxsav,nmixr)
        if (procid .eq. master) then
          readerror = .false.
          do  30  j = 1, k
            if (lbin) then
              if (.not. lddump(a(1,1,j+1,1),na,ifi)) readerror = .true.
              if (.not. lddump(a(1,1,j+1,2),na,ifi)) readerror = .true.
            else
              if (.not. lfdmp(a(1,1,j+1,1),na,ifi)) readerror = .true.
              if (.not. lfdmp(a(1,1,j+1,2),na,ifi)) readerror = .true.
            endif
            if (readerror) goto 313
            nmixr = j
   30     continue
        endif
  313   continue
        if (nmixr .gt. 0) then
          call MPI_BCAST(a,nda*nsp*(mxsav+2)*2,MPI_DOUBLE_PRECISION,
     .                   master,MPI_COMM_WORLD,ierr)
          if (mlog) then
            call gettime(datim)
            call awrit3(' mixrho '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' bcast (mix) a %i d.p. numbers',' ',
     .        256,lgunit(3),procid,numprocs,nda*nsp*(mxsav+2)*2)
          endif
        endif
   31   continue
#endif
        if (iprint() .ge. 20) then
          j = awrite('%x mixrho:  sought %i iter from file',
     .      outs,len(outs),0,nmix,0,0,0,0,0,0,0)
          outs(j+2:len(outs)) = fnam
          call awrit2('%a; read %i.  RMS DQ=%1,3;3e',outs,80,0,nmixr,
     .      rms2)
          if (rmsdel .ne. 0)
     .      call awrit1('%a  last it=%1,3;3e',outs,80,0,rmsdel)
          call awrit0('%a',outs,-len(outs),-stdo)
        endif
C       call prm('a after read',a,nda*nsp,nda*nsp,(mxsav+2)*2)

C --- Save iterations on disk ---
      else
#if MPI | MPIK
        if (procid .eq. master) then
#endif
        rewind (-ifi)
        if (lbin) write (-ifi) min(nmix+1,mxsav), nda*nsp
        if (.not. lbin) write (-ifi,*) min(nmix+1,mxsav), nda*nsp
        do  j = 1, min(nmix+1,mxsav)
          if (lbin) then
            call dpdump(a(1,1,j,1),nda*nsp,ifi)
            call dpdump(a(1,1,j,2),nda*nsp,ifi)
          else
            call dfdump(a(1,1,j,1),nda*nsp,ifi)
            call dfdump(a(1,1,j,2),nda*nsp,ifi)
          endif
        enddo
C       call prm('a after write',a,nda*nsp,nda*nsp,(mxsav+2)*2)
#if MPI | MPIK
        endif
#endif
      endif

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [orhold orhnew]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pvmix5 


      subroutine pvmix6(broy,nmix,mmix,mxsav,beta,wc,nda,a)

C- Mixing of the total density
C ------------------------------------------------------------------
Ci  broy   : 0 for Anderson mixing
Ci         : 1 for Broyden mixing
Ci  nmix   : nmix: number of iter to try and mix
Ci  mmix   : number of iterates available to mix
Ci  mxsav  : second dimension of a
Cio Inputs/Outputs
Cio  a     :On input:
Cio         (*,i,1)  output vector for prev. iteration i
Cio         (*,i,2)  input  vector for prev. iteration i
Cio        :On output:
Cio         (*,0,1) and (*,0,2) estimated new vector
Cio         (*,i,1), i>1 destroyed
Cr  Remarks
Cr    Broyden mixing adapted from Duane Johnson.  Notation:
Cr    x^(m): input vector for iteration m
Cr    F^(m): difference between output and input vector in iteration m
C ------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nda,nmix,mmix,mxsav,broy
      double precision beta,a(nda,0:mxsav+1,2)
C ... Local parameters
      integer j,iprint,amix
      double precision tjmax,tj(10),norm(100),kpvt(10),wc,wctrue,rms2
C     integer idamax
C     double precision dval,rmskm(100),rmsmx
      parameter (tjmax = 5d0)

C ... This is the result of linear mixing, q=0 smooth rho
C     sqmix = a(1,0,1)*beta + a(1,0,2)*(1-beta)
      real(8):: aaa(nda)

ccccccccccccccccccccccccccccccccccccccccc
c takao simple mixing test
c      print *,' vvv: takao simple mixing test'
c      aaa = a(:,0,2)
c      aaa = .2d0*a(:,0,1)+.8d0*a(:,0,2)
c      aaa = .40d0*a(:,0,1)+.6d0*a(:,0,2)

c      aaa = .9d0*a(:,0,1)+.1d0*a(:,0,2)
c      aaa = 1.d0*a(:,0,1)+.0d0*a(:,0,2)
c       aaa = .5d0*a(:,0,1)+.5d0*a(:,0,2)
c      aaa = .8d0*a(:,0,1)+.2d0*a(:,0,2)
c      call dcopy(nda,aaa,1,a(1,0,1),1)
c      call dcopy(nda,aaa,1,a(1,0,2),1)
c      return
cccccccccccccccccccccccccccccccccccccccccc


C --- Anderson mixing ---
      if (nmix .eq. 0 .or. nmix .eq. 1 .or. broy .eq. 0) then
C ... amix needs f-x for prior iterations
      do  30  j = 1, nmix
   30 call daxpy(nda,-1d0,a(1,j,2),1,a(1,j,1),1)
      call query('beta',4,beta)
      nmix = amix(nda,nmix,mxsav,0,beta,iprint(),tjmax,norm,kpvt,a,tj,
     .  rms2)
C      if (iprint() .ge. 20) write(stdo,101) strn,rms2,(tj(j), j=1,nmix)
C  101 format(1x,a4,':     rms diff',f10.6:'  tj=',5f8.3)
C     call prm('a out',a,nda,nda,(mxsav+2)*2)
C ... Copy x*=a(*,0,2) to a(*,0,1)
      call dcopy(nda,a(1,0,2),1,a(1,0,1),1)

C --- Broyden mixing, Duane Johnson's approach ---
      elseif (broy .eq. 1) then
        call pqmixb(nda,nmix,mmix,mxsav,beta,wc,rms2,a,wctrue)
C   ... Copy x*=a(*,0,2) to a(*,0,1)
        call dcopy(nda,a(1,0,2),1,a(1,0,1),1)
      else
        call rx('pvmix6: bad value for broy')
      endif

C ... Add the change relative to linear mixing
C     sqmix = (a(1,0,1) - sqmix)*vol

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pvmix6 


      subroutine pvmix7(ssite,sspec,nbas,nsp,nda,a,n1,n2,n3,k1,k2,k3,
     .  locmix,wt,k9,kmxr,nlm0,qkl,ng,ng0,kv,ips0,gv,crho,wk,orhold,
     .  smrho)
       
       use m_struc_def  !Cgetarg

C- Restore mixed density into specific arrays
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: a nr rmt lmxl rg
Ci     Stored:    *
Ci     Passed to: *
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nda   :leading dimension of a
Ci   a     :mixed density from pvmix6, smoothed + local densities
Ci         :For local parts, a contains SCALED rho1+rho2, rho1-rho2
Ci         :If locmix =  0,1, or 3, a contains spherical part on mesh
Ci         :If locmix =  2,3 a full rho, GkL expansion
Ci   n1..3 :FT mesh
Ci   k1..3 :dimensions smrho
Ci   locmix:mixing mode for local densities
Ci   ng    :number of FT G-vectors
Ci   ng0   :condensed number of G vector (excluding hermitian equiv)
Ci   kv    :indices for gather/scatter operations (gvlist.f)
Ci   ips0  :permutation array mapping ng0 vectors to ng vectors (lgstar)
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   crho  :FT coefficients of smrho(G)
Ci   wk    :complex work array of dimension (k1,k2,k3)
Ci   smrho :smooth density that generated the hamiltonian
Ci   orhold:local  density that generated the hamiltonian,
Ci         :a portion of which which has been linearly mixed (pvmix4)
Ci         :orhold(1) contains rho1+rho2;
Ci         :orhold(2) contains rho1-rho2; see pvmix9
Co Outputs
Co   smrho :overwritten by mixed smooth density
Co   orhold:mixed local density is returned, in rho1+rho2,rho10rho2 form
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nsp,n1,n2,n3,k1,k2,k3,ng,ng0,nda,na,nr,nbas,locmix,k9,
     .  kv(ng,3),ips0(ng),kmxr,nlm0
Cdelw1 do not change  [['call'], 'orhold', 'pvmix7']
      integer orhold(3,1)
Cgetarg       double precision ssite(1),sspec(1),gv(ng,3),a(nda,nsp),
Cgetarg      .  qkl(0:kmxr,nlm0,nsp,4,nbas),rf,wt(2)
       real(8):: gv(ng,3) , a(nda,nsp) , qkl(0:kmxr,nlm0,nsp,4,nbas) 
     . , rf , wt(2) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)

      double complex crho(ng,nsp),smrho(k1,k2,k3,nsp),wk(k1,k2,k3)
C ... Local parameters
Cdelw1       integer orofi,orwgt
       real(8) ,allocatable :: rv_w_orofi(:)
       real(8) ,allocatable :: rv_w_orwgt(:)

Cdelw1 do not change  [['array', 'wref', 'call'], 'orsm', 'pvmix7']
      integer ib,is,igetss,i,m,lmxl,nlml,off,np,k9l,orsm(2)
      double precision aat,rmt,rg,xx
C ... Heap
      integer w(1)
      common /w/ w

C ... Restore mixed smoothed density
      do  i = 1, nsp
C       Add difference of mixed, old to minimize perturbation
ccccccccccccccccc
c        print *,' vvv pvmix7 TEST!!! a=0 test'
c        a=0d0
ccccccccccccccc
        call dscal(2*ng0,-1d0,crho(1,i),1)
        call daxpy(2*ng0,1d0,a(1,i),1,crho(1,i),1)
        call lgstar(2,ng,1,gv,m,ips0,crho(1,i))
        call gvputf(ng,1,kv,k1,k2,k3,crho(1,i),wk)
        call fftz3(wk,n1,n2,n3,k1,k2,k3,1,0,1)
        call daxpy(2*k1*k2*k3,1d0,wk,1,smrho(1,1,1,i),1)
      enddo

C ... Update local densities
      na = 1 + 2*ng0
      k9l = k9 - mod(k9,10)
      do  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           real(8) :: aat , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg         call upack('spec a nr rmt lmxl',sspec,is,aat,nr,rmt,lmxl)
         
         aat=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt
         lmxl=sspec(is)%lmxl

        if (lmxl .eq. -1) cycle
        rf = rmt/3
Cgetarg...info...           real(8) :: rg , real(8) :: sspec%rg 1
Cgetarg         call upack2('spec rg',sspec,is,rg)
         
         rg=sspec(is)%rg

        nlml = (lmxl+1)**2
Cdelw1         call defrr(orofi, nr)
         allocate(rv_w_orofi(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) rv_w_orofi(:)=0.0d0

Cdelw1         call defrr(orwgt, nr)
         allocate(rv_w_orwgt(nr))
Ctakao_ZeroClear_NotRequiered          if (nr<0) rv_w_orwgt(:)=0.0d0

Cdelw1         call radmsh(rmt,aat,nr,w(orofi))
         call radmsh ( rmt , aat , nr , rv_w_orofi ) 

Cdelw1         call radwgt(rmt,aat,nr,w(orwgt))
         call radwgt ( rmt , aat , nr , rv_w_orwgt ) 


C       Scale rho1+rho2 to match scaled spherical parts and sm parts
C       call prrmsh('1 mix7a',w(orofi),w(orhold(1,ib)),nr,nr,nlml*nsp)
Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1         call pvmix9(k9l,0,nr,nlml*nsp,0,rf,w(orofi),
Cdelw1      .    w(orhold(1,ib)),w(orhold(2,ib)))
         call pvmix9 ( k9l , 0 , nr , nlml * nsp , 0 , rf , rv_w_orofi 
     .   , w ( orhold ( 1 , ib ) ) , w ( orhold ( 2 , ib ) ) ) 

C       call prrmsh('1 mix7b',w(orofi),w(orhold(1,ib)),nr,nr,nlml*nsp)

C       Overwrite sph. rhold with scaled mixed rho; do not unscale
        if (locmix .le. 1 .or. locmix .eq. 3) then
          do  i = 1, nsp
            off = 1 + nr*nlml*(i-1)
Cdelw1 1 is digit, dropped
            call dpscop(a(na+00,i),w(orhold(1,ib)),nr,1,off,1d0)
Cdelw1 2 is digit, dropped
            call dpscop(a(na+nr,i),w(orhold(2,ib)),nr,1,off,1d0)
C            call pvmix9(k9l,-1,nr,1,off-1,rf,w(orofi),w(orhold(1,ib)),
C     .        w(orhold(2,ib)))
C          call prrmsh('1 mix7c',w(orofi),w(orhold(1,1)),nr,nr,nlml*nsp)
C          call prrmsh('2 mix7c',w(orofi),w(orhold(2,1)),nr,nr,nlml*nsp)
          enddo
          na = na + 2*nr
        endif

C   ... Mixing mode 2,3: add mixed sm density to existing rhold
        if (locmix .ge. 2) then

C         call prrmsh('1 scale',w(orofi),w(orhold(1,ib)),nr,nr,nlml*nsp)
          np = (1+kmxr)*nlml
          do  m = 1, 2
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
            call defrr(orsm(m),nr*nlml*nsp)
            do  i = 1, nsp
              call dpscop(a(na+np*(m-1),i),qkl(0,1,i,m,ib),np,1,1,1d0)
            enddo
            i = 1001
            if (wt(1) .eq. 0) i = 11001
            if (wt(2) .eq. 0) i = 21001
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
Cdelw1             call pkl2ro(i,1,rg,kmxr,nr,nlml,nsp,w(orofi),w(orwgt),
Cdelw1      .        kmxr,nlm0,w,qkl(0,1,1,m,ib),w(orsm(m)),w,xx)
             call pkl2ro ( i , 1 , rg , kmxr , nr , nlml , nsp , rv_w_orofi 
     .       , rv_w_orwgt , kmxr , nlm0 , w , qkl ( 0 , 1 , 1 , m , ib ) , 
     .       w ( orsm ( m ) ) , w , xx ) 

          enddo
C         Don't undo scaling of smoothed density since rho1+rho2 scaled
C          call pvmix9(k9l,-1,nr,nlml*nsp,0,rf,w(orofi),
C     .      w(orsm(1)),w(orsm(2)))
C         Add to rhold(1,2)
Cdelw1 1 is digit, dropped
Cdelw1 1 is digit, dropped
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
          call daxpy(nr*nlml*nsp,1d0,w(orsm(1)),1,w(orhold(1,ib)),1)
Cdelw1 2 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
          call daxpy(nr*nlml*nsp,1d0,w(orsm(2)),1,w(orhold(2,ib)),1)
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
          call rlse(orsm(2))
Cdelw1 do not change orsm because of ['array', 'wref', 'call']
          call rlse(orsm(1))

          na = na + 2*np
        endif
C       call prrmsh('1 mix7d',w(orofi),w(orhold(1,ib)),nr,nr,nlml*nsp)

Cdelw1 1 is digit, dropped
Cdelw1 2 is digit, dropped
Cdelw1         call pvmix9(k9l,-1,nr,nlml*nsp,0,rf,w(orofi),w(orhold(1,ib)),
Cdelw1      .    w(orhold(2,ib)))
         call pvmix9 ( k9l , - 1 , nr , nlml * nsp , 0 , rf , rv_w_orofi 
     .   , w ( orhold ( 1 , ib ) ) , w ( orhold ( 2 , ib ) ) ) 


C      call prrmsh('rho1,mixed',w(orofi),w(orhold(1,ib)),nr,nr,nlml*nsp)
C      call prrmsh('rho2,mixed',w(orofi),w(orhold(2,ib)),nr,nr,nlml*nsp)

Cdelw1 rlse name= orofi old_list= orofi orwgt orsm 
Cdelw1 rlse name= orofi new_list= (None)
Cdelw1 not deallocate orsm because of [array wref call]
Cdelw1         call rlse(orofi)
         if (allocated(rv_w_orwgt)) deallocate(rv_w_orwgt)
         if (allocated(rv_w_orofi)) deallocate(rv_w_orofi)
         call rlse(orsm)

      enddo                     !Loop over sites
      na = na-1
      if (nda .ne. na) call rx('mixrho: bug in pvmix7')

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [orhold orsm]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pvmix7 


C      subroutine pvmix8(mode,nr,nlml,nsp,rho1,rho2)
CC- Handle possible constraints for change in smooth density
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode  :a compound of digits :
CCi         :10s digit
CCi         :  0 apply to rho1 only; rho2 is not touched
CCi         :  1 apply to both rho1 and rho2
CCi         :100s digit for spin polarized case
CCi         :  1 zero out charge, rho1+ + rho1- (and rho2+ + rho2-)
CCi         :  2 zero out spin, rho1+ - rho1- (and rho2+ + rho2-)
CCi   nr    :number of radial mesh points
CCi   nlml  :L-cutoff for rho1,rho2
CCi   nsp   :2 for spin-polarized case, otherwise 1
CCio Inputs/Outputs
CCio  rho1  :overwritten with constraints controlled by mode
CCio  rho2  :overwritten with constraints controlled by mode
CCr Remarks
CCu Updates
CCu   16 Dec 08 First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer mode,nr,nlml,nsp
C      double precision rho1(nr,nlml,nsp),rho2(nr,nlml,nsp)
CC ... Local parameters
C      integer np,mode0,mode1,mode2,i
C      double precision xx
C
CC ... Setup
C      mode0 = mod(mode,10)
C      mode1 = mod(mode/10,10)
C      mode2 = mod(mode/100,10)
C      np = nr*nlml
C      if (mode2 .eq. 0 .or. nsp .ne. 2) return
C
C      i = 20                   ! No core
C      if (mode1 .eq. 0) i = 30 ! No rho2
C      call splrho(i,nsp,nr,nlml,rho1,rho2,xx)
C      if (mode2 .eq. 1) then   ! Zero density
C        call dpzero(rho1(1,1,1),np)
C        if (mode1 .ne. 0) then ! Including rho2
C          call dpzero(rho2(1,1,1),np)
C        endif
C      endif
C      if (mode2 .eq. 2) then   ! Zero spin
C        call dpzero(rho1(1,1,nsp),np)
C        if (mode1 .ne. 0) then ! Including rho2
C          call dpzero(rho2(1,1,nsp),np)
C        endif
C      endif
C      call splrho(i+1,nsp,nr,nlml,rho1,rho2,xx)
C
C      end
      subroutine pvmix9(mode,linv,nr,nlml,off,rf,rofi,rho1,rho2)

C- Transformation of local densities rho1,rho2 for mixing
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :transformation mode; see Remarks
Ci         : mode  transform
Ci         : 1s digit
Ci         : <=0    do nothing
Ci         :  1     rho1 <= rho1+rho2  rho2 <= rho1-rho2
Ci         : 10s digit
Ci         : 10     scale rho1 by gaussian exp(-(r/rf)**2)
Ci         : 20     scale rho1 by gaussian 1/(1+(r/rf)**2)
Ci   linv  :0  forward transform
Ci         :-1 inverse transform
Ci   nr    :number of radial mesh points
Ci   nlml  :number of L channels
Ci   off   :offset to start of first point (0 for first point)
Ci   rf    :factor used in scaling; see mode
Cio Inputs/Outputs
Ci   rho1  :On input, local true density, or the transformation
Ci         :of it if the inverse is to be made
Ci         :On output, some linear combination of rho1,rho2.
Ci   rho2  :On input, local smoothed density, or the transformation
Ci         :of it if the inverse is to be made
Ci         :On output, some linear combination of rho1,rho2
Cr Remarks
Cr   This routine takes as input rho1,rho2 and returns in
Cr   rho1,rho2 some linear combination of them, depending on mode
Cr   If linv=1, the inverse transform is done.
Cu Updates
Cu   18 Nov 01 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nr,nlml,linv,off
      double precision rho1(nr,nlml),rho2(nr,nlml),rofi(nr),rf
C ... Local parameters
      integer i,ilm,mode0,mode1
      double precision ri

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (mode .le. 0) return
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      call isanrg(linv,-1,0,'pvmix9:','linv', .true.)
      call isanrg(mode0,0,1,'pvmix9:','1s digit mode', .true.)
      call isanrg(mode1,0,2,'pvmix9:','10s digit mode', .true.)

C ... Forward
      if (linv .eq. 0) then
        if (mode0 .eq. 1) then
          call dsumdf(nr*nlml,1d0,rho1,off,1,rho2,off,1)
        endif

        if (mode1 .eq. 1) then
          do  ilm = 1, nlml
          do  i = 1, nr
            ri = rofi(i)/rf
            rho1(i+off,ilm) = rho1(i+off,ilm) * exp(-ri**2)
          enddo
          enddo

        elseif (mode1 .eq. 2) then
          do  ilm = 1, nlml
          do  i = 1, nr
            ri = rofi(i)/rf
            rho1(i+off,ilm) = rho1(i+off,ilm) / (1 + ri**2)
          enddo
          enddo
        endif

C ... Inverse
      else
        if (mode1 .eq. 1) then
          do  ilm = 1, nlml
          do  i = 1, nr
            ri = rofi(i)/rf
            rho1(i+off,ilm) = rho1(i+off,ilm) * exp(ri**2)
          enddo
          enddo

        elseif (mode1 .eq. 2) then
          do  ilm = 1, nlml
          do  i = 1, nr
            ri = rofi(i)/rf
            rho1(i+off,ilm) = rho1(i+off,ilm) * (1 + ri**2)
          enddo
          enddo
        endif

        if (mode0 .eq. 1) then
          call dsumdf(nr*nlml,0.5d0,rho1,off,1,rho2,off,1)
        endif

      endif

C      call prrmsh('pvmix9 rho1',rofi,rho1,nr,nr,nlml)
C      call prrmsh('pvmix9,rho2',rofi,rho2,nr,nr,nlml)

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pvmix9 


