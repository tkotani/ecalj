Cgetarg...info...           structure ['lmfopb', 'spec', 'sspec']
Cgetarg...info...           structure ['lmfopb', 'lat', 'slat']
Cgetarg...info...           structure ['lmfop2', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['lmfop2', 'm_struc_def']
Cgetarg...info...           use_to_add ['lmfopb', 'm_struc_def']
      subroutine lmfopb(sopts,slat,sspec,ehf)
       
       use m_struc_def  !Cgetarg

C- Optimize lmf basis
C ----------------------------------------------------------------------
Ci Inputs
Ci   sopts
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxb rmt
Ci     Stored:    name
Ci     Passed to: spacks uspecb
Co Outputs
Ci   ir    :0 optimization completed
Cl Local variables
Cl   ivar  :index to current variable for minimization
Cl   ncall :number of function calls for which energy aleady calculated
Cl         :for variable being minimized.
Cl   iter  :number of function calls so far in search for min.
Cl         :When iter reaches ncall, lmfopb returns with ir=-1
Cl   vhist :table of energy values for prior 1..ncall points
Cr Remarks
Cr
Cu Updates
Cu   28 Aug 04 Interactive mode queries whether to restore
Cu   07 Sep 03 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character sopts*(*)
Cgetarg       double precision slat(1),sspec(1),ehf
       real(8):: ehf 
       type(s_lat)::slat
       type(s_spec)::sspec(*)

C ... Local parameters
      logical lrsm,leh,ioorbp,ix !takao ioorpb is a logical function
      integer mxcall,nvmax
      parameter (mxcall=20,nvmax=100)
      integer i,is,j1,j2,j3,nspec,nglob,llst(6),lmxb,nvar,iter,
     .  it,ncall,ivar,ifi,fopna,ir,lsort,lpass,ipr,l,lt,stdo,lgunit
      integer optlst(3,nvmax),optls0(3,nvmax),iwk(nvmax,5)
      character dc*1,slabl*8,strnl*6,pnam(2)*3
      double precision val,rmt,ax,bx,cx,xmin,xmax,fa,fb,fc,fopb
      double precision vhist(2,mxcall),wk(nvmax),estart
      common /lmfop/ vhist,iter,ncall
      external fopb
C ... for lmf basis
      integer nkap0,n0
      parameter (nkap0=3,n0=10)
      integer lh(nkap0),nkapi,ik
      double precision rsmh(n0,nkap0),eh(n0,nkap0),e,rsm,tol,gmax
      save ivar,xmin,xmax,tol,lpass,optlst,estart,gmax
      data ivar /0/ lpass /0/ estart /0d0/


C     call shstru('spec',sspec,1,1)

      pnam(1) = 'rsm'
      pnam(2) = 'eh'
      nspec = nglob('nspec')
      if (sopts .eq. ' ') return
      call getpr(ipr)
      stdo = lgunit(1)

C --- Parse arguments in sopts ---
C     ls = len(sopts)
      j1 = 1
      dc = sopts(j1:j1)
      j1 = j1+1
      nvar = 0

C ... Return here to resume parsing for arguments
   40 continue
      call nwordg(sopts,0,dc//' ',1,j1,j2)
      if (j2 .ge. j1) then
C       print *, sopts(j1:j2)
        if (sopts(j1:j2) .eq. 'atm')  then
          call rx('lmfopb needs implement atm')
        elseif (sopts(j1:j2) .eq. 'sort')  then
          lsort = 1
        elseif (sopts(j1:j1+4) .eq. 'spec=')  then
          call nwordg(sopts,0,dc//', ',1,j1,j3)
          call locase(sopts(j1+5:j3))
          do  is = 1, nspec
Cgetarg             call spacks(0,'spec name',sspec,slabl,is,is)
             do i_spacks=is,is 
             call spacks_copy('u',sspec(i_spacks)%name,is,is,slabl,i_spacks)
             enddo

            call locase(slabl)
C           Found a species: add orbitals to list
            if (slabl .eq. sopts(j1+5:j3)) then
              strnl = '012345'
              lrsm = .false.
              leh  = .false.

C         ... Return here to resume parsing for species arguments
  140         continue
              if (j3 .ge. j2) goto 150
              j1 = j3+2
              call nwordg(sopts,0,','//dc//' ',1,j1,j3)
              if (sopts(j1:j3) .eq. 'rs')  then
                lrsm = .true.
                goto 140
              elseif (sopts(j1:j3) .eq. 'e')  then
                leh = .true.
                goto 140
              elseif (sopts(j1:j1+1) .eq. 'l=')  then
                strnl = sopts(j1+2:j3)
                goto 140
              else
                goto 999
              endif

C         ... Add specification to list
  150         continue
              if (.not. lrsm .and. .not. leh) lrsm = .true.
              call strip(strnl,j1,j3)
              read(strnl,'(6I1)') (llst(j1), j1=1,j3)
Cgetarg...info...           integer :: lmxb , integer(8) :: sspec%lmxb 1
Cgetarg               call upack('spec lmxb',sspec,is,lmxb,0,0,0)
               
               lmxb=sspec(is)%lmxb

C             Exclude states l>lmxb
              do  j1 = 1, j3
                if (llst(j1) .gt. lmxb) llst(j1) = -1
              enddo
              call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkapi)
C             call lmfop2(0,is,sspec,rsmh,eh)
C             Exclude states rsm=0
              ik = 1
              do  j1 = 1, j3
                if (nkapi .eq. 0) llst(j1) = -1
                l = llst(j1)
                if (l .ge. 0) then
                  e = eh(l+1,ik)
                  rsm = rsmh(l+1,ik)
                  if (rsm .le. 0) llst(j1) = -1
                endif
              enddo

C         ... Add to list, smoothing radii
              if (lrsm) then
                do  j1 = 1, j3
                  if (llst(j1) .ge. 0) then
                    nvar = nvar+1
                    optls0(1,nvar) = is
                    optls0(2,nvar) = llst(j1)
                    optls0(3,nvar) = 1
                  endif
                enddo
              endif

C         ... Add to list, smoothing energies
              if (leh) then
                do  j1 = 1, j3
                  if (llst(j1) .ge. 0) then
                    nvar = nvar+1
                    optls0(1,nvar) = is
                    optls0(2,nvar) = llst(j1)
                    optls0(3,nvar) = 2
                  endif
                enddo
              endif

              j1 = j2+2
              goto 40
            endif
          enddo
          call rxs('lmfop: specified nonexistent species: ',
     .      sopts(j1:j3))
        else
          goto 999
        endif
        j1 = j2+2
        goto 40
      endif

C --- Order variables to be minimized ---
      if (lsort .eq. 1 .and. lpass .eq. 0) then
    
C       Order rs variables before e variables
        call icopy(nvar,optls0(3,1),3,iwk,1)
        call ivheap(1,nvar,iwk,iwk(1,2),101)
        call ivprm(3,nvar,optls0,iwk(1,3),iwk(1,2),.true.)
      
C       Count the number of rs variables
        j1 = 0
        do  i = 1, nvar
          if (optls0(3,i) .eq. 1) j1 = i
        enddo
      
C       Put rs variables in ascending order
        do  i = 1, j1
          is = optls0(1,i)
          l  = optls0(2,i)
          lt = optls0(3,i)
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg           call upack('spec rmt',sspec,is,rmt,0,0,0)
           
           rmt=sspec(is)%rmt

Cgetarg           call spacks(0,'spec name',sspec,slabl,is,is)
           do i_spacks=is,is 
           call spacks_copy('u',sspec(i_spacks)%name,is,is,slabl,i_spacks)
           enddo

          call lmfop2(0,is,sspec,rsmh,eh)
          wk(i) = min(rsmh(l+1,ik),rmt)
        enddo
        call dvheap(1,j1,wk,iwk(1,2),1d-3,101)
        call ivprm(3,j1,optls0,iwk(1,3),iwk(1,2),.true.)
      
C       Put e variables in ascending order
        do  i = j1+1, nvar
          is = optls0(1,i)
          l  = optls0(2,i)
          lt = optls0(3,i)
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg           call upack('spec rmt',sspec,is,rmt,0,0,0)
           
           rmt=sspec(is)%rmt

Cgetarg           call spacks(0,'spec name',sspec,slabl,is,is)
           do i_spacks=is,is 
           call spacks_copy('u',sspec(i_spacks)%name,is,is,slabl,i_spacks)
           enddo

          call lmfop2(0,is,sspec,rsmh,eh)
          wk(i-j1) = eh(l+1,ik)
        enddo
        call dvheap(1,nvar-j1,wk,iwk(1,2),1d-3,101)
        call ivprm(3,nvar-j1,optls0(1,j1+1),iwk(1,3),iwk(1,2),.true.)
      endif

C --- First pass: copy sorted table to permanent array ---
      if (lpass .eq. 0) then

C       Keep static record of these variables; reset each new pass
Cgetarg...info...           real(8) :: gmax , real(8) :: slat%gmax 1
Cgetarg         call upack1('lat gmax',slat,gmax)
         
         gmax=slat%gmax


        call icopy(3*nvar,optls0,1,optlst,1)
        lpass = 1
      else
C       Restore these variables
Cgetarg...info...           real(8) :: gmax , real(8) :: slat%gmax 1
Cgetarg         call pack2('lat gmax',slat,gmax,0)
         
         slat%gmax=gmax 

      endif

C --- Minimization procedure for next variable ---
   60 continue
      if (nvar .eq. 0) then
        call info0(10,1,0,
     .    ' LMFOPB (warning): no orbitals specified to optimize')
      endif

C ... Exit when optimization complete
      if (ncall .eq. 0 .and. ivar .eq. nvar) then
        call info2(20,1,0,' LMFOPB:  before optimization ehf=%;4d.'//
     .    '  After optimization ehf=%;4d',estart,vhist(2,iter))
        call rx0('Optimization complete')
      endif

C ... Setup for new variable, or re-entrance after generation of ehf
      if (ncall .eq. 0) then
        vhist(1,1) = 0
        vhist(1,2) = 0
        ivar = ivar+1
      else
        vhist(2,ncall) = ehf
        if (estart .eq. 0) estart = ehf
      endif

C --- Initial Setup ---
      if (ncall .eq. 0 .and. ivar .eq. 1 .and. ipr .ge. 20) then

        call info2(10,1,0,
     .    ' LMFOPB:  optimizing energy wrt %i parameters:%N'//
     .    '   spec%7fl  type     start',nvar,0)
        ik = 1
        do  j1 = 1, nvar
          is = optlst(1,j1)
          l  = optlst(2,j1)
          lt = optlst(3,j1)

Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg           call upack('spec rmt',sspec,is,rmt,0,0,0)
           
           rmt=sspec(is)%rmt

Cgetarg           call spacks(0,'spec name',sspec,slabl,is,is)
           do i_spacks=is,is 
           call spacks_copy('u',sspec(i_spacks)%name,is,is,slabl,i_spacks)
           enddo

C         call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkapi)
          call lmfop2(0,is,sspec,rsmh,eh)
C         Starting guesses; also set limit to range
          if (lt .eq. 1) then
            val = min(rsmh(l+1,ik),rmt)
            xmin = min(1d0,val-.1d0)
            xmax = rmt
          endif
          if (lt .eq. 2) then
            val = eh(l+1,ik)
            xmin = min(-1d0,val-.5d0)
            xmax = -.02d0
          endif
          write(stdo,333) is,slabl,l,pnam(lt),val
  333     format(i4,':',a8,i2,3x,a3,2x,f8.3)
        enddo
      endif

C ... Get parameters for this variable
      is = optlst(1,ivar)
      l  = optlst(2,ivar)
      lt = optlst(3,ivar)
C     call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkapi)
      call lmfop2(0,is,sspec,rsmh,eh)
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg       call upack('spec rmt',sspec,is,rmt,0,0,0)
       
       rmt=sspec(is)%rmt


C ... Starting guesses; also set limit to range
      if (ncall .eq. 0) then
        if (lt .eq. 1) then
          val = min(rsmh(l+1,ik),rmt)
          xmin = min(1d0,val-.1d0)
          xmax = rmt
          tol = .02d0
        endif
        if (lt .eq. 2) then
          val = eh(l+1,ik)
          xmin = min(-1d0,val-.5d0)
          xmax = -.02d0
          tol = .04d0
        endif
        ax = val
        if (val .lt. (xmin+xmax)/2) then
          cx = min(val+.2d0,xmax)
          if (abs(cx-ax) .lt. .1d0) ax = ax-.1d0
        else
          cx = max(val-.2d0,xmin)
          if (abs(cx-ax) .lt. .1d0) ax = ax+.1d0
        endif
        vhist(1,1) = ax
        vhist(1,2) = cx

Cgetarg         call spacks(0,'spec name',sspec,slabl,is,is)
         do i_spacks=is,is 
         call spacks_copy('u',sspec(i_spacks)%name,is,is,slabl,i_spacks)
         enddo

        call info5(20,1,0,' LMFOPB:  begin optimization of var #%i in'//
     .    ' species '//slabl//'%a: start=%d range=(%d,%d)',ivar,ax,
     .    xmin,xmax,0)
      else
Cgetarg         call spacks(0,'spec name',sspec,slabl,is,is)
         do i_spacks=is,is 
         call spacks_copy('u',sspec(i_spacks)%name,is,is,slabl,i_spacks)
         enddo

        call info2(20,1,0,' LMFOPB: continue optimization of var #%i,'//
     .    ' species '//slabl//'%a val=%d',ivar,vhist(1,ncall))
      endif

C ... Go through minimization steps until run out of energies, or converged
      iter = 0
      ax = vhist(1,1)
      cx = vhist(1,2)
      call mnbrak(ax,bx,cx,fa,fb,fc,fopb,xmin,xmax,1,mxcall,it,ir)
      if (ir .eq. -2) vhist(1,iter) = cx
      if (lt .eq. 1) then
        rsmh(l+1,ik) = vhist(1,iter)
      endif
      if (lt .eq. 2) then
        eh(l+1,ik) = vhist(1,iter)
      endif
      call lmfop2(1,is,sspec,rsmh,eh)
      if (ir .eq. -2) then
        call info5(20,1,0,' LMFOPB: var #%i min appears outside'//
     .    ' range (%d,%d) ... use %d : ehf=%;6d',ivar,xmin,xmax,
     .    vhist(1,iter),vhist(2,iter))
        ncall = 0
        ifi = fopna('basp',-1,0)
        rewind ifi
        ix= ioorbp(1,2,1,nspec,sspec,0,-ifi)
        call fclr(' ',ifi)
        j1 = 1
        call query('enter <ret> or 0<ret> to '//
     .    ' continue or restore starting value before continuing',2,j1)
        if (j1 .eq. 0) then
          if (lt .eq. 1) then
            rsmh(l+1,ik) = vhist(1,1)
          endif
          if (lt .eq. 2) then
            eh(l+1,ik) = vhist(1,1)
          endif
          call lmfop2(1,is,sspec,rsmh,eh)
          ifi = fopna('basp',-1,0)
          rewind ifi
          ix= ioorbp(1,2,1,nspec,sspec,0,-ifi)
          call fclr(' ',ifi)
        endif
        goto 60
      endif
C ... mnbrak wants another function call
C     if (iter .ge. ncall) then
      if (ir .eq. -1) then
        ncall = ncall+1
        return
      endif

C ... Minimize
      call brent(ax,bx,cx,fa,fb,fc,fopb,1,tol,.01d0,0,ir)

      if (ir .ge. 0) vhist(1,iter) = bx
      if (ir .ge. 0) vhist(2,iter) = fb

      if (lt .eq. 1) then
        rsmh(l+1,ik) = vhist(1,iter)
      endif
      if (lt .eq. 2) then
        eh(l+1,ik) = vhist(1,iter)
      endif
      call lmfop2(1,is,sspec,rsmh,eh)

      if (ir .ge. 0) then
        call info5(20,1,0,' LMFOPB: var #%i converged in %i iterations'
     .    //' to %;3d:  ehf=%;6d ... writing file basp',
     .    ivar,iter,bx,fb,0)
        ifi = fopna('basp',-1,0)
        rewind ifi
        ix= ioorbp(1,2,1,nspec,sspec,0,-ifi)
        call fclr(' ',ifi)
        j1 = 1
        call query('enter <ret> or 0<ret> to'//
     .    ' continue or restore starting value before continuing',2,j1)
        if (j1 .eq. 0) then
          if (lt .eq. 1) then
            rsmh(l+1,ik) = vhist(1,1)
          endif
          if (lt .eq. 2) then
            eh(l+1,ik) = vhist(1,1)
          endif
          call lmfop2(1,is,sspec,rsmh,eh)
          ifi = fopna('basp',-1,0)
          rewind ifi
          ix= ioorbp(1,2,1,nspec,sspec,0,-ifi)
          call fclr(' ',ifi)
        endif
      else
        call info2(20,1,0,' LMFOPB:  found=F   x=%d',vhist(1,iter),0)
      endif

C ... brent wants another function call
      if (ir .eq. -1) then
        ncall = ncall+1
        return
      endif

C     Done with this variable; start another
      ncall = 0
      goto 60

C     return
  999 call rxs('lmfopb: failed to parse switch: ',sopts(j1:j2))
Cgetarg       end
       end subroutine lmfopb 

      subroutine lmfop2(lpack,is,sspec,rsmh,eh)
       
       use m_struc_def  !Cgetarg

C- Pack or unpack species parameters
C     implicit none
      integer lpack,is
Cgetarg       double precision sspec(1)
       type(s_spec)::sspec(*)

      integer nkap0,n0
      parameter (nkap0=3,n0=10)
C     integer lh(nkap0),nkapi,npqni,ik
      double precision rsmh(n0,nkap0),eh(n0,nkap0)
C Local
      integer ikap
      double precision orbp(n0,2,nkap0)

      call dpzero(orbp,n0*nkap0*2)
Cgetarg...info...           real(8) :: orbp(n0,2,nkap0) , real(8) :: sspec%orbp 60
Cgetarg       call upack2('spec orbp',sspec,is,orbp)
       
       i_copy_size=size(sspec(is)%orbp) 
       call dcopy(i_copy_size,sspec(is)%orbp,1,orbp,1) 

C ... Unpack
      if (lpack .eq. 0) then
C       call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkapi)
        ikap = 1
        call dcopy(n0,orbp(1,1,ikap),1,rsmh(1,ikap),1)
        call dcopy(n0,orbp(1,2,ikap),1,eh(1,ikap),1)
      else
        ikap = 1
        call dcopy(n0,rsmh(1,ikap),1,orbp(1,1,ikap),1)
        call dcopy(n0,eh(1,ikap),1,orbp(1,2,ikap),1)
Cgetarg...info...           real(8) :: orbp(n0,2,nkap0) , real(8) :: sspec%orbp 60
Cgetarg         call pack2('spec orbp',sspec,is,orbp)
         
         i_copy_size=size(sspec(is)%orbp) 
         call dcopy(i_copy_size,orbp,1,sspec(is)%orbp,1) 

      endif
Cgetarg       end
       end subroutine lmfop2 


      double precision function fopb(x,nfcn0,ir)

      double precision x(1)
      integer nfcn0,ir

      integer mxcall,ncall,iter
      parameter (mxcall=20)
      double precision vhist(2,mxcall)
      common /lmfop/ vhist,iter,ncall

      iter = iter+1
      if (iter .le. ncall) then
        if (x(1) .ne. vhist(1,iter))
     .    call rxi('lmfopb: history corrupted iter %i .. aborting',iter)
        ir = iter
        fopb = vhist(2,iter)
      else
        ir = -1
        if (iter .gt. mxcall)
     .    call rxi('lmfopb: failed to converge in %i calls',mxcall)
        vhist(1,iter) = x(1)
        fopb = 0
      endif
Cgetarg       end
       end function fopb 


