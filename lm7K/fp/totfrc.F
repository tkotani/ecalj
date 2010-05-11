Cdelw1 var ['owk', '1', 'owk', '3*nbas', 'real(8)', 'rv_w_', 'totfrc']
Cdelw1 var ['og', '-1', 'undef', '0', 'undef', 'x', 'totfrc']
Cdelw1 var ['oistab', '-1', 'undef', '0', 'undef', 'x', 'totfrc']
Cdelw1 undel [['undef', 'wref'], 'og', 'totfrc']
Cdelw1 undel [['undef', 'wref'], 'oistab', 'totfrc']
Cgetarg...info...           structure ['totfrc', 'site', 'ssite']
Cgetarg...info...           structure ['totfrc', 'lat', 'slat']
Cgetarg...info...           use_to_add ['totfrc', 'm_struc_def']
      subroutine totfrc(nbas,ssite,slat,leks,fes1,fes2,dfhf,f)
       
       use m_struc_def  !Cgetarg

C- Add together and print contributions to force
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ssite :struct containing site-specific information
Ci   leks: :<0 no forces calculated
Ci         : 0 Harris forces only
Ci         :>0 also print HKS forces
Ci         :>1 use HKS forces instead of HF forces
Ci   fes1  :contribution to HF forces from estat + xc potential
Ci   fes2  :contribution to KS forces from estat + xc potential
Ci   dfhf  :2nd order corr. HF forces from ansatz density shift (dfrce)
Cio Inputs/Outputs
Cio     f  :On input, f is the contribution to force from eigval sum
Cio        :On output, f is the total force
Cr Remarks
Cu Updates
Cu   12 Apr 03 Prints out max correction to Harris force
Cu   30 May 00 Adapted from nfp totforce.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,leks
Cgetarg       double precision ssite(*),slat(*),f(3,nbas),fes1(3,nbas),
Cgetarg      .  fes2(3,nbas),dfhf(3,nbas)
       real(8):: f(3,nbas) , fes1(3,nbas) , fes2(3,nbas) , dfhf(3,nbas) 
       type(s_site)::ssite(*)
       type(s_lat)::slat

C ... Local parameters
Cdelw1 do not change  [['undef', 'wref'], 'oistab', 'totfrc']
Cdelw1 do not change  [['undef', 'wref'], 'og', 'totfrc']
Cdelw1       integer stdo,stdl,lgunit,ipr,ipl,ib,m,owk,ng,oistab,og,ibmx,isw
       integer:: stdo , stdl , lgunit , ipr , ipl , ib , m , ng , oistab 
     . , og , ibmx , isw 
       real(8) ,allocatable :: rv_w_owk(:)

      double precision fev1(3),fev2(3),fhar(3),fks(3),c,ddot,fmax,dfmax
C ... Heap
      integer w(1)
      common /w/ w

      stdo = lgunit(1)
      stdl = lgunit(2)
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [og oistab]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (leks .lt. 0) return
      c = 1000d0
      call getpr(ipr)
      ipl = 1
      fmax = -1
      dfmax = -1

      if (ipr .ge. 30) then
        if (ddot(3*nbas,dfhf,1,dfhf,1) .ne. 0) then
          write(stdo,'(/'' Forces, with eigenvalue correction'')')
        else
          write(stdo,'(/''Forces:'')')
        endif
        write(stdo,201)
  201   format('  ib',11x,'estatic',18x,'eigval',20x,'total')
      endif
      do  ib = 1, nbas
        do  m = 1, 3
          fev1(m) = f(m,ib) + dfhf(m,ib)
          fev2(m) = f(m,ib)
          fhar(m) = fev1(m) + fes1(m,ib)
          f(m,ib) = fhar(m)
        enddo
        if (dsqrt(ddot(3,fhar,1,fhar,1)) .gt. fmax) then
          ibmx = ib
          fmax = dsqrt(ddot(3,fhar,1,fhar,1))
        endif
        if (dsqrt(ddot(3,dfhf(1,ib),1,dfhf(1,ib),1)) .gt. dfmax) then
          ibmx = ib
          dfmax = dsqrt(ddot(3,dfhf(1,ib),1,dfhf(1,ib),1))
        endif
        if (ipr .ge. 30) write(stdo,200) ib,(c*fes1(m,ib),m=1,3),
     .    (c*fev1(m),m=1,3),(c*fhar(m),m=1,3)
  200   format(i4,3f8.2,1x,3f8.2,1x,3f8.2)
        if (leks. eq. 0 .and. ipl .gt. 0 .and. ipr .ge. 30)
     .      write(stdl,710) ib,(c*fhar(m),m=1,3)
Cgetarg...info...           real(8) :: f(3,nbas) , real(8) :: ssite%force 3
Cgetarg         call pack2('site force',ssite,ib,f(1,ib))
         
         i_copy_size=size(ssite(ib)%force) 
         call dcopy(i_copy_size,f(1,ib),1,ssite(ib)%force,1) 

  710   format('fp ib',i4,'  fh ',3f8.2,2x,3f8.2)

        if (leks .ge. 1) then
          do  m = 1, 3
            fks(m)  = fev2(m) + fes2(m,ib)
            if (leks .gt. 1) f(m,ib) = fks(m)
Cgetarg...info...           real(8) :: fks(3) , real(8) :: ssite%force 3
Cgetarg             if (leks .gt. 1) call pack2('site force',ssite,ib,fks)
             if ( leks .gt. 1 ) then 
             
             i_copy_size=size(ssite(ib)%force) 
             call dcopy(i_copy_size,fks,1,ssite(ib)%force,1) 
             endif 

          enddo
          if (ipr .gt. 40) write(stdo,210) (c*fes2(m,ib),m=1,3),
     .       (c*fev2(m),m=1,3),(c*fks(m),m=1,3)
  210     format('  KS',3f8.2,1x,3f8.2,1x,3f8.2)
          if (ipl .gt. 0 .and. ipr .ge.30)
     .       write (stdl,711) ib,(c*fhar(m),m=1,3),(c*fks(m),m=1,3)
  711     format('fp ib',i4,'  fh ',3f8.2,'   fks ',3f8.2)
        endif
      enddo
      call info5(10,0,0,' Maximum Harris force = %;3g mRy/au (site %i)'
     .  //'%?#n#  Max eval correction = %;1d##',
     .  c*fmax,ibmx,isw(dfmax.gt.0),c*dfmax,0)

C     Symmetrize forces to machine precision
Cgetarg...info...           integer :: ng , integer(8) :: slat%nsgrp 1
Cgetarg...info...           integer :: oistab , integer(8) :: slat%oistab 1
Cgetarg...info...           integer :: og , integer(8) :: slat%osymgr 1
Cgetarg       call upack('lat nsgrp oistab osymgr',slat,ng,oistab,og,0,0)
       
       ng=slat%nsgrp
Cdelw1 warning(1) , probably  oistab  is not defined yet at linenumber= 122
Cdelw1 do not change oistab because of ['undef', 'wref']
       oistab=slat%oistab
Cdelw1 warning(1) , probably  og  is not defined yet at linenumber= 123
Cdelw1 do not change og because of ['undef', 'wref']
       og=slat%osymgr

      if (ng .gt. 1) then
      call info(30,1,0,' Symmetrize forces ...',0,0)
Cdelw1       call defdr(owk,3*nbas)
       allocate(rv_w_owk(3*nbas))
Ctakao_ZeroClear_NotRequiered        if (3*nbas<0) rv_w_owk(:)=0.0d0

Cdelw1 warning(1) , probably  og  is not defined yet at linenumber= 128
Cdelw1 do not change og because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oistab  is not defined yet at linenumber= 128
Cdelw1 do not change oistab because of ['undef', 'wref']
Cdelw1 warning(1) , probably  og  is not defined yet at linenumber= 128
Cdelw1 do not change og because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oistab  is not defined yet at linenumber= 128
Cdelw1 do not change oistab because of ['undef', 'wref']
Cdelw1       call symfor(nbas,1,w(og),ng,w(oistab),w(owk),f)
       call symfor ( nbas , 1 , w ( og ) , ng , w ( oistab ) , rv_w_owk 
     . , f ) 

Cdelw1 rlse name= owk old_list= owk 
Cdelw1 rlse name= owk new_list= (None)
Cdelw1       call rlse(owk)
       if (allocated(rv_w_owk)) deallocate(rv_w_owk)

      endif

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [og oistab]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine totfrc 


