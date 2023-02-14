module m_gtv2
  use m_ftox
  use m_lgunit,only:stdo
  use m_MPItk,only:master_mpi
  public rval2,gtv2_setrcd
  private
    integer::nrecs,reclnr
    character(:),allocatable:: recrd(:)
contains
  subroutine gtv2_setrcd(recrdin)
    character(len=*):: recrdin(:)
    nrecs = size(recrdin)
    reclnr= len(recrdin(1))   !write(6,*)nrecs,reclnr
    allocate(character(reclnr)::recrd(nrecs))
    recrd=recrdin
  end subroutine gtv2_setrcd
  subroutine rval2(cattok, ch,rr,rv, nout,nreq,defa) !Choose one of nout,nreq,defa(default), or none
      integer:: NULL=-99999                          !Choose one of ch,rr,rv
      integer,optional:: nout,nreq          
      real(8),optional:: defa(:),rr
      real(8),optional,allocatable:: rv(:)
      character(*),optional:: ch
      character(*):: cattok
      integer::ncat,i,ncount,ndefa,lx
      real(8):: arr(1000),rvx(1000)
      logical:: nomode,nrmode,ndmode,debug=.false.
      character(4):: modec
      character(8):: xn
      character(256):: outx
      if(debug) write(stdo,*)'cccccccc rval2 ',trim(cattok)
      ncat=len(trim(cattok))
      if(present(ch)) then
         do i=1,nrecs
            if( recrd(i)(1:ncat)==trim(cattok) ) then
               read(recrd(i)(ncat+1:),"(a)") ch
               goto 1012
            endif
         enddo
         ch=''
1012     continue
         if(debug) write(stdo,*)'cccccccc ',trim(cattok),'ch=###'//trim(ch)//'###'
         return
      endif   
      nomode=.false.; nrmode=.false.; ndmode=.false.
      if(present(nout)) nomode=.True. !nout mode
      if(present(nreq)) nrmode=.True. !nreq mode
      if(present(defa)) ndmode=.True. !default mode
      if(count([nomode,nrmode,ndmode])>1) call rx('rval2: chooose one of nout nreq default or none')
      if(ndmode) ndefa =size(defa)
      do i=1,nrecs
         if( recrd(i)(1:ncat+1)==trim(cattok)//' ' ) then
!            print *,'BZ_N xxxxx',recrd(i)(1:ncat), 'xxxx', trim(recrd(i)(ncat+1:))
            call getdval(trim(recrd(i)(ncat+2:)),ncount,arr) !Read undefinit number of real(8) array
            exit
         endif
      enddo
      if(ndmode)then !default mode
         if(ncount>ndefa) ncount=ndefa !truncation
         if(ncount>0.and.ncount/=ndefa)call rx('rval2: '//trim(cattok)//' default mode. ncount='//xn(ncount)//'/=ndefa='//xn(ndefa))
      elseif(nrmode) then !nrequest mode
         if(ncount>nreq) ncount=nreq !truncation
         if(ncount/=nreq) call rx('rval2: '//trim(cattok)//' nreq mode. ncount='//xn(ncount)//'/=nreq='//xn(nreq))
      endif
      
      if(present(rv)) then !rv mode
        if(allocated(rv)) deallocate(rv)
        if(ncount==0.and.ndmode) then
           ncount=ndefa
           allocate(rv(ncount))
           rv(1:ncount)=defa(1:ncount)
        else   
           allocate(rv(ncount))
           rv(1:ncount)=arr(1:ncount)
        endif
        rvx(1:ncount)=rv(1:ncount)
      else !rr mode
        if(.not.present(rr)) call rx('rval2:neither rr nor rv given')
        if(ncount==0.and.ndmode) then
           ncount=1
           rr=defa(1)
        elseif(ncount==1) then   
           rr=arr(1)
        else
           rr=NULL
        endif
        rvx(1)=rr
      endif
      if(present(nout)) nout=ncount !nomode
!print
      modec='----'
      if(nomode) modec='----'
      if(nrmode) modec='requ'
      if(ndmode) modec='defa'
      if(debug.and.present(rv)) write(stdo,*)'cccccccrr rval2 mode ncount val ',&
           trim(cattok),nomode,nrmode,ndmode,ncount,rv(1:ncount)
      if(debug.and.present(rr)) write(stdo,*)'cccccccrv rval2 mode ncount val ',&
           trim(cattok),nomode,nrmode,ndmode,ncount,rr
      outx='rval2: '//trim(cattok)
      lx=len_trim(outx)
      if(master_mpi) write(stdo,ftox) trim(outx)//repeat(' ',30-lx),trim(modec),'n=',ncount,'val=',ftof(rvx(1:ncount),8)
 end subroutine
end module

logical function parmxp(strnin,lstrn,broy,nmix,wgt,beta,wc,killj)
  implicit none
  integer :: broy
  character strnin*(*)
  integer :: iter,lstrn,nmix
  double precision :: wgt(3),beta,wc,betv,rmserr,akillj !elind,
  integer :: i,j,np,it(5),parg,nit,nmixj,jp,kp,killj,nitj, &
       iprint,i1mach,a2vec,lstblk,lstitj,iblk,nbump,k,ia0,lbroy,iterm
  logical :: lpr,lagain,cmdopt
  character outs*100,fnam*8,num*10,strn*1000
  double precision :: bet,elin,wt(3),wcj,rmsc,errmin,xx
  parmxp = .true.
  if (strnin == ' ' .OR. lstrn <= 0) goto 9999
  lbroy = broy
  wcj = wc
  nmixj = nmix
  bet = beta
  wt = wgt
  if (wgt(3) == -9) wt(3) = 0
  strn=adjustl(strnin)
  np =0
  jp =0
  iblk = 1
  call chrps2(strn,'AaBbCc',6,np,jp,it)!Broyden or Anderson mixing
  if (it(1) == 0) goto 999
  lbroy = 0
  if (it(1) >= 3) lbroy = 1
  if (it(1) >= 5) lbroy = 2
  if (lbroy == 1) outs =' mixrho: mixing mode=A' ! call awrit0('%a%bB',outs,len(outs),0)
  if (lbroy == 2) outs =' mixrho: mixing mode=B' !call awrit0('%a%bC',outs,len(outs),0)
  ! ... Pick up nmix
  iterm = index(trim(strn),',')
  if(iterm>2) read(strn(2:iterm-1),*) nmixj
  iterm = index(trim(strn),' ')
  if(iterm>2) read(strn(2:iterm-1),*) nmixj
  ! ... Pick up rmsc
  rmsc = -1
  jp = np
  i = parg(',r<',strn,jp,lstrn,',; ',2,1,it,rmsc)
  if (i < 0) goto 999
  ! ... Pick up mixing wc
  jp = np
  if (lbroy == 1) then
     i = parg(',wc=',strn,jp,lstrn,',; ',2,1,it,wcj)
     if (i < 0) goto 999
  endif
  ! ... Pick up mixing beta
  jp = np
  i = parg(',b=',strn,jp,lstrn,',; ',2,1,it,bet)
  if (i < 0) goto 999
  ! ... Pick up weights
  jp = np
  i = parg(',w=',strn,jp,lstrn,',; ',2,2,it,wt)
  if (i < 0) goto 999
  jp = np
  j = parg(',wa=',strn,jp,lstrn,',; ',2,1,it,wt(3))
  if (j < 0) goto 999
  !...  Pick up iteration number for file kill
  akillj = -1
  jp = np
  i = parg(',k=',strn,jp,lstrn,',; ',2,1,it,akillj)
  killj=akillj
  if (i < 0) goto 999
  broy = lbroy
  nmix = nmixj
  wgt(1) = wt(1)
  wgt(2) = wt(2)
  if (wgt(3) == -9) wt(3) = 0
  wgt(3) = wt(3)
  beta = bet
  wc = wcj
  goto 9999
999 outs = 'parmxp: parse failed:'//strn(1:lstrn) !error exit
  parmxp = .false.
9999 continue ! --- Normal exit ---
end function parmxp
integer function parg(tok,strn,ip,lstr,sep,itrm,narg,it,res)
  use m_MPItk,only:procid
  !- Returns res= real(8) values from a string
  implicit none
  integer :: lstr,ip,cast,narg,itrm,it(1)
  character*(*) tok,sep,strn
  double precision :: res(narg)
  logical :: ldum,parstr
  character term*1
  character(100)::ddd
  integer :: jp,np,nsep,lentok,a2vec,ipx
  integer::itx(100),ixx
  character(8):: xt
  nsep = len(sep)
  lentok = len(tok)
  term = tok(lentok:lentok)
  ! --- Find end of string ---
  jp = ip
  it(1) = 0
  if (itrm <= nsep) call chrps2(strn,sep(itrm:nsep),nsep-itrm+1,lstr,jp,it)
  if (it(1) /= 0) then
     np = jp
  else
     np = lstr
  endif
  if (tok /= ' ') then
     if (narg == 0 .AND. np == lentok) then
        ip = np+1
        parg = 0
        if (strn(1:np) == tok(1:np)) parg = 1
        return
     elseif ( .NOT. parstr(strn,tok,np-lentok,lentok,term,ip,jp)) then
        parg = 0
        ip = np+1
        return
     endif
  else
     jp = ip
  endif
  ! --- Parse for vector of binary values to convert
  if (narg == 0) then
     ip = jp
     parg = 1
     return
  endif
  ip = jp
  cast=4
  ddd=trim(strn(ip+1:))
  call getdval(ddd, parg, res)
end function parg
