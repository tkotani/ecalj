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
