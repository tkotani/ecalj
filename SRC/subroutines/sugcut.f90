module m_sugcut
  public sugcut
  integer,allocatable,protected,public :: ngcut(:,:,:)
  private
contains  
subroutine sugcut(mode) 
  use m_lmfinit,only: nspec,alat=>lat_alat,tol=>lat_tolft,n0,nkap0,nkaphh,lhh,nkapii,sspec=>v_sspec,slabl
  use m_supot,only: gv=>rv_a_ogv,ng=>lat_ng
  use m_uspecb,only:uspecb
  use m_lgunit,only:stdo
  !- Find max recip for each spec and orbital block, store in struct.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1 make cutoffs for standard envelope functions
  !i         :2 make cutoffs for extended local orbitals
  !i         :3 combination 1+2
  !o Outputs
  !o    sspec%ngcut
  !u Updates
  !u   16 Aug 04 New mode for getting cutoffs, local orbs.
  !u             Changed argument list
  !u   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
  !u    9 May 00 Adapted from nfp su_gvcut.f
  ! ----------------------------------------------------------------------
  implicit none
  integer:: mode,lh(nkap0) !, ngcut(n0,nkap0)
  integer:: ipr,iprint,is,irep,icut,i,ik,l,lcut,nkapi,nkap1,nkap2
  real(8):: rsmh(n0,nkap0),eh(n0,nkap0),tpi,tpiba2,gg0,gg,e,rsm,gam,gmax,top
  character(8) :: spid
  character(1) :: ccc,ccl
  if (ng == 0) return
  if(mode==1) allocate(ngcut(n0,nkap0,nspec) )
  ipr = iprint()
  tpi = 8d0*datan(1d0)
  tpiba2 = (tpi/alat)**2
  if (ipr >= 20) then
     if (mode == 1) write(stdo,887) tol
     if (mode == 2) write(stdo,888) tol
     write(stdo,774)
887  format(/' sugcut:  make orbital-dependent reciprocal vector',' cutoffs for tol=',1p,e9.2)
888  format(/' sugcut:  orbital-dependent cutoffs for local',' orbitals, tol=',1p,e9.2)
  endif
  gg = -1
  if(mode==1) ngcut=0
  do  is = 1, nspec
     spid = slabl(is) !sspec(is)%name
     nkap1 = 1
     call uspecb(is,rsmh,eh)
     nkap2 = nkaphh(is)
     nkapi = nkapii(is)
     if (mode>1) then
        if(mode==2) nkap1 = nkapi+1
     else
        nkap2 = nkapi
     endif
!     if(mode==2 ) ngcut=sspec(is)%ngcut
     gg0 = gg
     do  ik = nkap1, nkap2
        lcut = -1
        do  l  = 0, lhh(ik,is)
           e = eh(l+1,ik)
           rsm = rsmh(l+1,ik)
           if (rsm /= 0) then
              if (l < lhh(ik,is) .AND. l > lcut) then
                 lcut = l-1
12               lcut = lcut+1
                 if (lcut < lhh(ik,is)) then
                    if (rsmh(lcut+2,ik) == rsm .AND. eh(lcut+2,ik) == e) goto 12
                 endif
              endif
              !     ... Get cutoff radius where exp(-gam*gmax)*gmax**l equals tol
              gam = rsm*rsm/4d0
              gmax = 1d0
              do  irep = 1, 10
                 gmax = dsqrt(-dlog(tol/gmax**l)/gam)
              enddo
              !     ... Find first longer vector, icut is one less
              icut = 1
              do  i = 1, ng
                 gg = tpiba2*(gv(i,1)**2+gv(i,2)**2+gv(i,3)**2)
                 if (gg > gmax*gmax) goto 90
                 icut = i
                 gg0 = gg
              enddo
90            continue
              top = dexp(-gam*gg0)*dsqrt(gg0)**l
              ccc = ' '
              if (icut == ng) ccc = '*'
              ccl = ' '
              if (l < lcut) ccl = '*'
              if (ipr >= 20) write(stdo,773) spid,l,ccl,rsm,e,gmax,top,icut,ccc
773           format(2x,a,i2,a1,f7.2,f7.2,f8.3,1p,e12.2,0p,i8,a)
774           format(' spec      l    rsm',4x,'eh',5x,'gmax',4x,'last term',4x,'cutoff')
              ngcut(l+1,ik,is) = icut
           endif
        enddo
     enddo
!     sspec(is)%ngcut=ngcut
  enddo
end subroutine sugcut

end module m_sugcut
