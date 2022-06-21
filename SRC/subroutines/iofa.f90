integer function iofa(spid,nxi0,nxi,exi,hfc,hfct,rsm,z,rmt,a,nr, &
     qc,ccof,ceh,stc,rho,rhoc,v,ifi)
  use m_lmfinit,only: nsp,lrel
  use m_lgunit,only:stdo
  !- I/O for free-atom data, one species
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   spid  :species label
  !i   nxi0  :leading dimension of hfc,hfct
  !i   ifi   :file logical unit, but >0 for read, <0 for write
  ! o File I/O
  ! o  nxi   :number of energies used to fit tail of valence density
  ! o  exi   :energies that fit tail of valence density
  ! o  hfc   :coefficients for fit of tail valence density
  ! o  hfct  :coefficients for fit of tail valence+core density
  ! o  rsm   :smoothing radius for fit of tail valence density
  ! o  z     :nuclear charge
  ! o  rmt   :muffin-tin radius, in a.u.
  ! o  a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
  ! o  nr    :number of radial mesh points
  ! o  qc    :Sphere core charge
  ! o  ccof  :coefficient to core density fit by unsmoothed Hankel
  ! o  ceh   :hankel function energy to core density fit
  ! o  stc   :core kinetic energy
  ! o  rho   :valence density
  ! o  rhoc  :core density
  ! o  v     :spherical potential
  !r Remarks
  !u Updates
  !u   10 Jun 00 spin polarized
  !u   20 May 00 adapted from nfp rw_fa.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: ifi,nr,nxi,nxi0
  character spid*8
  double precision :: a,ccof,ceh,qc,rmt,rsm,stc,z, &
       hfc(nxi0,2),rho(2*nr),v(2*nr),exi(nxi0),rhoc(2*nr),hfct(nxi0,2)
  integer :: i,jfi,lrel0,nglob,nsp0,ipr,iprint
  logical :: isanrg
  character msg*23
  ipr    = iprint()
  msg  = '         File mismatch:'
  iofa   = -1
  ! --- Input ---
  if (ifi > 0) then
     jfi = ifi
     read(jfi,201,end=998,err=998) spid
201  format(9x,a8/)
     read(jfi,102) z,a,nsp0,lrel0,nr,rmt,rsm
     if (isanrg(lrel0, lrel,lrel,msg,'lrel', .TRUE. )) stop
     if (isanrg(nsp0,  nsp,nsp,  msg,'nsp', .TRUE. )) stop
     read(jfi,103) nxi
     read(jfi,205) (exi(i),i=1,nxi)
     read(jfi,205) (hfc(i,1),i=1,nxi)
     read(jfi,205) (hfct(i,1),i=1,nxi)
205  format(5d16.8)
     if (nsp == 2) read(jfi,*) (hfc(i,2),i=1,nxi)
     if (nsp == 2) read(jfi,*) (hfct(i,2),i=1,nxi)
     read (jfi,210) qc,ccof,ceh,stc
210  format(5x,4f16.7)
     read(ifi,333) rho(1:nr)
     if (nsp == 2) read(ifi,333) rho(nr+1:2*nr)
     read(ifi,333) rhoc(1:nr)
     if (nsp == 2) read(ifi,333) rhoc(nr+1:2*nr)
     read(ifi,333) v(1:nr)
     if (nsp == 2) read(ifi,333) v(nr+1:2*nr)
  endif
  ! --- Output ---
  if (ifi < 0)  then
     jfi = -ifi
     write(jfi,101) spid
101  format('-------- ',a8,' ----------'/ &
          '     z       a      nsp lrel  nr     rmt          rsm')
     write(jfi,102) z,a,nsp,lrel,nr,rmt,rsm
102  format(f8.2,f9.4,3i5,2f13.7)
     write(jfi,103) nxi
103  format(i4)
     write(jfi,105) (exi(i),i=1,nxi)
     write(jfi,105) (hfc(i,1),i=1,nxi)
     write(jfi,105) (hfct(i,1),i=1,nxi)
     if (nsp == 2) write(jfi,105) (hfc(i,2),i=1,nxi)
     if (nsp == 2) write(jfi,105) (hfct(i,2),i=1,nxi)
105  format(1p,5d16.8)
     write (jfi,110) qc,ccof,ceh,stc
110  format(' core',4f16.7)
     write(jfi,333) rho(1:nr) !call dfdump(rho,nr,-jfi)
     if (nsp == 2) write(jfi,333) rho(nr+1:2*nr) !call dfdump(rho(1+nr),nr,-jfi)
     write(jfi,333) rhoc(1:nr) !call dfdump(rhoc,nr,-jfi)
     if (nsp == 2) write(jfi,333) rhoc(nr+1:2*nr) !call dfdump(rhoc(1+nr),nr,-jfi)
     write(jfi,333) v(1:nr)! call dfdump(v,nr,-jfi)
     if (nsp == 2) write(jfi,333) v(nr+1:2*nr) !call dfdump(v(1+nr),nr,-jfi)
     return
  endif
  iofa = 0
  return
998 if(ipr > 0) write(stdo,'('' iofa  : missing species id ... nothing read'')')
333 format(1p,4e20.13)
end function iofa
