subroutine efsimplef2ax ( legas, esmr, valn,ef)
  use m_READ_BZDATA,only: nqbz,nqibz,ginv,  qibz,wibz,qbz
  use m_genallcf_v3,only: nspin,z,natom,nclass,iclass,nl,konfig=>konf
  use m_readeigen, only: readeval
  use m_readhbe,only: nband
  use m_mpi, only: mpi__root
  use m_hamindex,only:   zbak
  implicit none
  intent(in)::            legas 
  intent(inout)::                esmr, valn !in for legas=T
  intent(out)::                             ef
  !!== Calculate efermi for discrete sum. (not for tetrahedron method) ==
  !! You need to call init_reaeigen before you call this.
  !! user readeval (readeigen.f) to get eigenvalues.
  ! nspin   = 1, paramagnetic
  !           2, ferromagnetic
  ! ef      = fermi level
  ! nband   = no. states
  ! nqbz    = no. k-points
  ! valn    = number of valence electron.
  ! -------------------
  !      e(iband) < efm : occupation is one
  ! efm< e(iband) < efp : occupation is wfacef.
  ! efp< e(iband)       : occupation is zero
  integer:: is,iq,ik,isig,kpx,ifev(2), ieaord(nband*nqibz*nspin),mbytes,mwords,iwksize, iqibz
  real(8)   :: ekt(nband, nqibz,nspin), ektx(nband*nqibz*nspin)
  real(8)   :: wgt(nband, nqibz,nspin), wgtx(nband*nqibz*nspin)
  real(8)   :: qx(3),qbas(3,3),wwg !qbzx(3), ,ginv(3,3),
  integer:: ncore,l,ia,ic ,ierr
  real(8)   :: valn,ef
  integer :: nbnqnsp,ix,ikx=-9999,ikini,nne
  real(8)    :: ew1,ew2,ein,valx,enumef_gauss,esmr, efini ,eee2,wwg2 ,enumef
  logical :: legas,autoew,GaussSmear=.true. !is external
  integer:: if8301,if8302 !nqbz,
  autoew =.false.
  if(GaussSmear) then; write(6,*)' efsimplef2(gaussian mode):start'
  else;                write(6,*)' efsimplef2(rectangular mode):start';  endif
  if(esmr<=0d0) autoew= .TRUE. 
  if(legas) then
     write(6,*)' efsimplef2: legas=T use given valn = ',valn
  else ! total valence charge
     valn    = 0d0
     do ia   = 1,natom
        ic    = iclass(ia)
        valn  = valn + z(ic)
        write(6,*)' ia z(ic)=',ia, z(ic)
        do    l = 0,nl-1
           write(6,*)' l (konfig(l+1,ic)-l-1) 2*(2l+1)=',l,(konfig(l+1,ic)-l-1),( 2*l +1)*2
           valn  = valn - (konfig(l+1,ic)-l-1) *( 2*l +1)*2
        end do
     end do
  endif
  valn=valn-zbak !subtruct back ground charge
  do is = 1,nspin
     do iq = 1,nqibz
        ekt(:,iq,is) = readeval(qibz(:,iq),is)
     enddo
  enddo
  if(abs(sum(wibz(1:nqibz))-2d0)>1d-10) then
     write(6,*) 'sum (wibz)=', sum(wibz(1:nqibz))
     call rx( 'efsimplef2: wibzsumerr')
  endif
  do is = 1,nspin
     do iq = 1,nqibz
        wgt(1:nband,iq,is) = wibz(iq)
        if(nspin==2) wgt(1:nband,iq,is) = wgt(1:nband,iq,is)/2d0
     enddo
  enddo
  ! ekt and wgt
  call dcopy ( nband*nqibz*nspin, ekt,1, ektx ,1)
  call dcopy ( nband*nqibz*nspin, wgt,1, wgtx ,1)
  call sortea( ektx,ieaord,nband*nqibz*nspin,isig)
  if (mpi__root) then
     open(newunit=if8301,file = "DOSACC.lda")
     open(newunit=if8302,file = "DOSACC2.lda")
  endif
  wwg = 0d0
  wwg2= 0d0
  eee2= -1d99
  ikini= 1
  ierr= 1
  nne = nband*nqibz*nspin
  do ik = 1, nne
     if(eee2 +1d-4< ektx(ieaord(ik)) .OR. ik==nne ) then        ! degeneracy check
        if (mpi__root) then
           if(ik/=1) write(if8302,"(2i6,2d23.15)") ikini,ik-1,eee2,wwg2
        endif
        wwg2 = wgtx(ieaord(ik))
        eee2 = ektx(ieaord(ik))
        ikini =ik
     else
        wwg2= wwg2 + wgtx(ieaord(ik))
     endif
     wwg = wwg + wgtx(ieaord(ik))
     if(wwg<valn+2d0) write(6,*) ik,ieaord(ik),ektx(ieaord(ik)),wwg
     if (mpi__root) then
        write(if8301,"(2i6,3d23.15)") &
             ik,ieaord(ik),ektx(ieaord(ik)),wwg,wgtx(ieaord(ik))
     endif
     if( wwg>valn-1d-8 .AND. ierr==1 ) then
        write(6,*)
        efini = .5*(ektx(ieaord(ik+1))+ ektx(ieaord(ik)))
        if(autoew) then
           if(ik<3) call rx( ' efsimplef2: ik<3')
           esmr  = ektx(ieaord(ik)) - ektx(ieaord(ik-1))
        endif
        ierr=0
     endif
     if( wwg > valn+1d0) ikx=ik
  enddo
  if(ierr==1) call rx( ' efsimplef2: ierr=1 given nval is too large')
  nbnqnsp = nband*nqibz*nspin
  if(GaussSmear) then
     valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp)) ,efini,esmr,nbnqnsp)
  else
     valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp))      ,efini,esmr,nbnqnsp)
  endif
  write(6,*) 'valx at efini=',efini,valx
  if(abs(valx-valn)<1d-8) then
     ef=efini
     goto 8891
  endif
  ew1= ektx(ieaord(1))-0.01d0
  ew2= ektx(ieaord(ikx))
  nbnqnsp =nband*nqibz*nspin
  do ix = 1,100
     ein = 0.5d0*(ew1+ew2)
     if(GaussSmear) then
        valx= enumef_gauss(wgtx(ieaord(1:nbnqnsp)) ,ektx(ieaord(1:nbnqnsp)),ein,esmr,nbnqnsp)
     else
        valx= enumef( wgtx(ieaord(1:nbnqnsp)),ektx(ieaord(1:nbnqnsp))     ,ein,esmr,nbnqnsp)
     endif
     if(valx>valn) ew2=ein
     if(valx<valn) ew1=ein
     if(abs(ew1-ew2)<1d-15) exit
  enddo
  ef = 0.5d0*(ew1+ew2)
8891 continue
  if (mpi__root) then
     write(if8301,*) " ef=",ef
     close(if8301)
     write(if8302,*) " ef=",ef
     close(if8302)
  endif
  if(GaussSmear) then
     write(6,*)' efsimplef2ax(gauss):end'
  else
     write(6,*)' efsimplef2ax:end'
  endif
end subroutine efsimplef2ax
real(8) function enumef_gauss( wgtx,ektx,ein,esmr,nbnqnsp)
  implicit none
  integer:: nbnqnsp,ik
  real(8) :: ektx(nbnqnsp),wgtx(nbnqnsp),wwg, derfcx,ein,esmr
  wwg = 0d0
  do ik = 1, nbnqnsp
     wwg= wwg + wgtx(ik) *0.5d0* erfc( -(ein-ektx(ik))/sqrt(2d0)/esmr )
  enddo
  enumef_gauss = wwg
END function enumef_gauss
real(8) function enumef( wgtx,ektx,ein,esmr,nbnqnsp)
  implicit real*8(a-h,o-z)
  integer:: nbnqnsp,ik
  real(8) :: ektx(nbnqnsp),wgtx(nbnqnsp)
  wwg = 0d0
  do ik = 1, nbnqnsp
     if    (  ektx(ik) + 0.5d0*esmr < ein ) then
        wwg  = wwg + wgtx(ik)
     elseif(  ektx(ik) - 0.5d0*esmr < ein ) then
        wwg  = wwg + wgtx(ik)*(ein- (ektx(ik)-0.5d0*esmr))/esmr
     endif
  enddo
  enumef = wwg
END function enumef
