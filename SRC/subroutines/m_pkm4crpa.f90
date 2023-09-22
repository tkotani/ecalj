!>read pkmU for cRPA calclaitons ===
module m_pkm4crpa 
  implicit none !!! subroutine wmaxloc in maxloc2.F write pkmU file.
  public:: Readpkm4crpa
  private
  real(8),allocatable:: pkmud(:,:,:),qvecud(:,:,:)
  integer:: nqbz,nsp,iko_ix,iko_fx,nwf
  logical:: init=.true.
contains
  function Readpkm4crpa(ib,q,isp) result(wpw)
    use m_iqindx_qtt,only: iqindx2_
    intent(in)::            ib,q,isp
    !      intent(out)::                    wpw
    integer:: ib,isp,iq,ibin,iqin,ifi,ispin
    real(8):: wpw,q(3),wpwin,qvec(3),qin(3)
    character*8::fname
    logical:: outofrange,debug=.false.
    if(debug)print *,' iiiiiiqqq readpkm4crpa',init
    fname= 'pkm4crpa'
    if(init) then
       open(newunit=ifi,file=fname,form='formatted',status='unknown')
       read(ifi,*)
       read(ifi,*)
       read(ifi,*) nqbz,nwf,nsp,iko_ix,iko_fx
       read(ifi,*)
       allocate(qvecud(3,nqbz,nsp),pkmud(iko_ix:iko_fx,nqbz,nsp))
       qvecud=0d0; pkmud=0d0
       do
          read(ifi,*,end=991) wpwin, ibin,iqin,ispin,qin
          pkmud(ibin,iqin,ispin)=wpwin
          qvecud(:,iqin,ispin)=qin
       enddo
991    continue
       init=.false.
       close(ifi)
    endif
    if(debug) print *,' qqqqqqqqqqqqq',nqbz,nwf,nsp,iko_ix,iko_fx
    if(debug) print *,' qqqqqqqq init',isp,init,ib
    outofrange = ib<iko_ix .or. iko_fx<ib
    if(outofrange) then
       wpw=0d0
       if(debug) print *,' outofrange',outofrange
       return
    endif
    !! find iq
    if(debug) print *,'qqqqqqqqqqq  q=',q
    call iqindx2_(q, iq, qvec) !qvec is used q. q-qu is a G vector.
    if(debug) print *,'qqqqqqqqqqq qvec iq=',qvec,iq
    wpw = pkmud(ib,iq,isp)
    if(debug) print *,'qqqqqqqq xxxxxxx qqqqqqqq',qvecud(:,iq,isp)
    if(sum(abs(qvecud(:,iq,isp)-qvec))>1d-5 ) then
       call rx('readpkm4crpa:readq in pkmu is not consistent with given q')
    endif
  end function Readpkm4crpa
end module m_pkm4crpa
