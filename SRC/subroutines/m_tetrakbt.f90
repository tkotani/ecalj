!>finite-temperature tetrahedron method: the temperature kbt (tetrakbt_init) and integtetn.
!! `[gw] tetrakbt=true` is served by lindtet6_kbt (in tetwt5.f90), the exact energy convolution of the T=0 lindtet6.
!! (2026-10-02 05:26) The old midpoint-factorization routines (tetrakbt, eaf_triangle, factri, ... ; they gave a wrong
!! chi0 and were no longer called) are removed; the diagnosis and how to get them back: ecaljdoc/MD/past_log.md §13.
module m_tetrakbt
  use m_keyvalue,only: getkeyvalue
  implicit none
  public:: tetrakbt_init, kbt,integtetn
  private
  real(8):: tt,kbt !tt: temperature [K]
  logical:: init = .true.
  real(8):: kb=8.6171d-5 ![eV/K]
contains
  !-----------------------------------------------------
  subroutine tetrakbt_init()
    use m_GWinput, only: gwinput_init, gwinput_loaded, tg_t_tetrakbt => t_tetrakbt
    real(8):: temperature, rydberg
    call gwinput_init()
    if (gwinput_loaded) then
       temperature = tg_t_tetrakbt
    else
       call rx('m_GWinput: legacy GWinput reader is disabled; ctrlg.<sname>.toml is required.')
!       call getkeyvalue("GWinput","t_tetrakbt",temperature,default=3d+2)
    endif
    tt = temperature+1d-12 !avoid 0
    kbt=kb*tt/rydberg()
    if (init) then
       write(6,"(' tetrakbt_init: T[K], kbt[Ry], kbt[eV]',I5,2E13.5)") int(tt),kbt,kbt*rydberg()
       init = .false.
    endif
  end subroutine tetrakbt_init
  !-----------------------------------------------------
  subroutine integtetn(e, ee, integb)  ! Calculate primitive integral of integb = 1/pi Imag[\int^ee dE' 1/(E' -e(k))] = \int^ee dE' S[E']
    !! \remark
    !!  S[E] : is area of the cross-section between the omega-constant plane and the tetrahedron. [here we assumee e1<e2<e3<e4].
    !!  Normalization is not considered! Rath&Freeman Integration of imaginary part on Eq.(17)
    implicit none
    real(8)::  e(1:4), ee, integb,a,b !,e1,e2,e3,e4 V1,V2,V3,V4 , ,D1,D2,D3,D4
    associate( e1=>e(1)-3d-8, e2=>e(2)-2d-8, e3=>e(3)-1d-8, e4=>e(4))
      associate(V1=>ee-e1, V2=>ee-e2, V3=>ee-e3,V4=>ee-e4)
        if(ee<e1) then
           integb=0d0
        elseif( e1<=ee .AND. ee<e2 ) then
           integb = V1**3/((e4-e1)*(e3-e1)*(e2-e1))
        elseif( e2<=ee .AND. ee<e3 ) then
           a  =  V1/ ((e4-e1)*(e3-e1)*(e2-e1))**(1d0/3d0)
           b  =  V2/(-(e4-e2)*(e3-e2)*(e1-e2))**(1d0/3d0)
           integb = (a-b) * (a**2+a*b+b**2)
        elseif( e3<=ee .AND. ee<e4 ) then
           integb = 1d0 - V4**3/((e1-e4)*(e2-e4)*(e3-e4))
        elseif( ee==e4 ) then
           integb = 1d0
        else
           call rx( ' integtetn: ee>e4')
        endif
      endassociate
    endassociate
  end subroutine integtetn
end module m_tetrakbt


