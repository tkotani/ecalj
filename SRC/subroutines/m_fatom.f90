!>Free atom density determined by lmfa. lmf stores data into rst.*. But unchanged.
! mpibc1_s_spec was here; every caller (iors, rdovfa) now reads __atm/
! rst on every rank itself so the broadcast is no longer needed.
module m_fatom
  integer,parameter::  n0=10
  type s_spec
     ! I think lmfa detemines all the following data and write to atm.* files
     ! The data is used for lmf-MPIK and copied into rst file (unchanged).
     real(8):: qc    !core charge
     real(8):: rsmfa !rsm to fit free atom density
     real(8):: ctail !coefficients to fit of free-atom core tail by unsm. Hankel
     real(8):: etail !energy to fit of free-atom core tail
     real(8):: stc   !core kinetic energy
     integer:: nxi    ! Number of energies in fit of free-atom density tails
     real(8):: exi(n0)    ! Hankel energies for fit to c.d.; fit to free-atom density tails.
     real(8):: chfa(n0,2) ! coefficients to fit of free-atom density tails
     real(8),allocatable :: rv_a_orhoc(:) !pointer to core density
  end type s_spec
  type(s_spec),allocatable:: sspec(:) !just allocated for iors and rdovfa. Not touched.
end module m_fatom
