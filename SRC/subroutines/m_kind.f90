!>To contol the computational size for dielectric function and G\times W
module m_kind 
#ifdef __MP
  integer,parameter:: kindrcxq=4 !or 4
  integer,parameter:: kindgw=4  !or 4
  integer,parameter:: kindzmel=4  !or 4
#else
  integer,parameter:: kindrcxq=8 !or 4
  integer,parameter:: kindgw=8  !or 4
  integer,parameter:: kindzmel=8  !or 4
#endif
endmodule m_kind

