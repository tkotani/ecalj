      module m_globalvariables
      implicit none
      type t_globalvariables
      integer:: avw=-999999
      integer:: l_avw=0   ! flag for avw, =1 if its value is set
      integer:: lrel=-999999
      integer:: l_lrel=0
      integer:: lxcf=-999999
      integer:: l_lxcf=0
      integer:: mxorb=-999999
      integer:: l_mxorb=0
      integer:: nat=-999999
      integer:: l_nat=0
      integer:: nbas=-999999
      integer:: l_nbas=0
      integer:: nbasp=-999999
      integer:: l_nbasp=0
      integer:: nkaph=-999999
      integer:: l_nkaph=0
      integer:: nl=-999999
      integer:: l_nl=0
      integer:: nlmax=-999999
      integer:: l_nlmax=0
      integer:: nsp=-999999
      integer:: l_nsp=0
      integer:: nspc=-999999
      integer:: l_nspc=0
      integer:: nspec=-999999
      integer:: l_nspec=0
      integer:: stde=-999999
      integer:: l_stde=0
      integer:: stdl=-999999
      integer:: l_stdl=0
      integer:: stdo=-999999
      integer:: l_stdo=0
      end type t_globalvariables
      type(t_globalvariables) :: globalvariables
      end module m_globalvariables

