      module m_work
!! These are used to pass data set.
!! It is not clear where they are generated and where cited.
!!     bz
      use m_struc_def,only: s_rv1
c      integer , allocatable ::  iv_a_oidtet (:) !moved to m_mkqp
c      integer , allocatable ::  iv_a_oipq (:)
c      real(8) , allocatable ::  rv_p_oqp (:) 
c      integer , allocatable ::  iv_a_ostar (:)
c      real(8) , allocatable ::  rv_a_owtkp (:)
!! s_array      
      integer , allocatable ::  iv_a_oips (:)
!! ham
      real(8) , allocatable ::  rv_a_ohrs (:)
      integer, allocatable :: iv_a_oiaxs(:)
      integer , allocatable ::  iv_a_oindxo (:)
      integer, allocatable :: iv_a_ontabs(:)
      integer , allocatable ::  iv_a_ooffH (:)
      real(8) , allocatable ::  rv_p_oqsig (:) 
!! lat
      complex(8) , allocatable ::  zv_a_obgv (:)
      real(8) , allocatable ::  rv_a_ocg (:)
      real(8) , allocatable ::  rv_a_ocy (:)
      real(8) , allocatable ::  rv_a_odlv (:)
      real(8) , allocatable ::  rv_a_ogv (:)
      integer , allocatable ::  iv_a_oidxcg (:)
      integer , allocatable ::  iv_a_oips0 (:)
      integer , allocatable ::  iv_a_ojcg (:)
      integer , allocatable ::  iv_a_okv (:)
      real(8), allocatable :: rv_a_opos(:)
      real(8) , allocatable ::  rv_a_oqlv (:)
!! spot
      real(8) , allocatable ::  rv_a_obxc (:)
      real(8) , allocatable ::  rv_a_ogrrme (:)
      real(8) , allocatable ::  rv_a_omad (:)
      type(s_rv1), allocatable :: sv_a_oorhat(:)
      real(8) , allocatable ::  rv_a_opmpol (:)
      real(8) , allocatable ::  rv_a_opnu (:)
      real(8) , allocatable ::  rv_a_opp (:)
      real(8) , allocatable ::  rv_a_opprel (:)
      real(8) , allocatable ::  rv_a_oqnu (:)
      real(8) , allocatable ::  rv_a_oqpp (:)
      real(8) , allocatable ::  rv_a_oqt (:)
      real(8) , allocatable ::  rv_a_orhrmx (:)
      complex(8) , allocatable ::  zv_a_osmpot (:)
      complex(8) , allocatable ::  zv_a_osmrho (:)
      real(8) , allocatable ::  rv_a_osop (:)
      real(8) , allocatable ::  rv_a_ovdif (:)
      real(8) , allocatable ::  rv_a_oves (:)
      real(8) , allocatable ::  rv_a_ovintr (:)
      real(8) , allocatable ::  rv_a_ovrmax (:)
      
      end module
