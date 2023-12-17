!>pointer arrays
module m_struc_def 
  public s_rv1, s_rv2, s_nv2, s_cv1, s_cv2,s_cv3,s_cv4, s_sblock,s_rv5,s_rv4,s_cv5
  private
  type s_i
     integer,allocatable:: i(:)
  end type s_i
  type s_rv1
     real(8),allocatable:: v(:)
  end type s_rv1
  type s_rv2
     real(8),allocatable:: v(:,:)
  end type s_rv2
  type s_rv4
     real(8),allocatable:: v(:,:,:,:)
  end type s_rv4
  type s_rv5
     real(8),allocatable:: v(:,:,:,:,:)
  end type s_rv5
  type s_cv5
     complex(8),allocatable:: cv(:,:,:,:,:)
  end type s_cv5
  type s_cv1
     complex(8),allocatable:: cv(:)
  end type s_cv1
  type s_sblock
     complex(8),allocatable:: sdiag(:,:,:,:,:) !(1,1) and (2,2) spin block as sdiag(:,:,:,:, isp)
     complex(8),allocatable:: soffd(:,:,:,:,:) !(1,2) and (2,1) spin block as soffd(:,:,:,:, isp)
  end type s_sblock
  type s_nv2
     integer,allocatable:: nv2(:,:)
  end type s_nv2
  type s_cv2
     complex(8),allocatable:: cv2(:,:)
  end type s_cv2
  type s_cv3
     complex(8),allocatable:: cv3(:,:,:)
  end type s_cv3
  type s_cv4
     complex(8),allocatable:: cv4(:,:,:,:)
  end type s_cv4
end module m_struc_def
