!> Single space-group operations (g,ag): products, equality, lattice vectors, the symbol of an operation.
!> Split from m_mksym_util (2026-10-02); used by m_symderive and m_symfind.
module m_symop_util
  use m_lgunit,only:stdo
  use m_ftox
  public spgcop,spgprd,spgeql,grpeql,latvec,asymop
  private
  real(8),parameter:: toll=1d-4,tiny=1d-4,epsr=1d-12
contains
  subroutine asymop(grpin,ag,asep,sg)  ! Generate the symbolic representation sg of a group operation 
    !i  grpin,ag :  space group rotation + translation matrix
    !i  asep: 
    !o  sg  :  symbolic representation of group op
    implicit none
    real(8) :: grp(3,3),ag(3),vecg(3),grpin(3,3),costbn,detop,ddet33,dnrm2,sinpb3,vfac,wk(9)
    character(*):: sg,asep
    integer :: nrot,ip,isw,i1,i2,fmtv,llen,i,idamax,j,in
    logical :: li
    real(8),parameter:: twopi = 8*datan(1d0)
    character(8):: xn
    grp=grpin
    call dinv33(grp,0,wk,detop)
    if(dabs(dabs(detop)-1d0)>tiny) call rx('Exit -1 ASYMOP: determinant of group op must be +/- 1, but is '//trim(ftof(detop)))
    detop = dsign(1d0,detop) !sign of determinant
    li = detop<0d0 !   ... li is T if to multiply by inversion
    grp= detop*grp !Multiply operation grp with detop to guarantee pure rotation 
    costbn = 0.5d0*(-1 + grp(1,1) + grp(2,2) + grp(3,3))
    if (dabs(costbn-1d0) < tiny) then
       nrot = 1
       vecg=0d0 
    else
       nrot = idnint(twopi/dacos(dmax1(-1d0,costbn)))
       if (nrot == 2) then
          vecg = 0.5d0*[((grp(i,i)+1.0d0),i=1,3)]
          j = idamax(3,vecg,1)
          if(vecg(j) < 0d0)call rx('ASYMOP: bad operation j='//trim(xn(j))//'. Diagonal element is '//ftof(grp(j,j)))
          vecg(j) = dsqrt(vecg(j))
          vfac = 0.5d0/vecg(j)
          do i = 1, 3
             if (i /= j) vecg(i) = vfac*grp(i,j)
          enddo
       else
          vecg=[grp(3,2)-grp(2,3), grp(1,3)-grp(3,1), grp(2,1)-grp(1,2)]
       endif
       sinpb3 = dsqrt(.75d0) 
       if (dabs((sinpb3-dabs(vecg(1)))*(sinpb3-dabs(vecg(2)))*(sinpb3-dabs(vecg(3)))) > tiny) then
          do  j = 3, 1,-1 !Renormalize at least one component to 1 to allow for abbreviations as 'D', 'X', 'Y' or 'Z'
             vfac = dabs(vecg(j))
             if(vfac > tiny) vecg=1d0/vfac*vecg
          enddo
       endif
    endif
    sg=''
!    write(stdo,ftox)'nrotnnnnn',nrot,li,'vecg',vecg,'ag=',ag
    if(nrot == 1) then ! Rotational part
       sg = merge('i','e',li) 
       ip=len(trim(sg))+1
    else
       if(li.and.nrot==2) then
          sg='m'
       else   
          sg=merge('i*','  ',li)//'r'//char(48+nrot)
       endif
       ip=len(trim(sg))+1
       call rxx(.not. parsvc2(.true.,sg,ip,vecg),'bug in asymop 2')!rotation axis
    endif
    sg=adjustl(sg)
    if(sum(abs(ag))>tiny) then !Translational part added
!       print *,'sg=',sg,'asep=',trim(asep)
       if(asep(1:1)/=' ') sg=trim(sg)//trim(asep)
       ip=len(trim(sg))+1
       call rxx(.not. parsvc2(.false.,sg,ip,ag),'bug in asymop 1')
     endif  
  end subroutine asymop
  logical function parsvc2(modex,t,ip,v)  
    implicit none
    logical:: modex
    integer :: ip
    real(8) :: v(3)
    character(*) :: t
    real(8) :: x,y,z,d
    character sout*50, add*1,soutx*50
    integer :: itrm,ix(3),ich,iopt,m,i,iz,id,mx !,awrite !,a2vec
    character(9),parameter:: rchr='(XxYyZzDd'
    parsvc2 = .true.
    t(ip:ip)=' '
    if(modex) then
       if( all(abs(v(:)-[1d0,1d0,1d0])<tiny) ) t(ip:ip)='d'
       if( all(abs(v(:)-[1d0,0d0,0d0])<tiny) ) t(ip:ip)='x'
       if( all(abs(v(:)-[0d0,1d0,0d0])<tiny) ) t(ip:ip)='y'
       if( all(abs(v(:)-[0d0,0d0,1d0])<tiny) ) t(ip:ip)='z'
    endif
    if(t(ip:ip)==' ') then
       t(ip:ip)='('
       do i = 1, 3
          write(sout,ftox) ftof(v(i))
          if(abs(nint(v(i))-v(i))<1d-6) write(sout,ftox) nint(v(i))
          sout = trim(adjustl(sout))//merge(')',',',i==3)
          m= len_trim(sout)
          t(ip+1:ip+m)=trim(sout) 
          ip = ip+m
       enddo
    endif
    ip = ip+1
  end function parsvc2
  subroutine spgcop(g,ag,h,ah)
    real(8):: h(9),g(9),ag(3),ah(3)
    h = merge(0d0, g, dabs(g) <1d-8)
    ah= merge(0d0,ag, dabs(ag)<1d-8)
  end subroutine spgcop
  subroutine spgprd(g1,a1,g2,a2,g,a)
    implicit none
    real(8) :: g1(3,3),g2(3,3),g(3,3),sum,a1(3),a2(3),a(3),h(3,3),ah(3)
    integer :: i,j,k
    h=matmul(g1,g2)
    g=h !tk does not know why g=matmul(g1,g2) fails for gfortran gcc9.4.0 2023march
    ah=a1+matmul(g1,a2)
    a=ah
  end subroutine spgprd
  logical function spgeql(g1,a1,g2,a2,qb) ! Determines whether space group op g1 is equal to g2
    implicit none !i      g1,a1 :first space group,  g2,a2 :second space group, qb:reciprocal lattice vectors
    integer :: m,iq,iac
    real(8) :: g1(9),g2(9),a1(3),a2(3),qb(3,3),adiff(3)
    adiff = matmul(a1-a2,qb)
    spgeql= all([dabs(g1-g2),abs(adiff-nint(adiff))]<toll)
  end function spgeql
  logical function grpeql(g1,g2)    !- Checks if G1 is equal to G2
    implicit none
    real(8):: g1(9),g2(9)
    grpeql = all(dabs(g1-g2)<toll)
  end function grpeql
  logical function latvec(n,tol,qlat,vec) ! Checks whether a set of vec(1:3,n) are lattice vectors
    implicit none
    integer:: n
    real(8):: qlat(3,3),vec(3,n),tol, vdiff(n,3)
    vdiff  = matmul(transpose(vec(:,:)),qlat(:,:))
    latvec = all(reshape(abs(vdiff-nint(vdiff)),[n*3]) < tol)
  end function latvec
end module m_symop_util
