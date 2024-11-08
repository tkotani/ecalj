!> Clebsch-Gordan coefficients
module m_scg
  public scg,scg_sizechk
  ! --- CG coefficienets. <LM3|lm1 lm2>
  ! inxcg = lm1(lm1-1)/2 + lm2 (lm1>lm2)
  ! Injcg = indxcg(inxcg) to indxcg(inxcg)-1
  ! cg(inxcg)  : = <lm3|lm1 lm2>
  ! jcg(lnjcg) : = lm3
contains
  subroutine scg(lmax,c,cindx,js)
    !   call scg(lmxcg,cg,indxcg,jcg)
    !- Computes Clebsch-Gordan coefficients
    ! ----------------------------------------------------------------
    !i Inputs
    !i   lmax
    !o Outputs
    !o   c,cindx,js
    !r Remarks
    !r   (FORMERLY S104 IN ASW)
    ! ----------------------------------------------------------------
    implicit none
    integer :: lmax
    integer :: cindx(*),js(*)
    double precision :: c(*)
    integer :: i,i1,i2,i3,i31,i32,ic,j1,j1s,j2,j2s,k2,l1,l2,l3,lmindx, &
         m1=99999,m2=99999,m3=99999,mb,n1,n2,n3,nl,nm3=99999,s1,s2,s3,t1,t2,t3
    double precision :: q1,sr2,t,srpi,fs !,f100,f102
    double precision :: fac(161)
    data srpi /1.772453851d0/
    fs(i) = 1 + 4*(i/2) - 2*i
    mb = 999999
    nl = lmax+1
    sr2 = dsqrt(2d0)
    fac(1) = 1d0
    do  11  i = 1, 160
      fac(i+1) = i*fac(i)
11  enddo
    ic = 0
    lmindx = 0
    do  104  i1 = 1, nl
      l1 = i1-1
      j1s = 2*l1+1
      do  103  j1 = 1, j1s
        m1 = j1-i1
        n1 = iabs(m1)
        s1 = 0
        if (m1 < 0) s1 = 1
        t1 = 0
        if (m1 == 0) t1 = 1
        do  102  i2 = 1, i1
          l2 = i2-1
          i31 = l1 - l2 + 1
          i32 = l1 + l2 + 1
          j2s = 2*l2 + 1
          k2 = j1s*j2s
          if (i2 == i1) j2s = j1
          do  101  j2 = 1, j2s
            lmindx = lmindx + 1
            cindx(lmindx) = ic+1
            m2 = j2-i2
            n2 = iabs(m2)
            s2 = 0
            if (m2 < 0) s2 = 1
            t2 = 0
            if (m2 == 0) t2 = 1
            !!              
            !  if (m1*m2) 2,3,4
            if(m1*m2<0) then
               m3 = -n1 - n2
               mb = -iabs(n1-n2)
               if(mb == 0) then
                  nm3 = 1
               else
                  nm3 = 2
               endif
            elseif(m1*m2==0) then
               m3 = m1+m2
               nm3 = 1
            else 
               m3 = n1+n2
               mb = iabs(n1-n2)
               nm3 = 2
            endif
5           continue              
            n3 = iabs(m3)
            s3 = 0
            if (m3 < 0) s3 = 1
            t3 = 0
            if (m3 == 0) t3 = 1
            q1 = dsqrt(dble(k2))*fs(n3+(s1+s2+s3)/2)/(2*sr2**(1+t1+t2+t3))
            do  6  i3 = i31, i32,2
              l3 = i3-1
              if (n3 > l3) goto 6
              t = 0d0
              if (n1+n2 == -n3) t = t + f102(fac,l1,l2,l3)
              if (n1+n2 == n3)  t = t + f100(fac,l1,l2,l3,n1,n2,n3)*fs(n3+s3)
              if (n1-n2 == -n3) t = t + f100(fac,l1,l2,l3,n1,-n2,-n3)*fs(n2+s2)
              if (n1-n2 == n3)  t = t + f100(fac,l1,l2,l3,-n1,n2,-n3)*fs(n1+s1)
              ic = ic+1
              c(ic) = q1*t*f102(fac,l1,l2,l3)/(srpi*dsqrt(dble(2*l3+1)))
              js(ic) = l3*(l3+1) + m3 + 1
6           enddo
            nm3 = nm3-1
            m3 = mb
            if (nm3 > 0) goto 5
101       enddo
102     enddo
103   enddo
104 enddo
    cindx(lmindx+1) = ic+1
    return
  end subroutine scg

  subroutine scg_sizechk(lmax,lnjcg,lnxcg) !(lmax,c,cindx,js)
    !  computes clebsch-gordan coefficients (formerly s104 in asw)
    !  but here all is doubleprecision
    implicit real*8 (a-h,p-z)
    implicit integer(o)
    integer :: s1,s2,s3,t1,t2,t3,lmax,lnjcg,lnxcg,i1,i2,i3,i31,i32,j1,j1s,j2,k1,l1,l2,j2s,l3 &
         ,lmindx,m1=99999,m2=99999,m3=99999,mb,ic,k2,n2,n3,nl,nm3=99999,n1
    lnjcg=0
    lnxcg=0
    mb = 999999
    srpi = dsqrt(4*datan(1d0))
    nl=lmax+1
    sr2=dsqrt(2.d0)
    ic=0
    lmindx=0
    do 1111 i1=1,nl
      l1=i1-1
      j1s=2*l1+1
      do 111 j1=1,j1s
        m1=j1-i1
        n1=iabs(m1)
        s1=0
        if(m1 < 0) s1=1
        t1=0
        if(m1 == 0) t1=1
        do 11 i2=1,i1
          l2=i2-1
          i31=l1-l2+1
          i32=l1+l2+1
          j2s=2*l2+1
          k2=j1s*j2s
          if(i2 == i1) j2s=j1
          do 1 j2=1,j2s
            lmindx=lmindx+1
            lnxcg=max(lnxcg,lmindx)
            m2=j2-i2
            n2=iabs(m2)
            s2=0
            if(m2 < 0) s2=1
            t2=0
            if(m2 == 0) t2=1
            if(m1*m2 <0) then !  if(m1*m2) 2,3,4
               m3=-n1-n2
               mb=-iabs(n1-n2)
               if(mb == 0) then
                  nm3=1
               else
                  nm3=2
               endif
            elseif(m1*m2==0) then
               m3=m1+m2
               nm3=1
            else
               m3=n1+n2
               mb=iabs(n1-n2)
               nm3=2
            endif
5           continue
            n3=iabs(m3)
            s3=0
            if(m3 < 0) s3=1
            t3=0
            if(m3 == 0) t3=1
            do 6 i3=i31,i32,2
              l3=i3-1
              if(n3 > l3) goto 6
              t=0.d0
              ic=ic+1
              lnjcg=max(ic,lnjcg)
6           enddo
            nm3=nm3-1
            m3=mb
            if(nm3 > 0) goto 5
1         enddo
11      enddo
111   enddo
1111 enddo
    lnxcg=max(lnxcg,lmindx+1)
    return
  end subroutine scg_sizechk

  double precision function f100(fac,j1,j2,j3,m1,m2,m3)
    !o Outputs
    !o   F100
    integer :: j1,j2,j3,m1,m2,m3
    double precision :: fac(50)
    integer :: m,n,n1,n2
    double precision :: t,t1
    if (m3 /= m1+m2) goto 2
    t = (2*j3+1)*fac(j1+j2-j3+1)*fac(j3+j1-j2+1)*fac(j3+j2-j1+1)/ &
         fac(j1+j2+j3+2)
    t = dsqrt(t*fac(j1+m1+1)*fac(j1-m1+1)*fac(j2+m2+1)*fac(j2-m2+1)* &
         fac(j3+m3+1)*fac(j3-m3+1))
    n1 = max0(j2-j3-m1,j1-j3+m2,0) + 1
    n2 = min0(j1+j2-j3,j1-m1,j2+m2) + 1
    if (n1 > n2) goto 2
    t1 = 0d0
    do  m = n1, n2
      n = m-1
      t1 = t1 + &
           dble(1+4*(n/2)-2*n)/(fac(m)*fac(j1+j2-j3-n+1)*fac(j1-m1-n+1)* &
           fac(j2+m2-n+1)*fac(j3-j2+m1+n+1)*fac(j3-j1-m2+n+1))
    enddo
    f100 = t*t1
    return
2   f100 = 0d0
    return
  END function f100

  double precision function f102(fac,l1,l2,l3)
    integer :: l1,l2,l3
    double precision :: fac(50)
    integer :: lt,p,x
    lt = l1 + l2 + l3
    p = lt/2
    if (2*p /= lt) goto 1
    f102 = dsqrt(dble(2*l3+1)/dble(lt+1))
    f102 = f102*fac(p+1)/dsqrt(fac(2*p+1))
    x = p-l1
    f102 = f102*dsqrt(fac(2*x+1))/fac(x+1)
    x = p-l2
    f102 = f102*dsqrt(fac(2*x+1))/fac(x+1)
    x = p-l3
    f102 = f102*dsqrt(fac(2*x+1))/fac(x+1)
    if (x > 2*(x/2)) f102 = -f102
    return
1   continue
    f102 = 0d0
  END function f102

endmodule m_scg
