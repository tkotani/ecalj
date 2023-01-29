module m_vxtrap !- Extrapolate potential, radial wave function outside MT boundary
  public rwftai
  private
contains
  subroutine rwftai(rmt,a,nrmt,nrbig,ribig,phi,dphi,tphi,l, ehl,rsml,g)
    use m_hansr,only :hansr
    !- Extend radial wave function outside MT boundary
    ! ----------------------------------------------------------------------
    !i Inputs 
    !i         :Compute radial wave function on extended mesh using
    !i         : rsml,ehl for tail, scale gz so that value
    !i         : matches envelope h(rsm,eh)
    !i   rmt   :augmentation radius, in a.u., by species
    !i   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
    !i   nrmt  :number of radial mesh points between origin and rmt
    !i   l     :quantum number for this wave function
    !i   ehl   :energy of smoothed Hankel tail for local orbital
    !i   rsml  :smoothing radius for smoothed Hankel tail of local orbital
    !i  nrbig  :number of radial mesh points on extended mesh
    !i  ribig  :points for extended radial mesh
    !o Outputs
    !o   g      :radial wave function extended to points nrmt..nrbig
    !l Local variables
    !l   rbig   :rmax of extended mesh
    !l   nxtn   :number of points between rmt and rbig
    implicit none
    integer :: mode,nrmt,nrbig,l
    double precision :: a,rmt,phi,dphi,tphi,ribig(*),g(nrbig,2),ehl,rsml
    integer :: nrx,idn,IPRTDB,nxtn,ir,info,mode0
    parameter (nrx=1501,IPRTDB=10)
    double precision :: rbig,rwgt(nrx),xi(nrx,0:l),r2(nrx)
    double precision :: fac1,fac2,rsm,ekin,eh,alfa,beta
    rbig = ribig(nrbig)
    nxtn = nrbig-nrmt+1
    call radwgt(rbig,a,nrbig,rwgt)
    rwgt(nrmt) = rwgt(nrmt)/2
    r2(1:nxtn) = [(ribig(ir-1+nrmt)**2,ir=1,nxtn)]
    call hansr(rsml,0,l,1,[l],[ehl],[r2],nrx,nxtn,[idn],11,xi)
    fac1 = g(nrmt,1)/rmt/xi(1,l) ! Hankel function scaled to match g at nrmt
    fac2 = g(nrmt,2)/rmt/xi(1,l)
    g(1:nrmt,:)=g(1:nrmt,:)/fac1 
    fac2 = fac2 / fac1
    fac1 = 1
    g(nrmt:nrmt+nxtn-1,1) = [(xi(ir,l)*ribig(ir-1+nrmt) * fac1,ir=1,nxtn)]
    g(nrmt:nrmt+nxtn-1,2) = [(xi(ir,l)*ribig(ir-1+nrmt) * fac2,ir=1,nxtn)]
  end subroutine rwftai
end module m_vxtrap
