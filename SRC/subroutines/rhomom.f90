!>Multipole moments Q_aL = Q_aL^Zc+ Q_aL^v (See.Eq.(28) JPSJ034702)
module m_rhomom
  public rhomom
  private
  contains
  !i   orhoat: vector of offsets containing site density
  !o   qmom  : Multipole moments
  !o         : = integral r^l Y_L (rho1-rho2) + smooth core 
  !o   vsum  : sum over all sites ib of difference vs1_ib - vs2_ib
  !o         : vs1_ib = integral in sphere of estat potential[true rho]
  !o         : vs2_ib = integral in sphere of estat potential[sm rho]
subroutine rhomom(sv_p_orhoat, qmom,vsum)
  use m_struc_def
  use m_lmfinit,only: z_i=>z,nr_i=>nr,lmxa_i=>lmxa,rmt_i=>rmt,lmxb_i=>lmxb,lmxl_i=>lmxl,spec_a
  use m_lmfinit,only: kmxt_i=>kmxt,lfoca_i=>lfoca,rfoca_i=>rfoca,rg_i=>rg, nsp,nbas,jnlml,ispec,nlmxlx
  use m_lgunit,only:stdo
!  use m_hansr,only:corprm
  implicit none
  intent(in) ::    sv_p_orhoat
  intent(out)::                 qmom,vsum
  type(s_rv1) :: sv_p_orhoat(3,nbas)
  integer:: ipr,iprint,j1,ib,is,igetss,lmxl,nr,nlml,ilm,j,lfoc
  real(8) ,allocatable :: rofi(:),rwgt(:), h_rv(:), v_rv(:)
  real(8):: qmom(nlmxlx,nbas) , vsum,vs1(nbas),vs2(nbas), z,qc,a,rmt,qcorg,qcorh,qsc,cofg,cofh,rg,ceh,rfoc
  real(8),parameter:: fpi  = 16d0*datan(1d0), y0 = 1d0/dsqrt(fpi),pi = 4d0*datan(1d0), srfpi = dsqrt(4d0*pi)
  ipr = iprint()
  if(ipr>=45) write(stdo,"(/' rhomom:   ib   ilm      qmom',8x,'Qval',7x, 'Qc',8x,'Z')")
  do  ib = 1, nbas
     is = ispec(ib) 
     lmxl=lmxl_i(is)
     if (lmxl == -1) cycle
     z  = z_i(is)
     a  = spec_a(is)
     nr = nr_i(is)
     rmt= rmt_i(is)
     rg = rg_i(is)
     call corprm(is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     nlml = (lmxl+1)**2
     call pvrhom(rmt,a,nlml,nr,nsp, sv_p_orhoat(1,ib)%v, sv_p_orhoat(2,ib)%v,cofg,z, qmom(:,ib),vs1(ib),vs2(ib))
     if(ipr>=45) then
        qc = qcorg+qcorh
        write(stdo,220) ib,1,qmom(1,ib),qmom(1,ib)/y0,qc,z
        do ilm = 2, nlml; if (dabs(qmom(ilm,ib)) > 1d-6) write(stdo,220) ib,ilm,qmom(ilm,ib); enddo
     endif
  enddo
  vsum=sum(vs1-vs2)
  220 format(i13,i6,f12.6,f12.6,2f9.2)
end subroutine rhomom
subroutine pvrhom(rmt,a,nlml,nr,nsp, rho1,rho2,cofg,z, qmomj,vsum1,vsum2)  !Multipole moments qmom = Q_L = Q_aL^Zc+ Q_aL^v (See.Eq.(28),(30) in JPSJ034702)
  use m_hansmr,only: hansmr,hansmronly
  use m_hansr,only:  hansr
  use m_ll,only:ll
  !i   nlml  : L-cutoff for charge (l+1)**2
  !i   rho1  :local true density*r**2, tabulated on a radial mesh
  !i   rho2  :local smoothed density*r**2, tabulated on a radial mesh
  !i   cofg  : Q_aL^c
  !i   z     :nuclear charge
  !o   qmomj :multipole moments = Q_aL^Zc + Q_aL^v (See.Eq.(28) JPSJ034702) = Q_aL^Zc + \integral r^l (rho1-rho2)
  implicit none
  integer :: nlml,nr,nsp,i,ilm,l,m,lmxl,isp
  real(8) :: cofg,z,rmt,rofi(nr),rwgt(nr),h(nr),qmomj(nlml), &
       rho1(nr,nlml,nsp),rho2(nr,nlml,nsp),a,vsum1,vsum2,vhrho,vsum,cg,af,q2,v(nr)
  real(8),parameter:: pi=4d0*datan(1d0), srfpi = dsqrt(4d0*pi),y0 = 1d0/srfpi,fpi = 4d0*pi
  call radmsh( rmt , a , nr , rofi )
  call radwgt( rmt , a , nr , rwgt )
  qmomj=0d0
  do isp=1,nsp
     do  ilm=1,nlml
        qmomj(ilm)= qmomj(ilm)+ sum(rwgt*rofi**ll(ilm)*(rho1(:,ilm,isp)-rho2(:,ilm,isp))) !, isp=1,nsp) ]),ilm=1,nlml)] !Q_aL^v Eq.(28)
     enddo
  enddo
! next line gave wrong results for nvfortran  
!  qmomj= [(sum([ (sum(rwgt*rofi**ll(ilm)*(rho1(:,ilm,isp)-rho2(:,ilm,isp))), isp=1,nsp) ]),ilm=1,nlml)] !Q_aL^v Eq.(28)
  qmomj(1) = qmomj(1) -y0*z + cofg ! Eq.(25). Add Q_al^Zc.  -y0*z = QaL[-Z_a\delta(\bfr)] ; cofg=QaL[ n^c_a(\bfr) - n^c_sH,a(\bfr)]
  call poiss0(z,  rofi, [(0d0,i=1,2*nr)],  nr,0d0, v, vhrho,vsum,1)
  vsum1 = sum(rwgt(2:nr)*rofi(2:nr)**2*(v(2:nr)-2d0*z/rofi(2:nr)))
  call poiss0(0d0,rofi,[(0d0,i=1,2*nr)],nr,0d0, v, vhrho,vsum,1)
  vsum2 = fpi*sum(rwgt*rofi**2*v)
end subroutine pvrhom
endmodule m_rhomom
