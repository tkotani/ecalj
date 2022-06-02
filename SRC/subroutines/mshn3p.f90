subroutine mshn3p ( nbas , ssite , sspec , lmet , lrout , lfrce &
     , qval , ef0 , def , sumqv , sumev , n1 , n2 , n3 , k1 , k2 , &
     k3 , smrho , sv_p_oqkkl , f , lrep )
  use m_struc_def
  use m_lmfinit,only: nkaph,nsp
  use m_lgunit,only:stdo,stdl


  !      use m_globalvariables
  !- Interpolate the density on the three-point energy mesh for metals
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nbas  :size of basis
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i     Stored:    *
  !i     Passed to: *
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: rsma lmxl lmxb lmxa kmxt
  !i     Stored:    *
  !i     Passed to: *
  !i   lmet  :0 semiconductor--do nothing
  !i         :>0 interpolate using 3-point sampling scheme; see Remarks
  !i   lrout :0 interpolate ef0,sumqv,sumev only
  !i         :1 also interpolate smrho,oqkkl,f
  !i   lfrce :>0 interpolate forces
  !i   qval  :total valence charge
  !i   def   :Fermi level window.
  !i   k1,k2,k3:dimensions smrho,smpot
  ! o Inputs/Outputs
  !i   ef0   :On input: trial  Fermi level; see Remarks
  !i         :On output: interpolated Fermi level
  !i   sumqv :on input, total charge at each of the trial Fermi levels
  !i         :second channel is magentic moment
  !i         :out output, interpolated charge and magnetic moment
  !i   sumev :on input, eigenvalue sum at each trial Fermi level
  !i         :out output, interpolated eigenvalue sum
  !o Outputs
  !o   smrho :smooth density on uniform mesh interpolated to true Fermi level
  !o   oqkkl :site density matrix interpolated to true Fermi level
  !o   f     :forces interpolated to true Fermi level
  !o   lrep  :set to nonzero if Fermi level outside acceptable bounds
  !r Remarks
  !r   Input quantities smrho,qkkl,f were accumulated input for three
  !r   separate fermi levels, ef0, ef0+def, ef0-def.  The true Fermi
  !r   level is calculated, and a quadratic interpolation is made
  !r   for these three quantities.
  !r
  !r   Interpolated data is put back into first slot (iq=1).
  !r   If Ef is crazy and the charges are unusable, set flag lrep and return.
  !u Updates
  !u   02 Jan 06 sumqv resolved by spin
  !u   21 Mar 01 Bug fix for case q1,q3 close to qval
  !u   23 Jan 01 Added lrout switch
  !u   17 Jun 00 spin polarized
  !u   22 May 00 Adapted from nfp ip_density.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer:: n1 , n2 , n3 , k1 , k2 , k3 , lmet , lfrce , lrout &
       , lrep , nbas,i_copy_size
  type(s_rv1) :: sv_p_oqkkl(3,1)

  double precision :: def,ef0,qval
  real(8):: sumev(2,3) , sumqv(3,2) , f(3,nbas,*)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)

  double complex smrho(k1,k2,k3,3,2)
  ! ... Local parameters
  integer :: ib,ipl,is,kmax,lmxa,lmxh,lmxl,m,nelt1,nelt2, &
       nelt3,nlma,nlmh,lgunit,ipr,iprint
  double precision :: a,b,c,def1,disc,eadd,efold,p(3),px,py,pz, &
       q1,q2,q3,a1,a2,a3,rsm,sen,sev,sqv,sav,w1,w2,w3,x,x1,x2,xlin
  equivalence (px,p(1)),(py,p(2)),(pz,p(3))
  ! ... Heap

  !      stdo = lgunit(1)
  !      stdl = lgunit(2)
  ipr = iprint()
  ipl = ipr
  lrep = 0
  ! angenglob      nsp  = nglob('nsp')
  !      nsp  = globalvariables%nsp
  ! angenglob      nkaph = nglob('nkaph')
  !      nkaph = globalvariables%nkaph

  ! ... Semiconductors: don't do anything
  if (lmet == 0) then
     sqv = sumqv(1,1) + sumqv(1,2)
     !       smm = sumqv(1,2) - sumqv(1,2)
     sev = sumev(1,1)
     if (ipr >= 20) write(stdo,800) sqv,sev
800  format(/' Semiconductor branch:  sumqv=',f12.6, &
          '   sumev=',f12.6)
     return
  endif

  ! --- Metals: make weights w1,w2,w3 for the three trial Ef's ---
  q1 = sumqv(1,1) + sumqv(1,2)
  q2 = sumqv(2,1) + sumqv(2,2)
  q3 = sumqv(3,1) + sumqv(3,2)
  a1 = sumqv(1,1) - sumqv(1,2)
  a2 = sumqv(2,1) - sumqv(2,2)
  a3 = sumqv(3,1) - sumqv(3,2)
  if (ipr >= 30) write(stdo,300) qval,ef0-def,ef0,ef0+def,q1,q2,q3
300 format(/' Interpolate to Fermi energy for qval=',f9.3 &
       /' trial Fermi energies: ',3f12.5 &
       /' band occupations:     ',3f12.5)
  if (ipr >= 30 .AND. nsp == 2) write(stdo,301) a1,a2,a3
301 format(' magnetic moment:      ',3f12.5)
  if (dabs(q3-q1) < 1d-8) then
     if (dabs(q2-qval) > 1d-8) then
        write(stdo,"('mshn3p: Input Fermi energy unusable, Q= ',f15.8)") q2
        sumev(1,1) = 0
        sumev(2,1) = 0
        lrep = 1
        return
     endif
     w1 = 0d0
     w2 = 1d0
     w3 = 0d0
     x = 0
     if (ipr >= 10) write(stdo,430) w1,w2,w3
430  format(' weights (converged):  ',3f12.5)
     goto 90
  endif
  if ((qval > q3 .OR. qval < q1) .AND. ipr > 0) &
       write(stdo,"('mshn3p: qval outside window; extrapolating')")
  a = (q1+q3-2d0*q2)/(2d0*def*def)
  b = (q3-q1)/(2d0*def)
  c = q2-qval
  disc = b*b - 4d0*a*c
  if (qval <= q2) xlin = -def*(qval-q2)/(q1-q2)
  if (qval > q2) xlin = def*(qval-q2)/(q3-q2)
  if (disc < 0d0) then
     x = xlin
     if (qval >= q2) then
        w1 = 0d0
        w2 = 1d0 - x/def
        w3 = x/def
     else
        w1 = -x/def
        w2 = 1d0 + x/def
        w3 = 0d0
     endif
     if (ipr >= 10) write(stdo,450) w1,w2,w3
450  format(' weights from lin fit: ',3f12.5)
  else
     x1 = (-b+dsqrt(disc))/(2d0*a)
     x2 = (-b-dsqrt(disc))/(2d0*a)
     x = x1
     if (dabs(x2-xlin) < dabs(x1-xlin)) x=x2
     if (x1*x2 < 0d0) then
        if (qval > q2) x = dmax1(x1,x2)
        if (qval <= q2) x = dmin1(x1,x2)
     endif
     w1 = 0.5d0*(x/def-1d0)*x/def
     w2 = 1d0-(x/def)**2
     w3 = 0.5d0*(x/def+1d0)*x/def
     if (ipr >= 10) write(stdo,460) w1,w2,w3
460  format(' weights from quad fit:',3f12.5)
  endif

  ! ... Interpolate sumev,sumqv
90 sqv = w1*sumqv(1,1) + w2*sumqv(2,1) + w3*sumqv(3,1) &
       + w1*sumqv(1,2) + w2*sumqv(2,2) + w3*sumqv(3,2)
  sav = w1*sumqv(1,1) + w2*sumqv(2,1) + w3*sumqv(3,1) &
       - w1*sumqv(1,2) - w2*sumqv(2,2) - w3*sumqv(3,2)
  sev = w1*sumev(1,1) + w2*sumev(1,2) + w3*sumev(1,3)
  sen = w1*sumev(2,1) + w2*sumev(2,2) + w3*sumev(2,3)
  if (ipr >= 20) write(stdo,470) &
       q1,q2,q3,sqv, &
       sumev(1,1),sumev(1,2),sumev(1,3),sev, &
       sumev(2,1),sumev(2,2),sumev(2,3),sen
470 format(' interpolated charge:  ',3f12.5,'  ->',f11.5 &
       /' eigenvalue sum:       ',3f12.5,'  ->',f11.5 &
       /' entropy term:         ',3f12.5,'  ->',f11.5)
  if (ipr >= 20 .AND. nsp == 2) write(stdo,471) a1,a2,a3,sav
471 format(' interpolated moment:  ',3f12.5,'  ->',f11.5)
  sumqv(1,1) = (sqv + sav)/nsp
  sumqv(1,2) = (sqv - sav)/nsp
  sumev(1,1) = sev
  sumev(2,1) = sen

  ! --- Interpolate coeffs of local density ---
  if (lrout > 0) then
     do  ib = 1, nbas

        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pos)
        call dcopy(i_copy_size,ssite(ib)%pos,1,p,1)


        rsm=sspec(is)%rsma
        lmxl=sspec(is)%lmxl
        lmxh=sspec(is)%lmxb


        lmxa=sspec(is)%lmxa
        kmax=sspec(is)%kmxt

        !       call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkapi)

        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2

        !   ... Case Pkl*Pkl, P*H, H*H
        if (lmxa > -1) then
           nelt1 = (kmax+1)*(kmax+1)*nlma*nlma
           nelt2 = (kmax+1)*nkaph*nlma*nlmh
           nelt3 = nkaph*nkaph*nlmh*nlmh
           call mshn31 ( w1 , w2 , w3 , nelt1 , nsp , sv_p_oqkkl( 1 , ib )%v &
                )

           call mshn31 ( w1 , w2 , w3 , nelt2 , nsp , sv_p_oqkkl( 2 , ib )%v &
                )

           call mshn31 ( w1 , w2 , w3 , nelt3 , nsp , sv_p_oqkkl( 3 , ib )%v &
                )

        endif

     enddo

     ! --- Interpolate smooth density ---
     call mshn32(w1,w2,w3,n1,n2,n3,k1,k2,k3,nsp,smrho)
     !     call zprm3('sm rho-out',0,smrho,k1,k2,k3*nsp)

     ! --- Interpolate forces ---
     if (lfrce /= 0) then
        do  ib = 1, nbas
           do  m = 1, 3
              f(m,ib,1) = w1*f(m,ib,1) + w2*f(m,ib,2) + w3*f(m,ib,3)
           enddo
        enddo
     endif
  endif

  ! --- Choose new energy window ---
  efold = ef0
  def1 = def
  eadd = dmin1(x,1.0d0)
  ef0 = ef0+eadd
  def = 2d0*dabs(eadd)
  def = dmin1(def,1.0d0,def1*4d0)
  def = dmax1(def,0.001d0,def1/4d0)
  if (ipr >= 20) write(stdo,455) efold,def1
455 format(' old energy window:   ef=',f10.5,'   plusminus',f9.5)
  if (ipr >= 10) write(stdo,456) ef0,def
456 format(' new energy window:   ef=',f10.5,'   plusminus',f9.5)
  if (ipl > 0) write(stdl,711) w1,w2,w3,ef0,def,sev,sen
711 format('nf w',3f7.3,'  ef',2f8.5,'   E0',f12.6,'  -sgS',f10.6)

end subroutine mshn3p


subroutine mshn31(w1,w2,w3,nelt,nsp,qkkl)

  !- Interpolate coeffs of local density
  !     implicit none
  ! ... Passed parameters
  integer :: nelt,nsp
  double precision :: w1,w2,w3,qkkl(nelt,3,nsp)
  ! ... Local parameters
  integer :: i,isp

  do  isp = 1, nsp
     do  i = 1, nelt
        qkkl(i,1,isp) = w1*qkkl(i,1,isp) + &
             w2*qkkl(i,2,isp) + &
             w3*qkkl(i,3,isp)
     enddo
  enddo
  if (nsp == 2) call dcopy(nelt,qkkl(1,1,2),1,qkkl(1,2,1),1)
end subroutine mshn31


subroutine mshn32(w1,w2,w3,n1,n2,n3,k1,k2,k3,nsp,smrho)

  !- Interpolate smooth mesh density
  !     implicit none
  ! ... Passed parameters
  integer :: n1,n2,n3,k1,k2,k3,nsp
  double precision :: w1,w2,w3
  double complex smrho(k1,k2,k3,3,nsp)
  ! ... Local parameters
  integer :: i1,i2,i3,isp

  do  isp = 1, nsp
     do  i3 = 1, n3
        do  i2 = 1, n2
           do  i1 = 1, n1
              smrho(i1,i2,i3,1,isp) = w1*smrho(i1,i2,i3,1,isp) &
                   + w2*smrho(i1,i2,i3,2,isp) &
                   + w3*smrho(i1,i2,i3,3,isp)
           enddo
        enddo
     enddo
  enddo
  if (nsp == 2) &
       call dcopy(2*k1*k2*k3,smrho(1,1,1,1,2),1,smrho(1,1,1,2,1),1)
end subroutine mshn32


