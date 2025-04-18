subroutine qpe1_sc (ifqpe,iftote,iftote2,itq,q, &
     elda,vxc,sex,sexcore, &
     rsec,csec,jin,deltaw,alat,ef, &! & sf..13may remove zfac after csec
     ntq,nq,is, &
     eshift0,eshift02,eshtlda)
  !     o                eqp,wqp)

  ! 92.04.18
  ! calculates the quasiparticle energies
  ! E(k,t) = e(k,t) + Z [SEx(k,t) + SEc(k,t) - xcLDA(k,t)]
  ! e(k,t) = LDA eigenvalue
  ! Z      = [1 - dSEc(e(k,t))/dw]^(-1)
  ! SEx(k,t)   = <psi(k,t)| SEx |psi(k,t)>
  ! SEc(k,t)   = <psi(k,t)| SEc |psi(k,t)>, SEc = GWc
  ! xcLDA(k,t) = <psi(k,t)| vxc |psi(k,t)>

  ! jin    0 make no shifts in eigenvalues
  !   999999 use shists as given
  !        >0 make shifts so eval(jin)=0 (both LDA,QP)
  !        <0 use shifts by EFERMI files

  ! ifqpe = unit output file
  ! itq   = state label
  ! q     = q-vectors
  ! elda  = e(k,t)
  ! vxc   = xcLDA
  ! sex   = SEx
  ! sec   = SEc
  ! zfac  = Z
  ! eshift0 = energy shift for quasiparticle (eV)
  ! ntq,nq,is = no. t,k and spin label

  ! eqp   = E(k,t) in eV
  ! wqp   = quasiparticle width

  implicit real*8 (a-h,o-z)
  dimension elda(ntq,nq),vxc(ntq,nq),sex(ntq,nq),sexcore(ntq,nq), &
       rsec(ntq,nq),csec(ntq,nq), zfac(ntq,nq),itq(ntq),q(3,ntq,nq)
  dimension eqp(ntq,nq),eqp2(ntq,nq),wqp(ntq,nq)
  logical :: legas
  integer:: ifqpe,iftote,iftote2,itq,ntq,jin,nq,is,if1111,ifelda,ifepq,iq,it, ifeflda,ifefqp,ifefqp1
  zfac=1d0 !sf..13may2002 zfac=1 to not change to much the rest of this code
  ! in SC GW
  if(jin==999999) then
     !         print *,' jin=999999 use given shift '
  elseif(jin>0) then
     iq = jin/ntq + 1
     it = jin - (iq-1)*ntq
     print *,' eshift is so as to be zero for UP of (it iq)= ',it,iq
     eshift  = zfac(it,iq) &
          * (sex(it,iq)+sexcore(it,iq)+rsec(it,iq)-vxc(it,iq)) !sf
     eshift0 = - elda(it,iq) - eshift
     eshift  = &
          (sex(it,iq)+sexcore(it,iq)+rsec(it,iq)-vxc(it,iq)) !sf
     eshift02 = - elda(it,iq) - eshift
     eshtlda  = - elda(it,iq)
  elseif (jin<0) then
     open(newunit=ifeflda  ,file='EFERMI')
     open(newunit=ifefqp   ,file='EFERMI.QP')
     open(newunit=ifefqp1  ,file='EFERMI.QPz=1')
     read(ifeflda,*) eftetlda
     read(ifefqp,* ) eftetqp
     read(ifefqp1,*) eftetqp1
     close(ifeflda)
     close(ifefqp)
     close(ifefqp1)
     eshift0  = (-eftetqp  + ef)*rydberg()
     eshift02 = (-eftetqp1 + ef)*rydberg()
     eshtlda  = (-eftetlda + ef)*rydberg()
  else
     eshift0 =0d0
     eshift02=0d0
     eshtlda =0d0
  endif

  if (is == 1) then
     write (ifqpe,*) &
          '==============================================================='
     write (ifqpe,*) ' quasiparticle energies MAJORITY'
     write (ifqpe,*) &
          '==============================================================='
     write(ifqpe,"('E_shift=', 3d24.16,' eV')")eshtlda,eshift0,eshift02
  endif

  if (is == 2) then
     write (ifqpe,*) &
          '==============================================================='
     write (ifqpe,*) ' quasiparticle energies MINORITY'
     write (ifqpe,*) &
          '==============================================================='
     write(ifqpe,"('E_shift=', 3d24.16,' eV')")eshtlda,eshift0,eshift02
  endif

  ! loop over q-vector
  write (iftote, *) nq,ntq,ef
  write (iftote2,"(2i9,4d24.16)") &
       nq,ntq, ef*rydberg(), eshtlda, eshift0, eshift02
  write (ifqpe,*)
  write (ifqpe,"(a)") &
       !     .'           q               state  SEx   SExcore SEc    vxc    dSE
       !     .  dSEnoZ  eLDA  eQPnoZ  eQPnoZ  eHF  Z=1  FWHM=2Z*Simg  ReS(elda)'
       ! '           q               state  SEx   SExcore SEc    vxc   ---' &
       ! // '   dSEnoZ  eQP(starting by lmf)'// &
       ! '  eHF  Z=1  FWHM=2Z*Simg ReS(elda)'
       '           q               state   SEx    SExcore   SEc     vxc     ---'// &
       '    dSEnoZ   eQP(starting by lmf)    eHF   Z=1  FWHM=2Z*Simg ReS(elda)'

  do      iq = 1,nq
     do      it = 1,ntq

        eshift   = zfac(it,iq) &
             * (sex(it,iq)+sexcore(it,iq)+rsec(it,iq)-vxc(it,iq)) !sf
        eqp(it,iq)  = elda(it,iq) + eshift + eshift0
        eshift2   = &
             (sex(it,iq)+sexcore(it,iq)+rsec(it,iq)-vxc(it,iq)) !sf
        eqp2(it,iq) = elda(it,iq) + eshift2 + eshift02
        fwhm  =  2d0*csec(it,iq) * zfac(it,iq)  !takao multiply zfac  !sf..13May2002
        ehf   =  elda(it,iq) + sex(it,iq)+ sexcore(it,iq) - vxc(it,iq)


        if(elda(it,iq)<1d20) then
           ehfx = ehf
           dsenoz = eshift2
           if(abs(sex(it,iq))+abs(sexcore(it,iq))+abs(rsec(it,iq))==0d0) then
              dsenoz= 0d0
              ehfx  = 0d0
           endif
           write(ifqpe,6100) q(1:3,it,iq),itq(it),sex(it,iq),sexcore(it,iq) ,rsec(it,iq),&
                vxc(it,iq), 0d0, dsenoz, elda(it,iq)+eshtlda, &
                0d0, 0d0, ehfx,zfac(it,iq),fwhm, &
                sex(it,iq)+sexcore(it,iq)+rsec(it,iq) 
6100       format (3f9.5,1x,i3,1x,10f8.3,f5.2,f10.5,3x,f10.5)
        endif

        eqp01= elda(it,iq) + eshift
        eqp02= elda(it,iq) + eshift2

        write(iftote,"(3f12.7,1x,2i4,1x,4d24.16)") &
             q(1:3,it,iq),itq(it),iq, elda(it,iq), eqp01, eqp02, zfac(it,iq)

        write(iftote2,"(3f12.7,1x,2i4,1x,4d24.16)") &
             q(1:3,it,iq),itq(it),iq, elda(it,iq)+eshtlda, &
             eqp(it,iq),eqp2(it,iq), zfac(it,iq)
     end do
     write (ifqpe,*)
  end do

  !------------------------------------------------------------------------
  INQUIRE (FILE = 'LEGAS', EXIST = legas)
  if(legas) then
     rydberginv=1d0/rydberg()
     hartreeinv=.5d0/rydberg()
     print *," EGAS mode eshift?"
     read(5,*) eshift0
     open(newunit=if1111,file='egas.rlt')
     eshift0 = eshift0/rydberginv
     print *,' deltaw=',deltaw
     print *,' alat  =',alat
     print *,' ef    =',ef

     pi         = 4.d0*datan(1.d0)
     tpia       = 2.d0*pi/alat

     do      iq = 1,nq
        do      it = 1,ntq
           sm1     =  (sex(it,iq)+ rsec(it,iq)-vxc(it,iq) + eshift0) !sf
           sm2     =  (            csec(it,iq)           )           !sf
           zinv1=1d0 !sf   1d0-hartreeinv*(rsec(3,it,iq) - rsec(1,it,iq))/2d0/deltaw
           zinv2=0d0 !sf    -hartreeinv*(csec(3,it,iq) - csec(1,it,iq))/2d0/deltaw

           write (if1111,6110) tpia*sqrt(sum(q(1:3,it,iq)**2)/ef) &
                ,q(1:3,it,iq),itq(it), &
                rydberginv*sm1,rydberginv*sm2,zinv1,zinv2, &
                rydberginv/dcmplx(zinv1,zinv2)*dcmplx(sm1,sm2)
6110       format (f10.5,2x, &
                3f9.5,1x,i2,'  M=',2f12.5'  Zinv=',2f12.5,'  E=',2f12.5)
        end do
        write (if1111,*)
     end do
  endif
  !------------------------------------------------------------------------

  ! formats
  ! 6000 format (1x,'q =',)
  return
end subroutine qpe1_sc


