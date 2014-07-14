      subroutine qpe1_mlw (ifqpe,iftote,iftote2,itq,q,
     i                elda_in,vxc_in,sex_in,sexcore_in,
     i                rsec_in,csec_in,jin,deltaw,alat,ef,omg,
     i                umn,iko_ix,iko_fx,iwini,iwfin,
     d                nband,ntq,ntw,nq,is,
     o                eshift0,eshift02,eshtlda)

c 92.04.18
c calculates the quasiparticle energies
c E(k,t) = e(k,t) + Z [SEx(k,t) + SEc(k,t) - xcLDA(k,t)]
c e(k,t) = LDA eigenvalue
c Z      = [1 - dSEc(e(k,t))/dw]^(-1)
c SEx(k,t)   = <psi(k,t)| SEx |psi(k,t)>
c SEc(k,t)   = <psi(k,t)| SEc |psi(k,t)>, SEc = GWc
c xcLDA(k,t) = <psi(k,t)| vxc |psi(k,t)>
      
C jin    0 make no shifts in eigenvalues
C   999999 use shists as given
C        >0 make shifts so eval(jin)=0 (both LDA,QP)
C        <0 use shifts by EFERMI files

c ifqpe = unit output file
c itq   = state label
c q     = q-vectors
c elda  = e(k,t)
c vxc   = xcLDA
c sex   = SEx
c sec   = SEc
c zfac  = Z
c eshift0 = energy shift for quasiparticle (eV)
c ntq,nq,is = no. t,k and spin label

c eqp   = E(k,t) in eV
c wqp   = quasiparticle width

      implicit real*8 (a-h,o-z)
      dimension elda_in(nband,nq),vxc_in(ntw,nq),
     i          sex_in(ntw,nq),sexcore_in(ntw,nq),
     i          rsec_in(iwini:iwfin,ntw,nq),csec_in(iwini:iwfin,ntw,nq),
     i          elda(ntq,nq),vxc(ntq,nq),sex(ntq,nq),sexcore(ntq,nq),
     i          rsec(iwini:iwfin,ntq,nq),csec(iwini:iwfin,ntq,nq),
     i          itq(ntq),q(3,nq)
      complex*8  umn(iko_ix:iko_fx,ntw,nq)

      dimension eqp(ntq,nq),eqp2(ntq,nq),wqp(ntq,nq),omg(iwini:iwfin)
      logical legas,lprev
c
!      if(jin==999999) then
!c         print *,' jin=999999 use given shift '
!      elseif(jin>0) then
!        iq = jin/ntq + 1
!        it = jin - (iq-1)*ntq
!        print *,' eshift is so as to be zero for UP of (it iq)= ',it,iq
!        eshift  = zfac(it,iq)
!     .           * (sex(it,iq)+sexcore(it,iq)+rsec(2,it,iq)-vxc(it,iq))
!        eshift0 = - elda(it,iq) - eshift
!        eshift  = 
!     .             (sex(it,iq)+sexcore(it,iq)+rsec(2,it,iq)-vxc(it,iq))
!        eshift02 = - elda(it,iq) - eshift
!        eshtlda  = - elda(it,iq)
!      elseif (jin<0) then
!        ifeflda  = iopenxx ('EFERMI')
!        ifefqp   = iopenxx ('EFERMI.QP')
!        ifefqp1  = iopenxx ('EFERMI.QPz=1')
!        read(ifeflda,*) eftetlda
!        read(ifefqp,* ) eftetqp
!        read(ifefqp1,*) eftetqp1
!        close(ifeflda)
!        close(ifefqp)
!        close(ifefqp1)
!        eshift0  = (-eftetqp  + ef)*rydberg()
!        eshift02 = (-eftetqp1 + ef)*rydberg()
!        eshtlda  = (-eftetlda + ef)*rydberg()
!      else
!       eshift0 =0d0
!       eshift02=0d0
!       eshtlda =0d0
!      endif
c

c LDA energy
      do ip = 1,nq
      do i  = 1,ntq
        elda(it,ip) = (elda_in(itq(i),ip)-ef)*rydberg()
      enddo  
      enddo  

c from Wannier to Kohn-Sham
c static part
      call wan2ks(sex_in,sex,umn,itq,iko_ix,iko_fx,ntw,ntq,nq)
      call wan2ks(sexcore_in,sexcore,umn,itq,iko_ix,iko_fx,ntw,ntq,nq)
      call wan2ks(vxc_in,vxc,umn,itq,iko_ix,iko_fx,ntw,ntq,nq)
      call wan2ks2(rsec_in,rsec,umn,itq,iko_ix,iko_fx,
     & iwini,iwfin,ntw,ntq,nq)
      call wan2ks2(csec_in,csec,umn,itq,iko_ix,iko_fx,
     & iwini,iwfin,ntw,ntq,nq)

      if (is .eq. 1) then
      write (ifqpe,*)
     .'==============================================================='
      write (ifqpe,*) ' quasiparticle energies MAJORITY'
      write (ifqpe,*)
     .'==============================================================='
!      write(ifqpe,"('E_shift=', 3d24.16,' eV')")eshtlda,eshift0,eshift02
      endif

      if (is .eq. 2) then
      write (ifqpe,*)
     .'==============================================================='
      write (ifqpe,*) ' quasiparticle energies MINORITY'
      write (ifqpe,*)
     .'==============================================================='
!      write(ifqpe,"('E_shift=', 3d24.16,' eV')")eshtlda,eshift0,eshift02
      endif

c loop over q-vector
      write (iftote, *) nq,ntq,ef
      write (iftote2,"(2i9,4d24.16)") 
     &  nq,ntq, ef*rydberg(), eshtlda, eshift0, eshift02 
      write (ifqpe,*) 
      write (ifqpe,"(a)")
     &'           q               state  SEx   SExcore SEc    vxc    dSE'
     &//
     &'  dSEnoZ  eLDA    eQP  eQPnoZ   eHF  Z    FWHM=2Z*Simg  ReS(elda)'
c
      do      iq = 1,nq
      do      it = 1,ntq
c
         do iw = iwini,iwfin
            omega = omg(iw)
            dw2 = elda(it,iq) -omega
     .          + (sex(it,iq)+sexcore(it,iq)+rsec(2,it,iq)-vxc(it,iq))
            if (dw2.le.0.0d0) then
              iw1 = iw-1
              iw2 = iw
              w1 = abs(dw2)/ (dw1+abs(dw2))
              w2 = abs(dw1)/ (dw1+abs(dw2))
              rsc = w1*rsec(iw1,it,iq) + w2*rsec(iw2,it,iq)
              csc = w1*csec(iw1,it,iq) + w2*csec(iw2,it,iq)
            endif
            dw1 = dw2
         enddo ! iw   

        eshift   =  sex(it,iq)+sexcore(it,iq)+sec-vxc(it,iq)
        eqp(it,iq)  = elda(it,iq) + eshift !+ eshift0
        eshift2   = eshift
        eqp2(it,iq) = elda(it,iq) + eshift2 !+ eshift02
        fwhm  =  2d0*csc
        ehf   =  elda(it,iq) + sex(it,iq)+ sexcore(it,iq) - vxc(it,iq)

        write(ifqpe,6100) q(1:3,iq),itq(it),sex(it,iq),sexcore(it,iq)
     &     ,rsec(2,it,iq),
     &     vxc(it,iq), eshift,eshift2,elda(it,iq),
     &      eqp(it,iq), eqp2(it,iq), ehf,zfac,fwhm, 
     .      sex(it,iq)+sexcore(it,iq)+rsec(2,it,iq)
 6100   format (3f9.5,1x,i2,1x,10f7.2,f5.2,f10.5,3x,f10.5)

        eqp01= elda(it,iq) !+ eshift
        eqp02= elda(it,iq) !+ eshift2
        write(iftote,"(3f12.7,1x,2i4,1x,4d24.16)") 
     &   q(1:3,iq),itq(it),iq, elda(it,iq), eqp01, eqp02, zfac

        write(iftote2,"(3f12.7,1x,2i4,1x,4d24.16)") 
     &   q(1:3,iq),itq(it),iq, elda(it,iq)+eshtlda, 
     &      eqp(it,iq),eqp2(it,iq), zfac
      end do
      write (ifqpe,*) 
      end do

c------------------------------------------------------------------------

c formats
c 6000 format (1x,'q =',)
      return
      end
c------------------------------------------------------------------------
      subroutine wan2ks(se_in,se,umn,itq,iko_ix,iko_fx,ntw,ntq,nq)
      implicit real*8(a-h,o-z)

      integer(4) :: itq(ntq)
      real(8) :: se_in(ntw,nq),se(ntq,nq)
      complex(8) :: umn(iko_ix:iko_fx,ntw,nq)

      se = 0d0
      do i  = 1,ntq
        it = itq(i)
        if (it.ge.iko_ix.and.it.le.iko_fx) then
          do ip = 1,nq
          do iwf = 1,ntw
            se(i,ip) = se(i,ip) + abs(umn(it,iwf,ip))**2*se_in(iwf,ip)
          enddo
          enddo ! ip
        endif
      enddo ! i

      return
      end
c------------------------------------------------------------------------
      subroutine wan2ks2(se_in,se,umn,itq,iko_ix,iko_fx,
     &  iwini,iwfin,ntw,ntq,nq)
      implicit real*8(a-h,o-z)

      integer(4) :: itq(ntq)
      real(8) :: se_in(iwini:iwfin,ntw,nq),se(iwini:iwfin,ntq,nq)
      complex(8) :: umn(iko_ix:iko_fx,ntw,nq)

      se = 0d0
      do i  = 1,ntq
       it = itq(i)
       if (it.ge.iko_ix.and.it.le.iko_fx) then
         do iw = iwini,iwfin
         do ip = 1,nq
         do iwf = 1,ntw
           se(iw,i,ip) = se(iw,i,ip)
     .    + abs(umn(it,iwf,ip))**2*se_in(iw,iwf,ip)
         enddo
         enddo ! ip
         enddo ! iw
       endif
      enddo ! i   

      return
      end
c------------------------------------------------------------------------


