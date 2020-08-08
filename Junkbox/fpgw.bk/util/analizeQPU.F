     program aQPU

c 15 aug 2002
c compare QPU files

      implicit real*8 (a-h,o-z)
      dimension elda(ntq,nq),vxc(ntq,nq),sex(ntq,nq),sexcore(ntq,nq),
     i          rsec(ntq,nq),csec(ntq,nq), !sf..13May2002
     i          zfac(ntq,nq),itq(ntq),q(3,ntq,nq)

      dimension eqp(ntq,nq),eqp2(ntq,nq),wqp(ntq,nq)
      logical legas
c


      eshift0 =0d0
      eshift02=0d0
      eshtlda =0d0

      open(unit=3888,
     . file='~/sirun/si.sc.k=4x4x4dw=.0025omgc=.06isig=1/TOTE1run'
     .   form='formatted',status='old')
      ifqpe1=3888
      open(unit=3889,
     . file='~/sirun/si.sc.k=4x4x4dw=.0025omgc=.06isig=1/TOTE2run'
     .   form='formatted',status='old')
      ifqpe2=3889

      read(ifqpe1,*)
      read(ifqpe2,*)

      read (ifqpe1,*)
      read (ifqpe1,*)
      read (ifqpe2,*)
      read (ifqpe2,*)

      read(ifqpe1,"(a, 3d24.16,a)")eshtlda,eshift0,eshift02
      read(ifqpe1,"(a, 3d24.16,a)")eshtlda2,eshift02,eshift022

c loop over q-vector
      write (iftote, *) nq,ntq,ef
      write (iftote2,"(2i9,4d24.16)") 
     &  nq,ntq, ef*rydberg(), eshtlda, eshift0, eshift02
      write (ifqpe,*)
      write (ifqpe,"(a)")
     &'           q               state  SEx   SExcore SEc    vxc    dSE
     &  dSEnoZ  eLDA  eQPnoZ  eQPnoZ  eHF  Z=1  FWHM=2Z*Simg  ReS(elda)'
c
      do      iq = 1,nq
        do      it = 1,ntq
c
          eshift   = zfac(it,iq)
     .           * (sex(it,iq)+sexcore(it,iq)+rsec(it,iq)-vxc(it,iq)) !sf
          eqp(it,iq)  = elda(it,iq) + eshift + eshift0
          eshift2   = 
     .             (sex(it,iq)+sexcore(it,iq)+rsec(it,iq)-vxc(it,iq)) !sf
          eqp2(it,iq) = elda(it,iq) + eshift2 + eshift02
          fwhm  =  2d0*csec(it,iq) * zfac(it,iq)  !takao multiply zfac  !sf..13May2002
          ehf   =  elda(it,iq) + sex(it,iq)+ sexcore(it,iq) - vxc(it,iq)

          write(ifqpe,6100) q(1:3,it,iq),itq(it),sex(it,iq),sexcore(it,iq)
     &     ,rsec(it,iq),                                 !sf
     &     vxc(it,iq), eshift,eshift2,elda(it,iq)+eshtlda,
     &      eqp(it,iq), eqp2(it,iq), ehf,zfac(it,iq),fwhm, 
     .      sex(it,iq)+sexcore(it,iq)+rsec(it,iq)        !sf
 6100     format (3f9.5,1x,i2,1x,10f7.2,f5.2,f10.5,3x,f10.5)

          eqp01= elda(it,iq) + eshift
          eqp02= elda(it,iq) + eshift2
c
          write(iftote,"(3f12.7,1x,2i4,1x,4d24.16)") 
     &   q(1:3,it,iq),itq(it),iq, elda(it,iq), eqp01, eqp02, zfac(it,iq)

          write(iftote2,"(3f12.7,1x,2i4,1x,4d24.16)") 
     &   q(1:3,it,iq),itq(it),iq, elda(it,iq)+eshtlda, 
     &      eqp(it,iq),eqp2(it,iq), zfac(it,iq)
        end do
        write (ifqpe,*)
      end do
c
c------------------------------------------------------------------------
      INQUIRE (FILE = 'LEGAS', EXIST = legas)
      if(legas) then
        rydberginv=1d0/rydberg()
        hartreeinv=.5d0/rydberg()
        print *," EGAS mode eshift?"
        read(5,*) eshift0
        open(1111,file='egas.rlt')
        eshift0 = eshift0/rydberginv
        print *,' deltaw=',deltaw
        print *,' alat  =',alat
        print *,' ef    =',ef
c
        pi         = 4.d0*datan(1.d0)
        tpia       = 2.d0*pi/alat
c
        do      iq = 1,nq
          do      it = 1,ntq
            sm1     =  (sex(it,iq)+ rsec(it,iq)-vxc(it,iq) + eshift0) !sf
            sm2     =  (            csec(it,iq)           )           !sf
            zinv1=1d0 !sf   1d0-hartreeinv*(rsec(3,it,iq) - rsec(1,it,iq))/2d0/deltaw
            zinv2=0d0 !sf    -hartreeinv*(csec(3,it,iq) - csec(1,it,iq))/2d0/deltaw

            write (1111,6110) tpia*sqrt(sum(q(1:3,it,iq)**2)/ef)
     &    ,q(1:3,it,iq),itq(it),
     &    rydberginv*sm1,rydberginv*sm2,zinv1,zinv2,
     &    rydberginv/dcmplx(zinv1,zinv2)*dcmplx(sm1,sm2)
 6110       format (f10.5,2x,
     &     3f9.5,1x,i2,'  M=',2f12.5'  Zinv=',2f12.5,'  E=',2f12.5)
          end do
          write (1111,*)
        end do
      endif
c------------------------------------------------------------------------

c formats
c 6000 format (1x,'q =',)
      return
      end


