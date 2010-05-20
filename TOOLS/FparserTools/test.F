Cgetarg...info...           structure ['bndfp', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['bndfp', 'ham', 'sham']
Cgetarg...info...           structure ['bndfp', 'pot', 'spot']
Cgetarg...info...           structure ['bndfp', 'lat', 'slat']
Cgetarg...info...           structure ['bndfp', 'bz', 'sbz']
Cgetarg...info...           structure ['bndfp', 'site', 'ssite']
Cgetarg...info...           structure ['bndfp', 'spec', 'sspec']
      subroutine bndfp(nbas,nsp,nlibu,lmaxu,lldau,
     .  ssite,sspec,slat,sctrl,sham,spot,sbz,sstrn,
     .  ndham,leks,lrout,lfrce,lpnu,dmxp,iter,maxit,evl,frc,
     .  dmatu,vorb,llmfgw)
      do 2010 iq = iqini, iqend

        do  isp = 1, nsp
          if (onesp .eq. 0 .or. isp .eq. onesp) then

cccccccccccccc
ctakao do we need to do shorbz? --> furthermore, shorbz is probably wrong.---> see new algolism shortn3.
            qpo=qp
            if(.not.noshorbz()) then
              if (isp .eq. 1 .or. isp .eq. onesp) call shorbz(qp,qp,qlat,plat)
            endif


ccccccccccccccccccccccccc
c        qp=(/-1d0/3d0,1d0/3d0,1d0/3d0/)
c        qp=(/.5d0,.5d0,-.5d0/)
c        qp=(/-1d0,0d0,0d0/)
ccccccccccccccccccccccccc



C   ... For this qp, G vectors for PW basis and hamiltonian dimension
            if (pwemax .gt. 0 .and. mod(pwmode,10) .gt. 0) then
              pwgmin = dsqrt(pwemin)
              pwgmax = dsqrt(pwemax)
              call pshpr(1)
              call dpzero(xv,3)
              if (mod(pwmode/10,10) .eq. 1) call dpcopy(qp,xv,1,3,1d0)
              !print *,'xxx xxx 1111'
              call gvlst2(alat,plat,xv,0,0,0,pwgmin,pwgmax,0,
     .      0,0,napw,dum,dum,dum,dum)
              call poppr
              call defi(oigv2,3*napw)
              call defi(okv,3*napw)
              call pshpr(iprint()-10)
              !print *,'xxx xxx 2222'
              call gvlst2(alat,plat,xv,0,0,0,pwgmin,pwgmax,0,
     .      2,0,napw,w(okv),dum,dum,w(oigv2))
              !print *,'xxx xxx 3333'
              call rlse(okv)
              call poppr
              ndimh = ldim + napw  !ndimh is iq-dependent.
              if (mod(pwmode,10) .eq. 2) ndimh = napw
              if (ndimh .gt. ndham) then
                call fexit2(-1,111,'%N Exit -1 : BNDFP: '//
     .        'ndimh=%i exceeds ndham=%i.  Try increasing '//
     .        'input NPWPAD',ndimh,ndham)
              endif
            else
              ndimh = ldim
              oigv2 = 1
            endif
            ndimhx = ndimh*nspc !this is iq-dependent.

C   ... Make Hamiltonian and overlap matrices
Cgetarg...info...           integer :: nqsig , integer(8) :: sham%nqsig 1
Cgetarg...info...           integer :: oqsig , integer(8) :: sham%oqsig 1
Cgetarg          call upack2('ham nqsig oqsig',sham,nqsig,oqsig)

            nqsig=sham%nqsig
            oqsig=sham%oqsig

            if (oqsig .eq. 0) oqsig = 1

            ispc = min(isp,nspc) !nspc=2 for non-colinier case
            !nspc=1


c ---
            if(isp==1) then
              if(allocated(hamm)) deallocate(hamm)
              if(allocated(ovlm)) deallocate(ovlm)
c      if(allocated(hamsom)) deallocate(hamsom)
              nnspc=1
              if(nspc==2) nnspc=4
              allocate(hamm(ndimh,ndimh,nnspc),ovlm(ndimh,ndimh,nnspc))
c      if(lso/=0) then
c        allocate(hamsom(ndimh,ndimh))
c      else
c        allocate(hamsom(1,1)) !dummy
c      endif
            endif
            if(ispc==1) then
              hamm=0d0
              ovlm=0d0
            endif


            i = lrsig*10
C Commented by takao lwsig should be renewed.
C C       lwsig=-1: sigm LDA -> orbital basis
C C       lwsig= 1: sigm orbital -> LDA basis
C C       lwsig= 2: sigm orbital -> LDA basis, high energy parts replaced
C C       lwsig= 3: Write evecs of LDA hamiltonian to file
C C       lwsig= 4: Write evecs of hamiltonian to file
C C       lwsig= 5: sigm orbital -> orbital, high energy parts replaced
C C       hambls returns with evecs,evals of the LDA hamiltonian in s,h
C C       Transform sigm to LDA basis: hambls returns sigm(LDA) in h
C         if (lwsig .eq. 1) then
C           i = i + 3000
C         elseif (lwsig .eq. 2) then
C           i = i + 4000
C         elseif (lwsig .eq. -1 .or. lwsig .eq. 3) then
C           i = i + 1000
C           if (nspc .eq. 2)
C      .      call rx('transform sigm not implemented in noncoll case')
C         elseif (lwsig .eq. -1 .or. lwsig .eq. 4) then
C           i = i + 6000
C         elseif (lwsig .eq. 5) then
C           i = i + 5000
C         endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ctakaox
c        if(iq==1) then
c              ifio = fopna('ovlmat',-1,4)
c        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *,' mmm: qp=',qp


            if(debug) print *,' mmm1: qp=',qp
            if(sigmamode) then
c --- Generate interpolation weight wgtq.
              allocate( wgtq(nbas,nbas,nk1,nk2,nk3) )
              pi2i = 2d0*4d0*datan(1d0)*(0d0,1d0)
              wgtq = 0d0
              !print *,' mmm2: qp=',qp
              do ib1=1,nbas
                do ib2=1,nbas
                  do ii = 1,npair(ib1,ib2)
                    rrr = matmul(plat, nlat(:,ii,ib1,ib2))
ccccccccccccccccccccccccccccccccccc
c          write(6,"(a, 3i5,2x,3f8.3)")'vvv1: qwgt=', ib1,ib2,ii
ccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccc
c          qpx=(/-0.50000D+00, -0.50000D+00, -0.50000D+00/)
ccccccccccccccccccccccccccccccccccccccccccc

c wgtq is interpolation weight for MTO.
                    do ik1=1,nk1
                      do ik2=1,nk2
                        do ik3=1,nk3
                       wgtq(ib1,ib2,ik1,ik2,ik3) = wgtq(ib1,ib2,ik1,ik2,ik3) 
     &     + qwgt(ii,ib1,ib2)*exp( -pi2i*sum( (qp-qsmesh(:,ik1,ik2,ik3))*rrr ) )
c     &      + qwgt(ii,ib1,ib2)*exp( -pi2i*sum( (qpx-qsmesh(:,ik1,ik2,ik3))*rrr ) )
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
 888              continue
                enddo
              enddo
              deallocate(igx)
 1019         continue
            endif 
 2010 continue
      end subroutine bndfp
