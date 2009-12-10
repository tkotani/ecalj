      subroutine makalp(nl,nda,nbas,nkap,e,hcr,io,ips,lmx,ac,alpha)
C- Generates screening parameters alpha
C ----------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1.
Ci         :alpha made for channels l=0..nl-1 for which hcr>0
Ci   nda   :leading dimension of alpha
Ci   nbas  :size of basis
Ci   nkap  :number of energies for which alpha is made
Ci   nkap  :number of kappa's
Ci   e     :Hankel energies, atomic units
Ci   hcr   :l-dependent hard core screening radii in atomic units
Ci   io:  1s digit:
Ci        0: use alpha's as passed by ac
Ci        1: make alpha fit to val or val,slo at hcr.
Ci           2-kappa case:
Ci                                  -1
Ci                   (  k_1   k_2 )    (  j_1   j_2 )
Ci           alpha = (            )    (            )
Ci                   ( k'_1  k'_2 )    ( j'_1  j'_2 )
Ci
Ci        2: make alpha to transform functions to val-slope functions
Ci           2-kappa case
Ci            (k_1  k_2 ) alpha =  ( val slo ) in 2-kap
Ci                                  -1
Ci                   (  k_1   k_2 )    (  k'_2  -k_2 )
Ci           alpha = (            )  = (             ) / (k_1 k'_2 - k'_1 k_2)
Ci                   ( k'_1  k'_2 )    ( -k'_1   k_1 )
Ci
Ci        Add 100's digit for OKA conventions
Ci
Ci   ac:  input alpha's to be transferred to alpha (iopt = 0 only)
Co Outputs
Co   alpha: vector of screening parameters.
Cl Local variables
Cr Remarks
Cr   For a single hankel function, alpha is stored as alpha(lm,site)
Cr   (single value for a given site and l channel)
Cr   For a double kappa hankel, alpha is stored as alpha(lm,site,i,j)
Cr   where i,j is a 2 by 2 matrix; Additionally the determinants
Cr   stored in 1,3 and 2,3
Cr
Cr   To generate the inverse of alpha, follow with a call to invalp.
Cu Updates
Cu   06 Aug 06 Redesigned
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,nda,nbas,nkap,io,ips(1),lmx(1)
      double precision e(2)
      double precision alpha(nda,nbas,nkap,nkap),
     .  ac(0:nl-1,1),hcr(nl,*)
C Local parameters
      integer lmax,ib,l,lm,m,is,ikap,jkap,iprint,stdo,
     .  io0,io2,nglob
      double precision ak(0:9,2),aj(0:9,2),dk(0:9,2),dj(0:9,2),detk
      double precision akl(0:9),ajl(0:9),dkl(0:9),djl(0:9)

      io0 = mod(io,10)
      io2 = mod(io/100,10)
      lmax = nl-1
      if (nl .eq. 5)
     .  print *, 'makalp: (warning) alphas need checking for nl=5'

C --- For each atom in the basis ... ---
      do  ib = 1, nbas
        is = ips(ib)

C   --- Get information for hcr, conventional Hankels ---
        if (io0 .ne. 0) then
          m = io2*10+0
          do  l = 0, lmax
          do  ikap = 1, nkap
            call radkj(e(ikap),hcr(l+1,is),lmax,akl,ajl,dkl,djl,m)
            ak(l,ikap) = akl(l)
            aj(l,ikap) = ajl(l)
            dk(l,ikap) = dkl(l)
            dj(l,ikap) = djl(l)
          enddo
          enddo
        endif

C   --- For each orbital make the appropriate alpha ---
        lm = 0
        do  l = 0, lmax
          do  m = -l, l
          lm = lm+1

          if (hcr(l+1,is) .ne. 0) then
          if (io0 .eq. 0 .and. nkap .eq. 1) then
            alpha(lm,ib,1,1) = ac(l,is)
          elseif (io0 .eq. 0) then
            call rxi('makalp not implemented for io=',io)
          else if (io0 .eq. 1 .and. nkap .eq. 1) then
            alpha(lm,ib,1,1) =  +aj(l,1)/ak(l,1)
          else if (io0 .eq. 1 .and. nkap .eq. 2) then
            detk = +ak(l,1)*dk(l,2) - dk(l,1)*ak(l,2)
            alpha(lm,ib,1,1) = (+aj(l,1)*dk(l,2) - dj(l,1)*ak(l,2))/detk
            alpha(lm,ib,2,1) = (-aj(l,1)*dk(l,1) + dj(l,1)*ak(l,1))/detk
            alpha(lm,ib,1,2) = (+aj(l,2)*dk(l,2) - dj(l,2)*ak(l,2))/detk
            alpha(lm,ib,2,2) = (-aj(l,2)*dk(l,1) + dj(l,2)*ak(l,1))/detk
          else if (io0 .eq. 2 .and. nkap .eq. 1) then
            alpha(lm,ib,1,1) = 1/ak(l,1)
          else if (io0 .eq. 2 .and. nkap .eq. 2) then
            detk = +ak(l,1)*dk(l,2) - dk(l,1)*ak(l,2)
            alpha(lm,ib,1,1) = +dk(l,2)/detk
            alpha(lm,ib,1,2) = -ak(l,2)/detk
            alpha(lm,ib,2,1) = -dk(l,1)/detk
            alpha(lm,ib,2,2) = +ak(l,1)/detk
          endif
C         Not ready to handle hcr=0 yet
          else
            call rx('makalp not ready for hcr=0')
          endif
          enddo
        enddo
      enddo

C  --- Output ---
      if (iprint() .ge. 50) then
        stdo = nglob('stdo')
        write(stdo,'(/'' MAKALP: alpha, mode'',i2/6x,''l=0 ...'')') io0
        do  ib = 1, nbas
          write(stdo,'(''# Site'',i4)') ib
          do  ikap = 1, nkap
          do  jkap = 1, nkap
            write(stdo,885) (alpha(lm**2,ib,ikap,jkap), lm=1,nl)
          enddo
          enddo
        enddo
      endif

  885 format(5f12.6)

      end

