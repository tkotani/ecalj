      subroutine vsl(job,rmtg,rsmg,kmax,lmaxg,pmax,lmaxh,k0,ndim,p0,
     .     sf,fval,fslo,flap,pkl,gpkl)
C- lm-resolved value, slope and Laplacian of function F_L at |r| = rmt 
C  derived from its P_kM expansion 
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :if 1s digit = 1, returns pkl and gpkl; otherwise takes them on input
Ci         :if 10s digit = 0, all radial functions (pkl,gpk,fval,fslo, and flap)
Ci         :are divided by r^l, otherwise not. 
Ci         :if 1s digit = 0 then 10s digit is ignored; whether or not the radial   
Ci         :functions divided by r^l is defined by input arrays pkl and gpkl  
Ci   rmtg  :l-dependent radii at which value, slope, and Laplacian of F_L are 
Ci         :evaluated. Must either be specified for 0..lmaxg, or rmtg(0) could
Ci         :be set to some negative number. 
Ci         :The latter case implies that rmtg(l) = - rmtg(0) for all l.
Ci   rsmg  :smoothing radius for polynomials
Ci   kmax  :polynomial cutoff
Ci   lmaxg :l-cutoff for P_kM expansion
Ci   pmax  :max power of the Laplace operator (lap)^p F_L := F_pL, must be >=1
Ci   lmaxh :l-cutoff for F_L 
Ci   k0,ndim,p0:leading dimensions of the array of coefficients sf
Ci   sf    :coefficients sf(k,M,p,L) of decomposition of head function F_pL
Ci         :into auxiliary polynomials P_kM 
Ci         :here M=1..(lmaxg+1)^2, L=1..(lmaxh+1)^2, k=0..kmax, p=0..pmax 
Co Outputs
Co   fval   : lm-resolved value of F_L       |
Co   fslo   : lm-resolved slope of F_L       | at rmt (see Remarks)
Co   flap   : lm-resolved Laplacian of F_L   |
Co   pkl    : radial part of polynomials P_kl / r^l (if job's 10s digit = 0)
Co   gpkl   : radial part of [d/dr (P_kl)] / r^l    (  -    "    -         )
Co
Cr Remarks
Cr   fval, fslo, and flap are defined so that: 
Cr             F_L(r-R) = \sum_M fval(M,L) * r^m * Y_M(r)  | if job's 10s digit = 0  
Cr      (d/dr) F_L(r-R) = \sum_M fslo(M,L) * r^m * Y_M(r)  | otherwise without r^m
Cr         lap F_L(r-R) = \sum_M flap(M,L) * r^m * Y_M(r)  |
Cr   are (approximately) fulfilled at |r| = rmt
Cr  
Cr   in other words, fval(M,L), fslo(M,L), and flap(M,L) are M-resolved    
Cr   value, slope, and laplacian at the expansion site of F_L at the head site 
Cr  
Cr   fval, fslo, and flap correspond to E^(0)_LL', E^(1)_LL', and 
Cr   E^(2)_LL' in MvS's notes on interstitial fits
Cr
Cb Bugs
Cu Updates
Cu   15 May 08 rmtg made l-dependent (S. Lozovoi)
Cu   26 May 07 first written (S. Lozovoi)
C ----------------------------------------------------------------------
C     implicit none
C Input parameters
      integer job,kmax,lmaxg,pmax,lmaxh,k0,ndim,p0
      double precision rsmg,rmtg(0:lmaxg)
      double precision sf(0:k0,ndim,0:p0,(lmaxh+1)**2)
C Output parameters
      double precision fval(ndim,(lmaxh+1)**2),
     .   fslo(ndim,(lmaxh+1)**2),flap(ndim,(lmaxh+1)**2)
      double precision pkl(0:k0,0:lmaxg),gpkl(0:k0,0:lmaxg)
C Local variables
      integer kmax0,n0,jvec      
      parameter (n0 = 10, kmax0 = 20)
      integer ilm,jlm,ik,il,im             
      double precision sval,sslo,slap
c     double precision rsmgl(0:n0),rsmg0,rsx
      double precision rmtgl(0:n0),rmtg0,rmx
      double precision wkl(1,0:kmax0,0:n0+1),gwkl(1,0:kmax0,0:n0),wk(1)
      double precision tol
      data tol/1.d-15/

      if (pmax .lt. 1) call rx('vsl: pmax should be at least 1')
      if (lmaxg .gt. n0) call rx('vsl: lmaxg gt n0')
      if (kmax .gt. kmax0) call rx('vsl: kmax gt kmax0')
            

C --- Handle negative rmt
c     if (rsmg(0) .lt. 0d0) then
c       call dvset(rsmgl(0),1,lmaxg+1,-rsmg(0))
c     else
c       call dcopy(lmaxg+1,rsmg(0),1,rsmgl(0),1)
c     endif
      if (rmtg(0) .lt. 0d0) then
        call dvset(rmtgl(0),1,lmaxg+1,-rmtg(0))
      else
        call dcopy(lmaxg+1,rmtg(0),1,rmtgl(0),1)
      endif

c --- If job=1, make polynomials pkl and their radial derivatives gpkl 
      if (mod(job,10) .eq. 1) then
        if (mod(job/10,10) .eq. 0) then
          jvec = 10
        else  
          jvec = 11
        endif
        rmx = -1d2
        do il = lmaxg, 0, -1
          rmtg0 = rmtgl(il)
          if (dabs(rmtg0-rmx) .gt. tol) then

c ... need lmax+1 for gradients
            call vecpkl(rmtg0,rsmg,1,kmax,il+1,1,kmax0,
     .        wk,jvec,wkl,gwkl)

            rmx = rmtg0
          endif
c         do ik = 0, kmax
c           pkl(ik,il) = wkl(1,ik,il)
c           gpkl(ik,il) = gwkl(1,ik,il)
c         enddo
          call dcopy(kmax+1,wkl(1,0,il),1,pkl(0,il),1)
          call dcopy(kmax+1,gwkl(1,0,il),1,gpkl(0,il),1)
        enddo
      endif 

c --- Sum up polynomials to make value, slope, and Laplacian
      do jlm= 1, (lmaxh+1)**2
        ilm = 0
        do il = 0, lmaxg
          do im = -il, il 
            ilm = ilm + 1
            sval = 0d0
            sslo = 0d0
            slap = 0d0
            do ik = 0, kmax
              sval = sval + sf(ik,ilm,0,jlm)*pkl(ik,il)
              sslo = sslo + sf(ik,ilm,0,jlm)*gpkl(ik,il)
              slap = slap + sf(ik,ilm,1,jlm)*pkl(ik,il)
            enddo
            fval(ilm,jlm) = sval
            fslo(ilm,jlm) = sslo
            flap(ilm,jlm) = slap
          enddo
        enddo
      enddo

      end


