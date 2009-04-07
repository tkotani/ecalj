      subroutine addtog(nds,nda,iax,npr,sgalph,nprs,sg)
C- Appends structure constant vector for single sphere to sg
C ----------------------------------------------------------------
Ci Inputs
Ci   nds    :leading dimensions of sg
Ci   nda    :2*nda is the leading dimension of sgalph
Ci   iax    :neighbor table containing pair information (pairc.f)
Ci   npr    :number of neighbors in current cluster
Ci   sgalph :structure constant vector for a given cluster
Ci   nprs   :total number of strux accumulated before this one
Co Outputs
Co   sg     :sgalph is appended to sg
Cb Bugs
Cb   No check is made to ensure that sg is sufficiently large.
Cr Remarks
Cr   Converts sg to format sg(L',L,i,j,R'), the indices corresponding
Cr   to sg_iR'L',jRL (index R is dropped since sg is only stored for 
Cr   the head node of each cluster). 
Cr
Cr   The meaning of indecies is as follows. Given two sets of functions
Cr   {F_iRL,i=0,1} associated with each site R and channel L, one 
Cr   converts these to the "value-Laplacian" set of functions U_jRL as
Cr     U_jRL = \sum iR'L' sg(L',L,i,j,R') * F_iR'L'  
Cr
Cr   Index j = 0,1 refers to the "value" and "Laplacian" parts of the
Cr   screened basis U_jRL (ie the value of U_0RL is zero unless it is 
Cr   channel L at sphere R, same for Laplacian of U_1RL).
Cr
Cu Updates
Cu   17 Jan 08 (S. Lozovoi) Adapted from addtos.f
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nds,nda,npr,nprs
      integer niax,nkap
      parameter (niax=10, nkap=2)
      integer iax(niax,npr)
      double precision sgalph(2*nda,2*nda),
     .  sg(nds,nds,nkap,nkap,nprs+npr)
C Local parameters
      integer ipr,ilm,klm,ii,nlma,nlmb
      integer ikap,jkap

      nlmb = iax(9,1)
      ii = 0
      do  ipr = 1, npr
        nlma = iax(9,ipr)
        do  ikap = 1, nkap
        do  jkap = 1, nkap
          do  klm = 1, nlmb
          do  ilm = 1, nlma
            sg(ilm,klm,ikap,jkap,ipr+nprs) = 
     .        sgalph(ii+ilm+nlma*(ikap-1),klm+nlmb*(jkap-1))
          enddo
          enddo
        enddo
        enddo
        ii = ii+2*nlma
      enddo

      end

