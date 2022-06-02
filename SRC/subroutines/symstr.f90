subroutine symstr(mode,nds,nsites,iax,nsp,nkap,sflg,s,sc,asym)
  !- Symmetrize structure constants
  ! ----------------------------------------------------------------
  !i Inputs
  !i   mode  :0 s is real
  !i         :1 s is complex
  !i         :2 s is complex; symmetrize s+_ij = s_ji
  !i   nds   :leading dimension of s
  !i   nsites:number of pairs
  !i   iax   :neighbor table containing pair information
  !i         :iax(*,i) contains information for pair i
  !i         :This routine uses:
  !i         :iax(6,i): index to conjugate (jb,ib) pair matching (ib,jb)
  !i         :iax(8,i): if nonzero, points to an pair which is to
  !i         :        : substitute for pair i.
  !i   nsp   :2 for spin-polarized case, otherwise 1
  ! o Inputs/Outputs
  ! o  sflg  :marks which pairs in s have been symmetrized.
  ! o        :On output, any pairs symmetrized by symstr are marked.
  ! o  s,sc  :real-space hamiltonian or structure-constant matrix.
  ! o        :Only s is used if mode=0 or sc is used if mode=1.
  ! o        :Any pair i which has a counterpart j is symmetrized.
  ! o        :*If neither pair i not pair j=iax(6,i) is symmetrized,
  ! o        : both j and i are symmetrized by averaging the two.
  ! o        :*If j=iax(6,i) is already symmetrized, i is copied from j.
  ! o        :*If i is already symmetrized, j is copied from i.
  !o Outputs
  !o   asym  :maximum asymmetry found when symmetrizing
  !r Remarks
  !u Updates
  !u   05 Aug 06 Passes parms for 2-kappa (just return for now)
  !u   30 Mar 03 complex case (mode>0) can symm. s_ij=s_ji or s+_ij=sji
  !u   17 Jul 02 Redesigned for more general cases (complex, nsp)
  !u             Altered argument list.
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: mode,nds,nsites,nsp,sflg(nsites),niax,nkap
  parameter (niax=10)
  integer :: iax(niax,nsites)
  double precision :: s(nds,nds,nsp,nkap,nkap,nsites),asym
  double complex sc(nds,nds,nsp,nkap,nkap,nsites)
  ! Local parameters
  integer :: i,j,lm1,lm2,isp,ii,jj,ik,jk
  double precision :: tmp
  double complex tmpc

  logical:: isanrg, l_dummy_isanrg

  ! Cludge for now ...
  asym = 0
  if (nkap /= 1) return
  ik = 1
  jk = 1

  ! ino isanrg is logical function,       call isanrg(mode,0,2,'symstr:','mode',.true.)
  l_dummy_isanrg=isanrg(mode,0,2,'symstr:','mode',.true.)

  do  10  i = 1, nsites
     j = iax(6,i)
     if (j == 0) goto 10

     ii = iax(8,i)
     if (ii == 0) ii = i
     jj = iax(8,j)
     if (jj == 0) jj = j

     !         print *, i,ii,' ',j,jj,' ',sflg(ii),sflg(jj)
     !     ... Symmetrize ii,jj as (ii+jj)/2 if neither was symmetrized
     if (sflg(ii) == 0 .AND. sflg(jj) == 0) then
        if (mode == 0) then
           do  302  isp = 1, nsp
              do  301  lm1 = 1, nds
                 do  30  lm2 = 1, nds
                    tmp = (s(lm1,lm2,isp,ik,jk,ii)+s(lm2,lm1,isp,ik,jk,jj))/2
                    asym = max(asym,abs(s(lm1,lm2,isp,ik,jk,ii)-tmp))
                    asym = max(asym,abs(s(lm2,lm1,isp,ik,jk,jj)-tmp))
                    !              if (asym .gt. 1d-3) then
                    !                print *, 'asym=',asym
                    !              endif
                    s(lm1,lm2,isp,ik,jk,ii) = tmp
                    s(lm2,lm1,isp,ik,jk,jj) = tmp
30               enddo
301           enddo
302        enddo
        else if (mode == 1) then
           do  312  isp = 1, nsp
              do  311  lm1 = 1, nds
                 do  31  lm2 = 1, nds
                    tmpc = (sc(lm1,lm2,isp,ik,jk,ii)+ &
                         sc(lm2,lm1,isp,ik,jk,jj))/2
                    asym = max(asym,abs(sc(lm1,lm2,isp,ik,jk,ii)-tmpc))
                    asym = max(asym,abs(sc(lm2,lm1,isp,ik,jk,jj)-tmpc))
                    !              if (asym .gt. .004) then
                    !                print *, i,lm1,lm2
                    !              endif
                    sc(lm1,lm2,isp,ik,jk,ii) = tmpc
                    sc(lm2,lm1,isp,ik,jk,jj) = tmpc
31               enddo
311           enddo
312        enddo
        else
           do  322  isp = 1, nsp
              do  321  lm1 = 1, nds
                 do  32  lm2 = 1, nds
                    tmpc = (sc(lm1,lm2,isp,ik,jk,ii)+ &
                         dconjg(sc(lm2,lm1,isp,ik,jk,jj)))/2
                    asym = max(asym,abs(sc(lm1,lm2,isp,ik,jk,ii)-tmpc))
                    asym =max(asym,abs(dconjg(sc(lm2,lm1,isp,ik,jk,jj))-tmpc))
                    sc(lm1,lm2,isp,ik,jk,ii) = tmpc
                    sc(lm2,lm1,isp,ik,jk,jj) = dconjg(tmpc)
32               enddo
321           enddo
322        enddo
        endif
        !       ... Flag this s(ii), s(jj) as symmetrized
        sflg(ii) = 1
        sflg(jj) = 1
        !     ... Symmetrize ii from jj if ii not symmetrized
     elseif (sflg(ii) == 0) then
        do  342  isp = 1, nsp
           do  341  lm1 = 1, nds
              do  34  lm2 = 1, nds
                 if (mode == 0) then
                    s(lm1,lm2,isp,ik,jk,ii) = s(lm2,lm1,isp,ik,jk,jj)
                 else if (mode == 1) then
                    sc(lm1,lm2,isp,ik,jk,ii) = sc(lm2,lm1,isp,ik,jk,jj)
                 else
                    sc(lm1,lm2,isp,ik,jk,ii) = &
                         dconjg(sc(lm2,lm1,isp,ik,jk,jj))
                 endif
34            enddo
341        enddo
342     enddo
        !       ... Flag this s(ii) as symmetrized
        sflg(ii) = 1
        !     ... Symmetrize jj from ii if jj not symmetrized
     elseif (sflg(jj) == 0) then
        do  352  isp = 1, nsp
           do  351  lm1 = 1, nds
              do  35  lm2 = 1, nds
                 if (mode == 0) then
                    s(lm2,lm1,isp,ik,jk,jj) = s(lm1,lm2,isp,ik,jk,ii)
                 else if (mode == 1) then
                    sc(lm2,lm1,isp,ik,jk,jj) = sc(lm1,lm2,isp,ik,jk,ii)
                 else
                    sc(lm2,lm1,isp,ik,jk,jj) = &
                         dconjg(sc(lm1,lm2,isp,ik,jk,ii))
                 endif
35            enddo
351        enddo
352     enddo
        !       ... Flag this s(jj) as symmetrized
        sflg(jj) = 1
     endif
     goto 10
10 enddo

  !     call info(60,0,0,' symstr: max asymmetry=%d',asym,0)
end subroutine symstr

