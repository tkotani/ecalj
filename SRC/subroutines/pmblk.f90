subroutine pmblk(nblk,ipm,ioffo,ioffn,lds,SRC,mode,alfa,nlma,ldd, &
     DST,nmax)
  !- Clips and optionally permutes subblocks of a double precision array
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nblk  :number of subblocks to permute
  !i   ipm   :permutation of subblocks; see Remarks
  !i   ioffo :ioffo(1:nblk+1) markers to the start of subblocks in SRC:
  !i          Subblock i starts at ioffo(i)+1 and ends at ioffo(i+1)
  !i          Thus SRC subblock has has size ioffo(i+1)-ioffo(i)
  !i   ioffn :ioffn(1:nblk+1) markers to the start of subblocks in DST:
  !i          Subblock i starts at ioffn(i)+1 and ends at ioffn(i+1)
  !i          Thus DST subblock has has size ioffn(i+1)-ioffn(i)
  !i          NB: ioffn is the same as ioffo, but for DST.
  !i   lds   :leading dimension of SRC
  !i   SRC   :source matrix
  !i   mode:  one's digit: 0 SRC offsets to block i computed by ioffo(i)
  !i                     : 1 SRC offsets are relative:  offset to block i
  !i                     :       computed by ioffo(i)-ioffo(1)
  !i                     : 2 SRC is not permuted at all: columns read
  !i                     :       sequentially starting from col 1
  !i                     :       In this case, ioffo is not used.
  !i                     : 3     same as 2, but cols are read starting
  !i                     :       at col ioffo(1).
  !i                     :       In this case, only ioffo(1) is used.
  !i                     : NB: only 1 or 10s digit may be set to 2 or 3
  !i          ten's digit: 0 DST offsets to block i computed by ioffn(i)
  !i                     : 1 DST offsets are relative:  offset to block i
  !i                     :       computed by ioffn(i)-ioffn(1)
  !i                     : 2 DST is not permuted at all: columns written
  !i                     :       sequentially starting from col ioffn(1).
  !i                     :       In this case, ioffn is not used.
  !i                     : 3     same as 2, but cols are written starting
  !i                     :       at col ioffn(1)
  !i                     :       In this case, only ioffn(1) is used.
  !i                     : NB: only 1 or 10s digit may be set to 2 or 3
  !i          100's digit: 0 no permutation of subblock order
  !i                         In this case, ipm is not used
  !i                       1 permute sequence of subblocks in SRC by ipm
  !i                         i.e. SRC subblock i endpoints are computed
  !i                         from ioffo(j) and ioffo(j+1), j=ipm(i)
  !i                         NB: not compatible with 1s digit >=2
  !i                       2 permute sequence of subblocks in DST by ipm
  !i                         i.e. DST subblock i endpoints are computed
  !i                         from ioffo(j) and ioffn(j+1), j=ipm(i)
  !i                         NB: not compatible with 10s digit >=2
  !i         1000's digit: 1 add into new array, not copy
  !i                       2 copy (add) alfa(j)*SRC(*,j) into DST(*,j')
  !i                       3 both 1 and 2
  !i                       4 add 4 to exchange roles of rows and columns,
  !i                         i.e. permute rows for each nlma columns
  !i                       8 permute subblocks in both rows and cols
  !i                        (see Remarks).
  !i                       9 like 8, but add into DST.
  !i   alfa  :copy or add alfa(j)*SRC(*,j) into DST(*,j'), if mode set
  !i   nlma  :number of rows for which to permute columns
  !i          (not used if rows are also permuted; see 1000s digit mode)
  !i   ldd   :leading dimension of DST
  !o Outputs
  !o   DST   :destination matrix
  !o   nmax  :largest column index in DST affected by call to pmblk.
  !r Remarks
  !r  Columns are grouped into a sequence of subblocks.  pmblk copies
  !r  entire subblocks of SRC: columns within a subblock retain their
  !r  order when copied.  The subblock order, however, may be permuted,
  !r  as specified by the 10s digit of mode.
  !r
  !r  * The permutations may occur by permuting the order in which the
  !r    subblocks are copied from SRC (if 1000s digit mode is 1) or by
  !r    permuting the order in which subblocks are copied to DST (if 1000s
  !r    digit mode is 2).  By default subblocks are copied in their
  !r    natural order; thus there is no permutation of subblocks if 1000s
  !r    digit mode is 0.
  !r
  !r  * If both SRC and DST are being permuted, the size of the ith
  !r    subblock copied is the smaller of sizes of the subblock i in SRC
  !r    and subblock i in DST.  Thus, if the size of DST subblock i is
  !r    smaller than SRC, only the first portion of the SRC subblock is
  !r    copied; if the size DST is larger than SRC, only the first columns
  !r    in the DST subblock will be affected.
  !r
  !r  * By default pmblk permutes columns, doing the permutation for each of
  !r    nlma rows.  If mode 1000's digit >=8, pmblk permutes both columns
  !r    and rows.
  !u Updates
  !u   16 May 03 Added options 2 and 3 to 1 and 10s digit mode
  ! ----------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: nblk,lds,ldd,nlma,mode,ioffn(nblk),ioffo(nblk), &
       ipm(nblk),nmax
  double precision :: dst(ldd,1),src(lds,1),alfa(1)
  ! ... Local parameters
  integer :: ilma,ib,iofo,ilm,iofn,nlmi,ipn,ipo,offo0,offn0, &
       nn,no,mode1,mode2,mode3,mode4,jb,jofo,jlm,jofn,nlmj,jpn,jpo
  logical:: l_dummy_isanrg, isanrg

  mode1 = mod(mode,10)
  mode2 = mod(mode/10,10)
  mode3 = mod(mode/100,10)
  mode4 = mod(mode/1000,10)

  !     if (mode3.ne.0) call awrit2('pm ipm  %n,3i',' ',80,6,nblk,ipm)
  !     call awrit2('pm iofo %n,3i',' ',80,6,nblk+1,ioffo)
  !     call awrit2('pm iofn %n,3i',' ',80,6,nblk+1,ioffn)

  offn0 = 0
  offo0 = 0
  if (mod(mode1,2) == 1) offo0 = ioffo(1)
  if (mod(mode2,2) == 1) offn0 = ioffn(1)
  if (mode1 >= 2 .AND. mode2 >= 2) &
       call rx('pmblk: must permute one of SRC or DST')
  if (mode1 >= 2 .AND. mode3 == 1) &
       call rxi('pmblk: imcompatible switches, mode=',mode)
  if (mode2 >= 2 .AND. mode3 == 2) &
       call rxi('pmblk: imcompatible switches, mode=',mode)
  ! ino isanrg is logical function,       call isanrg(mode3,0,2,'pmblk:','100s digit mode', .true.)
  l_dummy_isanrg=isanrg(mode3,0,2,'pmblk:','100s digit mode', .true.)

  nmax = -1
  nlmi = 0
  iofo = offo0
  iofn = offn0
  if (mode4 <= 7) then
     do  10  ib = 1, nblk
        ipo = ib
        if (mode3 == 1) ipo = ipm(ib)
        ipn = ib
        if (mode3 == 2) ipn = ipm(ib)
        !         sizes and offsets to SRC and DST blocks
        if (mode1 >= 2) then
           iofo = iofo + nlmi
           iofn = ioffn(ipn) - offn0
           nlmi = ioffn(ipn+1)-ioffn(ipn)
        elseif (mode2 >= 2) then
           iofn = iofn + nlmi
           iofo = ioffo(ipo) - offo0
           nlmi = ioffo(ipo+1)-ioffo(ipo)
        else
           iofn = ioffn(ipn) - offn0
           iofo = ioffo(ipo) - offo0
           nn = ioffn(ipn+1)-ioffn(ipn)
           no = ioffo(ipo+1)-ioffo(ipo)
           nlmi = min(nn,no)
        endif

        if (mode4 == 0) then
           do  20  ilm = 1, nlmi
              do  20  ilma = 1, nlma
                 dst(ilma,ilm+iofn) = src(ilma,ilm+iofo)
20            enddo
           elseif (mode4 == 1) then
              do  25  ilm = 1, nlmi
                 do  25  ilma = 1, nlma
                    dst(ilma,ilm+iofn) = dst(ilma,ilm+iofn) &
                         + src(ilma,ilm+iofo)
25               enddo
              elseif (mode4 == 2) then
                 do  30  ilm = 1, nlmi
                    do  30  ilma = 1, nlma
                       dst(ilma,ilm+iofn) = src(ilma,ilm+iofo)*alfa(ilm+iofn)
30                  enddo
                 elseif (mode4 == 3) then
                    do  35  ilm = 1, nlmi
                       do  35  ilma = 1, nlma
                          dst(ilma,ilm+iofn) = dst(ilma,ilm+iofn) &
                               + src(ilma,ilm+iofo)*alfa(ilm+iofn)
35                     enddo
                    elseif (mode4 == 4) then
                       do  40  ilma = 1, nlma
                          do  40  ilm = 1, nlmi
                             dst(ilm+iofn,ilma) = src(ilm+iofo,ilma)
40                        enddo
                       elseif (mode4 == 5) then
                          do  45  ilma = 1, nlma
                             do  45  ilm = 1, nlmi
                                dst(ilm+iofn,ilma) = dst(ilm+iofn,ilma) &
                                     + src(ilm+iofo,ilma)
45                           enddo
                          elseif (mode4 == 6) then
                             do  50  ilma = 1, nlma
                                do  50  ilm = 1, nlmi
                                   dst(ilm+iofn,ilma) = src(ilm+iofo,ilma)*alfa(ilm+iofn)
50                              enddo
                             else
                                do  55  ilma = 1, nlma
                                   do  55  ilm = 1, nlmi
                                      dst(ilm+iofn,ilma) = dst(ilm+iofn,ilma) &
                                           + src(ilm+iofo,ilma)*alfa(ilm+iofn)
55                                 enddo
                                endif
                                nmax = max(nn+iofn,nmax)
10                           enddo
                          else
                             do  100  ib = 1, nblk
                                ipo = ib
                                if (mode3 == 1) ipo = ipm(ib)
                                ipn = ib
                                if (mode3 == 2) ipn = ipm(ib)
                                !         sizes and offsets to SRC and DST blocks, col index
                                if (mode1 >= 2) then
                                   iofo = iofo + nlmi
                                   iofn = ioffn(ipn) - offn0
                                   nlmi = ioffn(ipn+1)-ioffn(ipn)
                                elseif (mode2 >= 2) then
                                   iofn = iofn + nlmi
                                   iofo = ioffo(ipo) - offo0
                                   nlmi = ioffo(ipo+1)-ioffo(ipo)
                                else
                                   iofn = ioffn(ipn) - offn0
                                   iofo = ioffo(ipo) - offo0
                                   nn = ioffn(ipn+1)-ioffn(ipn)
                                   no = ioffo(ipo+1)-ioffo(ipo)
                                   nlmi = min(nn,no)
                                endif
                                nlmj = 0
                                jofo = offo0
                                jofn = offn0
                                do  110  jb = 1, nblk
                                   jpo = jb
                                   if (mode3 == 1) jpo = ipm(jb)
                                   jpn = jb
                                   if (mode3 == 2) jpn = ipm(jb)
                                   !           sizes and offsets to SRC and DST blocks, row index
                                   if (mode1 >= 2) then
                                      jofo = jofo + nlmj
                                      jofn = ioffn(jpn) - offn0
                                      nlmj = ioffn(jpn+1)-ioffn(jpn)
                                   elseif (mode2 >= 2) then
                                      jofn = jofn + nlmj
                                      jofo = ioffo(jpo) - offo0
                                      nlmj = ioffo(jpo+1)-ioffo(jpo)
                                   else
                                      jofn = ioffn(jpn) - offn0
                                      jofo = ioffo(jpo) - offo0
                                      nn = ioffn(jpn+1)-ioffn(jpn)
                                      no = ioffo(jpo+1)-ioffo(jpo)
                                      nlmj = min(nn,no)
                                   endif
                                   if (mode4 == 8) then
                                      do  120  jlm = 1, nlmj
                                         do  120  ilm = 1, nlmi
                                            dst(ilm+iofn,jlm+jofn) = src(ilm+iofo,jlm+jofo)
120                                      enddo
                                      else
                                         do  125  jlm = 1, nlmj
                                            do  125  ilm = 1, nlmi
                                               dst(ilm+iofn,jlm+jofn) = dst(ilm+iofn,jlm+jofn) + &
                                                    src(ilm+iofo,jlm+jofo)
125                                         enddo
                                         endif
110                                   enddo
                                      nmax = max(nn+iofn,nmax)
100                                enddo
                                endif
                              end subroutine pmblk
! #if TEST
!                               subroutine fmain
!                                 implicit none
!                                 integer :: nblk,ioffs(5),ioffd(5),iprm(4)
!                                 double precision :: s(10),d(10)
!                                 integer :: nmax,i,mode

!                                 ioffs(1) = 0
!                                 ioffs(2) = 1
!                                 ioffs(3) = 3
!                                 ioffs(4) = 6
!                                 ioffs(5) = 10

!                                 ioffd(1) = 0
!                                 ioffd(2) = 2
!                                 ioffd(3) = 6
!                                 ioffd(4) = 9
!                                 ioffd(5) = 10

!                                 nblk = 4

!                                 iprm(1) = 2
!                                 iprm(2) = 4
!                                 iprm(3) = 3
!                                 iprm(4) = 1

!                                 do  i = 1, 10
!                                    s(i) = 1.1d0*i
!                                 enddo

!                                 print 332, ioffs
! 332                             format(10i3)
!                                 print 332, ioffd
!                                 print 332, iprm

!                                 !     Do the permutation
!                                 !     NB: 4000 should have no effect; just tests different branch.
!                                 call dpzero(d,10)
!                                 mode = 4100
!                                 call pmblk(nblk,iprm,ioffs,ioffd,1,s,mode,0d0,1,1,d,nmax)

!                                 print 333, (s(i), i=1,10)
!                                 print 333, (d(i), i=1,10)
! 333                             format(10f8.1)

!                                 !     Undo the permutation
!                                 mode = 0200
!                                 call dpzero(s,10)
!                                 call pmblk(nblk,iprm,ioffd,ioffs,1,d,mode,0d0,1,1,s,nmax)
!                                 print 333, (s(i), i=1,10)

!                                 !     Same permutation, but use mode for sequential DST order
!                                 call dpzero(d,10)
!                                 mode = 4120
!                                 call pmblk(nblk,iprm,ioffs,ioffd,1,s,mode,0d0,1,1,d,nmax)
!                                 print 333, (d(i), i=1,10)

!                               end subroutine fmain
! #endif

