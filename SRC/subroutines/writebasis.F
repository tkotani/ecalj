c$$$      subroutine writebasis(nspec,sspec)!,ifi)
c$$$      use m_struc_def,only: s_spec
c$$$      use m_ext,only: sname
c$$$c takao wrote this based on Mark's ioorbp below.
c$$$c This is to write a basis set file mtopara.* at the end of lmfa.
c$$$c Default setting of basis set is embedded in this.
c$$$      implicit none
c$$$      integer mode,nkapi,nspec,ifi,nkmax
c$$$      type(s_spec)::sspec(*)
c$$$
c$$$C ... Local parameters
c$$$      logical scat,rdstrn
c$$$      character spid*8,sout*256,szz*256,sss*256
c$$$      integer jfi,n0,parg,nkap0
c$$$      parameter (n0=10,nkap0=3)
c$$$      integer is,ks,lp1,jp,ikap,lmxb,lmxi,j1,j2,ls,i,j,lpz,
c$$$     .it(n0),idxdn(n0,nkap0)
c$$$      double precision rtab(n0,nkap0),etab(n0,nkap0),orbp(n0,2,nkap0),
c$$$     .dasum, zz
c$$$      integer:: i_spacks,i_copy_size
c$$$      open(newunit=jfi,file='mtopara.'//trim(sname))
c$$$      do  110  is = 1, nspec
c$$$        spid = sspec(is)%name
c$$$
c$$$        i_copy_size=size(sspec(is)%orbp)
c$$$        call dcopy(i_copy_size,sspec(is)%orbp,1,orbp,1)
c$$$        lmxb=sspec(is)%lmxb
c$$$        zz=sspec(is)%z
c$$$
c$$$        etab=0d0
c$$$        rtab=0d0
c$$$        sout = spid
c$$$        ikap = 1
c$$$        call dcopy(n0,orbp(1,1,ikap),1,rtab(1,ikap),1)
c$$$        call dcopy(n0,orbp(1,2,ikap),1,etab(1,ikap),1)
c$$$        if (dasum(n0,orbp(1,1,ikap),1) .ne. 0) then
c$$$          write(jfi,"(a,$)")trim(sout)//'@ RSMH='
c$$$          do i=1,lmxb+1
c$$$            if(rtab(i,ikap)>=0d0) then
c$$$              write(jfi,"(' ',f5.3,$)") rtab(i,ikap)
c$$$            endif
c$$$          enddo
c$$$          write(jfi,"(a,$)")'  EH='
c$$$          do i=1,lmxb+1
c$$$            if(etab(i,ikap)<0d0) then
c$$$              write(jfi,"(' ',f6.3,$)") etab(i,ikap)
c$$$            endif
c$$$          enddo
c$$$        endif
c$$$c            write(jfi,"(' KMXA=',$)")
c$$$
c$$$c Fractional part of P for empty sphere
c$$$c l      pfree
c$$$c 0     .500000
c$$$c 1     .250000
c$$$c 2     .147584
c$$$c 3     .102416
c$$$c 4     .077979
c$$$c 5     .062833
c$$$        szz=' '
c$$$        if(1.5d0<zz .and.  zz<3.5d0) then
c$$$c          write(sss,"(f8.3)") zz-2d0
c$$$          szz='  PZ=11.9 P=2.6' !//trim(adjustl(sss))
c$$$        endif
c$$$c p channel local orbital
c$$$        if(10.5d0<zz .and. zz<12.5d0) then !2p
c$$$          szz='  PZ=0,12.9 P=0,3.3'
c$$$        endif
c$$$        if(18.5d0<zz .and. zz<20.5d0) then !3p
c$$$          szz='  PZ=0,13.9 P=0,4.3'
c$$$        endif
c$$$        if(36.5d0<zz .and. zz<38.5d0) then !4p lo
c$$$          szz='  PZ=0,14.9 P=0,5.3'
c$$$        endif
c$$$        if(54.5d0<zz .and. zz<56.5d0) then !5p lo
c$$$          szz='  PZ=0,15.9 P=0,6.3'
c$$$        endif
c$$$        if(86.5d0<zz .and. zz<88.5d0) then !6p lo
c$$$          szz='  PZ=0,16.9 P=0,7.3'
c$$$        endif
c$$$c d channel local orbital
c$$$        if(29.5d0<zz .and. zz<32.5d0) then !3d
c$$$          szz='  PZ=0,0,13.9 P=0,0,4.2'
c$$$        endif
c$$$        if(47.5d0<zz .and. zz<51.5d0) then !4d
c$$$          szz='  PZ=0,0,14.9 P=0,0,5.2'
c$$$        endif
c$$$        if(79.5d0<zz .and. zz<83.5d0) then !5d
c$$$          szz='  PZ=0,0,15.9 P=0,0,6.2'
c$$$        endif
c$$$c f channel
c$$$        if(70.5d0<zz .and. zz<73.5d0) then !4f
c$$$          szz='  PZ=0,0,0,14.9 P=0,0,0,5.15'
c$$$        endif
c$$$        write(jfi,"(a)") trim(szz)
c$$$ 110  continue
c$$$      close(jfi)
c$$$      end subroutine writebasis
c$$$
