Cgetarg...info...           structure ['writebasis', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['writebasis', 'm_struc_def']
      subroutine writebasis(nspec,sspec,ifi)

       
       use m_struc_def  !Cgetarg

c takao wrote this based on Mark's ioorbp below. 
c This is to write a basis set file mtopara.* at the end of lmfa.
c Default setting of basis set is embedded in this.
      implicit none
      integer mode,nkapi,nspec,ifi,nkmax
Cgetarg       double precision sspec(1)
       type(s_spec)::sspec(*)

C ... Local parameters
      logical scat,rdstrn
      character spid*8,sout*256,szz*256,sss*256
      integer jfi,n0,parg,nkap0
      parameter (n0=10,nkap0=3)
      integer is,ks,lp1,jp,ikap,lmxb,lmxi,j1,j2,ls,i,j,lpz,
     .  it(n0),idxdn(n0,nkap0)
      double precision rtab(n0,nkap0),etab(n0,nkap0),orbp(n0,2,nkap0),
     .  dasum, zz
      integer:: i_spacks,i_copy_size
        jfi = ifi
          do  110  is = 1, nspec
Cgetarg           call spacks(0,'spec name',sspec,spid,is,is)
           do i_spacks=is,is 
           call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
           enddo

Cgetarg...info...           real(8) :: orbp(n0,2,nkap0) , real(8) :: sspec%orbp 60
Cgetarg...info...           integer :: lmxb , integer(8) :: sspec%lmxb 1
Cgetarg...info...           real(8) :: zz , real(8) :: sspec%z 1
Cgetarg           call upack('spec orbp lmxb z',sspec,is, orbp,lmxb,zz,0)
           
           i_copy_size=size(sspec(is)%orbp) 
           call dcopy(i_copy_size,sspec(is)%orbp,1,orbp,1) 
           lmxb=sspec(is)%lmxb
           zz=sspec(is)%z

          etab=0d0
          rtab=0d0
          sout = spid
          ikap = 1
            call dcopy(n0,orbp(1,1,ikap),1,rtab(1,ikap),1)
            call dcopy(n0,orbp(1,2,ikap),1,etab(1,ikap),1)
            if (dasum(n0,orbp(1,1,ikap),1) .ne. 0) then
            write(jfi,"(a,$)")trim(sout)//'@ RSMH='
            do i=1,lmxb+1
              if(rtab(i,ikap)>=0d0) then
                write(jfi,"(' ',f5.3,$)") rtab(i,ikap)
              endif
            enddo
            write(jfi,"(a,$)")'  EH='
            do i=1,lmxb+1
              if(etab(i,ikap)<0d0) then
                write(jfi,"(' ',f6.3,$)") etab(i,ikap)
              endif
            enddo
            endif
c            write(jfi,"(' KMXA=',$)")

c Fractional part of P for empty sphere
c l      pfree   
c 0     .500000  
c 1     .250000  
c 2     .147584  
c 3     .102416  
c 4     .077979  
c 5     .062833  
            szz=' '
            if(1.5d0<zz .and.  zz<3.5d0) then
              write(sss,"(f8.3)") zz-2d0   
              szz='  #PZ=11.9 #P=2.6 #Q=0,'//trim(adjustl(sss))
            endif
c p channel local orbital
            if(10.5d0<zz .and. zz<12.5d0) then !2p
              szz='  PZ=0,12.9 P=0,3.3'
            endif
            if(18.5d0<zz .and. zz<20.5d0) then !3p
              szz='  PZ=0,13.9 P=0,4.3'
            endif
            if(36.5d0<zz .and. zz<38.5d0) then !4p lo
              szz='  PZ=0,14.9 P=0,5.3'
            endif
            if(54.5d0<zz .and. zz<56.5d0) then !5p lo
              szz='  PZ=0,15.9 P=0,6.3'
            endif
            if(86.5d0<zz .and. zz<88.5d0) then !6p lo
              szz='  PZ=0,16.9 P=0,7.3'
            endif
c d channel local orbital
            if(29.5d0<zz .and. zz<32.5d0) then !3d
              szz='  PZ=0,0,13.9 P=0,0,4.2'
            endif
            if(47.5d0<zz .and. zz<51.5d0) then !4d
              szz='  PZ=0,0,14.9 P=0,0,5.2'
            endif
            if(79.5d0<zz .and. zz<83.5d0) then !4d
              szz='  PZ=0,0,15.9 P=0,0,6.2'
            endif
c f channel
            if(70.5d0<zz .and. zz<73.5d0) then !4f
              szz='  PZ=0,0,0,14.9 P=0,0,0,5.15'
            endif
            write(jfi,"(a)") trim(szz)
  110     continue
Cgetarg       end
       end subroutine writebasis 

