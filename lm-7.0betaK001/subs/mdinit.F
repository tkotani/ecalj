      subroutine mdinit
C- Open the MD-files and print the corresponding header
C--------------------------------------------------------------------      
C     implicit none
      integer fopn
      
      call awrit0('#Time [fs] Temp. [K] [-----'//
     x     'Total Momentum x,y,z-----]     zeta       <zeta>'//
     x     '      <zeta>^2',' ',300,fopn('MD'))
      call awrit0('#Time [fs]   Hcons       epot+ekin        etot'//
     x     '          ekin     [--------Total force x,y,z-------]',
     x     ' ',300,fopn('ENG'))
      call awrit0('#Time [fs]  Q_0      [----------- Q_1 -----------'//
     x     '-]   [----------------------- Q_2 ----------------------]',
     x     ' ',300,fopn('MPOL1'))
      call awrit0('#Time [fs]  Q_0      [----------- Q_1 -----------'//
     x     '-]   [----------------------- Q_2 ----------------------]',
     x     ' ',300,fopn('MPOL2'))

      end

