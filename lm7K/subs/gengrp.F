      subroutine gengrp(csym,dorig,gens,instr,iopt,iorig,isym,ngen,
     .                  plat,qlat)
C- Supplies information about the space-group
C ----------------------------------------------------------------------
Ci Inputs:
Ci   instr :string containing space-group number of symbol
Ci         :used only for iopt=-1
Ci   iopt  := -1 gengrp returns space group info: isym,csym,gens
Ci          =  1 gengrp to supply space group number from group info
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci         :used only for iopt=1
Ci   qlat  :primitive translation vectors in reciprocal space
Ci         :used only for iopt=1
Cio Inputs/Outputs:
Cio  gens  :symbol for generators
Cio  iorig :origin choice
Cio  isym  :numbers characterizing the symmetry of lattice and crystal
Co         :data is output if iopt=-1: else some data is input
Co          isym(7)=space-group number
Cio         isym(6)=point-group number
Co          isym(5)=crystal system
Co          isym(4)=crystal family
Cio         isym(3)=centring type
Cio         isym(2)=Bravais-lattice number
Cio         isym(1)=Bravais-system number
Cio  ngen  :number of symmetry generators
Co Outputs:
Co   csym  :names corresponding to integer indices isym
Co   dorig :new origin (iopt=1)
Cr Remarks:
Cu Updates
Cu   02 Nov 01 Adapted from Stuttgart lmto56 spcgrp.f
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer iopt,iorig,isym(7),ngen
      double precision dorig(3),plat(9,*),qlat(9)
      character*10 instr
      character*30 csym(7)
      character*72 gens
C Local variables:
      integer idum,nelts,isym1(0:230),isym2(7,7),isym3,isym4(0:7),
     .  isym5(0:32),isym6(0:230),isym7,isym8,iprint,lgunit,stdo
      logical a2bin
      parameter(nelts=499)
      character*1  csym3(0:7),corig
      character*10 csym7
      character*53 symgrp
      character*68 c(nelts)
C     integer info
C     character*72 gens1(4)
C ... External calls
      external chcase,gtcsym,isanrg,query,rx,rxs,words
C Data statements:
      data isym1/0,2*1,13*2,59*3,68*4,3*6,1*5,1*6,1*5,6*6,1*5,4*6,2*5,
     .           4*6,2*5,27*6,36*7/
      data isym2/1,0,0,0, 0, 0,0,    2,3,3,3,3,0, 0,    4,5,5,5,6,7,0,
     .           8,0,0,0, 9, 0,0,   10,0,0,0,0,0,10,   11,0,0,0,0,0,0,
     .          12,0,0,0,13,14,0/
      data csym3/' ','P','C','B','A','I','F','R'/
      data isym4/0,1,2,3,4,2*5,6/
      data isym5/0,2*1,3*2,3*3,7*4,5*5,7*6,5*7/
      data isym6/0,1*1,1*2,3*3,4*4,6*5,9*6,22*7,28*8,6*9,2*10,6*11,
     .          10*12,12*13,12*14,20*15,4*16,2*17,7*18,6*19,6*20,6*21,
     .          1*22,2*23,6*24,4*25,4*26,4*27,5*28,7*29,8*30,6*31,10*32/

      c(  1)='  1 P1         E'
      c(  2)='  2 P-1        I'
      c(  3)='  3 P2         R2Y'
      c(  4)='  3 P121       R2Y'
      c(  5)='  4 P2_1       R2Y::(0,1/2,0)'
      c(  6)='  4 P12_11     R2Y::(0,1/2,0)'
      c(  7)='  5 C2         R2Y'
      c(  8)='  5 A2         R2Y'
      c(  9)='  5 I2         R2Y'
      c( 10)='  5 C121       R2Y'
      c( 11)='  5 A121       R2Y'
      c( 12)='  5 I121       R2Y'
      c( 13)='  6 Pm         MY'
      c( 14)='  6 P1m1       MY'
      c( 15)='  7 Pc         MY::(0,0,1/2)'
      c( 16)='  7 Pa         MY::(1/2,0,0)'
      c( 17)='  7 Pn         MY::(1/2,0,1/2)'
      c( 18)='  7 P1c1       MY::(0,0,1/2)'
      c( 19)='  7 P1a1       MY::(1/2,0,0)'
      c( 20)='  7 P1n1       MY::(1/2,0,1/2)'
      c( 21)='  8 Cm         MY'
      c( 22)='  8 Am         MY'
      c( 23)='  8 Im         MY'
      c( 24)='  8 C1m1       MY'
      c( 25)='  8 A1m1       MY'
      c( 26)='  8 I1m1       MY'
      c( 27)='  9 Cc         MY::(0,0,1/2)'
      c( 28)='  9 Cn         MY::(1/2,0,1/2)'
      c( 29)='  9 Aa         MY::(1/2,0,0)'
      c( 30)='  9 An         MY::(1/2,1/2,0)'
      c( 31)='  9 Ia         MY::(1/2,0,0)'
      c( 32)='  9 Ic         MY::(0,0,1/2)'
      c( 33)='  9 C1c1       MY::(0,0,1/2)'
      c( 34)='  9 C1n1       MY::(1/2,0,1/2)'
      c( 35)='  9 A1a1       MY::(1/2,0,0)'
      c( 36)='  9 A1n1       MY::(1/2,1/2,0)'
      c( 37)='  9 I1a1       MY::(1/2,0,0)'
      c( 38)='  9 I1c1       MY::(0,0,1/2)'
      c( 39)=' 10 P2/m       I MY'
      c( 40)=' 10 P12/m1     I MY'
      c( 41)=' 11 P2_1/m     I MY::(0,1/2,0)'
      c( 42)=' 11 P12_1/m1   I MY::(0,1/2,0)'
      c( 43)=' 12 C2/m       I MY'
      c( 44)=' 12 A2/m       I MY'
      c( 45)=' 12 I2/m       I MY'
      c( 46)=' 12 C12/m1     I MY'
      c( 47)=' 12 A12/m1     I MY'
      c( 48)=' 12 I12/m1     I MY'
      c( 49)=' 13 P2/c       I MY::(0,0,1/2)'
      c( 50)=' 13 P2/a       I MY::(1/2,0,0)'
      c( 51)=' 13 P2/n       I MY::(1/2,0,1/2)'
      c( 52)=' 13 P12/c1     I MY::(0,0,1/2)'
      c( 53)=' 13 P12/a1     I MY::(1/2,0,0)'
      c( 54)=' 13 P12/n1     I MY::(1/2,0,1/2)'
      c( 55)=' 14 P2_1/c     I MY::(0,1/2,1/2)'
      c( 56)=' 14 P2_1/a     I MY::(1/2,1/2,0)'
      c( 57)=' 14 P2_1/n     I MY::(1/2,1/2,1/2)'
      c( 58)=' 14 P12_1/c1   I MY::(0,1/2,1/2)'
      c( 59)=' 14 P12_1/a1   I MY::(1/2,1/2,0)'
      c( 60)=' 14 P12_1/n1   I MY::(1/2,1/2,1/2)'
      c( 61)=' 15 C2/c       I MY::(0,0,1/2)'
      c( 62)=' 15 C2/n       I MY::(1/2,0,1/2)'
      c( 63)=' 15 A2/a       I MY::(1/2,0,0)'
      c( 64)=' 15 A2/n       I MY::(1/2,1/2,0)'
      c( 65)=' 15 I2/a       I MY::(1/2,0,0)'
      c( 66)=' 15 I2/c       I MY::(0,0,1/2)'
      c( 67)=' 15 C12/c1     I MY::(0,0,1/2)'
      c( 68)=' 15 C12/n1     I MY::(1/2,0,1/2)'
      c( 69)=' 15 A12/a1     I MY::(1/2,0,0)'
      c( 70)=' 15 A12/n1     I MY::(1/2,1/2,0)'
      c( 71)=' 15 I12/a1     I MY::(1/2,0,0)'
      c( 72)=' 15 I12/c1     I MY::(0,0,1/2)'
      c( 73)=' 16 P222       R2X R2Y'
      c( 74)=' 17 P222_1     R2X R2Y::(0,0,1/2)'
      c( 75)=' 17 P22_12     R2X R2Y::(0,1/2,0)'
      c( 76)=' 17 P2_122     R2X::(1/2,0,0) R2Y'
      c( 77)=' 18 P22_12_1   R2X R2Y::(0,1/2,1/2)'
      c( 78)=' 18 P2_122_1   R2X::(1/2,0,1/2) R2Y'
      c( 79)=' 18 P2_12_12   R2X::(1/2,1/2,0) R2Y::(1/2,1/2,0)'
      c( 80)=' 19 P2_12_12_1 R2X::(1/2,1/2,0) R2Y::(0,1/2,1/2)'
      c( 81)=' 20 C222_1     R2X R2Y::(0,0,1/2)'
      c( 82)=' 20 A2_122     R2X::(1/2,0,0) R2Y'
      c( 83)=' 20 B22_12     R2X R2Y::(0,1/2,0)'
      c( 84)=' 21 C222       R2X R2Y'
      c( 85)=' 21 A222       R2X R2Y'
      c( 86)=' 21 B222       R2X R2Y'
      c( 87)=' 22 F222       R2X R2Y'
      c( 88)=' 23 I222       R2X R2Y'
      c( 89)=' 24 I2_12_12_1 R2X::(0,0,1/2) R2Y::(1/2,0,0)'
      c( 90)=' 25 Pmm2       MX MY'
      c( 91)=' 25 Pm2m       MX MZ'
      c( 92)=' 25 P2mm       MY MZ'
      c( 93)=' 26 Pmc2_1     MX MY::(0,0,1/2)'
      c( 94)=' 26 Pcm2_1     MX::(0,0,1/2) MY'
      c( 95)=' 26 Pm2_1b     MX MZ::(0,1/2,0)'
      c( 96)=' 26 Pb2_1m     MX::(0,1/2,0) MZ'
      c( 97)=' 26 P2_1ma     MY MZ::(1/2,0,0)'
      c( 98)=' 26 P2_1am     MY::(1/2,0,0) MZ'
      c( 99)=' 27 Pcc2       MX::(0,0,1/2) MY::(0,0,1/2)'
      c(100)=' 27 Pb2b       MX::(0,1/2,0) MZ::(0,1/2,0)'
      c(101)=' 27 P2aa       MY::(1/2,0,0) MZ::(1/2,0,0)'
      c(102)=' 28 Pma2       MX::(1/2,0,0) MY::(1/2,0,0)'
      c(103)=' 28 Pbm2       MX::(0,1/2,0) MY::(0,1/2,0)'
      c(104)=' 28 Pm2a       MX::(1/2,0,0) MZ::(1/2,0,0)'
      c(105)=' 28 Pc2m       MX::(0,0,1/2) MZ::(0,0,1/2)'
      c(106)=' 28 P2mb       MY::(0,1/2,0) MZ::(0,1/2,0)'
      c(107)=' 28 P2cm       MY::(0,0,1/2) MZ::(0,0,1/2)'
      c(108)=' 29 Pca2_1     MX::(1/2,0,1/2) MY::(1/2,0,0)'
      c(109)=' 29 Pbc2_1     MX::(0,1/2,0) MY::(0,1/2,1/2)'
      c(110)=' 29 Pb2_1a     MX::(1/2,1/2,0) MZ::(1/2,0,0)'
      c(111)=' 29 Pc2_1b     MX::(0,0,1/2) MZ::(0,1/2,1/2)'
      c(112)=' 29 P2_1ab     MY::(1/2,1/2,0) MZ::(0,1/2,0)'
      c(113)=' 29 P2_1ca     MY::(0,0,1/2) MZ::(1/2,0,1/2)'
      c(114)=' 30 Pnc2       MX::(0,1/2,1/2) MY::(0,1/2,1/2)'
      c(115)=' 30 Pcn2       MX::(1/2,0,1/2) MY::(1/2,0,1/2)'
      c(116)=' 30 Pn2b       MX::(0,1/2,1/2) MZ::(0,1/2,1/2)'
      c(117)=' 30 Pb2n       MX::(1/2,1/2,0) MZ::(1/2,1/2,0)'
      c(118)=' 30 P2an       MY::(1/2,1/2,0) MZ::(1/2,1/2,0)'
      c(119)=' 30 P2na       MY::(1/2,0,1/2) MZ::(1/2,0,1/2)'
c     c(118)=' 30 P2an       MY::(1/2,1/2,0) MZ::(0,1/2,1/2)'
c     c(119)=' 30 P2na       MY::(0,1/2,1/2) MZ::(0,1/2,1/2)'
      c(120)=' 31 Pmn2_1     MX MY::(1/2,0,1/2)'
      c(121)=' 31 Pnm2_1     MX::(0,1/2,1/2) MY'
      c(122)=' 31 Pn2_1m     MX::(0,1/2,1/2) MZ'
      c(123)=' 31 Pm2_1n     MX MZ::(1/2,1/2,0)'
      c(124)=' 31 P2_1mn     MY MZ::(1/2,1/2,0)'
      c(125)=' 31 P2_1nm     MY::(1/2,0,1/2) MZ'
      c(126)=' 32 Pba2       MX::(1/2,1/2,0) MY::(1/2,1/2,0)'
      c(127)=' 32 Pc2a       MX::(1/2,0,1/2) MZ::(1/2,0,1/2)'
      c(128)=' 32 P2cb       MY::(0,1/2,1/2) MZ::(0,1/2,1/2)'
      c(129)=' 33 Pna2_1     MX::(1/2,1/2,1/2) MY::(1/2,1/2,0)'
      c(130)=' 33 Pbn2_1     MX::(1/2,1/2,0) MY::(1/2,1/2,1/2)'
      c(131)=' 33 Pc2_1n     MX::(1/2,0,1/2) MZ::(1/2,1/2,1/2)'
      c(132)=' 33 Pn2_1a     MX::(1/2,1/2,1/2) MZ::(1/2,0,1/2)'
      c(133)=' 33 P2_1nb     MY::(1/2,1/2,1/2) MZ::(0,1/2,1/2)'
      c(134)=' 33 P2_1cn     MY::(0,1/2,1/2) MZ::(1/2,1/2,1/2)'
      c(135)=' 34 Pnn2       MX::(1/2,1/2,1/2) MY::(1/2,1/2,1/2)'
      c(136)=' 34 Pn2n       MX::(1/2,1/2,1/2) MZ::(1/2,1/2,1/2)'
      c(137)=' 34 P2nn       MY::(1/2,1/2,1/2) MZ::(1/2,1/2,1/2)'
      c(138)=' 35 Cmm2       MX MY'
      c(139)=' 35 Bm2m       MX MZ'
      c(140)=' 35 A2mm       MY MZ'
      c(141)=' 36 Cmc2_1     MX MY::(0,0,1/2)'
      c(142)=' 36 Ccm2_1     MX::(0,0,1/2) MY'
      c(143)=' 36 Bb2_1m     MX::(0,1/2,0) MZ'
      c(144)=' 36 Bm2_1b     MX MZ::(0,1/2,0)'
      c(145)=' 36 A2_1ma     MY MZ::(1/2,0,0)'
      c(146)=' 36 A2_1am     MY::(1/2,0,0) MZ'
      c(147)=' 37 Ccc2       MX::(0,0,1/2) MY::(0,0,1/2)'
      c(148)=' 37 Bb2b       MX::(0,1/2,0) MZ::(0,1/2,0)'
      c(149)=' 37 A2aa       MY::(1/2,0,0) MZ::(1/2,0,0)'
      c(150)=' 38 Amm2       MX MY'
      c(151)=' 38 Bmm2       MX MY'
      c(152)=' 38 Am2m       MX MZ'
      c(153)=' 38 Cm2m       MX MZ'
      c(154)=' 38 B2mm       MY MZ'
      c(155)=' 38 C2mm       MY MZ'
      c(156)=' 39 Abm2       MX::(0,1/2,0) MY::(0,1/2,0)'
      c(157)=' 39 Bma2       MX::(1/2,0,0) MY::(1/2,0,0)'
      c(158)=' 39 Cm2a       MX::(1/2,0,0) MZ::(1/2,0,0)'
      c(159)=' 39 Ac2m       MX::(0,1/2,0) MZ::(0,1/2,0)'
      c(160)=' 39 B2cm       MY::(1/2,0,0) MZ::(1/2,0,0)'
      c(161)=' 39 C2mb       MY::(1/2,0,0) MZ::(1/2,0,0)'
      c(162)=' 40 Ama2       MX::(1/2,0,0) MY::(1/2,0,0)'
      c(163)=' 40 Bbm2       MX::(0,1/2,0) MY::(0,1/2,0)'
      c(164)=' 40 Cc2m       MX::(0,0,1/2) MZ::(0,0,1/2)'
      c(165)=' 40 Am2a       MX::(1/2,0,0) MZ::(1/2,0,0)'
      c(166)=' 40 B2mb       MY::(0,1/2,0) MZ::(0,1/2,0)'
      c(167)=' 40 C2cm       MY::(0,0,1/2) MZ::(0,0,1/2)'
      c(168)=' 41 Aba2       MX::(1/2,1/2,0) MY::(1/2,1/2,0)'
      c(169)=' 41 Bba2       MX::(1/2,1/2,0) MY::(1/2,1/2,0)'
      c(170)=' 41 Cc2a       MX::(1/2,0,1/2) MZ::(1/2,0,1/2)'
      c(171)=' 41 Ac2a       MX::(1/2,1/2,0) MZ::(1/2,1/2,0)'
      c(172)=' 41 B2cb       MY::(1/2,1/2,0) MZ::(1/2,1/2,0)'
      c(173)=' 41 C2cb       MY::(1/2,0,1/2) MZ::(1/2,0,1/2)'
      c(174)=' 42 Fmm2       MX MY'
      c(175)=' 42 Fm2m       MX MZ'
      c(176)=' 42 F2mm       MY MZ'
      c(177)=' 43 Fdd2       MX::(1/4,1/4,1/4) MY::(1/4,1/4,1/4)'
      c(178)=' 43 Fd2d       MX::(1/4,1/4,1/4) MZ::(1/4,1/4,1/4)'
      c(179)=' 43 F2dd       MY::(1/4,1/4,1/4) MZ::(1/4,1/4,1/4)'
      c(180)=' 44 Imm2       MX MY'
      c(181)=' 44 Im2m       MX MZ'
      c(182)=' 44 I2mm       MY MZ'
      c(183)=' 45 Iba2       MX::(0,0,1/2) MY::(0,0,1/2)'
      c(184)=' 45 Ic2a       MX::(0,1/2,0) MZ::(0,1/2,0)'
      c(185)=' 45 I2cb       MY::(1/2,0,0) MZ::(1/2,0,0)'
      c(186)=' 46 Ima2       MX::(1/2,0,0) MY::(1/2,0,0)'
      c(187)=' 46 Ibm2       MX::(0,1/2,0) MY::(0,1/2,0)'
      c(188)=' 46 Ic2m       MX::(0,0,1/2) MZ::(0,0,1/2)'
      c(189)=' 46 Im2a       MX::(1/2,0,0) MZ::(1/2,0,0)'
      c(190)=' 46 I2mb       MY::(0,1/2,0) MZ::(0,1/2,0)'
      c(191)=' 46 I2cm       MY::(0,0,1/2) MZ::(0,0,1/2)'
      c(192)=' 47 Pmmm       I MX MY'
      c(193)=' 481Pnnn       I::(1/2,1/2,1/2) MX::(1/2,1/2,1/2) '//
     .                      'MY::(1/2,1/2,1/2)'
      c(194)=' 482Pnnn       I MX::(0,1/2,1/2) MY::(1/2,0,1/2)'
      c(195)=' 49 Pccm       I MX::(0,0,1/2) MY::(0,0,1/2)'
      c(196)=' 49 Pmaa       I MX MY::(1/2,0,0)'
      c(197)=' 49 Pbmb       I MX::(0,1/2,0) MY'
      c(198)=' 501Pban       I::(1/2,1/2,0) MX::(1/2,1/2,0) '//
     .                      'MY::(1/2,1/2,0)'
      c(199)=' 502Pban       I MX::(0,1/2,0) MY::(1/2,0,0)'
      c(200)=' 501Pncb       I::(0,1/2,1/2) MX::(0,1/2,1/2) '//
     .                      'MY::(0,1/2,1/2)'
      c(201)=' 502Pncb       I MX::(0,1/2,1/2) MY::(0,0,1/2)'
      c(202)=' 501Pcna       I::(1/2,0,1/2) MX::(1/2,0,1/2) '//
     .                      'MY::(1/2,0,1/2)'
      c(203)=' 502Pcna       I MX::(0,0,1/2) MY::(1/2,0,1/2)'
      c(204)=' 51 Pmma       I MX::(1/2,0,0) MY'
      c(205)=' 51 Pmmb       I MX MY::(0,1/2,0)'
      c(206)=' 51 Pbmm       I MX::(0,1/2,0) MY::(0,1/2,0)'
      c(207)=' 51 Pcmm       I MX::(0,0,1/2) MY'
      c(208)=' 51 Pmcm       I MX MY::(0,0,1/2)'
      c(209)=' 51 Pmam       I MX::(1/2,0,0) MY::(1/2,0,0)'
      c(210)=' 52 Pnna       I MX::(0,1/2,1/2) MY::(1/2,1/2,1/2)'
      c(211)=' 52 Pnnb       I MX::(1/2,1/2,1/2) MY::(1/2,0,1/2)'
      c(212)=' 52 Pbnn       I MX::(0,1/2,0) MY::(1/2,0,1/2)'
      c(213)=' 52 Pcnn       I MX::(0,0,1/2) MY::(1/2,1/2,1/2)'
      c(214)=' 52 Pncn       I MX::(1/2,1/2,1/2) MY::(0,0,1/2)'
      c(215)=' 52 Pnan       I MX::(0,1/2,1/2) MY::(1/2,0,0)'
      c(216)=' 53 Pmna       I MX MY::(1/2,0,1/2)'
      c(217)=' 53 Pnmb       I MX::(0,1/2,1/2) MY'
      c(218)=' 53 Pbmn       I MX::(1/2,1/2,0) MY'
      c(219)=' 53 Pcnm       I MX::(1/2,0,1/2) MY::(1/2,0,1/2)'
      c(220)=' 53 Pncm       I MX::(0,1/2,1/2) MY::(0,1/2,1/2)'
      c(221)=' 53 Pman       I MX MY::(1/2,1/2,0)'
      c(222)=' 54 Pcca       I MX::(1/2,0,1/2) MY::(0,0,1/2)'
      c(223)=' 54 Pccb       I MX::(0,0,1/2) MY::(0,1/2,1/2)'
      c(224)=' 54 Pbaa       I MX::(0,1/2,0) MY::(1/2,1/2,0)'
      c(225)=' 54 Pcaa       I MX::(0,0,1/2) MY::(1/2,0,0)'
      c(226)=' 54 Pbcb       I MX::(0,1/2,0) MY::(0,0,1/2)'
      c(227)=' 54 Pbab       I MX::(1/2,1/2,0) MY::(1/2,0,0)'
      c(228)=' 55 Pbam       I MX::(1/2,1/2,0) MY::(1/2,1/2,0)'
      c(229)=' 55 Pmcb       I MX MY::(0,1/2,1/2)'
      c(230)=' 55 Pcma       I MX::(1/2,0,1/2) MY'
      c(231)=' 56 Pccn       I MX::(1/2,0,1/2) MY::(0,1/2,1/2)'
      c(232)=' 56 Pnaa       I MX::(0,1/2,1/2) MY::(1/2,1/2,0)'
      c(233)=' 56 Pbnb       I MX::(1/2,1/2,0) MY::(1/2,0,1/2)'
      c(234)=' 57 Pbcm       I MX::(0,1/2,0) MY::(0,1/2,1/2)'
      c(235)=' 57 Pcam       I MX::(1/2,0,1/2) MY::(1/2,0,0)'
      c(236)=' 57 Pmca       I MX::(1/2,0,0) MY::(0,0,1/2)'
      c(237)=' 57 Pmab       I MX::(1/2,0,0) MY::(1/2,1/2,0)'
      c(238)=' 57 Pbma       I MX::(1/2,1/2,0) MY::(0,1/2,0)'
      c(239)=' 57 Pcmb       I MX::(0,0,1/2) MY::(0,1/2,0)'
      c(240)=' 58 Pnnm       I MX::(1/2,1/2,1/2) MY::(1/2,1/2,1/2)'
      c(241)=' 58 Pmnn       I MX MY::(1/2,1/2,1/2)'
      c(242)=' 58 Pnmn       I MX::(1/2,1/2,1/2) MY'
      c(243)=' 591Pmmn       I::(1/2,1/2,0) MX MY'
      c(244)=' 592Pmmn       I MX::(1/2,0,0) MY::(0,1/2,0)'
      c(245)=' 591Pmnm       I::(1/2,0,1/2) MX MY::(1/2,0,1/2)'
      c(246)=' 592Pmnm       I MX::(1/2,0,0) MY::(1/2,0,1/2)'
      c(247)=' 591Pnmm       I::(0,1/2,1/2) MX::(0,1/2,1/2) MY'
      c(248)=' 592Pnmm       I MX::(0,1/2,1/2) MY::(0,1/2,0)'
      c(249)=' 60 Pbcn       I MX::(1/2,1/2,0) MY::(0,0,1/2)'
      c(250)=' 60 Pcan       I MX::(0,0,1/2) MY::(1/2,1/2,0)'
      c(251)=' 60 Pnca       I MX::(1/2,1/2,1/2) MY::(0,1/2,1/2)'
      c(252)=' 60 Pnab       I MX::(1/2,1/2,1/2) MY::(1/2,0,0)'
      c(253)=' 60 Pbna       I MX::(0,1/2,0) MY::(1/2,1/2,1/2)'
      c(254)=' 60 Pcnb       I MX::(1/2,0,1/2) MY::(1/2,1/2,1/2)'
      c(255)=' 61 Pbca       I MX::(1/2,1/2,0) MY::(0,1/2,1/2)'
      c(256)=' 61 Pcab       I MX::(1/2,0,1/2) MY::(1/2,1/2,0)'
      c(257)=' 62 Pnma       I MX::(1/2,1/2,1/2) MY::(0,1/2,0)'
      c(258)=' 62 Pmnb       I MX::(1/2,0,0) MY::(1/2,1/2,1/2)'
      c(259)=' 62 Pbnm       I MX::(1/2,1/2,0) MY::(1/2,1/2,1/2)'
      c(260)=' 62 Pcmn       I MX::(1/2,0,1/2) MY::(0,1/2,0)'
      c(261)=' 62 Pmcn       I MX::(1/2,0,0) MY::(0,1/2,1/2)'
      c(262)=' 62 Pnam       I MX::(1/2,1/2,1/2) MY::(1/2,1/2,0)'
      c(263)=' 63 Cmcm       I MX MY::(0,0,1/2)'
      c(264)=' 63 Ccmm       I MX::(0,0,1/2) MY'
      c(265)=' 63 Amma       I MX::(1/2,0,0) MY'
      c(266)=' 63 Amam       I MX::(1/2,0,0) MY::(1/2,0,0)'
      c(267)=' 63 Bbmm       I MX::(0,1/2,0) MY::(0,1/2,0)'
      c(268)=' 63 Bmmb       I MX MY::(0,1/2,0)'
      c(269)=' 64 Cmca       I MX MY::(1/2,0,1/2)'
      c(270)=' 64 Ccmb       I MX::(1/2,0,1/2) MY'
      c(271)=' 64 Abma       I MX::(1/2,1/2,0) MY'
      c(272)=' 64 Acam       I MX::(1/2,1/2,0) MY::(1/2,1/2,0)'
      c(273)=' 64 Bbcm       I MX::(1/2,1/2,0) MY::(1/2,1/2,0)'
      c(274)=' 64 Bmab       I MX MY::(1/2,1/2,0)'
      c(275)=' 65 Cmmm       I MX MY'
      c(276)=' 65 Ammm       I MX MY'
      c(277)=' 65 Bmmm       I MX MY'
      c(278)=' 66 Cccm       I MX::(0,0,1/2) MY::(0,0,1/2)'
      c(279)=' 66 Amaa       I MX MY::(1/2,0,0)'
      c(280)=' 66 Bbmb       I MX::(0,1/2,0) MY'
      c(281)=' 67 Cmma       I MX MY::(1/2,0,0)'
      c(282)=' 67 Cmmb       I MX::(1/2,0,0) MY'
      c(283)=' 67 Abmm       I MX::(0,1/2,0) MY'
      c(284)=' 67 Acmm       I MX::(0,1/2,0) MY::(0,1/2,0)'
      c(285)=' 67 Bmcm       I MX::(1/2,0,0) MY::(1/2,0,0)'
      c(286)=' 67 Bmam       I MX MY::(1/2,0,0)'
      c(287)=' 681Ccca       I::(1/2,0,1/2) MX::(1/2,0,1/2) '//
     .                      'MY::(1/2,0,1/2)'
      c(288)=' 682Ccca       I MX::(1/2,0,1/2) MY::(0,0,1/2)'
      c(289)=' 681Cccb       I::(1/2,0,1/2) MX::(1/2,0,1/2) '//
     .                      'MY::(1/2,0,1/2)'
      c(290)=' 682Cccb       I MX::(0,0,1/2) MY::(1/2,0,1/2)'
      c(291)=' 681Abaa       I::(1/2,1/2,0) MX::(1/2,1/2,0) '//
     .                      'MY::(1/2,1/2,0)'
      c(292)=' 682Abaa       I MX::(0,1/2,0) MY::(1/2,1/2,0)'
      c(293)=' 681Acaa       I::(1/2,1/2,0) MX::(1/2,1/2,0) '//
     .                      'MY::(1/2,1/2,0)'
      c(294)=' 682Acaa       I MX::(0,1/2,0) MY::(1/2,0,0)'
      c(295)=' 681Bbcb       I::(1/2,1/2,0) MX::(1/2,1/2,0) '//
     .                      'MY::(1/2,1/2,0)'
      c(296)=' 682Bbcb       I MX::(0,1/2,0) MY::(1/2,0,0)'
      c(297)=' 681Bbab       I::(1/2,1/2,0) MX::(1/2,1/2,0) '//
     .                      'MY::(1/2,1/2,0)'
      c(298)=' 682Bbab       I MX::(1/2,1/2,0) MY::(1/2,0,0)'
      c(299)=' 69 Fmmm       I MX MY'
      c(300)=' 701Fddd       I::(1/4,1/4,1/4) MX::(1/4,1/4,1/4) '//
     .                      'MY::(1/4,1/4,1/4)'
      c(301)=' 702Fddd       I MX::(0,1/4,1/4) MY::(1/4,0,1/4)'
      c(302)=' 71 Immm       I MX MY'
      c(303)=' 72 Ibam       I MX::(0,0,1/2) MY::(0,0,1/2)'
      c(304)=' 72 Imcb       I MX MY::(1/2,0,0)'
      c(305)=' 72 Icma       I MX::(0,1/2,0) MY'
      c(306)=' 73 Ibca       I MX::(0,0,1/2) MY::(1/2,0,0)'
      c(307)=' 73 Icab       I MX::(0,1/2,0) MY::(0,0,1/2)'
      c(308)=' 74 Imma       I MX MY::(0,1/2,0)'
      c(309)=' 74 Immb       I MX::(1/2,0,0) MY'
      c(310)=' 74 Ibmm       I MX::(0,0,1/2) MY'
      c(311)=' 74 Icmm       I MX::(0,1/2,0) MY::(0,1/2,0)'
      c(312)=' 74 Imcm       I MX::(1/2,0,0) MY::(1/2,0,0)'
      c(313)=' 74 Imam       I MX MY::(0,0,1/2)'
      c(314)=' 75 P4         R4Z'
      c(315)=' 76 P4_1       R4Z::(0,0,1/4)'
      c(316)=' 77 P4_2       R4Z::(0,0,1/2)'
      c(317)=' 78 P4_3       R4Z::(0,0,-1/4)'
      c(318)=' 79 I4         R4Z'
      c(319)=' 80 I4_1       R4Z::(1/2,0,-1/4)'
      c(320)=' 81 P-4        I*R4Z'
      c(321)=' 82 I-4        I*R4Z'
      c(322)=' 83 P4/m       I R4Z'
      c(323)=' 84 P4_2/m     I R4Z::(0,0,1/2)'
      c(324)=' 851P4/n       I::(1/2,1/2,0) R4Z::(1/2,1/2,0)'
      c(325)=' 852P4/n       I R4Z::(1/2,0,0)'
      c(326)=' 861P4_2/n     I::(1/2,1/2,1/2) R4Z::(1/2,1/2,1/2)'
      c(327)=' 862P4_2/n     I R4Z::(0,1/2,1/2)'
      c(328)=' 87 I4/m       I R4Z'
      c(329)=' 881I4_1/a     I::(1/2,0,-1/4) R4Z::(1/2,0,-1/4)'
      c(330)=' 882I4_1/a     I R4Z::(-1/4,1/4,1/4)'
      c(331)=' 89 P422       R4Z R2X'
      c(332)=' 90 P42_12     R4Z::(1/2,1/2,0) R2X::(1/2,1/2,0)'
      c(333)=' 91 P4_122     R4Z::(0,0,1/4) R2X::(0,0,1/2)'
      c(334)=' 92 P4_12_12   R4Z::(1/2,1/2,1/4) R2X::(1/2,1/2,-1/4)'
      c(335)=' 93 P4_222     R4Z::(0,0,1/2) R2X'
      c(336)=' 94 P4_22_12   R4Z::(1/2,1/2,1/2) R2X::(1/2,1/2,1/2)'
      c(337)=' 95 P4_322     R4Z::(0,0,-1/4) R2X::(0,0,1/2)'
      c(338)=' 96 P4_32_12   R4Z::(1/2,1/2,-1/4) R2X::(1/2,1/2,1/4)'
      c(339)=' 97 I422       R4Z R2X'
      c(340)=' 98 I4_122     R4Z::(1/2,0,-1/4) R2X::(1/2,0,-1/4)'
      c(341)=' 99 P4mm       R4Z MX'
      c(342)='100 P4bm       R4Z MX::(1/2,1/2,0)'
      c(343)='101 P4_2cm     R4Z::(0,0,1/2) MX::(0,0,1/2)'
      c(344)='102 P4_2nm     R4Z::(1/2,1/2,1/2) MX::(1/2,1/2,1/2)'
      c(345)='103 P4cc       R4Z MX::(0,0,1/2)'
      c(346)='104 P4nc       R4Z MX::(1/2,1/2,1/2)'
      c(347)='105 P4_2mc     R4Z::(0,0,1/2) MX'
      c(348)='106 P4_2bc     R4Z::(0,0,1/2) MX::(1/2,1/2,0)'
      c(349)='107 I4mm       R4Z MX'
      c(350)='108 I4cm       R4Z MX::(0,0,1/2)'
      c(351)='109 I4_1md     R4Z::(1/2,0,-1/4) MX'
      c(352)='110 I4_1cd     R4Z::(1/2,0,-1/4) MX::(0,0,1/2)'
      c(353)='111 P-42m      I*R4Z R2X'
      c(354)='112 P-42c      I*R4Z R2X::(0,0,1/2)'
      c(355)='113 P-42_1m    I*R4Z R2X::(1/2,1/2,0)'
      c(356)='114 P-42_1c    I*R4Z R2X::(1/2,1/2,1/2)'
      c(357)='115 P-4m2      I*R4Z MX'
      c(358)='116 P-4c2      I*R4Z MX::(0,0,1/2)'
      c(359)='117 P-4b2      I*R4Z MX::(1/2,1/2,0)'
      c(360)='118 P-4n2      I*R4Z MX::(1/2,1/2,1/2)'
      c(361)='119 I-4m2      I*R4Z MX'
      c(362)='120 I-4c2      I*R4Z MX::(0,0,1/2)'
      c(363)='121 I-42m      I*R4Z R2X'
      c(364)='122 I-42d      I*R4Z R2X::(1/2,0,-1/4)'
      c(365)='123 P4/mmm     I R4Z MX'
      c(366)='124 P4/mcc     I R4Z MX::(0,0,1/2)'
      c(367)='1251P4/nbm     I::(1/2,1/2,0) R4Z MX::(1/2,1/2,0)'
      c(368)='1252P4/nbm     I R4Z::(1/2,0,0) MX::(0,1/2,0)'
      c(369)='1261P4/nnc     I::(1/2,1/2,1/2) R4Z MX::(1/2,1/2,1/2)'
      c(370)='1262P4/nnc     I R4Z::(1/2,0,0) MX::(0,1/2,1/2)'
      c(371)='127 P4/mbm     I R4Z MX::(1/2,1/2,0)'
      c(372)='128 P4/mnc     I R4Z MX::(1/2,1/2,1/2)'
      c(373)='1291P4/nmm     I::(1/2,1/2,0) R4Z::(1/2,1/2,0) MX'
      c(374)='1292P4/nmm     I R4Z::(1/2,0,0) MX::(1/2,0,0)'
      c(375)='1301P4/ncc     I::(1/2,1/2,0) R4Z::(1/2,1/2,0) '//
     .                      'MX::(0,0,1/2)'
      c(376)='1302P4/ncc     I R4Z::(1/2,0,0) MX::(1/2,0,1/2)'
      c(377)='131 P4_2/mmc   I R4Z::(0,0,1/2) MX'
      c(378)='132 P4_2/mcm   I R4Z::(0,0,1/2) MX::(0,0,1/2)'
      c(379)='1331P4_2/nbc   I::(1/2,1/2,1/2) R4Z::(1/2,1/2,1/2) '//
     .                      'MX::(1/2,1/2,0)'
      c(380)='1332P4_2/nbc   I R4Z::(1/2,0,1/2) MX::(0,1/2,0)'
      c(381)='1341P4_2/nnm   I::(1/2,1/2,1/2) R4Z::(1/2,1/2,1/2) '//
     .                      'MX::(1/2,1/2,1/2)'
      c(382)='1342P4_2/nnm   I R4Z::(1/2,0,1/2) MX::(0,1/2,1/2)'
      c(383)='135 P4_2/mbc   I R4Z::(0,0,1/2) MX::(1/2,1/2,0)'
      c(384)='136 P4_2/mnm   I R4Z::(1/2,1/2,1/2) MX::(1/2,1/2,1/2)'
      c(385)='1371P4_2/nmc   I::(1/2,1/2,1/2) R4Z::(1/2,1/2,1/2) MX'
      c(386)='1372P4_2/nmc   I R4Z::(1/2,0,1/2) MX::(1/2,0,0)'
      c(387)='1381P4_2/ncm   I::(1/2,1/2,1/2) R4Z::(1/2,1/2,1/2) '//
     .                      'MX::(0,0,1/2)'
      c(388)='1382P4_2/ncm   I R4Z::(1/2,0,1/2) MX::(1/2,0,1/2)'
      c(389)='139 I4/mmm     I R4Z MX'
      c(390)='140 I4/mcm     I R4Z MX::(0,0,1/2)'
      c(391)='1411I4_1/amd   I::(1/2,0,-1/4) R4Z::(1/2,0,-1/4) MX'
      c(392)='1412I4_1/amd   I R4Z::(1/4,-1/4,1/4) MX'
      c(393)='1421I4_1/acd   I::(1/2,0,-1/4) R4Z::(1/2,0,-1/4) '//
     .                      'MX::(0,0,1/2)'
      c(394)='1422I4_1/acd   I R4Z::(1/4,-1/4,1/4) MX::(0,0,1/2)'
      c(395)='143 P3         R3Z'
      c(396)='143 P311       R3Z'
      c(397)='144 P3_1       R3Z::(0,0,1/3)'
      c(398)='144 P3_111     R3Z::(0,0,1/3)'
      c(399)='145 P3_2       R3Z::(0,0,-1/3)'
      c(400)='145 P3_211     R3Z::(0,0,-1/3)'
      c(401)='146 R3         R3Z'
      c(402)='146 P31        R3Z'
      c(403)='147 P-3        I R3Z'
      c(404)='147 P-311      I R3Z'
      c(405)='148 R-3        I R3Z'
      c(406)='148 P-31       I R3Z'
      c(407)='149 P312       R3Z R2X'
      c(408)='150 P321       R3Z R2Y'
      c(409)='151 P3_112     R3Z::(0,0,1/3) R2X'
      c(410)='152 P3_121     R3Z::(0,0,1/3) R2Y::(0,0,1/3)'
      c(411)='153 P3_212     R3Z::(0,0,-1/3) R2X'
      c(412)='154 P3_221     R3Z::(0,0,-1/3) R2Y::(0,0,-1/3)'
      c(413)='155 R32        R3Z R2Y'
      c(414)='155 P32        R3Z R2Y'
      c(415)='156 P3m1       R3Z MY'
      c(416)='157 P31m       R3Z MX'
      c(417)='158 P3c1       R3Z MY::(0,0,1/2)'
      c(418)='159 P31c       R3Z MX::(0,0,1/2)'
      c(419)='160 R3m        R3Z MY'
      c(420)='160 P3m        R3Z MY'
      c(421)='161 R3c        R3Z MY::(0,0,1/2)'
      c(422)='161 P3c        R3Z MY::(1/2,1/2,1/2)'
      c(423)='162 P-31m      I R3Z MX'
      c(424)='163 P-31c      I R3Z MX::(0,0,1/2)'
      c(425)='164 P-3m1      I R3Z MY'
      c(426)='165 P-3c1      I R3Z MY::(0,0,1/2)'
      c(427)='166 R-3m       I R3Z MY'
      c(428)='166 P-3m       I R3Z MY'
      c(429)='167 R-3c       I R3Z MY::(0,0,1/2)'
      c(430)='167 P-3c       I R3Z MY::(1/2,1/2,1/2)'
      c(431)='168 P6         R6Z'
      c(432)='169 P6_1       R6Z::(0,0,1/6)'
      c(433)='170 P6_5       R6Z::(0,0,-1/6)'
      c(434)='171 P6_2       R6Z::(0,0,1/3)'
      c(435)='172 P6_4       R6Z::(0,0,-1/3)'
      c(436)='173 P6_3       R6Z::(0,0,1/2)'
      c(437)='174 P-6        I*R6Z'
      c(438)='175 P6/m       I R6Z'
      c(439)='176 P6_3/m     I R6Z::(0,0,1/2)'
      c(440)='177 P622       R6Z R2X'
      c(441)='178 P6_122     R6Z::(0,0,1/6) R2X::(0,0,1/6)'
      c(442)='179 P6_522     R6Z::(0,0,-1/6) R2X::(0,0,-1/6)'
      c(443)='180 P6_222     R6Z::(0,0,1/3) R2X::(0,0,1/3)'
      c(444)='181 P6_422     R6Z::(0,0,-1/3) R2X::(0,0,-1/3)'
      c(445)='182 P6_322     R6Z::(0,0,1/2) R2X::(0,0,1/2)'
      c(446)='183 P6mm       R6Z MX'
      c(447)='184 P6cc       R6Z MX::(0,0,1/2)'
      c(448)='185 P6_3cm     R6Z::(0,0,1/2) MX'
      c(449)='186 P6_3mc     R6Z::(0,0,1/2) MX::(0,0,1/2)'
      c(450)='187 P-6m2      I*R6Z MY'
      c(451)='188 P-6c2      I*R6Z::(0,0,1/2) MY::(0,0,1/2)'
      c(452)='189 P-62m      I*R6Z MX'
      c(453)='190 P-62c      I*R6Z::(0,0,1/2) MX::(0,0,1/2)'
      c(454)='191 P6/mmm     I R6Z MX'
      c(455)='192 P6/mcc     I R6Z MX::(0,0,1/2)'
      c(456)='193 P6_3/mcm   I R6Z::(0,0,1/2) MX'
      c(457)='194 P6_3/mmc   I R6Z::(0,0,1/2) MX::(0,0,1/2)'
      c(458)='195 P23        R3D R2X'
      c(459)='196 F23        R3D R2X'
      c(460)='197 I23        R3D R2X'
      c(461)='198 P2_13      R3D R2X::(1/2,1/2,0)'
      c(462)='199 I2_13      R3D R2X::(0,0,1/2)'
      c(463)='200 Pm-3       I R3D MX'
      c(464)='2011Pn-3       I::(1/2,1/2,1/2) R3D MX::(1/2,1/2,1/2)'
      c(465)='2012Pn-3       I R3D MX::(0,1/2,1/2)'
      c(466)='202 Fm-3       I R3D MX'
      c(467)='2031Fd-3       I::(1/4,1/4,1/4) R3D MX::(1/4,1/4,1/4)'
      c(468)='2032Fd-3       I R3D MX::(0,1/4,1/4)'
      c(469)='204 Im-3       I R3D MX'
      c(470)='205 Pa-3       I R3D MX::(1/2,1/2,0)'
      c(471)='206 Ia-3       I R3D MX::(0,0,1/2)'
      c(472)='207 P432       R4X R3D'
      c(473)='208 P4_232     R4X::(1/2,1/2,1/2) R3D'
      c(474)='209 F432       R4X R3D'
      c(475)='210 F4_132     R4X::(1/4,1/4,1/4) R3D'
      c(476)='211 I432       R4X R3D'
      c(477)='212 P4_332     R4X::(-1/4,-1/4,1/4) R3D'
      c(478)='213 P4_132     R4X::(1/4,1/4,-1/4) R3D'
      c(479)='214 I4_132     R4X::(1/4,1/4,-1/4) R3D'
      c(480)='215 P-43m      I*R4X R3D'
      c(481)='216 F-43m      I*R4X R3D'
      c(482)='217 I-43m      I*R4X R3D'
      c(483)='218 P-43n      I*R4X::(1/2,1/2,1/2) R3D'
      c(484)='219 F-43c      I*R4X::(1/2,1/2,1/2) R3D'
      c(485)='220 I-43d      I*R4X::(1/4,1/4,-1/4) R3D'
      c(486)='221 Pm-3m      I R4X R3D'
      c(487)='2221Pn-3n      I::(1/2,1/2,1/2) R4X R3D'
      c(488)='2222Pn-3n      I R4X::(0,1/2,0) R3D'
      c(489)='223 Pm-3n      I R4X::(1/2,1/2,1/2) R3D'
      c(490)='2241Pn-3m      I::(1/2,1/2,1/2) R4X::(1/2,1/2,1/2) R3D'
      c(491)='2242Pn-3m      I R4X::(1/2,0,1/2) R3D'
      c(492)='225 Fm-3m      I R4X R3D'
      c(493)='226 Fm-3c      I R4X::(1/2,1/2,1/2) R3D'
      c(494)='2271Fd-3m      I::(1/4,1/4,1/4) R4X::(1/4,1/4,1/4) R3D'
      c(495)='2272Fd-3m      I R4X::(1/4,0,1/4) R3D'
      c(496)='2281Fd-3c      I::(-1/4,-1/4,-1/4) R4X::(1/4,1/4,1/4) R3D'
      c(497)='2282Fd-3c      I R4X::(1/4,1/2,1/4) R3D'
      c(498)='229 Im-3m      I R4X R3D'
      c(499)='230 Ia-3d      I R4X::(1/4,1/4,-1/4) R3D'

      stdo = lgunit(1)


C --- Supply generators ---
      if (iopt .eq. -1) then

        call chcase(1,1,instr)
        call chcase(-1,10,instr(2:))
        isym7=0
        csym7=' '
        idum=0
C       Parse for integer group id.  Failing that, assign to string
        if (.not. a2bin(instr,isym7,2,0,' ',idum,-1)) csym7 = instr

        if ((isym7.ge.1.and.isym7.le.230) .or. csym7.ne.' ') then
C     ... Loop over nelts possibilities until grp id matches
        do  isym8 = 1, nelts
          read(c(isym8),200) isym(7),corig,csym(7),symgrp

C         Group found if integer ID or string ID matches
          if (isym7.eq.isym(7) .or. csym7.eq.csym(7)) then

C           For some groups, there is an alternate origin
            if (corig .eq. '1') then
              call query('choose origin for this space group (1 or 2):',
     .          2,iorig)
              if (iorig .eq. 1) then
                continue
              elseif (iorig .eq. 2) then
                read(c(isym8+1),200) isym(7),corig,csym(7),symgrp
                if (corig.ne.'2') call rxs(' SPCGRP: bad corig,',corig)
              else
                call isanrg(iorig,1,2,'SPCGRP:','iorig',.true.)
              endif
            elseif (corig .ne. ' ') then
              call rxs(' SPCGRP: bad corig,',corig)
            endif

C       ... Get isym(6-1) and csym(6-1)
            isym(6) = isym6(isym(7))
            isym(5) = isym5(isym(6))
            isym(4) = isym4(isym(5))
            do  isym3 = 1, 7
              if (csym(7)(1:1) .eq. csym3(isym3)) isym(3)=isym3
            enddo
            isym(1) = isym1(isym(7))
            isym(2) = isym2(isym(3),isym(1))
            call gtcsym(csym,isym)

C       ... Printout
            write(stdo,300) csym(7),isym(7),csym(5)
            if (isym(5) .eq. 5) then
              if (isym(1).eq.5 .and. isym(3).eq.1) write(*,303)
              if (isym(1).eq.5 .and. isym(3).eq.7) write(*,304)
              if (isym(1).eq.6) write(*,305)
            endif

C       ... Copy symgrp into gens
            gens = symgrp
            call words(symgrp,ngen)
            write(stdo,302) symgrp
            return
          endif
        enddo
        endif
        isym(7) = 0
        write(stdo,301) instr

C --- Determine space group number ---
      elseif (iopt .eq. 1) then
        do  isym8 = 1, nelts
          read(c(isym8),200) isym(7),corig,csym(7),symgrp
          if(isym(6).eq.isym6(isym(7)) .and.
     .       isym(1).eq.isym1(isym(7)) .and.
     .       csym(3).eq.csym(7)(1:1)) then
            iorig = 1
            if (corig .eq. '2') iorig = 2
            if (iprint() .ge. 110) write(stdo,306) symgrp
            call rx('gengrp: iopt=1 not implemented')
C            call strcp2(53,72,ngen,symgrp,gens1)
C            call gtorig(dorig,gens,gens1,info,ngen,plat,qlat)
C            if (info .eq. 0) return
          endif
        enddo
        csym(7) = 'undefined'
        isym(7) = 0
      endif
  200 format(i3,a1,a10,1x,a53)

  300 format(5x,'Space group: ',a10, ' No.:',i3,
     .      /5x,'Crystal system: ',a)
  301 format(5x,'GENGRP: no such space group, ',a)
  302 format(5x,'Generators: ',a53)
  303 format(5x,'Bravais system: rhombohedral; rhombohedral axes')
  304 format(5x,'Bravais system: rhombohedral; hexagonal axes')
  305 format(5x,'Bravais system: hexagonal')
  306 format(' SPCGRP: try ',a53)
      end

      subroutine gtcsym(csym,isym)
C- Gets names characterizing the symmetry of lattice and crystal
C ----------------------------------------------------------------------
Ci Inputs/Outputs:
Ci  isym  :numbers characterizing the symmetry of lattice and crystal
Co         :any isym(1..6) outside allowable range is initialized to zero
Co Outputs:
Co   csym  :names characterizing the symmetry of lattice and crystal
Co         :csym(i) is a string labelling the integer index isym(i).
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer isym(6)
      character*30 csym(6)
C Local variables:
      integer i,n(6)
      character*12 csym1(0:7)
      character*30 csym2(0:14)
      character*1  csym3(0:7)
      character*25 csym4(0:6)
      character*12 csym5(0:7)
      character*5  csym6(0:32)
C Data statements:
      data n/7,14,7,6,7,32/
      data csym1/'undefinite','triclinic','monoclinic','orthorhombic',
     .          'tetragonal','rhombohedral','hexagonal','cubic'/
      data csym2/'undefinite','triclinic','monoclinic primitive',
     .          'monoclinic base-centred','orthorhombic primitive',
     .          'orthorhombic base-centred','orthorhombic body-centred',
     .          'orthorhombic face-centred','tetragonal primitive',
     .          'tetragonal body-centred','rhombohedral','hexagonal',
     .          'cubic primitive','cubic body-centred',
     .          'cubic face-centred'/
      data csym3/' ','P','C','B','A','I','F','R'/
      data csym4/'undefinite','triclinic','monoclinic','orthorhombic',
     .           'tetragonal','hexagonal','cubic'/
      data csym5/'undefinite','triclinic','monoclinic','orthorhombic',
     .           'tetragonal','trigonal','hexagonal','cubic'/
      data csym6/'undef','1','-1','2','m','2/m','222','mm2','mmm',
     .           '4','-4','4/m','422','4mm','-42m','4/mmm',
     .           '3','-3','32','3m','-3m',
     .           '6','-6','6/m','622','6mm','-62m','6/mmm',
     .           '23','m-3','432','-432','m-3m'/

      do  i = 1, 6
        if (isym(i).gt.n(i) .or. isym(i).lt.0) isym(i)=0
      enddo

      csym(1) = csym1(isym(1))
      csym(2) = csym2(isym(2))
      csym(3) = csym3(isym(3))
      csym(4) = csym4(isym(4))
      csym(5) = csym5(isym(5))
      csym(6) = csym6(isym(6))
      end
C      subroutine strcp2(lstrin,lstrou,nelts,strin,strout)
CC-Changes the format of the string containing the generators
CC ----------------------------------------------------------------------
CCi Inputs:
CCi   strin  :
CCi   lstrin :
CCo Outputs:
CCo   strout :
CCo   nelts  :
CCr Remarks:
CC ----------------------------------------------------------------------
C      implicit none
CC Passed variables:
C      integer nelts,lstrin,lstrou
C      character strin(lstrin)
C      character*(*) strout(*)
CC Local variables:
C      integer i,j
CC External calls:
C      external  skipbl,strcop
C
C      nelts=0
C      j=0
C      call skipbl(strin,lstrin,j)
C      j=j+1
C   10 if (j.lt.lstrin) then
C       nelts=nelts+1
C       strout(nelts)=' '
C       call strcop(strout(nelts),strin(j),min0(lstrou,lstrin-j+1),' ',i)
C       j=j+i-1
C       call skipbl(strin,lstrin,j)
C       j=j+1
C       goto 10
C      endif
C      end

