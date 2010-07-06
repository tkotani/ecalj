C ... tests a2vec
      subroutine fmain
      implicit none
      character s*80
      integer ix(10),ncof,a2vec,ich,i2,i1mach,k
      double precision fcof(10)

      ich = 0
      s = '3,4:5'
      k = 5
      ncof = 3
      i2 = a2vec(s,len(s),ich,4,',: ',3,3,ncof,ix,fcof)
      call awrit4('strn: '//s(1:k)//' : ncof=%i a2vec=%i ix=%3:1i  '//
     .  'fcof=%3:1d',' ',80,i1mach(2),ncof,i2,ix,fcof)

      ich = 0
      s = 'x=3 x==4,x+5'
      k = 12
      ncof = -3
      i2 = a2vec(s,len(s),ich,4,',: ',3,3,ncof,ix,fcof)
      call awrit4('strn: '//s(1:k)//' : ncof=%i a2vec=%i ix=%3:1i  '//
     .  'fcof=%3:1d',' ',80,i1mach(2),ncof,i2,ix,fcof)

      ich = 0
      s = '9  8.8   7'
      k = 10
      ncof = 3
      i2 = a2vec(s,len(s),ich,4,',: ',3,1,ncof,ix,fcof)
      call awrit4('strn: '//s(1:k)//' : ncof=%i a2vec=%i ix=%3:1i  '//
     .  'fcof=%3:1d',' ',80,i1mach(2),ncof,i2,ix,fcof)

      end
