      character*8 function charext(num)
      integer(4) ::num
      charext = char(48+mod(num,10))
      if(num>9)   charext= char(48+mod(num/10,10))//charext
      if(num>99)  charext= char(48+mod(num/100,10))//charext
      if(num>999) charext= char(48+mod(num/1000,10))//charext
      if(num>9999)charext= char(48+mod(num/10000,10))//charext
Cstop2rx 2013.08.09 kino      if(num >99999) stop ' charext:can not produce'
      if(num >99999) call rx( ' charext:can not produce')
      end
