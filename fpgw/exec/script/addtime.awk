BEGIN {
   start=START

   poplist=0
}
/^\!TIME0/{ printf("       call realtimediff(%i,'')\n",start); push(start); start+=2;next }
/^\!TIME1/{ gsub("[!]TIME1","");i=pop(); 
            printf("       call realtimediff(%i,%s)\n",i,$0);next }
/^\!TIMESHOW/ { 
            print  "       call print_realtimediff()" ; next}
/^\!KINO/{ gsub("\!KINO","     "); print; next }
{ print }
func push(x) {
   list[poplist]=x 
   poplist++
}
func pop() {
   poplist--;
   return list[poplist]
}
