BEGIN{
comment="Cstop2rx 2013.08.09 kino"
#comment="------------------\n*"
}
/^ .*[ )]stop[ '"]/{
  if (stopprocess()) { next }
}
/^ .*[ )]stop$/{
  if (stopprocess()) { next }
}
/^ *stop /{
  if (stopprocess()) { next }
}
{ 
print 
}


func stopprocess() {

    if (match($0,"^ *call")) {
        return  0
    }
    i1=match($0,"stop *'")
    i2=match($0,"stop *\"")
    i3=match($0,"stop *!")
    i4=match($0,"stop *$")
  #  print i1,i2, i3,i4
    if (i1+i2+i3+i4==0 ) {
        return  0
    }

    msg="''"
    tail=""
    beforestop=""
    afterstop=""
    line=$0
    #print"---------------------------"
    #print FILENAME
    print  comment $0
    i0=match($0,"stop")
    beforestop=substr($0,1,i0-1)
    #print "beforestop=",beforestop"stop"

# after stop
    afterstop=substr($0,i0+4)
    #print "afterstop", afterstop
    i1=match(afterstop,"'")
    i1d=match(afterstop,"\"")
#print i1,i1d
    readsingle=0
    readdouble=0
    if (i1==0 && i1d>0) {
        readdouble=1
    }
    if (i1d==0 && i1>0){
        readsingle=1
    }
    if(i1d>0 && i1>0) {
        if (i1d>i1) {
            readsingle=1
        }
        else {
            readdouble=1
        }
    }
    if (readsingle>0 && readdouble>0) {
        print "error double",readsingle,readdouble
        exit(10)
    }
    if (readsingle==0 && readdouble==0)  {
        tail=afterstop
        writerx();
        readuntilcontinue()
        #print "tail=",tail
        return  1
    }

    if (readsingle) {
        s=substr(afterstop,i1+1)
        i2=index(s,"'")
    }
    else {
        s=substr(afterstop,i1d+1)
        i2=index(s,"\"")
    }
    if (readsingle) {
#        print i1,i2
    }
    else {
#        print i1d,i2
        i1=i1d
    }
    if ( (readsingle>0 || readdouble>0 ) && i2==0) {
        print comment,"continue line"
        msg=afterstop
        writerx(0); print""
        readuntilclosed()
        return  1
    }
    msg=substr(afterstop,i1,i2+1)
    ok=0
    if (match(msg,"OK!")) {  ok=1 }
    tail=substr(afterstop,i1+i2+1)
    #print "msg=",ok,msg
    if (length(tail)>0) {
        #print "tail=",tail
        if (tail=="//") {
            print comment,"continue line(2)"
            writerx(0); print tail
            readuntilclosed()
            return  1
        }
    }
    else {
       writerx(1)
       return  1
    }
    
    writerx(1)
    return 1
}


func writerx(end){
title=""
if (ok) {
printf("%s%s%s %s", title,beforestop,"call rx0(",msg)
}else {
printf("%s%s%s %s", title,beforestop,"call rx(",msg)
}
if (end) { print ")"}
}


func readuntilclosed()
{

    while (1) {
        getline
        if (readsingle)  {  i3=match($0,"'") }
        else  {  i3=match($0,"\"") }
        if (i3>0) {
           printf("%s )\n",$0)
           return
        }
        else {
           print
        }
    }

}

func readuntilcontinue()
{

    while (1) {
        getline
        if (match($0,"^     [^ ]")) {
             print "//"
             printf("%s",$0)
        }
        else {
             print ")"
             print 
             return
        }
    }

}
