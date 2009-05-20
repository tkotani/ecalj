#  make dependency of makefile like
#  ../main/hsfp0.sc.m.o: ../gwsrc/readqg.o
#  
# > grep -i "^[ \t]*module" ../*/*.f > LIST_MODULE
#
# > grep -i "^[ \t]*use" ../*/*.f > LIST_USE
#
# > awk -f moduledep.awk 
#


BEGIN{

#   system("grep -i \"^[ \t]*module\" ../*/*.f > LIST_MODULE");
#   system("grep -i \"^[ \t]*use\" ../*/*.f > LIST_USE");

   make_list_module();
   make_list_use();

   make_makefiledependency(); 
}


func make_list_module()
{
 infile="LIST_MODULE";

 imodule=0;
 ret=1;
 while (ret!=0) {
  ret = getline <infile;

  sub(":"," ");  # delete : 

  if ($2=="module") {
  imodule++;
  modulename[imodule]=tolower($3)
  modulefile[imodule]=$1;
  }
  else {
    print "$2 is not use. Fix",infile;
    exit(10);
  }


 }
 nmodule=imodule

# test
 
# for (i=1;i<=nmodule;i++) {
#    print modulename[i], modulefile[i];
# }

}


func find_modulefile(name, i)
{
   for (i=1;i<=nmodule;i++) {
     if ( modulename[i]==name ) {
        return modulefile[i];
     }
   }
   print name, " can not  be found in module list"
    exit(10);

}

func make_list_use()
{

  infile="LIST_USE";

  nuse=0;

  ret=1;
  while (ret!=0) {
   ret=getline <infile;

   sub(":"," ");
   gsub(","," ");

   if ($2=="use") {
    inuse++;
    add_use_list( $1, tolower($3) );
#    print $1, $3;
   }
   else {
     print "$2 is not use. Fix",infile;
     print $0;
     exit(10);
   } 

  }

#  for (i=1;i<=nuse;i++) {
#     print i,  uselist[i], nusename[i];
#     for (j=1;j<=nusename[i];j++) {
#       print "    ",j,usefile[i,j] #, find_modulefile( usefile[i,j] );
#     }
#  }

}


func add_use_list(name0, name1)
{
  found=0;
  for (i=1;i<=nuse;i++) {
    if (uselist[i]==name0) {
       found=1;
       break;
    }
  }
  uselistid=i;

  if (found) {
    for (i=1;i<=nusename[uselistid];i++) {
       if ( usefile[uselistid,i]==name1 ) {
          return;
       }
    }
    nusename[uselistid]++;
    usefile[uselistid,nusename[uselistid]]=  name1;

  }
  else {
    nuse++;
    uselist[nuse]=name0
    nusename[nuse]=1;
    usefile[nuse,nusename[nuse]]=name1;
  }

}


func make_makefiledependency()
{


  for (i=1;i<=nuse;i++) {
#     print i,  uselist[i], nusename[i];

     for (j=1;j<=nusename[i];j++) {
#       print "    ",j,usefile[i,j], find_modulefile( usefile[i,j] );
        name0= uselist[i];
        name1= find_modulefile( usefile[i,j] );
        sub(".f$",".o",name0);
        sub(".f$",".o",name1);
        printf("%s: %s\n", name0, name1);

     }
  }

}

