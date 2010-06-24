BEGIN {
 IGNORECASE=1
 g_nlist=0
 read_struc_def()
}
{ 
 type="undef"
 for (i=1;i<=NF;i++) {
	if (match($i,"^'%',") ) {
		element= $(i+1)
		gsub("[',]","",element)
		gsub("\\]","",element)
	}
	if (match($i,"^'type\\(") ){ 
		structype=$i
		gsub("[',]","",structype)
		gsub("\\]","",structype)
		gsub("type\\(","",structype)
		gsub("[)]","",structype)
	}
	if (match($i,"^'real") || match($i,"^'complex") ||  match($i,"^'integer") ) {
		type=$i
		gsub("[']","",type)
		gsub(",","",type)
	}
 }
 list_push( structype,element,type)
}
END{
  list_print()
}


func list_print(    i,struc,element,type,comment)
{
    for (i=1;i<=g_nlist;i++) {
        struc=list[i,0]
        element=list[i,1]
        type= struc_element_type( struc,element)
        comment=""; comment2=""
        if ( match(element, "^o")  && type!="integer(8)" ) { comment="unchange"; comment2=type }
	print struc,element,list[i,2],"#",comment,comment2
    }

}

func list_push( structype,element,type,   i)
{
    for (i=0;i<=g_nlist;i++) {
        if (list[i,0]==structype && list[i,1]==element) {
		if (type=="undef") { return}
		if (list[i,2]=="undef")  {
			list[i,2]=type	
			return
		}
		else if (list[i,2]!=type) {
                	print "push another type",structype,element,list[i,2], "add type",type 
		} 
        }
    }
    g_nlist++
    i=g_nlist
    list[i,0]=structype
    list[i,1]=element
    list[i,2]=type
}

func read_struc_def( id,ret,struc,n,nlist)
{
 
 file="/home/kino/kit/GW/7K/ecalj_2010_0513/lm7K/subs/m_struc_def.F"
 id=0
 while (1) {
    ret=getline < file
    if (ret<=0) { print "error"; break; }
    if ( match($0,"^[^ ]") ) { continue;}
    if (NF==2 && $1=="type") {
        id=id+1
        struc=$2
        strucid[id]=struc
        while (1) {
            ret=getline < file   
            if (ret<=0) { break; }
	    if ( match($0,"^[^ ]") ) { continue;}
            if (NF>=2 && $1=="end" && $2=="type") { break }
            gsub("::"," :: ")
            n=split($0,a)
            nlist=list[id,0,0]
            nlist=nlist+1
	    list[id,nlist,1]=a[1] # type
	    list[id,nlist,3]=a[3] # name(size)
            name=a[3]
            gsub("\\(.*\\)","",name)
	    list[id,nlist,2]=name  # name without size
	    list[id,nlist,4]=a[4]  # comment
	    list[id,0,0]=nlist
        }
    }
    if ($1=="subroutine") {
	break;
    }
 } 
 close(file)

 g_nid=id
}

func struc_element_type(struc,element, id,n){
   id=strucname[struc]
   for (id=1;id<=g_nid;id++) {
        if (strucid[id]==struc) {
        break 
        }
   }
   if ( id>g_nid ) {
   print "ERROR in struc_element_type, struc unknown ",struc,element
   exit(10)
   }
   n=list[id,0,0]
   for (i=1;i<=n;i++) {
        if (list[id,i,2]==element) {
            return list[id,i,1]
        }
   }
   print "ERROR in struc_element_type, element unknown ",struc,element
   exit(10)
}



