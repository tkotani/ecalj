/^     [^ ]/ {
   if (length(str)>0) {
   printf("      %s",substr(str,7));
   str="";
   }
}
{
   if (length(str)>0) { printf("      %s\n",substr(str,7));}
   str=$0;
}
END{ if (length(str)>0) { printf("      %s\n",substr(str,7));}
}
/^[^ ]/{ print str; str="";  }
