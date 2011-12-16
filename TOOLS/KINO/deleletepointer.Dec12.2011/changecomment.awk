/comment *=.*ckino/{
  addbegin=0
  if (match($0,"BEGIN")) { addbegin=1}
  if (addbegin) { print "BEGIN{" }
  print " str=strftime(\"%b.%d.%Y\");"
  print " comment=\"ckino \" str \": \";"
  if (addbegin) { print "}" }
  next;
}
{print;}
