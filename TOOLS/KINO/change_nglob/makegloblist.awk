BEGIN{ IGNORECASE=1 }
/^ .* subroutine /{
  processname=$0
  showprocessname=0
  next
}
/^ .*nglob *\(/{
  if (showprocessname==0) {print ""; print FILENAME;print processname ; showprocessname=1}
  print
  next
}
/^ .*dglob *\(/{
  if (showprocessname==0) {print "" ;print FILENAME; print processname; showprocessname=1}
  print
  next
}
/^ .* function /{
  begin=0
  if ($1=="function") {
    begin=1
  }
  if ($2=="function") {
    if ($1=="real" || $1=="integer" || $1=="logical" || $1=="complex" || $1=="real(8)" ) {
       begin=1
    }
  }
  if ($3=="function") {
    if ($1=="double" && $2=="precision") { begin=1}
    if ($1=="double" && $2=="complex") { begin=1}
  }
  processname=$0
  showprocessname=0
   next
}
