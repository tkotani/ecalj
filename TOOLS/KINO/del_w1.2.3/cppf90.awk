/#if/{ flag=do_macro($0); if (flag==0) { print; } else if ( }
/#else/{ }
/#elif/{ flag=do_macro($0) }
/#endif/{ flag=do_macro($0) }

func do_macro(line)
{
    a=match(line,"F90")
    b=match(line,"!")

    if (a==0) {
       return  0
    }

    if (a>0 && b==0) {
	return 1
    }
    if (a>0 && b>0) {
	return -1
    }

}

