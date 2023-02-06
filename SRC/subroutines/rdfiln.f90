module m_rdfiln! preprocessor. ctrl.* is conveter to ctrl_preprocessd.*
  public M_rdfiln_init
  private
contains
  subroutine m_rdfiln_init() !Read recln,recrd,nrecs,mxrecs for m_lmfinit_init
    !     -vfoobar replaced simplified contents of ctrl into recrd
    use m_lgunit,only:stdo
    use m_ext,only: sname
    integer:: i
    character(512):: aaa,cmdl,argv
    logical:: fileexist
    inquire(file='ctrl.'//trim(sname),exist=fileexist)
    if( .NOT. fileexist) call rx("No ctrl file found! ctrl."//trim(sname))
    aaa=''
    do i = 1, iargc()
       call getarg( i, argv )
       aaa=trim(aaa)//' '//trim(argv)
    enddo
    cmdl='ctrl2ctrlp.py '//trim(aaa)//'<ctrl.'//trim(sname)//' >ctrlp.'//trim(sname) !core is in python
    write(stdo,*)trim(cmdl)
    call system(cmdl)
  end subroutine m_rdfiln_init
endmodule m_rdfiln
