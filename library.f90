
subroutine CON_stop(instr)

  character (len=*), intent(in) :: instr

  write(*,*) instr
  stop
  
end subroutine CON_stop
