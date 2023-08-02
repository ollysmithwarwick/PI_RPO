module PI_RPO_Guess
  use parameters
  use, intrinsic :: iso_c_binding
  implicit none 
  character(len =20) :: mainfmt=''

contains

  subroutine getGuess(size, guessInfoFile, vec)
    integer, intent(in) :: size
    character (len =400), intent(in) :: guessInfoFile
    real(dp), dimension(size) :: vec

    integer :: guessMode = 0, lineN = 1
    character(len = 400) :: guessFile
    namelist /guessInfo/ guessMode, guessFile, lineN

    ! guessMode
    !  0 Read in full guess from a file (i.e 6Nx+3 size)
    !  1 Read in only positional guess from file (i.e 6*Nx+3)
    !      but set parameter guesses to 0 (to be dealt with elsewhere)
    !      Useful for reading in bisection output for example

    open(21, file = guessInfoFile)
    read(21, nml = guessInfo)
    close(21)
    vec = 0d0
    if (guessMode == 0) then
       call guessFile2Vec(size, guessFile, lineN, vec)
    else if (guessMode == 1) then
       call guessFile2Vec(size - 3, guessFile, lineN, vec(4:))
    end if
  end subroutine getGuess
    
  subroutine guessFile2Vec(size, fileName, lineN, vec_)
    integer, intent(in) :: size
    character(len = 100), intent(in) :: fileName
    real(dp), dimension(size) :: vec_
    integer, intent(in) :: lineN
    integer :: i1

    write(mainfmt, '(A1,I8,A11)') '(', size,'(E30.15E3))'
    open(12, file = fileName)
    do i1 = 1, lineN
       read(12, fmt = mainfmt) vec_(:)
    end do
    close(12)
  end subroutine guessFile2Vec
end module PI_RPO_Guess
