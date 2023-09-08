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

    real(dp), dimension(size) :: vec2
    integer :: guessMode = 0, lineN = 1, lineN2
    real(dp) :: interpolationFactor = 0d0, tFactor = 1d0
    character(len = 400) :: guessFile, guessFile2=''
    namelist /guessInfo/ guessMode, guessFile, lineN, guessFile2, lineN2, interpolationFactor, tFactor

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
       vec(1) = vec(1) * tFactor
       vec(2:3) = vec(2:3)*tFactor

    else if (guessMode == 1) then
       call guessFile2Vec(size - 3, guessFile, lineN, vec(4:))

    else if (guessMode == 2) then !Interpolate mode
       call guessFile2Vec(size,  guessFile,  lineN,  vec(:))
       call guessFile2Vec(size, guessFile2, lineN2, vec2(:))
       vec(4:) = vec2(4:) + (vec2(4:) - vec(4:))*interpolationFactor
       vec(1) = (vec2(1) + (vec2(1) - vec(1))*interpolationFactor) * tFactor
       vec(2:3) = (vec2(2:3) + (vec2(2:3) - vec(2:3))*interpolationFactor) * tFactor
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
       read(12, *) vec_(:)
    end do
    close(12)
  end subroutine guessFile2Vec
end module PI_RPO_Guess
