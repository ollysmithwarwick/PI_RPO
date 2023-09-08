program main
  use POPI, only : file2POPIRun, POPIRun, POPI_main
  use PI_IO
  use RPOSolve, only : RPONewtonParams, file2RPORun
  use parameters

  implicit none

  type(RPONewtonParams) :: rpoParams
  type(POPIRun)         :: POPIParams
  real(dp), dimension(:,:), allocatable :: states, guesses, fluxs, amps, avgFluxs, avgAmps, pkFluxs, pkAmps, J
  real(dp), dimension(:), allocatable :: tVals
  complex(dp), dimension(:,:), allocatable :: eVecs
  complex(dp), dimension(:), allocatable :: eVals
  complex(dp), dimension(:,:), allocatable :: eVals2D
  character(len = 400) :: SetupFile = 'POPISetup.in', rpoFile, POPIFile, StatesFile, TFile, OutDir, GuessFile, JFile='', eValFile='', eVecFile=''
  integer, dimension(2) :: shp
  integer :: info
  integer :: nframes, nits
  integer :: jt, jf
  logical :: jmode = .false.
  namelist /POPISetup/ rpoFile, POPIFile, OutDir, StatesFile, TFile, GuessFile, JFile, eValFile, eVecFile

  ! Setup file contains names of files containing PI model params and RPO solver params
  write(*,*) 'Reading setup file'
  open(21, file = SetupFile)
  read(21, nml = POPISetup)
  close(21)


  write(*,*) 'Setup read done'
  write(*,*) rpoFile
  write(*,*) POPIFile

  ! Convert POPIFile into POPIRun type (contains all info about PI specific implementation)
  call file2POPIRun(POPIFile,POPIParams)
  call file2RPORun(RPOFile, RPOParams)


  ! **************************************************************************************!
  !                    MOVE INTO POPI.f90?                                               !

  ! POPI Always uses the fixed NDTS mode of RPOSolve. Can find an ndts that gives close
  ! to target timestep or just select ndts. POPI Params will overwrite RPOParams

  if (POPIParams%ndts > 0) then    ! Overwrite RPO ndts if needed
     RPOParams%ndts = POPIParams%ndts
  else if (POPIParams%dt > 0) then
     RPOParams%ndts = floor(POPIParams%n_int *2d0 * pi/(POPIParams%L * POPIParams%S * POPIParams%dt))  ! Determine ndts from target dt (approx)
     POPIParams%ndts = floor(POPIParams%n_int *2d0 * pi/(POPIParams%L * POPIParams%S * POPIParams%dt))  ! Determine ndts from target dt (approx)
  else
     ! Stick with using ndts from RPO file (this defaults to 1000 in RPOSolver if ndts < 0)
  end if
  ! *************************************************************************************!

  write(*,*) 'POPI: Nx = ', POPIParams%Nx
  write(*,*) 'POPI: spf = ', POPIParams%spf
  nframes = int(POPIParams%ndts/POPIParams%spf)
  write(*,*) 'POPI: nf = ', nframes
  allocate(states(6*POPIParams%Nx, nframes + 2))
  allocate(J(6*POPIParams%Nx,6*POPIParams%Nx))
  allocate(tVals(nframes + 1))
  allocate(eVecs(6*POPIParams%Nx,6*POPIParams%Nx))
  allocate(eVals(6*POPIParams%Nx))
  allocate(eVals2D(1,6*POPIParams%Nx))


  ! Convert POPIFile into POPIRun type (contains all info about PI specific implementation)
!  call file2POPIRun(POPIFile,POPIParams)
  write(*,*) RPOParams
  if (.not. (JFile == '')) then
     jmode = .true.
     call POPI_main(rpoParams, POPIParams, states, guesses, info, J, eVals, eVecs)
  else
     call POPI_main(rpoParams, POPIParams, states, guesses, info)
  end if

!  write(*,*) evecs
  evals2D(1,:) = evals
  shp = shape(guesses)
  nits = shp(2)-1
  jf = 0
  tVals(0) = 0d0
  do jt = 1, POPIParams%ndts
    if (mod(jt,POPIParams%spf) == 0) then
       jf = jf + 1
       tVals(jf + 1) = jt * POPIParams%dt
    end if
  end do

  

  call system('mkdir -p '//trim(adjustl(OutDir)))
  call writeFile( trim(adjustl(OutDir))//trim(adjustl(StatesFile)), states)
  if (jmode) then
     call writeFile( trim(adjustl(OutDir))//trim(adjustl(JFile)), J)
     call writeFileComplex( trim(adjustl(OutDir))//trim(adjustl(eVecFile)), eVecs)
     call writeFileComplex(trim(adjustl(OutDir))//trim(adjustl(eValFile)), eVals2D)
  end if
  call writeFile1D( trim(adjustl(OutDir))//trim(adjustl(GuessFile)), guesses(:,nits+1))
  call writeFile1D( trim(adjustl(OutDir))//trim(adjustl(TFile))     ,  tvals)
  call system('cp '//trim(adjustl(POPIFile))//' '//trim(adjustl(OutDir)))
  call system('cp '//trim(adjustl(rpoFile))//' '//trim(adjustl(OutDir)))

end program main
