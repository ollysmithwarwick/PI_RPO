!*************************************************************************
!  Example with one extra parameter T.
!  Can put parameter and constraint =0 when not required.
!- - - - - - - - - - - - - - -
!  Newton vector:
!    x(1:3)   = (/S, shift, phase/)
!    x(4: 6* N_x+ 3) = vector containing real and imag parts (for tilde quantities) of each quantity at each point
!
!  Extra constraints:
!    (F(x)-x). dx/dt = 0 .
!       no update along direction of trajectory.
!    
!    (F_(shift)(x) - x).dx/dX (X being the coordinate, x is the vector)
!       no update in the X direction (the solver is solving for a length translation anyway)
!
!    (F_(phase)(x) - x).dx/d(phi)
!       no update in a phase rotatin direction (solver solves for global phase transition)
!
!
!  Jacobian approximation:
!    dF(x_n)/dx . dx = (F(x_n+eps.dx)-F(x_n))/eps Used to quickly approximate dF/dx, made more accurate by reducing eps
!    
!
!  File guess.in and guess.out:
!     S		initial guess for shear
!     shift     initial guess for translation length
!     phase     initial guess for phase shift
!     n_period  number of boundary condition periods (1 bc period = 2*pi/LS)  
!     ndts	number of timesteps taken in one period - if set to -1 this will split the total time into 0.01 unit chunks,
!               though as the total time changes with each guess (since S changes -> T = n_period*2*pi/L*S changes) the size
!               of each step changes - need to keep an eye to make sure it doesn't get too big
!     info_out  Prints out the info integer in guess.out and is ignored in guess.in (so you can copy the out file to the in file without changing anything)
!               0 If Converged
!               1 If the code is running - this should never be printed out, since it should only print once it is complete
!               2 If nits > nits_max
!               3 If trust region gets too small
!               4 If there was a memory error
!     nits      Number of newton iterations used to reach convergence (or whatever other outcome)
!               An approximate measure of how far the guess was from the solution, but depends on many factors
!*************************************************************************
!*************************************************************************
 module POPI
   use parameters
   use PIRunModule
   use RPOSolve
   integer :: Nx
   real(dp) :: dx, L
   integer :: n_int
   logical :: fixedS=.false., zeroS=.false.


   type POPIRun      
      double precision, dimension(:), allocatable :: initGuess
      real(dp)         :: L
      integer          :: Nx
      integer          :: spf
      integer          :: n_int, info
      logical          :: fixedS=.false.
      real(dp)         :: S = 1d0
      real(dp)         :: dt = -1d0
      integer          :: ndts = -1
   end type POPIRun
   
 contains
   subroutine assignPOPIRun(POPIRunIn, POPIRunOut)
     ! Copies one POPI run to another. Probably a better way to do this with operator overloading
     type(POPIRun), intent(in) :: POPIRunIn
     type(POPIRun), intent(out) :: POPIRunOut

     if (allocated(POPIRunOut%initGuess)) then
        deallocate(POPIRunOut%initGuess)
     end if

     allocate(POPIRunOut%initGuess(size(POPIRunIn%initGuess)))

     POPIRunOut%n_int      = POPIRunIn%n_int
     POPIRunOut%initGuess  = POPIRunIn%initGuess
     POPIRunOut%info       = POPIRunIn%info
     POPIRunOut%fixedS     = POPIRunIn%fixedS
     POPIRunOut%S          = POPIRunIn%S     
     POPIRunOut%dt         = POPIRunIn%dt
     POPIRunOut%ndts       = POPIRunIn%ndts

   end subroutine assignPOPIRun

   subroutine file2POPIRun(fileName, POPIRunOut)
     ! Uses namelist file to set up a POPIRun
     use PI_RPO_Guess, only : getGuess
     character(len = 400), intent(in)  :: fileName
     type(POPIRun) , intent(out) :: POPIRunOut
     character(len = 400) :: guessInfoFile

     real(dp)         :: L
     integer          :: Nx
     integer          :: spf
     integer          :: n_int, info
     logical          :: fixedS=.false.    ! If false, first element of init solves for S
     real(dp)         :: dt = -1d0
     integer          :: ndts = -1
     logical          :: overwriteParams = .false., overwriteS = .false. ! Will only overwite params in guess if true
     real(dp)         :: transGuess=0d0, phaseGuess=0d0, periodGuess =-1d0 ! Only changes params if overwrite is true
     real(dp)         :: S = 1d0           ! Only required if fixedS is true - otherwise uses first element in init
     real(dp), dimension(:), allocatable :: vec
          
     namelist /POPIParams/ L, Nx, spf, n_int, info, fixedS, S, dt, ndts, &
          guessInfoFile, overwriteParams, transGuess, PhaseGuess, periodGuess, &
          overwriteS

     write(*,*) 'POPI: Reading '
     
     open(12, file = fileName)
     read(12, nml = POPIParams)
     close(12)
     
     POPIRunOut%L           = L
     POPIRunOut%Nx          = Nx
     POPIRunOut%spf         = spf
     POPIRunOut%n_int       = n_int
     POPIRunOut%info        = info
     POPIRunOut%fixedS      = fixedS
     POPIRunOut%S           = S
     POPIRunOut%ndts        = ndts
     POPIRunOut%dt          = dt

     allocate(vec(6*Nx))
     allocate(POPIRunOut%initGuess(6*Nx + 3))
     vec = 0d0
     POPIRunOut%initGuess = 0d0
     call getGuess(6*Nx+3, guessInfoFile, vec)
     !vec(4:Nx+3) = 0d0
     POPIRunOut%initGuess(:)   = vec
     
     ! If in fixedS mode, parameter 1 contains period, otherwise it is shear
     if (fixedS) then
        if (overwriteParams .or. POPIRunOut%initGuess(1)==0d0) then 
           if (periodGuess <0) then
              POPIRunOut%initGuess(1) = 2d0 * pi * n_int / (L*S)
           else
              POPIRunOut%initGuess(1) = periodGuess
           end if          
        end if
     else
        if (overwriteS) then
           POPIRunOut%initGuess(1) = S
        end if
     end if


     if (overwriteParams) then
        POPIRunOut%initGuess(2) = transGuess
        POPIRunOut%initGuess(3) = PhaseGuess
     end if
        
   end subroutine file2POPIRun

   subroutine POPIRun2file(POPIRunOut, fileName, initInfoFileIn)

     
     character(len = 100), intent(in)  :: fileName
     character(len = 100), intent(in), optional :: initInfoFileIn
     character(len = 100) :: initInfoFile = ''
     type(POPIRun) , intent(out) :: POPIRunOut
     
     real(dp)         :: L
     integer          :: Nx
     integer          :: spf
     integer          :: n_int, info
     logical          :: fixedS=.false.
     real(dp)         :: S = 1d0
     real(dp)         :: transGuess, phaseGuess, periodGuess = -1d0
     real(dp)         :: dt
     integer          :: ndts

          
     namelist /POPIParams/ L, Nx, spf, n_int, info, fixedS, S, ndts, dt, initInfoFile, transGuess, phaseGuess, periodGuess
     
     
     L  =   POPIRunOut%L
     Nx  =   POPIRunOut%Nx
     spf  =   POPIRunOut%spf
     n_int =    POPIRunOut%n_int
     fixedS =    POPIRunOut%fixedS
     ndts = POPIRunOut%ndts
     dt = POPIRunOut%dt

     if (fixedS) then
        S    = POPIRunOut%S
        periodGuess = POPIRunOut%initGuess(1)
     else
        S = POPIRunOut%initGuess(1)
        periodGuess = -1d0
     end if

     if (present(initInfoFileIn)) then
        initInfoFile = initInfoFileIn
     end if
     transGuess = POPIRunOut%initGuess(2)
     phaseGuess = POPIRunOut%initGuess(3)

     open(12, file = fileName)
     write(12, nml = POPIParams)
     close(12)
   end subroutine POPIRun2file

   subroutine POPI_main(rpoParams, POPIParams, states_, guesses_, infoOut, J, evals, evecs)
     !*************************************************************************
     use RPOSolve, only : RPONewtonSetTimestep, RPONewtonSetDotprod, RPONewtonSetSaveGuess, stepOrbit, RPONewtonVecSetup, RPONewtonParams, RPONewtonSetup, RPOSolveNewton, dotprod
     implicit none

     double precision :: d
     integer :: i1, i2, i3
     type(RPONewtonParams) :: rpoParams
     type(POPIRun) :: POPIParams

     real(dp), dimension(:), allocatable :: initGuess
     real(dp), dimension(1:,1:) :: states_
     real(dp), dimension(:,:), allocatable :: guessesTmp_
     real(dp), dimension(:,:), allocatable :: guesses_
     real(dp), dimension(:,:), optional :: J
     complex(dp), dimension(:), optional :: eVals
     complex(dp), dimension(:,:), optional :: eVecs

     complex(dp), dimension(:), allocatable :: eVals_
     complex(dp), dimension(:,:), allocatable :: eVecs_

!     real(dp), dimension(:,:), allocatable, intent(out), optional :: Jac_
     integer, optional :: infoOut
     logical :: jacOut = .false., evalOut = .false.
     integer :: info, size, nframes

     integer, dimension(2) :: shp = (/0, 0/)
     logical :: saveGuessFlag = .true.
     integer :: nits = 0
     double complex, dimension(:), allocatable :: final_f
     external :: getrhs, saveorbit

     !   character(4) :: cnum

     states_ = 0d0
!     Jac_ = 0d0

     L   = POPIParams%L
     Nx  = POPIParams%Nx
     n_int = POPIParams%n_int

     allocate(initGuess(6*Nx))
     initGuess   = POPIParams%initGuess
     n_int       = POPIParams%n_int
     fixedS      = POPIParams%fixedS


     allocate(guessesTmp_(6*Nx + 3, rpoParams%nits + 1))
     
     ! If S=0 then the boundary conditions blow up - fix S at zero and use fixedS mode
     if (fixedS) then
        if (POPIParams%S==0d0) then
           zeroS = .true.
        end if
        if (rpoParams%ndts == -1) then
          rpoParams%ndts = 1000
        end if
        if (initGuess(1)==0.0) then !! Require non-zero period so if fixedS mode need period != 0
           error stop 'Need non zero period!'
        end if
     end if

     size  = 6 * Nx + 3

     if (present(J)) then
        JacOut = .true.
     end  if

     if(present(evals) .or. present(evecs)) then
        evalOut = .true.
        allocate(evals_(size - 3))
        allocate(evecs_(size - 3, size - 3))
     end if


     write(*,*) 'evals_', shape(evals_)
     write(*,*) 'evecs_', shape(evecs_)
     evals_=0d0
     evecs_=0d0

 !    if (present(Jac_)) then
 !       JacOut = .true.
 !       allocate(J(n-3, n-3))
 !       J = 0d0
 !       allocate(Jac_(n-3, n-3))
 !    end if

 !    write(*,*) 'allocated Jac_'

     guessesTmp_ = 0d0
     shp = shape(states_)
     nframes = shp(2) - 2

     !   call writeplotinfo                   ! Writes data to plotinfo.in for plotting later

!*****
!   LOAD STATE TO  new_x(3:) - this is your initial guess. The data structure can be determined from the allocation here, it's inverse is simple too
      !*****
!     call makepos(N_x)

     ! scale params by norm of state - i.e all errors are relative to the INITIAL GUESS not the final state
     
     !call store_traj(n_int*2*pi*L_inv/S, ndts) ! Stores the initial guess' trajectory in out.dat - if an error occurs or the code is exited early this can be helpful
     
     call storeTrajectory(initGuess, states_)

     if (initGuess(2) == 0) then
        write(*,*) 'POPI_Main: Calling getlshift'
        initGuess(2) = -bestshift(states_(:,1), states_(:,nframes+2), Nx)*L/Nx
        write(*,*) 'bestshift = ', initGuess(2)
        call storeTrajectory(initGuess, states_)
     end if
     
    if (initGuess(3) == 0) then
       write(*,*) 'POPI_Main: Calling getphase'
       initGuess(3) = bestphase(states_(:,1), states_(:,nframes+2), Nx)
       write(*,*) 'bestshift = ', initGuess(3)
       call storeTrajectory(initGuess, states_)
    end if


     write(*,*) 'Initial guess trajectory stored.'
!     info = 0
     if (rpoParams%nits > 0) then
        write(*,*) 'newton_PI: Calling newtonhook'
        call RPONewtonSetTimestep( timestepPI )
        call RPONewtonSetDotprod( dotprodPI )
        call RPONewtonSetSaveGuess( saveGuessPI )
        call RPONewtonVecSetup( (/ 3, 1, 2 /), (/ .false., .false., .false. /), info)
        call RPONewtonSetup(rpoParams)
        call RPOSolveNewton(size, initGuess, info)

!        write(*,*) 'NDTS = ', ndts
!        write(*,*) '   S = ', new_x(1)
        !     guesses_ = guesses
        !     deallocate(guesses)
        
        allocate(guesses_(size, nits+1))
        guesses_ = guessesTmp_(:,1: nits+1)
!        write(*,*) 'POPI_Main: Final guess  = ', guessesTmp_(:,nits+1)
        call storeTrajectory(guesses_(:,nits+1), states_)
     end if
     deallocate(guessesTmp_)

!      if (fixedS) then
!         rpoParams%SGuess = SGuess !Unchanged
!         PIParams%dt = new_x(1)/dble(ndts)
!      else
!         rpoParams%SGuess = new_x(1)
!         PIParams_%dt = n_int*2*pi /(new_x(1) * PIParams%L * dble(ndts)) ! n_int is number of intervals
!      end if
!      rpoParams%ShiftGuess = new_x(2)
!      rpoParams%PhaseGuess = new_x(3)
!      rpoParams%initGuess = new_x(4:)
!      rpoParams%info = info

!      PIParams_%S = rpoParams%SGuess

!      PIParams_%init = new_x(4:)
     
!      PIParams_%N_t = ndts

!      if (present(infoOut)) then
!         infoOut = info
!      end if

     if (JacOut .or. evalOut) then
        write(*,*) 'Getting full matrix'
        call RPOGetMatrix(size,3,J)
        write(*,*) 'Full matrix found'
        if (evalOut) then
           write(*,*) 'Getting Eigenvalues'
           call RPOGetEigenvalues(size - 3, J, eVals_, eVecs_)
           write(*,*) 'Getting Eigenvalues'
        end if
     end if
     
     if (present(eVals)) then
        eVals = eVals_
     end if
     if (present(eVecs)) then
        eVecs = eVecs_
     end if
     
!      if (present(fluxs)) then
!         call GetStats(states_, amps, fluxs, avgAmps, avgFluxs, pkAmps, pkFluxs)
!      end if

!      if (allocated(guesses)) then
! !        guesses_ = guesses
!         deallocate(guesses)
!      end if
!      deallocate(new_x)
!      deallocate(new_fx)
!      deallocate(initGuess)
!      write(*,*) 'newton_PI: Done'
     
   contains
     function timestepPI(n, in_, ndts, dt)
       use PIModule, only : PI_main 
       use PIRunModule
       implicit none
       real(dp), dimension(:), intent(in) :: in_
       integer, intent(in) :: ndts, n
       integer :: ndts_
       real(dp), intent(in) :: dt
       real(dp) :: dt_, duration
       real(dp), dimension(:), allocatable :: timestepPI
       real(dp), dimension(:,:), allocatable :: states
       
       type(PIRun) :: run0
       integer :: j
       allocate(timestepPI(n))
       allocate(states(n-3,2))
       timestepPI = 0d0

       
       if (ndts == 1) then
          duration = dt
          ndts_ = 1
       else if (ndts == 0) then
          duration = 0d0
          ndts_ = 1
       else if (fixedS) then
          duration = in_(1)
          ndts_ = ndts
       else
          duration = 2d0*pi*n_int/(in_(1)*POPIParams%L)
          ndts_ = ndts
       end if
       
       dt_ = duration/ndts_

      ! WRITE(*,*) 'TimestepPI: n    =', n
      ! WRITE(*,*) 'TimestepPI: ndts =', ndts
      ! WRITE(*,*) 'TimestepPI: dt   =', dt_
      ! WRITE(*,*) 'TimestepPI: in   =', in_
      ! WRITE(*,*) 'TimestepPI: Duration = ', duration

       run0%shift = in_(2)
       run0%phase = in_(3)
       run0%N_t   = ndts_
       run0%spf   = ndts_ + 1
       run0%dt    = dt_
       run0%L     = POPIParams%L
       run0%N_x   = POPIParams%Nx
       if (fixedS) then
          run0%S = POPIParams%S
       else
          run0%S = in_(1)
       end if
       run0%init = in_(4:)
       call PI_main(run0, states)
!       WRITE(*,*)  'TimestepPI: Returned'
       timestepPI = 0d0
       timestepPI(4:) = states(:,2)
       
     end function timestepPI

     double precision function dotprodPI(n_,a,b)
       implicit none
       integer,          intent(in) :: n_
       double precision,       dimension(:), intent(in) ::  a
       double precision,       dimension(:), intent(in) ::  b
       integer :: n1
       
       n1 = 1
       if (n_==-1) n1 = 4
       dotprodPI = dble(dot_product(a(n1:),b(n1:)))       
     end function dotprodPI

     subroutine saveGuessPI()
       !Called once per newton iteration
       use newton, only: new_x, new_nits
       implicit none
       if (saveGuessFlag) then
         guessesTmp_(:, new_nits + 1) = new_x(:)
       end if
       write(*,*) 'shape(new_x) = ', shape(new_x)
       write(*,*) 'shape(new_x) = ', shape(guessesTmp_)
       nits = new_nits
     end subroutine saveGuessPI


     subroutine storeTrajectory(x_, states_)
       ! use runtype
       use parameters
       !use orbit, only: n_int
       use PIModule, only: PI_main

       real(dp) :: dt_
       real(dp), dimension(1:) :: x_
       real(dp), intent(out), dimension(1:,1:) :: states_
       type(PIRun) :: runIn
       write(*,*) 'Calculating dt_'

       dt_ = n_int*2.0_dp*pi/(x_(1) * L * dble(RPOParams%ndts)) ! n_int is number of Boundary Condition intervals

       write(*,*) 'Calling AssignRun. x(1) = ', x_(1) 
       write(*,*) 'Calling AssignRun. L= ', L
       write(*,*) 'Calling AssignRun. dt = ', dt_
       write(*,*) 'Calling AssignRun. n_int = ', n_int
       write(*,*) 'Calling AssignRun. pi = ', pi
       write(*,*) 'Calling AssignRun. spf = ', POPIParams%spf
       !call AssignRun(PIParams, runIn)
       runIn%init = x_(4:)
       if (.not. fixedS) then
          runIn%S    = x_(1)
          runIn%dt   = dt_
       else
          runIn%dt=x_(1)/RPOParams%ndts
          runIn%S=POPIParams%S
       end if

       runIn%L    = L
       runIn%N_t  = RPOParams%ndts
       runIn%N_x   = POPIParams%Nx
       runIn%spf   = POPIParams%spf
       runIn%nonlin = 1
       runIn%shift = x_(2)
       runIn%phase = x_(3)

       write(*,*) 'Calling PI_main. dt = ', runIn%dt
       call PI_main(runIn, states_)
       !     write(*,*) 'states = ' , states_

       write(*,*) 'storeTraj complete'
     end subroutine storeTrajectory

     ! function getLShift(x)
     !   implicit none
     !   real(dp), intent(out) :: getlshift
     !   real(dp), dimension(n) :: y_
     !   y_(:) = 0d0
     !   call steporbit(n, ndts, 0d0, 0d0, PIParams, new_x, y_)
     !   new_x(2) = -bestshift(new_x(4:), y_(4:))*PIParams%L/N_x
     ! end subroutine getLShift
     
     ! subroutine getphase
     !   implicit none
     !   real(dp), dimension(n) :: y_
     !   call steporbit(n, ndts, new_x(2), 0d0, PIParams, new_x, y_)
     !   new_x(3) = bestphase(new_x(4:), y_(4:))
     ! end subroutine getphase
   !*************************************************************************
   end subroutine POPI_main
   !*************************************************************************

!-------------------------------------------------------------------------
!  determine best guess for Lshift
!-------------------------------------------------------------------------


   function bestShift(stateStart, stateEnd, N_x)
     real(dp), dimension(1:), intent(in) :: stateStart, stateEnd
     integer :: N_x
     integer:: bestShift
     complex(dp), dimension(N_x) :: nstart, phistart, nend, phiend
     complex(dp), dimension(N_x) :: fluxStart, fluxEnd
     double precision :: prod(N_x)
     integer :: i1, i2
     real(dp) :: s
     do i1 = 1, N_x
        phistart(i1) = stateStart(2*i1 + 2*N_x - 1) + im * stateStart(2*i1 + 2*N_x)
        nstart(i1)   = stateStart(2*i1 + 4*N_x - 1) + im * stateStart(2*i1 + 4*N_x)
        phiend(i1)   =   stateEnd(2*i1 + 2*N_x - 1) + im *   stateEnd(2*i1 + 2*N_x)
        nend(i1)     =   stateEnd(2*i1 + 4*N_x - 1) + im *   stateEnd(2*i1 + 4*N_x)
     end do
     
     fluxStart(:) = abs(phiStart(:)*conjg(nStart(:)) - nstart(:)*conjg(phiStart(:)))
     fluxEnd(:)   = abs(  phiEnd(:)*conjg(  nEnd(:)) -   nEnd(:)*conjg(  phiEnd(:)))
     prod(:) = 0.0
     do i1 = 1, N_x
        do i2 = 1, N_x
           s = (fluxStart(modulo(i2+i1-2, N_x)+1) - fluxEnd(i2))
           prod(i1) = prod(i1) + s * s
        end do
     end do
     bestshift = minloc(prod(:), dim = 1)
     bestshift = modulo((bestshift - 1 + N_x/2), N_x) - N_x/2
!     bestshift = bestshift-1
   end function bestShift
   
   function bestPhase(stateStart, stateEnd, N_x)
     use parameters
     real(dp), dimension(1:), intent(in) :: stateStart, stateEnd
     real(dp) :: bestPhase
     integer :: N_x
     complex(dp), dimension(N_x) :: nstart, phistart, nend, phiend
     complex(dp), dimension(N_x) :: fluxStart, fluxEnd
     double complex :: s
     double precision :: prod(6283)
     integer :: i1, i2
     
     do i1 = 1, N_x
        phistart(i1) = stateStart(2*i1 + 2*N_x - 1) + im * stateStart(2*i1 + 2*N_x)
        nstart(i1)   = stateStart(2*i1 + 4*N_x - 1) + im * stateStart(2*i1 + 4*N_x)
        phiend(i1)   =   stateEnd(2*i1 + 2*N_x - 1) + im *   stateEnd(2*i1 + 2*N_x)
        nend(i1)     =   stateEnd(2*i1 + 4*N_x - 1) + im *   stateEnd(2*i1 + 4*N_x)
     end do

     prod(:) = 0d0
     do i1 = 1, 6283
        do i2 = 1, N_x
           s = -phiStart(i2) + phiEnd(i2)*exp(im*(i1-1)*0.001d0)
           prod(i1) = prod(i1) + dble(s*conjg(s))
           s = -nStart(i2) + nEnd(i2)*exp(im*(i1-1)*0.001d0)
           prod(i1) = prod(i1) + dble(s*conjg(s))
        end do
     end do
!     write(*,*) 'bestshift: prod =  ', prod

     bestPhase = (minloc(prod(:), dim = 1) - 1) * 0.001d0
    end function bestPhase




   
   ! subroutine file2POPIRun(fileName, POPIRunOut)
   !   character(len = 100), intent(in)  :: fileName
   !   type(POPIRun) , intent(out) :: POPIRunOut
     
   !   integer  :: mgmres, nits, n_int, info
   !   real(dp) :: rel_err
   !   real(dp) :: del, mndl, mxdl
   !   real(dp) :: gtol, epsJ
   !   real(dp) :: SGuess, ShiftGuess = 0.0_dp, PhaseGuess = 0.0_dp, periodGuess=0.0_dp
   !   logical  :: fixedS = .false.
     
   !   namelist /POPIParams/ mgmres, nits, n_int, rel_err, del, mndl, mxdl, gtol, epsJ, SGuess, ShiftGuess, PhaseGuess, info, fixedS, periodGuess
     
   !   open(12, file = fileName)
   !   read(12, nml = POPIParams)
   !   close(12)
     
   !   POPIRunOut%mgmres     = mgmres
   !   POPIRunOut%nits       = nits
   !   POPIRunOut%n_int      = n_int
   !   POPIRunOut%rel_err    = rel_err
   !   POPIRunOut%del        = del
   !   POPIRunOut%mndl       = mndl
   !   POPIRunOut%mxdl       = mxdl
   !   POPIRunOut%gtol       = gtol
   !   POPIRunOut%epsJ       = epsJ
   !   POPIRunOut%SGuess     = SGuess
   !   POPIRunOut%ShiftGuess = ShiftGuess
   !   POPIRunOut%PhaseGuess = PhaseGuess
   !   POPIRunOut%info       = info
   !   POPIRunOut%fixedS     = fixedS
   !   POPIRunOut%periodGuess=periodGuess
   ! end subroutine file2POPIRun

   ! subroutine POPIRun2file(POPIRunOut, fileName)
   !   character(len = 100)  :: fileName
   !   type(POPIRun) , intent(in) :: POPIRunOut
     
   !   integer  :: mgmres, nits, n_int, info
   !   real(dp) :: rel_err
   !   real(dp) :: del, mndl, mxdl
   !   real(dp) :: gtol, epsJ
   !   real(dp) :: SGuess, ShiftGuess = 0.0_dp, PhaseGuess = 0.0_dp, PeriodGuess=0.0_dp
   !   logical  :: fixedS = .false.
     
   !   namelist /POPIParams/ mgmres, nits, n_int, rel_err, del, mndl, mxdl, gtol, epsJ, SGuess, ShiftGuess, PhaseGuess, info, fixedS, periodGuess
     
   !   mgmres      = POPIRunOut%mgmres
   !   nits        = POPIRunOut%nits  
   !   n_int       = POPIRunOut%n_int 
   !   rel_err     = POPIRunOut%rel_err
   !   del         = POPIRunOut%del    
   !   mndl        = POPIRunOut%mndl   
   !   mxdl        = POPIRunOut%mxdl   
   !   gtol        = POPIRunOut%gtol   
   !   epsJ        = POPIRunOut%epsJ   
   !   SGuess      = POPIRunOut%SGuess 
   !   ShiftGuess  = POPIRunOut%ShiftGuess
   !   PhaseGuess  = POPIRunOut%PhaseGuess
   !   info        = POPIRunOut%info
   !   fixedS      = POPIRunOut%fixedS
   !   periodGuess = POPIRunOut%periodGuess

   !   open(12, file = fileName)
   !   write(12, nml = POPIParams)
   !   close(12)
     

   ! end subroutine POPIRun2file

 end module POPI
