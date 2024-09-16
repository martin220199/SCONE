module criticalKineticPhysicsPackage_class

  use numPrecision
  use universalVariables
  use endfConstants
  use genericProcedures,              only : fatalError, printFishLineR, numToChar, rotateVector
  use hashFunctions_func,             only : FNV_1
  use dictionary_class,               only : dictionary
  use outputFile_class,               only : outputFile

  ! Timers
  use timer_mod,                      only : registerTimer, timerStart, timerStop, &
                                             timerTime, timerReset, secToChar

  ! Particle classes and Random number generator
  use particle_class,                 only : particle, P_NEUTRON, particleState
  use particleDungeon_class,          only : particleDungeon
  use RNG_class,                      only : RNG

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry
  use geometry_inter,                 only : geometry
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx, &
                                             gr_fieldIdx => fieldIdx, gr_fieldPtr => fieldPtr
  use geometryFactory_func,           only : new_geometry, killGeom

  ! Fields
  use field_inter,                    only : field
  use uniFissSitesField_class,        only : uniFissSitesField, uniFissSitesField_TptrCast
  use fieldFactory_func,              only : new_field

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat           => nMat
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_activate    => activate ,&
                                             ndReg_display     => display, &
                                             ndReg_kill        => kill, &
                                             ndReg_get         => get ,&
                                             ndReg_getMatNames => getMatNames
  use nuclearDatabase_inter,          only : nuclearDatabase

  ! Sources
  use source_inter,                   only : source
  use sourceFactory_func,             only : new_source

  ! Operators
  use collisionOperator_class,        only : collisionOperator
  use transportOperator_inter,        only : transportOperator

  ! Tallies
  use tallyCodes
  use tallyAdmin_class,               only : tallyAdmin
  use tallyResult_class,              only : tallyResult
  use keffAnalogClerk_class,          only : keffResult

  ! Factories
  use transportOperatorFactory_func,  only : new_transportOperator

  ! Visualisation
  use visualiser_class,               only : visualiser

  implicit none
  private

  !!
  !! Physics Package for eigenvalue calculations
  !!
  type, public,extends(physicsPackage) :: criticalKineticPhysicsPackage
    private
    ! Building blocks
    class(nuclearDatabase), pointer        :: nucData       => null()
    class(geometry), pointer               :: geom          => null()
    integer(shortInt)                      :: geomIdx       = 0
    type(collisionOperator)                :: collOpCritical, collOpKinetic
    class(transportOperator), allocatable  :: transOpCritical, transOpKinetic
    class(source), allocatable             :: initSource
    class(RNG), pointer                    :: pRNG          => null()
    type(tallyAdmin),pointer               :: inactiveTally => null()
    type(tallyAdmin),pointer               :: activeTally   => null()
    type(tallyAdmin),pointer               :: inactiveAtch  => null()
    type(tallyAdmin),pointer               :: activeAtch    => null()
    class(uniFissSitesField),pointer       :: ufsField      => null()
    type(tallyAdmin),pointer               :: tally         => null()

    ! Settings criticality
    integer(shortInt)  :: N_inactive
    integer(shortInt)  :: N_active
    integer(shortInt)  :: pop
    character(pathLen) :: outputFileCritical
    character(nameLen) :: outputFormatCritical
    integer(shortInt)  :: printSource = 0
    integer(shortInt)  :: particleType
    real(defReal)      :: keff_0
    integer(shortInt)  :: bufferSize
    logical(defBool)   :: UFS = .false.
    integer(shortInt)  :: N_stepBacks

    ! Calculation components criticality
    type(particleDungeon), pointer :: thisCycle    => null()
    type(particleDungeon), pointer :: nextCycle    => null()
    type(particleDungeon), pointer :: temp_dungeon => null()

    ! Settings kinetic
    character(pathLen) :: outputFileKinetic
    character(nameLen) :: outputFormatKinetic
    integer(shortInt)  :: N_cycles
    integer(shortInt)  :: N_timeBins
    real(defReal)      :: timeIncrement
    logical(defBool)   :: criticalWgtNorm
    logical(defBool)   :: useCombing
    logical(defBool)   :: usePrecursors
    logical(defBool)   :: useForcedPrecursorDecay
    integer(shortInt)  :: bufferShift

    real(defReal) :: minWgt = 0.25
    real(defReal) :: maxWgt = 1.25
    real(defReal) :: avWgt = 0.5

    ! Calculation components kinetic
    type(particleDungeon), pointer, dimension(:) :: currentTime       => null()
    type(particleDungeon), pointer, dimension(:) :: nextTime          => null()
    type(particleDungeon), pointer, dimension(:) :: tempTime          => null()
    type(particleDungeon), pointer, dimension(:) :: precursorDungeons => null() 

    ! Kinetic geometry
    class(dictionary), pointer, dimension(:)     :: kineticGeomDict  => null()
    integer(shortInt), dimension(:), allocatable :: kineticGeomTimeSteps

    ! Temporary storage of dictionary for memory efficiency
    class(dictionary), pointer                   :: dict             => null()

    ! Timer bins
    integer(shortInt) :: timerMain
    real (defReal)    :: time_transport = 0.0
    real (defReal)    :: CPU_time_start
    real (defReal)    :: CPU_time_end

  contains
    procedure :: init
    procedure :: run
    procedure :: printSettingsCritical
    procedure :: generateInitialState
    procedure :: cyclesCritical
    procedure :: collectResultsCritical

    procedure :: initKinetic
    procedure :: printSettingsKinetic
    procedure :: coupleKinetic
    procedure :: normCriticalWgts
    procedure :: cyclesKinetic
    procedure :: collectResultsKinetic


    procedure :: kill
    procedure :: russianRoulette
  end type criticalKineticPhysicsPackage

contains

  !!
  !! Initialise from individual components and dictionaries for inactive and active tally
  !!
  subroutine init(self, dict)
    class(criticalKineticPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)                    :: dict
    class(dictionary),pointer                           :: tempDict
    character(:), allocatable                           :: string
    type(dictionary)                                    :: locDict1, locDict2
    integer(shortInt)                                   :: seed_temp
    integer(longInt)                                    :: seed
    character(10)                                       :: time
    character(8)                                        :: date
    character(nameLen)                                  :: nucData, energy, geomName
    type(outputFile)                                    :: test_out_critical
    type(visualiser)                                    :: viz
    class(field), pointer                               :: field
    character(100), parameter :: Here ='init (criticalKineticPhysicsPackage_class.f90)'

    call cpu_time(self % CPU_time_start)

    ! Initialise general + criticality calcs first

    ! Read calculation settings
    call dict % get( self % pop,'pop')
    call dict % get( self % N_inactive,'inactive')
    call dict % get( self % N_active,'active')
    call dict % get( nucData, 'XSdata')
    call dict % get( energy, 'dataType')

    ! Parallel buffer size
    call dict % getOrDefault( self % bufferSize, 'buffer', 1000)

    ! Process type of data
    select case(energy)
      case('ce')
        self % particleType = P_NEUTRON_CE
      case default
        call fatalError(Here,"dataType must be 'mg' or 'ce'.")
    end select

    ! Read outputfile path
    call dict % getOrDefault(self % outputFileCritical,'outputFile','./outputCritical')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be cought early)
    call dict % getOrDefault(self % outputFormatCritical, 'outputFormat', 'asciiMATLAB')
    call test_out_critical % init(self % outputFormatCritical)

    ! Register timer
    self % timerMain = registerTimer('transportTime')

    ! Initialise RNG
    allocate(self % pRNG)

    ! *** It is a bit silly but dictionary cannot store longInt for now
    !     so seeds are limited to 32 bits (can be -ve)
    if( dict % isPresent('seed')) then
      call dict % get(seed_temp,'seed')

    else
      ! Obtain time string and hash it to obtain random seed
      call date_and_time(date, time)
      string = date // time
      call FNV_1(string,seed_temp)

    end if
    seed = seed_temp
    call self % pRNG % init(seed)

    ! Initial k_effective guess
    call dict % getOrDefault(self % keff_0,'keff_0', ONE)

    ! Read whether to print particle source per cycle
    call dict % getOrDefault(self % printSource, 'printSource', 0)

    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'eigenGeom'
    call new_geometry(tempDict, geomName)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Activate Nuclear Data *** All materials are active
    call ndReg_activate(self % particleType, nucData, self % geom % activeMats())
    self % nucData => ndReg_get(self % particleType)

    ! Call visualisation
    if (dict % isPresent('viz')) then
      print *, "Initialising visualiser"
      tempDict => dict % getDictPtr('viz')
      call viz % init(self % geom, tempDict)
      print *, "Constructing visualisation"
      call viz % makeViz()
      call viz % kill()
    endif

    ! Read uniform fission site option as a geometry field
    if (dict % isPresent('uniformFissionSites')) then
      self % ufs = .true.
      ! Build and initialise
      tempDict => dict % getDictPtr('uniformFissionSites')
      call new_field(tempDict, nameUFS)
      ! Save UFS field
      field => gr_fieldPtr(gr_fieldIdx(nameUFS))
      self % ufsField => uniFissSitesField_TptrCast(field)
      ! Initialise
      call self % ufsField % estimateVol(self % geom, self % pRNG, self % particleType)
    end if

    ! Read variance reduction option as a geometry field
    if (dict % isPresent('varianceReduction')) then
      ! Build and initialise
      tempDict => dict % getDictPtr('varianceReduction')
      call new_field(tempDict, nameWW)
    end if

    ! Build collision operator
    tempDict => dict % getDictPtr('collisionOperatorCritical')
    call self % collOpCritical % init(tempDict)

    ! Build transport operator
    tempDict => dict % getDictPtr('transportOperatorCritical')
    call new_transportOperator(self % transOpCritical, tempDict)

    ! Initialise active & inactive tally Admins
    tempDict => dict % getDictPtr('inactiveTally')
    allocate(self % inactiveTally)
    call self % inactiveTally % init(tempDict)

    tempDict => dict % getDictPtr('activeTally')
    allocate(self % activeTally)
    call self % activeTally % init(tempDict)

    ! Load Initial source
    if (dict % isPresent('source')) then ! Load definition from file
      call new_source(self % initSource, dict % getDictPtr('source'), self % geom)

    else
      call locDict1 % init(3)
      call locDict1 % store('type', 'fissionSource')
      call locDict1 % store('data', trim(energy))
      call new_source(self % initSource, locDict1, self % geom)
      call locDict1 % kill()

    end if

    ! Initialise active and inactive tally attachments
    ! Inactive tally attachment
    call locDict1 % init(2)
    call locDict2 % init(2)

    call locDict2 % store('type','keffAnalogClerk')
    call locDict1 % store('keff', locDict2)
    call locDict1 % store('display',['keff'])

    allocate(self % inactiveAtch)
    call self % inactiveAtch % init(locDict1)

    call locDict2 % kill()
    call locDict1 % kill()

    ! Active tally attachment
    call locDict1 % init(2)
    call locDict2 % init(2)

    call locDict2 % store('type','keffImplicitClerk')
    call locDict1 % store('keff', locDict2)
    call locDict1 % store('display',['keff'])

    allocate(self % activeAtch)
    call self % activeAtch % init(locDict1)

    call locDict2 % kill()
    call locDict1 % kill()

    ! Attach attachments to result tallies
    call self % inactiveTally % push(self % inactiveAtch)
    call self % activeTally % push(self % activeAtch)

    call dict % get( self % N_stepBacks,'decorrelationCycles')

    allocate(self % dict)
    self % dict = dict

    call self % printSettingsCritical()

  end subroutine init


  subroutine run(self)
    class(criticalKineticPhysicsPackage), intent(inout) :: self

    print *, repeat("<>",50)
    print *, "/\/\ EIGENVALUE CALCULATION /\/\"

    call self % generateInitialState()
    call self % cyclesCritical(self % inactiveTally, self % inactiveAtch, self % N_inactive)
    call self % cyclesCritical(self % activeTally, self % activeAtch, self % N_active - 1)
    call self % collectResultsCritical()

    print *
    print *, "\/\/ END OF EIGENVALUE CALCULATION \/\/"
    print *

    print *, repeat("<>",50)
    print *, "/\/\ INITIALISING TIME DEPENDENT /\/\"
    call self % initKinetic()
    print *
    print *, "\/\/ END OF TIME DEPENDENT INITIALISATION \/\/"
    print *

    print *, repeat("<>",50)
    print *, "/\/\ COUPLING TO TIME DEPENDENT /\/\"
    call self % coupleKinetic(self % inactiveTally, self % inactiveAtch, self % N_stepBacks)
    print *
    print *, "\/\/ END OF TIME DEPENDENT COUPLING \/\/"
    print *

    print *, repeat("<>",50)
    print *, "/\/\ TIME DEPENDENT CALCULATION /\/\"

    call self % cyclesKinetic(self % tally, self % N_cycles, self % N_timeBins, self % timeIncrement)
    call self % collectResultsKinetic()

    print *
    print *, "\/\/ END OF TIME DEPENDENT CALCULATION \/\/"
    print *
  end subroutine

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettingsCritical(self)
    class(criticalKineticPhysicsPackage), intent(in) :: self

    print *, repeat("<>",50)
    print *, "/\/\ EIGENVALUE CALCULATION WITH POWER ITERATION METHOD /\/\"
    print *, "Inactive Cycles:    ", numToChar(self % N_inactive)
    print *, "Active Cycles:      ", numToChar(self % N_active)
    print *, "Neutron Population: ", numToChar(self % pop)
    print *, "Initial RNG Seed:   ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettingsCritical

  !!
  !!
  !!
  subroutine generateInitialState(self)
    class(criticalKineticPhysicsPackage), intent(inout) :: self
    character(100), parameter :: Here =' generateInitialState( criticalKineticPhysicsPackage_class.f90)'

    ! Allocate and initialise particle Dungeons
    allocate(self % thisCycle)
    allocate(self % nextCycle)
    call self % thisCycle % init(10 * self % pop)
    call self % nextCycle % init(10 * self % pop)

    ! Generate initial surce
    print *, "GENERATING INITIAL FISSION SOURCE"
    call self % initSource % generate(self % thisCycle, self % pop, self % pRNG)
    print *, "DONE!"

  end subroutine generateInitialState

  !!
  !! Critical calc
  !!
  subroutine cyclesCritical(self, tally, tallyAtch, N_cycles)
    class(criticalKineticPhysicsPackage), intent(inout) :: self
    type(tallyAdmin), pointer,intent(inout)             :: tally
    type(tallyAdmin), pointer,intent(inout)             :: tallyAtch
    integer(shortInt), intent(in)                       :: N_cycles
    type(particleDungeon), save                         :: buffer
    integer(shortInt)                                   :: i, n, Nstart, Nend, nParticles
    class(tallyResult),allocatable                      :: res
    type(collisionOperator), save                       :: collOpCritical
    class(transportOperator),allocatable,save           :: transOpCritical
    type(RNG), target, save                             :: pRNG
    type(particle), save                                :: neutron
    real(defReal)                                       :: k_old, k_new
    real(defReal)                                       :: elapsed_T, end_T, T_toEnd
    character(100),parameter :: Here ='cycles (criticalKineticPhysicsPackage_class.f90)'
    !$omp threadprivate(neutron, buffer, collOpCritical, transOpCritical, pRNG)

    !$omp parallel
    ! Create particle buffer
    call buffer % init(self % bufferSize)

    ! Initialise neutron
    neutron % geomIdx = self % geomIdx

    ! Create a collision + transport operator which can be made thread private
    collOpCritical = self % collOpCritical
    transOpCritical = self % transOpCritical
    !$omp end parallel

    ! Set initial k-eff
    k_new = self % keff_0

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    do i=1,N_cycles

      ! Send start of cycle report
      Nstart = self % thisCycle % popSize()
      call tally % reportCycleStart(self % thisCycle)

      nParticles = self % thisCycle % popSize()

      !$omp parallel do schedule(dynamic)
      gen: do n = 1, nParticles

        ! TODO: Further work to ensure reproducibility!
        ! Create RNG which can be thread private
        pRNG = self % pRNG
        neutron % pRNG => pRNG
        call neutron % pRNG % stride(n)

        ! Obtain particle current cycle dungeon
        call self % thisCycle % copy(neutron, n)

        bufferLoop: do
          call self % geom % placeCoord(neutron % coords)

          ! Set k-eff for normalisation in the particle
          neutron % k_eff = k_new

          ! Save state
          call neutron % savePreHistory()

          ! Transport particle untill its death
          history: do
            call transOpCritical % transport(neutron, tally, buffer, self % nextCycle)
            if(neutron % isDead) exit history

            call collOpCritical % collide(neutron, tally, buffer, self % nextCycle)
            if(neutron % isDead) exit history
          end do history

          ! Clear out buffer
          if (buffer % isEmpty()) then
            exit bufferLoop
          else
            call buffer % release(neutron)
          end if

        end do bufferLoop

      end do gen
      !$omp end parallel do

      call self % thisCycle % cleanPop()

      ! Update RNG
      call self % pRNG % stride(self % pop + 1)

      ! Send end of cycle report
      Nend = self % nextCycle % popSize()
      call tally % reportCycleEnd(self % nextCycle)

      if (self % UFS) then
        call self % ufsField % updateMap()
      end if

      ! Normalise population
      call self % nextCycle % normSize(self % pop, pRNG)

      ! Flip cycle dungeons
      self % temp_dungeon => self % nextCycle
      self % nextCycle    => self % thisCycle
      self % thisCycle    => self % temp_dungeon

      ! Obtain estimate of k_eff
      call tallyAtch % getResult(res,'keff')

      select type(res)
        class is(keffResult)
          k_new = res % keff(1)

        class default
          call fatalError(Here, 'Invalid result has been returned')

      end select

      ! Load new k-eff estimate into next cycle dungeon
      k_old = self % nextCycle % k_eff
      self % nextCycle % k_eff = k_new

      ! Used to normalise fission source of the first active cycle
      self % keff_0 = k_new

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(N_cycles,defReal) * elapsed_T / i
      T_toEnd = max(ZERO, end_T - elapsed_T)


      ! Display progress
      call printFishLineR(i)
      print *
      print *, 'Cycle: ', numToChar(i), ' of ', numToChar(N_cycles)
      print *, 'Pop: ', numToChar(Nstart) , ' -> ', numToChar(Nend)
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))
      call tally % display()
    end do

    ! Load elapsed time
    self % time_transport = self % time_transport + elapsed_T

  end subroutine cyclesCritical

  !!
  !! Print calculation criticality results to file
  !!
  subroutine collectResultsCritical(self)
    class(criticalKineticPhysicsPackage), intent(inout) :: self
    type(outputFile)                          :: out
    character(nameLen)                        :: name

    call out % init(self % outputFormatCritical)

    name = 'seed'
    call out % printValue(self % pRNG % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Inactive_Cycles'
    call out % printValue(self % N_inactive,name)

    name = 'Active_Cycles'
    call out % printValue(self % N_active,name)

    call cpu_time(self % CPU_time_end)
    name = 'Total_CPU_Time'
    call out % printValue((self % CPU_time_end - self % CPU_time_start),name)

    name = 'Total_Transport_Time'
    call out % printValue(self % time_transport,name)

    ! Print Inactive tally
    name = 'inactive'
    call out % startBlock(name)
    call self % inactiveTally % print(out)
    call out % endBlock()

    ! Print Active attachment
    ! Is printed into the root block
    call self % activeAtch % print(out)

    name = 'active'
    call out % startBlock(name)
    call self % activeTally % print(out)
    call out % endBlock()

    call out % writeToFile(self % outputFileCritical)

  end subroutine collectResultsCritical

  !!
  !! Initialise from individual components and dictionaries for inactive and active tally
  !!
  subroutine initKinetic(self)
    class(criticalKineticPhysicsPackage), intent(inout) :: self
    class(dictionary),pointer                           :: tempDict
    integer(shortInt)                                   :: i
    character(len=10)                                   :: string
    type(outputFile)                                    :: test_out_kinetic
    logical(defBool)                                    :: useKineticGeom

    ! Read outputfile path
    call self % dict % getOrDefault(self % outputFileKinetic,'outputFile','./outputKinetic')

    ! Get output format and verify
    ! Initialise output file before calculation (so mistake in format will be cought early)
    call self % dict % getOrDefault(self % outputFormatKinetic, 'outputFormat', 'asciiMATLAB')
    call test_out_kinetic % init(self % outputFormatKinetic)

    call self % dict % get( self % N_cycles,'cycles')
    call self % dict % get( self % N_timeBins,'timeSteps')
    call self % dict % get( self % timeIncrement, 'timeIncrement')

    ! Whether to use implicit for weight normalisation (default = no)
    call self % dict % getOrDefault(self % criticalWgtNorm, 'criticalWgtNorm', .false.)

    ! Whether to use combing (default = no)
    call self % dict % getOrDefault(self % useCombing, 'combing', .false.)

    ! Whether to implement precursors (default = yes)
    call self % dict % getOrDefault(self % usePrecursors, 'precursors', .false.)

    ! Whether to use analog or implicit kinetic (default = Analog)
    call self % dict % getOrDefault(self % useForcedPrecursorDecay, 'useForcedPrecursorDecay', .false.)

    ! Build collision operator
    tempDict => self % dict % getDictPtr('collisionOperatorKinetic')
    call self % collOpKinetic % init(tempDict)

    ! Build transport operator
    tempDict => self % dict % getDictPtr('transportOperatorKinetic')
    call new_transportOperator(self % transOpKinetic, tempDict)

    ! Initialise tally Admin for kinetic
    tempDict => self % dict % getDictPtr('tally')
    allocate(self % tally)
    call self % tally % init(tempDict)

    ! Initialise kinetic geometry
    tempDict => self % dict % getDictPtr('kineticGeometry')
    call tempDict % getOrDefault(useKineticGeom, 'enable', .false.)
    if (useKineticGeom .eqv. .true.) then
      call tempDict % get(self % kineticGeomTimeSteps, 'timeSteps')
      allocate(self % kineticGeomDict(size(self % kineticGeomTimeSteps)))
      do i=1, size(self % kineticGeomTimeSteps)
        write(string, '(I10)') i
        self % kineticGeomDict(i) = self % dict % getDictPtr('kineticGeom' // trim(adjustl(string)))
      end do
    end if

    ! deallocate dict
    call self % dict % kill()
    deallocate(self % dict)

    ! Size particle dungeon
    allocate(self % currentTime(self % N_cycles))
    allocate(self % nextTime(self % N_cycles))

    do i = 1, self % N_cycles
      call self % currentTime(i) % init(5*self % bufferSize)
      call self % nextTime(i) % init(5*self % bufferSize)
    end do

    ! Size precursor dungeon
    if (self % usePrecursors) then
      allocate(self % precursorDungeons(self % N_cycles))
      do i = 1, self % N_cycles
        call self % precursorDungeons(i) % init(3*self % bufferSize)
      end do
    end if

    call self % printSettingsKinetic()

  end subroutine initKinetic

  !!
  !! Print settings of the physics package
  !!
  subroutine printSettingsKinetic(self)
    class(criticalKineticPhysicsPackage), intent(in) :: self
    real(defReal)                                    :: TStart, Tstop, Tincrement

    TStart = 0.0
    Tstop = self % timeIncrement * self % N_timeBins
    Tincrement = self % timeIncrement
    print *, repeat("<>",50)
    print *, "/\/\ TIME DEPENDENT CALCULATION /\/\"
    print *, "Time grid [start, stop, increment]: ", numToChar(TStart), numToChar(Tstop), numToChar(Tincrement)
    print *, "Source batches:                     ", numToChar(self % N_cycles)
    print *, "Initial Population per batch:       ", numToChar(self % pop)
    print *, "Initial RNG Seed:                   ", numToChar(self % pRNG % getSeed())
    print *
    print *, repeat("<>",50)
  end subroutine printSettingsKinetic

  subroutine coupleKinetic(self, tally, tallyAtch, N_stepBacks)
    class(criticalKineticPhysicsPackage), intent(inout) :: self
    integer(shortInt), intent(in)                       :: N_stepBacks
    type(tallyAdmin), pointer,intent(inout)             :: tally
    type(tallyAdmin), pointer,intent(inout)             :: tallyAtch
    type(particleDungeon), save                         :: buffer
    integer(shortInt)                                   :: i, n, Nstart, Nend, nParticles, Ncycles, batchIdx
    integer(shortInt)                                   :: batchTracker
    class(tallyResult),allocatable                      :: res
    type(collisionOperator), save                       :: collOpCritical
    class(transportOperator),allocatable,save           :: transOpCritical
    type(RNG), target, save                             :: pRNG
    type(particle), save                                :: neutron
    real(defReal)                                       :: k_old, k_new
    real(defReal)                                       :: elapsed_T, end_T, T_toEnd
    character(100),parameter :: Here ='coupleKinetic (criticalKineticPhysicsPackage_class.f90)'
    !$omp threadprivate(neutron, buffer, collOpCritical, transOpCritical, pRNG)

    !$omp parallel
    ! Initialise neutron
    neutron % geomIdx = self % geomIdx

    ! Create a collision + transport operator which can be made thread private
    collOpCritical = self % collOpCritical
    transOpCritical = self % transOpCritical
    !$omp end parallel

    ! Set initial k-eff
    k_new = self % keff_0

    Ncycles = self % N_cycles + N_stepBacks * (self % N_cycles - 1)
    do i = 1, Ncycles
      
      batchIdx = int((i - 1) / (self % N_stepBacks + 1)) + 1

      ! Send start of cycle report
      Nstart = self % thisCycle % popSize()
      call tally % reportCycleStart(self % thisCycle)

      nParticles = self % thisCycle % popSize()

      !$omp parallel do schedule(dynamic)
      gen: do n = 1, nParticles

        ! TODO: Further work to ensure reproducibility!
        ! Create RNG which can be thread private
        pRNG = self % pRNG
        neutron % pRNG => pRNG
        call neutron % pRNG % stride(n)

        ! Obtain particle current cycle dungeon
        call self % thisCycle % copy(neutron, n)
        neutron % criticalSource = .false.

        bufferLoop: do
          call self % geom % placeCoord(neutron % coords)

          ! Set k-eff for normalisation in the particle
          neutron % k_eff = k_new

          ! Save state
          call neutron % savePreHistory()

          ! Transport particle untill its death
          history: do
            call transOpCritical % transport(neutron, tally, buffer, self % nextCycle)
            if(neutron % isDead) exit history

            if (mod(i-1, self % N_stepBacks + 1) == 0) then
              neutron % criticalSource = .true.
              call collOpCritical % collide(neutron, tally, self % precursorDungeons(batchIdx), self % currentTime(batchIdx))
            else
              call collOpCritical % collide(neutron, tally, buffer, self % nextCycle)
            end if
            if(neutron % isDead) exit history
          end do history

          ! Clear out buffer
          if (buffer % isEmpty()) then
            exit bufferLoop
          else
            call buffer % release(neutron)
          end if

        end do bufferLoop

      end do gen
      !$omp end parallel do

      ! Update RNG
      call self % pRNG % stride(self % pop + 1)

      ! Send end of cycle report
      Nend = self % nextCycle % popSize()
      call tally % reportCycleEnd(self % nextCycle)


      if (mod(i-1, self % N_stepBacks + 1) /= 0) then

        call self % thisCycle % cleanPop()

        if (self % UFS) then
          call self % ufsField % updateMap()
        end if

        ! Normalise population
        call self % nextCycle % normSize(self % pop, pRNG)

        ! Flip cycle dungeons
        self % temp_dungeon => self % nextCycle
        self % nextCycle    => self % thisCycle
        self % thisCycle    => self % temp_dungeon

        ! Obtain estimate of k_eff
        call tallyAtch % getResult(res,'keff')

        select type(res)
          class is(keffResult)
            k_new = res % keff(1)

          class default
            call fatalError(Here, 'Invalid result has been returned')

        end select

        ! Load new k-eff estimate into next cycle dungeon
        k_old = self % nextCycle % k_eff
        self % nextCycle % k_eff = k_new

        ! Used to normalise fission source of the first active cycle
        self % keff_0 = k_new
      else
        batchTracker = batchTracker + 1
        print *, 'Generated critical batch: ', batchTracker, '/', self % N_cycles
      end if
  
    end do

    ! Normalise neutron and precursor weights to ensure dimension-less (only if implicit)
    if (self % criticalWgtNorm .eqv. .true.) call self % normCriticalWgts()

    ! Free memory
    allocate(self % thisCycle)
    call self % thisCycle % kill()
    deallocate(self % thisCycle)

    allocate(self % nextCycle)
    call self % nextCycle % kill()
    deallocate(self % nextCycle)

    allocate(self % temp_dungeon)
    call self % temp_dungeon % kill()
    deallocate(self % temp_dungeon)

    call self % inactiveTally % kill()
    call self % inactiveAtch % kill()
    call self % activeAtch % kill()
    call self % transOpCritical % kill()
    if (allocated(self % transOpCritical)) deallocate(self % transOpCritical)
    call self % collOpCritical % kill()

  end subroutine coupleKinetic

  !!
  !!
  !!
  subroutine cyclesKinetic(self, tally, N_cycles, N_timeBins, timeIncrement)
    class(criticalKineticPhysicsPackage), intent(inout) :: self
    type(tallyAdmin), pointer,intent(inout)             :: tally
    integer(shortInt), intent(in)                       :: N_timeBins, N_cycles
    integer(shortInt)                                   :: i, t, n
    integer(shortInt)                                   :: nParticles, nDelayedParticles, nPrecuCount
    type(particle), save                                :: p, p_d
    type(particleDungeon), save                         :: buffer
    type(collisionOperator), save                       :: collOpKinetic
    class(transportOperator),allocatable,save           :: transOpKinetic
    type(RNG), target, save                             :: pRNG
    real(defReal) , save                                :: decay_T
    real(defReal)                                       :: elapsed_T, end_T, T_toEnd
    real(defReal), intent(in)                           :: timeIncrement
    character(nameLen)                                  :: geomName
    integer(shortInt), dimension(:), allocatable        :: kineticGeomLocs
    integer(shortInt)                                   :: kineticGeomIdx
    character(100),parameter :: Here ='cycles (timeDependentPhysicsPackage_class.f90)'
    !$omp threadprivate(p, p_d, buffer, collOpKinetic, transOpKinetic, pRNG, decay_T)

    !$omp parallel
    ! Create particle buffer
    call buffer % init(self % bufferSize)

    ! Initialise neutron
    p % geomIdx = self % geomIdx
    p % k_eff = ONE

    ! Create a collision + transport operator which can be made thread private
    collOpKinetic = self % collOpKinetic
    transOpKinetic = self % transOpKinetic
    !$omp end parallel

    ! Number of particles in each batch
    nParticles = self % pop

    ! Reset and start timer
    call timerReset(self % timerMain)
    call timerStart(self % timerMain)

    do t = 1, N_timeBins

      ! Handle kinetic geometry
      if (allocated(self % kineticGeomTimeSteps)) then
        if (ANY(self % kineticGeomTimeSteps == t)) then
          call self % geom % kill()
          call killGeom()
          kineticGeomLocs = findloc(self % kineticGeomTimeSteps, t)
          kineticGeomIdx = kineticGeomLocs(1)
          geomName = 'kineticGeom'
          call new_geometry(self % kineticGeomDict(kineticGeomIdx), geomName)
          self % geomIdx = gr_geomIdx(geomName)
          self % geom    => gr_geomPtr(self % geomIdx)
          p % geomIdx = self % geomIdx
        end if
      end if

      do i = 1, N_cycles

        if ((t == 1) .and. (self % useCombing .eqv. .true.)) call self % currentTime(i) % combing(self % bufferSize, pRNG)
        nParticles = self % currentTime(i) % popSize()

        if ((self % usePrecursors .eqv. .true.) .and. (self % useForcedPrecursorDecay .eqv. .true.)) then
          if (t == 1) then
            nPrecuCount = 1
          else 
            nPrecuCount = self % precursorDungeons(i) % popSize() + 1
          end if
        end if

        call tally % reportCycleStart(self % currentTime(i))

        !$omp parallel do schedule(dynamic)
        gen: do n = 1, nParticles
          pRNG = self % pRNG
          p % pRNG => pRNG
          call p % pRNG % stride(n)
          call self % currentTime(i) % copy(p, n)
          p % timeMax = t * timeIncrement

          bufferLoop: do

            if ((p % fate == aged_FATE) .or. (p % fate == no_FATE)) then
              p % fate = no_FATE
              p % isdead = .false.
            else
              p % isdead = .true.
            end if

            call self % geom % placeCoord(p % coords)
            call p % savePreHistory()

            ! Transport particle untill its death
            history: do
              if(p % isDead) exit history
              call transOpKinetic % transport(p, tally, buffer, buffer)
              if(p % isDead) exit history
              if(p % fate == AGED_FATE) then
                call self % nextTime(i) % detain(p)
                exit history
              endif
              if (self % usePrecursors) then
                call collOpKinetic % collide(p, tally, self % precursorDungeons(i), buffer)
              else
                call collOpKinetic % collide(p, tally, buffer, buffer)
              end if
              if(p % isDead) exit history
            end do history

            ! Clear out buffer
            if (buffer % isEmpty()) then
              exit bufferLoop
            else
              call buffer % release(p)
            end if

          end do bufferLoop
        end do gen
        !$omp end parallel do

        if ((self % usePrecursors .eqv. .true.) .and. (self % useForcedPrecursorDecay .eqv. .false.)) then

          ! Analog delayed neutron handling
          nDelayedParticles = self % precursorDungeons(i) % popSize()
          nPrecuCount = 1
          if (nDelayedParticles > 0) then
            superGenDelayed: do
            !$omp parallel do schedule(dynamic)
            genDelayed: do n = nPrecuCount, nDelayedParticles
              call self % precursorDungeons(i) % copy(p, n)

              if ((p % time <= t*timeIncrement) .and. (p % time > (t-1)*timeIncrement)) then
                p % type = P_NEUTRON
                pRNG = self % pRNG
                p % pRNG => pRNG
                call p % pRNG % stride(n)

                bufferLoopDelayed: do

                  if ((p % fate == aged_FATE) .or. (p % fate == no_FATE)) then
                    p % fate = no_FATE
                    p % isdead = .false.
                  else
                    p % isdead = .true.
                  end if

                  call self % geom % placeCoord(p % coords)
                  p % timeMax = t*timeIncrement
                  call p % savePreHistory()

                  ! Transport particle until its death
                  historyDelayed: do
                    if(p % isDead) exit historyDelayed
                    call transOpKinetic % transport(p, tally, buffer, buffer)
                    if(p % isDead) exit historyDelayed
                    if(p % fate == AGED_FATE) then
                      call self % nextTime(i) % detain(p)
                      exit historyDelayed
                    endif
                    call collOpKinetic % collide(p, tally, self % precursorDungeons(i), buffer)
                    if(p % isDead) exit historyDelayed
                  end do historyDelayed

                  ! Clear out buffer
                  if (buffer % isEmpty()) then
                    exit bufferLoopDelayed
                  else
                    call buffer % release(p)
                  end if

                end do bufferLoopDelayed
              end if
            end do genDelayed
            !$omp end parallel do
            if (nDelayedParticles .eq. self % precursorDungeons(i) % popSize()) then
               exit superGenDelayed
            else
              nPrecuCount = nDelayedParticles+1
              nDelayedParticles = self % precursorDungeons(i) % popSize()
            end if
          end do superGenDelayed
          end if

        else if ((self % usePrecursors .eqv. .true.) .and. (self % useForcedPrecursorDecay .eqv. .true.)) then

          ! Impliciy delayed neutron handling
          nDelayedParticles = self % precursorDungeons(i) % popSize()

          if (nDelayedParticles > nPrecuCount - 1) then
            superGenDelayedImp: do
            !$omp parallel do schedule(dynamic)
            genDelayedImp: do n = nPrecuCount, nDelayedParticles

              ! Force decay at this time interval if generated in this time interval
              call self % precursorDungeons(i) % copy(p, n)

              decay_T = p % time + pRNG % get() * (t*timeIncrement - p % time)

              ! Weight adjustment
              p % w = p % forcedPrecursorDecayWgt(decay_T, t*timeIncrement - p % time)

              pRNG = self % pRNG
              p % pRNG => pRNG
              call p % pRNG % stride(n)

              ! Update parameters
              p % type = P_NEUTRON
              p % time = decay_T
              p % fate = no_FATE

              bufferLoopDelayedImp: do

                if ((p % fate == aged_FATE) .or. (p % fate == no_FATE)) then
                  p % fate = no_FATE
                  p % isdead = .false.
                else
                  p % isdead = .true.
                end if

                call self % geom % placeCoord(p % coords)
                p % timeMax = t*timeIncrement
                call p % savePreHistory()

                ! Transport particle until its death
                historyDelayedImp: do
                  if(p % isDead) exit historyDelayedImp
                  call transOpKinetic % transport(p, tally, buffer, buffer)
                  if(p % isDead) exit historyDelayedImp
                  if(p % fate == AGED_FATE) then
                    call self % nextTime(i) % detain(p)
                    exit historyDelayedImp
                  endif
                  call collOpKinetic % collide(p, tally, self % precursorDungeons(i), buffer)
                  if(p % isDead) exit historyDelayedImp
                end do historyDelayedImp

                ! Clear out buffer
                if (buffer % isEmpty()) then
                  exit bufferLoopDelayedImp
                else
                  call buffer % release(p)
                end if

              end do bufferLoopDelayedImp

            end do genDelayedImp
            !$omp end parallel do

            if (nDelayedParticles .eq. self % precursorDungeons(i) % popSize()) then
              exit superGenDelayedImp
            else
              nPrecuCount = nDelayedParticles+1
              nDelayedParticles = self % precursorDungeons(i) % popSize()
            end if
            end do superGenDelayedImp

          end if

          nDelayedParticles = self % precursorDungeons(i) % popSize()

          if (nDelayedParticles > 0) then

            ! Precursor population control
            if (nDelayedParticles > self % bufferSize) then
              call self % precursorDungeons(i) % precursorCombing(self % bufferSize, pRNG, timeIncrement*t)
            end if

            nDelayedParticles = self % precursorDungeons(i) % popSize()

            !$omp parallel do schedule(dynamic)
            genDelayedImpNext: do n = 1, nDelayedParticles
              ! pass to next time interval for Forced Precursor Decay

              call self % precursorDungeons(i) % copy(p, n)

              ! Sample decay time
              decay_T = timeIncrement * (t+pRNG % get())

              ! Weight adjustment
              p % w = p % forcedPrecursorDecayWgt(decay_T, timeIncrement)

              pRNG = self % pRNG
              p % pRNG => pRNG
              call p % pRNG % stride(n)

              ! Update parameters
              p % type = P_NEUTRON
              p % time = decay_T
              p % fate = no_FATE

              ! Add to current dungeon
              call self % nextTime(i) % detain(p)

            end do genDelayedImpNext
            !$omp end parallel do
          end if

        end if

        ! Update RNG
        call self % pRNG % stride(nParticles + 1)
        call tally % reportCycleEnd(self % currentTime(i))
        call self % currentTime(i) % cleanPop()

        ! Neutron population control
        if (self % useCombing) then
          call self % nextTime(i) % combing(self % bufferSize, pRNG)
        else if ((self % usePrecursors .eqv. .true.) .and. (self % useForcedPrecursorDecay .eqv. .true.)) then
          call self % nextTime(i) % combing(self % bufferSize, pRNG)
        end if

      end do

      self % tempTime  => self % nextTime
      self % nextTime  => self % currentTime
      self % currentTime => self % tempTime

      ! Calculate times
      call timerStop(self % timerMain)
      elapsed_T = timerTime(self % timerMain)

      ! Predict time to end
      end_T = real(N_timeBins,defReal) * elapsed_T / t
      T_toEnd = max(ZERO, end_T - elapsed_T)

      ! Display progress
      call printFishLineR(t)
      print *
      print *, 'Time step: ', numToChar(t), ' of ', numToChar(N_timeBins)
      print *, 'Elapsed time: ', trim(secToChar(elapsed_T))
      print *, 'End time:     ', trim(secToChar(end_T))
      print *, 'Time to end:  ', trim(secToChar(T_toEnd))
      call tally % display()

      call tally % setNumBatchesPerTimeStep(N_cycles)
    end do

  end subroutine cyclesKinetic

  !!
  !! Print calculation criticality results to file
  !!
  subroutine collectResultsKinetic(self)
    class(criticalKineticPhysicsPackage), intent(inout) :: self
    type(outputFile)                          :: out
    character(nameLen)                        :: name

    call out % init(self % outputFormatKinetic)

    name = 'seed'
    call out % printValue(self % pRNG % getSeed(),name)

    name = 'pop'
    call out % printValue(self % pop,name)

    name = 'Source_batches'
    call out % printValue(self % N_cycles,name)

    name = 'Time_increment'
    call out % printValue(self % timeIncrement,name)

    name = 'Time_bins'
    call out % printValue(self % N_timeBins,name)

    call cpu_time(self % CPU_time_end)
    name = 'Total_CPU_Time'
    call out % printValue((self % CPU_time_end - self % CPU_time_start),name)

    name = 'Transport_time'
    call out % printValue(timerTime(self % timerMain),name)

    ! Print tally
    call self % tally % print(out)

    call out % writeToFile(self % outputFileKinetic)

  end subroutine collectResultsKinetic

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(criticalKineticPhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

  !!
  !! Normalise particle weighs pre-kinetic calculations
  !!
  subroutine normCriticalWgts(self)
    class(criticalKineticPhysicsPackage), intent(inout) :: self
    integer(shortInt)                                   :: i
    real(defReal), save                                 :: nWgt, pWgt, totWgt, normFactor
    !$omp threadprivate(nWgt, pWgt, totWgt, normFactor)

    !$omp parallel do schedule(dynamic)
    do i = 1, self % N_cycles
      if (self % usePrecursors .eqv. .true.) then
        nWgt = self % currentTime(i) % popWeight()
        pWgt = self % precursorDungeons(i) % popWeight()
        totWgt = nWgt + pWgt
        normFactor = totWgt / (self % precursorDungeons(i) % popSize() + self % currentTime(i) % popSize())
        call self % currentTime(i) % scaleWeight(ONE / normFactor)
        call self % precursorDungeons(i) % scaleWeight(ONE / normFactor)
      else 
        nWgt = self % currentTime(i) % popWeight()
        totWgt = nWgt
        normFactor = totWgt / (self % currentTime(i) % popSize())
        call self % currentTime(i) % scaleWeight(ONE / normFactor)
      end if
    end do
    !$omp end parallel do

  end subroutine normCriticalWgts

  subroutine russianRoulette(self, p, avWgt)
    class(criticalKineticPhysicsPackage), intent(inout) :: self
    class(particle), intent(inout)     :: p
    real(defReal), intent(in)          :: avWgt

    if (p % pRNG % get() < (ONE - p % w/avWgt)) then
      p % isDead = .true.
    else
      p % w = avWgt
    end if

  end subroutine russianRoulette

end module criticalKineticPhysicsPackage_class