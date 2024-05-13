module nuAnalogClerk_class

  use numPrecision
  use tallyCodes
  use endfConstants
  use dictionary_class,      only : dictionary
  use particle_class,        only : particle, particleState, P_PRECURSOR, P_NEUTRON
  use particleDungeon_class, only : particleDungeon
  use outputFile_class,      only : outputFile

  use scoreMemory_class,     only : scoreMemory
  use tallyClerk_inter,      only : tallyClerk, kill_super => kill

  ! Tally Maps
  use tallyMap_inter,             only : tallyMap
  use tallyMapFactory_func,       only : new_tallyMap

  ! Nuclear Data Interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase

  implicit none
  private

  !! Locations of diffrent bins wrt memory address of the clerk
  integer(shortInt), parameter :: MEM_SIZE = 3
  integer(longInt), parameter  :: NU_FISS = 0, &  ! Number of neutrons generated by fission
                                  EVENT   = 1, &  ! Number of events
                                  NU_MEAN = 2     ! Averaged value

  !!
  !! Analog estimator for nu prompt or delayed
  !!
  !! Private Members:
  !!   delayed -> flag that indicated whether nu prompt or nu delayed is requested
  !!
  !! SAMPLE DICTIOANRY INPUT:
  !!
  !! myClerk {
  !!   type nuAnalogClerk; #delayed = 1;#
  !!   #map { <tallyMap definition> }#
  !! }
  !!
  type, public,extends(tallyClerk) :: nuAnalogClerk
    private
    logical(defBool) :: delayed = .false.
    class(tallyMap), allocatable :: map
  contains
    ! Procedures used during build
    procedure :: init
    procedure :: kill
    procedure :: validReports
    procedure :: getSize

    ! File reports and check status -> run-time procedures
    procedure :: reportSpawn
    procedure :: reportOutColl
    procedure :: reportCycleEnd

    ! Output procedures
    procedure  :: display
    procedure  :: print

  end type nuAnalogClerk

contains

  !!
  !! Initialise from dictionary and name
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine init(self, dict, name)
    class(nuAnalogClerk), intent(inout) :: self
    class(dictionary), intent(in)       :: dict
    character(nameLen), intent(in)      :: name

    ! Needs no settings, just load name
    call self % setName(name)

    ! Initialise delayed flag
    call dict % getOrDefault(self % delayed, 'delayed', .false.)

    ! Load map
    if (dict % isPresent('map')) then
      call new_tallyMap(self % map, dict % getDictPtr('map'))
    end if

  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(nuAnalogClerk), intent(inout) :: self

    ! Superclass
    call kill_super(self)

    ! Kill and deallocate map
    if (allocated(self % map)) then
      call self % map % kill()
      deallocate(self % map)
    end if

    self % delayed = .false.

  end subroutine kill

  !!
  !! Returns array of codes that represent diffrent reports
  !!
  !! See tallyClerk_inter for details
  !!
  function validReports(self) result(validCodes)
    class(nuAnalogClerk),intent(in)            :: self
    integer(shortInt),dimension(:),allocatable :: validCodes

    validCodes = [ spawn_CODE, outColl_CODE, cycleEnd_CODE ]

  end function validReports

  !!
  !! Return memory size of the clerk
  !!
  !! See tallyClerk_inter for details
  !!
  elemental function getSize(self) result(S)
    class(nuAnalogClerk), intent(in) :: self
    integer(shortInt)                :: S

    S = MEM_SIZE
    if (allocated(self % map)) S = S * self % map % bins(0)

  end function getSize

  !!
  !! Process fission report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportSpawn(self, MT, pOld, pNew, xsData, mem)
    class(nuAnalogClerk), intent(inout)   :: self
    integer(shortInt), intent(in)         :: MT
    class(particle), intent(in)           :: pOld
    class(particleState), intent(in)      :: pNew
    class(nuclearDatabase), intent(inout) :: xsData
    type(scoreMemory), intent(inout)      :: mem
    type(particleState)                   :: state
    integer(shortInt)                     :: binIdx
    integer(longInt)                      :: addr

    ! Get current particle state
    state = pOld

    ! Find bin index
    if (allocated(self % map)) then
      binIdx = self % map % map(state)
    else
      binIdx = 1
    end if

    ! Return if invalid bin index
    if (binIdx == 0) return

    ! Calculate bin address
    addr = self % getMemAddress() + MEM_SIZE * (binIdx - 1)

    if ((self % delayed .and. pNew % type == P_PRECURSOR) .or. &
        & (.not. self % delayed .and. pNew % type == P_NEUTRON)) then

      ! Score nu
      call mem % score(ONE, addr + NU_FISS)

    end if

  end subroutine reportSpawn

  !!
  !! Process outgoing collision report
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportOutColl(self, p, MT, muL, xsData, mem)
    class(nuAnalogClerk), intent(inout)  :: self
    class(particle), intent(in)          :: p
    integer(shortInt), intent(in)        :: MT
    real(defReal), intent(in)            :: muL
    class(nuclearDatabase),intent(inout) :: xsData
    type(scoreMemory), intent(inout)     :: mem
    type(particleState)                  :: state
    integer(shortInt)                    :: binIdx
    integer(longInt)                     :: addr

    ! Score analog event if a fission happened
    if (MT == N_FISSION) then

      ! Get current particle state
      state = p

      ! Find bin index
      if (allocated(self % map)) then
        binIdx = self % map % map(state)
      else
        binIdx = 1
      end if

      ! Return if invalid bin index
      if (binIdx == 0) return

      ! Calculate bin address
      addr = self % getMemAddress() + MEM_SIZE * (binIdx - 1)

      call mem % score(ONE, addr + EVENT)

    end if

  end subroutine reportOutColl

  !!
  !! Process end of the cycle
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine reportCycleEnd(self, end, mem)
    class(nuAnalogClerk), intent(inout) :: self
    class(particleDungeon), intent(in)  :: end
    type(scoreMemory), intent(inout)    :: mem
    integer(longInt)                    :: addr
    integer(shortInt)                   :: i, N
    real(defReal)                       :: nu, events, nu_score

    ! Score average number of neutrons produced by fission
    if (mem % lastCycle()) then

      if (allocated(self % map)) then
        N = self % map % bins(0)
      else
        N = 1
      end if

      ! Loop over map bins
      do i = 1, N

        ! Calculate bin address
        addr = self % getMemAddress() + MEM_SIZE * (i - 1)

        nu     = mem % getScore(addr + NU_FISS)
        events = mem % getScore(addr + EVENT)

        ! Calculate average nu
        if (events == ZERO) then
          nu_score = ZERO
        else
          nu_score = nu/events
        end if

        call mem % accumulate(nu_score, addr + NU_MEAN)

      end do

    end if

  end subroutine reportCycleEnd

  !!
  !! Display convergance progress on the console
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine display(self, mem)
    class(nuAnalogClerk), intent(in)  :: self
    type(scoreMemory), intent(in)     :: mem

    ! Does nothing

  end subroutine display

  !!
  !! Write contents of the clerk in the slot to output file
  !!
  !! See tallyClerk_inter for details
  !!
  subroutine print(self, outFile, mem)
    class(nuAnalogClerk), intent(in) :: self
    class(outputFile), intent(inout) :: outFile
    type(scoreMemory), intent(in)    :: mem
    real(defReal)                    :: nu, STD
    character(nameLen)               :: name
    integer(shortInt)                :: i
    integer(longInt)                 :: addr
    integer(shortInt),dimension(:),allocatable :: resArrayShape

    ! Print to output file
    call outFile % startBlock(self % getName())

    ! Print out map info
    if (allocated(self % map)) then
      call self % map % print(outFile)
      resArrayShape = self % map % binArrayShape()
    else
      resArrayShape = [1]
    end if

    name = 'Res'

    call outFile % startArray(name, resArrayShape)

    ! Print results to the file
    do i = 1, product(resArrayShape)
      addr = self % getMemAddress() + MEM_SIZE * (i - 1)
      call mem % getResult(nu, STD, addr + NU_MEAN)
      call outFile % addResult(nu, STD)
    end do

    call outFile % endArray()

    call outFile % endBlock()

  end subroutine print

end module nuAnalogClerk_class
