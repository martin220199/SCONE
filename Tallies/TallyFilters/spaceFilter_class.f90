module spaceFilter_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use particle_class,    only : particleState
  use dictionary_class,  only : dictionary
  use tallyFilter_inter, only : tallyFilter

  implicit none
  private

  !!
  !! Filter that tests if energy of a state is in a closed, energy (or group) interval.
  !!
  !! Private members:
  !!   Emin -> minimum value of energy
  !!   Emax -> maximum value of energy
  !!   Gtop -> Top energy group
  !!   Glow -> Bottom of allowable energy range
  !!
  !! Interface:
  !!   tallyFilter Interface
  !!   build -> build filter from components
  !!
  !! NOTE: Energy is not enforced to be +ve. Could be useful in debugging
  !!
  !! Sample Dictionary Input:
  !!
  !!   CEfilter {
  !!     type energyFilter;
  !!     Emin 1.0;
  !!     Emax 2.0;
  !!   }
  !!
  !!   MGfilter {
  !!     type energyFilter;
  !!     Gtop 3;
  !!     Glow 17;
  !!  }
  !!
  !!  filter {
  !!     type energyFilter;
  !!     Emin 1.0;
  !!     Emax 2.0;
  !!     Gtop 3;
  !!     Glow 17;
  !! }
  !!
  type, public,extends(tallyFilter) :: spaceFilter
    private
    real(defReal)     :: Xmin
    real(defReal)     :: Xmax
    real(defReal)     :: Ymin
    real(defReal)     :: Ymax
  contains
    procedure :: init
    procedure :: isPass

    !! Instance specific procedures
    generic :: build => build_CE, build_MG, build_CEMG
    procedure :: build_CE
    procedure :: build_MG
    procedure :: build_CEMG

  end type spaceFilter

contains

  !!
  !! Initialise energyFilter from dictionary
  !!
  subroutine init(self,dict)
    class(spaceFilter), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    real(defReal)                      :: X1, X2, Y1, Y2
    logical(defBool)                   :: hasMG, hasCE

    ! Detect which case
    hasMG = dict % isPresent('Glow') .or. dict % isPresent('Gtop')
    hasCE = dict % isPresent('Emin') .or. dict % isPresent('Emax')

    if(hasMG .and. hasCE) then
      ! CE-MG case
      call dict % get(X1,'Xmin')
      call dict % get(X2,'Xmax')
      call dict % get(Y1,'Ymin')
      call dict % get(Y2,'Ymax')

      call self % build(X1, X2, Y1, Y2)

    else if(hasMG) then
      ! MG Case
      call dict % get(X1,'Xmin')
      call dict % get(X2,'Xmax')
      call dict % get(Y1,'Ymin')
      call dict % get(Y2,'Ymax')

      call self % build(X1, X2, Y1, Y2)
    else
      ! CE Case
      call dict % get(X1,'Xmin')
      call dict % get(X2,'Xmax')
      call dict % get(Y1,'Ymin')
      call dict % get(Y2,'Ymax')

      call self % build(X1, X2, Y1, Y2)
    end if

  end subroutine init

  !!
  !! Returns true if energy value is between specified bounds
  !!
  elemental function isPass(self,state) result(passed)
    class(spaceFilter), intent(in)  :: self
    class(particleState), intent(in) :: state
    logical(defBool)                 :: passed
    real(defReal)                    :: X, Y

    ! MG particle
    if(state % isMG) then
      ! CE paricle
      X = state % r(1)
      Y = state % r(2)
      passed = (self % Xmin <= X) .and. (X <= self % Xmax) .and. (self % Ymin <= Y) .and. (Y <= self % Ymax)

    else
      ! CE paricle
      X = state % r(1)
      Y = state % r(2)
      passed = (self % Xmin <= X) .and. (X <= self % Xmax) .and. (self % Ymin <= Y) .and. (Y <= self % Ymax)

    end if
  end function isPass

  !!
  !! Build energyFilter for CE particles only from components
  !!
  !! Args:
  !!   Emin [in] -> minimum energy [MeV]
  !!   Emax [in] -> maximum energy [MeV]
  !!
  !! Errors:
  !!   fatalError if Emin > Emax
  !!
  subroutine build_CE(self, Emin, Emax)
    class(spaceFilter), intent(inout) :: self
    real(defReal), intent(in)          :: Emin
    real(defReal), intent(in)          :: Emax
    character(100), parameter :: Here = 'build_CE (spaceFilter_class.f90)'


  end subroutine build_CE
    
  !!
  !! Build energyFilter for MG particles only from components
  !!
  !! Args:
  !!   Gtop [in] -> maximum energy (lowest index) energy group
  !!   Glow [in] -> minimum energy (highest index) enegy group
  !!
  !! Errors:
  !!   fatalError if Gtop > Elow
  !!
  subroutine build_MG(self, Gtop, Glow)
    class(spaceFilter), intent(inout) :: self
    integer(shortInt), intent(in)      :: Gtop
    integer(shortInt), intent(in)      :: Glow
    character(100), parameter :: Here = 'build_MG (spaceFilter_class.f90)'


  end subroutine build_MG

  !!
  !! Build energyFilter for MG and CE particles from components
  !!
  !! Args:
  !!   Emin [in] -> minimum energy [MeV]
  !!   Emax [in] -> maximum energy [MeV]
  !!   Gtop [in] -> maximum energy (lowest index) energy group
  !!   Glow [in] -> minimum energy (highest index) enegy group
  !!
  !! Errors:
  !!   fatalError if Emin > Emax
  !!   fatalError if Gtop > Elow
  !!
  subroutine build_CEMG(self, Xmin, Xmax, Ymin, Ymax)
    class(spaceFilter), intent(inout) :: self
    real(defReal), intent(in)          :: Xmin
    real(defReal), intent(in)          :: Xmax
    real(defReal), intent(in)          :: Ymin
    real(defReal), intent(in)          :: Ymax
    character(100), parameter :: Here = 'build_CEMG (spaceFilter_class.f90)'



    self % Xmin = Xmin
    self % Xmax = Xmax

    self % Ymin = Ymin
    self % Ymax = Ymax

    ! Verify bounds
    if( self % Xmax <= self % Xmin) then
      call fatalError(Here,'Emin='// numToChar(self % Xmin) //' is larger or equal to Emax=' // numToChar(self % Xmax))
    end if

    if( self % Ymax <= self % Ymin) then
      call fatalError(Here,'Emin='// numToChar(self % Ymin) //' is larger or equal to Emax=' // numToChar(self % Ymax))
    end if

  end subroutine build_CEMG


end module spaceFilter_class
