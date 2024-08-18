module trackClerk_test

  use numPrecision
  use genericProcedures,              only : numToChar
  use trackClerk_class,               only :trackClerk
  use particle_class,                 only : particle
  use dictionary_class,               only : dictionary
  use scoreMemory_class,              only : scoreMemory
  use testNeutronDatabase_class,      only : testNeutronDatabase
  use outputFile_class,               only : outputFile
  use pFUnit_mod

  implicit none

!@testParameter(constructor = new_testNumber)
  type, extends(AbstractTestParameter) :: testNumber
    integer(shortInt) :: i
  contains
    procedure :: toString
  end type testNumber

!@testCase(constructor=newTest)
  type, extends(ParameterizedTestCase) :: test_trackClerk
    private
    logical(defBool)                            :: hasFilter
    logical(defBool)                            :: hasMap
    logical(defBool)                            :: has2Res
    integer(longInt), dimension(:), allocatable :: bins
    real(defReal), dimension(:), allocatable    :: results
  end type test_trackClerk

contains

  !!
  !! Build new test parameter form integer
  !!
  function new_testNumber(i) result (tstNum)
    integer(shortInt) :: i
    type(testNumber)  :: tstNum

    tstNum % i = i

  end function new_testNumber

  !!
  !! Write test parameter to string
  !!
  function toString(this) result(string)
    class(testNumber), intent(in) :: this
    character(:), allocatable :: string
    character(nameLen)        :: str

    write (str,*) this % i
    string = str

  end function toString

  !!
  !! Construct test case
  !!
  function newTest(testParam) result(tst)
    type(testNumber), intent(in) :: testParam
    type(test_trackClerk)        :: tst
    integer(shortInt)            :: testNum, Nbins, i
    real(defReal)                :: score1, score2

    testNum = testParam % i - 1

    ! Set test parameters
    tst % hasFilter = mod(testNum, 2)   == 1
    tst % hasMap    = mod(testNum/2, 2) == 1
    tst % has2Res   = mod(testNum/4, 2) == 1

    ! Load expected results

    ! Determine number of bins
    Nbins = 1
    if(tst % has2Res) Nbins = Nbins * 2
    if(tst % hasMap)  Nbins = Nbins * 7

    ! Allocate result arrays
    tst % bins    = [(int(i,longInt), i=1,Nbins)]
    tst % results = [(ZERO, i=1,Nbins)]

    ! Set approperiate results (wgt * L)
    score1 = 0.7_defReal * 0.3_defReal
    score2 = 1.3_defReal * 0.3_defReal

    select case(testParam % i)
      case(1) ! Single Bin, both score to the same
        tst % results(1) = score1 + score2

      case(2) ! Single bin with filter
        tst % results(1) = score1

      case(3) ! Map without filter
        tst % results(1) = score1
        tst % results(6) = score2

      case(4) ! Map with filter
        tst % results(1) = score1

      case(5) ! Single bin with 2nd Response
        tst % results(1) = score1 + score2
        tst % results(2) = score1 * 1.3_defReal + score2 * 1.3_defReal

      case(6) ! Single bin with filter and 2nd Response
        tst % results(1) = score1
        tst % results(2) = score1 * 1.3_defReal

      case(7) ! Map with 2nd response
        tst % results(1)  = score1
        tst % results(2)  = score1 * 1.3_defReal
        tst % results(11) = score2
        tst % results(12) = score2 * 1.3_defReal

      case(8) ! Map with 2nd response and filter
        tst % results(1)  = score1
        tst % results(2)  = score1 * 1.3_defReal
    end select
  end function newTest

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Scoring test
  !!
!@Test(cases=[1,2,3,4,5,6,7,8])
  subroutine testScoring(this)
    class(test_trackClerk), intent(inout)     :: this
    logical(defBool)                          :: hasFilter, hasMap, has2Res
    character(:),allocatable                  :: case
    type(trackClerk)                          :: clerk
    type(scoreMemory)                         :: mem
    type(particle)                            :: p
    type(testNeutronDatabase)                 :: nucData
    type(outputFile)                          :: outF
    type(dictionary)                          :: filterDict, mapDict, res1Dict, res2Dict, clerkDict
    character(nameLen)                        :: res1Name, res2Name, clerkName
    integer(shortInt)                         :: i
    real(defReal)                             :: res, L
    real(defReal), parameter :: TOL = 1.0E-9

    ! Copy test settings
    hasFilter = this % hasFilter
    hasMap    = this % hasMap
    has2Res   = this % has2Res

    ! Build case description
    case = 'Vanila case with: '
    if(hasFilter) case = case // ' Filter '
    if(hasMap)    case = case // ' Map '
    if(has2Res)   case = case // ' 2nd Response '

    ! Define filter dictionary
    call filterDict % init(3)
    call filterDict % store('type','testFilter')
    call filterDict % store('minIdx',0)
    call filterDict % store('maxIdx',5)

    ! Define Map dictionary
    call mapDict % init(2)
    call mapDict % store('type','testMap')
    call mapDict % store('maxIdx',7)

    ! Define 1st response dictionary and name
    res1Name = 'flux'
    call res1Dict % init(1)
    call res1Dict % store('type','fluxResponse')

    ! Define 2nd response dictionary and name
    res2Name ='testResponse'
    call res2Dict % init(2)
    call res2Dict % store('type','testResponse')
    call res2Dict % store('value', 1.3_defReal)

    ! Configure dictionary for the clerk
    call clerkDict % init(6)
    call clerkDict % store('type','trackClerk')
    call clerkDict % store(res1Name, res1Dict)
    call clerkDict % store(res2Name, res2Dict)

    ! Store filter or map
    if(hasFilter) call clerkDict % store('filter', filterDict)
    if(hasMap)    call clerkDict % store('map', mapDict)

    ! Store responses used
    if(has2Res) then
      call clerkDict % store('response', [res1Name, res2Name])
    else
      call clerkDict % store('response', [res1Name])
    end if

    ! Build Clerk
    clerkName ='myClerk'
    call clerk % init(clerkDict, clerkName)

    ! Create score memory
    call mem % init(int(clerk % getSize(), longInt) , 1)
    call clerk % setMemAddress(1_longInt)

    ! Set track lenght
    L = 0.3_defReal

    ! Perform scoring
    p % prePath % matIdx = 1
    p % w = 0.7_defReal
    call clerk % reportPath(p, L, nucData, mem)

    p % prePath % matIdx = 6
    p % w = 1.3_defReal
    call clerk % reportPath(p, L, nucData, mem)

    call mem % closeCycle(ONE)

    ! Verify results of scoring
    do i=1,size(this % bins)
      call mem % getResult(res, this % bins(i))
#line 221 "/home/mskrette/SCONE_cambridge_fork/SCONE/Tallies/TallyClerks/Tests/trackClerk_test.f90"
  call assertEqual(this % results(i), res, TOL, case // 'BIN : ' //numToChar(i) , &
 & location=SourceLocation( &
 & 'trackClerk_test.f90', &
 & 221) )
  if (anyExceptions()) return
#line 222 "/home/mskrette/SCONE_cambridge_fork/SCONE/Tallies/TallyClerks/Tests/trackClerk_test.f90"
    end do

    ! Verify that size of memory returned is correct
#line 225 "/home/mskrette/SCONE_cambridge_fork/SCONE/Tallies/TallyClerks/Tests/trackClerk_test.f90"
  call assertEqual(size(this % bins), clerk % getSize(), case // 'Memory size test:', &
 & location=SourceLocation( &
 & 'trackClerk_test.f90', &
 & 225) )
  if (anyExceptions()) return
#line 226 "/home/mskrette/SCONE_cambridge_fork/SCONE/Tallies/TallyClerks/Tests/trackClerk_test.f90"

    ! Verify that output calls are correct
    call outF % init('dummyPrinter', fatalErrors = .false.)
    call clerk % print (outF, mem)

#line 231 "/home/mskrette/SCONE_cambridge_fork/SCONE/Tallies/TallyClerks/Tests/trackClerk_test.f90"
  call assertTrue(outF % isValid(), case, &
 & location=SourceLocation( &
 & 'trackClerk_test.f90', &
 & 231) )
  if (anyExceptions()) return
#line 232 "/home/mskrette/SCONE_cambridge_fork/SCONE/Tallies/TallyClerks/Tests/trackClerk_test.f90"

    ! Clean up
    call nucData % kill()
    call clerkDict % kill()
    call filterDict % kill()
    call mapDict % kill()
    call res1Dict % kill()
    call res2Dict % kill()

  end subroutine testScoring

end module trackClerk_test

module WraptrackClerk_test
   use pFUnit_mod
   use trackClerk_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_trackClerk) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use trackClerk_test
        class (test_trackClerk), intent(inout) :: this
     end subroutine userTestMethod
   end interface

contains

   subroutine runMethod(this)
      class (WrapUserTestCase), intent(inout) :: this

      call this%testMethodPtr(this)
   end subroutine runMethod

   function makeCustomTest(methodName, testMethod, testParameter) result(aTest)
#ifdef INTEL_13
      use pfunit_mod, only: testCase
#endif
      type (WrapUserTestCase) :: aTest
#ifdef INTEL_13
      target :: aTest
      class (WrapUserTestCase), pointer :: p
#endif
      character(len=*), intent(in) :: methodName
      procedure(userTestMethod) :: testMethod
      type (testNumber), intent(in) :: testParameter
      aTest%test_trackClerk = newTest(testParameter)

      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
      call aTest%setTestParameter(testParameter)
   end function makeCustomTest

end module WraptrackClerk_test

function trackClerk_test_suite() result(suite)
   use pFUnit_mod
   use trackClerk_test
   use WraptrackClerk_test
   type (TestSuite) :: suite

   type (testNumber), allocatable :: testParameters(:)
   type (testNumber) :: testParameter
   integer :: iParam 
   integer, allocatable :: cases(:) 
 
   suite = newTestSuite('trackClerk_test_suite')

   cases = [1,2,3,4,5,6,7,8]
   testParameters = [(new_testNumber(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testScoring', testScoring, testParameter))
   end do


end function trackClerk_test_suite

