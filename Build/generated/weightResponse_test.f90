module weightResponse_test

  use numPrecision
  use endfConstants
  use weightResponse_class,           only : weightResponse
  use particle_class,                 only : particle, P_NEUTRON
  use dictionary_class,               only : dictionary
  use testNeutronDatabase_class,      only : testNeutronDatabase
  use pFUnit_mod

  implicit none

!@testCase
  type, extends(TestCase) :: test_weightResponse
    private
    type(weightResponse)        :: response_weight_m0
    type(weightResponse)        :: response_weight_m2
    type(testNeutronDatabase)   :: xsData
  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_weightResponse


contains

  !!
  !! Sets up test_macroResponse object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_weightResponse), intent(inout) :: this
    type(dictionary)                          :: tempDict

    ! Cross-sections:         Total
    call this % xsData % build(4.0_defReal)

    ! Set up weight response
    call tempDict % init(1)
    call tempDict % store('moment', 0)
    call this % response_weight_m0 % init(tempDict)
    call tempDict % kill()

    call tempDict % init(1)
    call tempDict % store('moment', 2)
    call this % response_weight_m2 % init(tempDict)
    call tempDict % kill()

  end subroutine setUp

  !!
  !! Kills test_weightResponse object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_weightResponse), intent(inout) :: this

    ! Kill and deallocate testTransportNuclearData
    call this % xsData % kill()

  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test correct behaviour of the filter
  !!
!@Test
  subroutine testGettingResponse(this)
    class(test_weightResponse), intent(inout) :: this
    type(particle)                            :: p
    real(defReal), parameter                  :: TOL = 1.0E-9

    p % type = P_NEUTRON
    p % w = 2.0_defReal

    ! Test response values
#line 78 "/home/mskrette/SCONE_cambridge_fork/SCONE/Tallies/TallyResponses/Tests/weightResponse_test.f90"
  call assertEqual(2.0_defReal, this % response_weight_m0 % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'weightResponse_test.f90', &
 & 78) )
  if (anyExceptions()) return
#line 79 "/home/mskrette/SCONE_cambridge_fork/SCONE/Tallies/TallyResponses/Tests/weightResponse_test.f90"
#line 79 "/home/mskrette/SCONE_cambridge_fork/SCONE/Tallies/TallyResponses/Tests/weightResponse_test.f90"
  call assertEqual(8.0_defReal, this % response_weight_m2 % get(p, this % xsData), TOL, &
 & location=SourceLocation( &
 & 'weightResponse_test.f90', &
 & 79) )
  if (anyExceptions()) return
#line 80 "/home/mskrette/SCONE_cambridge_fork/SCONE/Tallies/TallyResponses/Tests/weightResponse_test.f90"

  end subroutine testGettingResponse

end module weightResponse_test

module WrapweightResponse_test
   use pFUnit_mod
   use weightResponse_test
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(test_weightResponse) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use weightResponse_test
        class (test_weightResponse), intent(inout) :: this
     end subroutine userTestMethod
   end interface

contains

   subroutine runMethod(this)
      class (WrapUserTestCase), intent(inout) :: this

      call this%testMethodPtr(this)
   end subroutine runMethod

   function makeCustomTest(methodName, testMethod) result(aTest)
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
      aTest%testMethodPtr => testMethod
#ifdef INTEL_13
      p => aTest
      call p%setName(methodName)
#else
      call aTest%setName(methodName)
#endif
   end function makeCustomTest

end module WrapweightResponse_test

function weightResponse_test_suite() result(suite)
   use pFUnit_mod
   use weightResponse_test
   use WrapweightResponse_test
   type (TestSuite) :: suite

   suite = newTestSuite('weightResponse_test_suite')

   call suite%addTest(makeCustomTest('testGettingResponse', testGettingResponse))


end function weightResponse_test_suite

