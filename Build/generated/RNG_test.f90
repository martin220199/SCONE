module RNG_test

  use numPrecision
  use RNG_class, only : RNG
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test random numbers
  !!
!@Test
  subroutine testRN()
    integer(shortInt), parameter :: N = 1000
    type(RNG)                    :: pRNG
    real(defReal), dimension(N)  :: rand
    integer(shortInt)            :: i

    call pRNG % init(int(z'5c3a84c9', longInt))

    do i=1,N
      rand(i) = pRNG % get()
    end do

    ! Check correcness
#line 28 "/home/mskrette/SCONE_cambridge_fork/SCONE/RandomNumbers/Tests/RNG_test.f90"
  call assertGreaterThanOrEqual(ONE, rand, &
 & location=SourceLocation( &
 & 'RNG_test.f90', &
 & 28) )
  if (anyExceptions()) return
#line 29 "/home/mskrette/SCONE_cambridge_fork/SCONE/RandomNumbers/Tests/RNG_test.f90"
#line 29 "/home/mskrette/SCONE_cambridge_fork/SCONE/RandomNumbers/Tests/RNG_test.f90"
  call assertLessThanOrEqual(ZERO, rand, &
 & location=SourceLocation( &
 & 'RNG_test.f90', &
 & 29) )
  if (anyExceptions()) return
#line 30 "/home/mskrette/SCONE_cambridge_fork/SCONE/RandomNumbers/Tests/RNG_test.f90"

  end subroutine testRN

  !!
  !! Test skip forward and backwards
  !!
!@Test
  subroutine testSkip()
    type(RNG)         :: rand1
    type(RNG)         :: rand2
    integer(longInt)  :: seed
    real(defReal)     :: r_start, r2_start, r_end, r2_end
    integer(shortInt) :: i, N

    !! Initialise both RNGs to a nice number
    seed = int(z'5c3a84c9', longInt)
    call rand1 % init(seed)
    call rand2 % init(seed)

    !! Get initial random number
    r_start = rand1 % get()

    !! Move forward by 13456757 steps
    N = 13456757
    do i=1,N
      r_end = rand1 % get()
    end do

    ! Skip 2nd generator forward
    call rand2 % skip(int(N, longInt))
    r2_end = rand2 % get()
    
    ! Skip 2nd generator backwards. Must be 1 more becouse we drew a RN from generator
    call rand2 % skip(-int(N + 1, longInt))
    r2_start = rand2 % get()

    ! Verify values
#line 67 "/home/mskrette/SCONE_cambridge_fork/SCONE/RandomNumbers/Tests/RNG_test.f90"
  call assertEqual(r_end, r2_end, &
 & location=SourceLocation( &
 & 'RNG_test.f90', &
 & 67) )
  if (anyExceptions()) return
#line 68 "/home/mskrette/SCONE_cambridge_fork/SCONE/RandomNumbers/Tests/RNG_test.f90"
#line 68 "/home/mskrette/SCONE_cambridge_fork/SCONE/RandomNumbers/Tests/RNG_test.f90"
  call assertEqual(r_start, r2_start, &
 & location=SourceLocation( &
 & 'RNG_test.f90', &
 & 68) )
  if (anyExceptions()) return
#line 69 "/home/mskrette/SCONE_cambridge_fork/SCONE/RandomNumbers/Tests/RNG_test.f90"

  end subroutine testSkip



end module RNG_test

module WrapRNG_test
   use pFUnit_mod
   use RNG_test
   implicit none
   private

contains


end module WrapRNG_test

function RNG_test_suite() result(suite)
   use pFUnit_mod
   use RNG_test
   use WrapRNG_test
   type (TestSuite) :: suite

   suite = newTestSuite('RNG_test_suite')

   call suite%addTest(newTestMethod('testRN', testRN))

   call suite%addTest(newTestMethod('testSkip', testSkip))


end function RNG_test_suite

