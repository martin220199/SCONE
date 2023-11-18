module multipleCorrelatedLaws_test

  use numPrecision
  use endfConstants
  use RNG_class,                    only : RNG
  use correlatedLawENDF_inter,      only : correlatedLawENDF
  use multipleCorrelatedLaws_class, only : multipleCorrelatedLaws
  use testCorrelatedLaw_class,      only : testCorrelatedLaw
  use pFUnit_mod
  implicit none


  ! Fractional probability data taken from Fluorine-19 MT=16 JEFF 3.1.1
  ! Constant 0.5 fraction. Same for each law
  real(defReal), dimension(*), parameter     :: eGrid = [10.9870_defReal, 20.0_defReal]
  real(defReal), dimension(*), parameter     :: prob  = [0.5_defReal, 0.5_defReal]
  integer(shortInt), dimension(*), parameter :: bounds = [2]
  integer(shortInt), dimension(*), parameter :: inter  = [histogramInterpolation]

  ! Test laws settings
  real(defReal), dimension(2), parameter :: E_outs = [8.7_defReal, 0.3_defReal]
  real(defReal), dimension(2), parameter :: mu_outs = [0.77_defReal, -0.21_defReal]

  ! Test objects
  type(multipleCorrelatedLaws) :: law

contains


  !!
  !! Setup test enviroment
  !!
!@Before
  subroutine setUp()
    type(testCorrelatedLaw)               :: tempLaw
    class(correlatedLawENDF), allocatable :: corrLaw

    ! Allocate space
    call law % init(2)

    ! Load correlated laws
    tempLaw % E_out = E_outs(1)
    tempLaw % mu    = mu_outs(1)
    allocate(corrLaw, source = tempLaw)
    call law % addLaw(corrLaw, eGrid, prob, bounds, inter)

    tempLaw % E_out = E_outs(2)
    tempLaw % mu    = mu_outs(2)
    allocate(corrLaw, source = tempLaw)
    call law % addLaw(corrLaw, eGrid, prob, bounds, inter)

  end subroutine setUp

  !!
  !! Clean test enviroment
  !!
!@After
  subroutine cleanUp()

    ! Clean
    call law % kill()

  end subroutine cleanUp


  !!
  !! Test getting probability of sample
  !!
!@Test
  subroutine test_probabilityOf()
    real(defReal), parameter :: TOL = 1.0E-7_defReal

    ! Immpossible combination
#line 74 "/home/mskrette/SCONE_cambridge_fork/SCONE/NuclearData/emissionENDF/correlatedLawENDF/Tests/multipleCorrelatedLaws_test.f90"
  call assertEqual(ZERO, law % probabilityOf(mu_outs(1), E_outs(2), 11.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'multipleCorrelatedLaws_test.f90', &
 & 74) )
  if (anyExceptions()) return
#line 75 "/home/mskrette/SCONE_cambridge_fork/SCONE/NuclearData/emissionENDF/correlatedLawENDF/Tests/multipleCorrelatedLaws_test.f90"

    ! Probability with energy slightly below the range
#line 77 "/home/mskrette/SCONE_cambridge_fork/SCONE/NuclearData/emissionENDF/correlatedLawENDF/Tests/multipleCorrelatedLaws_test.f90"
  call assertEqual(0.5_defReal, law % probabilityOf(mu_outs(1), E_outs(1), eGrid(1)*0.999_defReal), TOL, &
 & location=SourceLocation( &
 & 'multipleCorrelatedLaws_test.f90', &
 & 77) )
  if (anyExceptions()) return
#line 78 "/home/mskrette/SCONE_cambridge_fork/SCONE/NuclearData/emissionENDF/correlatedLawENDF/Tests/multipleCorrelatedLaws_test.f90"

    ! Probability with incident energy a bit above the range
#line 80 "/home/mskrette/SCONE_cambridge_fork/SCONE/NuclearData/emissionENDF/correlatedLawENDF/Tests/multipleCorrelatedLaws_test.f90"
  call assertEqual(0.5_defReal, law % probabilityOf(mu_outs(2), E_outs(2), eGrid(2)*1.001_defReal), TOL, &
 & location=SourceLocation( &
 & 'multipleCorrelatedLaws_test.f90', &
 & 80) )
  if (anyExceptions()) return
#line 81 "/home/mskrette/SCONE_cambridge_fork/SCONE/NuclearData/emissionENDF/correlatedLawENDF/Tests/multipleCorrelatedLaws_test.f90"

    ! Normal enquiry
#line 83 "/home/mskrette/SCONE_cambridge_fork/SCONE/NuclearData/emissionENDF/correlatedLawENDF/Tests/multipleCorrelatedLaws_test.f90"
  call assertEqual(0.5_defReal, law % probabilityOf(mu_outs(2), E_outs(2), 17.0_defReal), TOL, &
 & location=SourceLocation( &
 & 'multipleCorrelatedLaws_test.f90', &
 & 83) )
  if (anyExceptions()) return
#line 84 "/home/mskrette/SCONE_cambridge_fork/SCONE/NuclearData/emissionENDF/correlatedLawENDF/Tests/multipleCorrelatedLaws_test.f90"

  end subroutine test_probabilityOf

  !!
  !! Test sampling
  !!
  !! Easy to do since it is just a binary choice
  !! USING THE FACT THAT LAW PROBABILITY IS UNIFORM IN INCIDENT ENERGY
  !!
!@Test
  subroutine test_sample()
    real(defReal), dimension(1200,2) :: samples
    type(RNG)                        :: rand
    real(defReal)                    :: E_in
    integer(shortInt)                :: i
    real(defReal)                    :: C1, C2, Chi

    ! Initialise RNG
    call rand % init(1945_longInt)

    ! Draw samples with energy slightly below
    E_in =  eGrid(1)*0.999_defReal
    do i = 1, 400
      call law % sample(samples(i,1), samples(i,2), E_in, rand)
    end do

    ! Draw slightly above the range
    E_in = eGrid(2)*1.001_defReal
    do i = 401, 800
      call law % sample(samples(i,1), samples(i,2), E_in, rand)
    end do

    ! Drow from the middle of the range
    E_in = 17.0_defReal
    do i = 801, 1200
      call law % sample(samples(i,1), samples(i,2), E_in, rand)
    end do

    ! Count occurences of each law
    C1 = real(count(samples(:,1) == mu_outs(1) .and. samples(:,2) == E_outs(1)), defReal)
    C2 = real(count(samples(:,1) == mu_outs(2) .and. samples(:,2) == E_outs(2)), defReal)

    Chi = (C1 - 600.0_defReal)**2/600_defReal + (C2 - 600.0_defReal)**2/600_defReal

    ! Value of 5.02 was chosen from 1-degree of freedom Chi-sq distribution to obtain
    ! 97.5% confidence interval
#line 130 "/home/mskrette/SCONE_cambridge_fork/SCONE/NuclearData/emissionENDF/correlatedLawENDF/Tests/multipleCorrelatedLaws_test.f90"
  call assertGreaterThan(5.02_defReal, Chi, &
 & location=SourceLocation( &
 & 'multipleCorrelatedLaws_test.f90', &
 & 130) )
  if (anyExceptions()) return
#line 131 "/home/mskrette/SCONE_cambridge_fork/SCONE/NuclearData/emissionENDF/correlatedLawENDF/Tests/multipleCorrelatedLaws_test.f90"

  end subroutine test_sample


end module multipleCorrelatedLaws_test

module WrapmultipleCorrelatedLaws_test
   use pFUnit_mod
   use multipleCorrelatedLaws_test
   implicit none
   private

contains


end module WrapmultipleCorrelatedLaws_test

function multipleCorrelatedLaws_test_suite() result(suite)
   use pFUnit_mod
   use multipleCorrelatedLaws_test
   use WrapmultipleCorrelatedLaws_test
   type (TestSuite) :: suite

   suite = newTestSuite('multipleCorrelatedLaws_test_suite')

   call suite%addTest(newTestMethod('test_probabilityOf', test_probabilityOf, setUp, cleanUp))

   call suite%addTest(newTestMethod('test_sample', test_sample, setUp, cleanUp))


end function multipleCorrelatedLaws_test_suite

