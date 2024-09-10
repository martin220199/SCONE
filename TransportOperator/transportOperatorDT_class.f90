!!
!! Transport operator for delta tracking
!!
module transportOperatorDT_class
  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar
  use particle_class,             only : particle, particleState
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary
  use rng_class,                  only : rng

  ! Superclass
  use transportOperator_inter,    only : transportOperator

  ! Geometry interfaces
  use geometry_inter,             only : geometry

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase

  implicit none
  private

  !!
  !! Transport operator that moves a particle with delta tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorDT
  contains
    procedure :: transit => deltaTracking
    procedure :: processParticlePrediction
  end type transportOperatorDT

contains

  subroutine deltaTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorDT), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon), intent(inout)     :: thisCycle
    class(particleDungeon), intent(inout)     :: nextCycle
    real(defReal)                             :: majorant_inv, sigmaT, distance
    character(100), parameter :: Here = 'deltaTracking (transportOperatorDT_class.f90)'

    ! Get majornat XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getMajorantXS(p)

    ! Should never happen! Prevents Inf distances
    if (abs(majorant_inv) > huge(majorant_inv)) call fatalError(Here, "Majorant is 0")

    DTLoop:do
      distance = -log( p% pRNG % get() ) * majorant_inv


      if (p % time + distance / p % getSpeed() > p % timeMax) then
        distance = distance * (p % timeMax - p % time)/(distance / p % getSpeed())
        p % fate = AGED_FATE
        p % time = p % timeMax
        call self % geom % teleport(p % coords, distance)
        return
      endif

      ! Move partice in the geometry
      call self % geom % teleport(p % coords, distance)

      ! Update time
      p % time = p % time + distance / p % getSpeed()

      ! If particle has leaked, exit
      if (p % matIdx() == OUTSIDE_FILL) then
        p % fate = LEAK_FATE
        p % isDead = .true.
        return
      end if

      ! Check for void
      if(p % matIdx() == VOID_MAT) then
        call tally % reportInColl(p, .true.)
        cycle DTLoop
      end if

      ! Give error if the particle somehow ended in an undefined material
      if (p % matIdx() == UNDEF_MAT) then
        print *, p % rGlobal()
        call fatalError(Here, "Particle is in undefined material")
      end if

      ! Obtain the local cross-section
      sigmaT = self % xsData % getTransMatXS(p, p % matIdx())

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real, report collision if virtual
      if (p % pRNG % get() < sigmaT*majorant_inv) then
        exit DTLoop
      else
        call tally % reportInColl(p, .true.)
      end if

    end do DTLoop

    call tally % reportTrans(p)
  end subroutine deltaTracking

  subroutine processParticlePrediction(self, p, p_p, maxT)
    class(transportOperatorDT), intent(inout) :: self
    type(particle), intent(in)                :: p
    type(particle), intent(out)               :: p_p
    real(defReal), intent(in)                 :: maxT
    type(particleState)                       :: temp_p
    real(defReal)                             :: distance
    character(100), parameter :: Here = 'processParticlePrediction (transportOperatorDT_class.f90)'

    temp_p = p
    p_p = temp_p
    distance = ONE / self % xsData % getMajorantXS(p) !self % xsData % getTotalMatXS(p, p % matIdx()) !
    !print *, self % xsData % getMajorantXS(p)

    !if (p % time + distance / p % getSpeed() > maxT) p_p % fate = LEAK_FATE

    call self % geom % teleport(p_p % coords, distance)
    if (p_p % matIdx() == OUTSIDE_FILL) p_p % fate = LEAK_FATE

    if (p_p % matIdx() == UNDEF_MAT) call fatalError(Here, "Particle is in undefined material")

  end subroutine processParticlePrediction


end module transportOperatorDT_class
