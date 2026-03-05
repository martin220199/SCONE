module neutronMGkineticstd_class

  use numPrecision
  use endfConstants
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG

  ! Particle types
  use particle_class,                only : particle, particleState, printType, P_NEUTRON
  use particleDungeon_class,         only : particleDungeon

  ! Abstract interface
  use collisionProcessor_inter,      only : collisionProcessor, collisionData ,init_super => init

  ! Nuclear Data Interface
  use nuclearDataReg_mod,            only : ndReg_getNeutronMG => getNeutronMG
  use nuclearDatabase_inter,         only : nuclearDatabase
  use mgNeutronDatabase_inter,       only : mgNeutronDatabase
  use mgNeutronMaterial_inter,       only : mgNeutronMaterial, mgNeutronMaterial_CptrCast
  use reactionHandle_inter,          only : reactionHandle
  use multiScatterMG_class,          only : multiScatterMG, multiScatterMG_CptrCast
  use fissionMG_class,               only : fissionMG, fissionMG_TptrCast

  ! Cross section packages
  use neutronXsPackages_class,       only : neutronMacroXSs


  ! Nuclear Data
  !use nuclearData_inter,              only : nuclearData
  !use perMaterialNuclearDataMG_inter, only : perMaterialNuclearDataMG

  ! Cross-section packages to interface with nuclear data
  !use xsMacroSet_class,               only : xsMacroSet, xsMacroSet_ptr

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for MG neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  NONE
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type            neutronMGkineticstd;
  !!   }
  !!
  type, public, extends(collisionProcessor) :: neutronMGkineticstd
    private
    class(mgNeutronDatabase), pointer, public :: xsData => null()
    class(mgNeutronMaterial), pointer, public :: mat    => null()
    class(mgNeutronMaterial), pointer, public :: mat2    => null()

    ! Precursors
    logical(defBool) :: usePrecursors
  contains
    ! Initialisation procedure
    procedure :: init

    ! Implementation of customisable procedures
    procedure :: sampleCollision
    procedure :: implicit
    procedure :: elastic
    procedure :: inelastic
    procedure :: capture
    procedure :: fission
    procedure :: cutoffs
  end type neutronMGkineticstd

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronMGkineticstd), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    character(100), parameter :: Here = 'init (neutronMGkineticstd_class.f90)'

    ! Call superclass
    call init_super(self, dict)

    call dict % getOrDefault(self % usePrecursors, 'precursors', .false.)

  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGkineticstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(neutronMacroXSs)                :: macroXSs, macroXSs2
    real(defReal)                        :: r
    character(100),parameter :: Here =' sampleCollision (neutronMGkineticstd_class.f90)'

    ! Verify that particle is MG neutron
    if( .not. p % isMG .or. p % type /= P_NEUTRON) then
      call fatalError(Here, 'Supports only MG Neutron. Was given CE '//printType(p % type))
    end if

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getNeutronMG()
    if(.not.associated(self % xsData)) call fatalError(Here, "Failed to get active database for MG Neutron")

    !call self % xsData % activateT(p)

    ! Get and verify material pointer
    self % mat => mgNeutronMaterial_CptrCast( self % xsData % getMaterial( p % matIdx()))
    if(.not.associated(self % mat)) call fatalError(Here, "Failed to get MG Neutron Material")



    ! Select Main reaction channel
    call self % mat % getMacroXSs(macroXSs, p % G, p % pRNG)
    !check material index, set macroXSS
    if (p % matIdx() == 9999999) then !1
      self % mat2 => mgNeutronMaterial_CptrCast(self % xsData % getMaterial(2))
      call self % mat2 % getMacroXSsT(macroXSs2, p % G, p % pRNG, p)

      if (p % time <= 1e-3) then
        macroXSs % capture = macroXSs % capture * (1e-3 - p % time) / 1e-3 + macroXSs2 % capture * p % time / 1e-3
        macroXSs % inelasticScatter=macroXSs%inelasticScatter*(1e-3-p%time)/1e-3+macroXSs2%inelasticScatter*p%time/1e-3
      else
        macroXSs % capture = macroXSs2 % capture
        macroXSs % inelasticScatter = macroXSs2 % inelasticScatter
      end if
      macroXSs % total = macroXSs % elasticScatter + macroXSs % inelasticScatter + macroXSs % capture +  macroXSs % fission
      !print *, 'tot', macroXSs % total, p % G

    end if

    r = p % pRNG % get()

    collDat % MT = macroXSs % invert(r)

  end subroutine sampleCollision

  !!
  !! Preform implicit treatment
  !!
  subroutine implicit(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGkineticstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(neutronMacroXSs)                :: macroXSs, macroXSs2
    type(fissionMG),pointer              :: fission
    type(particleState)                  :: pTemp
    real(defReal),dimension(3)           :: r, dir
    integer(shortInt)                    :: G_out, n, i
    real(defReal)                        :: wgt, w0, rand1, mu, phi
    real(defReal)                        :: sig_tot, k_eff, sig_nufiss
    character(100),parameter :: Here = 'implicit (neutronMGkineticstd_class.f90)'

    if ( self % mat % isFissile()) then
      ! Obtain required data
      wgt   = p % w                ! Current weight
      w0    = p % preHistory % wgt ! Starting weight
      k_eff = p % k_eff            ! k_eff for normalisation
      rand1 = p % pRNG % get()     ! Random number to sample sites

      call self % mat % getMacroXSs(macroXSs, p % G, p % pRNG)

      sig_tot    = macroXSs % total
      sig_nuFiss = macroXSs % nuFission

      if (p % matIdx() == 9999999) then !1
        call self % mat2 % getMacroXSsT(macroXSs2, p % G, p % pRNG, p)
        if (p % time <= 1e-3) then
          macroXSs % capture = macroXSs % capture * (1e-3 - p % time) / 1e-3 + macroXSs2 % capture * p % time / 1e-3
          macroXSs % inelasticScatter=macroXSs%inelasticScatter*(1e-3-p%time)/1e-3+macroXSs2%inelasticScatter*p%time/1e-3
        else
          macroXSs % capture = macroXSs2 % capture
          macroXSs % inelasticScatter = macroXSs2 % inelasticScatter
        end if
        macroXSs % total = macroXSs % elasticScatter + macroXSs % inelasticScatter + macroXSs % capture +  macroXSs % fission
        sig_tot = macroXSs % total
      end if

      ! Sample number of fission sites generated
      !n = int(wgt * sig_nuFiss/(sig_tot*k_eff) + r1, shortInt)
      n = int(abs( (wgt * sig_nuFiss) / (w0 * sig_tot * k_eff)) + rand1, shortInt)

      ! Shortcut if no particles were samples
      if (n < 1) return

      ! Get Fission reaction object
      fission => fissionMG_TptrCast( self % xsData % getReaction(macroFission, collDat % matIdx))
      if (.not.associated(fission)) call fatalError(Here, 'Failed to getrive fissionMG reaction object')

      ! Store new sites in the next cycle dungeon
      wgt =  sign(w0, wgt)
      r   = p % rGlobal()

      do i=1,n
        call fission % sampleOut(mu, phi, G_out, p % G, p % pRNG)
        dir = rotateVector(p % dirGlobal(), mu, phi)

        ! Copy extra detail from parent particle (i.e. time, flags ect.)
        pTemp       = p

        ! Overwrite position, direction, energy group and weight
        pTemp % r   = r
        pTemp % dir = dir
        pTemp % G   = G_out
        pTemp % wgt = wgt
        pTemp % time = p % time

        call nextCycle % detain(pTemp)
      end do

    end if

  end subroutine implicit

  !!
  !! Elastic Scattering
  !!
  subroutine elastic(self, p , collDat, thisCycle, nextCycle)
    class(neutronMGkineticstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    ! Do nothing. Should not be called

  end subroutine elastic

  !!
  !! Preform scattering
  !!
  subroutine inelastic(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGkineticstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    class(multiScatterMG),pointer        :: scatter, scatter2
    integer(shortInt)                    :: G_out   ! Post-collision energy group
    real(defReal)                        :: phi     ! Azimuthal scatter angle
    real(defReal)                        :: w_mul   ! Weight multiplier
    real(defReal)                        :: rem, interScat, interP0
    character(100),parameter :: Here = "inelastic (neutronMGkineticstd_class.f90)"

    ! Assign MT number
    collDat % MT = macroIEscatter

    if (collDat % matIdx /= 999999) then !1
      ! Get Scatter object
      scatter => multiScatterMG_CptrCast( self % xsData % getReaction(macroIEscatter, collDat % matIdx))
      if(.not.associated(scatter)) call fatalError(Here, "Failed to get scattering reaction object for MG neutron")

      ! Sample Mu and G_out
      call scatter % sampleOut(collDat % muL, phi, G_out, p % G, p % pRNG)

      ! Read scattering multiplicity
      w_mul = scatter % production(p % G, G_out)

      ! Update neutron state
      p % G = G_out
      p % w = p % w * w_mul
      call p % rotate(collDat % muL, phi)
  
    else 

      ! Get Scatter object
      scatter => multiScatterMG_CptrCast( self % xsData % getReaction(macroIEscatter, collDat % matIdx))
      if(.not.associated(scatter)) call fatalError(Here, "Failed to get scattering reaction object for MG neutron")

      scatter2 => multiScatterMG_CptrCast( self % xsData % getReaction(macroIEscatter, 2))
      if(.not.associated(scatter2)) call fatalError(Here, "Failed to get scattering reaction object for MG neutron")


      if (p % G < 0 .or. p % G > size(scatter % scatterXSs)) then
        call fatalError(Here, 'Invalid incident group number: '//numToChar(p % G))
      end if
      if (p % G < 0 .or. p % G > size(scatter2 % scatterXSs)) then
        call fatalError(Here, 'Invalid incident group number: '//numToChar(p % G))
      end if


      if (p % time <= 1e-3) then
        ! Perform sampling
        interScat = scatter % scatterXSs(p % G) * (1e-3 - p % time) / 1e-3 + scatter2 % scatterXSs(p % G) * p % time / 1e-3
        rem = p % pRNG % get() * interScat

        do G_out = 1,size(scatter2 % scatterXSs)
          interP0 = scatter % P0(G_out, p % G) * (1e-3 - p % time) / 1e-3 + scatter2 % P0(G_out, p % G) * p % time / 1e-3
          rem = rem - interP0
          if(rem < ZERO) exit
        end do
      else 
        ! Perform sampling
        rem = p % pRNG % get() * scatter2 % scatterXSs(p % G)
        do G_out = 1,size(scatter2 % scatterXSs)
          rem = rem - scatter2 % P0(G_out, p % G)
          if(rem < ZERO) exit
        end do
      end if

      collDat % muL  = TWO * p % pRNG % get() - ONE
      phi = TWO_PI * p % pRNG % get()

      ! Read scattering multiplicity
      w_mul = scatter % production(p % G, G_out)

      ! Update neutron state
      p % G = G_out
      p % w = p % w * w_mul
      call p % rotate(collDat % muL, phi)
    end if


  end subroutine inelastic

  !!
  !! Preform capture
  !!
  subroutine capture(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGkineticstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead = .true.

  end subroutine capture

  !!
  !! Preform fission
  !!
  subroutine fission(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGkineticstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    p % isDead = .true.

  end subroutine fission

  !!
  !! Applay cutoffs or post-collision implicit treatment
  !!
  subroutine cutoffs(self, p, collDat, thisCycle, nextCycle)
    class(neutronMGkineticstd), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle

    ! Do nothing

  end subroutine cutoffs

end module neutronMGkineticstd_class
