module neutronCEkineticimp_class

  use numPrecision
  use endfConstants
  use universalVariables,            only : precursorGroups
  use genericProcedures,             only : fatalError, rotateVector, numToChar
  use dictionary_class,              only : dictionary
  use RNG_class,                     only : RNG

  ! Particle types
  use particle_class,                only : particle, particleState, printType, P_NEUTRON, P_PRECURSOR
  use particleDungeon_class,         only : particleDungeon

  ! Abstarct interface
  use collisionProcessor_inter,      only : collisionProcessor, collisionData ,init_super => init

  ! Nuclear Data Interfaces
  use nuclearDataReg_mod,            only : ndReg_getNeutronCE => getNeutronCE
  use nuclearDatabase_inter,         only : nuclearDatabase
  use ceNeutronDatabase_inter,       only : ceNeutronDatabase
  use ceNeutronMaterial_class,       only : ceNeutronMaterial, ceNeutronMaterial_CptrCast
  use ceNeutronNuclide_inter,        only : ceNeutronNuclide, ceNeutronNuclide_CptrCast

  ! Nuclear reactions
  use reactionHandle_inter,          only : reactionHandle
  use uncorrelatedReactionCE_inter,  only : uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use neutronScatter_class,          only : neutronScatter, neutronScatter_TptrCast
  use fissionCE_class,               only : fissionCE, fissionCE_TptrCast

  ! Cross-Section Packages
  use neutronXsPackages_class,       only : neutronMicroXSs

  ! Scattering procedures
  use scatteringKernels_func, only : asymptoticScatter, targetVelocity_constXS, &
                                     asymptoticInelasticScatter

  ! Tally interfaces
  use tallyAdmin_class,       only : tallyAdmin

  implicit none
  private

  !!
  !! Standard (default) scalar collision processor for CE neutrons
  !!   -> Preforms implicit fission site generation
  !!   -> Preforms analog capture
  !!   -> Treats fission as capture (only implicit generation of 2nd-ary neutrons)
  !!   -> Does not create secondary non-neutron projectiles
  !!
  !! Settings:
  !!  minE       -> minimum energy cut-off [MeV] (default = 1.0E-11)
  !!  maxE       -> maximum energy. Higher energies are set to maximum (not re-rolled) [MeV]
  !!             (default = 20.0)
  !!  thresh_E   -> Energy threshold for explicit treatment of target nuclide movement [-].
  !!              Target movment is sampled if neutron energy E < kT * thresh_E where
  !!              kT is target material temperature in [MeV]. (default = 400.0)
  !!  thresh_A   -> Mass threshold for explicit tratment of target nuclide movement [Mn].
  !!              Target movment is sampled if target mass A < thresh_A. (default = 1.0)
  !!  minWgt     -> minimum particle weight for rouletting (optional)
  !!  maxWgt     -> maximum particle weight for splitting (optional)
  !!  avgWgt     -> weight of a particle on surviving rouletting (optional)
  !!  splitting  -> splits particles above certain weight (off by default)
  !!  roulette   -> roulettes particles below certain weight (off by default)
  !!  branchless -> branchless collision
  !!
  !! Sample dictionary input:
  !!   collProcName {
  !!   type             neutronCEkineticimp;
  !!   #minEnergy       <real>;#
  !!   #maxEnergy       <real>;#
  !!   #energyThreshold <real>;#
  !!   #massThreshold   <real>;#
  !!   #splitting       <logical>;#
  !!   #roulette        <logical>;#
  !!   #minWgt          <real>;#
  !!   #maxWgt          <real>;#
  !!   #avgWgt          <real>;#
  !!   #precursors      <logical>;#
  !!   #branchless      <logical>;#
  !!   }
  !!
  type, public, extends(collisionProcessor) :: neutronCEkineticimp
    private
    !! Nuclear Data block pointer -> public so it can be used by subclasses (protected member)
    class(ceNeutronDatabase), pointer, public :: xsData => null()
    class(ceNeutronMaterial), pointer, public :: mat    => null()
    class(ceNeutronNuclide),  pointer, public :: nuc    => null()

    !! Settings - private
    real(defReal)    :: minE
    real(defReal)    :: maxE
    real(defReal)    :: thresh_E
    real(defReal)    :: thresh_A
    real(defReal)    :: minWgt
    real(defReal)    :: maxWgt
    real(defReal)    :: avWgt
    logical(defBool) :: splitting
    logical(defBool) :: roulette
    logical(defBool) :: branchless

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

    ! Local procedures
    procedure,private :: scatterFromFixed
    procedure,private :: scatterFromMoving
    procedure,private :: scatterInLAB

    ! Variance reduction procedures
    procedure, private :: split
    procedure, private :: russianRoulette
  end type neutronCEkineticimp

contains

  !!
  !! Initialise from dictionary
  !!
  subroutine init(self, dict)
    class(neutronCEkineticimp), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    character(100), parameter :: Here = 'init (neutronCEkineticimp_class.f90)'

    ! Call superclass
    call init_super(self, dict)

    ! Read settings for neutronCEkineticimp
    ! Maximum and minimum energy
    call dict % getOrDefault(self % minE,'minEnergy',1.0E-11_defReal)
    call dict % getOrDefault(self % maxE,'maxEnergy',20.0_defReal)

    ! Thermal scattering kernel thresholds
    call dict % getOrDefault(self % thresh_E, 'energyThreshold', 400.0_defReal)
    call dict % getOrDefault(self % thresh_A, 'massThreshold', 1.0_defReal)

    ! Verify settings
    if( self % minE < ZERO ) call fatalError(Here,'-ve minEnergy')
    if( self % maxE < ZERO ) call fatalError(Here,'-ve maxEnergy')
    if( self % minE >= self % maxE) call fatalError(Here,'minEnergy >= maxEnergy')
    if( self % thresh_E < 0) call fatalError(Here,' -ve energyThreshold')
    if( self % thresh_A < 0) call fatalError(Here,' -ve massThreshold')

    ! Obtain precursor settings
    call dict % getOrDefault(self % usePrecursors, 'precursors', .false.)

    ! Obtain settings for variance reduction
    call dict % getOrDefault(self % splitting,'split', .false.)
    call dict % getOrDefault(self % roulette,'roulette', .false.)
    call dict % getOrDefault(self % minWgt,'minWgt', 0.4_defReal)
    call dict % getOrDefault(self % maxWgt,'maxWgt', 2.0_defReal)
    call dict % getOrDefault(self % avWgt,'avWgt', 1.0_defReal)
    call dict % getOrDefault(self % branchless,'branchless', .false.)

    ! Ensure compatability
    if (self % splitting .eqv. .true.) then
      if (self % maxWgt < 2 * self % minWgt) call fatalError(Here,&
              'Upper weight bound must be at least twice the lower weight bound')
    end if

    if (self % branchless .eqv. .true.) then
      if (self % usePrecursors .eqv. .false.) then
        self % splitting = .true.
        self % roulette = .true.
      else
        self % splitting = .false.
      end if
    end if
  end subroutine init

  !!
  !! Samples collision without any implicit treatment
  !!
  subroutine sampleCollision(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEkineticimp), intent(inout)   :: self
    class(particle), intent(inout)       :: p
    type(tallyAdmin), intent(inout)      :: tally
    type(collisionData), intent(inout)   :: collDat
    class(particleDungeon),intent(inout) :: thisCycle
    class(particleDungeon),intent(inout) :: nextCycle
    type(neutronMicroXSs)                :: microXSs
    real(defReal)                        :: r
    character(100),parameter :: Here = 'sampleCollision (neutronCEkineticimp_class.f90)'

    ! Verify that particle is CE neutron
    if(p % isMG .or. p % type /= P_NEUTRON) then
      call fatalError(Here, 'Supports only CE Neutron. Was given MG '//printType(p % type))
    end if

    ! Verify and load nuclear data pointer
    self % xsData => ndReg_getNeutronCE()
    if(.not.associated(self % xsData)) call fatalError(Here, 'There is no active Neutron CE data!')

    ! Verify and load material pointer
    self % mat => ceNeutronMaterial_CptrCast( self % xsData % getMaterial( p % matIdx()))
    if(.not.associated(self % mat)) call fatalError(Here, 'Material is not ceNeutronMaterial')

    ! Select collision nuclide
    collDat % nucIdx = self % mat % sampleNuclide(p % E, p % pRNG)

    self % nuc => ceNeutronNuclide_CptrCast(self % xsData % getNuclide(collDat % nucIdx))
    if(.not.associated(self % mat)) call fatalError(Here, 'Failed to retive CE Neutron Nuclide')

    ! Select Main reaction channel
    call self % nuc % getMicroXSs(microXss, p % E, p % pRNG)
    r = p % pRNG % get()

  if (self % branchless .eqv. .true.) then
    collDat % MT = microXss % invertBranchless(r)
  else
    collDat % MT = microXss % invert(r)
  end if

  end subroutine sampleCollision

  !!
  !! Perform implicit treatment
  !!
  subroutine implicit(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEkineticimp), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    type(collisionData), intent(inout)        :: collDat
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle
    type(fissionCE), pointer                  :: fiss
    type(neutronMicroXSs)                     :: microXSs
    type(particleState)                       :: pTemp
    real(defReal),dimension(3)                :: r, dir
    integer(shortInt)                         :: n, i
    real(defReal)                             :: wgt, w0, E_out, mu, phi, lambda, wgtFactor
    real(defReal)                             :: sig_nuPromptfiss, sig_nuDelayedfiss, sig_tot, k_eff
    character(100),parameter                  :: Here = 'implicit (neutronCEkineticimp_class.f90)'

    if (self % branchless .eqv. .true.) then
      ! Obtain micro cross-sections
      call self % nuc % getMicroXSs(microXSs, p % E, p % pRNG)

      ! Compute weight multiplier when applying branchless on isotope
      wgtFactor = (microXSs % nuFission + microXSs % elasticScatter + microXSs % inelasticScatter) &
                  / microXSs % total

      ! Modify weight at each collision
      p % w = p % w * wgtFactor

    ! Generate fission sites if nuclide is fissile
    else if ( self % nuc % isFissile()) then
      ! Obtain required data
      wgt   = p % w                ! Current weight
      w0    = p % preHistory % wgt ! Starting weight
      k_eff = p % k_eff            ! k_eff for normalisation

      call self % nuc % getMicroXSs(microXSs, p % E, p % pRNG)
      sig_tot    = microXSs % total

      ! Get fission Reaction
      fiss => fissionCE_TptrCast(self % xsData % getReaction(N_FISSION, collDat % nucIdx))
      if(.not.associated(fiss)) call fatalError(Here, "Failed to get fissionCE")
      sig_nuPromptfiss = fiss % releasePrompt(p % E) * microXSs % fission


      ! Sample number of prompt fission sites generated
      ! Support -ve weight particles
      n = int(abs( (wgt * sig_nuPromptfiss) / (w0 * sig_tot * k_eff)) + p % pRNG % get(), shortInt)

      if (n >= 1) then
        ! Store new sites in the next cycle dungeon
        wgt =  sign(w0, wgt)
        r   = p % rGlobal()

        do i = 1, n
          call fiss % samplePrompt(mu, phi, E_out, p % E, p % pRNG)
          dir = rotateVector(p % dirGlobal(), mu, phi)

          if (E_out > self % maxE) E_out = self % maxE

          ! Copy extra detail from parent particle (i.e. time, flags ect.)
          pTemp = p

          ! Overwrite position, direction, energy and weight
          pTemp % r   = r
          pTemp % dir = dir
          pTemp % E   = E_out
          pTemp % wgt = wgt
          pTemp % type      = P_NEUTRON
          pTemp % timeBirth = p % time

          call nextCycle % detain(pTemp)

          ! Report birth of new particle
          call tally % reportSpawn(N_FISSION, p, pTemp)

        end do
      end if

      if (self % usePrecursors) then

        ! Handle delayed neutrons using Forced Decay
        sig_nuDelayedfiss = fiss % releaseDelayed(p % E) * microXSs % fission
        n = int(abs( (wgt * sig_nuDelayedfiss) / (w0 * sig_tot * k_eff)) + p % pRNG % get(), shortInt)

        if (n >= 1) then

          wgt =  sign(w0, wgt)
          r   = p % rGlobal()

          do i = 1, n
            call fiss % sampleDelayed(mu, phi, E_out, p % E, p % pRNG, lambda)

            dir = rotateVector(p % dirGlobal(), mu, phi)

            if (E_out > self % maxE) E_out = self % maxE

            ! Copy extra detail from parent particle (i.e. time, flags ect.)
            pTemp       = p

            ! Overwrite particle attributes
            pTemp % r   = r
            pTemp % dir = dir
            pTemp % E   = E_out
            pTemp % wgt = wgt
            pTemp % type      = P_PRECURSOR
            pTemp % lambda    = lambda
            pTemp % timeBirth = p % time

            call thisCycle % detain(pTemp)

            ! Report birth of new particle
            call tally % reportSpawn(N_FISSION, p, pTemp)

          end do
        end if
      end if
    end if

  end subroutine implicit

  !!
  !! Process capture reaction
  !!
  subroutine capture(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEkineticimp), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    type(collisionData), intent(inout)        :: collDat
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle

    p % isDead = .true.

  end subroutine capture

  !!
  !! Process fission reaction
  !!
  subroutine fission(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEkineticimp), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    type(collisionData), intent(inout)        :: collDat
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle
    type(fissionCE), pointer                  :: fiss
    type(neutronMicroXSs)                     :: microXSs
    type(particleState)                       :: pTemp
    real(defReal),dimension(3)                :: r, dir
    real(defReal)                             :: wgt, w0, E_out, mu, phi, lambda, rand, probabilityOfPrompt, k_eff
    character(100),parameter                  :: Here = 'fission (neutronCEkineticimp_class.f90)'

    if ((self % branchless .eqv. .true.) .and. (self % nuc % isFissile() .eqv. .true.)) then

      wgt   = p % w                ! Current weight
      w0    = p % preHistory % wgt ! Starting weight
      k_eff = p % k_eff            ! k_eff for normalisation

      ! Get fission Reaction
      fiss => fissionCE_TptrCast(self % xsData % getReaction(N_FISSION, collDat % nucIdx))
      if(.not.associated(fiss)) call fatalError(Here, "Failed to get fissionCE")
      call self % nuc % getMicroXSs(microXSs, p % E, p % pRNG)

      probabilityOfPrompt = fiss % releasePrompt(p % E) / fiss % release(p % E)
      rand = p % pRNG % get()

      if (rand <= probabilityOfPrompt) then

        r   = p % rGlobal()
        call fiss % samplePrompt(mu, phi, E_out, p % E, p % pRNG)
        dir = rotateVector(p % dirGlobal(), mu, phi)

        if (E_out > self % maxE) E_out = self % maxE

        ! Copy extra detail from parent particle (i.e. time, flags ect.)
        pTemp       = p

        ! Overwrite position, direction, energy and weight
        pTemp % r   = r
        pTemp % dir = dir
        pTemp % E   = E_out
        pTemp % timeBirth = p % time

        call nextCycle % detain(pTemp)

        ! Report birth of new particle
        call tally % reportSpawn(N_FISSION, p, pTemp)

      else if (self % usePrecursors) then

        r   = p % rGlobal()
        call fiss % sampleDelayed(mu, phi, E_out, p % E, p % pRNG, lambda)
        if (E_out > self % maxE) E_out = self % maxE

        dir = rotateVector(p % dirGlobal(), mu, phi)

        ! Copy extra detail from parent particle (i.e. time, flags ect.)
        pTemp       = p

        ! Overwrite particle attributes
        pTemp % r   = r
        pTemp % dir = dir
        pTemp % E   = E_out
        pTemp % type      = P_PRECURSOR
        pTemp % lambda    = lambda
        pTemp % timeBirth = p % time

        call thisCycle % detain(pTemp)

        ! Report birth of new particle
        call tally % reportSpawn(N_FISSION, p, pTemp)

      end if
    end if

    p % isDead = .true.

  end subroutine fission

  !!
  !! Process elastic scattering
  !!
  !! All CE elastic scattering happens in the CM frame
  !!
  subroutine elastic(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEkineticimp), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    type(collisionData), intent(inout)        :: collDat
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle
    class(uncorrelatedReactionCE), pointer    :: reac
    logical(defBool)                          :: isFixed
    character(100),parameter :: Here = 'elastic (neutronCEkineticimp_class.f90)'

    ! Get reaction
    reac => uncorrelatedReactionCE_CptrCast( self % xsData % getReaction(collDat % MT, collDat % nucIdx))
    if(.not.associated(reac)) call fatalError(Here,'Failed to get elastic neutron scatter')

    ! Scatter particle
    collDat % A =  self % nuc % getMass()
    collDat % kT = self % nuc % getkT()

    isFixed = (p % E > collDat % kT * self % thresh_E) .and. (collDat % A > self % thresh_A)

    ! Apply criterion for Free-Gas vs Fixed Target scattering
    if (.not. reac % inCMFrame()) then
      call self % scatterInLAB(p, collDat, reac)
    elseif (isFixed) then
      call self % scatterFromFixed(p, collDat, reac)
    else
      call self % scatterFromMoving(p, collDat, reac)
    end if

  end subroutine elastic

  !!
  !! Process inelastic scattering
  !!
  subroutine inelastic(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEkineticimp), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    type(collisionData), intent(inout)        :: collDat
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle
    class(uncorrelatedReactionCE), pointer    :: reac
    character(100),parameter  :: Here =' inelastic (neutronCEkineticimp_class.f90)'

    ! Invert inelastic scattering and Get reaction
    collDat % MT = self % nuc % invertInelastic(p % E, p % pRNG)
    reac => uncorrelatedReactionCE_CptrCast( self % xsData % getReaction(collDat % MT, collDat % nucIdx))
    if(.not.associated(reac)) call fatalError(Here, "Failed to get scattering reaction")

    ! Scatter particle
    if (reac % inCMFrame()) then
      collDat % A =  self % nuc % getMass()
      call self % scatterFromFixed(p, collDat, reac)
    else
      call self % scatterInLAB(p, collDat, reac)
    end if

    ! Apply weigth change
    p % w = p % w * reac % release(p % E)

  end subroutine inelastic

  !!
  !! Apply cutoffs
  !!
  subroutine cutoffs(self, p, tally, collDat, thisCycle, nextCycle)
    class(neutronCEkineticimp), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    type(collisionData), intent(inout)        :: collDat
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle

    if (p % E < self % minE ) p % isDead = .true.
    if (p % isDead .eqv. .true.) return

    ! Splitting with fixed threshold
    if ((p % w > self % maxWgt) .and. (self % splitting .eqv. .true.)) then
      call self % split(p, thisCycle, self % maxWgt)
    ! Roulette with fixed threshold and survival weight
    elseif ((p % w < self % minWgt) .and. (self % roulette .eqv. .true.)) then
      call self % russianRoulette(p, self % avWgt)
    end if

  end subroutine cutoffs

  !!
  !! Perform Russian roulette on a particle
  !!
  subroutine russianRoulette(self, p, avWgt)
    class(neutronCEkineticimp), intent(inout) :: self
    class(particle), intent(inout)            :: p
    real(defReal), intent(in)                 :: avWgt

    if (p % pRNG % get() < (ONE - p % w/avWgt)) then
      p % isDead = .true.
    else
      p % w = avWgt
    end if

  end subroutine russianRoulette

  !!
  !! Split particle which has too large a weight
  !!
  subroutine split(self, p, nextCycle, maxWgt)
    class(neutronCEkineticimp), intent(inout) :: self
    class(particle), intent(inout)            :: p
    class(particleDungeon), intent(inout)     :: nextCycle
    real(defReal), intent(in)                 :: maxWgt
    integer(shortInt)                         :: mult, i

    ! This value must be at least 2
    mult = ceiling(p % w/maxWgt)

    ! Decrease weight
    p % w = p % w/mult

    ! Add split particle's to the dungeon
    do i = 1,mult-1
      call nextCycle % detain(p)
    end do

  end subroutine split

  !!
  !! Subroutine to perform scattering in LAB frame
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterInLAB(self, p, collDat, reac)
    class(neutronCEkineticimp), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(collisionData), intent(inout)        :: collDat
    class(uncorrelatedReactionCE), intent(in) :: reac
    real(defReal)                             :: phi    ! Azimuthal scatter angle
    real(defReal)                             :: E_out, mu

    ! Sample scattering angles and post-collision energy
    call reac % sampleOut(mu, phi, E_out, p % E, p % pRNG)

    ! Update neutron state
    p % E = E_out
    call p % rotate(mu, phi)
    collDat % muL = mu

  end subroutine scatterInLAB

  !!
  !! Subroutine to perform scattering from stationary target.
  !! Returns mu -> cos of deflection angle in LAB frame
  !!
  subroutine scatterFromFixed(self, p, collDat, reac)
    class(neutronCEkineticimp), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(collisionData), intent(inout)        :: collDat
    class(uncorrelatedReactionCE), intent(in) :: reac
    real(defReal)                             :: phi
    real(defReal)                             :: E_out
    real(defReal)                             :: E_outCM, mu
    integer(shortInt)                         :: MT

    ! Read data
    MT     = collDat % MT

    ! Sample mu , phi and outgoing energy
    call reac % sampleOut(mu, phi, E_outCM, p % E, p % pRNG)

    ! Save incident energy
    E_out = p % E

    if( MT == N_N_elastic) then
      call asymptoticScatter(E_out, mu, collDat % A)
    else
      call asymptoticInelasticScatter(E_out, mu, E_outCM, collDat % A)
    end if

    ! Update particle state
    call p % rotate(mu, phi)
    p % E = E_out
    collDat % muL = mu

  end subroutine scatterFromFixed

  !!
  !! Subroutine to perform scattering from moving target
  !! Supports only elastic collisions
  !!
  subroutine scatterFromMoving(self, p, collDat, reac)
    class(neutronCEkineticimp), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(collisionData),intent(inout)         :: collDat
    class(uncorrelatedReactionCE), intent(in) :: reac
    real(defReal)                             :: A, kT, mu
    real(defReal),dimension(3)                :: V_n           ! Neutron velocity (vector)
    real(defReal)                             :: U_n           ! Neutron speed (scalar)
    real(defReal),dimension(3)                :: dir_pre       ! Pre-collision direction
    real(defReal),dimension(3)                :: dir_post      ! Post-collicion direction
    real(defReal),dimension(3)                :: V_t, V_cm     ! Target and CM velocity
    real(defReal)                             :: phi, dummy

    ! Read data
    A      = collDat % A
    kT     = collDat % kT

    ! Get neutron direction and velocity
    dir_pre = p % dirGlobal()
    V_n     = dir_pre * sqrt(p % E)

    ! Sample velocity of target
    V_t = targetVelocity_constXS(p % E, dir_pre, A, kT, p % pRNG)

    ! Calculate Centre-of-Mass velocity
    V_cm = (V_n + V_t *A)/(A+1)

    ! Move Neutron velocity to CM frame, store speed and calculate new normalised direction
    V_n = V_n - V_cm
    U_n = norm2(V_n)
    V_n = V_n / U_n

    ! Sample mu and phi in CM frame
    call reac % sampleOut(mu, phi, dummy, p % E, p % pRNG)

    ! Obtain post collision speed
    V_n = rotateVector(V_n, mu, phi) * U_n

    ! Return to LAB frame
    V_n = V_n + V_cm

    ! Calculate new neutron speed and direction
    U_n = norm2(V_n)
    dir_post = V_n / U_n

    ! Update particle state and calculate mu in LAB frame
    p % E = U_n * U_n
    call p % point(dir_post)
    collDat % muL = dot_product(dir_pre, dir_post)

  end subroutine scatterFromMoving


end module neutronCEkineticimp_class
