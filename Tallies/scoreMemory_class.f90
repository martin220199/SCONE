module scoreMemory_class

  use numPrecision
  use universalVariables, only : array_pad
  use genericProcedures,  only : fatalError, numToChar
  use openmp_func,        only : ompGetMaxThreads, ompGetThreadNum

  implicit none
  private

  !! Parameters for indexes of per cycle SCORE, Cumulative Sum and Cumulative Sum of squares
  integer(shortInt), parameter :: CSUM  = 1, &
                                  CSUM2 = 2

  !! Size of the 2nd Dimension of bins
  integer(shortInt), parameter :: DIM2 = 2


  !!
  !! scoreMemory is a class that stores space for scores from tallies.
  !! It is separate from tallyClerks and individual responses to allow:
  !!   -> Easy writing and (later) reading from file for archivisation of results
  !!   -> Easy possibility of extention to tally higher moments of result
  !!   -> Possibility of extension to tally covariance of selected tally bins
  !!   -> Easy copying and recombination of results for OpenMP shared memory parallelism
  !!   -> Easy, output format-independent way to perform regression tests
  !!   -> Easy handling of different batch sizes
  !!
  !! For every bin index there are two positions: CSUM, CSUM2. All are initialised to 0.
  !! For scoring, an array is created with dimension (Nbins,nThreads) to mitigate false sharing.
  !! On accumulation, this array adds to the normal bin array.
  !!
  !! Interface:
  !!     init(N,idBS): Initialise with integer size N and integer id. Optional integer Batch Size.
  !!
  !!     kill(): Elemental. Return to uninitialised state.
  !!
  !!     score(score,idx): Score in the bin under idx. FatalError if idx is outside bounds. Score
  !!         is defReal, shortInt or longInt
  !!
  !!     accumulate(score,idx): Accumulate result in cumulative sums in bin under idx. FatalError
  !!         if idx is outside bounds. Score is defReal, shortInt or longInt.
  !!
  !!     getResult(mean, STD, idx, samples): Retrieve mean value and standard deviation of the
  !!         estimate under idx. Use optional samples to specify number of estimates used to
  !!         evaluate mean and STD from default, which is number of batches in score memory.
  !!         STD is optional.
  !!
  !!     getScore(idx): Return current value of score from bin under idx. FatalError if idx is
  !!         outside bounds.
  !!
  !!     closeBin(normFactor,idx): Multiplies score under bin by normFactor and accumulates it in
  !!         cumulative sums. Then sets the bin to zero.
  !!
  !!     closeCycle(normFactor): Multiplies all scores by normFactor and accumulates them in
  !!         cumulative sums. Sets all scors to zero.
  !!
  !!     lastCycle(): Return true if the next call to closeCycle will close a batch.
  !!
  !!     getBatchSize(): Returns number of cycles that constitute a single batch.
  !!
  !!     setNumBatchesPerTimeStep(batchN): Sets batchN member of scoreMemory for kinetic treatment
  !!
  !! Example use case:
  !!
  !!  do batches=1,20
  !!    do hist=1,10
  !!      call scoreMem % score(hist,1)        ! Score hist (1,10) in bin 1
  !!      call scoreMem % accumulate(hist,2)   ! Accumulate hist in CSUMs of bin 2
  !!    end do
  !!    call scoreMem % closeCycle(ONE)        ! Close batch without normalisation (factor = ONE)
  !!  end do
  !!
  !!  call scoreMem % getResult(mean,STD,1) ! Get result from bin 1 with STD
  !!  call scoreMem % getResult(mean,2,200) ! Get mean from bin 2 assuming 200 samples
  !!
  !! NOTE:  Following indexing is used in bins class member
  !!        bins(binIndex,binType) binType is CSUM/CSUM2
  !! NOTE2: If batch size is not a denominator of cycles scored results accumulated
  !!        in extra cycles are discarded in current implementation
  !!
  type, public :: scoreMemory
      !private
      real(defReal),dimension(:,:,:),allocatable :: bins, bins_b          !! Space for storing cumul data (2nd dim size is always 2!)
      real(defReal),dimension(:,:,:),allocatable :: parallelBins, parallelBins_b  !! Space for scoring for different threads
      integer(longInt)                         :: N = 0         !! Size of memory (number of bins)
      integer(shortInt)                        :: nThreads = 0  !! Number of threads used for parallelBins
      integer(shortInt)                        :: id            !! Id of the tally
      integer(shortInt)                        :: batchN = 0    !! Number of Batches
      integer(shortInt)                        :: cycles = 0    !! Cycles counter
      integer(shortInt)                        :: batchSize = 1 !! Batch interval size (in cycles)
      integer(shortInt)                        :: basis
      integer(shortInt), dimension(:), allocatable :: maxFetOrder
      real(defReal), dimension(:), allocatable :: piecewise
      real(defReal)                              :: a, b
      real(defReal), dimension(:), allocatable   :: FET_evalPoints, minT, maxT, deltaT  
      procedure(scoreFET_signature), pointer     :: scoreFET => null()
      procedure(getFETResult_signature), pointer :: getFETResult => null()
      logical(defBool)                           :: useFET = .false.
  contains
    ! Interface procedures
    procedure :: init
    procedure :: kill
    generic   :: score      => score_defReal, score_shortInt, score_longInt
    generic   :: accumulate => accumulate_defReal, accumulate_shortInt, accumulate_longInt
    generic   :: getResult  => getResult_withSTD, getResult_withoutSTD
    procedure :: getScore
    procedure :: closeCycle
    procedure :: closeBin
    procedure :: lastCycle
    procedure :: getBatchSize
    procedure :: setNumBatchesPerTimeStep
    procedure :: mapPiecewise
    procedure :: getResult_FET

    ! Public procedures required for Fourier FET
    procedure :: getResult_Fourier_b
    procedure :: getFETResult_Fourier

    ! Private procedures
    procedure, private :: score_defReal
    procedure, private :: score_shortInt
    procedure, private :: score_longInt
    procedure, private :: accumulate_defReal
    procedure, private :: accumulate_shortInt
    procedure, private :: accumulate_longInt
    procedure, private :: getResult_withSTD
    procedure, private :: getResult_withoutSTD

    ! Private Legendre Procedures
    procedure, private :: scoreFET_Legendre
    procedure, private :: transDomain_Legendre
    procedure, private :: weight_Legendre
    procedure, private :: orthonormalisationTransformed_Legendre
    procedure, private :: orthonormalisationStandard_Legendre
    procedure, private :: getFETResult_Legendre

    ! Private Chebyshev1 Procedures
    procedure, private :: scoreFET_Chebyshev1
    procedure, private :: transDomain_Chebyshev1
    procedure, private :: weight_Chebyshev1
    procedure, private :: orthonormalisationTransformed_Chebyshev1
    procedure, private :: orthonormalisationStandard_Chebyshev1
    procedure, private :: getFETResult_Chebyshev1

    ! Private Chebyshev2 Procedures
    procedure, private :: scoreFET_Chebyshev2
    procedure, private :: transDomain_Chebyshev2
    procedure, private :: weight_Chebyshev2
    procedure, private :: orthonormalisationTransformed_Chebyshev2
    procedure, private :: orthonormalisationStandard_Chebyshev2
    procedure, private :: getFETResult_Chebyshev2


    ! Private Fourier Procedures
    procedure, private :: scoreFET_Fourier
    procedure, private :: transDomain_Fourier
    procedure, private :: weight_Fourier
    procedure, private :: orthonormalisationTransformed_Fourier
    procedure, private :: orthonormalisationStandard_Fourier
    procedure, private :: orthonormalisationTransformed_b_Fourier
    procedure, private :: orthonormalisationStandard_b_Fourier

    ! Private Jacobi Procedures
    procedure, private :: scoreFET_Jacobi
    procedure, private :: transDomain_Jacobi
    procedure, private :: weight_Jacobi
    procedure, private :: orthonormalisationTransformed_Jacobi
    procedure, private :: orthonormalisationStandard_Jacobi
    procedure, private :: getFETResult_Jacobi

  end type scoreMemory

contains

  !!
  !! Allocate space for the bins given number of bins N
  !! Optionaly change batchSize from 1 to any +ve number
  !!
  subroutine init(self, N, id, batchSize, maxFetOrder, minT, maxT, FET_evalPoints, basisFlag, a, b, piecewise)
    class(scoreMemory),intent(inout)       :: self
    integer(longInt),intent(in)            :: N
    integer(shortInt),intent(in)           :: id
    integer(shortInt),optional,intent(in)  :: batchSize
    integer(shortInt), dimension(:), optional, intent(in) :: maxFetOrder
    real(defReal), dimension(:), optional, intent(in)    :: minT, maxT
    integer(shortInt), optional,intent(in) :: FET_evalPoints
    integer(shortInt), optional,intent(in) :: basisFlag
    real(defReal), optional,intent(in)     :: a, b
    real(defReal), dimension(:), optional, intent(in) :: piecewise
    integer(shortInt)                      :: i
    real(defReal)                          :: offset, deltaEval
    character(100),parameter :: Here = 'Init (scoreMemory_class.f90)'

    self % nThreads = ompGetMaxThreads()
    if (present(maxFetOrder)) then
      self % useFET = .true.
      self % basis = basisFlag
      select case(basisFlag)
        case(0)
          self % scoreFET => scoreFET_Legendre
          self % getFETResult => getFETResult_Legendre
        case(1)
          self % scoreFET => scoreFET_Chebyshev1
          self % getFETResult => getFETResult_Chebyshev1
        case(2)
          self % scoreFET => scoreFET_Chebyshev2
          self % getFETResult => getFETResult_Chebyshev2
        case(5)
          self % scoreFET => scoreFET_Fourier
          allocate( self % parallelBins_b(maxval(maxFetOrder) + 1, self % nThreads, size(maxFetOrder)))
          self % parallelBins_b = ZERO
          allocate( self % bins_b(maxval(maxFetOrder) + 1, DIM2, size(maxFetOrder)))
          self % bins_b = ZERO
        case(6)
          self % scoreFET => scoreFET_Jacobi
          self % getFETResult => getFETResult_Jacobi
          self % a = a
          self % b = b
        case default
          call fatalError(Here, 'Need to define the basis function')
      end select

      self % maxFetOrder = maxFetOrder
      self % piecewise = piecewise
      self % deltaT = maxT - minT

      allocate(self % minT(size(maxFetOrder)))
      allocate(self % maxT(size(maxFetOrder)))
      self % minT = minT
      self % maxT = maxT

      allocate(self % FET_evalPoints(FET_evalPoints))
      deltaEval = (maxval(maxT) - minval(minT)) / real(FET_evalPoints)
      offset = minval(minT) + deltaEval / TWO
      do i = 1, FET_evalPoints
        self % FET_evalPoints(i) = offset
        offset = offset + deltaEval
      end do

      ! Allocate space and zero all bins
      allocate( self % bins(maxval(maxFetOrder) + 1, DIM2, size(maxFetOrder)))
      self % bins = ZERO

      ! Note the array padding to avoid false sharing
      allocate( self % parallelBins(maxval(maxFetOrder) + 1, self % nThreads, size(maxFetOrder)))
      self % parallelBins = ZERO

      ! Save size of memory
      self % N = N

      ! Assign memory id
      self % id = id

      ! Set batchN, cycles and batchSize to default values
      self % batchN    = 0
      self % cycles    = 0
      self % batchSize = 1

      if(present(batchSize)) then
        if(batchSize > 0) then
          self % batchSize = batchSize
        else
          call fatalError(Here,'Batch Size of: '// numToChar(batchSize) //' is invalid')
        end if
      end if

    else
      ! Allocate space and zero all bins
      allocate( self % bins(N, DIM2,1))
      self % bins = ZERO

      ! Note the array padding to avoid false sharing
      allocate( self % parallelBins(N + array_pad, self % nThreads,1))
      self % parallelBins = ZERO

      ! Save size of memory
      self % N = N

      ! Assign memory id
      self % id = id

      ! Set batchN, cycles and batchSize to default values
      self % batchN    = 0
      self % cycles    = 0
      self % batchSize = 1

      if(present(batchSize)) then
        if(batchSize > 0) then
          self % batchSize = batchSize
        else
          call fatalError(Here,'Batch Size of: '// numToChar(batchSize) //' is invalid')
        end if
      end if
    end if

  end subroutine init

  !!
  !! Deallocate memory and return to uninitialised state
  !!
  subroutine kill(self)
   class(scoreMemory), intent(inout) :: self

   if(allocated(self % bins)) deallocate(self % bins)
   if(allocated(self % parallelBins)) deallocate(self % parallelBins)
   self % N = 0
   self % nThreads = 0
   self % batchN = 0

  end subroutine kill

  !!
  !! Score a result on a given single bin under idx
  !!
  subroutine score_defReal(self, score, idx, t)
    class(scoreMemory), intent(inout) :: self
    real(defReal), intent(in)         :: score
    integer(longInt), intent(in)      :: idx
    real(defReal), intent(in), optional  :: t
    integer(shortInt)                 :: thread_idx
    character(100),parameter :: Here = 'score_defReal (scoreMemory_class.f90)'

    if (self % useFET .eqv. .true.) then
      call self % scoreFET(score, t)
    else
      ! Verify bounds for the index
      if( idx < 0_longInt .or. idx > self % N) then
        call fatalError(Here,'Index '//numToChar(idx)//' is outside bounds of &
                              & memory with size '//numToChar(self % N))
      end if

      ! Add the score
      thread_idx = ompGetThreadNum() + 1
      self % parallelBins(idx, thread_idx,1) = &
              self % parallelBins(idx, thread_idx,1) + score
    end if

  end subroutine score_defReal

  !!
  !! Score a result with shortInt on a given bin under idx
  !!
  subroutine score_shortInt(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    integer(shortInt), intent(in)     :: score
    integer(longInt), intent(in)      :: idx

    call self % score_defReal(real(score, defReal), idx)

  end subroutine score_shortInt

  !!
  !! Score a result with longInt on a given bin under idx
  !!
  subroutine score_longInt(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    integer(longInt), intent(in)      :: score
    integer(longInt), intent(in)     :: idx

    call self % score_defReal(real(score, defReal), idx)

  end subroutine score_longInt

  !!
  !! Increment the result directly on cumulative sums
  !!
  subroutine accumulate_defReal(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    real(defReal), intent(in)         :: score
    integer(longInt), intent(in)      :: idx
    character(100),parameter :: Here = 'accumulate_defReal (scoreMemory_class.f90)'

    ! Verify bounds for the index
    if( idx < 0_longInt .or. idx > self % N) then
      call fatalError(Here,'Index '//numToChar(idx)//' is outside bounds of &
                            & memory with size '//numToChar(self % N))
    end if

    ! Add the score
    self % bins(idx, CSUM,1)  = self % bins(idx, CSUM,1)  + score
    self % bins(idx, CSUM2,1) = self % bins(idx, CSUM2,1) + score * score

  end subroutine accumulate_defReal

  !!
  !! Increment the result directly on cumulative sums with shortInt score
  !!
  subroutine accumulate_shortInt(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    integer(shortInt), intent(in)     :: score
    integer(longInt), intent(in)     :: idx

    call self % accumulate_defReal(real(score, defReal), idx)

  end subroutine accumulate_shortInt

  !!
  !! Increment the result directly on cumulative sums with longInt score
  !!
  subroutine accumulate_longInt(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    integer(longInt), intent(in)      :: score
    integer(longInt), intent(in)      :: idx

    call self % accumulate_defReal(real(score, defReal), idx)

  end subroutine accumulate_longInt

  !!
  !! Close Cycle
  !! Increments cycle counter and detects end-of-batch
  !! When batch finishes it normalises all scores by the factor and moves them to CSUMs
  !!
  subroutine closeCycle(self, normFactor)
    class(scoreMemory), intent(inout) :: self
    real(defReal),intent(in)          :: normFactor
    integer(longInt)                  :: i
    integer(shortInt)                 :: k, j
    real(defReal), save               :: res
    !$omp threadprivate(res)

    if (self % useFET .eqv. .true.) then

      ! Increment Cycle Counter
      self % cycles = self % cycles + 1

      if(mod(self % cycles, self % batchSize) == 0) then ! Close Batch

        do j = 1, size(self % maxFetOrder)
          !$omp parallel do
          do k = 1, self % maxFetOrder(j) + 1
            
            ! Normalise scores
            !self % parallelBins(k,:) = self % parallelBins(k,:) * normFactor
            res = sum(self % parallelBins(k,:,j))
            
            ! Zero all score bins
            self % parallelBins(k,:,j) = ZERO
          
            ! Increment cumulative sums 
            self % bins(k,CSUM,j)  = self % bins(k, CSUM,j) + res
            self % bins(k,CSUM2,j) = self % bins(k, CSUM2,j) + res * res

            if (self % basis == 5) then
              res = sum(self % parallelBins_b(k,:,j))

              ! Zero all score bins
              self % parallelBins_b(k,:,j) = ZERO

              ! Increment cumulative sums 
              self % bins_b(k,CSUM,j)  = self % bins_b(k, CSUM,j) + res
              self % bins_b(k,CSUM2,j) = self % bins_b(k, CSUM2,j) + res * res
            end if

          end do
          !$omp end parallel do
          ! Increment batch counter
          !self % batchN = self % batchN + 1
        end do

      end if

    else
      ! Increment Cycle Counter
      self % cycles = self % cycles + 1

      if(mod(self % cycles, self % batchSize) == 0) then ! Close Batch
        
        !$omp parallel do
        do i = 1, self % N
          
          ! Normalise scores
          self % parallelBins(i,:,1) = self % parallelBins(i,:,1) * normFactor
          res = sum(self % parallelBins(i,:,1))
          
          ! Zero all score bins
          self % parallelBins(i,:,1) = ZERO
        
          ! Increment cumulative sums 
          self % bins(i,CSUM,1)  = self % bins(i,CSUM,1) + res
          self % bins(i,CSUM2,1) = self % bins(i,CSUM2,1) + res * res

        end do
        !$omp end parallel do

        ! Increment batch counter
        self % batchN = self % batchN + 1

      end if
    end if

  end subroutine closeCycle

  !!
  !! Close Cycle
  !! Multiplies score in bin under idx by normFactor, accumulates it and sets it to zero
  !!
  subroutine closeBin(self, normFactor, idx)
    class(scoreMemory), intent(inout) :: self
    real(defReal),intent(in)          :: normFactor
    integer(longInt), intent(in)      :: idx
    real(defReal)                     :: res
    character(100),parameter :: Here = 'closeBin (scoreMemory_class.f90)'

    ! Verify bounds for the index
    if( idx < 0_longInt .or. idx > self % N) then
      call fatalError(Here,'Index '//numToChar(idx)//' is outside bounds of &
                            & memory with size '//numToChar(self % N))
    end if

    ! Normalise score
    self % parallelBins(idx, :, 1) = self % parallelBins(idx, :, 1) * normFactor

    ! Increment cumulative sum
    res = sum(self % parallelBins(idx,:, 1))
    self % bins(idx,CSUM, 1)  = self % bins(idx,CSUM, 1) + res
    self % bins(idx,CSUM2, 1) = self % bins(idx,CSUM2, 1) + res * res

    ! Zero the score
    self % parallelBins(idx,:,1) = ZERO

  end subroutine closeBin


  !!
  !! Return true if next closeCycle will close a batch
  !!
  function lastCycle(self) result(isIt)
    class(scoreMemory), intent(in) :: self
    logical(defBool)               :: isIt

    isIt =  mod(self % cycles + 1, self % batchSize) == 0

  end function lastCycle

  !!
  !! Return batchSize
  !!
  pure function getBatchSize(self) result(S)
    class(scoreMemory), intent(in) :: self
    integer(shortInt)              :: S

    S = self % batchSize

  end function getBatchSize

  !!
  !! Set number of batches used per time step in ScoreMemory
  !!
  !! Args:
  !!   batchN
  !!
  !! Errors:
  !!   None
  !!
  subroutine setNumBatchesPerTimeStep(self, batchN) 
    class(scoreMemory), intent(inout) :: self
    integer(shortInt), intent(in)     :: batchN

    self % batchN = batchN

  end subroutine setNumBatchesPerTimeStep

  ! FET LEGENDRE FUNCTIONS

  subroutine scoreFET_signature(self, score, t)
    class(scoreMemory), intent(inout) :: self
    real(defReal), intent(in)         :: score, t
  end subroutine scoreFET_signature

  subroutine scoreFET_Legendre(self, score, t)
    class(scoreMemory),intent(inout)  :: self
    real(defReal), intent(in)         :: score, t
    integer(shortInt)                 :: k
    integer(shortInt)                 :: thread_idx, piecewise_idx
    real(defReal)                     :: t_trans, phi = ONE
    real(defReal)                     :: p_k, p_prev, p_curr, p_next
    character(100),parameter :: Here = 'scoreFET_Legendre (scoreMemory_class.f90)'

    thread_idx = ompGetThreadNum() + 1
    piecewise_idx = self % mapPiecewise(t)

    !time transform
    t_trans = self % transDomain_Legendre(t, piecewise_idx)

    !basis weight
    call self % weight_Legendre(t_trans, phi)
    do k = 0, self % maxFetOrder(piecewise_idx)

      ! Handle base cases
      if (k == 0) then
        p_k = 1.0_defReal
        
      else if (k == 1) then
        p_k = t_trans

        p_prev = 1.0_defReal
        p_curr = t_trans
        
      else
        p_next = ((2.0_defReal * real(k) - 1.0_defReal) * t_trans * p_curr - (real(k) - 1.0_defReal) * p_prev) / real(k)
        p_prev = p_curr
        p_curr = p_next
        p_k = p_curr

      end if

      ! Add the score
      self % parallelBins(k+1, thread_idx, piecewise_idx) = &
              self % parallelBins(k+1, thread_idx, piecewise_idx) + score * phi * p_k
    end do

  end subroutine scoreFET_Legendre

  function transDomain_Legendre(self, t, piecewise_idx) result(t_trans)
    class(scoreMemory),intent(in) :: self
    real(defReal), intent(in)     :: t
    integer(shortInt), intent(in) :: piecewise_idx
    real(defReal)                 :: t_trans

    t_trans = TWO * ((t - self % minT(piecewise_idx)) / self % deltaT(piecewise_idx)) - ONE

  end function transDomain_Legendre

  subroutine weight_Legendre(self, t_trans, phi) 
    class(scoreMemory), intent(in) :: self
    real(defReal), intent(in)      :: t_trans
    real(defReal), intent(inout)   :: phi

    phi = ONE

  end subroutine weight_Legendre

  function orthonormalisationTransformed_Legendre(self, k, piecewise_idx) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k, piecewise_idx
    real(defReal)                  :: km

    km = (TWO * k + ONE) / self % deltaT(piecewise_idx)

  end function orthonormalisationTransformed_Legendre

  function orthonormalisationStandard_Legendre(self, k) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k
    real(defReal)                  :: km

    km = (TWO * k + ONE) / TWO

  end function orthonormalisationStandard_Legendre

  subroutine getFETResult_signature(self, mean, STD, time, fet_coeff_arr, fet_coeff_std_arr)
    class(scoreMemory), intent(in)                           :: self
    real(defReal), intent(out)                               :: mean
    real(defReal), intent(out)                               :: STD
    real(defReal), intent(in)                                :: time
    real(defReal), dimension(maxval(self% maxFetOrder) + 1, size(self% maxFetOrder)), intent(in) :: fet_coeff_arr, fet_coeff_std_arr 
  end subroutine getFETResult_signature

  subroutine getFETResult_Legendre(self, mean, STD, time, fet_coeff_arr, fet_coeff_std_arr)
    class(scoreMemory), intent(in)                           :: self
    real(defReal), intent(out)                               :: mean
    real(defReal), intent(out)                               :: STD
    real(defReal), intent(in)                                :: time
    real(defReal), dimension(maxval(self% maxFetOrder) + 1, size(self% maxFetOrder)), intent(in) :: fet_coeff_arr, fet_coeff_std_arr 
    real(defReal)                                            :: t_trans
    real(defReal)                                            :: factor
    integer(shortInt)                                        :: k, piecewise_idx
    real(defReal)                  :: p_k, p_prev, p_curr, p_next
    character(100),parameter :: Here = 'getFETResult_Legendre (scoreMemory_class.f90)'

    piecewise_idx = self % mapPiecewise(time)
    t_trans = self % transDomain_Legendre(time, piecewise_idx)
    mean = ZERO
    STD = ZERO

    do k = 0, self % maxFetOrder(piecewise_idx)
      factor = fet_coeff_arr(k + 1, piecewise_idx) * self % orthonormalisationTransformed_Legendre(k, piecewise_idx)

      ! Handle base cases
      if (k == 0) then
        p_k = 1.0_defReal
        
      else if (k == 1) then
        p_k = t_trans

        p_prev = 1.0_defReal
        p_curr = t_trans

      else
        p_next = ((2.0_defReal * real(k) - 1.0_defReal) * t_trans * p_curr - (real(k) - 1.0_defReal) * p_prev) / real(k)
        p_prev = p_curr
        p_curr = p_next
        p_k = p_curr

      end if

      mean = mean + factor * p_k
      STD = STD + (fet_coeff_std_arr(k + 1, piecewise_idx)**TWO) * self % orthonormalisationStandard_Legendre(k)
    end do

    STD = SQRT(STD)

  end subroutine getFETResult_Legendre

  ! FET CHEBYSHEV1 FUNCTIONS

  subroutine scoreFET_Chebyshev1(self, score, t)
    class(scoreMemory),intent(inout) :: self
    real(defReal), intent(in)        :: score, t
    integer(shortInt)                :: k
    integer(shortInt)                :: thread_idx, piecewise_idx
    real(defReal)                    :: t_trans, phi = ONE
    real(defReal)                    :: p_k, p_prev, p_curr, p_next
    character(100),parameter         :: Here = 'scoreFET_Chebyshev1 (scoreMemory_class.f90)'

    thread_idx = ompGetThreadNum() + 1
    piecewise_idx = self % mapPiecewise(t)

    !time transform
    t_trans = self % transDomain_Chebyshev1(t, piecewise_idx)

    !basis weight
    call self % weight_Chebyshev1(t_trans, phi)

    do k = 0, self % maxFetOrder(piecewise_idx)

      ! Handle base cases
      if (k == 0) then
        p_k = 1.0_defReal
        
      else if (k == 1) then
        p_k = t_trans

        p_prev = 1.0_defReal
        p_curr = t_trans
        
      else
        p_next = 2.0_defReal * t_trans * p_curr - p_prev
        p_prev = p_curr
        p_curr = p_next
        p_k = p_curr

      end if

      ! Add the score
      self % parallelBins(k+1, thread_idx, piecewise_idx) = &
              self % parallelBins(k+1, thread_idx, piecewise_idx) + score * phi * p_k
    end do

  end subroutine scoreFET_Chebyshev1

  function transDomain_Chebyshev1(self, t, piecewise_idx) result(t_trans)
    class(scoreMemory),intent(in) :: self
    real(defReal), intent(in)     :: t
    integer(shortInt), intent(in)     :: piecewise_idx
    real(defReal)                 :: t_trans

    t_trans = TWO * ((t - self % minT(piecewise_idx)) / self % deltaT(piecewise_idx)) - ONE

  end function transDomain_Chebyshev1

  subroutine weight_Chebyshev1(self, t_trans, phi) 
    class(scoreMemory), intent(in) :: self
    real(defReal), intent(in)      :: t_trans
    real(defReal), intent(inout)   :: phi

    phi = ONE / sqrt(ONE - t_trans**TWO)

  end subroutine weight_Chebyshev1

  function orthonormalisationTransformed_Chebyshev1(self, k, piecewise_idx) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k
    integer(shortInt), intent(in)     :: piecewise_idx
    real(defReal)                  :: km

    if (k == 0) then
      km = TWO / (PI * self % deltaT(piecewise_idx))
    else
      km = TWO * TWO / (PI * self % deltaT(piecewise_idx))
    end if

  end function orthonormalisationTransformed_Chebyshev1

  function orthonormalisationStandard_Chebyshev1(self, k) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k
    real(defReal)                  :: km

    if (k == 0) then
      km = ONE / PI
    else
      km = TWO / PI
    end if

  end function orthonormalisationStandard_Chebyshev1

  subroutine getFETResult_Chebyshev1(self, mean, STD, time, fet_coeff_arr, fet_coeff_std_arr)
    class(scoreMemory), intent(in)                           :: self
    real(defReal), intent(out)                               :: mean
    real(defReal), intent(out)                               :: STD
    real(defReal), intent(in)                                :: time
    real(defReal), dimension(maxval(self% maxFetOrder + 1), size(self% maxFetOrder)), intent(in) :: fet_coeff_arr, fet_coeff_std_arr 
    real(defReal)                                            :: t_trans
    real(defReal)                                            :: factor
    integer(shortInt)                                        :: k, piecewise_idx
    real(defReal)                                            :: p_k, p_prev, p_curr, p_next
    character(100),parameter :: Here = 'getFETResult_Chebyshev1 (scoreMemory_class.f90)'

    piecewise_idx = self % mapPiecewise(time)
    t_trans = self % transDomain_Chebyshev1(time, piecewise_idx)
    mean = ZERO
    STD = ZERO
    do k = 0, self % maxFetOrder(piecewise_idx)
      factor = fet_coeff_arr(k + 1, piecewise_idx) * self % orthonormalisationTransformed_Chebyshev1(k, piecewise_idx)

      ! Handle base cases
      if (k == 0) then
        p_k = 1.0_defReal
        
      else if (k == 1) then
        p_k = t_trans

        p_prev = 1.0_defReal
        p_curr = t_trans

      else
        p_next = 2.0_defReal * t_trans * p_curr - p_prev
        p_prev = p_curr
        p_curr = p_next
        p_k = p_curr

      end if

      mean = mean + factor * p_k
      STD = STD + (fet_coeff_std_arr(k + 1, piecewise_idx)**TWO) * self % orthonormalisationStandard_Chebyshev1(k)
    end do

    STD = SQRT(STD)

  end subroutine getFETResult_Chebyshev1

  ! FET CHEBYSHEV2 FUNCTIONS

  subroutine scoreFET_Chebyshev2(self, score, t)
    class(scoreMemory),intent(inout)  :: self
    real(defReal), intent(in)         :: score, t
    integer(shortInt)                 :: k
    integer(shortInt)                 :: thread_idx, piecewise_idx
    real(defReal)                     :: t_trans, phi = ONE
    real(defReal)                     :: p_k, p_prev, p_curr, p_next
    character(100),parameter :: Here = 'scoreFET_Chebyshev2 (scoreMemory_class.f90)'

    thread_idx = ompGetThreadNum() + 1
    piecewise_idx = self % mapPiecewise(t)

    !time transform
    t_trans = self % transDomain_Chebyshev2(t, piecewise_idx)

    !basis weight
    call self % weight_Chebyshev2(t_trans, phi)

    do k = 0, self % maxFetOrder(piecewise_idx)

      ! Handle base cases
      if (k == 0) then
        p_k = 1.0_defReal
        
      else if (k == 1) then
        p_k = TWO * t_trans

        p_prev = 1.0_defReal
        p_curr = TWO * t_trans
        
      else
        p_next = 2.0_defReal * t_trans * p_curr - p_prev
        p_prev = p_curr
        p_curr = p_next
        p_k = p_curr

      end if

      ! Add the score
      self % parallelBins(k+1, thread_idx, piecewise_idx) = &
              self % parallelBins(k+1, thread_idx, piecewise_idx) + score * phi * p_k
    end do

  end subroutine scoreFET_Chebyshev2

  function transDomain_Chebyshev2(self, t, piecewise_idx) result(t_trans)
    class(scoreMemory),intent(in) :: self
    real(defReal), intent(in)     :: t
    integer(shortInt), intent(in)     :: piecewise_idx
    real(defReal)                 :: t_trans

    t_trans = TWO * ((t - self % minT(piecewise_idx)) / self % deltaT(piecewise_idx)) - ONE

  end function transDomain_Chebyshev2

  subroutine weight_Chebyshev2(self, t_trans, phi) 
    class(scoreMemory), intent(in) :: self
    real(defReal), intent(in)      :: t_trans
    real(defReal), intent(inout)   :: phi

    phi = sqrt(ONE - t_trans**TWO)

  end subroutine weight_Chebyshev2

  function orthonormalisationTransformed_Chebyshev2(self, k, piecewise_idx) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k
    integer(shortInt), intent(in)     :: piecewise_idx
    real(defReal)                  :: km

    km = TWO * TWO / (PI * self % deltaT(piecewise_idx))

  end function orthonormalisationTransformed_Chebyshev2

  function orthonormalisationStandard_Chebyshev2(self, k) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k
    real(defReal)                  :: km

    km = TWO / PI

  end function orthonormalisationStandard_Chebyshev2

  subroutine getFETResult_Chebyshev2(self, mean, STD, time, fet_coeff_arr, fet_coeff_std_arr)
    class(scoreMemory), intent(in)                           :: self
    real(defReal), intent(out)                               :: mean
    real(defReal), intent(out)                               :: STD
    real(defReal), intent(in)                                :: time
    real(defReal), dimension(maxval(self% maxFetOrder) + 1, size(self% maxFetOrder)), intent(in) :: fet_coeff_arr, fet_coeff_std_arr 
    real(defReal)                                            :: t_trans
    real(defReal)                                            :: factor
    integer(shortInt)                                        :: k, piecewise_idx
    real(defReal)                  :: p_k, p_prev, p_curr, p_next
    character(100),parameter :: Here = 'getFETResult_Chebyshev1 (scoreMemory_class.f90)'

    piecewise_idx = self % mapPiecewise(time)
    t_trans = self % transDomain_Chebyshev2(time, piecewise_idx)
    mean = ZERO
    STD = ZERO

    do k = 0, self % maxFetOrder(piecewise_idx)
      factor = fet_coeff_arr(k + 1, piecewise_idx) * self % orthonormalisationTransformed_Chebyshev2(k, piecewise_idx)

      ! Handle base cases
      if (k == 0) then
        p_k = 1.0_defReal
        
      else if (k == 1) then
        p_k = TWO * t_trans

        p_prev = 1.0_defReal
        p_curr = TWO * t_trans

      else
        p_next = 2.0_defReal * t_trans * p_curr - p_prev
        p_prev = p_curr
        p_curr = p_next
        p_k = p_curr

      end if

      mean = mean + factor * p_k
      STD = STD + (fet_coeff_std_arr(k + 1, piecewise_idx)**TWO) * self % orthonormalisationStandard_Chebyshev2(k)
    end do

    STD = SQRT(STD)

  end subroutine getFETResult_Chebyshev2

  ! FET FOURIER FUNCTIONS

  subroutine scoreFET_Fourier(self, score, t)
    class(scoreMemory),intent(inout)  :: self
    real(defReal), intent(in)         :: score, t
    integer(shortInt)                 :: k
    integer(shortInt)                 :: thread_idx, piecewise_idx
    real(defReal)                     :: t_trans, phi = ONE
    real(defReal)                     :: p_k, p_prev, p_curr, p_next
    character(100),parameter :: Here = 'scoreFET_Fourier (scoreMemory_class.f90)'

    thread_idx = ompGetThreadNum() + 1
    piecewise_idx = self % mapPiecewise(t)

    !time transform
    t_trans = self % transDomain_Fourier(t, piecewise_idx)

    !basis weight
    call self % weight_Fourier(t_trans, phi)

    do k = 0, self % maxFetOrder(piecewise_idx)

      p_k = cos(k * t_trans)
      ! Add the score
      self % parallelBins(k+1, thread_idx, piecewise_idx) = &
              self % parallelBins(k+1, thread_idx, piecewise_idx) + score * phi * p_k

      p_k = sin(k * t_trans)
      ! Add the score
      self % parallelBins_b(k+1, thread_idx, piecewise_idx) = &
              self % parallelBins_b(k+1, thread_idx, piecewise_idx) + score * phi * p_k
    end do

  end subroutine scoreFET_Fourier

  function transDomain_Fourier(self, t, piecewise_idx) result(t_trans)
    class(scoreMemory),intent(in) :: self
    real(defReal), intent(in)     :: t
    integer(shortInt), intent(in)     :: piecewise_idx
    real(defReal)                 :: t_trans

    t_trans = TWO_PI*((t - self % minT(piecewise_idx)) / self % deltaT(piecewise_idx)) - PI

  end function transDomain_Fourier

  subroutine weight_Fourier(self, t_trans, phi) 
    class(scoreMemory), intent(in) :: self
    real(defReal), intent(in)      :: t_trans
    real(defReal), intent(inout)   :: phi
    real(defReal)                  :: t

    phi = ONE

  end subroutine weight_Fourier

  function orthonormalisationTransformed_Fourier(self, k, piecewise_idx) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k, piecewise_idx
    real(defReal)                  :: km

    if (k == 0) then
      km = ONE / (self % maxT(piecewise_idx) - self % minT(piecewise_idx))
    else
      km = TWO / (self % maxT(piecewise_idx) - self % minT(piecewise_idx))
    end if

  end function orthonormalisationTransformed_Fourier

  function orthonormalisationStandard_Fourier(self, k) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k
    real(defReal)                  :: km

    if (k == 0) then
      km = ONE / TWO_PI
    else
      km = ONE / PI
    end if

  end function orthonormalisationStandard_Fourier

  function orthonormalisationTransformed_b_Fourier(self, k, piecewise_idx) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k, piecewise_idx
    real(defReal)                  :: km
    character(100),parameter :: Here = 'orthonormalisationTransformed_b_Fourier (scoreMemory_class.f90)'

    if (k == 0) then
      call fatalError(Here, 'k0 not defined for sine basis')
    else
      km = TWO / (self % maxT(piecewise_idx) - self % minT(piecewise_idx))
    end if

  end function orthonormalisationTransformed_b_Fourier

  function orthonormalisationStandard_b_Fourier(self, k) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k
    real(defReal)                  :: km
    character(100),parameter :: Here = 'orthonormalisationStandard_b_Fourier (scoreMemory_class.f90)'

    if (k == 0) then
      call fatalError(Here, 'k0 not defined for sine basis')
    else
      km = ONE / PI
    end if

  end function orthonormalisationStandard_b_Fourier

  subroutine getFETResult_Fourier(self, mean, STD, time, fet_coeff_arr, fet_coeff_arr_b, fet_coeff_std_arr, fet_coeff_std_arr_b)
    class(scoreMemory), intent(in)                           :: self
    real(defReal), intent(out)                               :: mean
    real(defReal), intent(out)                               :: STD
    real(defReal), intent(in)                                :: time
    real(defReal), dimension(maxval(self% maxFetOrder) + 1,&
         size(self% maxFetOrder)), intent(in) :: fet_coeff_arr, fet_coeff_std_arr 
    real(defReal), dimension(maxval(self% maxFetOrder) + 1,&
         size(self% maxFetOrder)), intent(in) :: fet_coeff_arr_b, fet_coeff_std_arr_b
    real(defReal)                                            :: t_trans
    real(defReal)                                            :: factor
    integer(shortInt)                                        :: k, piecewise_idx
    real(defReal)                  :: p_k, p_prev, p_curr, p_next
    character(100),parameter :: Here = 'getFETResult_Fourier (scoreMemory_class.f90)'

    piecewise_idx = self % mapPiecewise(time)
    t_trans = self % transDomain_Fourier(time, piecewise_idx)

    mean = ZERO
    STD = ZERO

    do k = 0, self % maxFetOrder(piecewise_idx)
      factor = fet_coeff_arr(k + 1, piecewise_idx) * self % orthonormalisationTransformed_Fourier(k, piecewise_idx)
      p_k = cos(k * t_trans)
      mean = mean + factor * p_k
      STD = STD + (fet_coeff_std_arr(k + 1, piecewise_idx)**TWO) * self % orthonormalisationStandard_Fourier(k)

      if (k > 0) then
        factor = fet_coeff_arr_b(k + 1, piecewise_idx) * self % orthonormalisationTransformed_b_Fourier(k, piecewise_idx)
        p_k = sin(k * t_trans)
        mean = mean + factor * p_k
        STD = STD + (fet_coeff_std_arr_b(k + 1, piecewise_idx)**TWO) * self % orthonormalisationStandard_b_Fourier(k)
      end if

    end do

    STD = SQRT(STD)

  end subroutine getFETResult_Fourier

  ! FET JACOBI FUNCTIONS

  subroutine scoreFET_Jacobi(self, score, t)
    class(scoreMemory),intent(inout)  :: self
    real(defReal), intent(in)         :: score, t
    integer(shortInt)                 :: k
    integer(shortInt)                 :: thread_idx, piecewise_idx
    real(defReal)                     :: t_trans, phi = ONE
    real(defReal)                     :: p_k, p_prev, p_curr, p_next
    real(defReal)                     :: an, bn, cn, dn
    character(100),parameter :: Here = 'scoreFET_Jacobi (scoreMemory_class.f90)'

    thread_idx = ompGetThreadNum() + 1
    piecewise_idx = self % mapPiecewise(t)

    !time transform
    t_trans = self % transDomain_Jacobi(t) 

    !basis weight
    call self % weight_Jacobi(t_trans, phi)

    do k = 0, self % maxFetOrder(piecewise_idx)

      ! Handle base cases
      if (k == 0) then
        p_k = ONE
        
      else if (k == 1) then
        p_k = HALF * (self % a + self % b + TWO)*t_trans + HALF*(self%a-self%b)
        p_prev = ONE
        p_curr = HALF * (self % a + self % b + TWO)*t_trans + HALF*(self%a-self%b)

      else
        an = TWO * real(k) * (real(k) + self % a + self % b) * (TWO*real(k) - TWO + self % a + self % b)
        bn = (TWO*real(k) - ONE + self % a + self % b) * &
             ((TWO*real(k) - TWO + self % a + self % b)*(TWO*real(k) + self % a + self % b)*t_trans + self % a**TWO - self % b**TWO)
        cn = TWO * (real(k) - ONE + self % a) * (real(k) - ONE + self % b) * (TWO * real(k) + self % a + self % b)

        p_next = (bn * p_curr - cn * p_prev) / an

        !p_next = ( (TWO*real(k)-ONE)*((TWO*real(k)-TWO)*TWO*real(k)*t_trans)*p_curr &
        !          - TWO*(real(k)-ONE)*(real(k)-ONE)*(TWO*real(k))*p_prev  ) / (TWO*real(k)*real(k)*(TWO*real(k)-TWO))

        p_prev = p_curr
        p_curr = p_next
        p_k = p_curr

      end if

      ! Add the score
      self % parallelBins(k+1, thread_idx, piecewise_idx) = &
              self % parallelBins(k+1, thread_idx, piecewise_idx) + score * phi * p_k
    end do

  end subroutine scoreFET_Jacobi

  function transDomain_Jacobi(self, t) result(t_trans)
    class(scoreMemory),intent(in) :: self
    real(defReal), intent(in)     :: t
    real(defReal)                 :: t_trans
    integer(shortInt)             :: piecewise_idx

    piecewise_idx = self % mapPiecewise(t)
    t_trans = TWO * ((t - self % minT(piecewise_idx)) / self % deltaT(piecewise_idx)) - ONE

  end function transDomain_Jacobi

  subroutine weight_Jacobi(self, t_trans, phi) 
    class(scoreMemory), intent(in) :: self
    real(defReal), intent(in)      :: t_trans
    real(defReal), intent(inout)   :: phi

    phi = ((ONE-t_trans)**self%a)*((ONE+t_trans)**self%b)

  end subroutine weight_Jacobi

  function orthonormalisationTransformed_Jacobi(self, k, piecewise_idx) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k, piecewise_idx
    real(defReal)                  :: km

    km = self % orthonormalisationStandard_Jacobi(k) * TWO / self % deltaT(piecewise_idx)

  end function orthonormalisationTransformed_Jacobi

  function orthonormalisationStandard_Jacobi(self, k) result(km)
    class(scoreMemory), intent(in) :: self
    integer(shortInt), intent(in)  :: k
    real(defReal)                  :: km, num, denom
    character(100),parameter :: Here = 'orthonormalisationStandard_Jacobi (scoreMemory_class.f90)'

    if ((k == 0) .and. (k + self % a + self % b <= -1)) then

      if ((self % a == -0.5) .and. (self % b == -0.5)) then
        km = ONE / PI
        return
      else
        call fatalError(Here, 'Undefined Gamma! Choose other a,b for Jacobian')
      end if

    else
      num = (TWO**(self % a + self % b + ONE)) * gamma(k + self % a + ONE) * gamma(k + self % b + ONE)
      denom = (TWO * k + self % a + self % b + ONE) * gamma(k + ONE) * gamma(k + self % a + self % b + ONE)
    end if

    km = denom / num

  end function orthonormalisationStandard_Jacobi

  subroutine getFETResult_Jacobi(self, mean, STD, time, fet_coeff_arr, fet_coeff_std_arr)
    class(scoreMemory), intent(in)                           :: self
    real(defReal), intent(out)                               :: mean
    real(defReal), intent(out)                               :: STD
    real(defReal), intent(in)                                :: time
    real(defReal), dimension(maxval(self% maxFetOrder) + 1, size(self% maxFetOrder)), intent(in) :: fet_coeff_arr, fet_coeff_std_arr
    real(defReal)                                            :: t_trans
    real(defReal)                                            :: factor
    integer(shortInt)                                        :: k, piecewise_idx
    real(defReal)                  :: p_k, p_prev, p_curr, p_next
    real(defReal)                     :: an, bn, cn, dn
    character(100),parameter :: Here = 'getFETResult_Jacobi (scoreMemory_class.f90)'

    piecewise_idx = self % mapPiecewise(time)
    t_trans = self % transDomain_Jacobi(time)

    mean = ZERO
    STD = ZERO

    do k = 0, self % maxFetOrder(piecewise_idx)
      factor = fet_coeff_arr(k + 1, piecewise_idx) * self % orthonormalisationTransformed_Jacobi(k, piecewise_idx)

      ! Handle base cases
      if (k == 0) then
        p_k = ONE
        
      else if (k == 1) then
        p_k = HALF * (self % a + self % b + TWO)*t_trans + HALF*(self%a-self%b)
        p_prev = ONE
        p_curr = HALF * (self % a + self % b + TWO)*t_trans + HALF*(self%a-self%b)

      else
        an = TWO * real(k) * (real(k) + self % a + self % b) * (TWO*real(k) - TWO + self % a + self % b)
        bn = (TWO*real(k) - ONE + self % a + self % b) * &
             ((TWO*real(k) - TWO + self % a + self % b)*(TWO*real(k) + self % a + self % b)*t_trans + self % a**TWO - self % b**TWO)
        cn = TWO * (real(k) - ONE + self % a) * (real(k) - ONE + self % b) * (TWO * real(k) + self % a + self % b)

        p_next = (bn * p_curr - cn * p_prev) / an

        !p_next = ( (TWO*real(k)-ONE)*((TWO*real(k)-TWO)*TWO*real(k)*t_trans)*p_curr &
        !          - TWO*(real(k)-ONE)*(real(k)-ONE)*(TWO*real(k))*p_prev  ) / (TWO*real(k)*real(k)*(TWO*real(k)-TWO))

        p_prev = p_curr
        p_curr = p_next
        p_k = p_curr

      end if

      mean = mean + factor * p_k
      STD = STD + (fet_coeff_std_arr(k + 1, piecewise_idx)**TWO) * self % orthonormalisationStandard_Jacobi(k)
    end do

    STD = SQRT(STD)

  end subroutine getFETResult_Jacobi

  !!
  !! Load mean result and Standard deviation into provided arguments
  !! Load from bin indicated by idx
  !! Returns 0 if index is invalid
  !!
  elemental subroutine getResult_withSTD(self, mean, STD, idx, samples, piecewise_idx)
    class(scoreMemory), intent(in)         :: self
    real(defReal), intent(out)             :: mean
    real(defReal),intent(out)              :: STD
    integer(longInt), intent(in)           :: idx
    integer(shortInt), intent(in),optional :: samples
    integer(shortInt), intent(in),optional :: piecewise_idx
    integer(shortInt)                      :: N
    real(defReal)                          :: inv_N, inv_Nm1

    ! Check if # of samples is provided
    if( present(samples)) then
      N = samples
    else
      N = self % batchN
    end if

    ! Calculate mean
    mean = self % bins(idx, CSUM, 1) / N

    ! Calculate STD
    inv_N   = ONE / N
    if( N /= 1) then
      inv_Nm1 = ONE / (N - 1)
    else
      inv_Nm1 = ONE
    end if
    STD = self % bins(idx, CSUM2, 1) *inv_N * inv_Nm1 - mean * mean * inv_Nm1
    STD = sqrt(STD)


  end subroutine getResult_withSTD


  !!
  !! Load mean result and Standard deviation into provided arguments
  !! Load from bin indicated by idx
  !! Returns 0 if index is invalid
  !!
  elemental subroutine getResult_FET(self, mean, STD, idx, samples, piecewise_idx)
    class(scoreMemory), intent(in)         :: self
    real(defReal), intent(out)             :: mean
    real(defReal),intent(out)              :: STD
    integer(longInt), intent(in)           :: idx
    integer(shortInt), intent(in),optional :: samples
    integer(longInt), intent(in),optional :: piecewise_idx
    integer(shortInt)                      :: N
    real(defReal)                          :: inv_N, inv_Nm1

    ! Check if # of samples is provided
    if( present(samples)) then
      N = samples
    else
      N = self % batchN
    end if

    ! Calculate mean
    mean = self % bins(idx, CSUM, piecewise_idx) / N

    ! Calculate STD
    inv_N   = ONE / N
    if( N /= 1) then
      inv_Nm1 = ONE / (N - 1)
    else
      inv_Nm1 = ONE
    end if
    STD = self % bins(idx, CSUM2, piecewise_idx) *inv_N * inv_Nm1 - mean * mean * inv_Nm1
    STD = sqrt(STD)


  end subroutine getResult_FET

  elemental subroutine getResult_Fourier_b(self, mean, STD, idx, samples, piecewise_idx)
    class(scoreMemory), intent(in)         :: self
    real(defReal), intent(out)             :: mean
    real(defReal),intent(out)              :: STD
    integer(longInt), intent(in)           :: idx
    integer(shortInt), intent(in),optional :: samples
    integer(longInt), intent(in),optional :: piecewise_idx
    integer(shortInt)                      :: N
    real(defReal)                          :: inv_N, inv_Nm1

    !! Verify index. Return 0 if not present
    !if( idx - 1 < 0_longInt .or. idx - 1 > self % maxFetOrder) then
    !  mean = ZERO
    !  STD = ZERO
    !  return
    !end if

    ! Check if # of samples is provided
    if( present(samples)) then
      N = samples
    else
      N = self % batchN
    end if

    ! Calculate mean
    mean = self % bins_b(idx, CSUM, piecewise_idx) / N

    ! Calculate STD
    inv_N   = ONE / N
    if( N /= 1) then
      inv_Nm1 = ONE / (N - 1)
    else
      inv_Nm1 = ONE
    end if
    STD = self % bins_b(idx, CSUM2, piecewise_idx) *inv_N * inv_Nm1 - mean * mean * inv_Nm1
    STD = sqrt(STD)

  end subroutine getResult_Fourier_b

  !!
  !! Load mean result provided argument
  !! Load from bin indicated by idx
  !! Returns 0 if index is invalid
  !!
  elemental subroutine getResult_withoutSTD(self, mean, idx, samples)
    class(scoreMemory), intent(in)         :: self
    real(defReal), intent(out)             :: mean
    integer(longInt), intent(in)           :: idx
    integer(shortInt), intent(in),optional :: samples
    integer(shortInt)                      :: N

    !! Verify index. Return 0 if not present
    if( idx < 0_longInt .or. idx > self % N) then
      mean = ZERO
      return
    end if

    ! Check if # of samples is provided
    if( present(samples)) then
      N = samples
    else
      N = self % batchN
    end if

    ! Calculate mean
    mean = self % bins(idx, CSUM,1) / N

  end subroutine getResult_withoutSTD

  !!
  !! Obtain value of a score in a bin
  !! Return ZERO for invalid bin address (idx)
  !!
  elemental function getScore(self, idx) result (score)
    class(scoreMemory), intent(in) :: self
    integer(longInt), intent(in)   :: idx
    real(defReal)                  :: score

    if(idx <= 0_longInt .or. idx > self % N) then
      score = ZERO
    else
      score = sum(self % parallelBins(idx, :, 1))
    end if

  end function getScore

  function mapPiecewise(self, t) result(piecewise_idx)
    class(scoreMemory), intent(in) :: self
    real(defReal), intent(in)      :: t
    integer(shortInt)              :: piecewise_idx

    if (t <= self % piecewise(1)) then
      piecewise_idx = 1_shortInt
    else if (t > self % piecewise(1) .and. t <= self % piecewise(2)) then 
      piecewise_idx = 2_shortInt
    else if (t > self % piecewise(2) .and. t <= self % piecewise(3)) then 
      piecewise_idx = 3_shortInt
    else if (t > self % piecewise(3)) then 
      piecewise_idx = 4_shortInt
    else
      piecewise_idx = 1_shortInt
    end if

  end function mapPiecewise

end module scoreMemory_class
