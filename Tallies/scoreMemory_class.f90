module scoreMemory_class

  use numPrecision
  use universalVariables, only : array_pad
  use genericProcedures,  only : fatalError, numToChar
  use openmp_func,        only : ompGetMaxThreads, ompGetThreadNum
  use RNG_class,                      only : RNG

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
      real(defReal),dimension(:,:),allocatable     :: bins              !! Space for storing cumul data (2nd dim size is always 2!)
      real(defReal),dimension(:,:),allocatable     :: parallelBins      !! Space for scoring for different threads
      integer(longInt)                             :: N = 0             !! Size of memory (number of bins)
      integer(shortInt)                            :: nThreads = 0      !! Number of threads used for parallelBins
      integer(shortInt)                            :: id                !! Id of the tally
      integer(shortInt)                            :: batchN = 0        !! Number of Batches
      integer(shortInt)                            :: cycles = 0        !! Cycles counter
      integer(shortInt)                            :: batchSize = 1     !! Batch interval size (in cycles)

      !real(defReal), dimension(:,:), allocatable   :: bootstrapBins
      !real(defReal), dimension(:), allocatable     :: plugInMean
      !real(defReal), dimension(:), allocatable     :: plugInVar
      integer(shortInt)                            :: nBootstraps
      real(defReal), dimension(:,:), allocatable   :: plugInSamples

      !real(defReal), dimension(:,:), allocatable   :: plugInSamples, plugInSamplesMean, plugInSamplesVar
      !real(defReal), dimension(:), allocatable     :: bootstrapAccumulator, unbiasedMeans, unbiasedVars
      !real(defReal), dimension(:), allocatable     ::  biasedMeans, biasedVars, normBias
      !integer(shortInt), dimension(:), allocatable :: scoreBinTracker
      integer(shortInt)                            :: nTimeBins, cyclesPerTime
      integer(shortInt)                            :: Ntallies
      integer(shortInt)                            :: bootstrapV = 0
      integer(longInt)                             :: timeN

      real(defReal), dimension(:), allocatable   :: bootstrapMean, bootstrapVar


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
    procedure :: reportTimeEnd
    procedure :: bootstrapV1
    procedure :: bootstrapV2
    procedure :: bootstrapV3
    ! Private procedures
    procedure, private :: score_defReal
    procedure, private :: score_shortInt
    procedure, private :: score_longInt
    procedure, private :: accumulate_defReal
    procedure, private :: accumulate_shortInt
    procedure, private :: accumulate_longInt
    procedure, private :: getResult_withSTD
    procedure, private :: getResult_withoutSTD

  end type scoreMemory

contains

  !!
  !! Allocate space for the bins given number of bins N
  !! Optionaly change batchSize from 1 to any +ve number
  !!
  subroutine init(self, N, id, batchSize, timeSteps, CyclesPerTime, nBootstraps, bootstrapV)
    class(scoreMemory),intent(inout)      :: self
    integer(longInt),intent(in)           :: N
    integer(shortInt),intent(in)          :: id
    integer(shortInt),optional,intent(in) :: batchSize, timeSteps, CyclesPerTime, nBootstraps, bootstrapV
    character(100), parameter :: Here= 'init (scoreMemory_class.f90)'

    self % nThreads = ompGetMaxThreads()

    ! Note the array padding to avoid false sharing
    allocate( self % parallelBins(N + array_pad, self % nThreads))
    self % parallelBins = ZERO

    ! Save size of memory
    self % N = N

    ! Assign memory id
    self % id = id

    ! Set batchN, cycles and batchSize to default values
    self % batchN    = 1
    self % cycles    = 0
    self % batchSize = 1

    if(present(batchSize)) then
      if(batchSize > 0) then
        self % batchSize = batchSize
      else
        call fatalError(Here,'Batch Size of: '// numToChar(batchSize) //' is invalid')
      end if
    end if

    self % nTimeBins = timeSteps
    self % Ntallies = N / self % nTimeBins
    self % TimeN = 1_longInt
    self % CyclesPerTime = CyclesPerTime


    !TODO: not include for bootstrap
    ! Allocate space and zero all bins
    allocate(self % bins(N, DIM2))
    self % bins = ZERO

    if (present(bootstrapV)) then !bootstrap score
      self % bootstrapV = bootstrapV
      self % nBootstraps = nBootstraps
    
      allocate(self % plugInSamples(self % Ntallies, CyclesPerTime))
      self % plugInSamples = ZERO

      allocate(self % bootstrapMean(self % N))
      allocate(self % bootstrapVar(self % N))
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
  subroutine score_defReal(self, score, idx)
    class(scoreMemory), intent(inout) :: self
    real(defReal), intent(in)         :: score
    integer(longInt), intent(in)      :: idx 
    integer(shortInt),save                 :: thread_idx
    real(defReal), save                     :: ratio
    character(100),parameter :: Here = 'score_defReal (scoreMemory_class.f90)'
    !$omp threadprivate(ratio)

    ! Verify bounds for the index
    if( idx < 0_longInt .or. idx > self % N) then
      call fatalError(Here,'Index '//numToChar(idx)//' is outside bounds of &
                            & memory with size '//numToChar(self % N))
    end if

    thread_idx = ompGetThreadNum() + 1

    ! Add the score
    self % parallelBins(idx, thread_idx) = &
            self % parallelBins(idx, thread_idx) + score

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
    self % bins(idx, CSUM)  = self % bins(idx, CSUM)  + score
    self % bins(idx, CSUM2) = self % bins(idx, CSUM2) + score * score

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
    class(scoreMemory), intent(inout)       :: self
    real(defReal),intent(in)                :: normFactor
    integer(longInt)                        :: i
    integer(longInt)                        :: loc
    real(defReal), save                     :: res
    !$omp threadprivate(res)

    ! Increment Cycle Counter
    self % cycles = self % cycles + 1

    loc = (self % timeN - 1_longInt) * self % Ntallies
    !$omp parallel do
    do i = 1, self % Ntallies

      ! Normalise scores
      self % parallelBins(loc + i,:) = self % parallelBins(loc + i,:) * normFactor
      res = sum(self % parallelBins(loc + i,:))

      ! Zero all score bins
      self % parallelBins(loc + i,:) = ZERO

      !TODO: not include for bootstrap
      ! Increment cumulative sums
      self % bins(loc + i,CSUM)  = self % bins(loc + i,CSUM) + res
      self % bins(loc + i,CSUM2) = self % bins(loc + i,CSUM2) + res * res

      if (self % bootstrapV > 0) then
        self % plugInSamples(i, self % batchN) = res
      end if
    end do
    !$omp end parallel do

    ! Increment batch counter
    self % batchN = self % batchN + 1
    !self % batchSize = self % batchSize + 1


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
    self % parallelBins(idx, :) = self % parallelBins(idx, :) * normFactor

    ! Increment cumulative sum
    res = sum(self % parallelBins(idx,:))
    self % bins(idx,CSUM)  = self % bins(idx,CSUM) + res
    self % bins(idx,CSUM2) = self % bins(idx,CSUM2) + res * res

    ! Zero the score
    self % parallelBins(idx,:) = ZERO

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

    !self % batchN = batchN

  end subroutine setNumBatchesPerTimeStep

  !!
  !! Load mean result and Standard deviation into provided arguments
  !! Load from bin indicated by idx
  !! Returns 0 if index is invalid
  !!
  elemental subroutine getResult_withSTD(self, mean, STD, idx, samples)
    class(scoreMemory), intent(in)         :: self
    real(defReal), intent(out)             :: mean
    real(defReal),intent(out)              :: STD
    integer(longInt), intent(in)           :: idx
    integer(shortInt), intent(in),optional :: samples
    integer(shortInt)                      :: N, i
    real(defReal)                          :: inv_N, inv_Nm1

    !! Verify index. Return 0 if not present
    if( idx < 0_longInt .or. idx > self % N) then
      mean = ZERO
      STD = ZERO
      return
    end if

    ! Check if # of samples is provided
    if( present(samples)) then
      N = samples
    else
      N = self % cyclesPerTime !self % batchN
    end if

    ! Calculate mean
    mean = (self % bins(idx, CSUM) / N)
    !
    !! Calculate STD
    inv_N   = ONE / N
    if( N /= 1) then
      inv_Nm1 = ONE / (N - 1)
    else
      inv_Nm1 = ONE
    end if
    STD = self % bins(idx, CSUM2) * inv_N * inv_Nm1 - mean * mean * inv_Nm1
    STD = sqrt(STD)
  end subroutine getResult_withSTD

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
      N = self % cyclesPerTime !self % batchN
    end if

    ! Calculate mean
    mean = self % bins(idx, CSUM) / N

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
      score = sum(self % parallelBins(idx, :))
    end if

  end function getScore

  subroutine reportTimeEnd(self, rand)
    class(scoreMemory), intent(inout) :: self
    class(RNG), intent(inout)         :: rand

    !perform bootstrapping based on version.
    if (self % bootstrapV == 1) then
      call self % bootstrapV1(rand)
      self % plugInSamples = ZERO

    else if (self % bootstrapV == 2) then
      call self % bootstrapV2(rand)
      self % plugInSamples = ZERO

    else if (self % bootstrapV == 3) then
      call self % bootstrapV3(rand)
      self % plugInSamples = ZERO
    end if

    self % batchN = 1
    self % timeN = self % timeN + 1_longInt
  end subroutine reportTimeEnd


  !Bootstrap means
  subroutine bootstrapV1(self, rand)
    class(scoreMemory), intent(inout)        :: self
    class(RNG), intent(inout)                :: rand
    integer(shortInt)                        :: i, N, j, r, B
    integer(longInt)                         :: loc
    real(defReal)                            :: inv_Bm1
    real(defReal), save                      :: res, plugInMean
    integer(shortInt), save                  :: idx
    real(defReal), save                      :: bootstrapAccumulator
    real(defReal), dimension(:), allocatable, save :: bootstrapMean, bootstrapVar
    real(defReal), save                      :: mean, var

    N = self % cyclesPerTime
    B = self % nBootstraps

    if (B /= 1) then
      inv_Bm1 = ONE / (B - 1)
    else
      inv_Bm1 = ONE
    end if

    allocate(bootstrapMean(self % Ntallies))
    allocate(bootstrapVar(self % Ntallies))
    bootstrapAccumulator = ZERO
    bootstrapMean = ZERO
    bootstrapVar = ZERO

    loc = (self % timeN - 1_longInt) * self % Ntallies

    do i = 1, 1 !self % Ntallies
      !$omp parallel do
      do j = 1, self % CyclesPerTime
        write(10, '(F24.6, ",")') self % plugInSamples(i,j)
      end do 
      !$omp end parallel do
    end do

    !$omp parallel do private(mean, var, bootstrapAccumulator, plugInMean, bootstrapMean, bootstrapVar)
    do i = 1, self % Ntallies

      plugInMean = sum(self % plugInSamples(i,:)) / N
      bootstrapMean(i) = plugInMean
      bootstrapVar(i) = plugInMean * plugInMean

      !$omp parallel do private(res, bootstrapAccumulator)
      do r = 1, self % nBootstraps - 1
        !$omp parallel do private(idx)
        do j = 1, N
          call rand % stride(j)
          validSample: do
            idx = int(N * rand % get()) + 1
            if (idx <= N) then
              bootstrapAccumulator = bootstrapAccumulator + self % plugInSamples(i,idx)
              exit validSample
            else
              cycle validSample
            end if

          end do validSample
        end do
        !$omp end parallel do

        res = bootstrapAccumulator
        bootstrapAccumulator = ZERO
        bootstrapMean(i) = bootstrapMean(i) + res / N
        bootstrapVar(i) = bootstrapVar(i) + (res / N) * (res / N)
      end do
      !$omp end parallel do

      mean = bootstrapMean(i) / self % nBootstraps
      var = bootstrapVar(i) * inv_Bm1 - mean * mean * inv_Bm1 * B

      self % bootstrapMean(loc + i) = mean
      self % bootstrapVar(loc + i) = var

    end do
    !$omp end parallel do

    deallocate(bootstrapMean)
    deallocate(bootstrapVar)

  end subroutine bootstrapV1


  !Bootstrap unbiased var
  subroutine bootstrapV2(self, rand)
    class(scoreMemory), intent(inout)        :: self
    class(RNG), intent(inout)                :: rand
    integer(shortInt)                        :: i, N, j, r, B
    integer(longInt)                         :: loc
    real(defReal)                            :: inv_Bm1, inv_N, inv_Nm1
    real(defReal), save                      :: res, res_sq, pluginVar
    integer(shortInt), save                  :: idx
    real(defReal), save                      :: bootstrapAccumulator, bootstrapAccumulator_sq
    real(defReal), dimension(:), allocatable, save :: bootstrapMean
    real(defReal), save                      :: mean, mean_biasAdj

    N = self % cyclesPerTime
    B = self % nBootstraps

    if (B /= 1) then
      inv_Bm1 = ONE / (B - 1)
    else
      inv_Bm1 = ONE
    end if

    inv_N   = ONE / N
    if (N /= 1) then
      inv_Nm1 = ONE / (N - 1)
    else
      inv_Nm1 = ONE
    end if


    allocate(bootstrapMean(self % Ntallies))
    bootstrapAccumulator = ZERO
    bootstrapAccumulator_sq = ZERO
    bootstrapMean = ZERO

    loc = (self % timeN - 1_longInt) * self % Ntallies


    !$omp parallel do 
    do i=1, size(self % plugInSamples(1,:))
      write(10, '(F24.16, ",")') self % plugInSamples(1,i)
    end do
    !$omp end parallel do 

    !$omp parallel do private(mean, mean_biasAdj, bootstrapAccumulator, bootstrapAccumulator_sq, plugInVar, bootstrapMean)
    do i = 1, self % Ntallies

      plugInVar = self % bins(loc + i, CSUM2) * inv_N * inv_Nm1 &
      - (self % bins(loc + i, CSUM) / N) * (self % bins(loc + i, CSUM) / N) * inv_Nm1

      bootstrapMean(i) = plugInVar

      !$omp parallel do private(res, res_sq, bootstrapAccumulator, bootstrapAccumulator_sq)
      do r = 1, self % nBootstraps - 1
        !$omp parallel do private(idx)
        do j = 1, N
          call rand % stride(j)
          validSample: do
            idx = int(N * rand % get()) + 1
            if (idx <= N) then
              bootstrapAccumulator = bootstrapAccumulator + self % plugInSamples(i,idx)
              bootstrapAccumulator_sq = bootstrapAccumulator_sq + self % plugInSamples(i,idx) * self % plugInSamples(i,idx)
              exit validSample
            else
              cycle validSample
            end if

          end do validSample
        end do
        !$omp end parallel do

        res = bootstrapAccumulator
        res_sq = bootstrapAccumulator_sq
        bootstrapAccumulator = ZERO
        bootstrapAccumulator_sq = ZERO

        bootstrapMean(i) = bootstrapMean(i) & 
        +  res_sq * inv_N * inv_Nm1 - (res / N) * (res / N) * inv_Nm1

      end do
      !$omp end parallel do

      mean = bootstrapMean(i) / self % nBootstraps
      mean_biasAdj = TWO * plugInVar - mean

      self % bootstrapMean(loc + i) = mean        !biased
      self % bootstrapVar(loc + i) = mean_biasAdj !unbiased

    end do
    !$omp end parallel do

    deallocate(bootstrapMean)
  end subroutine bootstrapV2

 !Bootstrap maximum likelihood var
  subroutine bootstrapV3(self, rand)
    class(scoreMemory), intent(inout)        :: self
    class(RNG), intent(inout)                :: rand
    integer(shortInt)                        :: i, N, j, r, B
    integer(longInt)                         :: loc
    real(defReal)                            :: inv_Bm1, inv_N
    real(defReal), save                      :: res, res_sq, plugInVar
    integer(shortInt), save                  :: idx
    real(defReal), save                      :: bootstrapAccumulator, bootstrapAccumulator_sq
    real(defReal), dimension(:), allocatable, save :: bootstrapMean
    real(defReal), save                      :: mean, mean_biasAdj

    N = self % cyclesPerTime
    B = self % nBootstraps

    if (B /= 1) then
      inv_Bm1 = ONE / (B - 1)
    else
      inv_Bm1 = ONE
    end if

    inv_N   = ONE / N

    allocate(bootstrapMean(self % Ntallies))
    bootstrapAccumulator = ZERO
    bootstrapAccumulator_sq = ZERO
    bootstrapMean = ZERO

    loc = (self % timeN - 1_longInt) * self % Ntallies


    !$omp parallel do 
    do i=1, size(self % plugInSamples(1,:))
      write(10, '(F24.16, ",")') self % plugInSamples(1,i)
    end do
    !$omp end parallel do 

    !$omp parallel do private(mean, mean_biasAdj, bootstrapAccumulator, bootstrapAccumulator_sq, plugInVar, bootstrapMean)
    do i = 1, self % Ntallies

      plugInVar = self % bins(loc + i, CSUM2) * inv_N * inv_N &
      - (self % bins(loc + i, CSUM) / N) * (self % bins(loc + i, CSUM) / N) * inv_N

      bootstrapMean(i) = plugInVar

      !$omp parallel do private(res, res_sq, bootstrapAccumulator, bootstrapAccumulator_sq)
      do r = 1, self % nBootstraps - 1
        !$omp parallel do private(idx)
        do j = 1, N
          call rand % stride(j)
          validSample: do
            idx = int(N * rand % get()) + 1
            if (idx <= N) then
              bootstrapAccumulator = bootstrapAccumulator + self % plugInSamples(i,idx)
              bootstrapAccumulator_sq = bootstrapAccumulator_sq + self % plugInSamples(i,idx) * self % plugInSamples(i,idx)
              exit validSample
            else
              cycle validSample
            end if

          end do validSample
        end do
        !$omp end parallel do

        res = bootstrapAccumulator
        res_sq = bootstrapAccumulator_sq
        bootstrapAccumulator = ZERO
        bootstrapAccumulator_sq = ZERO

        bootstrapMean(i) = bootstrapMean(i) & 
        +  res_sq * inv_N * inv_N - (res / N) * (res / N) * inv_N

      end do
      !$omp end parallel do

      mean = bootstrapMean(i) / self % nBootstraps
      mean_biasAdj = TWO * plugInVar - mean

      self % bootstrapMean(loc + i) = mean        !biased
      self % bootstrapVar(loc + i) = mean_biasAdj !unbiased

    end do
    !$omp end parallel do

    deallocate(bootstrapMean)
  end subroutine bootstrapV3

end module scoreMemory_class