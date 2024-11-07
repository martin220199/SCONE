program scone

  use numPrecision
  use genericProcedures,          only : printStart
  use openmp_func,                only : ompSetNumThreads
  use commandLineUI,              only : getInputFile, clOptionIsPresent, addClOption, getFromCL
  use dictionary_class,           only : dictionary
  use dictParser_func,            only : fileToDict
  use physicsPackage_inter,       only : physicsPackage
  use physicsPackageFactory_func, only : new_physicsPackage
  use vizPhysicsPackage_class,    only : vizPhysicsPackage
  use timer_mod                 , only : registerTimer, timerStart, timerStop, timerTime, secToChar

  implicit none
  type(dictionary)                  :: input
  class(physicsPackage),allocatable :: core
  character(:),allocatable          :: inputPath
  integer(shortInt)                 :: timerIdx
  integer(shortInt)                 :: cores

  integer(shortInt)                 :: k
  character(len=20)                 :: filename_1, filename_2, filename_3, filename_4, filename_5, filename_6
  integer                           :: file_unit_1, file_unit_2, file_unit_3, file_unit_4, file_unit_5, file_unit_6

  ! Add command line options here
  call addClOption('--plot', 0, ['int'],&
          'Executes geometry plotting specified by a viz dict in the input file')
#ifdef _OPENMP
  call addClOption('--omp', 1, ['int'], &
          'Number of OpenMP threads in a parallel calculation')
#endif

  ! Get path to input file
  call getInputFile(inputPath)

  ! Set Number of threads
  if (clOptionIsPresent('--omp')) then
    call getFromCL(cores, '--omp', 1)
  else
    cores = 1
  end if
  call ompSetNumThreads(cores)

  ! Register timer
  timerIdx = registerTimer('Main Timer')

  call printStart()

  call fileToDict(input, inputPath)

  call timerStart(timerIdx)

  ! Define the filename and assign a file unit number
  filename_1 = 'samples.txt'
  file_unit_1 = 10
  open(unit=file_unit_1, file=filename_1, status='replace')

  filename_2 = 'bs_std.txt'
  file_unit_2 = 11
  open(unit=file_unit_2, file=filename_2, status='replace')

  filename_3 = 'mc_std.txt'
  file_unit_3 = 12
  open(unit=file_unit_3, file=filename_3, status='replace')

  filename_4 = 'bs_time.txt'
  file_unit_4 = 13
  open(unit=file_unit_4, file=filename_4, status='replace')

  filename_5 = 'mc_time.txt'
  file_unit_5 = 14
  open(unit=file_unit_5, file=filename_5, status='replace')

  filename_6 = 'bs_mean.txt'
  file_unit_6 = 15
  open(unit=file_unit_6, file=filename_6, status='replace')


  if (clOptionIsPresent('--plot')) then
    allocate(vizPhysicsPackage :: core)
    call core % init(input)
    call core % run()
  
  else

    do k = 1, 1

      print *, '------- k', k

      if (allocated(core)) then
        call core % kill()
        deallocate(core)
      end if

      allocate( core, source = new_physicsPackage(input))


      call core % run()

    end do

  end if


  ! Close the file
  close(file_unit_1)
  close(file_unit_2)
  close(file_unit_3)
  close(file_unit_4)
  close(file_unit_5)
  close(file_unit_6)

  call timerStop(timerIdx)
  print *, 'Total calculation time: ', trim(secToChar(timerTime(timerIdx)))
  print *, 'Have a good day and enjoy your result analysis!'
end program scone
