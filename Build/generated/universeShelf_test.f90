module universeShelf_test

  use numPrecision
  use universalVariables,  only : OUTSIDE_MAT
  use dictionary_class,    only : dictionary
  use dictParser_func,     only : charToDict
  use charMap_class,       only : charMap
  use surfaceShelf_class,  only : surfaceShelf
  use cellShelf_class,     only : cellShelf
  use universe_inter,      only : universe
  use uniFills_class,      only : uniFills
  use universeShelf_class, only : universeShelf
  use pfUnit_mod

  implicit none

  ! Parameters
  character(*), parameter :: SURFS_DEF = " sph {id 1; type sphere; origin (0.0 0.0 0.0); radius 9.0;}"
  character(*), parameter :: CELLS_DEF = " "
  character(*), parameter :: UNIS_DEF = &
  " root {id 1; type rootUniverse; border 1; fill u<10>;} &
  & lat { id 10; type latUniverse; shape (2 2 0); pitch (1.0 1.0 0.0); padMat void; map (20 2 2 2);} &
  & pin { id 2; type pinUniverse; radii (1.0 0.0); fills (fuel water);} &
  & pin2 {id 20; type pinUniverse; radii (1.0 1.5 0.0); fills (u<21> fuel water);} &
  & fuel {id 21; type pinUniverse; radii (999.0 0.0); fills (fuel fuel);}"

  ! Variables
  type(charMap)       :: mats
  type(surfaceShelf)  :: surfs
  type(cellShelf)     :: cells
  type(universeShelf) :: unis

contains

  !!
  !! Build test variables
  !!
!@Before
  subroutine set_up()
    type(uniFills)     :: fills
    type(dictionary)   :: dict
    character(nameLen) :: name
    integer(shortInt)  :: idx, idx2, idx20

    ! Build materials
    name = 'water'
    call mats % add(name, 1)

    name = 'fuel'
    call mats % add(name, 2)

    name = 'void'
    call mats % add(name, 3)

    ! Build surfaces
    call charToDict(dict, SURFS_DEF)
    call surfs % init(dict)
    call dict % kill()

    ! Build cells
    call charToDict(dict, CELLS_DEF)
    call cells % init(dict, surfs, mats)
    call dict % kill()

    ! Build universes
    call charToDict(dict, UNIS_DEF)
    call unis % init(fills, dict, cells, surfs, mats)

    ! Checks fills
    ! Size
#line 71 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
  call assertEqual(5, size(fills % uni), &
 & location=SourceLocation( &
 & 'universeShelf_test.f90', &
 & 71) )
  if (anyExceptions()) return
#line 72 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"

    ! Fills in lattice universe
    idx = unis % getIdx(10)
    idx2 = unis % getIdx(2)
    idx20 = unis % getIdx(20)
#line 77 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
  call assertEqual([-idx2, -idx2, -idx20, -idx2, 3], fills % uni(idx) % fill , &
 & location=SourceLocation( &
 & 'universeShelf_test.f90', &
 & 77) )
  if (anyExceptions()) return
#line 78 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"

  end subroutine set_up

  !!
  !! Clean after test
  !!
!@After
  subroutine clean_up()

    call mats  % kill()
    call surfs % kill()
    call cells % kill()
    call unis  % kill()

  end subroutine clean_up

  !!
  !! Test getting ptr/idx/id
  !!
!@Test
  subroutine test_get()
    integer(shortInt) :: idx
    class(universe), pointer :: ptr

    ! Universe ID 10
    idx = unis % getIdx(10)
    ptr => unis % getPtr(idx)
#line 105 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
  call assertEqual(10, ptr % id(), &
 & location=SourceLocation( &
 & 'universeShelf_test.f90', &
 & 105) )
  if (anyExceptions()) return
#line 106 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
#line 106 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
  call assertEqual(10, unis % getID(idx), &
 & location=SourceLocation( &
 & 'universeShelf_test.f90', &
 & 106) )
  if (anyExceptions()) return
#line 107 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"

    ! Universe ID 21 -> With the fast access 
    idx = unis % getIdx(21)
    ptr => unis % getPtr_fast(idx)
#line 111 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
  call assertEqual(21, ptr % id(), &
 & location=SourceLocation( &
 & 'universeShelf_test.f90', &
 & 111) )
  if (anyExceptions()) return
#line 112 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
#line 112 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
  call assertEqual(21, unis % getID(idx), &
 & location=SourceLocation( &
 & 'universeShelf_test.f90', &
 & 112) )
  if (anyExceptions()) return
#line 113 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"

    ! Universe ID 1
    idx = unis % getIdx(1)
    ptr => unis % getPtr(idx)
#line 117 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
  call assertEqual(1, ptr % id(), &
 & location=SourceLocation( &
 & 'universeShelf_test.f90', &
 & 117) )
  if (anyExceptions()) return
#line 118 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
#line 118 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
  call assertEqual(1, unis % getID(idx), &
 & location=SourceLocation( &
 & 'universeShelf_test.f90', &
 & 118) )
  if (anyExceptions()) return
#line 119 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"

    ! Test size
#line 121 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"
  call assertEqual(5, unis % getSize(), &
 & location=SourceLocation( &
 & 'universeShelf_test.f90', &
 & 121) )
  if (anyExceptions()) return
#line 122 "/home/mskrette/SCONE_cambridge_fork/SCONE/Geometry/Universes/Tests/universeShelf_test.f90"

  end subroutine test_get


end module universeShelf_test

module WrapuniverseShelf_test
   use pFUnit_mod
   use universeShelf_test
   implicit none
   private

contains


end module WrapuniverseShelf_test

function universeShelf_test_suite() result(suite)
   use pFUnit_mod
   use universeShelf_test
   use WrapuniverseShelf_test
   type (TestSuite) :: suite

   suite = newTestSuite('universeShelf_test_suite')

   call suite%addTest(newTestMethod('test_get', test_get, set_up, clean_up))


end function universeShelf_test_suite

