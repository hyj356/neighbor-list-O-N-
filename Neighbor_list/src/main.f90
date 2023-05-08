program main
  use global_var, only:WP, stdout
  use fileIO, only: ReadData, Model, ReadEamAlloyFile, eamAlloyFile
  use neighbor, only: find_neighbor, BuildNeighbors
  use computePE, only: Energy
  implicit none
  type(Model) :: atoms
  type(eamAlloyFile) :: potential 
  real(WP), allocatable :: f(:, :), pe(:)
  real(WP) :: start, ends
  integer :: i

  call ReadData("../model/CuTa_huge.lmp", atoms)
  call ReadEamAlloyFile("../model/CuTa.eam.alloy", potential)
  write(stdout, *) "lammps result: -173868.420235853 ev"

  call cpu_time(start)
  call find_neighbor(atoms, potential%cutoff)
  call Energy(potential, atoms, f, pe)
  call cpu_time(ends)
  write(stdout, *) sum(pe), 'ev'
  write(stdout, *) "Cell-list method cost ", ends - start, 's'

  call cpu_time(start)
  call BuildNeighbors(atoms, potential%cutoff2)
  call Energy(potential, atoms, f, pe)
  call cpu_time(ends)
  write(stdout, *) sum(pe), 'ev'
  write(stdout, *) "Verlet-list method cost ", ends - start, 's'
  
  open(999, file='tmp.txt', action='write')
  do i = 1, atoms%nAtoms
    write(999, '(*(F0.6, 1x))') f(:, i), pe(i)
  end do
  ! write(stdout, *) sum(pe), 'ev'
  ! write(stdout, *) sum(f(1, :)), 'ev/A'
  ! write(stdout, *) sum(f(2, :)), 'ev/A'
  ! write(stdout, *) sum(f(3, :)), 'ev/A'
  ! Cu50Ta50: -173868.420235853 ev
  
end program main