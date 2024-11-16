program lambda_phi4_2d

  use parameters
  use field
  use pbc
  use dynamics
  use observables
  implicit none

  integer :: i, save_unit, term_unit, obs_unit
  real(dp) :: msq(1), msqi = -10.0_dp, msqf = 10.0_dp
  print"(a)", "2 dimensional lambda phi^4 theory"

  msq = [1.0_dp]![(msqi + ((msqf - msqi)/(size(msq) - 1))*i, i = 0, size(msq) - 1 )]
  
  call read_input()
  call allocate_pbc(Lx)
  allocate(magnetization(Nmeasurements))
  allocate(phi(Lx,Lx))

  call hot_start(phi)
  !call cold_start(phi)
  
  open(newunit = save_unit, file = "data/configurations.dat")
  open(newunit = term_unit,  file = "data/cold_data.dat")
  open(newunit = obs_unit,  file = "data/data.dat")

  do i = 1, size(msq)
     call save_configuration(phi,save_unit)
     call thermalization(phi,Nthermalization,algorithm,[msq(i),lambda],save_unit,term_unit)
     call take_measurements(phi, Nmeasurements, Nskip, algorithm, [msq(i),lambda], obs_unit)
     print*, msq(i), sum(abs(magnetization))/size(magnetization)
  end do
  
end program lambda_phi4_2d
