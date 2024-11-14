module parameters

  use precision
  implicit none

  integer(i4) :: Lx, Ly
  real(dp) :: lambda
  real(dp) :: epsilon

  namelist /lattice/ Lx, Ly
  namelist /pars/ lambda, epsilon
  
contains

  subroutine read_input()

    character(100) :: parameters_file
    integer(i4) :: inunit
    write(*,"(a)") "Please enter the parameters file:"
    read(*,"(a)") parameters_file
    write(*,"(2a)") "User typed: ", trim(parameters_file)
    open(newunit = inunit, file = trim(parameters_file))

    read(inunit, nml = lattice)
    read(inunit, nml = pars)

    close(inunit)

    write(*, nml = lattice)
    write(*, nml = pars)
    
  end subroutine read_input

    
end module parameters
