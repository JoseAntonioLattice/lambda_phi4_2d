module dynamics

  use precision
  implicit none

contains

  subroutine save_configuration(phi, save_unit)
    real(dp), intent(in) :: phi(:,:)
    integer(i4), intent(in) :: save_unit
    integer(i4) :: Lx, Ly, i

    Lx = size(phi(:,1))
    Ly = size(phi(1,:))
    do i = 1, Lx
       write(save_unit,"("//int2str(Ly)//"(f10.4,2x))") phi(i,:)
    end do
    write(save_unit,*) " "
  end subroutine save_configuration

  function int2str(i) result(k)
    integer(i4), intent(in) :: i
    character(:), allocatable :: k
    character(20) :: l

    write(l,*) i
    k = trim(l)
    
  end function int2str
  
  subroutine hot_start(phi)
    use parameters, only : epsilon
    real(dp), intent(out) :: phi(:,:)

    call random_number(phi)

    phi = epsilon * (2*phi - 1.0_dp)

  end subroutine hot_start

  subroutine cold_start(phi)
    real(dp), intent(out) :: phi(:,:)

    phi = 0.0_dp

  end subroutine cold_start

  
  subroutine thermalization(phi,Nthermalization,algorithm,args,save_unit,obs_unit)
    use functions, only : action
    real(dp), intent(inout) :: phi(:,:)
    integer(i4), intent(in) :: Nthermalization
    character(*), intent(in) :: algorithm
    real(dp), intent(in) :: args(2)
    integer(i4), intent(in) :: save_unit, obs_unit
    integer(i4) :: i,Lx, Ly, Vol

    Lx = size(phi(:,1))
    Ly = size(phi(1,:))
    Vol = Lx*Ly
    
    do i = 1, Nthermalization
       call sweeps(phi,trim(algorithm),args)
       call save_configuration(phi, save_unit)
       write(obs_unit,*) action(phi,args(1),args(2))/vol, sum(phi)/vol
    end do
    
  end subroutine thermalization

  subroutine take_measurements(phi, Nmeasurements, Nskip, algorithm, args, obs_unit, conf_unit)
    use functions, only : action
    use observables
    real(dp), intent(inout) :: phi(:,:)
    integer(i4), intent(in) :: Nmeasurements, Nskip
    character(*), intent(in) :: algorithm
    real(dp), intent(in) :: args(2)
    integer(i4), intent(in) :: obs_unit, conf_unit
    integer(i4) :: i, j, k, Lx, Ly, Vol

    Lx = size(phi(:,1))
    Ly = size(phi(1,:))
    Vol = Lx*Ly
    
    do i = 1, Nmeasurements
       do j = 1, Nskip
          call sweeps(phi,trim(algorithm),args)
       end do
       magnetization(i) = sum(phi)
       mag1(i) = sum(phi(1,:))/Lx
       do k = 1, Lx
          correlation(i,k) = mag1(i)*sum(phi(k,:))/Lx**2
       end do
       write(obs_unit,*) action(phi,args(1),args(2))/Vol, magnetization(i)/Vol, abs(magnetization(i))/Vol, mag1(i), correlation(i,:)
       call save_configuration(phi, conf_unit)
    end do
    
  end subroutine take_measurements

  subroutine sweeps(phi,algorithm,args)
    real(dp), intent(inout) :: phi(:,:)
    character(*), intent(in) :: algorithm
    real(dp), dimension(2), intent(in) :: args
    integer(i4) :: Lx, Ly, i, j
    real(dp) :: acceptance_rate, msq, lambda
    Lx = size(phi(:,1))
    Ly = size(phi(1,:))

    msq = args(1)
    lambda = args(2)
    do i = 1, Lx
       do j = 1, Ly
          call metropolis(phi,[i,j],acceptance_rate,msq,lambda)
       end do
    end do
    
  end subroutine sweeps

  subroutine metropolis(phi,x,acceptance_rate,msq,lambda) 
    use functions, only : DS, DS2
    use parameters, only : epsilon
    real(dp), intent(inout) :: phi(:,:)
    integer(i4), intent(in) :: x(2)
    real(dp), intent(out) :: acceptance_rate
    real(dp), intent(in) :: msq, lambda
    real(dp) :: r,phi_p

    call random_number(phi_p)
    phi_p = phi(x(1),x(2)) + epsilon * (2*phi_p - 1.0_dp)
    
    call random_number(r)
    acceptance_rate = min(1.0_dp, exp(-DS(phi,phi_p,x,msq,lambda)))

    if (r <= acceptance_rate) phi(x(1),x(2)) = phi_p
  end subroutine metropolis
 

  
end module dynamics
