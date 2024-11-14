module functions

  use precision
  implicit none

contains

  function lagrangian(phi,x,msq,lambda)
    real(dp), intent(in) :: phi(:,:)
    integer(i4), intent(in) :: x(2)
    real(dp), intent(in) :: msq, lambda

    real(dp) :: lagrangian

    lagrangian = 0.5 * ( ( phi(ip(x(1)),x(2)) - phi(x(1),x(2)) ) ** 2  &
                       + ( phi(x(1),ip(x(2))) - phi(x(1),x(2)) ) ** 2 )&
               + 0.5 * msq * phi(x(1),x(2)) ** 2 &
               + 0.25 * lambda * phi(x(1),x(2)) ** 4 
    
  end function lagrangian
  
  function action(phi,msq,lambda)
    real(dp), intent(in) :: phi(:,:)
    real(dp), intent(in) :: msq, lambda

    real(dp) :: action

    action  =  0.0_dp

    do x1 = 1, size(phi(:,1))
       do x2 = 1, size(phi(1,:))
          action = action + lagrangian(phi,[x1,x2],msq,lambda)
       end do
    end do

  end function action

  function DS(phi,msq,lambda)
    real(dp), intent(in) :: phi(:,:)
    real(dp), intent(in) :: msq, lambda
    real(dp) :: DS, phi_p, r

    call random_number(r)
    phi_p = r
    DS = 2 * (phi_p**2 - phi(x(1),x(2))**2) &
    + (phi_p - phi(x(1),x(2))) &
    * ( phi(ip(x(1)),x(2)) + phi(x(1),ip(x(2))) + phi(im(x(1)),x(2)) + phi(x(1),im(x(2))) )&
    + 0.5 * msq * ( phi_p**2 - phi(x(1),x(2))**2 ) &
    + 0.25 * lambda * ( phi_p**4 - phi(x(1),x(2))**4 )
    
  end function DS
  
end module functions
