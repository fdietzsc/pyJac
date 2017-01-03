module pyjac_mod
  use iso_c_binding
  interface
    subroutine compute_hi(T,h) bind(c,name='eval_h')
      use iso_c_binding
      ! value keyword important, because usually
      ! fortran passes by reference, but pyjac
      ! requires value to be passed
      real(kind=c_double), intent(in   ), value :: T
      real(kind=c_double) :: h(*)
    end subroutine compute_hi
  end interface

  interface
    subroutine compute_cpi(T,cp) bind(c,name='eval_cp')
      use iso_c_binding
      ! value keyword important, because usually
      ! fortran passes by reference, but pyjac
      ! requires value to be passed
      real(kind=c_double), intent(in   ), value :: T
      real(kind=c_double) :: cp(*)
    end subroutine compute_cpi
  end interface

  interface
    subroutine convert_x2y(X,Y) bind(c,name='mole2mass')
      use iso_c_binding
      ! value keyword important, because usually
      ! fortran passes by reference, but pyjac
      ! requires value to be passed
      real(kind=c_double), intent(in   ) :: X(*)
      real(kind=c_double) :: Y(*)
    end subroutine convert_x2y
  end interface

  interface
    subroutine convert_y2x(Y,X) bind(c,name='mass2mole')
      use iso_c_binding
      ! value keyword important, because usually
      ! fortran passes by reference, but pyjac
      ! requires value to be passed
      real(kind=c_double), intent(in   ) :: Y(*)
      real(kind=c_double) :: X(*)
    end subroutine convert_y2x
  end interface

  interface
    real(c_double) function computeDensity(T,P,X) bind(c,name='getDensity')
      use iso_c_binding
      ! value keyword important, because usually
      ! fortran passes by reference, but pyjac
      ! requires value to be passed
      real(kind=c_double), intent(in   ), value :: T
      real(kind=c_double), intent(in   ), value :: P
      real(kind=c_double), intent(in   ) :: X(*)
    end function computeDensity
  end interface

  interface
    subroutine reaction_rates(T,P,Y,omega) bind(c,name='dydt')
      use iso_c_binding
      ! value keyword important, because usually
      ! fortran passes by reference, but pyjac
      ! requires value to be passed
      real(kind=c_double), intent(in   ), value :: T
      real(kind=c_double), intent(in   ), value :: P
      real(kind=c_double), intent(in   ) :: Y(*)
      ! remember omega(1) is source term for temperature
      real(kind=c_double) :: omega(*)
    end subroutine reaction_rates
  end interface

  interface
    subroutine jacobian(T,P,Y,jac) bind(c,name='eval_jacob')
      use iso_c_binding
      ! value keyword important, because usually
      ! fortran passes by reference, but pyjac
      ! requires value to be passed
      real(kind=c_double), intent(in   ), value :: T
      real(kind=c_double), intent(in   ), value :: P
      real(kind=c_double), intent(in   ) :: Y(*)
      ! remember omega(1) is source term for temperature
      real(kind=c_double) :: jac(*)
    end subroutine jacobian
  end interface

end module pyjac_mod
