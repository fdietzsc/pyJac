module pyjac_mod
  use iso_c_binding

  interface
    subroutine fcinterface_create(T,h) bind(c,name='eval_h')
      use iso_c_binding
      real(kind=c_double), intent(in   ) :: T
      real(kind=c_double), intent(inout) :: h(*)
    end subroutine fcinterface_create
  end interface

end module pyjac_mod
