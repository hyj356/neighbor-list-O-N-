module spline
  use global_var, only: WP, STDOUT
  implicit none

  private

  public :: simple_spline_t
  type simple_spline_t

     integer            :: n = -1   !< number of entries
     real(WP)           :: x0       !< x数组的起点
     real(WP)           :: cut      !< end point
     real(WP)           :: dx       !< 相邻x的差值

     logical            :: associated = .false.

     real(WP), pointer  :: y(:)        => NULL()   ! value tables
     real(WP), pointer  :: d2y(:)      => NULL()   ! second derivatives table

     real(WP), pointer  :: coeff1(:)   => NULL()  ! spline coefficients
     real(WP), pointer  :: coeff2(:)   => NULL()  ! spline coefficients
     real(WP), pointer  :: coeff3(:)   => NULL()  ! spline coefficients
     real(WP), pointer  :: dcoeff1(:)  => NULL()  ! spline coefficients
     real(WP), pointer  :: dcoeff2(:)  => NULL()  ! spline coefficients
     real(WP), pointer  :: dcoeff3(:)  => NULL()  ! spline coefficients

  endtype simple_spline_t

  public :: init
  interface init
     module procedure simple_spline_init
  endinterface

  public :: del
  interface del
     module procedure simple_spline_del
  endinterface

  public :: associate
  interface associate
     module procedure simple_spline_associate
  endinterface

  public :: func
  interface func
     module procedure simple_spline_f, simple_spline_f_array
  endinterface

  public :: dfunc
  interface dfunc
     module procedure simple_spline_df
  endinterface

  public :: f_and_df
  interface f_and_df
     module procedure simple_spline_f_and_df
  endinterface

contains

  !>
  !! 
  !!  根据传入参数初始化Spline线条的系数表格
  !!
  !<
  pure subroutine simple_spline_init(this, n, x0, dx, y)
    implicit none

    type(simple_spline_t), intent(out)  :: this
    integer, intent(in)                 :: n
    real(WP), intent(in)                :: x0
    real(WP), intent(in)                :: dx
    real(WP), intent(in)                :: y(n)

    ! ---

    real(WP), parameter  :: sig = 0.5_WP

    integer  :: i, k
    real(WP) :: p, qn, un, u(n)

    ! ---

    this%associated  = .false.

    this%n     = n
    this%x0    = x0
    this%dx    = dx

    allocate(this%y(n))
    allocate(this%d2y(n))
    allocate(this%coeff1(n-1))
    allocate(this%coeff2(n-1))
    allocate(this%coeff3(n-1))
    allocate(this%dcoeff1(n-1))
    allocate(this%dcoeff2(n-1))
    allocate(this%dcoeff3(n-1))

    this%y(:)   = y(:)

    this%d2y(1) = 0.0_WP ! natural simple_spline, i.e. second derivatives vanishes
    u(1)        = 0.0_WP

    do i = 2, n-1
      p        = sig*this%d2y(i-1) + 2
      this%d2y(i) = (sig-1)/p             !-p/2
      u(i)     = (6.0_WP*((this%y(i+1)-this%y(i))/dx &
            -(this%y(i)-this%y(i-1))/dx)/(2*dx)-sig*u(i-1))/p
    enddo

    qn = 0.0_WP ! 自然边界条件
    un = 0.0_WP
    !    qn  = 0.5_WP
    !    un  = 3.0_WP*(0.0_WP-(this%y(n)-this%y(n-1)))/(dx**2)

    this%d2y(n) = (un-qn*u(n-1))/(qn*this%d2y(n-1)+1.0_WP)

    do k = n-1, 1, -1
      this%d2y(k) = this%d2y(k)*this%d2y(k+1) + u(k)
    enddo

    do k = 1, n-1
      this%coeff1(k)   = this%y(k+1) - this%y(k) - ( 2*this%d2y(k) + this%d2y(k+1) )*this%dx**2/6
      this%coeff2(k)   = this%d2y(k) * this%dx**2/2
      this%coeff3(k)   = ( this%d2y(k+1) - this%d2y(k) )*this%dx**2/6
    enddo

    this%dcoeff1  = this%coeff1/this%dx
    this%dcoeff2  = 2*this%coeff2/this%dx
    this%dcoeff3  = 3*this%coeff3/this%dx

    this%cut = this%x0+dx*(n-1)

  endsubroutine simple_spline_init

  !>
  !! 
  !!  删除一个spline的系数表格
  !!
  !<
  elemental subroutine simple_spline_del(this)
    implicit none

    type(simple_spline_t), intent(inout)  :: this

    ! ---

    if (.not. this%associated) then
       if (associated(this%y))        deallocate(this%y)
       if (associated(this%d2y))      deallocate(this%d2y)
       if (associated(this%coeff1))   deallocate(this%coeff1)
       if (associated(this%coeff2))   deallocate(this%coeff2)
       if (associated(this%coeff3))   deallocate(this%coeff3)
       if (associated(this%dcoeff1))  deallocate(this%dcoeff1)
       if (associated(this%dcoeff2))  deallocate(this%dcoeff2)
       if (associated(this%dcoeff3))  deallocate(this%dcoeff3)

       this%y       => NULL()
       this%d2y     => NULL()
       this%coeff1  => NULL()
       this%coeff2  => NULL()
       this%coeff3  => NULL()
       this%dcoeff1 => NULL()
       this%dcoeff2 => NULL()
       this%dcoeff3 => NULL()
    endif

  end subroutine simple_spline_del

  !>
  !!  将传入的this的spline指向另一个spline类型的变量w
  !!  这个操作类似于python中的浅拷贝, 仅仅只是改变指针
  !!  的指向, 并不是将数据直接拷贝过来
  !! 
  !<
  subroutine simple_spline_associate(this, w)
    implicit none

    type(simple_spline_t), intent(inout) :: this  !< 准备链接的spline变量
    type(simple_spline_t), intent(in)    :: w     !< 被链接的spline变量

    ! --

    this%associated = .true.

    this%n        = w%n
    this%x0       = w%x0
    this%cut      = w%cut
    this%dx       = w%dx

    this%y        => w%y
    this%d2y      => w%d2y

    this%coeff1   => w%coeff1
    this%coeff2   => w%coeff2
    this%coeff3   => w%coeff3

    this%dcoeff1  => w%dcoeff1
    this%dcoeff2  => w%dcoeff2
    this%dcoeff3  => w%dcoeff3

  end subroutine simple_spline_associate

  !>
  !! 根据传入的spline系数表插值计算函数的值
  !<
  function simple_spline_f(this, x, ierror)
    implicit none

    type(simple_spline_t), intent(in)  :: this
    real(WP), intent(in)               :: x     
    real(WP)                           :: simple_spline_f

    ! ---
    real(WP)  :: dx, B, xf
    integer   :: i
    integer,intent(inout) :: ierror
    ! ---

    dx  = this%dx

    if (x == this%cut) then
       xf        = this%n
       i         = this%n-1
    else
       xf        = (x-this%x0)/dx+1
       i         = floor(xf)
    endif

    if (i < 1 .or. i >= this%n) then
       simple_spline_f  = 1.0_WP
       ierror = 1       ! 如果发现待插的数值x超出范围, 给出一个不为0的值提升插值错误
       !write(STDOUT, *) "x = " // x // " outside of the defined interval."
       Return
    endif

    B                = xf-i
    simple_spline_f  = this%y(i) + B*(this%coeff1(i) + B*(this%coeff2(i) + B*this%coeff3(i)))
    ierror = 0          ! 如果通过前面的检测插值成功, 那么返回一个为0的值
    
  end function simple_spline_f

  !>
  !! 
  !!  作用大体类似于上述的函数, 但是传入的x是一个一维的参数
  !! 
  !<
  function simple_spline_f_array(this, x, ierror)
    implicit none

    type(simple_spline_t), intent(in)  :: this
    real(WP), intent(in)               :: x(:)
    real(WP)                           :: simple_spline_f_array(lbound(x, 1):ubound(x, 1))
    integer, intent(inout)             :: ierror(:)

    ! ---

    integer  :: i

    ! ---

    do i = lbound(x, 1), ubound(x, 1)
       simple_spline_f_array(i)  = simple_spline_f(this, x(i), ierror(i))
    enddo

  end function simple_spline_f_array

  !>
  !! 
  !!  返回函数的在插值点的一阶导数
  !! 
  !<
  function simple_spline_df(this, x, ierror)
    implicit none

    type(simple_spline_t), intent(in)  :: this
    real(WP), intent(in)               :: x
    real(WP)                           :: simple_spline_df
    integer, intent(inout), optional   :: ierror

    ! ---
    real(WP)  :: dx, B, xf
    integer   :: i

    ! ---

    dx  = this%dx

    if (x == this%cut) then
       xf        = this%n
       i         = this%n-1
    else
       xf        = (x-this%x0)/dx+1
       i         = floor(xf)
    endif

    if (i < 1 .or. i >= this%n) then
       simple_spline_df  = 1.0_WP
       ierror = 1       ! 如果发现待插的数值x超出范围, 给出一个不为0的值提升插值错误
       !RAISE_ERROR("x = " // x // " outside of the defined interval.", ierror)
    endif

    B                 = xf-i
    simple_spline_df  = this%dcoeff1(i) + B*(this%dcoeff2(i) + B*this%dcoeff3(i))
    ierror = 0          ! 如果通过前面的检测插值成功, 那么返回一个为0的值

  endfunction simple_spline_df

  !>
  !! 
  !!  返回函数的在插值点的函数值和一阶导数
  !! 
  !<
  subroutine simple_spline_f_and_df(this, x, f, df, extrapolate, ierror)
    implicit none

    type(simple_spline_t), intent(in)    :: this
    real(WP),              intent(in)    :: x
    real(WP),              intent(out)   :: f
    real(WP),              intent(out)   :: df
    logical,     optional, intent(in)    :: extrapolate
    integer,     optional, intent(inout) :: ierror

    ! ---
    real(WP)  :: dx, B, xf
    integer   :: i

    ! ---

    dx = this%dx

    if (present(extrapolate) .and. extrapolate) then
       xf = (x-this%x0)/dx+1
       i  = floor(xf)

       if (i < 1) then
          i = 1
       else if (i >= this%n) then
          i = this%n-1
       endif
    else
       if (x == this%cut) then
          xf = this%n
          i  = this%n-1
       else
          xf = (x-this%x0)/dx+1
          i  = floor(xf)
       endif

       if (i < 1 .or. i >= this%n) then
          f  = 1.0_WP
          df = 1.0_WP
          ierror = 1       ! 如果发现待插的数值x超出范围, 给出一个不为0的值提示插值错误
          !RAISE_ERROR("x = " // x // " outside of the defined interval.", ierror)
       endif
    endif

    B   = xf-i
    f   = this%y(i) + B*(this%coeff1(i) + B*(this%coeff2(i) + B*this%coeff3(i)))
    df  = this%dcoeff1(i) + B*(this%dcoeff2(i) + B*this%dcoeff3(i))
    ierror = 0          ! 如果通过前面的检测插值成功, 那么返回一个为0的值

  endsubroutine simple_spline_f_and_df

end module spline