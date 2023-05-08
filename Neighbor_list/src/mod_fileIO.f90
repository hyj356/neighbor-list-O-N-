module fileIO
  use global_var, only: wp, stdout, Bohr, Hartree, si, li
  use spline, only: simple_spline_t, init, associate
  implicit none
  private
  public :: ReadData, Model, &
  AllocMem, eamAlloyFile, ReadEamAlloyFile, WriteData,    &
  WriteLog

  type eamAlloyFile
    !! 一个用于存储eam/alloy类势函数文件中所有有效信息的类, 使用type进行管理,方便调用
    character(len=:), allocatable :: filename   !< 势函数文件名称
    real(wp), allocatable :: F(:, :)            !< 多体势, 也称为嵌入势的数组
    real(wp), allocatable :: Z(:, :, :)         !< 两体势数组(pair interaction)
    real(wp), allocatable :: Rho(:, :)          !< 电子密度, 与距离有关, 用于计算多体势
    ! real(wp), allocatable :: rx(:)              !< 原子距离在截断半径上的等差数组
    ! real(wp), allocatable :: rhox(:)            !< 电荷密度的等差数组
    real(wp), allocatable :: Amass(:)           !< 所有元素的原子质量
    real(wp), allocatable :: lattice(:)         !< 所有元素的晶格常数
    integer, allocatable  :: AtomicNumber(:)    !< 所有元素的原子序数
    character(len=3), allocatable :: latType(:)     !< 所有元素的晶格类型, 对于金属来说, 一般只有FCC, BCC, HCP三种晶格类型
    character(len=2), allocatable :: elementName(:) !< 所有元素的符号名称
    real(wp) :: dr                 !< r的取点间隔
    real(wp) :: dri                !< r的取点间隔的导数, 避免反复的除法运算
    real(wp) :: drho               !< 电子密度rho的取点间隔
    real(wp) :: drhoi              !< 电子密度rho的取点间隔
    real(wp) :: cutoff             !< eam势函数的截断半径
    real(wp) :: cutoff2            !< eam势函数的截断半径的平方
    integer :: atomTypes           !< eam/alloy文件中一共有几种金属元素
    integer  :: Nr                 !< r的取点个数
    integer  :: Nrho               !< 电子密度rho的的数据点个数
    type(simple_spline_t), allocatable  :: fF(:)      !< Embedding function
    type(simple_spline_t), allocatable  :: fphi(:, :) !< phi(r) - repulsive part
    type(simple_spline_t), allocatable  :: frho(:)    !< rho(r) - embedding density
  end type eamAlloyFile

  type Neighbors
    integer, allocatable :: ID(:)     !< 对应每个原子在截断半径以内的ID
    integer :: neb_n = 0              !< 对应原子的邻居原子个数
  end type

  type Model
    integer, allocatable :: typeId(:) !< 原子的类型ID
    real(wp), allocatable :: r(:, :)  !< 原子的坐标
    real(wp) :: xlo, xhi              !< 原子在x方向上的上下限
    real(wp) :: ylo, yhi              !< 原子在y方向上的上下限
    real(wp) :: zlo, zhi              !< 原子在z方向上的上下限
    real(wp) :: lx, ly, lz            !< 原子在三个方向上的长度
    integer :: nAtoms                 !< 模型一共有几个原子
    integer :: nTypes                 !< 模型一共有几种原子
    type(Neighbors), allocatable :: IDs(:)    !< 对应每个原子的邻居
  contains 
    procedure, pass(self), public :: InitialList
    procedure, pass(self), public :: DeleteList
  end type

  interface AllocMem    !! 对AllocMem进行重载, 以使其可以针对不同数据类型不同维度的数组分配内存并在出错时给出合理报错
    module procedure :: AllocMem_1d_dp, AllocMem_2d_dp, AllocMem_3d_dp
    module procedure :: AllocMem_1d_si, AllocMem_2d_si, AllocMem_3d_si
    module procedure :: AllocMem_1d_li, AllocMem_2d_li, AllocMem_3d_li
  end interface

contains

  subroutine ReadEamAlloyFile(fileName, eamPotential)
    !! 将eam/alloy文件读取到eampotential这个派生变量里面
    character(len=*), intent(in) :: filename
    type(eamAlloyFile), intent(out) :: eamPotential  !< 势函数文件的数据集合
    integer :: fileId             !< 势函数文件对应的通道ID
    integer :: io_Stat            !< 检测文件读取是否成功的整数
    logical :: isExist            !< 检测文件是否存在逻辑变量
    integer :: i, j, k            !< 循环变量

    !! 查询文件是否存在, 如果不存在就报错并终止程序
    inquire(file=filename, exist=isExist)
    if (.not. isExist) then
      write(stdout, '(A, A)') filename, " can't be found in this path."
      stop
    end if
    eamPotential%filename = filename
    open(newunit=fileId, file=filename, action='read')

    !! 跳过前三行的注释行
    do i = 1, 3
      read(fileId, *)
    end do
    read(fileId, '(I5)', advance='no', iostat=io_Stat) eamPotential%atomTypes   !! 读取势函数里面有几种元素
    if (io_Stat /= 0) then
      write(stdout, '(A)') 'Error occuring when reading the number of elements.'
      stop
    end if

    !! 根据读取到的数据分配内存
    call AllocMem(eamPotential%Amass, eamPotential%atomTypes)   !! 所有元素的原子质量
    call AllocMem(eamPotential%lattice, eamPotential%atomTypes) !! 所有元素的晶格常数
    allocate(eamPotential%latType(eamPotential%atomTypes))      !! 所有元素的晶格类型
    allocate(eamPotential%elementName(eamPotential%atomTypes))  !! 所有元素的元素符号名称
    call AllocMem(eamPotential%AtomicNumber, eamPotential%atomTypes)  !! 所有元素的原子序号

    !! 读取所有元素的元素符号
    read(fileId, *, iostat=io_Stat) eamPotential%elementName
    if (io_Stat /= 0) then
      write(stdout, '(A)') 'Error occuring when reading the name of elements.'
      stop
    end if

    !! 第五行包含了势函数中非常重要的一些信息, 具体含义在type定义时已经说明
    read(fileId, *, iostat=io_Stat) eamPotential%Nrho, eamPotential%drho, eamPotential%Nr, eamPotential%dr, eamPotential%cutoff
    if (io_Stat /= 0) then
      write(stdout, '(A)') 'Error occuring when reading 5th line of eam/alloy file.'
      stop
    end if

    !! 根据读取到的数据开始分配内存
    call AllocMem(eamPotential%F, eamPotential%Nrho, eamPotential%atomTypes)    ! 嵌入能数组, 即F(ρ)
    call AllocMem(eamPotential%Z, eamPotential%Nr, eamPotential%atomTypes, eamPotential%atomTypes) ! 对势数组, 即φ(r)
    call AllocMem(eamPotential%Rho, eamPotential%Nr, eamPotential%atomTypes)    ! 电子密度数组, 即ρ(r)
    ! call InitialArray(eamPotential%Rhox, eamPotential%rx, &
    !   eamPotential%dRho, eamPotential%dR, &
    !   eamPotential%Nrho, eamPotential%Nr)   !! 原子之间距离的等差数组和电子密度的等差数组
    allocate(eamPotential%fF(eamPotential%atomTypes))
    allocate(eamPotential%frho(eamPotential%atomTypes))
    allocate(eamPotential%fphi(eamPotential%atomTypes, eamPotential%atomTypes))

    !! 接下来开始正式读取嵌入能函数F(ρ)和电子密度函数ρ(r)的具体数值, N种元素共有N个数据块
    do i = 1, eamPotential%atomTypes
      !! 第i种元素的原子序数, 原子质量, 晶格常数, 晶格类型
      read(fileId, *, iostat=io_Stat) eamPotential%AtomicNumber(i), eamPotential%Amass(i), &
                      eamPotential%lattice(i), eamPotential%latType(i)                  
      read(fileId, *, iostat=io_Stat) (eamPotential%F(j, i), j = 1, eamPotential%Nrho)  !! 第i种元素的嵌入能数组
      read(fileId, *, iostat=io_Stat) (eamPotential%Rho(j, i), j = 1, eamPotential%Nr)  !! 第i种元素的电子密度数组
      if (io_Stat /= 0) then
        write(stdout, '(A)') 'Error occuring when reading arrary of embedding energy or electron density.'
        stop
      end if
      ! 根据读取到的嵌入能数组和电子密度数组初始化spline系数
      call init(eamPotential%fF(i), eamPotential%Nrho, 0.0_WP, eamPotential%drho, eamPotential%F(:, i))
      call init(eamPotential%frho(i), eamPotential%Nr, 0.0_WP, eamPotential%dr, eamPotential%Rho(:, i))
    end do

    !! 之后读取对势能数组φ(r), 以CuTa合金为例, 应该有Cu-Cu, Cu-Ta, Ta-Ta对势能, Cu-Ta和Ta的具体数值完全一样
    do i = 1, eamPotential%atomTypes
      do j = 1, i
        !! 这里要注意到Fortran的列优先准则, 就可以提高速度
        read(fileID, *, iostat=io_Stat) (eamPotential%Z(k, i, j), k = 1, eamPotential%Nr)     
        if (io_Stat /= 0) then
          write(stdout, '(A)') 'Error occuring when reading arrary of pair interactions.'
          stop
        end if
        !! 根据读取到的数据计算φ(r)数组的spline插值系数
        call init(eamPotential%fphi(i, j), eamPotential%Nr, 0.0_WP, eamPotential%dr, eamPotential%Z(:, i, j))
        if (i /= j) then
          call associate(eamPotential%fphi(j, i), eamPotential%fphi(i, j))
        endif
      end do
    end do
    close(fileId) !! 文件读取结束

    !! 对读取到的数据做一些处理
    do i = 1, eamPotential%atomTypes
      do j = 1, i
        if (i /= j) then
          do k = 1, eamPotential%Nr
            eamPotential%Z(k, j, i) = eamPotential%Z(k, i, j)   !! 将Z改成对称矩阵, 同样注意列优先准则
          end do
        end if
      end do
    end do

    !! 为减小计算量, 提前处理好一些数据
    eamPotential%dri = 1.d0 / eamPotential%dr
    eamPotential%drhoi = 1.d0 / eamPotential%drho
    eamPotential%cutoff2 = eamPotential%cutoff * eamPotential%cutoff  

  end subroutine ReadEamAlloyFile

  ! subroutine InitialArray(Rhox, rx, dRho, dR, Nrho, Nr)
  !   !! 初始化对于Rho和根据cutoff和dR计算的数学公式
  !   real(wp), allocatable, dimension(:), intent(out) :: Rhox, rx  !< 分别对应电子密度的均分和原子距离的均分数组
  !   real(wp), intent(in) ::dRho, dR         !< 上述数组相邻元素之间的步长
  !   integer, intent(in) :: Nrho, Nr         !< 数组的大小
  !   integer :: i                            !< 循环变量

  !   call AllocMem(Rhox, Nrho)
  !   call AllocMem(rx, Nr)

  !   do concurrent(i = 0 : Nrho - 1)
  !     Rhox(i + 1) = i * dRho
  !   end do

  !   do concurrent(i = 0 : Nrho - 1)
  !     Rx(i + 1) = i * dR
  !   end do

  ! end subroutine InitialArray

  !! 此子程序用于初始化每个原子的邻居列表的内存
  pure subroutine InitialList(self)
    class(Model), intent(inout) :: self
    integer :: i 

    if (.not. allocated(self%IDs)) then
      allocate(self%IDs(self%nAtoms))
      do i = 1, self%nAtoms
        allocate(self%IDs(i)%ID(0))
      end do
    end if

  end subroutine InitialList

  !! 此子程序用于删除原来每个原子的邻居列表的内存
  pure subroutine DeleteList(self)
    class(Model), intent(inout) :: self
    integer :: i 

    if (allocated(self%IDs)) then
      do i = 1, self%nAtoms
        deallocate(self%IDs(i)%ID)
      end do
      deallocate(self%IDs)
    end if

  end subroutine DeleteList

  !! 此子程序用于读取ATOMSK输出的atomic格式的原子模型, 想要让此程序运行
  !! 模型格式必须严格按照ATOMSK输出的.lmp文件的格式
  subroutine ReadData(filename, Atoms)
    character(len=*), intent(in) :: filename      !< 模型文件的名称
    type(Model), intent(out) :: Atoms           !< 存储原子坐标和邻居列表的ATOM类数组
    integer :: fileId             !< 模型文件对应的通道ID
    integer :: i                  !< 循环变量
    logical :: isExist            !< 检测文件是否存在的逻辑变量
    integer :: tempi              !< 用来跳过第一列原子ID的数据
    integer :: jump               !< 根据原子数目确定跳过几行

    inquire(file=filename, exist=isExist)
    if (.not. isExist) then
      write(stdout, '(A, A)') filename, " can't be found in this path."
      stop
    end if
    open(newunit=fileId, file=filename, action='read')

    !! 跳过前2行注释行
    read(fileId, *)
    read(fileId, *)

    !! 读取有多少个原子和多少种元素, 并分配内存
    read(fileId, *) Atoms%nAtoms
    read(fileId, *) Atoms%nTypes
    read(fileId, *)
    allocate(Atoms%r(Atoms%nAtoms, 3))
    allocate(Atoms%typeID(Atoms%nAtoms))

    !! 读取盒子的几个极限
    read(fileId, *) Atoms%xlo, Atoms%xhi
    read(fileId, *) Atoms%ylo, Atoms%yhi
    read(fileId, *) Atoms%zlo, Atoms%zhi
    Atoms%lx = Atoms%xhi - Atoms%xlo
    Atoms%ly = Atoms%yhi - Atoms%ylo
    Atoms%lz = Atoms%zhi - Atoms%zlo
    jump = 6 + Atoms%nTypes

    !! 连续跳过jump行
    do i = 1, jump
      read(fileId, *)
    end do

    !! 读取原子坐标
    do i = 1, Atoms%nAtoms
      read(fileId, *) tempi, Atoms%typeId(i), Atoms%r(i, :)
    end do

    close(fileId) !! 关闭文件

    call Atoms%InitialList()    !! 初始化内存

  end subroutine ReadData

  !! 此子程序仿照ATOMSK格式的输出文件进行输出, 使用此子程序
  !! 写出的模型文件必定可以被上面的ReadData读取
  subroutine WriteData(filename, eampot, Atoms)
    !! 将迭代完成之后的模型以lammps和ovito可以读取的格式输出到文件中
    character(len=*), intent(in) :: filename            !< 模型文件的名称
    type(eamAlloyFile), intent(in) :: eampot            !< 势函数文件名称
    type(model), intent(in) :: Atoms     !< 存储原子坐标和邻居列表的ATOM类数组
    integer :: fileId                   !< 模型文件对应的通道ID
    integer :: i                        !< 循环变量

    !write(stdout, *) nAtoms, nTypes
    open(newunit=fileId, file=filename, action='write')
    write(fileId, '(A)') '# This is the file after MC interation.'
    write(fileId, *)
    write(fileId, '(T10, I0, 1X, A)') Atoms%nAtoms, ' atoms'
    write(fileId, '(T10, I0, 1X, A)') Atoms%nTypes, ' atom types'

    write(fileId, *)
    write(fileId, *) Atoms%xlo, Atoms%xhi, 'xlo xhi'
    write(fileId, *) Atoms%ylo, Atoms%yhi, 'ylo yhi'
    write(fileId, *) Atoms%zlo, Atoms%zhi, 'zlo zhi'
    write(fileId, *)

    write(fileId, '(A)') 'Masses'
    write(fileId, *)

    do i = 1, Atoms%nTypes
      write(fileId, '(T10, I0, 1X, G0, 1X, A, 1x, A)') i, eampot%Amass(i), '# ', eampot%elementName(i)
    end do
    write(fileId, *)

    write(fileId, *) 'Atoms # atomic'
    write(fileId, *)
    do i = 1, Atoms%nAtoms
      write(fileId, *) i, Atoms%typeId, Atoms%r(i, :)
    end do
    close(fileID)
    
  end subroutine WriteData

  subroutine WriteLog(fileID, step, energy, time)
    !! 将计算过程中的数据输出到log文件中
    integer, intent(in) :: fileID
    integer, intent(in) :: step       !< 当前仿真步数
    real(wp), intent(in) :: energy    !< 当前模型势能
    real(wp), intent(in) :: time      !< 计算耗时
  
    write(fileID, *) step, energy, time

  end subroutine WriteLog

  subroutine AllocMem_1d_dp(Array, dim)
    !! 对传入的双精度数组进行分配内存, 如果报错可以给出报错信息
    real(wp), intent(inout), allocatable :: Array(:)  !< 传入的双精度数组
    integer, intent(in) :: dim                        !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(dim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_1d_dp

  subroutine AllocMem_2d_dp(Array, xdim, ydim)
    !! 对传入���二维双精度数组进行分配内存, ��果报错可以给出报错信息
    real(wp), intent(inout), allocatable :: Array(:, :)  !< 传入的双精度数组
    integer, intent(in) :: xdim, ydim                    !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_2d_dp

  subroutine AllocMem_3d_dp(Array, xdim, ydim, zdim)
    !! 对传入的三维双精度数组进行分配内存, 如果报错可以给出报错信息
    real(wp), intent(inout), allocatable :: Array(:, :, :)  !< 传入的双精度数组
    integer, intent(in) :: xdim, ydim, zdim                 !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim, zdim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_3d_dp

  subroutine AllocMem_1d_si(Array, dim)
    !! 对传入的双精度数组进行分配内存, 如果报错可以给出报错信息
    integer(si), intent(inout), allocatable :: Array(:)  !< 传入的双精度数组
    integer, intent(in) :: dim                        !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(dim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_1d_si

  subroutine AllocMem_2d_si(Array, xdim, ydim)
    !! 对传入的二维双精度数组进行分配内存, 如果报错可以给出报错信息
    integer(si), intent(inout), allocatable :: Array(:, :)  !< 传入的双精度数组
    integer, intent(in) :: xdim, ydim                    !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_2d_si

  subroutine AllocMem_3d_si(Array, xdim, ydim, zdim)
    !! 对传���的三维双精度数组进行分配��存, 如果报错可以给出报错信��
    integer(si), intent(inout), allocatable :: Array(:, :, :)  !< 传入的��精度数组
    integer, intent(in) :: xdim, ydim, zdim                 !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim, zdim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_3d_si

  subroutine AllocMem_1d_li(Array, dim)
    !! 对传入的双精度数组进行分配内存, 如果报错可以给出报错信息
    integer(li), intent(inout), allocatable :: Array(:)  !< 传入的双精度数组
    integer, intent(in) :: dim                        !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(dim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_1d_li

  subroutine AllocMem_2d_li(Array, xdim, ydim)
    !! 对传入���二维双精度数组进行��配内存, 如果报错可以给出报错����息
    integer(li), intent(inout), allocatable :: Array(:, :)  !< 传入的双精度数���
    integer, intent(in) :: xdim, ydim                    !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_2d_li

  subroutine AllocMem_3d_li(Array, xdim, ydim, zdim)
    !! 对传入的三维双精度数组进行分配内存, 如果报错可以给出报错信息
    integer(li), intent(inout), allocatable :: Array(:, :, :)  !< 传入的双精度数组
    integer, intent(in) :: xdim, ydim, zdim                 !< 维度大小
    integer :: allo_flag                !< 判断内存分配是否成功的整数
    character(len=256) :: err_Message   !< 编译器传出的错误信息

    if (allocated(Array)) then
      deallocate(Array, stat=allo_flag, errmsg=err_Message)
      if (allo_flag /= 0) then
        write(stdout, '(A)') "Error occuring when deallocate array. The error message is:"
        write(stdout, '(A)') trim(err_Message)
        stop
      end if
    end if

    allocate(Array(xdim, ydim, zdim), stat=allo_flag)
    if (allo_flag /= 0) then
      write(stdout, '(A)') "Error occuring when allocating array. The error message is:"
      write(stdout, '(A)') trim(err_Message)
      stop
    end if

  end subroutine AllocMem_3d_li

end module fileIO
