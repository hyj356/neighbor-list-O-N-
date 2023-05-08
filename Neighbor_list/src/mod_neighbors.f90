module neighbor
  use global_var, only: wp, stdout
  use fileIO,     only: Model
  implicit none
  
  private
  public :: find_neighbor, BuildNeighbors

contains
  !>
  !! 考虑周期性边界条件计算2个原子之间的距离
  !<
  pure function Distance(Atoms, i, j) result(res)
    type(Model), intent(in) :: Atoms     !< 存储原子坐标和邻居列表的ATOM类数组
    integer, intent(in) :: i, j          !< 2个原子之间的下标
    real(wp) :: dx, dy, dz               !< xyz方向上的差距, 考虑周期性边界条件
    real(wp) :: res                      !< 考虑周期性边界条件下第i个与第j个原子之间的距离

    dx = Atoms%r(i, 1) - Atoms%r(j, 1)
    dx = dx - nint(dx / Atoms%lx) * Atoms%lx
    dy = Atoms%r(i, 2) - Atoms%r(j, 2)
    dy = dy - nint(dy / Atoms%ly) * Atoms%ly
    dz = Atoms%r(i, 3) - Atoms%r(j, 3)
    dz = dz - nint(dz / Atoms%lz) * Atoms%lz

    res = dx * dx + dy * dy + dz * dz

  end function Distance

  !>
  !!  此函数用于计算第i个原子属于哪个盒子
  !<
  pure function find_cell(atoms, i, icutoff, numCell) result(res)
    type(Model), intent(in) :: atoms      !< 包含盒子和每个原子的xyz坐标信息的自定义类
    integer, intent(in)     :: i          !< 第i个原子归属于那个盒子里面
    integer, intent(in)     :: numCell(4) !< 盒子在xyz上可以被分成多少段, 第四个元素表示总共有几个盒子
    real(WP), intent(in)    :: icutoff    !< 盒子的边长的倒数, 这里我们默认盒子3个边长完全一样
    real(WP) :: dr
    integer, dimension(4)   :: res        !< 返回盒子在xyz上被划分成几份

    ! 计算第i个原子在xyz上分布属于第几个盒子, 向下取整
    ! 因为我们默认盒子的id是从1开始的
    dr = atoms%r(i, 1) - atoms%xlo  ! 获取某一个坐标轴距离下限的距离
    res(1) = floor(dr * icutoff)

    ! 这是为了防止有部分原子刚好在边界上, 此时dr = 0
    ! 而ceiling(0) = 0, 这是不合理的, 所以要对这种情况进行修正, 以下同理
    if(res(1) < 0)           res(1) = res(1) + numCell(1)    
    if(res(1) >= numCell(1)) res(1) = res(1) - numCell(1)
    
    dr = atoms%r(i, 2) - atoms%ylo  
    res(2) = floor(dr * icutoff)
    if(res(2) < 0)           res(2) = res(2) + numCell(2)    
    if(res(2) >= numCell(2)) res(2) = res(2) - numCell(2)

    dr = atoms%r(i, 3) - atoms%zlo  
    res(3) = floor(dr * icutoff)
    if(res(3) < 0)           res(3) = res(3) + numCell(3)    
    if(res(3) >= numCell(3)) res(3) = res(3) - numCell(3)

    ! 计算第i个原子归属于哪个cell
    res(4) = res(1) + numCell(1)  * (res(2) + numCell(2) * res(3))

  end function find_cell

  !>
  !!  此子程序用于构建每个原子的邻居ID列表
  !<
  subroutine find_neighbor(atoms, cutoff)
    type(Model), intent(inout) :: atoms   !< 包含盒子和每个原子的xyz坐标信息的自定义类
    real(WP), intent(in)    :: cutoff     !< 小盒子的边长
    real(WP) :: icutoff                   !< cutoff的导数
    real(WP) :: cutoff2                   !< cutoff的平方
    integer  :: numcells(4)
    integer  :: cell(4)
    integer  :: tmp(3)
    integer  :: i, j, k, n1, m, n2        !< 循环下标
    integer  :: neighbor_cell             !< 邻居cell的id
    integer, allocatable :: cell_count(:) !< 记录每个cell里面一共有多少个原子以确定循环次数
    integer, allocatable :: cell_count_sum(:) !< cell_count_sum的前缀和数组
    integer, allocatable :: cell_contents(:)  !< 记录每个原子对应的cellid

    ! 数据处理
    icutoff = 1.0_WP / cutoff
    cutoff2 = cutoff * cutoff

    ! 计算xyz方向上共有几个小盒子, 向下取整
    numcells(1) = floor(atoms%lx * icutoff)
    numcells(2) = floor(atoms%ly * icutoff)
    numcells(3) = floor(atoms%lz * icutoff)
    numcells(4) = numcells(1) * numcells(2) * numcells(3)
    !write(stdout, *) numcells(4)

    ! 分配内存并初始化
    allocate(cell_count(0:numcells(4)-1), cell_count_sum(0:numcells(4)-1), cell_contents(0:atoms%nAtoms-1))
    call InitialList(atoms)
    cell_count = 0; cell_count_sum = 0; cell_contents = 0

    !$omp parallel num_threads(2) default(firstprivate) &
    !$omp& shared(atoms, cutoff2, cell_count, cell_count_sum, cell_contents) &
    !$omp& private(n2, neighbor_cell, cell, tmp)
    ! 循环每个原子找到其对应的盒子ID
    !$omp do
    do i = 1, atoms%nAtoms
      cell = find_cell(atoms, i, icutoff, numcells)
      cell_count(cell(4)) = cell_count(cell(4)) + 1
    end do
    !$omp end do

    !$omp single
    ! 构建前缀和数组
    do i = 1, numcells(4) - 1
      cell_count_sum(i) = cell_count_sum(i - 1) + cell_count(i - 1)
    end do

    ! 清零列表
    cell_count = 0
    !$omp end single

    ! 重新累加列表
    !$omp do
    do i = 1, atoms%nAtoms
      cell = find_cell(atoms, i, icutoff, numcells)
      cell_contents(cell_count_sum(cell(4)) + cell_count(cell(4))) = i
      cell_count(cell(4)) = cell_count(cell(4)) + 1
    end do
    !$omp end do

    ! 接着开始循环构建邻居列表
    !$omp do
    do n1 = 1, atoms%nAtoms
      cell = find_cell(atoms, n1, icutoff, numcells) ! 找到当前原子对应的cell的id
      do k = -1, 1
        do j = -1, 1
          do i = -1, 1
            ! 考虑周期性条件对xyz的id进行修正
            tmp(1) = i; tmp(2) = j; tmp(3) = k
            tmp = cell(1:3) + tmp
            tmp = MODULO(tmp, numcells(1:3))    ! 对xyz的cellid进行求余操作
            neighbor_cell = tmp(1) + numCells(1) * (tmp(2) + numCells(2) * tmp(3))
            do m = 1, cell_count(neighbor_cell)
              ! 这里减一是因为原子的下标是从1开始的
              n2 = cell_contents(cell_count_sum(neighbor_cell) + m - 1)   ! 邻居原子的id
              ! 这里运用了verlet列表的思想减小计算量
              if (n1 < n2) then
                ! 考虑周期性边界条件计算模型, 如果符合条件, 就加入到邻居列表里面
                if (Distance(atoms, n1, n2) < cutoff2) then
                  !$omp critical
                  ! 对可分配数组的增减操作同一时间只能有一个线程进行操作
                  atoms%IDs(n1)%ID = [atoms%IDs(n1)%ID, n2]
                  atoms%IDs(n2)%ID = [atoms%IDs(n2)%ID, n1]
                  !$omp end critical
                  Atoms%IDs(n1)%neb_n = Atoms%IDs(n1)%neb_n + 1 ! 邻居数加一
                  Atoms%IDs(n2)%neb_n = Atoms%IDs(n2)%neb_n + 1 ! 邻居数加一
                end if
              end if
            end do
          end do
        end do 
      end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine find_neighbor

    !! 根据传入的截断半径的平方构建每个原子在截断半径内的近邻原子
  subroutine BuildNeighbors(Atoms, cutoff2)
    type(Model), intent(inout) :: Atoms    ! 原子模型
    real(WP), intent(in) :: cutoff2     ! 寻找邻居时的截断半径的平方
    real(WP) :: rij2                    ! 第i和第j个原子之间距离的平方
    integer  :: i, j

    call InitialList(Atoms)
    do i = 1, Atoms%nAtoms - 1
      do j = i + 1,  Atoms%nAtoms
        rij2 = Distance(Atoms, i, j)
        if (rij2 < cutoff2) then
          Atoms%IDs(i)%ID = [Atoms%IDs(i)%ID, j]    !! Fortran2003特性, 类似于python的list类型的append, 但因涉及反复申请内存, 效率极低
          Atoms%IDs(j)%ID = [Atoms%IDs(j)%ID, i]    !! 比较好的折中办法是使用指针和派生数据类型搭建链表,
          Atoms%IDs(i)%neb_n = Atoms%IDs(i)%neb_n + 1   !! 第i和第j个原子的邻居数量加一
          Atoms%IDs(j)%neb_n = Atoms%IDs(j)%neb_n + 1
        end if
      end do
    end do

  end subroutine BuildNeighbors

  !! 此子程序用于初始化每个原子的邻居列表的内存
  pure subroutine InitialList(self)
    type(Model), intent(inout) :: self
    integer :: i 

    if (.not. allocated(self%IDs)) then
      ! 如果没有分配内存, 就给每个邻居列表分配一个0内存
      allocate(self%IDs(self%nAtoms))
      do concurrent(i = 1:self%nAtoms)
        allocate(self%IDs(i)%ID(0))
      end do
    else
      ! 如果已经分配了内存, 就释放内存再分配0内存, 并将邻居数量归零
      do concurrent(i = 1:self%nAtoms)
        deallocate(self%IDs(i)%ID)
        allocate(self%IDs(i)%ID(0))
        self%IDs(i)%neb_n = 0
      end do
    end if

  end subroutine InitialList
end module neighbor