module computePE
  use global_var, only: wp, stdout
  use fileIO, only: model, eamAlloyFile
  use spline, only: func, dfunc
  implicit none
  private
  public :: Energy

  !logical, private, save :: isNeighbored = .false. !! 用来判断当前模型的近邻列表是否已经建立, 避免反复建立邻居列表增加运行时间

  interface Energy
    module procedure ::  EamAlloyEnergy
  end interface Energy

contains

  !! 考虑周期性边界条件计算2个原子之间的距离
  pure function Distance(Atoms, i, j) result(res)
    type(Model), intent(in) :: Atoms     !< 存储原子坐标和邻居列表的ATOM类数组
    integer, intent(in) :: i, j                       !< 2个原子之间的下标
    real(wp) :: dx, dy, dz                            !< xyz方向上的差距, 考虑周期性边界条件
    real(wp), dimension(3) :: res                     !< 考虑周期性边界条件下第i个与第j个原子之间的距离

    dx = Atoms%r(j, 1) - Atoms%r(i, 1)
    dx = dx - nint(dx / Atoms%lx) * Atoms%lx
    dy = Atoms%r(j, 2) - Atoms%r(i, 2)
    dy = dy - nint(dy / Atoms%ly) * Atoms%ly
    dz = Atoms%r(j, 3) - Atoms%r(i, 3)
    dz = dz - nint(dz / Atoms%lz) * Atoms%lz

    res(1) = dx; res(2) = dy; res(3) = dz

  end function Distance

  subroutine EamAlloyEnergy(eampot, Atoms, f, pe)
    !! 根据eam势函数文件与模型文件计算模型的势能
    type(eamAlloyFile), intent(in) :: eampot            !< 一个包含了势函数所有数据和处理好数据之后的类
    type(Model), intent(in):: Atoms                      !< 存储原子坐标和邻居列表的ATOM类数组
    real(wp), allocatable, intent(inout) :: f(:, :)                      !< 输出的模型的能量
    real(WP), allocatable, intent(inout)  :: pe(:)  !< 存储每个原子在xyz方向上面的受力
    real(wp), dimension(3) :: rij     !< 两个原子之间的距离的向量rij
    real(WP), allocatable  :: Fp(:)
    real(WP) :: s_fx, s_fy, s_fz      !< 计算每个原子的受力
    real(WP) :: f12_x, f21_x, f12_y, f21_y, f12_z, f21_z
    real(WP) :: s_pe                  !< 计算每个原子的势能
    real(WP) :: Fp1, Fp2              !< 第i个和第j个原子之间的嵌入能对电子密度的导数
    real(WP) :: rhop1, rhop2          !< 第i个和第j个原子之间的电子密度对距离的导数
    real(wp) :: dr                    !< 两个原子之间的距离
    real(WP) :: rho                   !< 电子密度
    real(WP) :: phi, phip             !< 对势φ对r的导数
    real(WP) :: inv_dr              !< 变量dR的倒数   
    real(WP) :: phi_temp
    integer :: cen_type !< 对应中心原子的type
    integer :: neb      !< 第i个原子的邻居原子数目
    integer :: neb_id   !< 第i个原子的邻居原子下标
    integer :: neb_type !< 第i个原子的邻居原子类型
    integer :: ierror   !< 通过查询这个参数查询插值是否成功
    integer :: i, j

    !! 初始化每一个变量和分配内存
    if (.not. allocated(f)) allocate(f(3, Atoms%nAtoms))
    if (.not. allocated(pe)) allocate(pe(Atoms%nAtoms))
    allocate(Fp(Atoms%nAtoms))
    !! 第一层循环遍历每一个原子, 第二层循环遍历每一个原子的近邻
    !$omp parallel num_threads(2) default(private) shared(atoms, eampot, pe, Fp, f)
    !$omp do
    do i = 1, Atoms%nAtoms
      cen_type = Atoms%typeID(i)   !! 获取第i个原子的类型
      neb = Atoms%IDs(i)%neb_n     !! 获取第i个原子的邻居数量
      rho = 0.0_WP
      !! 计算嵌入能和嵌入能F对电子密度ρ的导师
      do j = 1, neb
        neb_id = Atoms%IDs(i)%ID(j)
        neb_type = Atoms%typeId(neb_id)                     !! 获取第i个原子的第j个邻居原子的类型
        rij = Distance(Atoms, i, neb_id)                    !! 考虑周期性边界条件计算原子之间向量rij
        dr = norm2(rij)                                     !! 计算2个粒子之间的绝对距离
        rho = rho + func(eampot%frho(neb_type), dr, ierror) !! 计算2个原子之间的电子密度
      end do
      Fp(i) = dfunc(eampot%fF(cen_type), rho, ierror)       !! 计算出每个原子嵌入能F对电子密度ρ的导数
      pe(i) = func(eampot%fF(cen_type), rho, ierror)        !! 将嵌入能放入到pe数组中
    end do
    !$omp end do

    !! 开始第二次遍历整个模型, 计算每个原子的受力
    !$omp do
    do i = 1, Atoms%nAtoms
      cen_type = Atoms%typeID(i)   !! 获取第i个原子的类型
      neb = Atoms%IDs(i)%neb_n     !! 获取第i个原子的邻居数量
      Fp1 = Fp(i)                   
      s_fx = 0.0_WP;  s_fy = 0.0_WP;  
      s_fz = 0.0_WP;  s_pe = 0.0_WP;
      !if (neb == 0) f(1, i)=0.0_WP; f(2, i)=0.0_WP; f(3, i)=0.0_WP; pe(i)=0.0_WP; cycle
      do j = 1, neb
        neb_id = Atoms%IDs(i)%ID(j)
        neb_type = Atoms%typeId(neb_id)     !! 获取第i个原子的第j个邻居原子的类型
        Fp2 = Fp(neb_id)
        rij = Distance(Atoms, i, neb_id)    !! 获取rij的三个分量
        dr = norm2(rij)      
        inv_dr = 1.0_WP / dr  
        !! 计算第i和第j个原子之间的两体势φ(r)和对应的导数
        phi_temp = func(eampot%fphi(cen_type, neb_type), dr, ierror)
        phi = phi_temp * inv_dr  
        !! 这里注意从势函数文件中读取到的是Z(r) = φ(r)*r
        !! 所以φ`(r) = (Z`(r)*r - Z(r))/r**2
        phip = (dfunc(eampot%fphi(cen_type, neb_type), dr, ierror) * dr -  &
                phi_temp ) * inv_dr**2 * 0.5_WP
        !! 计算第i个原子和第j个原子之间互相由于对方产生的电子密度的导数   
        rhop1 = dfunc(eampot%frho(cen_type), dr, ierror)  
        rhop2 = dfunc(eampot%frho(neb_type), dr, ierror)  
        !! 进行数据处理
        phip = phip * inv_dr
        rhop1 = rhop1 * inv_dr
        rhop2 = rhop2 * inv_dr
        f12_x = rij(1) * (phip + Fp1 * rhop2)
        f12_y = rij(2) * (phip + Fp1 * rhop2)
        f12_z = rij(3) * (phip + Fp1 * rhop2)
        f21_x = -rij(1) * (phip + Fp2 * rhop1)
        f21_y = -rij(2) * (phip + Fp2 * rhop1)
        f21_z = -rij(3) * (phip + Fp2 * rhop1)
        
        s_pe = s_pe + phi
        s_fx = s_fx + f12_x - f21_x
        s_fy = s_fy + f12_y - f21_y
        s_fz = s_fz + f12_z - f21_z
      end do
      !! 数据归总到统一的数组里面
      f(1, i) = s_fx; f(2, i) = s_fy; 
      f(3, i) = s_fz; pe(i) = pe(i) + s_pe * 0.5_WP
    end do
    !$omp end do
    !$omp end parallel    

  end subroutine EamAlloyEnergy
end module computePE
