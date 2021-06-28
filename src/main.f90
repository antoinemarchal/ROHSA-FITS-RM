!! Regularized Optimization for Hypers-spectral Analysis (ROHSA)
program ROHSA
! USE OMP_LIB

  use mod_constants
  use mod_start
  use mod_inout
  use mod_rohsa
  use mod_fits
  
  implicit none

  logical :: noise           !! if false --> STD map computed by ROHSA with lstd and ustd (if true given by the user)
  logical :: regul           !! if true --> activate regulation
  logical :: descent         !! if true --> activate hierarchical descent to initiate the optimization
  logical :: save_grid       !! save grid of fitted parameters at each step of the multiresolution process
  logical :: lym             !! if true --> activate 2-Gaussian decomposition for Lyman alpha nebula emission
  logical :: init_spec       !! if true --> use params mean spectrum with input
  logical :: init_grid       !! if true --> use fileinit to give the initialization of the last grid
  logical :: norm_var        !! if true --> normalize the var sig energy term                                                 
  
  integer :: n_gauss         !! number of gaussian to fit
  integer :: n_gauss_add     !! number of gaussian to add at each step
  integer :: m               !! number of corrections used in the limited memory matrix by LBFGS-B
  integer :: lstd            !! lower bound to compute the standard deviation map of the cube (if noise .eq. false)
  integer :: ustd            !! upper bound to compute the standrad deviation map of the cube (if noise .eq. false)
  integer :: iprint          !! print option 
  integer :: iprint_init     !! print option init
  integer :: maxiter         !! max iteration for L-BFGS-B alogorithm
  integer :: maxiter_init    !! max iteration for L-BFGS-B alogorithm (init mean spectrum)

  real(xp) :: lambda_amp     !! lambda for amplitude parameter
  real(xp) :: lambda_mu      !! lamnda for mean position parameter
  real(xp) :: lambda_sig     !! lambda for dispersion parameter

  real(xp) :: lambda_var_amp !! lambda for variance amplitude parameter
  real(xp) :: lambda_var_mu  !! lambda for variance mean position parameter
  real(xp) :: lambda_var_sig !! lambda for variance dispersion parameter

  real(xp) :: lambda_lym_sig !! lambda for variance dispersion parameter

  real(xp) :: amp_fact_init  !! times max amplitude of additional Gaussian
  real(xp) :: sig_init       !! dispersion of additional Gaussian
  real(xp) :: lb_sig_init    !! lower bound sigma init
  real(xp) :: ub_sig_init    !! upper bound sigma init
  real(xp) :: lb_sig         !! lower bound sigma
  real(xp) :: ub_sig         !! upper bound sigma

  real(xp) :: sig_rmsf !!fixed dispersion of the rmsf

  character(len=512) :: filename_parameters !! name of the parameters file (default parameters.txt)
  character(len=512) :: filename            !! name of the data file
  character(len=512) :: fileout             !! name of the output result
  character(len=512) :: timeout             !! name of the output result
  character(len=512) :: filename_noise      !! name of the file with STD map (if noise .eq. true)
  character(len=512) :: filename_init_spec
  character(len=512) :: fileinit           !! name of the file with init last grid
  character(len=8)   :: init_option !!Init ROHSA with the mean or the std spectrum    

  real(xp) :: start, finish

  real(kind=4), allocatable, dimension(:,:,:) :: array
  real(kind=4), dimension(:,:), allocatable :: std_data    !! standard deviation map fo the cube is given by the user 
  real(xp), dimension(:,:,:), allocatable :: data        !! initial fits data
  real(xp), dimension(:,:,:), allocatable :: data_init   !! initial fits data init full grid
  real(xp), dimension(:,:), allocatable   :: std_cube    !! standard deviation map fo the cube is given by the user 
  real(xp), dimension(:), allocatable     :: params_init     !! standard deviation map fo the cube is given by the user 
  real(xp), dimension(:,:,:), allocatable :: grid_params !! parameters to optimize at final step (dim of initial cube)
  real(kind=4), allocatable, dimension(:,:,:) :: grid_fits
  integer, dimension(3) :: dim_data !! dimension of reshape cube

  integer :: thread_id, i

  ! !$OMP PARALLEL PRIVATE(thread_id)
  ! thread_id = OMP_GET_THREAD_NUM()

  ! DO i=0,OMP_GET_MAX_THREADS()
  !    IF (i == thread_id) THEN
  !       PRINT *, "Hello from process: ", thread_id
  !    END IF
  !    !$OMP BARRIER
  ! END DO
  ! !$OMP END PARALLEL
     
  call cpu_time(start)

  !Print header and get filename in argument
  call get_command_argument(1, filename_parameters)
    
  !Default user parameters
  n_gauss = 6

  lambda_amp = 1._xp
  lambda_mu = 1._xp
  lambda_sig = 1._xp
  lambda_var_sig = 1._xp
  lambda_lym_sig = 0._xp

  amp_fact_init = 2._xp/3._xp
  sig_init = 5._xp
  lb_sig_init = 0.001_xp
  ub_sig_init = 100._xp
  lb_sig = 0.001_xp
  ub_sig = 100._xp
  maxiter_init = 15000
  maxiter = 800
  m = 10
  noise = .false.
  regul = .true.
  descent = .false.
  lstd = 1; ustd = 20
  init_option = "mean"
  iprint = -1
  iprint_init = -1
  save_grid = .true.
  lym = .false.
  init_grid = .false.
  init_spec = .false.
  norm_var = .false.

  sig_rmsf = 1._xp
 
  !Read parameters
  call read_parameters(filename_parameters, filename, fileout, timeout, filename_noise, filename_init_spec, &
       n_gauss, lambda_amp, lambda_mu, lambda_sig, lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, &
       amp_fact_init, sig_init, lb_sig_init, ub_sig_init, lb_sig, ub_sig, init_option, maxiter_init, maxiter, &
       m, noise, regul, descent, lstd, ustd, iprint, iprint_init, save_grid, lym, init_grid, fileinit, init_spec, &
       norm_var, sig_rmsf)

  !Call header
  call header()  

  print*, "filename = '",trim(filename),"'"

  !Load data
  call read_fits(filename, array)
  call roll_fits(array,data)

  if (init_spec .eqv. .true.) then
     call read_array(filename_init_spec, params_init)
  end if

  !Load grid init
  if (init_grid .eqv. .true.) then
     call read_output(fileinit, data_init, n_gauss, shape(data))
  end if

  if (noise .eqv. .true.) then
     if (filename_noise == " ") then
        print*, "--> noise = .true. (no input rms map)"
     end if
     call read_fits2D(filename_noise, std_data)
     call roll_fits2D(std_data,std_cube)
  end if

  !Call ROHSA subroutine
  dim_data = shape(data)
  allocate(grid_params(3*n_gauss, dim_data(2), dim_data(3)))
  allocate(grid_fits(dim_data(3), dim_data(2),3*n_gauss))

  call main_rohsa(data, std_cube, grid_params, fileout, timeout, n_gauss, lambda_amp, lambda_mu, lambda_sig, &
       lambda_var_amp, lambda_var_mu, lambda_var_sig, lambda_lym_sig, amp_fact_init, sig_init, lb_sig_init, &
       ub_sig_init, lb_sig, ub_sig, maxiter_init, maxiter, m, noise, regul, descent, lstd, ustd, init_option, &
       iprint, iprint_init, save_grid, lym, init_grid, fileinit, data_init, params_init, init_spec, norm_var, &
       sig_rmsf)  

  do i=1, n_gauss
     grid_params(2+(3*(i-1)),:,:) = grid_params(2+(3*(i-1)),:,:) - 1._xp
  end do

  !Write output fits file
  print*, " "
  print*, "_____ Write output file _____"
  print*, " "
  call unroll_fits(grid_params, grid_fits)
  call writefits3D(fileout,grid_fits,dim_data(3),dim_data(2),3*n_gauss)

  call ender()

  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start
   
end program ROHSA
