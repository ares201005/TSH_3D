program surfacehopping
  use dynamik
  implicit none
  
  integer :: trajectory
  ! parameters read from input file
  integer :: N_traj
  integer :: state_init
  double precision :: r_init(3), p_init(3)
  double precision :: mass
  double precision :: dE, angle_dp
  double precision :: upper_limit
  double precision :: dt
  double precision :: t_snapshot, old_time, new_time
  ! variables for statistics
  integer :: timestep_counter,fname
  integer :: regular_traj, irregular_traj, end_state(2)
  double precision :: deltaE_max, deltaNorm_max, population(2)


  open(100,file='trajinfo.out',status='replace')

  ! read input parameters
  call read_input_parameters(N_traj, state_init, r_init, p_init, mass, dE, angle_dp, upper_limit, dt, t_snapshot)

  ! read potentials
  call read_orig_pes

  ! set random seed
  call set_nran
  
  ! initialize statistics variables
  regular_traj = 0
  end_state(:) = 0
  population(:) = 0.d0
  deltaE_max = 0.d0
  deltaNorm_max = 0.d0

  all_trajs: do trajectory = 1,N_traj
     fname=200
     timestep_counter = 0
     call set_initial_cond(state_init, r_init, p_init, mass, dE, angle_dp, dt)
     single_traj: do

        old_time = ret_actual_time()
        call timestep
        new_time = ret_actual_time()
        !projectile leaves region of interest
        if ((abs(ret_actual_posX()).gt.abs(upper_limit)).OR.(abs(ret_actual_posY()).gt.abs(upper_limit)).OR.(abs(ret_actual_posZ()).gt.abs(upper_limit))) then
           regular_traj = regular_traj + 1
           end_state(ret_actual_state()) = end_state(ret_actual_state()) + 1
           population(:) = population(:) + ret_actual_population()
           exit
        end if
!        fname=340
!        if ( (old_time.le.t_snapshot).AND.(new_time.gt.t_snapshot) ) call print_coords(fname,trajectory)
!        timestep_counter = timestep_counter + 1
!        if ( (timestep_counter.eq.100).AND.(fname.lt.1000) )then
!           if (fname.eq.340) call print_coords(fname,trajectory)
!!           call print_coords(fname,trajectory)
!           fname = fname + 1
!           timestep_counter = 0
!        end if

     end do single_traj

     call print_traj_info()
     if (deltaE_max.lt.ret_deltaE_max()) deltaE_max = ret_deltaE_max()
     if (deltaNorm_max.lt.ret_deltaNorm_max()) deltaNorm_max = ret_deltaNorm_max()
     
     if ((trajectory.gt.(N_traj/10.d0)).AND.(regular_traj.eq.0)) exit

  end do all_trajs
  
  !print statistics
  write(100,*) '* INPUT******************************************************'
  write(100,*) '* ', N_traj
  write(100,*) '* ', state_init
  write(100,*) '* ', sngl(r_init)
  write(100,*) '* ', sngl(p_init)
  write(100,*) '* ', sngl(mass)
  write(100,*) '* ', sngl(dE), sngl(angle_dp)
  write(100,*) '* ', sngl(upper_limit)
  write(100,*) '* ', dt
  write(100,*) '* SUMMARY****************************************************'
  write(100,*) '* stopped after ', trajectory, 'trajectories'
  write(100,*) '* deltaE_total_max ', deltaE_max
  write(100,*) '* deltaNorm_total_max ', deltaNorm_max
  write(100,*) '* regular trajectories ', regular_traj
  write(100,*) '* final states', end_state(:)
  if (regular_traj.ne.0) then
     write(100,*) '* final percentage', sngl(dble(end_state(:))/dble(regular_traj))
     write(100,*) '* final population', sngl(population(:)/dble(regular_traj))
  else
     write(100,*) '* final percentage', sngl(dble(end_state(:)))
     write(100,*) '* final population', sngl(dble(population(:)))
  end if
  write(100,*) '*************************************************************'
  close(100)

end program surfacehopping
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!----------read initial conditions
subroutine read_input_parameters(Ntraj,state_init,r_init,p_init,mass,dE,angle_dp,upper_limit,dt,t_snapshot)
  implicit none

  integer :: Ntraj                                      ! number of trajectories per initial conditions
  integer :: state_init                                 ! initial state
  double precision :: r_init(3), p_init(3)              ! initial position and initial momentum
  double precision :: mass                              ! mass of particle
  double precision :: dE, angle_dp                      ! energy width of beam and beam divergence (angle!)
  double precision :: upper_limit                       ! boarder at which trajectory is stopped
  double precision :: dt                                ! time step for integration of cc equations and equation of motion
  double precision :: t_snapshot                        ! time at which position of particle is written

  open(99,file='input.inp',status='old')
  read(99,*) Ntraj
  read(99,*) state_init
  read(99,*) r_init(:)
  read(99,*) p_init(:)
  read(99,*) mass
  read(99,*) dE
  read(99,*) angle_dp
  read(99,*) upper_limit
  read(99,*) dt
  read(99,*) t_snapshot
  close(99)

  return
end subroutine read_input_parameters
