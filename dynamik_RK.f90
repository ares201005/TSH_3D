module dynamik
  use read_data
  implicit none
  save

  ! public subroutines/functions from module read_data
  public :: read_orig_pes
  ! public subroutines/functions concerning the dynamics
  public :: set_nran, set_initial_cond, timestep,ret_actual_posX,ret_actual_posY,ret_actual_posZ,ret_actual_state,ret_deltaE_max,ret_deltaNorm_max,ret_actual_population,ret_actual_time
  ! public subroutines/functions for creating ouput
  public :: printtraj, printstat, print_traj_info, print_Pmat, print_coords

  private
  integer :: nran ! argument of the random number generator
  double precision :: mue ! masses of part 1 and 2, total mass, reduced mass
  double precision :: stepwidth ! length of time step
  
  ! all the following quantities are updated after every timestep (in subroutine timestep) !!!
  integer :: actual_state,jump_counter ! defines in which state the particle is in
  double precision :: actual_pos(3) ! position of particle
  double precision :: start_pos(3)  ! initial position of particle
  double precision :: actual_vel(3) ! velocity of particle
  double precision :: start_vel(3)  ! initial velocity of particle
  double precision :: V_mat(N_states,N_states) ! potential matrix in adiabatic basis
  double precision :: P_mat(N_states,N_states,3) ! matrix with the non-adiabatic coupling elements
  double precision :: actual_time
  double precision :: b_mat(N_states,N_states), g_mat(N_states,N_states)
  double precision :: max_E_total, min_E_total, max_Norm, min_Norm
  double precision :: actual_phase, actual_phase_v2
  complex*16 :: a_mat(N_states,N_states)
  
contains
  
!---------- set the initial value for random number generation
  subroutine set_nran
    implicit none
    
    integer :: time

    nran = -int(time())
    !set argument for random number generator
!    nran = -325456  

  end subroutine set_nran


!---------- set initial conditions for the dynamic of the problem
  subroutine set_initial_cond(state_init, r_init, p_init, mass, dX, dP, dt)
    implicit none

    integer :: state_init
    double precision :: r_init(3), p_init(3), mass, dX, dP, dt
    ! variables for energy and direction spread
    double precision :: energy_old, energy_new, dE_new, direction(3)

    ! set mass
    mue = mass

    ! set position
    actual_pos(:) = r_init(:)
    actual_pos(1) = actual_pos(1) + gasdev(nran,dX)
    actual_pos(2) = actual_pos(2) + gasdev(nran,dX)
    start_pos(:) = actual_pos(:)

    ! set velocity
    actual_vel(:) = p_init(:)/mue
!    dP = 2.7726d0/dX ! comes from 4*ln(2)
!    actual_vel(1) = (p_init(1) + gasdev(nran,dP))/mue
!    actual_vel(2) = (p_init(2) + gasdev(nran,dP))/mue
!!    energy_old = 0.5d0 * dot_product(p_init(:),p_init(:))/mue
!!    energy_new = energy_old + gasdev(nran) * (dE / 27.211383d0)/2.35482d0
!!    direction(:) = p_init(:) / sqrt(dot_product(p_init(:),p_init(:)))
!!    actual_vel(:) = direction(:) * sqrt(2.d0 * energy_new / mue)
!!    start_vel(:) = actual_vel(:)
!!    !alternatively set momentum uncertainty
!!    actual_vel(1) = (p_init(1) + gasdev(nran,angle_dp))/mue
!!    actual_vel(2) = (p_init(2) + gasdev(nran,angle_dp))/mue
!!    actual_vel(3) = p_init(3)/mue
    start_vel(:) = actual_vel(:)

    ! set initial coefficients
    actual_state = state_init
    a_mat = (0.d0,0.d0)
    a_mat(actual_state,actual_state) = (1.d0,0.d0)

    ! set initial time
    actual_time=0.d0

    ! set width of timestep
    stepwidth=dt

    ! set the jump counter to zero
    jump_counter = 0
    
    ! set actual phase
    actual_phase = 0.d0
    actual_phase_v2 = 0.d0

    ! set potential matrix
    call set_V_mat(actual_pos)

    ! set coupling matrix
    call set_P_mat(actual_pos)

    ! set initial values for energy and norm conservation
    max_E_total = actual_total_energy(); min_E_total = actual_total_energy()
    max_Norm = actual_norm(); min_Norm = actual_norm()

  end subroutine set_initial_cond



!---------- calculate force in X,Y,Z direction
  function forceX(help,state,index)
    implicit none

    integer :: state,index
    double precision :: help(3),xplus(3),xminus(3),forceX,h

    xplus = help
    xminus = help

    h=1.d-5
    xplus(index) = xplus(index)+h
    xminus(index) = xminus(index)-h

    forceX = -(interpolate_energy(xplus,state) - interpolate_energy(xminus,state))/(2.d0*h)

  end function forceX

!---------- calculate V matrix at given distance r
  subroutine set_V_mat(x)
    implicit none

    integer :: i,j
    double precision :: x(3)
    
    V_mat = 0.d0

    do i=1,N_states
       V_mat(i,i) = interpolate_energy(x,i)
    end do

  end subroutine set_V_mat



!---------- calculate P matrix at given distance r
  subroutine set_P_mat(x)
    implicit none

    integer :: i,j,index
    double precision :: x(3)

    P_mat = 0.d0

    index = 0
    do i=1,N_states
       do j=i+1,N_states
          index = index + 1
          P_mat(i,j,:) = interpolate_nacv(x,index)
          P_mat(j,i,:) = -P_mat(i,j,:)
       end do
    end do

  end subroutine set_P_mat



!---------- driver routine for ONE TIMESTEP
  subroutine timestep
    implicit none

    double precision :: motion(6)
    double precision :: old_pos(3)

    !make sure the initial conditions are set

    !timestep in coupled channel equations
    call rk4_complex(a_mat,N_states,actual_time,stepwidth,a_mat,aDot)
    !timestep in coupled channel equations with phase correction
!    call rk4_complex(a_mat,N_states,actual_time,stepwidth,a_mat,aDot_phase)

    old_pos(:) = actual_pos(:)

    !timestep in equations of motion
    motion(1) = actual_vel(1)
    motion(2) = actual_vel(2)
    motion(3) = actual_vel(3)
    motion(4) = actual_pos(1)
    motion(5) = actual_pos(2)
    motion(6) = actual_pos(3)
    call rk4(motion,6,actual_time,stepwidth,motion,vrDot)
    actual_vel(1) = motion(1)
    actual_vel(2) = motion(2)
    actual_vel(3) = motion(3)
    actual_pos(1) = motion(4)
    actual_pos(2) = motion(5)
    actual_pos(3) = motion(6)

    !update V-matrix
    call set_V_mat(actual_pos)

    !update P-matrix
    call set_P_mat(actual_pos)

    !update time
    actual_time=actual_time+stepwidth

    !update phase
    actual_phase = actual_phase + mue * sqrt( dot_product(actual_vel(:),actual_vel(:)) ) * sqrt( dot_product((actual_pos(:)-old_pos(:)),(actual_pos(:)-old_pos(:))) )
    actual_phase_v2 = actual_phase_v2 + dot_product(actual_vel(:),(actual_pos(:)-old_pos(:))) * mue

    !decide whether a jump from one state to another is performed
    call jump

    !check total energy and norm
    if (actual_total_energy().gt.max_E_total) max_E_total = actual_total_energy()
    if (actual_total_energy().lt.min_E_total) min_E_total = actual_total_energy()
    if (actual_norm().gt.max_Norm) max_Norm = actual_norm()
    if (actual_norm().lt.min_Norm) min_Norm = actual_norm()

    !apply decoherence corrections to density matrix a_mat
!    call decoherence

    return
  end subroutine timestep



!---------- performs timestep of matrix a_mat
  subroutine aDot(x,y,dydx)
    implicit none

    integer :: k,j,l
    double precision :: x,help_pos
    complex*16 :: y(N_states,N_states),dydx(N_states,N_states),imag

!    help_pos = actual_pos + actual_vel * (x-actual_time)
!    call set_V_mat(help_pos)
    imag = cmplx(0.d0,1.d0)

    do k=1,N_states
       do j=1,N_states
          
          dydx(k,j) = (0.d0,0.d0)

          do l=1,N_states
             dydx(k,j) = dydx(k,j) - imag * ( y(l,j)*(V_mat(k,l) - imag*dot_product(actual_vel(:),P_mat(k,l,:)) ) &
                                            - y(k,l)*(V_mat(l,j) - imag*dot_product(actual_vel(:),P_mat(l,j,:)) ))
          end do

       end do
    end do

  end subroutine aDot



!---------- performs timestep of matrix a_mat with phase correction according to
!           N. Shenvi et al., J. Chem. Phys. 135, 024101 (2011)
    subroutine aDot_phase(x,y,dydx)
    implicit none

    integer :: k,j,l
    double precision :: x, tot_energy, kin_energy(N_states), mom_vec(N_states), V_mat_phase(N_states,N_states)
    complex*16 :: y(N_states,N_states), dydx(N_states,N_states), imag

    imag = cmplx(0.d0,1.d0)

    tot_energy = actual_total_energy()
    V_mat_phase(:,:) = 0.d0

    ! set mom_vec and V_mat_phase
    do k=1, N_states
       kin_energy(k) = tot_energy - V_mat(k,k)
       if (kin_energy(k).lt.0.d0) kin_energy(k) = 0.d0
       mom_vec(k) = sqrt(2.d0 * mue * kin_energy(k))
       V_mat_phase(k,k) = - sqrt(dot_product(actual_vel(:),actual_vel(:))) * mom_vec(k)   ! actual_vel and mom_v are pos.
    end do

!    write(125,*) actual_time,V_mat_phase(1,1),V_mat_phase(2,2)

    do k=1,N_states
       do j=1,N_states
          
          dydx(k,j) = (0.d0,0.d0)

          do l=1,N_states
             dydx(k,j) = dydx(k,j) - imag * ( y(l,j)*(V_mat_phase(k,l) - imag*dot_product(actual_vel(:),P_mat(k,l,:)) ) &
                                            - y(k,l)*(V_mat_phase(l,j) - imag*dot_product(actual_vel(:),P_mat(l,j,:)) ))
          end do

       end do
    end do

  end subroutine aDot_phase




!---------- performs timestep of x and v
  subroutine vrDot(x,y,dydx)
    implicit none
    
    double precision x,y(6),dydx(6)
    
    dydx(1) = forceX(y(4:6),actual_state,1)/mue
    dydx(2) = forceX(y(4:6),actual_state,2)/mue
    dydx(3) = forceX(y(4:6),actual_state,3)/mue
    dydx(4) = y(1)
    dydx(5) = y(2)
    dydx(6) = y(3)
    
  end subroutine vrDot



!---------- decides if a jump to another state is performed
! Tully, J Chem Phys 93, 1061 (1990)
  subroutine jump
    implicit none

    integer :: i, j, ii, jj, old_state, new_state, nacv_index
    double precision :: rand, sum, ekin, deltaE

!   set b_mat
    do i = 1, N_states
       do j = 1, N_states

          b_mat(i,j) = 2.d0 * aimag( conjg(a_mat(i,j)) * V_mat(i,j) ) - 2.d0 * real( conjg(a_mat(i,j) )&
                                                                             * dot_product(actual_vel(:),P_mat(i,j,:)))
       end do
    end do

!   calculate the probability to jump from state i,j (stored in g_mat(i,j))

    do i = 1, N_states
       do j = 1, N_states

          if (b_mat(j,i).eq.0.d0) then
             g_mat(i,j) = 0.d0
          else
             g_mat(i,j) = stepwidth * b_mat(j,i) / real(a_mat(i,i))          
          end if

          if (g_mat(i,j).lt.0.d0) g_mat(i,j)=0.d0

          ! check if g(i,i) = 0
          if (i.eq.j) then
             if (g_mat(i,j).ne.0.d0) then
                write(*,*) 'ACHTUNG'
!                g_mat(i,j) = 0.d0
             end if
          end if

       end do
    end do

    sum = 0.d0
    rand = ran2(nran)

    do j=1,N_states
       sum = sum + g_mat(actual_state,j)
       
       jump_probable: if (sum.ge.rand) then
          ekin = 0.5d0 * mue * dot_product(actual_vel,actual_vel)
          ! deltaE corresponds to the change of potential energy due to the jump -> deltaE_kin = -deltaE
          deltaE = V_mat(j,j) - V_mat(actual_state,actual_state)
!          deltaE = interpolate_energy(actual_pos,j) - interpolate_energy(actual_pos,actual_state)

          ! perform jump only if it is classically allowed
          jump_allowed: if (ekin.gt.deltaE) then
             
             !   set old_state
             old_state = actual_state

             ! set new actual state
             actual_state = j

             ! calculate number of coupling component responsible for the jump
             new_state = actual_state

             ! scale velocity
!             call scale_velocity_VEL(deltaE)
             call scale_velocity_NACV(deltaE,old_state,new_state)

             jump_counter = jump_counter + 1

          end if jump_allowed

          ! when the jump is performed (or not allowed) leave this subroutine and go on with the propagation
          return

       end if jump_probable

    end do

  end subroutine jump



!---------- scales the velocity vector ALONG actual_vel
!           so that it yields the new kinetic energy e_kin_new = e_kin - deltaE
  subroutine scale_velocity_VEL(deltaE)
    implicit none

    double precision :: ekin_new, ekin_old, v_new, deltaE, norm
    double precision :: direction(3)
    
    norm = dsqrt(dot_product(actual_vel,actual_vel))
    direction = actual_vel/norm
    ekin_old = 0.5d0 * mue * norm**2
    ekin_new = ekin_old - deltaE
    actual_vel = direction * dsqrt(2.d0*ekin_new/mue)
    
  end subroutine scale_velocity_VEL



!---------- scales the velocity vector ALONG the non-adiabatic coupling vector leading to the jump
!           so that it yields the new kinetic energy e_kin_new = e_kin - deltaE
  subroutine scale_velocity_NACV(deltaE,old_state,new_state)
    implicit none

    integer :: old_state,new_state
    double precision :: ekin_old, deltaE, PinV, alpha_plus, alpha_minus, discriminant
    double precision :: nacv(3), nacv_norm
    double precision :: vn1(3), vn2(3), dot1, dot2

    ekin_old = 0.5d0 * mue * dot_product(actual_vel,actual_vel)
    nacv(:)=P_mat(old_state,new_state,:)
    nacv_norm = dsqrt(dot_product(nacv,nacv))
        
    PinV = dot_product(nacv,actual_vel)
    discriminant = (PinV/nacv_norm**2)**2 - 2.d0*deltaE/(mue*nacv_norm**2)
    ! if there is no solution scale velocity vector along actual_vel
    if (discriminant.lt.0.d0) then
       call scale_velocity_VEL(deltaE)
       return
    end if

    ! there are two solutions for the new velocity v_new = v_old + alpha*nacv
    alpha_plus =  -PinV/nacv_norm**2 + sqrt(discriminant)
    alpha_minus = -PinV/nacv_norm**2 - sqrt(discriminant)

    ! select the solution with greater inner product
    vn1(:) = actual_vel(:) + alpha_plus * nacv(:)
    vn2(:) = actual_vel(:) + alpha_minus * nacv(:)
    dot1 = dot_product(vn1,actual_vel)
    dot2 = dot_product(vn2,actual_vel)

    actual_vel = vn1
    if (dot2.gt.dot1) actual_vel = vn2

  end subroutine scale_velocity_NACV



!---------- include decoherence corrections on a_mat as proposed in
!           Granucci et al. in J. Chem. Phys. 126, 134114, 2007
  subroutine decoherence()
    implicit none
    
    integer :: k,j,l
    double precision, parameter :: alpha=0.1d0
    double precision :: e_kin, sum_a, ll_factor
    double precision :: correction_mat(N_states,N_states)
    
    e_kin = 0.5d0 * dot_product(actual_vel,actual_vel) * mue
    
    correction_mat(:,:) = 1.d0
    
    sum_a = 0.d0
    
    ! set diagonal elements of correction matrix and calculate sum(a'(k,k) for all k ne actual_state) and the ll_factor
    do k=1, N_states
       if (k.ne.actual_state) then 
          correction_mat(k,k) = exp(-stepwidth * (e_kin/(e_kin + alpha))&
               * 2.d0 * abs(V_mat(k,k)-V_mat(actual_state,actual_state)))
          sum_a = sum_a + real(a_mat(k,k)) * correction_mat(k,k)
       end if
    end do
    
    ll_factor = (1.d0 - sum_a) / real(a_mat(actual_state,actual_state))
    correction_mat(actual_state,actual_state) = ll_factor
    
    
    ! set off-diagonal elements of correction matrix
    do k=1, N_states-1
       do j=k+1, N_states
          
          if ((k.ne.actual_state).AND.(j.ne.actual_state)) correction_mat(k,j)=&
               exp(-stepwidth * (e_kin/(e_kin + alpha))&
               * (abs(V_mat(k,k)-V_mat(actual_state,actual_state))&
               + abs(V_mat(j,j)-V_mat(actual_state,actual_state))) )
          
          if (k.eq.actual_state) correction_mat(k,j)=&
               exp(-stepwidth * (e_kin/(e_kin + alpha))&
               * abs(V_mat(j,j)-V_mat(actual_state,actual_state)) )&
               * sqrt(ll_factor)
          
          if (j.eq.actual_state) correction_mat(k,j)=&
               exp(-stepwidth * (e_kin/(e_kin + alpha))&
               * abs(V_mat(k,k)-V_mat(actual_state,actual_state)) )&
               * sqrt(ll_factor)
          
          correction_mat(j,k) = correction_mat(k,j)
          
       end do
    end do
    
    !apply corrections to the density matrix a_mat
    do k=1, N_states
       do j=1, N_states
          
          a_mat(k,j) = a_mat(k,j) * correction_mat(k,j)
          
       end do
    end do
    
  end subroutine decoherence
  
  

!---------- calculate the actual total energy
  function actual_total_energy()
    implicit none

    double precision :: actual_total_energy

    actual_total_energy = V_mat(actual_state,actual_state) + 0.5d0 * mue * dot_product(actual_vel,actual_vel)

    return
  end function actual_total_energy



!---------- calculate the actual norm
  function actual_norm()
    implicit none

    integer :: i
    double precision :: actual_norm

    actual_norm = 0.d0
    
    do i = 1, N_states
       actual_norm = actual_norm + real(a_mat(i,i))
    end do

    return
  end function actual_norm



!---------- return actual state
  function ret_actual_state()
    implicit none

    integer :: ret_actual_state

    ret_actual_state = actual_state

  end function ret_actual_state



!---------- return actual time
  function ret_actual_time()
    implicit none

    double precision :: ret_actual_time

    ret_actual_time = actual_time

  end function ret_actual_time



!---------- return actual diagonal elements of a_matrix
  function ret_actual_population
    implicit none

    integer :: i
    double precision :: ret_actual_population(N_states)

    do i=1, N_states
       ret_actual_population(i) = real(a_mat(i,i))
    end do

  end function ret_actual_population



!---------- return the x value of actual_pos
  function ret_actual_posX()
    implicit none
    
    double precision :: ret_actual_posX

    ret_actual_posX = actual_pos(1)

    return
  end function ret_actual_posX


!---------- return the y value of actual_pos
  function ret_actual_posY()
    implicit none
    
    double precision :: ret_actual_posY

    ret_actual_posY = actual_pos(2)

    return
  end function ret_actual_posY



!---------- return the z value of actual_pos
  function ret_actual_posZ()
    implicit none
    
    double precision :: ret_actual_posZ

    ret_actual_posZ = actual_pos(3)

    return
  end function ret_actual_posZ


!---------- return the maximum change of the total energy
  function ret_deltaE_max()
    implicit none

    double precision :: ret_deltaE_max

    ret_deltaE_max = max_E_total-min_E_total

    return
  end function ret_deltaE_max



!---------- return the maximum change in the total norm
  function ret_deltaNorm_max()
    implicit none

    double precision :: ret_deltaNorm_max

    ret_deltaNorm_max = max_Norm-min_Norm

    return
  end function ret_deltaNorm_max



!---------- print trajectory information after trajectory ended
  subroutine print_coords(fname,trajectory)
    implicit none
    
    integer :: fname, trajectory

!    if (trajectory.eq.1) then
!       open(fname,status='new')
!    else
       open(fname,status='old',position='append')
!    end if

!   if (fname.eq.340) write(*,*) actual_time

2000 format(4F20.14,I3)
!    write(fname,*) actual_pos(1), actual_pos(2), actual_state
    write(fname,2000) actual_pos(1), actual_pos(2), actual_phase, actual_phase_v2, actual_state

    close(fname,status='keep')

    return
  end subroutine print_coords



!---------- print trajectory information after trajectory ended
  subroutine print_traj_info()
    implicit none
    
1000 format(2I3,14F20.14)
    write(100,1000) actual_state, jump_counter, start_pos, actual_pos, start_vel*mue, actual_vel*mue,&
               max_E_total-min_E_total, max_Norm-min_Norm

    return
  end subroutine print_traj_info



!---------- prints the actual position and velocity of particle 1 and 2
  subroutine printtraj()
    implicit none

    write(99,*) actual_time,sngl(actual_pos(:))
    write(98,*) actual_time,sngl(actual_vel(:))

  end subroutine printtraj



!---------- prints the actual coupling matrix P_mat
  subroutine print_Pmat()
    implicit none
    
123 format (25(F15.10,X))
    write(90,123) actual_time,P_mat(1,2,:),P_mat(1,3,:),P_mat(1,4,:),P_mat(2,3,:),P_mat(2,4,:),P_mat(3,4,:)

    return
  end subroutine print_Pmat



!---------- prints the status of the calculation
  subroutine printstat
    implicit none

    double precision :: kinetic_energy,kin_pot_energy

123 format (10(F15.10,X))

    kinetic_energy = 0.5d0*mue*dot_product(actual_vel,actual_vel)
    kin_pot_energy = kinetic_energy + V_mat(actual_state,actual_state)

    write(97,*) actual_time,jump_counter,actual_state
    write(96,123) actual_time,V_mat(1,1),V_mat(2,2),V_mat(3,3),V_mat(4,4),V_mat(actual_state,actual_state)
    write(95,*) actual_time,kin_pot_energy
    write(94,123) actual_time,real(a_mat(1,1)),aimag(a_mat(1,1)),real(a_mat(2,2)),aimag(a_mat(2,2)),real(a_mat(3,3)),aimag(a_mat(3,3)),real(a_mat(4,4)),aimag(a_mat(4,4))

  end subroutine printstat

!---------- print norm
  subroutine print_norm
    implicit none
    
    write(*,*) real(a_mat(1,1)+a_mat(2,2))

  end subroutine print_norm

!---------- print different quantities for debugging
  subroutine printstuff
    implicit none

124 format (5(F15.10,X))

    write(93,124) actual_time,g_mat(1,1),g_mat(1,2),g_mat(2,1),g_mat(2,2)
    write(92,124) actual_time,sngl(V_mat(1,1)),sngl(V_mat(1,2)),sngl(V_mat(2,1)),sngl(V_mat(2,2))
    write(91,124) actual_time,sngl(b_mat(1,1)),sngl(b_mat(1,2)),sngl(b_mat(2,1)),sngl(b_mat(2,2))

  end subroutine printstuff

!---------- random number generator
      FUNCTION ran2(idum)
      INTEGER :: idum
      INTEGER, parameter :: IM1=2147483563,IM2=2147483399,IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668
      INTEGER, parameter :: IQ2=52774,IR1=12211,IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB
      DOUBLE PRECISION ran2
      DOUBLE PRECISION, parameter :: AM=1.d0/IM1,EPS=1.2d-7,RNMX=1.d0-EPS
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
        end do
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END FUNCTION ran2



!---------- returns random number with gaussian distribution with a full width at half maximum of fwhm
      FUNCTION gasdev(idum,fwhm)
      INTEGER idum
      DOUBLE PRECISION gasdev
      INTEGER iset
      DOUBLE PRECISION fac,gset,rsq,v1,v2,fwhm
      SAVE iset,gset
      DATA iset/0/
!      if (iset.eq.0) then
1       v1=2.d0*ran2(idum)-1.d0
        v2=2.d0*ran2(idum)-1.d0
        rsq=v1**2+v2**2
        if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1
        fac=sqrt(-2.d0*log(rsq)/rsq)*fwhm/2.35482d0
        gset=v1*fac
        gasdev=v2*fac
        iset=1
!      else
!        gasdev=gset
!        iset=0
!      endif
      return
     END FUNCTION gasdev



!---------- performs a Runge-Kutta step for a NxN complex variable field depending on a real variable
      SUBROUTINE rk4_complex(y,n,x,h,yout,derivs)
      INTEGER n
      DOUBLE PRECISION h,x
      COMPLEX*16 dydx(n,n),y(n,n),yout(n,n)
      EXTERNAL derivs
      INTEGER i,j
      DOUBLE PRECISION h6,hh,xh
      COMPLEX*16 dym(n,n),dyt(n,n),yt(n,n)

      hh=h*0.5d0
      h6=h/6.d0
      xh=x+hh
      call derivs(x,y,dydx)
      do i=1,n
         do j=1,n
            yt(i,j)=y(i,j)+hh*dydx(i,j)
         end do
      end do
      call derivs(xh,yt,dyt)
      do i=1,n
         do j=1,n
            yt(i,j)=y(i,j)+hh*dyt(i,j)
         end do
      end do
      call derivs(xh,yt,dym)
      do i=1,n
         do j=1,n
            yt(i,j)=y(i,j)+h*dym(i,j)
            dym(i,j)=dyt(i,j)+dym(i,j)
         end do
      end do
      call derivs(x+h,yt,dyt)
      do i=1,n
         do j=1,n
            yout(i,j)=y(i,j)+h6*(dydx(i,j)+dyt(i,j)+2.d0*dym(i,j))
         end do
      end do
      return

     end SUBROUTINE rk4_complex

!---------- performs a Runge-Kutta step for a real variable field
      SUBROUTINE rk4(y,n,x,h,yout,derivs)
      INTEGER n
      DOUBLE PRECISION h,x,dydx(n),y(n),yout(n)
      EXTERNAL derivs
      INTEGER i
      DOUBLE PRECISION h6,hh,xh,dym(n),dyt(n),yt(n)

      hh=h*0.5d0
      h6=h/6.d0
      xh=x+hh
      call derivs(x,y,dydx)
      do i=1,n
        yt(i)=y(i)+hh*dydx(i)
      end do
      call derivs(xh,yt,dyt)
      do i=1,n
        yt(i)=y(i)+hh*dyt(i)
      end do
      call derivs(xh,yt,dym)
      do i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
      end do
      call derivs(x+h,yt,dyt)
      do i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))
      end do
      return

      end SUBROUTINE rk4

end module dynamik
