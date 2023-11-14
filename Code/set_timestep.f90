      
      subroutine set_timestep(av,g,bcs)

!     This subroutine sets a single value for the time step based on the 
!     stagnation speed of sound and the minimum length scale of any element

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(in) :: g
      type(t_bconds), intent(in) :: bcs
      real :: astag, v_max

!     Calculate the stagnation speed of sound from the inlet stagnation
!     temperature and gas constants
!     INSERT
      astag = sqrt(av%gam*av%rgas*bcs%tstag)

!     Assume that the maximum flow speed is also equal to "astag". This will 
!     be pessimistic for subsonic flows but may be optimistic for supersonic 
!     flows. In the latter case the length of the time step as determined by 
!     may need to be reduced by improving this routine or varying the CFL number
!     INSERT
      v_max = astag

!     Calculate the timestep using the CFL number and store it in "av%dt"
!     INSERT
      av%dt = av%cfl*g%l_min/(v_max+astag)

!     Print the calculated timestep and some intermediate values
!     INSERT
      write(6,*) 'Time step Calculations'
      write(6,*) 'dt =' , av%dt, '      astag =', astag, '      v_max=', v_max
      write(6,*)

      end subroutine set_timestep


