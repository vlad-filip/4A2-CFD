      
      subroutine calc_areas(g)

!     Calculate the area of the quadrilateral cells and the lengths of the side
!     facets

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_grid), intent(inout) :: g
      integer :: ni, nj

!     Declare integers or any extra variables you need here
!     INSERT
      real :: l_min, li_min, lj_min
      real, dimension(g%ni-1,g%nj-1) :: a1,a2,b1,b2

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Calculate the areas of the cells and store in g%area. The area of any
!     quadrilateral is half of the magnitude of the cross product of the two
!     vectors that form the diagonals. Check the order of your product so that
!     the values come out positive! You can do this using two nested loops in
!     the i and j-directions or in a vectorised way by indexing the coordinate
!     arrays with lists of indices
!     INSERT
      a1 = g%x(1:ni-1,1:nj-1)-g%x(2:ni,2:nj)
      a2 = g%y(1:ni-1,1:nj-1)-g%y(2:ni,2:nj)
      b1 = g%x(2:ni,1:nj-1)-g%x(1:ni-1,2:nj)
      b2 = g%y(2:ni,1:nj-1)-g%y(1:ni-1,2:nj)
      g%area = (a1*b2-a2*b1)/2

!     Calculate the projected lengths in the x and y-directions on all of the
!     "i = const" facets and store them in g%lx_i and g%ly_i. When combined
!     together these two components define a vector that is normal to the facet,
!     pointing inwards towards the centre of the cell. This is only the case for
!     the left hand side of the cell, the vector stored in position i,j points
!     towards the centre of the i,j cell
!     INSERT
      g%lx_i = -g%y(:,1:nj-1)+g%y(:,2:nj)
      g%ly_i = g%x(:,1:nj-1)-g%x(:,2:nj)
      

!     Now repeat the calculation for the project lengths on the "j=const"
!     facets. 
!     INSERT
      g%lx_j = g%y(1:ni-1,:)-g%y(2:ni,:)
      g%ly_j = -g%x(1:ni-1,:)+g%x(2:ni,:)

!     Find the minimum length scale in the mesh, this is defined as the length
!     of the shortest side of all the cells. Call this length "l_min", it is used
!     to set the timestep from the CFL number. Start by calculating the lengths
!     of the i and j facets by using the intrinsic function "hypot", this avoids
!     underflow and overflow errors. Then find the overal minimum value using
!     both the "min" and "minval" functions.
!     INSERT
      li_min = minval(hypot(g%lx_i,g%ly_i))
      lj_min = minval(hypot(g%lx_j,g%ly_j))
      g%l_min = min(li_min,lj_min)
      
!     Print the overall minimum length size that has been calculated
      write(6,*) 'Calculated cell areas and facet lengths'
      write(6,*) '  Overall minimum element size = ', g%l_min
      write(6,*)  

      end subroutine calc_areas
