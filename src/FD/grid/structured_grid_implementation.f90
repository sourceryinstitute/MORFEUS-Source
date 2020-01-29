!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(structured_grid_interface) structured_grid_implementation
  use assertions_interface, only : assert,assertions
  implicit none

  integer, dimension(:), allocatable :: global_block_shape

contains

    module procedure set_global_block_shape
      global_block_shape = shape_array
    end procedure

    module procedure get_global_block_shape
      shape_array = global_block_shape
    end procedure

    module procedure clone
      this%nodal_values = original%nodal_values
      this%global_bounds = original%global_bounds
      this%metadata = original%metadata
    end procedure

    module procedure num_cells
      associate( grid_shape => shape(this%nodal_values) )
        cell_count = (grid_shape(1)-1) * (grid_shape(2)-1) * (grid_shape(3)-1)
        !! in each coordinatee direction, there are one fewer cells than grid points
      end associate
    end procedure

    module procedure get_scalar
      call assert(this%space_dimension()==3 .and. this%free_tensor_indices()==0, "structured_grid%get_scalar: scalar on 3D grid")
      call assert(this%num_time_stamps()==1, "structured_grid%get_scalar: single snapshot of scalar data at 3D grid vertices")
      scalar_values = this%nodal_values(:,:,:,1,1,1)
    end procedure

    module procedure vectors
      call assert(this%space_dimension()==3 .and. this%free_tensor_indices()==1 .and. this%num_time_stamps()==1, &
        "structured_grid%vectors: single snapshot of 3D vector data")
      vectors3D = this%nodal_values(:,:,:,:,1,1)
    end procedure

    module procedure diffusion_coefficient
      coefficient = 1._r8k
      !coefficient = this%metadata%diffusion_coefficent
    end procedure

    module procedure write_formatted
      integer xloc, yloc, zloc

      associate( num_tensor_components => size(this%nodal_values,4)*size(this%nodal_values,5) )
        associate( num_time_stamps => size(this%nodal_values,6) )
          call assert(num_tensor_components==3, "structured_grid%write_formatted: vector data")
          call assert(num_time_stamps==1,"structured_grid%write_formatted: single time instant")
      end associate; end associate

      do xloc=1,lbound(this%nodal_values,1),ubound(this%nodal_values,1)
        do yloc=1,lbound(this%nodal_values,2),ubound(this%nodal_values,2)
          do zloc=1,lbound(this%nodal_values,3),ubound(this%nodal_values,3)
            write(unit,'(3(e12.4,","),3x,a)') this%nodal_values(xloc,yloc,zloc,:,1,1), "1"
            write(unit,*) new_line('a')
      end do; end do; end do

    end procedure write_formatted

    module procedure space_dimension
      integer, parameter :: max_space_dimension=3
      logical finite_extents(max_space_dimension)

      ! Requires
      call assert(allocated(this%nodal_values),"structured_grid%space_dimension: allocated(this%nodal_values)")

      finite_extents = [size(this%nodal_values,1)>1,size(this%nodal_values,2)>1,size(this%nodal_values,3)>1]
      num_dimensions = count(finite_extents)

    end procedure

    module procedure free_tensor_indices
      integer, parameter :: max_free_indices=2
      logical free_indices(max_free_indices)

      ! Requires
      call assert(allocated(this%nodal_values),"nodal values allocated")

      free_indices = [size(this%nodal_values,4)>1,size(this%nodal_values,5)>1]
      num_free_indices = count(free_indices)
    end procedure

    module procedure num_time_stamps
      call assert(allocated(this%nodal_values),"nodal values allocated")
      num_times = size(this%nodal_values,6)
    end procedure

    module procedure get_tag
      this_tag = this%metadata%get_tag()
    end procedure

    module procedure set_metadata
      this%metadata = metadata
    end procedure

    module procedure set_vector_components

      integer alloc_status
      integer, parameter :: components=3, time_stamps=1, dyad_components=1

      associate( x_shape=>shape(x_nodes),  y_shape=>shape(y_nodes), z_shape=>shape(z_nodes) )

      ! Requires
      if (assertions) &
        call assert( all(x_shape==y_shape) .and. all(y_shape==z_shape), "set_vector_components: inconsistent nodal values arrays" )

        allocate( this%nodal_values(x_shape(1),x_shape(2),x_shape(3),components,dyad_components,time_stamps), stat=alloc_status )
        call assert( alloc_status==0, "set_vector_components: allocation successful" )

        this%nodal_values(:,:,:,1,dyad_components,time_stamps) = x_nodes(:,:,:)
        this%nodal_values(:,:,:,2,dyad_components,time_stamps) = y_nodes(:,:,:)
        this%nodal_values(:,:,:,3,dyad_components,time_stamps)= z_nodes(:,:,:)

      end associate

    end procedure

    module procedure set_scalar

      integer alloc_status
      integer, parameter :: time_stamps=1

      if (allocated(this%nodal_values)) deallocate(this%nodal_values)
      allocate( this%nodal_values(size(scalar,1),size(scalar,2),size(scalar,3),1,1,time_stamps), stat=alloc_status )
      call assert( alloc_status==0, "set_scalar: allocation successful" )

      this%nodal_values(:,:,:,1,1,time_stamps) = scalar

    end procedure

    module procedure increment_scalar

      integer, parameter :: time_stamps=1

      this%nodal_values(:,:,:,1,1,time_stamps) = this%nodal_values(:,:,:,1,1,time_stamps) + scalar

    end procedure

    module procedure subtract
      error stop "structured_grid%subtract: this implementation generates a runtime invalid memory reference with GCC 8.3"
      call assert (allocated(this%nodal_values) .and. allocated(rhs%nodal_values), "structured_grid%subtract: operands allocated")
      call assert (shape(this%nodal_values) == shape(rhs%nodal_values), "structured_grid%subtract: operands conform")
      call assert (same_type_as(this, rhs), "structured_grid%subtract: operand types match")
      difference%nodal_values = this%nodal_values - rhs%nodal_values
    end procedure

    module procedure compare
      call assert (allocated(this%nodal_values).and.allocated(reference%nodal_values),"structured_grid%compare: operands allocated")
      call assert (shape(this%nodal_values) == shape(reference%nodal_values), "structured_grid%subtract: operands conform")
      call assert (same_type_as(this, reference), "structured_grid%subtract: operand types match")
      associate(L_infinity_norm => maxval(abs(this%nodal_values - reference%nodal_values)))
        call assert(L_infinity_norm <= tolerance, "structured_grid%compare: L_infinity_norm <= tolerance")
      end associate
    end procedure

end submodule
