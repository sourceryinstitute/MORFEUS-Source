!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
program main
  use assertions_interface, only : assert, max_errmsg_len
  use problem_discretization_interface, only :  problem_discretization
  use cartesian_grid_interface, only : cartesian_grid
  use inbox_interface, only : inbox, face_index
  implicit none

  type(inbox), allocatable, dimension(:,:), codimension[:] :: halo
    !! container for data from neighboring blocks: dimension=(my_blocks_padded, size(face_index)), codimension=(num_images)
  type(problem_discretization) block_structured_grid
    !! encapsulate the global grid structure
  type(cartesian_grid) prototype
    !! pass the cartesian_grid type
  integer, parameter :: num_structured_grids(*) = [3,3,3]
    !! number of subdomains in each coordinate direction
  integer, parameter :: lower=1, upper=2, success=0, nsteps=1
    !! my_blocks bounds, number of halo interfaces, allocation status success value, number of time steps
  integer block_, face, alloc_stat, step
    !! block loop counter, allocatiaon status result, time step number
  character(len=max_errmsg_len) error_message
    !! allocation error message

  associate( me => this_image(), ni => num_images(), num_blocks => product(num_structured_grids) )

    call assert( ni<=num_blocks, "test-halo-exchange: enough blocks to distribute to images")

    call block_structured_grid%partition( num_structured_grids, prototype )
      !! partition the block-structured grid into subdomains with connectivity implied by the supplied shape array

    associate( block_bounds => block_structured_grid%my_subdomains() )
      associate( my_blocks => [(block_, block_ = block_bounds(lower), block_bounds(upper))] )
        associate( my_blocks_padded => size(my_blocks) + merge(0,1,mod(num_blocks,ni)==0) )
           !! pad to ensure standard-conforming coarray allocation: all images must allocate same bounds & cobounds

          allocate( halo(my_blocks_padded, size(face_index))[*], stat=alloc_stat, errmsg=error_message )
          call assert(alloc_stat==success, "test-halo-exchange: halo allocation", error_message)

          loop_over_time_steps: &
          do step=1,nsteps

              loop_over_blocks: &
              do block_ = block_bounds(lower), block_bounds(upper)

                loop_over_neighbors: &
                do face = lbound(face_index,1), ubound(face_index,1)

                 !associate( neighbor_id => block_structured_grid%my_neighbor(block_, face_index(face)) )
                 !  associate( neighbor_image => block_structured_grid%neighbor_image(neighbor_id) )
                 !     associate( my_blocks_array_b_index=> findloc(block_, my_blocks_array) )
                 !       if (step>1) call halo(my_blocks_array_b_index, face_index(face))%event_wait() ! atomically decrement counter
                 !     end associate
                 !   end associate
                 ! end associate

                end do loop_over_neighbors
              end do loop_over_blocks

          end do loop_over_time_steps

        end associate
      end associate
    end associate
  end associate

  print *,"Test passed"
end program

       ! call halo(,)[]%set_sender_block_id(b)
       ! call halo(,)[]%set_step(step)

       ! write(greeting,*) "Hello from image ",me," of ",ni," on step ",step
       ! ! Signal image 1 that a new greeting is ready for pickup:
       ! event post(halo%message_ready(me)[1])              ! Atomically increments my event greeting-ready counter on image 1

       ! write(greeting,*) "Hello from image ",me," of ",ni," on step ",step
       ! print *,greeting

       ! spin_query_work: block
       !   integer :: image,ready_count
       !   logical, dimension(2:ni) :: message_not_read

       !   message_not_read=.true.

       !   spin: do while( any( message_not_read  ) ) ! Loop until all greetings have been printed
       !     query: do image=2,ni             ! Atomically access each event's counter
       !       if (message_not_read(image)) then      ! Print greetings that have not been printed during this step
       !         call event_query( message_ready(image), ready_count)
       !         work_if_ready: select case(ready_count)
       !           case(0) ! keep spinning until greeting is ready
       !           case(1) ! event posted so get and print greeting
       !             event wait(message_ready(image))
       !             print *,greeting[image]
       !             event post(message_read[image])
       !             message_not_read(image)=.false.
       !           case default
       !             if (ready_count<0) error stop "compiler bug: negative event_query count"
       !             error stop "multiple events happened since the last event query"
       !         end select work_if_ready
       !       end if
       !     end do query
       !   end do spin

       ! end block spin_query_work
