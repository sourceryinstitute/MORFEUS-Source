!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
submodule(cartesian_grid_interface) cartesian_grid_implementation
  !! author: Damian Rouson and Karla Morris
  use kind_parameters, only : i4k, r8k
  use assertions_interface,only : assert, max_errmsg_len
  implicit none

contains

  module procedure div_scalar_flux

    integer(i4k) i, j, k, alloc_stat
    integer(i4k), parameter :: success=0
    real(r8k), parameter :: half=0.5_r8k
    character(len=max_errmsg_len) :: alloc_error
    real(r8k), allocatable, dimension(:,:,:) :: div_flux_x, div_flux_y, div_flux_z

    call assert( same_type_as(this, vertices), "div_scalar_flux: consistent types" )

    allocate(div_flux, mold=this, stat=alloc_stat, errmsg=alloc_error )
    call assert( alloc_stat==success, "div_scalar_flux (cartesian): result allocation fails with message '"//alloc_error//"'" )

    associate( positions => vertices%vectors(), s=>this%get_scalar() )
      associate( npoints => shape(positions(:,:,:,1)) )

        allocate(div_flux_x, div_flux_y, div_flux_z, mold=positions(:,:,:,1), stat=alloc_stat, errmsg=alloc_error )
        call assert( alloc_stat==success, "div_scalar_flux (cartesian): allocate(div_flux_{x,y,z}) (error: "//alloc_error//")" )

        associate( x=>positions(:,:,:,1) )
          do concurrent(k=1:npoints(3), j=1:npoints(2), i=2:npoints(1)-1)
            associate( &
              dx_m => half*(x(i+1,j,k) - x(i-1,j,k)), &
              dx_f =>       x(i+1,j,k) - x(i,j,k), &
              dx_b =>       x(i,j,k)   - x(i-1,j,k), &
              s_f => half*( s(i+1,j,k) + s(i,j,k)  ), &
              s_b => half*( s(i,j,k)   + s(i-1,j,k) ))
              associate( &
                D_f => this%diffusion_coefficient( s_f ), &
                D_b => this%diffusion_coefficient( s_b) )

                div_flux_x(i,j,k) = &
                  D_f*(s(i+1,j,k) - s(i,j,k)  )/(dx_f*dx_m) - &
                  D_b*(s(i,j,k)   - s(i-1,j,k))/(dx_b*dx_m)
              end associate
            end associate
          end do
        end associate

        hardwire_known_boundary_values: &
        block
          associate( y=>positions(:,:,:,2))
              div_flux_x(1,:,:) = 2*y(1,:,:)
              div_flux_x(npoints(1),:,:) = 2*y(npoints(1),:,:)
          end associate
        end block hardwire_known_boundary_values

        associate( y=>positions(:,:,:,2))
          do concurrent(k=1:npoints(3), j=2:npoints(2)-1, i=1:npoints(1))
            associate( &
              dy_m => (y(i,j+1,k) - y(i,j-1,k))*half, &
              dy_f =>  y(i,j+1,k) - y(i,j,k ), &
              dy_b =>  y(i,j,k)   - y(i,j-1,k), &
              s_f => half*(y(i,j+1,k) + y(i,j,k)), &
              s_b =>  half*(y(i,j,k) + y(i,j-1,k)) )

              associate( &
                D_f => this%diffusion_coefficient( s_f ), &
                D_b => this%diffusion_coefficient( s_b) )

                !div_flux_y(i,j,k) = &
                !  D_f*(s(i,j+1,k) - s(i,j,k))/(dy_f*dy_m) - &
                !  D_b*(s(i,j,k) - s(i,j-1,k))/(dy_b*dy_m)
             end associate
           end associate
         end do
       end associate

       associate( z=>positions(:,:,:,3))
         do concurrent(k=2:npoints(3)-1, j=1:npoints(2), i=1:npoints(1))
           associate( &
             dz_m => (z(i,j,k+1) - z(i,j,k-1))*half, &
             dz_f =>  z(i,j,k+1) - z(i,j,k), &
             dz_b =>  z(i,j,k)   - z(i,j,k-1), &
             s_f => half*(z(i,j,k+1)+z(i,j,k)), &
             s_b => half*(z(i,j,k)+z(i,j,k-1)) )
             associate( &
               D_f => this%diffusion_coefficient( s_f ), &
               D_b => this%diffusion_coefficient( s_b) )

               !div_flux_z(i,j,k) = &
               !  D_f*(s(i,j,k+1) - s(i,j,k))/(dz_f*dz_m) - &
               !  D_b*(s(i,j,k) - s(i,j,k-1))/(dz_b*dz_m)
             end associate
           end associate
         end do
       end associate

       call div_flux%set_scalar( div_flux_x )
       !call div_flux%set_scalar( div_flux_x + div_flux_y + div_flux_z )

      end associate
    end associate
  end procedure div_scalar_flux

end submodule
