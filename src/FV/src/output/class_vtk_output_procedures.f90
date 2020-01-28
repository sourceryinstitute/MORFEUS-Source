!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
SUBMODULE (class_vtk_output) class_vtk_output_procedures
    USE iso_fortran_env, ONLY : r8k => real64, output_unit
    USE class_psblas,    ONLY : psb_dpk_
    USE class_mesh,      ONLY : mesh
    IMPLICIT NONE
    !! author: Ian Porter, NRC
    !! date: 01/23/2019
    !!
    !! This submodule implements the routines necessary to interface morfeus w/ vtkmofo
    !!
    !! Always write cell ID & processor ID
    CHARACTER(LEN=10), DIMENSION(*), PARAMETER :: fixed_cell_data  = &
        & [ 'cellIDs   ', 'procIDs   ' ]
    INTEGER, PARAMETER :: n_base_cell_values = SIZE(fixed_cell_data,DIM=1)

CONTAINS

    MODULE PROCEDURE write_vtk_morfeus
        USE class_psblas,   ONLY : mypnum_, nprocs_
        USE vtk_attributes, ONLY : attributes, scalar, vector
        USE vtk_cells,      ONLY : vtkcell, vtkcell_list, set_cell_type
        USE vtk_datasets,   ONLY : unstruct_grid
        USE vtk,            ONLY : vtk_serial_write
        IMPLICIT NONE
        INTEGER :: i, j, ibegin, itype
        INTEGER, DIMENSION(:), ALLOCATABLE :: cell_ids
        INTEGER, DIMENSION(:), ALLOCATABLE :: v2cconn
        INTEGER, DIMENSION(:), ALLOCATABLE :: icverts
        INTEGER, DIMENSION(:), ALLOCATABLE :: iproc
        REAL(psb_dpk_), DIMENSION(:),     ALLOCATABLE :: scalar_val
        REAL(psb_dpk_), DIMENSION(:,:),   ALLOCATABLE :: points
        REAL(psb_dpk_), DIMENSION(:,:,:), ALLOCATABLE :: vector_vals
        CLASS(vtkcell), ALLOCATABLE :: dummy
        TYPE(unstruct_grid)         :: vtk_mesh
        TYPE(attributes),   DIMENSION(:), ALLOCATABLE :: cell_vals_to_write
        TYPE(attributes),   DIMENSION(:), ALLOCATABLE :: point_vals_to_write
        TYPE(vtkcell_list), DIMENSION(:), ALLOCATABLE :: morfeus_cells

        ! Sets output path (bypass b/c vtkmofo handles this)
        !IF(PRESENT(iter)) CALL out%set_output_path(iter)

        CALL write_vtk_mesh (msh, points, cell_ids, v2cconn, icverts, iproc)

!WRITE(output_unit,*) 'Before allocate field_vals in write_vtk_morfeus on image: ',mypnum_()
!write(output_unit,*) size(points), size(cell_ids), size(v2cconn), size(icverts)

        ALLOCATE(vector_vals(1:SIZE(vfield),1:SIZE(cell_ids),1:3), source=0.0_psb_dpk_)
        DO i = 1, SIZE(vfield)
            !! This collapses all images to one
!            vector_vals(i,1:3,:) = out%get_vector_field(vfield(i))
            !! This keeps things local
            BLOCK
                USE class_vector, ONLY : vector, ASSIGNMENT(=)
                TYPE(vector), DIMENSION(:), ALLOCATABLE :: x_loc
                ASSOCIATE (my_field => vfield(i))
                    CALL my_field%get_x(x_loc)
                    DO j = 1, SIZE(x_loc)
                        vector_vals(i,j,1:3) = x_loc(j)
                    END DO
                END ASSOCIATE
            END BLOCK
        END DO

        ALLOCATE(morfeus_cells(1:SIZE(cell_ids))); ibegin = 1
        DO i = 1, SIZE(morfeus_cells)
            SELECT CASE (icverts(i))
                !! Take the # of vertices and convert it to the type of cell
            CASE (3)
                itype = 5         !! Triangle
            CASE (4)
                IF (msh%ncd>2) THEN
                    itype = 10    !! Tetra
                ELSE
                    itype = 9     !! Quad
                END IF
            CASE (5)
                itype = 14        !! Pyramid
            CASE (6)
                itype = 13        !! Wedge
            CASE (8)
                itype = 12        !! Hexahedron
            CASE DEFAULT
                ERROR STOP 'Unsupported # of vertices'
            END SELECT

            ASSOCIATE (connectivity => v2cconn(ibegin:ibegin+icverts(i)-1))
                dummy = set_cell_type(itype)                              !! Sets the proper type of cell based on value of icverts
                CALL dummy%setup(points=connectivity-1)                   !! Need to fill in an array of connectivity
                ALLOCATE(morfeus_cells(i)%cell,source=dummy)              !! Set cell into array of cells
                DEALLOCATE(dummy)
            END ASSOCIATE

            ibegin = ibegin + icverts(i)
        END DO

        CALL vtk_mesh%init (points=points, cell_list=morfeus_cells, datatype='double')
                               !! Set up the geometry using a vtk structured grid
!        CALL vtk_parallel_write (geometry=vtk_mesh, image=mypnum_(), filename=TRIM(out%path_()), multiple_io=.TRUE.) !!
        CALL vtk_serial_write (geometry=vtk_mesh, filename=TRIM(out%path_()), multiple_io=.TRUE.) !!

        ASSOCIATE (n_scalars => SIZE(fixed_cell_data) + SIZE(sfield), &
            &      n_vectors => SIZE(vfield))
            ALLOCATE(cell_vals_to_write(1:n_scalars+n_vectors))
            DO i = 1, SIZE(cell_vals_to_write)
                IF (i <= SIZE(fixed_cell_data)) THEN
                    !! Set cell values
                    IF (.NOT. ALLOCATED(cell_vals_to_write(i)%attribute)) THEN
                        ALLOCATE(scalar::cell_vals_to_write(i)%attribute)
                    END IF
                    SELECT CASE (i)
                    CASE (1)
                        CALL cell_vals_to_write(i)%attribute%init (TRIM(fixed_cell_data(i)), numcomp=1, int1d=cell_ids)
                    CASE (2)
                        CALL cell_vals_to_write(i)%attribute%init (TRIM(fixed_cell_data(i)), numcomp=1, int1d=iproc)
                    END SELECT
                ELSE IF (i <= n_scalars) THEN
                    !! Program supplied scalars
                    IF (.NOT. ALLOCATED(cell_vals_to_write(i)%attribute)) THEN
                        ALLOCATE(scalar::cell_vals_to_write(i)%attribute)
                    END IF
                    ASSOCIATE (my_id => i-SIZE(fixed_cell_data))
                        IF (ALLOCATED(scalar_val)) DEALLOCATE(scalar_val)
                        CALL sfield(my_id)%get_x(scalar_val)
                        print*,size(scalar_val)
                        CALL cell_vals_to_write(i)%attribute%init (sfield(my_id)%name_(), numcomp=1, real1d=scalar_val)
                    END ASSOCIATE
                ELSE
                    !! Program supplied vectors
                    ASSOCIATE (cnt => i - n_scalars)
                        !! Set cell values
                        IF (.NOT. ALLOCATED(cell_vals_to_write(i)%attribute)) THEN
                            ALLOCATE(vector::cell_vals_to_write(i)%attribute)
                        END IF
                        !! Supplied values
                        CALL cell_vals_to_write(i)%attribute%init (vfield(cnt)%name_(), numcomp=1, real2d=vector_vals(cnt,:,:))
                    END ASSOCIATE
                END IF
            END DO
        END ASSOCIATE

        CALL vtk_serial_write (celldatasets=cell_vals_to_write)  ! pointdatasets=point_vals_to_write

        CALL vtk_serial_write (finished=.TRUE.)                 !! Finish writing the serial file for each image

!        IF (mypnum_() == 0) CALL vtk_parallel_write(nprocs_())    !! Only do this on image 0

    END PROCEDURE write_vtk_morfeus

    SUBROUTINE write_vtk_mesh(msh, points, cell_ids, v2cconn, icverts, iproc)
        USE class_psblas
        USE class_cell
        USE class_connectivity
        USE class_face
        USE class_iterating
        USE class_mesh
        USE class_output
        USE class_vertex
        USE tools_output_basics
        IMPLICIT NONE
        !
        TYPE(mesh),      INTENT(IN) :: msh
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: cell_ids
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: v2cconn
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: icverts
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: iproc
        REAL(psb_dpk_), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: points
        !
        INTEGER :: info, err_act
        INTEGER :: icontxt, mypnum
        INTEGER :: i, ic, ig, ncells, ngc, ngroups
        INTEGER, POINTER :: ic2g(:) => NULL()
        INTEGER, ALLOCATABLE  :: igroup(:)
        INTEGER, ALLOCATABLE :: i_loc(:)
        TYPE(cell), ALLOCATABLE :: cells(:)
!        TYPE(face), ALLOCATABLE :: faces(:)
        TYPE(vertex), ALLOCATABLE :: verts(:)
        TYPE(connectivity) :: v2c!, f2c, v2f

        ! Sets error handling for PSBLAS-2 routines
        CALL psb_erractionsave(err_act)

        mypnum  = mypnum_()
        icontxt = icontxt_()

        ! Global number of cells
        ncells = psb_cd_get_global_cols(msh%desc_c)

        CALL psb_geall(i_loc,msh%desc_c,info)
        CALL psb_check_error(info,'write_mesh','psb_geall',icontxt)

        ALLOCATE(igroup(ncells),iproc(ncells),stat=info)
        IF(info /= 0) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Gathers MESH components on P0
!!!        CALL l2g_vertex(msh%verts,verts,msh%desc_v)
!        CALL l2g_face(msh%faces,faces,msh%desc_f,msh%desc_c)
!!!        CALL l2g_cell(msh%cells,cells,msh%desc_c)
!        CALL l2g_conn(msh%v2f,v2f,msh%desc_v,msh%desc_f)
!!!        CALL l2g_conn(msh%v2c,v2c,msh%desc_v,msh%desc_c)
!        CALL l2g_conn(msh%f2c,f2c,msh%desc_f,msh%desc_c)
         verts = msh%verts
         cells = msh%cells
         v2c = msh%v2c

        ! Gathers processor IDs
        i_loc(:) = mypnum
        CALL psb_gather(iproc,i_loc,msh%desc_c,info,root=0)
        CALL psb_check_error(info,'write_mesh','psb_gather',icontxt)

        ! Gathers group IDs
        i_loc(:) = 0
        ngroups = msh%c2g%nel_()

        DO ig = 1, ngroups
            CALL msh%c2g%get_ith_conn(ic2g,ig)
            ngc = SIZE(ic2g) ! number of group cells
            DO i = 1, ngc
                ic = ic2g(i)
                i_loc(ic) = ig
            END DO
        END DO

        CALL psb_gather(igroup,i_loc,msh%desc_c,info,root=0)
        CALL psb_check_error(info,'write_mesh','psb_gather',icontxt)

        CALL wr_vtk_mesh(msh%ncd, verts, cells, v2c, points, cell_ids, v2cconn, icverts)

        ! Frees Memory
        NULLIFY(ic2g)
        DEALLOCATE(igroup)

        CALL psb_gefree(i_loc,msh%desc_c,info)
        CALL psb_check_error(info,'write_mesh','psb_gefree',icontxt)

!        CALL free_conn(f2c)
        CALL free_conn(v2c)
!        CALL free_conn(v2f)
        IF(ALLOCATED(cells)) CALL free_cell(cells)
!        IF(ALLOCATED(faces)) CALL free_face(faces)
        IF(ALLOCATED(verts)) CALL free_vertex(verts)

        ! ----- Normal termination -----
        CALL psb_erractionrestore(err_act)

100     FORMAT(' ERROR! Memory allocation failure in WRITE_MESH')
!200     FORMAT(' ERROR! Unsupported output format in WRITE_MESH')

    END SUBROUTINE write_vtk_mesh

    SUBROUTINE wr_vtk_mesh(ncd, verts, cells, v2c, points, cell_ids, v2cconn, icverts)
        USE class_psblas
        USE class_cell
        USE class_connectivity
        USE class_vertex
        IMPLICIT NONE

        ! Input parameters
        INTEGER,            INTENT(IN) :: ncd
        TYPE(vertex),       INTENT(IN) :: verts(:)
        TYPE(cell),         INTENT(IN) :: cells(:)
        TYPE(connectivity), INTENT(IN) :: v2c
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: cell_ids
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: v2cconn
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: icverts
        REAL(psb_dpk_), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: points

        ! Local variables
        INTEGER :: i, info
        INTEGER, ALLOCATABLE :: v2clookup(:)             ! CSR data for v2c connectivity
        INTEGER :: v2cnconn                              ! number of connections to faces
        INTEGER :: nverts, ncells                        ! number of vertices & cells
        REAL(psb_dpk_), ALLOCATABLE :: xpos(:), ypos(:), zpos(:)

        nverts = SIZE(verts)
        ncells = SIZE(cells)

!        IF(mypnum_() /= 0) RETURN ! this need only be done by one processor

        ! pack up the data into arrays and pass it off to a C function to handle the IO
        ALLOCATE(xpos(nverts),ypos(nverts),zpos(nverts),icverts(ncells),stat=info)
        IF(info /= 0) THEN
            WRITE(*,300)
            CALL abort_psblas
        END IF

        ! use the getters to pack position arrays
        xpos = verts(:)%x_()
        ypos = verts(:)%y_()
        zpos = verts(:)%z_()

        ALLOCATE(points(3,1:SIZE(xpos)),source=0.0_r8k)
        points(1,:) = xpos; points(2,:) = ypos; points(3,:) = zpos  !! Set x,y,z positions

        ! VTK format needs v2c connectivity
        CALL v2c%get_conn_csr(v2clookup,v2cconn)
        v2cnconn = SIZE(v2cconn)

        icverts(:) = cells%nv_()

        ! passing:
        ! (0) nverts
        ! (1) vertex positions
        ! (2) ncells
        ! (3) size of v2c array (see next)
        ! (4) connectivity array for v2c
        ! (5) array listing number of vertices for each cell

        ! Set the cells IDs
        ALLOCATE(cell_ids(ncells),stat=info)
        IF(info /= 0) THEN
            WRITE(*,300)
            CALL abort_psblas
        ELSE
            cell_ids=[ (i,i=1,ncells) ]
        END IF

300     FORMAT(' ERROR! Memory allocation failure in WR_VTK_MESH')

    END SUBROUTINE wr_vtk_mesh

END SUBMODULE class_vtk_output_procedures
