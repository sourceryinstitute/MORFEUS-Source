!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
SUBMODULE (vtkmofo_io) vtkmofo_io_implementation
    USE iso_fortran_env, ONLY : r8k => real64
    USE class_psblas,    ONLY : psb_dpk_
    USE class_mesh, ONLY : mesh
    IMPLICIT NONE
    !! author: Ian Porter, NRC
    !! date: 01/23/2019
    !!
    !! This submodule implements the routines necessary to interface morfeus w/ vtkmofo
    !!
    CHARACTER(len=*), PARAMETER :: filename = 'heat_pipe.vtk'
    CHARACTER(len=*), PARAMETER :: title = 'heatpipe example'
    INTEGER, PARAMETER          :: unit = 10
    INTEGER, PARAMETER          :: n_cell_values     = 3
    INTEGER, PARAMETER          :: n_point_values    = 1
    INTEGER, PARAMETER          :: n_params_to_write = n_cell_values + n_point_values
    CHARACTER(LEN=10), DIMENSION(n_cell_values-1),  PARAMETER :: cell_dataname  = &
        & [ 'cellIDs   ', 'procIDs   ' ]
    CHARACTER(LEN=15), DIMENSION(n_point_values), PARAMETER :: point_dataname = &
        & [ 'Temperature(K) ' ]

    CONTAINS

    MODULE PROCEDURE write_vtk_morfeus
        IMPLICIT NONE
        INTEGER :: i
        CHARACTER(LEN=100) :: file_name
        INTEGER, DIMENSION(:), ALLOCATABLE :: cell_ids
        INTEGER, DIMENSION(:), ALLOCATABLE :: v2cconn
        INTEGER, DIMENSION(:), ALLOCATABLE :: icverts
        INTEGER, DIMENSION(:), ALLOCATABLE :: iproc
        REAL(psb_dpk_), DIMENSION(:,:), ALLOCATABLE :: points
        REAL(psb_dpk_), DIMENSION(:,:), ALLOCATABLE :: field_vals

        CALL write_vtk_mesh (msh, points, cell_ids, v2cconn, icverts, iproc, out, iter, file_name)

        ALLOCATE (field_vals(1:SIZE(sfield), 1:SIZE(cell_ids)))
        DO i = 1, SIZE(sfield)
            field_vals(i,:) = set_scalar_field(sfield(i), sname(i),out,iter)
        END DO

        !        DO i = 1, SIZE(vfield)
        !            call write_vector_field(vfield(i),vname(i),out,iter)
        !        END DO

        CALL write_vtkmofo (msh%ncd, points, cell_ids, v2cconn, icverts, iproc, field_vals(1,:), sname, TRIM(file_name))

  END PROCEDURE write_vtk_morfeus

  SUBROUTINE write_vtk_mesh(msh, points, cell_ids, v2cconn, icverts, iproc, out, iter, filename)
      USE class_psblas
      USE class_cell
      USE class_connectivity
      USE class_face
      USE class_iterating
      USE class_mesh
      USE class_output
      USE class_vertex
      !
      USE tools_output_basics

      IMPLICIT NONE
      !
      TYPE(mesh),      INTENT(IN) :: msh
      TYPE(output),    INTENT(INOUT) :: out
      TYPE(iterating), INTENT(IN), OPTIONAL :: iter
      CHARACTER(len=*), INTENT(OUT) :: filename
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
      TYPE(face), ALLOCATABLE :: faces(:)
      TYPE(vertex), ALLOCATABLE :: verts(:)
      TYPE(connectivity) :: v2f, v2c, f2c
      CHARACTER(len=32) :: path

      ! Sets error handling for PSBLAS-2 routines
      CALL psb_erractionsave(err_act)

      mypnum  = mypnum_()
      icontxt = icontxt_()

      ! Sets output path
      IF(PRESENT(iter)) CALL out%set_output_path(iter)
      path = out%path_()
      filename = TRIM(path)
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
      CALL l2g_vertex(msh%verts,verts,msh%desc_v)
      CALL l2g_face(msh%faces,faces,msh%desc_f,msh%desc_c)
      CALL l2g_cell(msh%cells,cells,msh%desc_c)
      CALL l2g_conn(msh%v2f,v2f,msh%desc_v,msh%desc_f)
      CALL l2g_conn(msh%v2c,v2c,msh%desc_v,msh%desc_c)
      CALL l2g_conn(msh%f2c,f2c,msh%desc_f,msh%desc_c)

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

      IF(mypnum == 0) THEN
          CALL wr_vtk_mesh(msh%ncd,verts,cells,v2c, points, cell_ids, v2cconn, icverts)
      END IF

      ! Frees Memory
      NULLIFY(ic2g)
      DEALLOCATE(igroup)

      CALL psb_gefree(i_loc,msh%desc_c,info)
      CALL psb_check_error(info,'write_mesh','psb_gefree',icontxt)

      CALL free_conn(f2c)
      CALL free_conn(v2c)
      CALL free_conn(v2f)
      IF(ALLOCATED(cells)) CALL free_cell(cells)
      IF(ALLOCATED(faces)) CALL free_face(faces)
      IF(ALLOCATED(verts)) CALL free_vertex(verts)


      ! ----- Normal termination -----
      CALL psb_erractionrestore(err_act)

100   FORMAT(' ERROR! Memory allocation failure in WRITE_MESH')
200   FORMAT(' ERROR! Unsupported output format in WRITE_MESH')

  END SUBROUTINE write_vtk_mesh

  SUBROUTINE wr_vtk_mesh(ncd,verts,cells,v2c, points, cell_ids, v2cconn, icverts)
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

      !  character(len=40) :: file_name

      nverts = SIZE(verts)
      ncells = SIZE(cells)

      IF(mypnum_() /= 0) RETURN ! this need only be done by one processor

      ! remember to add null to the end of the string for C compatibility
      !  file_name = trim(path)//char(0)

      ! make sure that our int variable is the size that C expects
      IF(bit_SIZE(info) /= 32) THEN
          WRITE(6,200)
          CALL abort_psblas
      ENDIF

      ! pack up the data into arrays and pass it off to a C function to handle
      ! the IO

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


      ! VTK_MESH is a C function that does the actual writing
      ! returns an integer value where 0 indicates success in I/O and
      ! a non-zero value indicates an error.

      ! Write cells IDs
      ALLOCATE(cell_ids(ncells),stat=info)
      IF(info /= 0) THEN
          WRITE(*,300)
          CALL abort_psblas
      END IF

      cell_ids=(/ (i,i=1,ncells) /)

100   FORMAT(' ERROR! VTK_MESH currently only supports 3D output')
200   FORMAT(' ERROR! Wrong size for integer variables in WR_VTK_MESH')
300   FORMAT(' ERROR! Memory allocation failure in WR_VTK_MESH')
400   FORMAT(' Error in writing VTK mesh. Code =', i3)
  END SUBROUTINE wr_vtk_mesh

  FUNCTION set_scalar_field(fld,field,out,iter) RESULT (x_glob)
      USE class_psblas
      USE class_cell
      USE class_iterating
      USE class_mesh
      USE class_output
      USE class_scalar_field
      !
      USE tools_output_basics

      IMPLICIT NONE
      !
      TYPE(scalar_field), INTENT(IN) :: fld
      CHARACTER(len=*),   INTENT(IN) :: field
      TYPE(output),       INTENT(INOUT) :: out
      TYPE(iterating),    INTENT(IN), OPTIONAL :: iter
      !
      INTEGER :: err_act, icontxt, info, mypnum, ncells_glob
      REAL(psb_dpk_), ALLOCATABLE :: x_loc(:)
      REAL(psb_dpk_), ALLOCATABLE :: x_glob(:)
      CHARACTER(len=32) :: path
      TYPE(cell), ALLOCATABLE :: cells_glob(:)
      TYPE(mesh), POINTER :: msh => NULL()


      ! Sets error handling for PSBLAS-2 routines
      CALL psb_erractionsave(err_act)

      mypnum  = mypnum_()
      icontxt = icontxt_()

      ! Sets output path
      IF(PRESENT(iter)) CALL out%set_output_path(iter)
      path = out%path_()

!!$  msh   => msh_(fld)
      CALL fld%get_mesh(msh)
      CALL fld%get_x(x_loc)

      ! Is FLD cell-centered?
      IF(fld%on_faces_()) THEN
          WRITE(*,100)
          CALL abort_psblas
      END IF

      ! Global number of cells
      ncells_glob = psb_cd_get_global_cols(msh%desc_c)

      ALLOCATE(x_glob(ncells_glob),stat=info)
      IF(info /= 0) THEN
          WRITE(*,200)
          CALL abort_psblas
      END IF

      ! Gathers cell-centered values
      CALL psb_gather(x_glob,x_loc,msh%desc_c,info,root=0)
      CALL psb_check_error(info,'set_scalar_field','psb_gather',icontxt)

      ! Local-to-global reallocation of CELL object
      CALL l2g_cell(msh%cells,cells_glob,msh%desc_c)

      IF(mypnum == 0) THEN
          SELECT CASE(out%fmt_())
          CASE(vtk_)
              ! call wr_vtk_field(field,x_glob,msh%ncd,path)
              ! Only need value for x_glob
          CASE default
              WRITE(*,300)
              CALL abort_psblas
          END SELECT
          WRITE(*,'()')
      END IF


      IF(ALLOCATED(cells_glob)) CALL free_cell(cells_glob)
      DEALLOCATE(x_loc)
      NULLIFY(msh)

      ! ----- Normal Termination -----
      CALL psb_erractionrestore(err_act)

100   FORMAT(' ERROR! Face-centered field in set_scalar_field')
200   FORMAT(' ERROR! Memory allocation failure in set_scalar_field')
300   FORMAT(' ERROR! Unsupported output format in set_scalar_field')

  END FUNCTION set_scalar_field


  SUBROUTINE set_vector_field(fld,field,out,iter)
      USE class_psblas
      USE class_cell
      USE class_vector
      USE class_iterating
      USE class_mesh
      USE class_output
      USE class_vector_field
      !
      USE tools_output_basics

      IMPLICIT NONE
      !
      TYPE(vector_field), INTENT(IN) :: fld
      CHARACTER(len=*),   INTENT(IN) :: field
      TYPE(output),       INTENT(INOUT) :: out
      TYPE(iterating),    INTENT(IN), OPTIONAL :: iter
      !
      INTEGER :: err_act, icontxt, info, mypnum, ncells_glob
      TYPE(vector), ALLOCATABLE :: x_loc(:)
      TYPE(vector), ALLOCATABLE :: x_glob(:)
      CHARACTER(len=32) :: path
      TYPE(cell), ALLOCATABLE :: cells_glob(:)
      TYPE(mesh), POINTER :: msh => NULL()

      ! Sets error handling for PSBLAS-2 routines
      CALL psb_erractionsave(err_act)

      mypnum  = mypnum_()
      icontxt = icontxt_()

      ! Sets output path
      IF(PRESENT(iter)) CALL out%set_output_path(iter)
      path = out%path_()

!!$  msh   => msh_(fld)
      CALL fld%get_mesh(msh)
      CALL fld%get_x(x_loc)

      ! Is FLD cell-centered?
      IF(fld%on_faces_()) THEN
          WRITE(*,100)
          CALL abort_psblas
      END IF

      ! Global number of cells
      ncells_glob = psb_cd_get_global_cols(msh%desc_c)

      ALLOCATE(x_glob(ncells_glob),stat=info)
      IF(info /= 0) THEN
          WRITE(*,200)
          CALL abort_psblas
      END IF

      ! Gathers cell-centered values
      CALL l2g_vector(x_loc,x_glob,msh%desc_c)
!!$  call psb_check_error(info,'set_vector_field','psb_gather',icontxt)

      ! Local-to-global reallocation of CELL object
      CALL l2g_cell(msh%cells,cells_glob,msh%desc_c)

      !if(mypnum == 0) then
      !   call wr_vtk_field(field,x_glob,msh%ncd,path)
      !end if

      IF(ALLOCATED(cells_glob)) CALL free_cell(cells_glob)
      DEALLOCATE(x_glob,x_loc)
      NULLIFY(msh)


      ! ----- Normal Termination -----
      CALL psb_erractionrestore(err_act)

100   FORMAT(' ERROR! Face-centered field in set_vector_field')
200   FORMAT(' ERROR! Memory allocation failure in set_vector_field')
300   FORMAT(' ERROR! Unsupported output format in set_vector_field')

  END SUBROUTINE set_vector_field

  MODULE PROCEDURE write_vtkmofo
      USE vtk_attributes, ONLY : attributes, scalar
      USE vtk_cells,      ONLY : vtkcell, vtkcell_list, set_cell_type
      USE vtk_datasets,   ONLY : unstruct_grid
      USE vtk,            ONLY : vtk_legacy_write
      !        USE vtk_vars
      !! author: Ian Porter, NRC
      !! date: 01/23/2019
      !!
      !! This subroutine translates morfeus data structure
      !! into a vtkmofo data structure and writes the vtk file
      !!
      !! For an unstructured grid, the legacy vtk file should look like this:
      !!
      !! DATASET UNSTRUCTURED_GRID
      !! POINTS n dataType
      !! p0x p0y p0z
      !! p1x p1y p1z
      !! ...
      !! p(n-1)x p(n-1)y p(n-1)z
      !! CELLS m size
      !! numPoints0, i, j, k, l, ...
      !! numPoints1, i, j, k, l, ...
      !! numPoints2, i, j, k, l, ...
      !! ...
      !! numPointsm-1, i, j, k, l, ...
      !! CELL_TYPES m
      !! type0
      !! type1
      !! type2
      !! ...
      !! typem-1
      !!
      INTEGER :: i, ibegin, itype
      CLASS(vtkcell), ALLOCATABLE :: dummy
      TYPE (unstruct_grid)        :: vtk_mesh
      TYPE (attributes), DIMENSION(n_cell_values)  :: cell_vals_to_write
      TYPE (attributes), DIMENSION(n_point_values) :: point_vals_to_write
      TYPE(vtkcell_list), DIMENSION(:), ALLOCATABLE   :: morfeus_cells

      ALLOCATE(morfeus_cells(1:SIZE(cell_ids))); ibegin = 1
      DO i = 1, SIZE(morfeus_cells)
          SELECT CASE (icverts(i))
              !! Take the # of vertices and convert it to the type of cell
          CASE (3)
              itype = 5         !! Triangle
          CASE (4)
              IF (ncd>2) THEN
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
      DO i = 1, n_cell_values
          !! Cell values
          IF (.NOT. ALLOCATED(cell_vals_to_write(i)%attribute)) THEN
              ALLOCATE(scalar::cell_vals_to_write(i)%attribute)
          END IF
          SELECT CASE (i)
          CASE (1)  !! Cell IDs
              CALL cell_vals_to_write(i)%attribute%init (TRIM(cell_dataname(i)), numcomp=1, int1d=cell_ids)
          CASE (2)  !! Processor IDs
              CALL cell_vals_to_write(i)%attribute%init (TRIM(cell_dataname(i)), numcomp=1, int1d=iproc)
          CASE (3)  !! Temperature
              CALL cell_vals_to_write(i)%attribute%init (TRIM(sname(1)), numcomp=1, real1d=sfield)
          END SELECT
      END DO

      !        DO i = 1, n_point_values
      !            !! Point values
      !            CALL point_vals_to_write(i)%attribute%init (TRIM(point_dataname(i)), numcomp=1, real1d=point_vals(:,i))
      !        END DO

      CALL vtk_legacy_write (unit=15, geometry=vtk_mesh, filename=file_name, celldatasets=cell_vals_to_write, &
          &                    multiple_io=.FALSE., title=title)
      !          &                    pointdatasets=point_vals_to_write, multiple_io=.TRUE., title=title)

    END PROCEDURE

END SUBMODULE vtkmofo_io_implementation
