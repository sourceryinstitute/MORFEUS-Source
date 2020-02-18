!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!

MODULE write_exodus

  USE iso_fortran_env, ONLY : INT64, REAL64
  USE class_psblas
  USE class_cell
  USE class_face
  USE class_vertex
  USE class_mesh
  USE class_connectivity
  USE class_iterating
  USE class_scalar_field
  USE class_vector_field
  USE class_vector
  USE class_output

  IMPLICIT NONE
  INCLUDE 'exodusII.inc'

  !! author: Hari Radhakrishnan, GSE
  !! date 12/12/2019
  !!
  !! This module implements the routines for writing an exodus file with mesh and results

CONTAINS

  SUBROUTINE write_exo_morfeus(msh, scalars, vectors, iter)

    TYPE(mesh), INTENT(IN) :: msh
    TYPE(iterating), INTENT(IN), OPTIONAL :: iter
    TYPE(scalar_field), INTENT(IN), DIMENSION(:), OPTIONAL :: scalars
    TYPE(vector_field), INTENT(IN), DIMENSION(:), OPTIONAL :: vectors
    !
    CHARACTER(len=100) :: filename
    INTEGER, DIMENSION(:), ALLOCATABLE :: cell_ids
    INTEGER, DIMENSION(:), ALLOCATABLE :: v2cconn
    INTEGER, DIMENSION(:), ALLOCATABLE :: icverts
    INTEGER, DIMENSION(:), ALLOCATABLE :: iproc
    REAL(psb_dpk_), DIMENSION(:,:), ALLOCATABLE :: points

    INTEGER(kind = int64) :: ierr, titlelen
    INTEGER(kind = int64) :: exodus_file_id, num_dims, num_nodes, num_elems,  &
      num_elem_blks, num_node_sets, num_side_sets
    INTEGER(kind = int64) :: cpu_ws, io_ws, num_props, prop_value
    REAL :: vers
    REAL(kind = real64), ALLOCATABLE :: xx(:), yy(:), zz(:)
    INTEGER :: info, err_act
    INTEGER :: icontxt, mypnum
    INTEGER :: i, ig, ic, nc2g, nv2c, nconn, ncells, ngc, ngroups
    INTEGER(kind = int64), ALLOCATABLE :: conn(:)
    CHARACTER(len = MXSTLN) :: coord_names(3), ctype
    INTEGER, POINTER :: ic2g(:) => NULL()
    INTEGER, POINTER :: iv2c(:) => NULL()
    INTEGER, ALLOCATABLE  :: igroup(:)
    INTEGER, ALLOCATABLE :: i_loc(:)
    TYPE(cell), ALLOCATABLE :: cells(:)
    TYPE(face), ALLOCATABLE :: faces(:)
    TYPE(vertex), ALLOCATABLE :: verts(:)
    TYPE(connectivity) :: v2f, v2c, f2c
    CHARACTER(len=32) :: path
    REAL(psb_dpk_), ALLOCATABLE :: scalar_local(:), scalar_global(:), f_s(:)
    TYPE(vector), ALLOCATABLE :: vector_local(:), vector_global(:)
    REAL(psb_dpk_), ALLOCATABLE :: f_x(:), f_y(:), f_z(:)
    ! Sets error handling for PSBLAS-2 routines
    CALL psb_erractionsave(err_act)
    mypnum  = mypnum_()
    icontxt = icontxt_()

    ! Sets output path
    ! IF(PRESENT(iter)) CALL out%set_output_path(iter)
    ! path = out%path_()
    filename = "out_nemo.e" !TRIM(path)
    print *, filename
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

    ! Open exodusII file for writing
    cpu_ws = 0
    io_ws = 8
    call exopts (EXVRBS, IERR) !Verbose error reporting of EXODUS routines
    exodus_file_id = excre("output.exo", EXCLOB, cpu_ws, io_ws, ierr)

    num_dims = msh%ncd
    num_elem_blks = ngroups
    num_node_sets = 0
    num_side_sets = msh%nbc
    num_nodes = size(verts)
    num_elems = size(cells)

    ! Write the database parameters to the exodus file
    CALL expini(exodus_file_id, "MORFEUS Output", num_dims, num_nodes, &
      num_elems, num_elem_blks, num_node_sets, num_side_sets, ierr)

    ! Write the mesh coordinates
    ALLOCATE(xx(num_nodes), yy(num_nodes), zz(num_nodes))

    xx(:) = verts(:)%x_()
    yy(:) = verts(:)%y_()
    zz(:) = verts(:)%z_()

    CALL expcor(exodus_file_id, xx, yy, zz, ierr)

    ! Write the mesh coordinate names
    coord_names(1) = "X"
    coord_names(2) = "Y"
    coord_names(3) = "Z"
    CALL expcon(exodus_file_id, coord_names, ierr)

    ! Write the element block parameters
    DO ig = 1, num_elem_blks !Loop over groups
      CALL msh%c2g%get_ith_conn(ic2g,ig)
      nc2g = SIZE(ic2g) !Number of cells in group
      CALL msh%v2c%get_ith_conn(iv2c,1) ! First cell
      nv2c = SIZE(iv2c) ! Number of nodes per cell
      ctype = "HEX8"
      CALL expelb(exodus_file_id, INT(ig,INT64), ctype, INT(nc2g,INT64), INT(nv2c,INT64), 1_INT64, ierr)

      nconn = nv2c*nc2g !Connectivity vector size

      print *, ig, nconn
      IF (ALLOCATED(conn)) DEALLOCATE(conn)
      ALLOCATE(conn(nconn))
      i = 0
      DO ic = 1, nc2g
        CALL msh%v2c%get_ith_conn(iv2c,ic)
        !print *, ic, i, iv2c
        conn(i+1:i+nv2c) = iv2c
        i = i + nv2c
      END DO
      call expelc(exodus_file_id, INT(ig,INT64), conn, ierr)
    END DO

    ! Loop over scalar fields and gather them on the root processor
    DO i = 1, size(scalars)
      ! Is the field cell-centered?
      IF(scalars(i)%on_faces_()) THEN
          WRITE(*,100)
          CYCLE
      END IF
      CALL scalars(i)%get_x(scalar_local)


      ALLOCATE(scalar_global(ncells),stat=info)
      IF(info /= 0) THEN
          WRITE(*,200)
          CALL abort_psblas
      END IF

      ! Gathers cell-centered values
      CALL psb_gather(scalar_global,scalar_local,msh%desc_c,info,root=0)

      IF(mypnum == 0) THEN
        DO ig = 1, num_elem_blks
          CALL msh%c2g%get_ith_conn(ic2g,ig)
          nc2g = SIZE(ic2g) !Number of cells in group
          IF(ALLOCATED(f_s)) DEALLOCATE(f_s)
          ALLOCATE(f_s(nc2g))
          f_s(1:nc2g) = scalar_global(ic2g(1:nc2g))
          CALL expev(exodus_file_id, 0.0, i, INT(ig, INT64), INT(nc2g, INT64), f_s, ierr)
        END DO
      END IF
    END DO

    ! Loop over vector fields and gather them on the root processor
    DO i = 1, size(vectors)
      ! Is the field cell-centered?
      IF(vectors(i)%on_faces_()) THEN
          WRITE(*,100)
          CYCLE
      END IF
      CALL vectors(i)%get_x(vector_local)


      ALLOCATE(vector_global(ncells),stat=info)
      IF(info /= 0) THEN
          WRITE(*,200)
          CALL abort_psblas
      END IF

      ! Gathers cell-centered values
      CALL l2g_vector(vector_global,vector_local,msh%desc_c)

      IF(mypnum == 0) THEN
        DO ig = 1, num_elem_blks
          CALL msh%c2g%get_ith_conn(ic2g,ig)
          nc2g = SIZE(ic2g) !Number of cells in group
          IF(ALLOCATED(f_x)) DEALLOCATE(f_x, f_y, f_z)
          ALLOCATE(f_x(nc2g), f_y(nc2g), f_z(nc2g))
          f_x(1:nc2g) = vector_global(ic2g(1:nc2g))%x_()
          f_y(1:nc2g) = vector_global(ic2g(1:nc2g))%y_()
          f_z(1:nc2g) = vector_global(ic2g(1:nc2g))%z_()
          CALL expev(exodus_file_id, 0.0, i, INT(ig, INT64), INT(nc2g, INT64), f_s, ierr)
        END DO
      END IF
    END DO


    ! DO i = 1, size(scalars) + 3*size(vectors)
    !   IF (i <= size(scalars)) THEN
    !     !write scalars to exodus file
    !   ELSE

    !   END IF
    ! END DO
    ! Close the exodus file
    CALL exclos(exodus_file_id, ierr)
    STOP

100 FORMAT(' WARNING! Face-centered field in WRITE_EXODUS. Field will not be output.')
200 FORMAT(' ERROR! Memory allocation failure in WRITE_EXODUS')

  END SUBROUTINE write_exo_morfeus

END MODULE
