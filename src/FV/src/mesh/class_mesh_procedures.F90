!
!     (c) 2019 Guide Star Engineering, LLC
!     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
!     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
!     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
!
SUBMODULE(class_mesh) class_mesh_procedures
    USE psb_base_mod
    USE class_psblas
    USE class_cell
    USE class_connectivity
    USE class_face
    USE class_vector
    USE class_keytable
    USE class_surface
    USE class_vertex
    IMPLICIT NONE

CONTAINS

    ! ----- Constructors -----

    MODULE PROCEDURE create_mesh
        !! Global constructor
        USE tools_mesh,          ONLY : rd_inp_mesh, rd_gambit_mesh, rd_cgns_mesh, cmp_mesh_c2c, cmp_mesh_desc,  &
            &                           cmp_mesh_f2b, cmp_mesh_f2f, cmp_mesh_part, cmp_mesh_renum, cmp_mesh_v2b, &
            &                           cmp_mesh_v2ve, cmp_moving_surf, supplement_v2c, supplement_v2f
        USE tools_mesh_basics,   ONLY : geom_face, geom_cell, geom_diff
        USE tools_output_basics, ONLY : wr_mtx_pattern
        USE class_least_squares, ONLY : set_least_squares

        LOGICAL :: mtx_pat
        INTEGER :: icontxt, mypnum, nprocs
        INTEGER :: irenum, ipart, nswpref
        INTEGER :: nverts, nfaces, ncells
        REAL(psb_dpk_) :: scale
        CHARACTER(len=nlen) :: mesh_file

        icontxt = icontxt_()
        mypnum = mypnum_()
        nprocs = nprocs_()

        ! Preliminary check
        IF(msh%set) THEN
            WRITE(*,100)
            CALL abort_psblas
        END IF

        ! Reads mesh-related input parameters
        CALL rd_inp_mesh(input_file,sec,&
            & mesh_file,scale,irenum,ipart,nswpref,mtx_pat)

        ! Imports the mesh on P0 (verts, faces, cells, v2c, v2f, f2c,c2g ...)
        CALL import_mesh(msh,mesh_file)

#ifdef DEBUGMSH
        CALL nemo_mesh_size(msh)
#endif

        ! Applies the scaling factor
        IF (mypnum==0) THEN
            msh%verts = scale * msh%verts

            nverts = SIZE(msh%verts)
            nfaces = SIZE(msh%faces)
            ncells = SIZE(msh%cells)
        ENDIF

!!$ call cmp_mesh_v2v(nverts,msh%v2c,msh%v2v)         ! Computes MESH%V2V
        CALL cmp_mesh_v2ve(msh%ncd,msh%v2f,msh%v2v)          !
        CALL cmp_mesh_f2f(nfaces,msh%f2c,msh%f2f)            ! Computes MESH%F2F
        CALL cmp_mesh_c2c(msh%faces,msh%f2c,msh%c2c)         ! Computes MESH%C2C
        CALL cmp_mesh_v2b(msh%v2f,msh%faces,msh%nbc,msh%v2b) ! Computes MESH%V2B
        IF(mypnum == 0) WRITE(*,'()')

        ! REMARK: CMP_MESH_V2B must be called BEFORE the G2L_MESH!

#ifdef DEBUGMSH
        CALL nemo_mesh_size(msh)
#endif

        ! Dump mesh-related sparsity pattern before renumbering
        IF(mtx_pat) THEN
            CALL wr_mtx_pattern(msh%c2c,'pattern_cell_1.mtx')
!!$       call wr_mtx_pattern(msh%v2v,'pattern_vert_1.mtx')
        END IF

        ! Computes and applies matrix renumbering
        CALL cmp_mesh_renum(irenum,msh%cells,msh%faces,msh%c2c,&
            & msh%f2c,msh%v2c,msh%c2g)

        ! Dumps mesh-related sparsity pattern after renumbering
        IF(mtx_pat .AND. irenum > 0) THEN
            CALL wr_mtx_pattern(msh%c2c,'pattern_cell_2.mtx')
!!$       call wr_mtx_pattern(msh%v2v,'pattern_vert_2.mtx')
        END IF

        ! Detects conceptual surfaces, needed for surface-constrained
        ! vertex motion. Conceptual surfaces cannot exist in 2d.
        IF (msh%ncd == 3) THEN
            CALL cmp_moving_surf(msh%nbc,msh%v2b,msh%verts,msh%surf)
!!!   to be done bcast_surf
        ENDIF

        ! Calculates the partitioning according to IPART value
        CALL cmp_mesh_part(ipart,nswpref,msh%c2c)

        ! Allocates and assemblies the MSH%DESC descriptor
        CALL cmp_mesh_desc(msh%v2v,msh%v2c,msh%f2f,msh%f2c,msh%c2c,&
            &msh%desc_v,msh%desc_f,msh%desc_c)


        ! Store supplemental connectivity info for smoothing the mesh
        CALL supplement_v2c(msh%v2c,msh%desc_v,msh%ov2c_sup,msh%c2ov_sup)
        CALL bcast_face(msh%faces)
        CALL supplement_v2f(msh%v2f,msh%faces, msh%desc_v,msh%ov2f_sup,msh%f2ov_sup)

#ifdef DEBUGMSH
        CALL  nemo_mesh_size(msh)

        IF (mypnum ==0 ) CALL pr_mesh_size(msh)
#endif

        ! ----- GLOBAL 2 LOCAL REALLCATION ----
        CALL g2l_mesh(msh)

#ifdef DEBUGMSH
        CALL  nemo_mesh_size(msh)

        IF (mypnum ==0 ) CALL pr_mesh_size(msh)
#endif

        ! ----- From here on data are LOCAL -----

        CALL cmp_mesh_f2b(msh%faces,msh%nbc,msh%f2b) ! Computes MESH%F2B
        IF(mypnum == 0) WRITE(*,'()')

        IF(mypnum == 0) THEN
            WRITE(*,*) 'Computing mesh metrics'
            WRITE(*,*)
        END IF

        CALL sw_geo%tic()

        ! Computes face-related metrics members MSH%FACE_CNTR, MSH%AF, MSH%AREA
        CALL geom_face(msh%verts,msh%v2f,msh%ncd,&
            & msh%face_cntr,msh%af,msh%area)

        ! Computes cell-related metrics members MSH%CELL_CNTR, MSH%VOL
        CALL geom_cell(msh%verts,msh%faces,msh%cells,msh%v2f,msh%v2c,msh%f2c,&
            & msh%ncd,msh%cell_cntr,msh%vol)

        ! Computes face-related metrics members MSH%DF, MSH%DIST, MSH%INTERP
        CALL geom_diff(msh%faces,msh%f2b,msh%face_cntr,msh%af,msh%cell_cntr, &
            & msh%df,msh%dist,msh%interp)

        ! Computes metrics for cell-centered least squares regression
        CALL set_least_squares(msh%lsr,msh%ncd,msh%desc_c,msh%c2c,msh%f2b, &
            & msh%faces,msh%cell_cntr,msh%face_cntr)

        CALL sw_geo%toc()

        ! Sets MSH%SET logical flag on .true.
        msh%set = .TRUE.
        msh%ngp = msh%c2g%nel_()

        IF(mypnum == 0) THEN
            WRITE(*,200) 'Mesh summary'
            WRITE(*,300) '- meshfile:', TRIM(ADJUSTL(mesh_file))
            WRITE(*,400) '- Number of cells:   ', ncells
            WRITE(*,400) '- Number of faces:   ', nfaces
            WRITE(*,400) '- Number of vertices:', nverts
            WRITE(*,500) '- Number of bc:      ', msh%nbc
            WRITE(*,500) '- Number of groups:  ', msh%ngp
            WRITE(*,*)
        END IF

100     FORMAT(' ERROR! MESH object in CREATE_MESH has been already set')
200     FORMAT(1x,a)
300     FORMAT(1x,a,1x,a)
400     FORMAT(1x,a,i8)
500     FORMAT(1x,a,i5)

    END PROCEDURE create_mesh


    SUBROUTINE import_mesh(msh,mesh_file)
        USE tools_mesh
        IMPLICIT NONE
        !
        TYPE(mesh), INTENT(INOUT) :: msh
        CHARACTER(len=*), INTENT(INOUT) :: mesh_file
        !
        INTEGER :: i, intbuf(2)
        INTEGER :: icontxt, mypnum
        CHARACTER(len=4) :: fmt


        icontxt = icontxt_()
        mypnum  = mypnum_()

        CALL sw_msh%tic()

        ! Reads the mesh object on P0...
        IF(mypnum == 0) THEN

            ! Sets mesh format
            mesh_file = ADJUSTR(mesh_file)
            i = SCAN(mesh_file,'.',back=.TRUE.)
            IF(i == 0) THEN
                WRITE(*,100)
                WRITE(*,'(a24,a80)')'Expected to find file: ',mesh_file
                CALL abort_psblas
            END IF
            i = i + 1
            fmt = mesh_file(i:LEN(mesh_file))

            SELECT CASE (fmt)
            CASE('cgns')
#ifdef HAVE_CGNS
                CALL rd_cgns_mesh(mesh_file, msh%id, msh%nbc, msh%ncd, &
                    & msh%verts, msh%faces, msh%cells, &
                    & msh%v2f, msh%v2c, msh%f2c, msh%c2g)
#else
                error STOP "MORFEUS built without CNGS support! Unable to continue."
#endif
            CASE('e')
#ifdef HAVE_EXODUS
                CALL rd_exodus_mesh(mesh_file, msh%id, msh%nbc, msh%ncd, &
                    & msh%verts, msh%faces, msh%cells, &
                    & msh%v2f, msh%v2c, msh%f2c, msh%c2g)
#else
                error STOP "MORFEUS built without exodus support! Unable to continue."
#endif
            CASE('neu')
                CALL rd_gambit_mesh(mesh_file, msh%id, msh%nbc, msh%ncd, &
                    & msh%verts, msh%faces, msh%cells, &
                    & msh%v2f, msh%v2c, msh%f2c, msh%c2g)
            CASE default
                WRITE(*,200)
                CALL abort_psblas
            END SELECT

        END IF

        ! ... then bcast msh%(id,nbc,ncd)
        IF(mypnum == 0) THEN

            CALL psb_bcast(icontxt,msh%id)

            intbuf(1) = msh%nbc
            intbuf(2) = msh%ncd
            CALL psb_bcast(icontxt,intbuf)
        ELSE
            CALL psb_bcast(icontxt,msh%id)

            CALL psb_bcast(icontxt,intbuf)
            msh%nbc = intbuf(1)
            msh%ncd = intbuf(2)
        END IF

        CALL sw_msh%toc()

100     FORMAT('ERROR! Unsupported mesh filename in IMPORT_MESH')
200     FORMAT('ERROR! Unsupported mesh format in IMPORT_MESH')

    END SUBROUTINE import_mesh

    ! ----- Global To Local -----

    SUBROUTINE g2l_mesh(msh)
        IMPLICIT NONE
        !!
        TYPE(mesh), INTENT(INOUT) :: msh

        CALL sw_g2l%tic()

        IF(mypnum_() == 0) THEN
            WRITE(*,*) 'Global to local reallocation of MESH object'
            WRITE(*,*)
        END IF

        CALL g2l_conn(msh%v2f,msh%desc_v,msh%desc_f)    ! Reallocates MSH%V2F
        CALL g2l_face(msh%faces,msh%desc_f,msh%desc_c)  ! Reallocates MSH%FACES
        CALL bcast_cell(msh%cells)                      ! Broadcast   MSH%CELLS
        CALL g2l_cell(msh%cells,msh%desc_c)             ! Reallocates MSH%CELLS
        CALL bcast_vertex(msh%verts)                    ! Broadcast   MSH%VERTS
        CALL g2l_vertex(msh%verts,msh%desc_v)           ! Reallocates MSH%VERTS

        CALL g2l_conn(msh%v2c,msh%desc_v,msh%desc_c)    ! Reallocates MSH%V2C
        CALL bcast_conn(msh%c2g)                        ! Broadcast   MSH%C2G
        CALL g2l_conn(msh%c2g,msh%desc_c)               ! Reallocates MSH%C2G
        !

        CALL bcast_conn(msh%v2b)                        ! Broadcast   MSH%V2B
        CALL g2l_conn(msh%v2b,msh%desc_v)               ! Reallocates MSH%V2B

        CALL sw_g2l%toc()

    END SUBROUTINE g2l_mesh

    MODULE PROCEDURE free_mesh
    !! ----- Destructor -----
        USE psb_base_mod
        USE class_least_squares, ONLY : free_least_squares
        IMPLICIT NONE
        !
        INTEGER :: err_act, icontxt, info
        INTEGER :: ib

        ! Sets error handling for PSBLAS-2 routines
        info = 0
        CALL psb_erractionsave(err_act)

        icontxt = icontxt_()

        IF(mypnum_() == 0) WRITE(*,*) 'Deallocating mesh object'

        CALL free_conn(msh%v2c)
        CALL free_conn(msh%v2f)
        CALL free_conn(msh%f2c)
        CALL free_conn(msh%c2g)
        !
        CALL free_conn(msh%v2v)
        CALL free_conn(msh%f2f)
        CALL free_conn(msh%c2c)
        !
        CALL free_conn(msh%v2b)
        CALL free_conn(msh%f2b)

        IF ( ALLOCATED(msh%surf) ) THEN
            DO ib = 1, msh%nbc
                CALL msh%surf(ib)%free_surface()
            END DO
            DEALLOCATE(msh%surf)
        ENDIF

        ! these supplemental info keytables are only filled for
        ! parallel runs, so check before freeing

        IF ( msh%ov2c_sup%exists() ) CALL msh%ov2c_sup%free_keytable()
        IF ( msh%c2ov_sup%exists() ) CALL msh%c2ov_sup%free_keytable()

        IF ( msh%ov2f_sup%exists() ) CALL msh%ov2f_sup%free_keytable()
        IF ( msh%f2ov_sup%exists() ) CALL msh%f2ov_sup%free_keytable()

        CALL free_vertex(msh%verts)
        CALL free_face(msh%faces)
        CALL free_cell(msh%cells)

        CALL psb_cdfree(msh%desc_v,info)
        CALL psb_check_error(info,'free_mesh%desc_v','psb_cdfree',icontxt)
        !
        CALL psb_cdfree(msh%desc_f,info)
        CALL psb_check_error(info,'free_mesh%desc_f','psb_cdfree',icontxt)
        !
        CALL psb_cdfree(msh%desc_c,info)
        CALL psb_check_error(info,'free_mesh%desc_c','psb_cdfree',icontxt)

        ! Face-related metrics
        DEALLOCATE(msh%area)
        DEALLOCATE(msh%dist)
        DEALLOCATE(msh%interp)
        CALL free_vector(msh%face_cntr)
        CALL free_vector(msh%af)
        CALL free_vector(msh%df)

        ! Cell-related metrics
        DEALLOCATE(msh%vol)
        CALL free_vector(msh%cell_cntr)

        ! Least Squares metrics
        CALL free_least_squares(msh%lsr)

        msh%set = .FALSE.

        IF(mypnum_() == 0) WRITE(*,*)

        ! ----- Normal Termination -----
        CALL psb_erractionrestore(err_act)

    END PROCEDURE free_mesh

    MODULE PROCEDURE check_mesh_consistency
        IMPLICIT NONE
        !! Checks the consistency of two meshes: MSH1 and MSH2

        IF(.NOT. ASSOCIATED(msh1, msh2)) THEN
            WRITE(*,100) TRIM(WHERE)
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! In ', a,': mesh pointers are not consistent')

    END PROCEDURE check_mesh_consistency

    ! ----- Check Operations -----

    MODULE PROCEDURE check_mesh_unused_el
        IMPLICIT NONE
        ! Scans through numerous connectivities, ensuring that in each one,
        ! all faces & cells are referenced at least once.

        INTEGER :: mypnum
        INTEGER :: failures

        mypnum = mypnum_()

        IF(mypnum == 0) WRITE(*,*) ' Checking C2C connectivity...'
        failures = msh%c2c%unused_elements()
        IF (failures > 0) THEN
            WRITE(*,100) failures, 'cells', 'MSH%C2C', mypnum
            CALL abort_psblas
        END IF

        IF(mypnum == 0) WRITE(*,*) ' Checking F2C connectivity...'
        failures = msh%f2c%unused_elements()
        IF (failures > 0) THEN
            WRITE(*,100) failures, 'faces', 'MSH%F2C', mypnum
            CALL abort_psblas
        END IF

        IF(mypnum == 0) WRITE(*,*) ' Checking V2C connectivity...'
        failures = msh%v2c%unused_elements()
        IF (failures > 0) THEN
            WRITE(*,100) failures, 'vertices', 'MSH%V2C', mypnum
            CALL abort_psblas
        END IF

        IF(mypnum == 0) WRITE(*,*) ' Checking V2F connectivity...'
        failures = msh%v2f%unused_elements()
        IF (failures > 0) THEN
            WRITE(*,100) failures, 'vertices', 'MSH%V2F', mypnum
            CALL abort_psblas
        END IF

100     FORMAT(' ERROR! ',i3,' unused ',a,' in ',a,' connectivity detected',&
            & ' on process ',i2)

    END PROCEDURE check_mesh_unused_el

    SUBROUTINE nemo_mesh_size(msh)
        IMPLICIT NONE
        TYPE(mesh), INTENT(IN) :: msh
        INTEGER(kind=nemo_int_long_)   ::isz
        REAL:: sz,sz1
        INTEGER :: icontxt, mypnum, nprocs,ip
        icontxt=icontxt_()
        mypnum = mypnum_()
        nprocs = nprocs_()
        isz = msh%nemo_sizeof()
        sz=isz
        IF(mypnum == 0) THEN
            WRITE(*,'()')
            WRITE(*,*) 'On process',mypnum, 'mesh size',sz
        END IF
        DO ip = 1, nprocs - 1
            IF(mypnum == ip) THEN
                CALL psb_snd(icontxt,sz,0)
            ELSEIF(mypnum == 0) THEN
                CALL psb_rcv(icontxt,sz1,ip)
                WRITE(*,*) 'On process',ip, 'mesh size', sz1
            END IF
        END DO
        CALL psb_sum(icontxt,isz)
        IF (mypnum==0) THEN
            WRITE(*,*) 'Total size of mesh   :',isz
            WRITE(*,*)
        ENDIF
    END SUBROUTINE nemo_mesh_size

    MODULE PROCEDURE nemo_mesh_sizeof
        USE psb_base_mod
        IMPLICIT NONE

        INTEGER(kind=nemo_int_long_)   :: val

        val = 3 * nemo_sizeof_int + LEN(msh%id)
        val = val + msh%v2c%nemo_sizeof() &
            & + msh%v2f%nemo_sizeof()&
            & + msh%f2c%nemo_sizeof()&
            & + msh%c2g%nemo_sizeof()&
            & + msh%v2v%nemo_sizeof()&
            & + msh%f2f%nemo_sizeof()&
            & + msh%c2c%nemo_sizeof()&
            & + msh%v2b%nemo_sizeof()&
            & + msh%f2b%nemo_sizeof()

        val = val + msh%ov2c_sup%nemo_sizeof() &
        & + msh%c2ov_sup%nemo_sizeof() &
        & + msh%ov2f_sup%nemo_sizeof() &
        & + msh%f2ov_sup%nemo_sizeof()

        IF (ALLOCATED(msh%verts))&
            & val = val + SUM(msh%verts%nemo_sizeof())
        IF (ALLOCATED(msh%faces))&
            & val = val + SUM(msh%faces%nemo_sizeof())
        IF (ALLOCATED(msh%cells))&
            & val = val + SUM(msh%cells%nemo_sizeof())

        val = val + psb_sizeof(msh%desc_v)&
            & + psb_sizeof(msh%desc_f)&
            & + psb_sizeof(msh%desc_c)

        IF (ALLOCATED(msh%area))&
            & val = val + nemo_sizeof_dp * SIZE(msh%area)
        IF (ALLOCATED(msh%dist))&
            & val = val + nemo_sizeof_dp * SIZE(msh%dist)
        IF (ALLOCATED(msh%interp))&
            & val = val + nemo_sizeof_dp * SIZE(msh%interp)

        IF (ALLOCATED(msh%face_cntr))&
            & val = val + SUM(msh%face_cntr%nemo_sizeof())
        IF (ALLOCATED(msh%af))&
            & val = val + SUM(msh%af%nemo_sizeof())
        IF (ALLOCATED(msh%df))&
            & val = val + SUM(msh%df%nemo_sizeof())

        IF (ALLOCATED(msh%vol))&
            & val = val + nemo_sizeof_dp * SIZE(msh%vol)
        IF (ALLOCATED(msh%cell_cntr))&
            & val = val + SUM(msh%cell_cntr%nemo_sizeof())

        IF (ALLOCATED(msh%lsr))&
            & val = val + SUM(msh%lsr%nemo_sizeof())
        IF (ALLOCATED(msh%surf))&
            & val = val + SUM(msh%surf%nemo_sizeof())

        nemo_mesh_sizeof = val

    END PROCEDURE nemo_mesh_sizeof

    SUBROUTINE pr_mesh_size(msh)
        USE psb_base_mod
        IMPLICIT NONE

        TYPE(mesh), INTENT(IN) :: msh

        WRITE(*,*) 'Size of v2c: ', msh%v2c%nemo_sizeof()
        WRITE(*,*) 'Size of v2f:' , msh%v2f%nemo_sizeof()
        WRITE(*,*) 'Size of f2c:' , msh%f2c%nemo_sizeof()
        WRITE(*,*) 'Size of c2g:' , msh%c2g%nemo_sizeof()
        WRITE(*,*) 'Size of v2v:' , msh%v2v%nemo_sizeof()
        WRITE(*,*) 'Size of f2f:' , msh%f2f%nemo_sizeof()
        WRITE(*,*) 'Size of c2c:' , msh%c2c%nemo_sizeof()
        WRITE(*,*) 'Size of v2b:' , msh%v2b%nemo_sizeof()
        WRITE(*,*) 'Size of f2b:' , msh%f2b%nemo_sizeof()

        WRITE(*,*) 'Size of v2c_sup:' , msh%ov2c_sup%nemo_sizeof()
        WRITE(*,*) 'Size of c2ov_sup:', msh%c2ov_sup%nemo_sizeof()
        WRITE(*,*) 'Size of v2f_sup:' , msh%ov2f_sup%nemo_sizeof()
        WRITE(*,*) 'Size of f2ov_sup:' , msh%f2ov_sup%nemo_sizeof()

        IF (ALLOCATED(msh%verts))&
            & WRITE(*,*) 'Size of verts:' , SUM(msh%verts%nemo_sizeof())
        IF (ALLOCATED(msh%faces))&
            &  WRITE(*,*) 'Size of faces:' , SUM(msh%faces%nemo_sizeof())
        IF (ALLOCATED(msh%cells))&
            & WRITE(*,*) 'Size of cells:' , SUM(msh%cells%nemo_sizeof())

        WRITE(*,*) 'Size of desc_v:' , psb_sizeof(msh%desc_v)
        WRITE(*,*) 'Size of desc_f:' , psb_sizeof(msh%desc_f)
        WRITE(*,*) 'Size of desc_c:' , psb_sizeof(msh%desc_c)

        IF (ALLOCATED(msh%area))&
            & WRITE(*,*) 'Size of msh%area:' , nemo_sizeof_dp * SIZE(msh%area)
        IF (ALLOCATED(msh%dist))&
            & WRITE(*,*) 'Size of msh%dist:' , nemo_sizeof_dp * SIZE(msh%dist)
        IF (ALLOCATED(msh%interp))&
            & WRITE(*,*) 'Size of mesh&interp:' , nemo_sizeof_dp * SIZE(msh%interp)

        IF (ALLOCATED(msh%face_cntr))&
            & WRITE(*,*) 'Size of msh&face_cntr:' , SUM(msh%face_cntr%nemo_sizeof())
        IF (ALLOCATED(msh%af))&
            & WRITE(*,*) 'Size of msh%af:' , SUM(msh%af%nemo_sizeof())
        IF (ALLOCATED(msh%df))&
            & WRITE(*,*) 'Size of msh%df:' , SUM(msh%df%nemo_sizeof())

        IF (ALLOCATED(msh%vol))&
            & WRITE(*,*) 'Size of msh%vol:' , nemo_sizeof_dp * SIZE(msh%vol)
        IF (ALLOCATED(msh%cell_cntr))&
            & WRITE(*,*) 'Size of msh%cell_cntr:' , SUM(msh%cell_cntr%nemo_sizeof())

        IF (ALLOCATED(msh%lsr))&
            & WRITE(*,*) 'Size of msh%lsr:' , SUM(msh%lsr%nemo_sizeof())
        IF (ALLOCATED(msh%surf))&
            & WRITE(*,*) 'Size of msh%surf:' , SUM(msh%surf%nemo_sizeof())

    END SUBROUTINE pr_mesh_size

END SUBMODULE class_mesh_procedures
