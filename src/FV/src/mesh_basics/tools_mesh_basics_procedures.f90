SUBMODULE(tools_mesh_basics) tools_mesh_basics_procedures
  IMPLICIT NONE
CONTAINS
  include "geom_face.f90"
  include "geom_cell.f90"
  include "geom_diff.f90"
  include "geom_tet_center.f90"
  include "geom_tet_dihedral_angle.f90"
  include "geom_hex_dihedral_angle.f90"
  include "geom_tet_quality.f90"
  include "geom_hex_quality.f90"
  include "geom_tet_volume.f90"
END SUBMODULE tools_mesh_basics_procedures
