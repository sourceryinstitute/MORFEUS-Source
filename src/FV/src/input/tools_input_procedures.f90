SUBMODULE(tools_input) tools_input_procedures
  IMPLICIT NONE
CONTAINS
  include "get_par_l.f90"
  include "get_par_i.f90"
  include "get_par_d.f90"
  include "get_par_h.f90"
  include "get_par_v.f90"
  include "read_par_l.f90"
  include "read_par_i.f90"
  include "read_par_d.f90"
  include "read_par_h.f90"
  include "read_par_v.f90"
  include "find_section.f90"
END SUBMODULE tools_input_procedures
