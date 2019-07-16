/*
  !
  !     (c) 2019 Guide Star Engineering, LLC
  !     This Software was developed for the US Nuclear Regulatory Commission (US NRC)
  !     under contract "Multi-Dimensional Physics Implementation into Fuel Analysis under
  !     Steady-state and Transients (FAST)", contract # NRC-HQ-60-17-C-0007
  !
*/

int call_smooth2d(int *num_incident_vtx, int *num_incident_tri,
                  double *free_pos,double incident_vtx[][2],
                   int vtx_connectivity[][2],int *tangled);
