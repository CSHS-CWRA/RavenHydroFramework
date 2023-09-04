#ifndef MFUSGPP_H
#define MFUSGPP_H
#ifdef _MODFLOW_USG_
using namespace std;

/* MODFLOW-USG Basic Model Interface (BMI) Fortran-imported functions
*
* All non-PBJ functions/routines are defined in Extern.f90
* PBJ routines are defined in rvn2pbju1.f90
*
* If you're ever getting a unresolved externals error, one thing to check
* is if you've declared a function with CAPITALS. The exported Fortran
* functions/subroutines are always lowercase.
*/
namespace MFUSG
{
	// Function declarations
	extern "C"
	{
		// Routines to run MFUSG
    void   init_mfusg            (char *fname, int *fnlen, int *export_nodes, int *export_nlay,
			                            int *export_nper, int *export_mxiter, int *ni_unit);
		void   init_stress_period    (int *export_kper);
		void   init_time_step        (int *export_kstp, int *export_IDOFLOW);
		void   formulate_boundaries  (int *export_kiter);
		void   solve_matrix          (int *export_ibflag, bool *export_icnvg);
		void   post_solution_updates ();
		void   calc_boundary_budgets ();
		void   print_save_data       ();
		void   mf_shutdown           ();

		// Routines for getting data from MFUSG
		void   get_mf_memory         (int *export_nodlay, int *export_nstp, double *export_perlen);
		double get_node_head         (int *node);
		double get_node_area         (int *node);
		int    get_node_ibound       (int *node);

		// Routines to work with MFUSG data
		void   mask_iunit            (int *export_iunit); // See Extern.f90 in MFUSGLib for iunit/package table (CUNIT)
		void   add_to_flow_eq        (int *node, double *AMAT_change, double *RHS_change);
	    void   add_to_flow_budget    (int *nnodes, int *rvn_nodes, double *rates, char *pckg_name);
        void   first_active_below    (int *node);
		void   set_irvncb            (int *unit_no);
		void   set_nrchop            (int *new_NRCHOP);

		// PBJ Routines
		void   init_pbj_settings     (int *rvn_nseg, int *rvn_mode, int *rvn_condtype);
		void   add_pbj_segment       (int *iseg, int *rvn_n1,     int *rvn_n2,     int *rvn_n3,
									                        double *rvn_wa1, double *rvn_wa2, double *rvn_wa3,
			                                    double *rvn_wb1, double *rvn_wb2, double *rvn_wb3,
									                        double *rvn_elev1, double *rvn_elev2,
			                                    double *rvn_cond1, double *rvn_cond2,
			                                    double *rvn_len);
		void   init_pbj_values       ();
		void   update_pbj_by_basin   (int *rvn_nseg, int *rvn_seg, double *rvn_basin_river_depth);
		double get_pbj_segment_flows (double *q);
	}

}
#endif
#endif
