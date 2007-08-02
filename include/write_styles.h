      integer, parameter ::
     &     number_of_columns = 80
      integer, parameter ::
     &     wst_left   = 0,
     &     wst_center = 1,
     &     wst_right  = 2,
     &     wst_left_wide = 3,
     &     wst_center_wide = 4,
     &     wst_right_wide  = 5,
     &     wst_noframe = 0,
     &     wst_uline_single = 10,
     &     wst_uline_double = 20,
     &     wst_lines_single = 30,
     &     wst_lines_double = 40,
     &     wst_around_single = 50,
     &     wst_around_double = 60
      integer, parameter ::
     &     wst_dbg_subr = wst_lines_single+wst_left,
     &     wst_dbg_func = wst_lines_double+wst_left,
     &     wst_section  = wst_center_wide+wst_around_double,
     &     wst_subsection = wst_center+wst_uline_single,
     &     wst_title = wst_left+wst_uline_single
