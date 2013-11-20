*----------------------------------------------------------------------*
      subroutine write_title(lulog,style,title)
*----------------------------------------------------------------------*
      implicit none

      include 'write_styles.h'

      integer, intent(in) ::
     &     lulog, style
      character, intent(in) ::
     &     title*(*)
      
      integer ::
     &     lentitle, nxspace
      character(256) ::
     &     fmttit, fmtline

      lentitle = len_trim(title)
      select case(mod(style,10))
      case(wst_left,wst_left_wide)
        nxspace = 0
      case(wst_center,wst_center_wide)
        nxspace = max(0,(number_of_columns-lentitle)/2)
      case(wst_right,wst_right_wide)
        nxspace = max(0,number_of_columns-lentitle-1)
      end select

      select case(style-mod(style,10))
      case(wst_noframe)
        write(fmttit,'("(x,a)")')
      case(wst_uline_single,wst_lines_single)
        write(fmtline,'("(x,",i3,"(""-""))")') lentitle+2
        write(fmttit,'("(2x,a)")')
      case(wst_uline_double,wst_lines_double)
        write(fmtline,'("(x,",i3,"(""=""))")') lentitle+2
        write(fmttit,'("(2x,a)")')
      case(wst_around_single)
        write(fmtline,'("(x,""+"",",i3,"(""-""),""+"")")') lentitle+2
        write(fmttit,'("(x,""|"",x,a,x,""|"")")')
      case(wst_around_double)
        write(fmtline,'("(x,""+"",",i3,"(""=""),""+"")")') lentitle+2
        write(fmttit,'("(x,""|"",x,a,x,""|"")")')
      end select

      select case(style-mod(style,10))
      case(wst_noframe)
        write(lulog,fmttit) trim(title)
      case(wst_uline_single,wst_uline_double)
        write(lulog,fmttit) trim(title)
        write(lulog,fmtline)
      case(wst_lines_single,wst_lines_double,
     &     wst_around_single,wst_around_double)
        write(lulog,fmtline)
        write(lulog,fmttit) trim(title)
        write(lulog,fmtline)
      end select

      return
      end
