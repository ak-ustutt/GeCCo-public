*----------------------------------------------------------------------*
      subroutine write_title(luout,style,title)
*----------------------------------------------------------------------*
      implicit none

      include 'write_styles.h'

      integer, intent(in) ::
     &     luout, style
      character, intent(in) ::
     &     title*(*)
      
      integer ::
     &     lentitle, nxspace
      character ::
     &     fmttit*128, fmtline*128

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
        write(luout,fmttit) trim(title)
      case(wst_uline_single,wst_uline_double)
        write(luout,fmttit) trim(title)
        write(luout,fmtline)
      case(wst_lines_single,wst_lines_double,
     &     wst_around_single,wst_around_double)
        write(luout,fmtline)
        write(luout,fmttit) trim(title)
        write(luout,fmtline)
      end select

      return
      end
