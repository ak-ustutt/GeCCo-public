*----------------------------------------------------------------------*
      subroutine arg_inp(l_exit,l_infile,l_logfile,
     &                        name_infile,name_logfile)
*----------------------------------------------------------------------*
*     process the command line argument list to GeCCo
*----------------------------------------------------------------------*

      implicit none
      
      logical, intent(out) ::
     &     l_exit, l_infile, l_logfile
      character(len=*), intent(out) ::
     &     name_infile, name_logfile
      
      integer ::
     &     iarg, nargs
      logical ::
     &     l_arg_ml

      character ::
     &     argstr*256

      l_exit = .false.
      l_infile = .false.
      l_logfile = .false.
      l_arg_ml = .false.

      nargs = command_argument_count()

      do iarg = 1, nargs
        call getarg(iarg,argstr)

        select case(trim(argstr))
        case('-h','--help')
          call help()
          l_exit = .true.

        case('-l','--logfile')
          l_arg_ml = .true.

        case default
          if (l_arg_ml) then
            name_logfile=argstr
            l_logfile = .true.
            l_arg_ml = .false.
          else
            name_infile=argstr
            l_infile = .true.
            ! issue error, if this was not the last argument
            if (iarg.ne.nargs) then
              call help()
              call quit(0,'arg_inp','arguments after input file?')
            end if
          end if

        end select

      end do

      return

      end

