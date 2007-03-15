*----------------------------------------------------------------------*
      subroutine arg_inp(l_exit,l_infile,name_infile)
*----------------------------------------------------------------------*
*     process the command line argument list to GeCCo
*----------------------------------------------------------------------*

      implicit none
      
      logical, intent(out) ::
     &     l_exit, l_infile
      character, intent(out) ::
     &     name_infile*(*)
      
      integer ::
     &     iarg, nargs

      character ::
     &     argstr*256

      l_exit = .false.
      l_infile = .false.

      nargs = iargc()

      do iarg = 1, nargs
        call getarg(iarg,argstr)

        select case(trim(argstr))
        case('-h','--help')
          call help()
          l_exit = .true.

        case default
          name_infile=argstr
          l_infile = .true.
          ! issue error, if this was not the last argument
          if (iarg.ne.nargs) then
            call help()
            call quit(0,'arg_inp','arguments after input file?')
          end if

        end select

      end do

      return

      end

