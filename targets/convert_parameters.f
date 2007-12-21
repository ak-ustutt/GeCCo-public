*----------------------------------------------------------------------*
*     a set of routines to convert parameters to "actions" into strings
*     and vice versa
*
*     the following routines are contained
*
*       hop_parameters
*       xop_parameters
*       dens_parameters
*       cloneop_parameters
*       form_parameters
*       me_list_parameters
*       import_parameters
*       solve_parameters
*
*----------------------------------------------------------------------*
      subroutine hop_parameters(rw,parameters,
     &     min_rank,max_rank,iformal)

      implicit none
      
      integer, intent(inout) ::
     &     rw,min_rank,max_rank,iformal
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(8(i5,x))')
     &        min_rank,max_rank,iformal
      else
        read(parameters,'(8(i5,x))')
     &       min_rank,max_rank,iformal
      end if

      return
      end

      subroutine xop_parameters(rw,parameters,
     &     dagger,min_rank,max_rank,ncadiff,iformal)

      implicit none
      
      logical, intent(inout) ::
     &     dagger
      integer, intent(inout) ::
     &     rw,min_rank,max_rank,ncadiff,iformal
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(l,x,9(i5,x))')
     &       dagger,min_rank,max_rank,ncadiff,iformal
      else
        read(parameters,'(l,x,9(i5,x))')
     &       dagger,min_rank,max_rank,ncadiff,iformal
      end if

      return
      end

      subroutine dens_parameters(rw,parameters,
     &     min_rank,max_rank,iformal)

      implicit none
      
      integer, intent(inout) ::
     &     rw,min_rank,max_rank,iformal
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(8(i5,x))')
     &        min_rank,max_rank,iformal
      else
        read(parameters,'(8(i5,x))')
     &        min_rank,max_rank,iformal
      end if

      return
      end

      subroutine cloneop_parameters(rw,parameters,name,dagger)

      implicit none
      
      integer, intent(in) ::
     &     rw
      logical, intent(inout) ::
     &     dagger
      character, intent(inout) ::
     &     parameters*(*), name*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(l,x,a)')
     &       dagger,name
      else
        read (parameters,'(l,x,a)')
     &       dagger,name
      end if

      return
      end

      subroutine form_parameters(rw,
     &     parameters,n_par_str,title)

      implicit none
      
      integer, intent(in) ::
     &     rw, n_par_str
      character*(*), intent(inout) ::
     &     parameters(n_par_str),
     &     title

      if (rw.lt.0) then
        write(parameters(1),'(a)') title
      else
        read(parameters(1),'(a)') title
      end if

      return
      end

      subroutine opt_parameters(rw,parameters,ncat,nint)

      implicit none
      
      integer, intent(in) ::
     &     rw
      integer, intent(inout) ::
     &     ncat,nint
      character(*), intent(inout) ::
     &     parameters

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(2i4,a)') ncat,nint
      else
        read(parameters,'(2i4,a)') ncat,nint
      end if

      return
      end

      subroutine me_list_parameters(rw,parameters,
     &     absym,casym,gamma,s2,ms)

      implicit none
      
      integer, intent(inout) ::
     &     rw,absym,casym,gamma,s2,ms
      character, intent(inout) ::
     &     parameters*(*)

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(8(i5,x))')
     &        absym,casym,gamma,s2,ms
      else
        read(parameters,'(8(i5,x))')
     &       absym,casym,gamma,s2,ms
      end if

      return
      end

      subroutine import_parameters(rw,parameters,
     &     env_type)

      implicit none
      
      integer, intent(in) ::
     &     rw
      character(*), intent(inout) ::
     &     env_type
      character(*), intent(inout) ::
     &     parameters

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(a)') env_type
      else
        read(parameters,'(a)') env_type
      end if

      return
      end

      subroutine solve_parameters(rw,parameters,
     &     nopt,nroots)

      implicit none
      
      integer, intent(in) ::
     &     rw
      integer, intent(inout) ::
     &     nopt,nroots
      character(*), intent(inout) ::
     &     parameters

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(2i4)') nopt,nroots
      else
        read(parameters,'(2i4)') nopt,nroots
      end if

      return
      end

      subroutine evalprop_parameters(rw,parameters,
     &     ndens,rank,env_type)

      implicit none
      
      integer, intent(in) ::
     &     rw
      integer, intent(inout) ::
     &     ndens,rank
      character(*), intent(inout) ::
     &     parameters, env_type

      if (rw.lt.0) then
        parameters(1:len(parameters)) = ' '
        write(parameters,'(2i4,a)') ndens,rank,env_type
      else
        read(parameters,'(2i4,a)') ndens,rank,env_type
      end if

      return
      end

      
