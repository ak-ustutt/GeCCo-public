*----------------------------------------------------------------------*
      subroutine get_argument_value(
     &     context,argkey,keycount,argcount,
     &     ival,iarr,lval,larr,xval,xarr,string)
*----------------------------------------------------------------------*
*     return dimension (and type) of argument for keyword given 
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*     the keyword must be active (status > 0)
*     the first appearance in history is evaluated, unless count is set
*----------------------------------------------------------------------*
      use keyword_trees,only : 
     &     inp_arg_from_context,reg_arg_from_context,
     &     Node,
     &     atr_val,atr_name,atr_len,
     &     getAttribute
      use parse_input, only :get_argument_dimension_core
      use FoX_common, only : rts

      implicit none
      include 'par_vtypes.h'
      include 'stdunit.h'

      integer,parameter::
     &     ntest=00
      character(len=18),parameter ::
     &     i_am="get_argument_value"

      character, intent(in) ::
     &     context*(*), argkey*(*)
      integer, intent(in), optional ::
     &     keycount,argcount
      logical, intent(out), optional ::
     &     lval, larr(*)
      integer, intent(out), optional ::
     &     ival, iarr(*)
      real(8), intent(out), optional ::
     &     xval, xarr(*)
      character, intent(out), optional ::
     &     string*(*)

      type(Node),pointer::
     &     curarg

      integer ::
     &     iargcount, iargcount_target,
     &     ikeycount, ikeycount_target,
     &     idx, itype, dim, ex, num
      logical ::
     &     succ, dummy

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         if (present(keycount)) write (lulog,*) "keycount:",keycount
         if (present(argcount)) write (lulog,*) "argcount:",argcount
         write (lulog,*) "looking for ",argkey," in context:",context
      end if

      succ=.false.

      ikeycount = 1
      if (present(keycount)) ikeycount = keycount
      ikeycount_target=ikeycount
      
      iargcount = 1
      if (present(argcount)) iargcount = argcount
      iargcount_target=iargcount


      curarg=>inp_arg_from_context(context,argkey,.false.,
     &     ikeycount_target,iargcount_target)

      if (ntest.ge.100) then
         if (associated(curarg))
     &       write (lulog,*)  "active_node in input found"
         if (.not.associated(curarg)) 
     &        write (lulog,*) "active_node not in input found"
      end if

      ikeycount_target=ikeycount
      iargcount_target=iargcount
      if(.not.associated(curarg))
     &     curarg=>reg_arg_from_context(context,argkey,.false.,
     &     ikeycount_target,iargcount_target)


      if (associated(curarg))then
         call get_argument_dimension_core(curarg,dim,
     &        itype,dummy)
         if (ntest.ge.100) then 
            write(lulog,'(" argument has ",i3," elements")') dim 
         end if
      end if 
      select case(itype)
      case (vtyp_log)
         if (.not.(present(lval).or.present(larr)))
     &        call quit(1,i_am,
     &        trim(context)//'->'//trim(argkey)//
     &        'no l-value array present')
         if (present(lval)) 
     &        call rts(getAttribute(curarg,atr_val),lval,
     &        iostat=ex)
         if (present(larr))  
     &        call rts(getAttribute(curarg,atr_val),
     &        larr(1:dim),iostat=ex,num=num) 
         if (ex.le. 0) succ = .true.
      case (vtyp_int)
         if (.not.(present(ival).or.present(iarr)))
     &        call quit(1,i_am,
     &        trim(context)//'->'//trim(argkey)//
     &        'no i-value array present')
         if (present(ival)) then
            call rts(getAttribute(curarg,atr_val),
     &           ival,iostat=ex,num=num)
         end if 
         if (present(iarr)) 
     &        call rts(trim(getAttribute(curarg,atr_val)),
     &        iarr(1:dim),iostat=ex,num=num)
         if (ex.le. 0) succ = .true.
         
      case (vtyp_rl8)
         if (.not.(present(xval).or.present(xarr)))
     &        call quit(1,i_am,
     &        trim(context)//'->'//trim(argkey)//
     &        'no r-value array present')
         if (present(xval)) 
     &        call rts(getAttribute(curarg,atr_val),xval
     &        ,iostat=ex)
         if (present(xarr)) 
     &        call rts(getAttribute(curarg,atr_val),
     &        xarr(1:dim),iostat=ex,num=num) 
         if (ex.le. 0) succ = .true.
      case (vtyp_str)
         if (.not.(present(string)))
     &        call quit(1,'get_argument_value',
     &        trim(context)//'->'//trim(argkey)//
     &        'no r-value array present')
         string = trim(getAttribute(curarg,atr_val))
         succ = .true.
      end select
     

      if (.not.succ)
     &     call quit(1,i_am,
     &     'Could not provide any value for '//trim(context)//
     &     '.'//trim(argkey))


      return
      end
