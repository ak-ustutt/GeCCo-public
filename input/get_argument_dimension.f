


*----------------------------------------------------------------------*
      subroutine get_argument_dimension(dim,
     &     context,argkey,keycount,argcount,type)
*----------------------------------------------------------------------*
*     return dimension (and type) of argument for keyword given 
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*     the keyword must be active (status > 0)
*     the first appearance in history is evaluated, unless count is set
*----------------------------------------------------------------------*

      use keyword_trees,only :tree_get_arg_from_context,
     &     fetch_input_keyword_tree,fetch_registry_keyword_tree,
     &     Node,tree_t
      use parse_input, only :get_argument_dimension_core

      implicit none
      include 'par_vtypes.h'
      include 'stdunit.h'

      integer,parameter::
     &     ntest=00
      character(len=22),parameter ::
     &     i_am="get_argument_dimension"
      character, intent(in) ::
     &     context*(*), argkey*(*)
      integer, intent(in), optional ::
     &     keycount,argcount
      integer, intent(out),optional::
     &     type


      type(Node), pointer ::
     &     curarg
      type(tree_t)::
     &     registry,input
      character ::
     &     curcontext*1024
      integer ::
     &     iargcount, iargcount_target,
     &     ikeycount, ikeycount_target,
     &     idx, itype, dim
      logical ::
     &     succ



      input=fetch_input_keyword_tree()
      registry=fetch_registry_keyword_tree()


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


      curarg=>tree_get_arg_from_context(input,context,argkey,.false.,
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
     &     curarg=>tree_get_arg_from_context(registry,context,
     &     argkey,.false.,
     &     ikeycount_target,iargcount_target)

      if (associated(curarg))then
         call get_argument_dimension_core(curarg,dim,
     &        itype,succ)
         if (ntest.ge.100) then 
            write(lulog,'(" argument has ",i3," elements")') dim 
         end if
      end if 
            
      if (.not.succ)
     &     call quit(1,i_am,
     &     'Could not provide any value for '//trim(context)//
     &     '.'//trim(argkey))

      if (present(type))type=itype

      return
      end
