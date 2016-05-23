*------------------------------------------------------------------------*
!>    sets blocks of a matrix element list to val
!!
!!    
!!    @param[inout] mel matrix element list to be worked on
!!    @param[in] occ_descr string which defines the blocks to be set to val from
!!    @param[in] val  use is dependend on mode
!!    @param[in] mode string that defines the action
!!                 "SET"  sets all elements of the specified blocks to val
!!                 "SCALE" scales all elements of the specified blocks by val
!!                 "SHIFT" shifts all elements of the specified blocks to val
      subroutine modify_blks(mel,occ_descr,val,mode)
*------------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'


      integer, parameter ::
     &     ntest = 00
      integer ::
     &     iprint

      character(len=11) ::
     &     i_am = 'modify_blks'
      
      type(me_list), intent(inout) ::
     &     mel
      real(8),intent(in) ::
     &     val
      character(len=*), intent(in) ::
     &     occ_descr,mode

      real(8),pointer::
     &     buffer(:)
      logical ::
     &     closeit
      integer ::
     &     ifree
      integer ::
     &     nblkmax,idisc_off,idxst,nbuff,idxnd,nblk
      integer ::
     &     iocc,nocc
!, len_op, nblk, nblkmax, nbuff, idxst, idxnd,
!     &     idum, iblk, idisc_off, iset, idx

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop

!     variables for process_occ_descr
!     maxdef: maximum number of possible occupation classes
!     maxlist:  number of blocks
!     nocc_cmp: number of created blocks(blocks or occ_cls?)
!     occ_def_cmp: the ca_strings in touple form
        integer, parameter ::
     &     maxdef = 512
        integer ::
     &       maxlist,nocc_cmp,iocc_cmp
        integer ::
     &     occ_def_cmp(ngastp,2,maxdef)


!     number of joined vertices
        integer::
     &       njoined,len_op

!     indices for the original operator
!     idx_blk_* are starting and ending index of the relevant part in the buffer
!     blk_* are starting and ending index of the relevant occupaiton classes
!     blk_no is the number current block; nblk ist the total number of blocks
        integer::
     &       idx_blk_nd,idx_blk_st,idx,
     &       blk_nd,blk_st,blk_no,occ 


!     pointer to the block definitions of the original operator
        integer,pointer ::
     &       occ_def(:,:,:)

!     the function to check wether two blocks are identical
        logical, external ::
     &       iocc_equal_n






      iprint=max(iprlvl,ntest)



      if (iprint.ge.100) then
        call write_title(lulog,wst_dbg_subr,'modify_blocks')
        write(lulog,*) 'occ_descr:',occ_descr
        write(lulog,*) 'val :',val
      end if


      ifree = mem_setmark('modify_blocks')

      ffop => mel%fhand
      op => mel%op

      if (iprint.ge.100) then
         write(lulog,*) 'mel :',trim(mel%label)
         write(lulog,*) 'op :',trim(op%name)
      end if


      if (.not.associated(ffop))
     &     call quit(1,i_am,'No file assigned to list: '//
     &     trim(mel%label))

      if (ffop%unit.le.0) then
        call file_open(ffop)
        closeit = .true.
      else
        closeit = .false.
      end if

!     maxnumber of records fitting into buffer
      nblkmax = ifree/ffop%reclen


      if (nblkmax.le.0) then
        write(lulog,*) 'free memory (words):  ',ifree
        write(lulog,*) 'block length (words): ',ffop%reclen
        call quit(1,i_am,'not even 1 record fits into memory?')
      end if

!     offset on file (if more than one instance of operator exists(e.g.: multiple roots were calculated))
      idisc_off = ffop%length_of_record*(ffop%current_record-1)

!     index-start right after the offset
      idxst = idisc_off+1


        occ_def=>op%ihpvca_occ

        njoined=op%njoined


        maxlist = maxdef/njoined

!     Convert the string of occupation classes into an array 
        call process_occ_descr(occ_def_cmp,nocc_cmp,
     &       occ_descr,njoined,maxlist)

 

      if (.not.ffop%buffered) then

        len_op = mel%len_op

        nblk = min((len_op-1)/ffop%reclen + 1,nblkmax)

!     Bufferlength is either the full operator or a length that can hold the maximum amount of 
!      records fitting into the free space
        nbuff = min(len_op,nblk*ffop%reclen)

        if (nbuff .lt. len_op) 
     &       call quit(1,i_am,"not prepared ")

        ifree = mem_alloc_real(buffer,nbuff,'buffer')
        

        buffer(1:nbuff) = 0d0


!     taktik: 

!     while idxst (starting index of block on file)
!     is lower than end of operator:
!     go on
!        while (idxst.lt.idisc_off-1+len_op)

!     idxnd=endindex on file.
!     either offset+length of op 
!     or (in case of iteration over blocks ) start of this block plus  
        idxnd = min(idisc_off-1+len_op,
     &               idxst-1+nbuff)


! nocc number of full (all vertices) occupation classes
        nocc=op%n_occ_cls

      if (iprint.ge.100) then
        write(lulog,*) 'nbuff:  ',nbuff
        write(lulog,*) 'idxst:  ',idxst
        write(lulog,*) 'idxnd:  ',idxnd   
        write(lulog,*) 'nocc:  ',nocc
        write(lulog,*) 'nocc_cmp:  ',nocc_cmp   
      end if
        call get_vec(ffop,buffer,idxst,idxnd)
!     Looping over all blocks of the MEL 
        do iocc=1,nocc

!     Looping over all blocks of the description
           cmp_loop:do iocc_cmp=1,nocc_cmp

              if (iocc_equal_n(
     &             occ_def_cmp(1,1,(iocc_cmp-1)*njoined+1),.false.,
     &             occ_def(1,1,(iocc-1)*njoined+1),.false.,njoined))
     &             then
                idx_blk_st=mel%off_op_occ(iocc)+1
                idx_blk_nd=idx_blk_st+mel%len_op_occ(iocc)-1

                if (iprint.ge.100) then
                   write(lulog,*) "setting the following block to val"                   
                   write(lulog,*) 'idx_blk_st:  ',idx_blk_st
                   write(lulog,*) 'idx_blk_nd:  ',idx_blk_nd
                   write(lulog,*) 'MODE:',mode
                end if

                select case(trim(mode))
                case("SET")
                   buffer(idx_blk_st:idx_blk_nd)=val
                case("SHIFT")
                   buffer(idx_blk_st:idx_blk_nd)=
     &                  buffer(idx_blk_st:idx_blk_nd)+val
                case("SCALE")
                   buffer(idx_blk_st:idx_blk_nd)=
     &                  buffer(idx_blk_st:idx_blk_nd)*val
                case default
                   call quit(1,i_am,'unknown mode')
                end select
                
                exit cmp_loop ! assumes any occupation only occurs once in the 
                              ! comparison loop
                              ! so after match we can search for the next occupation
              end if
           end do cmp_loop
        end do
        call put_vec(ffop,buffer,idxst,idxnd)
      else
         call  quit(1,i_am,'Not prepared for buffered vector')
      end if
!      end while
      call touch_file_rec(ffop)

      if (closeit)
     &     call file_close_keep(ffop)

      ifree = mem_flushmark()
      return

      end




!-----------------------------------------------------------------
!> Depreciated, retained for backward-compatibility
!
      subroutine set_blks(mel,occ_descr,val)
!-----------------------------------------------------------------
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'

      type(me_list), intent(inout) ::
     &     mel
      real(8),intent(in) ::
     &     val
      character(len=*), intent(in) ::
     &     occ_descr
      call modify_blks(mel,occ_descr,val,"SET")

      end 
