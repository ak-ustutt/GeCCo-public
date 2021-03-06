*------------------------------------------------------------------------*
*	info about file
*------------------------------------------------------------------------*
      integer, parameter ::
     &     ftyp_da_unf = 1,
     &     ftyp_sq_unf = 2,
     &     ftyp_sq_frm = 3,
     &     ftyp_st_unf = 4
*------------------------------------------------------------------------*
      integer, parameter ::
     &     maxfilnam = 512 
      type filinf
         character*(maxfilnam) ::
     &     name
         integer ::
     &     unit 
         integer ::
     &     type, reclen
         logical ::
     &     buffered
         integer ::
     &     nbuffer
         integer, pointer ::
     &     idxrec(:)
         integer, pointer ::
     &     incore(:)
         real(8), pointer ::
     &     buffer(:)
         ! alternative buffering mechanism: we need an ID
         ! cf. memman.f (mem_buffers)
         integer ::
     &     buf_id
         ! logical "records" for repeated data-structures,
	 ! see assign_file_to_op() for an example
	 integer ::
     &       active_records(2),
     &       current_record,
     &       length_of_record  ! length of that log. record in r8-words
	 ! offset in a "superfile", if several file_handle are
	 ! assigned the same unit (not yet used)
	 integer ::
     &       recoff_superfile
         ! modification time
	 integer(8), pointer ::
     &        last_mod(:)
      end type filinf
*------------------------------------------------------------------------*
*	name, unit: self-explaining, unit == -1 means not assigned
*       type: 1 direct access
*             2 sequential unformatted
*             3 sequential formatted
*       reclen: record length for direct access
*       buffered: set to .false. by default, if set to .true., the
*             behaviour is "user-defined", i.e. some of the below
*             array should be set to the appropriate values
*	idxrec: an index array for direct access files ("user defined")
*       incore: an array for indicating that certain parts of the
*               file are incore ("user defined", e.g. blocks of operator
*               or records of the file ...)
*       buffer: a buffer for incore parts of the file contents
*------------------------------------------------------------------------*
*	a note on buffering for "operators":
*	incore(iocc_cls) holds info whether the occupation iocc_cls
*	is held incore (in the buffer) (if > 0), the number is set
*       equal the length of the occ. class
*       idxrec will then contain the offset of the current block
*       
*	as the length increases with iocc_cls, we currently will hold 
*       the first few occupations (say 1 to 4 for the Hamilton) incore, 
*       so the addressing array on the operator structure can be used 
*       for "buffer" as well
*------------------------------------------------------------------------*
