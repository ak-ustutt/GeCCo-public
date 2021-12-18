#/bin/bash
function print_msg {
          token=$1
          message=$2
          shift 2
          if [ -f $1 ] 
	  then
              files=$(grep -l $token $@ 2>/dev/null )
	      if [ $? -eq 0 ]
	      then
		  echo "============================================="; 
		  echo $message ; 
		  echo "---------------------------------------------"; 
		  for file in $files
		  do
		      echo ${file%.err} ;
		  done
		  echo "============================================="; 
	      fi
	  fi
}

function print_all {
	print_msg "check_OK" " Successful tests:" $@ 
	print_msg "check_failed" "Failed tests:" $@
	print_msg "check_buggy" "The following tests have buggy checks:" $@
	print_msg "setup_buggy" "The following tests had buggy inputs:" $@



}

function print_summary {
    if [ -f $1 ] 
    then
	if [ $(grep -e "failed" -e "buggy" $@ &> /dev/null ; echo $?) -eq 1 ] ; then 
	    echo "A lucky day: all tests seem to run correctly!" 
	else  
	    echo "Oh no: some tests give errors, see above!" 
	    echo "!!! ERRORS !!!" 
	fi
    fi  #else it is probably just *.err
}
