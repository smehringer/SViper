add_test(NAME sviper_invoke_help_page 
	 COMMAND ../sviper -h)
add_test(NAME sviper_basic_usage 
	COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/sviper_test.sh" "../" "/tmp/")
