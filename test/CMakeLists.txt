# Function to add new executables
function(add_new_executable ARG)
  add_executable(example_${ARG} main_${ARG}.cpp)
  target_link_libraries(example_${ARG} ${LIB_NAME})
endfunction()

if(${LIB_NAME}_BATCH_MODE)
  add_new_executable(6dof)
  add_new_executable(homography)
  add_new_executable(rotation)
  add_new_executable(translation2d)

  add_new_executable(test_sequence)
endif()

if(${LIB_NAME}_INCREMENTAL_MODE)
  add_new_executable(incremental_6dof)
  add_new_executable(incremental_rotation)
  add_new_executable(incremental_translation2d)

  add_new_executable(incremental_test_sequence)
endif()
