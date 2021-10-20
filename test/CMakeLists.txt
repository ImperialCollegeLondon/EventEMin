add_executable(example_6dof main_6dof.cpp)
add_executable(example_homography main_homography.cpp)
add_executable(example_rotation main_rotation.cpp)
add_executable(example_translation2d main_translation2d.cpp)

add_executable(example_test_sequence main_test_sequence.cpp)

target_link_libraries(example_6dof ${LIB_NAME}_lib)
target_link_libraries(example_homography ${LIB_NAME}_lib)
target_link_libraries(example_rotation ${LIB_NAME}_lib)
target_link_libraries(example_translation2d ${LIB_NAME}_lib)

target_link_libraries(example_test_sequence ${LIB_NAME}_lib)
