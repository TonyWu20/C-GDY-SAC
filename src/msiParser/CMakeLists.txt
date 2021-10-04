# src/msiParser/CMakeLists.txt

add_library(msiParser parse_mol.c)
set (CMAKE_C_FLAGS "-lpcre2 -I/opt/homebrew/include")
target_include_directories(msiParser PUBLIC $(CMAKE_CURRENT_SOURCE_DIR))
