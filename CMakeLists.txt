project (C-GDY-SAC)

add_executable(main.o src/main.c src/file.c src/parse_mol.c)
set (CMAKE_C_FLAGS "-Wall -lpcre2-8 -I/opt/homebrew/include -g")
set (CC gcc-11)
