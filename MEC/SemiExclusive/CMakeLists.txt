cmake_minimum_required(VERSION 3.15)

# Project name and language
project(em_evgen_project LANGUAGES Fortran)

# Use mpif90 as the Fortran compiler (same as your script)
set(CMAKE_Fortran_COMPILER mpif90)

# Add your executable and its source files
add_executable(semiexc
    progress_bar.f90
    mathtool.f90
    mympi.f90
    SpectralFunction.f90
    event_module.f90
    isospin_operators.f90
    currents_opt_v3_an1.f90
    integrated_semidiff_xsec_v5_an1.f90
    em_evgen_an1.f90
)

# Add the same compile flags as in compile_em.sh
target_compile_options(semiexc PRIVATE
    -O3
    -fallow-argument-mismatch
    -march=native
)

# --- Data files that should live next to the executable in the build dir ---
set(DATA_FILES
    C12_rho2b_qQ.txt
    n2b_c12_new_fmt.dat
    input_xsec_tot.in
    rho_1.dat
)

set(COPIED_FILES)

foreach(f IN LISTS DATA_FILES)
    set(src "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
    set(dst "${CMAKE_CURRENT_BINARY_DIR}/${f}")

    add_custom_command(
        OUTPUT ${dst}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${src} ${dst}
        DEPENDS ${src}
        COMMENT "Copying ${f}"
    )

    list(APPEND COPIED_FILES ${dst})
endforeach()

add_custom_target(copy_data ALL
    DEPENDS ${COPIED_FILES}
)

add_dependencies(semiexc copy_data)
