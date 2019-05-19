# Tunguska - 64k Demo by Team210 at Solskogen 2019
# Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)

project(gpn19-visuals)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(WIN32)
    find_program(KKRUNCHY_K7 kkrunchy_k7.exe HINTS ${PROJECT_BINARY_DIR})
    
    if(${CTRL_ALT_TEST_MINIFIER})
        find_program(MINIFIER shader_minifier.exe HINTS ${PROJECT_BINARY_DIR})
        add_custom_target(minification
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        )
    else()
        add_custom_target(minification
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/minification/
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/load.h" "${PROJECT_SOURCE_DIR}/gfx/load.frag" --no-minification
            COMMAND python shader-compressor.py -o "${PROJECT_SOURCE_DIR}/gfx/post.h" "${PROJECT_SOURCE_DIR}/gfx/post.frag" --no-minification
        )
    endif()
    
    add_custom_target(texture_packing
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/font/
        COMMAND python pack.py
    )
    
    add_custom_target(symbolization
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/gfx/
        COMMAND python symbolize.py decayingfactory.frag fogforest.frag -o symbols.h
    )
    
    add_definitions(/Gy /O1 /fp:fast /GR- /GS- /MT)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ENTRY:demo /SUBSYSTEM:Windows  /OPT:REF  /OPT:ICF=10 /VERBOSE /NODEFAULTLIB") 
    add_executable(gpn19-visuals t.c)
    
    # for debug output: User32.lib ucrt.lib msvcrt.lib
    target_link_libraries(gpn19-visuals opengl32.lib Winmm.lib User32.lib ucrt.lib msvcrt.lib)
    
    add_dependencies(gpn19-visuals minification)
    add_dependencies(gpn19-visuals texture_packing)
    add_dependencies(gpn19-visuals symbolization)
    add_custom_command(TARGET gpn19-visuals COMMAND ${KKRUNCHY_K7} --best ${CMAKE_CURRENT_BINARY_DIR}/Release/gpn19-visuals.exe --out ${CMAKE_CURRENT_BINARY_DIR}/Release/gpn19-visualsc.exe POST_BUILD)
endif()

