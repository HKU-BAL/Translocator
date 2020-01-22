cmake_minimum_required(VERSION 2.8)
project(Translocator)

find_package (ZLIB)
#find_package (SeqAn REQUIRED)


#configure_file( Version.h.in ${CMAKE_SOURCE_DIR}/src/Version.h )

#set(SEQAN_INCLUDE_DIRS "/data2/junwenwang/m204333/Project/seqan-v2.2.0/include")

include_directories (../lib/bamtools-2.3.0/src)
include_directories(../lib/tclap-1.2.1/include)
include_directories(../lib/seqan-library-2.2.0/include)

find_package (Threads)

IF(STATIC)
	SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
	SET(BUILD_SHARED_LIBRARIES OFF)
	SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
ENDIF()


# Add CXX flags found by find_package (SeqAn).
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

add_library(libminimap2
		../lib/minimap2/kthread.c
		../lib/minimap2/kalloc.c
		../lib/minimap2/ksw2_extz2_sse.c
		../lib/minimap2/ksw2_extd2_sse.c
		../lib/minimap2/ksw2_exts2_sse.c
		../lib/minimap2/misc.c
		../lib/minimap2/bseq.c
		../lib/minimap2/sketch.c
		../lib/minimap2/sdust.c
		../lib/minimap2/options.c
		../lib/minimap2/index.c
		../lib/minimap2/chain.c
		../lib/minimap2/align.c
		../lib/minimap2/hit.c
		../lib/minimap2/map.c
		../lib/minimap2/format.c
		../lib/minimap2/pe.c
		../lib/minimap2/esterr.c
		../lib/minimap2/ksw2_ll_sse.c
		../lib/minimap2/splitidx.c
		)

target_link_libraries(libminimap2 zlibstatic)
target_link_libraries (libminimap2 ${CMAKE_THREAD_LIBS_INIT})

add_executable(translocator
                    tree/Breakpoint_Tree.cpp
                    Genotyper/Genotyper.cpp
					Alignment.cpp
					BamParser.cpp
		Translocator.cpp
					Ignore_Regions.cpp
					tree/Intervall_bed.cpp
					sub/Detect_Breakpoints.cpp
					sub/Breakpoint.cpp
					tree/IntervallTree.cpp
					tree/IntervallList.cpp
					realign/SWCPU.cpp
					realign/Realign.cpp
					print/VCFPrinter.cpp
					print/BedpePrinter.cpp
					print/IPrinter.cpp
					tree/BinTree.cpp
					print/NGMPrinter.cpp
					plane-sweep/PlaneSweep_slim.cpp
					cluster/Cluster_SVs.cpp
					force_calling/Force_calling.cpp
					force_calling/VCF_parser.cpp
        			sub/Realign.cpp)

TARGET_LINK_LIBRARIES(translocator BamTools-static)
TARGET_LINK_LIBRARIES(translocator zlibstatic)
TARGET_LINK_LIBRARIES(translocator libminimap2)

add_executable(translocator-debug
                    tree/Breakpoint_Tree.cpp
                    Genotyper/Genotyper.cpp
					Alignment.cpp
					BamParser.cpp
		Translocator.cpp
					Ignore_Regions.cpp
					tree/Intervall_bed.cpp
					tree/IntervallList.cpp
					sub/Detect_Breakpoints.cpp
					sub/Breakpoint.cpp
					tree/IntervallTree.cpp
					realign/SWCPU.cpp
					realign/Realign.cpp
					print/VCFPrinter.cpp
					print/BedpePrinter.cpp
					print/IPrinter.cpp
					tree/BinTree.cpp
					print/NGMPrinter.cpp
					plane-sweep/PlaneSweep_slim.cpp
					cluster/Cluster_SVs.cpp
					force_calling/Force_calling.cpp
					force_calling/VCF_parser.cpp
        			sub/Realign.cpp)

SET_TARGET_PROPERTIES(translocator-debug PROPERTIES COMPILE_FLAGS "-g3 -O0")



#target_link_libraries(translocator-debug pthread)
TARGET_LINK_LIBRARIES(translocator-debug BamTools-static)
TARGET_LINK_LIBRARIES(translocator-debug zlibstatic)
TARGET_LINK_LIBRARIES(translocator-debug libminimap2)





