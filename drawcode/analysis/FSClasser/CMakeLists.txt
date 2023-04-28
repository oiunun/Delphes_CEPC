########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( FSClasser )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )
#SET( ${FastJet_HOME} /besfs/groups/higgs/users/lig/analysis/LCFIplus/fastjet/install )
#SET( ${FastJet_LIBRARIES} ${FastJet_HOME}/lib) 
#SET( ${FastJet_INCLUDE_DIRS} ${FastJet_HOME}/include) 


### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )
FIND_PACKAGE( ROOT REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )



FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
FIND_PACKAGE( MarlinUtil 1.4.0 REQUIRED ) # minimum required Marlin version
FIND_PACKAGE( MarlinKinfit 0.1.2 REQUIRED ) # minimum required Marlin version
#FIND_PACKAGE( LCFIPlus 0.5.2 REQUIRED ) # minimum required Marlin version
#FIND_PACKAGE( FastJet 2.4.2 REQUIRED ) # minimum required Marlin version
#FIND_PACKAGE( FastJet 3.3.4 REQUIRED ) # minimum required Marlin version
#
#INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIRS} )
#INCLUDE_DIRECTORIES( ${MarlinUtil_INCLUDE} )
#INCLUDE_DIRECTORIES( ${MarlinKinfit_INCLUDE_DIRS} )
#INCLUDE_DIRECTORIES( ${FastJet_INCLUDE_DIRS} )
#INCLUDE_DIRECTORIES( ${LCFIPlus_INCLUDE} )
#INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )


##
FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )
#FOREACH( pkg Marlin MarlinUtil MarlinKinfit LCIO FastJet ROOT GSL )
FOREACH( pkg Marlin MarlinUtil MarlinKinfit LCIO ROOT GSL )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()
LINK_LIBRARIES(     /besfs5/groups/higgs/users/lig/analysis/LCFIplus/fastjet/install/lib/libfastjet.so  )
LINK_LIBRARIES(     /besfs5/groups/higgs/users/lig/analysis/LCFIplus/fastjet/install/lib/libsiscone.so  )
LINK_LIBRARIES(     /besfs5/groups/higgs/users/lig/analysis/LCFIplus/fastjet/install/lib/libfastjetplugins.so )

##

#LINK_LIBRARIES( ${Marlin_LIBRARIES} )
#ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

#LINK_LIBRARIES( ${MarlinKinfit_LIBRARIES} )
LINK_LIBRARIES( ${FastJet_LIBRARIES} )


# optional package
#FIND_PACKAGE( AIDA )
#IF( AIDA_FOUND )
#    INCLUDE_DIRECTORIES( ${AIDA_INCLUDE_DIRS} )
#    LINK_LIBRARIES( ${AIDA_LIBRARIES} )
#    ADD_DEFINITIONS( "-DMARLIN_USE_AIDA" )
#    MESSAGE( STATUS "AIDA -- found" )
#ELSE()
#    MESSAGE( STATUS "AIDA -- not found" )
#ENDIF()

#FIND_PACKAGE( GSL 1.12 REQUIRED )
#INCLUDE_DIRECTORIES( ${GSL_INCLUDE_DIRS} )
#LINK_LIBRARIES( ${GSL_LIBRARIES} )



### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic -Wno-unused-variable -Wwrite-strings -Wno-unused-but-set-variable " )
ADD_DEFINITIONS( "-Wno-long-long -fvar-tracking-assignments-toggle -std=c++11" )

# include directories
INCLUDE_DIRECTORIES( ./include )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

#SET( Marlin_DEPENDS_INCLUDE_DIRS /cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6/cepcsoft/0.1.0-b7/Framework/LCIO/02-04-03/include;/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6/cepcsoft/0.1.0-b7/Simulation/GEAR/01-04/include;/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6/cepcsoft/0.1.0-b7/Simulation/GEAR/01-04/include/gear;/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6/external/CLHEP/2.1.3.1/install/include;/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6/cepcsoft/0.1.0-b7/Framework/ILCUTIL/01-01/include;/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6/cepcsoft/0.1.0-b7/Analysis/MarlinKinfit/00-01-02/include;/cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6/cepcsoft/0.1.0-b7/Service/FastJet/2.4.2/install/include )
