cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

project(Ensitsp)
enable_testing()
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall -Werror -O3 -std=gnu99")

#########
# Gestion des logins
#########
# Vous devez editer la ligne suivante pour y mettre vos logins
#########
set(LOGINS login1 login2 login3)
list(SORT LOGINS)

#####
# Détection de la bibliothèque de threads et de la bibliothèque de math
#####
find_package(Threads)
find_library(MISC_LIBRARY NAMES m libm rt librt)


##
# Si vous utilisez d'autres fichiers .c, en plus de ceux listés,
# il faut les ajouter ici
##
add_executable(ensitsp src/tsp-genmap.c  src/tsp-job.c  src/tsp-main-para.c  src/tsp-print.c  src/tsp-tsp.c src/tsp-lp.c)
target_link_libraries(ensitsp ${CMAKE_THREAD_LIBS_INIT} m rt)

##
# Programme de test
##
add_test(Tests ../tests/validite.pl)
add_custom_target(check ../tests/validite.pl)

#add_test(MesuresDePerf ../mesures/mesures.pl)


##
# Construction de l'archive
##
string(REPLACE ";" "-" LOGINS_SANS_POINTVIRGULE "${LOGINS}")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH ${LOGINS_SANS_POINTVIRGULE})
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES
  "~$"
  "\\\\.o$"
  "^${PROJECT_SOURCE_DIR}/build/"
  "^${PROJECT_SOURCE_DIR}/.git/"
)
include(CPack)
