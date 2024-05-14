<<<<<<< HEAD
<<<<<<< HEAD
# Install script for directory: C:/Users/diana_9kbxrpt/Downloads/uni/2ano_2semestre/DA/projeto2/TravellingSalesman
=======
# Install script for directory: C:/Users/teres/Downloads/TravellingSalesman
>>>>>>> 54538a3cad605198209c6b11da5db88c7e04d2a2
=======
# Install script for directory: C:/Users/teres/Downloads/TravellingSalesman
>>>>>>> 4cc1d21c3d11c064b06c972fe624693abc6edffb

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/TravellingSalesman")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "C:/Program Files/JetBrains/CLion 2023.2.2/bin/mingw/bin/objdump.exe")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
<<<<<<< HEAD
<<<<<<< HEAD
file(WRITE "C:/Users/diana_9kbxrpt/Downloads/uni/2ano_2semestre/DA/projeto2/TravellingSalesman/cmake-build-debug/${CMAKE_INSTALL_MANIFEST}"
=======
file(WRITE "C:/Users/teres/Downloads/TravellingSalesman/cmake-build-debug/${CMAKE_INSTALL_MANIFEST}"
>>>>>>> 54538a3cad605198209c6b11da5db88c7e04d2a2
=======
file(WRITE "C:/Users/teres/Downloads/TravellingSalesman/cmake-build-debug/${CMAKE_INSTALL_MANIFEST}"
>>>>>>> 4cc1d21c3d11c064b06c972fe624693abc6edffb
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
