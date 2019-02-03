#!/usr/bin/python

import commands, os, SConstructLib 

common_opts = " -W  -Wall "

cpp_flags  = " -O3 " + common_opts 

system_infos = os.uname()
current_os   = system_infos[0]
host_name    = system_infos[1]
current_dir  = os.getcwd()
rosetta_home = ''

#for bragg
rosetta_home = ### ROSETTA_INSTALL_PATH 
cppunit_top  = ### CPPUNIT_PATH 

rosetta_shared_lib_path = SConstructLib.get_rosetta_shared_lib_path(rosetta_home, "single")
rosetta_platform_path   = SConstructLib.get_platform_path(rosetta_home, current_os)

rosetta_lib=os.path.join(rosetta_home, "/rosetta_source/build/src/release/linux/2.6/64/x86/gcc")

VariantDir('build', 'src', duplicate=0)

src_files = [
              './src/FileIO.cc',
              './src/LocalProfile.cc', 
              './src/init_rosetta.cc',
              './src/RosettaFragment.cc',
              './src/FragmentCluster.cc',
              './src/DistRange.cc',
              './src/DistMatrix.cc',
              './src/Singleton.cc',
              './src/rmsd.cc',
              './src/qcprot.cc',
              './src/SimpPDB.cc',
              './src/Stru.cc',
              './src/Residue.cc', 
              './src/FragmentManager.cc',
              './src/SegmentProfile.cc'
             ]

libraries = [ 'ObjexxFCL', 'core', 'utility', 'protocols', 'numeric', 'devel', 
              'dl', 'cppunit']

header_files = [ current_dir + '/src/',
                 current_dir + '/rosetta_wrapper/',
                 current_dir + '/unittest/',
                 rosetta_home + '/src',
                 rosetta_home + '/external/include',
                 rosetta_home + '/external/boost_1_38_0',
                 rosetta_home + rosetta_platform_path,
                 cppunit_top + "/include/"
               ]

libpaths  = [ rosetta_lib, current_dir + '/src/', cppunit_top + "/lib/" ]

env = Environment(CXX = "g++",
                  CCFLAGS = cpp_flags,
                  LIBS    = libraries,
                  CPPPATH = header_files,
                  LIBPATH = libpaths)  

# env.Program(target = "SegmentProfiler", source  = src_files)

del(src_files[len(src_files) - 1])
src_files.append("./src/FragmentGenerator.cc")
env.Program(target = "NEFILIM", source  = src_files)
