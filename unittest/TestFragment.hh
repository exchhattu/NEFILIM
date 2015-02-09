#ifndef CPPUNITTEST_FRAGMENT_HH
#define CPPUNITTEST_FRAGMENT_HH

#include<vector>
#include<iostream>

#include<cppunit/TestCase.h>
#include<cppunit/TestAssert.h>
#include<cppunit/TestFixture.h>
#include<cppunit/extensions/HelperMacros.h>
#include<cppunit/ui/text/TextTestRunner.h>
#include<cppunit/TestResult.h>
#include<cppunit/TestResultCollector.h>
#include<cppunit/TestRunner.h>

#include<FileIO.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>


using namespace std;
typedef core::pose::Pose Pose;

class TestBruteFragment: public CppUnit::TestFixture {

  private:

    Pose model_pose;
    Pose ref_pose;

    string model_model_id;
    string ref_model_id; 

    vector<float> carmsds_rsd_one; 

  public:

    void setup(const string modelpath, const string refpath); 

    void testFileOpen();

    void testLocalRMSDCal(); 

};

#endif
