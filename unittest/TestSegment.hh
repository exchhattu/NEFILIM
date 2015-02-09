#ifndef CPPUNITTEST_TESTSEGMENT_HH
#define CPPUNITTEST_TESTSEGMENT_HH

#include<iostream>
#include<string>
#include<vector>

#include<cppunit/TestCase.h>
#include<cppunit/TestFixture.h>
#include<cppunit/extensions/HelperMacros.h>
#include<cppunit/ui/text/TextTestRunner.h>
#include<cppunit/TestResult.h>
#include<cppunit/TestResultCollector.h>
#include<cppunit/TestRunner.h>


// include libraries which have to be tested

#include<FileIO.hh>
#include<Model.hh>
#include<SubModel.hh>
#include<ProfileVector.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

// using namespace CppUnit;
using namespace std;


class TestSegmentProfile: public CppUnit::TestFixture {
  private:
    size_t frag_size_;
    float bin_size_;

    typedef core::pose::Pose Pose;
    typedef ProfileVector<int> IProfileVector;

    Pose pose_;
    string model_id_;
    vector<SubModel> submodels_;
    vector<IProfileVector> model_profiles_;

  public:
    void setup();
    void testEquality();
    void testProfileVectors();

};

#endif

