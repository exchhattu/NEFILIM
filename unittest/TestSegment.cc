#include<iostream>
#include<string>
#include<vector>

// cpp unit test library
// #include<cppunit/TestCase.h>
#include<cppunit/TestFixture.h>
#include<cppunit/extensions/HelperMacros.h>
// #include<cppunit/ui/text/TextTestRunner.h>
// #include<cppunit/TestResult.h>
// #include<cppunit/TestResultCollector.h>
// #include<cppunit/TestRunner.h>


// include libraries which have to be tested

#include<FileIO.hh>
#include<Model.hh>
#include<SubModel.hh>
#include<ProfileVector.hh>

#include<TestSegment.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

using namespace std;


void TestSegmentProfile::setup() {
  FileIO fileio;
  string pdb_file_path = "./testdata/1ctf_atom.pdb"; 
  fileio.read_pose_no_stout(pdb_file_path, model_id_, pose_);
  
  Model model(model_id_, pose_, true);
  model.generate_overlapping_sub_models(pose_, submodels_, 9, 1.0, 0, 30);
  
  for(vector<SubModel>::iterator vit = submodels_.begin(); vit != submodels_.end(); vit++)  {
    model_profiles_.push_back((*vit).get_prof_vect());
  }
}

void TestSegmentProfile::testEquality() {
  CPPUNIT_ASSERT(model_id_=="1ctf_atom");
  CPPUNIT_ASSERT(submodels_.size() > 0);
}


void TestSegmentProfile::testProfileVectors() {
  CPPUNIT_ASSERT(submodels_.size() > 0);
}
