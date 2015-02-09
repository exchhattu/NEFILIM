#ifndef CPPUNITTEST_TESTSUBSTRUCTURE_HH
#define CPPUNITTEST_TESTSUBSTRUCTURE_HH

#include<iostream>
#include<string>
#include<vector>
#include<map>

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
#include<ProfileVectorImpl.hh>
#include<Distance.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

// using namespace CppUnit;
using namespace std;


class TestSubStructure: public CppUnit::TestFixture {
  private:
    size_t frag_size_;
    Distance distance_;
    string models_path_;
    string ref_coord_path_;
    string models_coords_path_;

    ProfileVectorImpl prof_vector_impl_;

    map<string, float> lsqkab_rmsd_;
    map<string, float> computed_rmsd_;
    map<string, float> known_submodels_;

    map<string, map<string, float> > found_fragments_;

    vector<SubModel> submodels_;

    void read_carmsd(const char* input_file_path, map<string, float> &ca_rmsds);

  public:
    void setup();

    void test_CARMSD_Equality(); 

    void setup_search(); 

    void test_query_to_ref_submodels(); 

    void test_found_fragments(); 

    void test_q_score(); 

    void compute_carmsd_query_to_database(); 
};

#endif

