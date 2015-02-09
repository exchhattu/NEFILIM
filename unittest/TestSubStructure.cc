#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>

// cpp unit test library
#include<cppunit/TestFixture.h>
#include<cppunit/extensions/HelperMacros.h>

// include libraries which have to be tested
#include<FileIO.hh>
#include<Model.hh>
#include<SubModel.hh>
#include<ProfileVector.hh>
#include<ProfileVectorImpl.hh>
#include<Distance.hh>

#include<TestSubStructure.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

using namespace std;

void TestSubStructure::setup() {
  models_path_ = "./models.lst";
  ProfileVectorImpl prof_vector_impl(models_path_, 9, 0.5, false);
  prof_vector_impl_ = prof_vector_impl;
  prof_vector_impl_.filter_intra_submodels(0.0);

  ref_coord_path_ = "./reference.out";
  distance_.set_reference_submodels(ref_coord_path_.c_str());
  // since the reference sets are already in sorted order.
  // models_coords_path_ = "./coordinates.out";
  // distance_.compute_and_dump_carmsd(models_coords_path_);
  read_carmsd("./3chy_a_15_1_lsqkab.out", lsqkab_rmsd_); 
  read_carmsd("./3chy_a_15_1.out", computed_rmsd_); 
}

void TestSubStructure::read_carmsd(const char* input_file_path, map<string, float> &ca_rmsds) {
  ifstream input_file(input_file_path, ifstream::in); 
  string line; 
  try {
    if(input_file.is_open()) {
      while(getline(input_file, line)) {
        stringstream sstream_(line);
        string submodel_id; 
        float ca_rmsd = 0.0;
        sstream_>>submodel_id>>ca_rmsd;
        ca_rmsds[submodel_id] = ca_rmsd;
      }
      input_file.close();
    }
    else {
      cout<<"Fatal:"<<input_file_path<<" cannot be opened."<<endl;
      cout<<__FILE__<<" line number: "<<__LINE__<<endl;
      exit(1);
    }
  }
  catch(exception &e) {
    cout<<"Fatal: "<<e.what()<<endl;
    cout<<__FILE__<<" line number: "<<__LINE__<<endl;
    exit(1);
  }
}

void TestSubStructure::test_CARMSD_Equality() {
  CPPUNIT_ASSERT(computed_rmsd_.size() > 0 && lsqkab_rmsd_.size() > 0); 
  CPPUNIT_ASSERT(lsqkab_rmsd_.size() == computed_rmsd_.size());
  map<string, float>::iterator mit1 = lsqkab_rmsd_.begin();
  map<string, float>::iterator mit2 = computed_rmsd_.begin();
  while( mit1 != lsqkab_rmsd_.end() && mit2 != computed_rmsd_.end()) {
    float carmsd1 = mit1->second, carmsd2 = mit2->second;
    string submodel_id1 = mit1->first, submodel_id2 = mit2->first;
    transform(submodel_id1.begin(), submodel_id1.end(), submodel_id1.begin(), ::tolower);
    transform(submodel_id2.begin(), submodel_id2.end(), submodel_id2.begin(), ::tolower);
    CPPUNIT_ASSERT(submodel_id1==submodel_id2);
    CPPUNIT_ASSERT(abs(carmsd1-carmsd2) < 0.010);
    mit1++; mit2++;
  }
}

void TestSubStructure::setup_search() {
  string query_model_path = "./3chy_a_2_1_query.pdb";
  prof_vector_impl_.filter_intra_submodels(submodels_, query_model_path, 0.0); 
  distance_.compute_distance_to_references_n_search(submodels_, found_fragments_, 0.20); 
  distance_.show_found_fragments(found_fragments_);
  known_submodels_["1ctf_a_38_1"] = 1.0;
  known_submodels_["1ctf_a_42_1"] = 1.0;
  known_submodels_["3chy_a_1_1"]  = 1.0;
  known_submodels_["3chy_a_44_1"] = 1.0;
}

void TestSubStructure::test_query_to_ref_submodels() { 
  using namespace std;
  map<string, float> carmsd_to_refs;
  distance_.compute_distance_to_references(submodels_, carmsd_to_refs); 
  map<string, float>::iterator mit = carmsd_to_refs.begin();
  while(mit != carmsd_to_refs.end()) {
    string submodel_id = mit->first;
    transform(submodel_id.begin(), submodel_id.end(), submodel_id.begin(), ::tolower);
    cout<<submodel_id<<endl;
    if(submodel_id== "1ctf_a_32_1")
      CPPUNIT_ASSERT(mit->second==3.50);
    else if(submodel_id == "1ctf_a_40_1")
      CPPUNIT_ASSERT(mit->second==3.29);
    else if(submodel_id == "3chy_a_15_1")
      CPPUNIT_ASSERT(mit->second==3.73);
    mit++;
  }
}

void TestSubStructure::test_found_fragments() {
  map<string, map<string,float> >::iterator vit = found_fragments_.begin();
  cout<<"This is with threshold 0.20"<<endl;
  while(vit != found_fragments_.end()) {
    map<string, float> submodel_ids   = vit->second;  
    map<string, float>::iterator mit1 = submodel_ids.begin();
    map<string, float>::iterator mit2 = known_submodels_.begin();
    while(mit1 != submodel_ids.end()) {
      string submodel_id = mit1->first;
      transform(submodel_id.begin(), submodel_id.end(), submodel_id.begin(), ::tolower);
      CPPUNIT_ASSERT(submodel_id == mit2->first);
      mit1++; mit2++;
    }
    vit++;
  }
}

void TestSubStructure::test_q_score() {
  CPPUNIT_ASSERT(found_fragments_.size()>0);
  CPPUNIT_ASSERT(submodels_.size()>0);
  map<string, SubModel> selected_submodels;
  Distance distance;
  distance.get_selected_submodels(found_fragments_, "./coordinates.out", selected_submodels);
  // cout<<"Submodel size: "<<selected_submodels.size()<<endl;
  // cout<<"Number of selected submodels: "<<selected_submodels.size()<<endl;
  distance.compute_and_dump_cosine(submodels_, selected_submodels, found_fragments_, "./cosine.out");
}

void TestSubStructure::compute_carmsd_query_to_database() {
  Distance distance;
  distance.set_reference_submodels("./3chy_a_2_1_query.pdb", 9);
  distance.compute_and_dump_carmsd("./coordinates.out");
}
