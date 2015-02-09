#include<TestFragment.hh>
#include<FragmentGeneratorMain.hh>


void TestBruteFragment::setup(const string modelpath, const string refpath) {
  FileIO fileio;
  cout<<refpath<<endl;
  cout<<modelpath<<endl;
  Pose pose;
  cout<<"before"<<endl;
  core::io::pdb::centroid_pose_from_pdb(pose, "/home/rojan/src/ResiduesProfiles/ResiduesProfile/unittest/testdata/S_00001846.pdb");
  cout<<"why"<<endl;
  fileio.read_pose_no_stout(refpath, ref_model_id, ref_pose);

  // fileio.parse_file_path(modelpath, model_model_id);
  // cout<<"helllo"<<endl;
  
  carmsds_rsd_one.push_back(4.229); // for index 0 
  carmsds_rsd_one.push_back(1.061); // for index 12 
  carmsds_rsd_one.push_back(2.326); // for index 18 
  carmsds_rsd_one.push_back(3.794); // for index 58 

  carmsds_rsd_one.push_back(3.922); // for index 59 
}

void TestBruteFragment::testFileOpen() {
  CPPUNIT_ASSERT(!model_pose.empty());
  CPPUNIT_ASSERT(!ref_pose.empty());
  CPPUNIT_ASSERT(model_model_id=="S_00001846");
  CPPUNIT_ASSERT(ref_model_id=="1CTF_A");
}

void TestBruteFragment::testLocalRMSDCal() {
  // size_t rsd_pos = 1; 
  // const size_t frag_size = 9;
  // vector<float> ca_rmsds = compute_fragment_carmsd(ref_pose, model_pose, rsd_pos, frag_size); 
  // assert(ca_rmsds.size() == 60);
  // CPPUNIT_ASSERT(ca_rmsds[0] == carmsds_rsd_one[0]);
  // CPPUNIT_ASSERT(ca_rmsds[12] == carmsds_rsd_one[1]);
  // CPPUNIT_ASSERT(ca_rmsds[18] == carmsds_rsd_one[2]);
  // CPPUNIT_ASSERT(ca_rmsds[58] == carmsds_rsd_one[3]);

  // rsd_pos = 24; ca_rmsds.clear();
  // ca_rmsds = compute_fragment_carmsd(ref_pose, model_pose, rsd_pos, frag_size); 
  // CPPUNIT_ASSERT(ca_rmsds[59] == carmsds_rsd_one[5]);
}


