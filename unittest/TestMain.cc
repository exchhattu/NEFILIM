#include<iostream>
#include<init_rosetta.hh>

#include<string>

#include<cppunit/TestFixture.h>

// #include<TestSegment.hh>
// #include<TestFragment.hh>
#include<TestSubStructure.hh>

using namespace std;


Rosetta rosetta_option_wrapper(int argc, char* argv[]) {
  if(argc < 2) return Rosetta(0);
  char** argvs = new char*[2];
  const int new_argc = 2;
  
  for(int i = 0; i < argc ; i++) {
    if(strcmp(argv[i],"-r") == 0) {
      argvs[0] = argv[0];
      argvs[1] = argv[i+1];
      return Rosetta(new_argc, argvs);
    }
  }
  return Rosetta(0);
}


int main(int argc, char* argv[]) {
  cout<<"Unit testing ..."<<endl;
  Rosetta rosetta = rosetta_option_wrapper(argc, argv);
  if(rosetta.get_init() == 0) 
    cout<<"Fatal: rosetta database has to be set. "<<endl;

  TestSubStructure testSubStructure;
  testSubStructure.setup();
  testSubStructure.test_CARMSD_Equality();

  testSubStructure.setup_search(); 
  testSubStructure.test_query_to_ref_submodels(); 
  testSubStructure.test_found_fragments(); 

  testSubStructure.test_q_score(); 
  // testSubStructure.compute_carmsd_query_to_database(); 

  return 1;
}

