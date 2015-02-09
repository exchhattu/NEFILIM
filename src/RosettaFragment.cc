
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<functional>
#include<numeric>

#include <core/scoring/rms_util.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include<RosettaFragment.hh>
#include<rmsd.hh>

using namespace std;
using namespace core::scoring;
using namespace core;


RosettaFragment::RosettaFragment(const int rsd_pos, const string line, 
                                  const int window_size, const bool is_in_top25) { 
  if(!line.empty()) { 
    rsd_pos_ = rsd_pos;
    string aminoacid, secstruct, pdb_id, chain_id;
    float phi, psi, omega;
    stringstream streams(line);
    streams>>pdb_id>>chain_id>>start_pos_>>aminoacid>>secstruct>>phi>>psi>>omega;
    end_pos_ = start_pos_ + (window_size - 1);
    transform(pdb_id.begin(), pdb_id.end(), pdb_id.begin(),::toupper);
    if(chain_id=="-" || chain_id=="_") 
      chain_id_ = "A";
    else
      chain_id_ = chain_id;
    pdb_id_ = pdb_id;
    subseq_ = aminoacid; subsec_ = secstruct;
    phis_.push_back(phi); psis_.push_back(psi); omegas_.push_back(omega);
    window_size_ = window_size;
    is_in_top25_ = is_in_top25;
  }
} 

void RosettaFragment::parse_line(const string line) {
  if(line.empty()) return; 
  int cur_rsd_pos;
  float phi, psi, omega;
  string aminoacid, secstruct, pdb_id, chain_id;
  stringstream streams(line);
  streams>>pdb_id>>chain_id>>cur_rsd_pos>>aminoacid>>secstruct>>phi>>psi>>omega;

  subseq_ += aminoacid; subsec_ += secstruct;
  phis_.push_back(phi); psis_.push_back(psi); omegas_.push_back(omega);
  if(phis_.size() == (size_t)window_size_ && 
     psis_.size() == (size_t)window_size_ && 
     omegas_.size() == (size_t)window_size_) {
    Pose pose;
    chemical::make_pose_from_sequence(pose, subseq_,
             *(chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID )));
    assert(pose.total_residue() == (size_t)window_size_);
    for(size_t pos = 1; pos <= pose.total_residue(); pos++) {
      if (!pose.residue(pos).is_protein() ) continue;
      pose.set_phi(pos, phis_[pos-1]); 
      pose.set_psi(pos, psis_[pos-1]); 
      pose.set_omega(pos, omegas_[pos-1]); 
    }
    pose_ = pose;
  }
}

void RosettaFragment::print2stdoutput() {
  cout<<"pdb_id: "<<pdb_id_<<" "
      <<chain_id_ <<" "<<rsd_pos_<<" "
      <<start_pos_<<" "<<end_pos_<<" "
      <<window_size_<<" "<<subseq_<<" "
      <<subsec_<<" "<<carmsd_<<"  "
      <<is_in_top25_<<endl;
  // int total_rsd = pose_.total_residue();
  // for(int i = 1; i <= total_rsd; i++) {
  //   cout<<"pdb id "<<i<<" "<<phis_[i-1]<<" "<<psis_[i-1]<<" "<<omegas_[i-1]<<endl;
  //   cout<<"pdb id "<<i<<" "<<pose_.phi(i)<<" "<<pose_.psi(i)<<" "<<pose_.omega(i)<<endl;
  // } 
}

void RosettaFragment::carmsd(Pose segmented_native_pose) {
  carmsd_ = CA_rmsd(segmented_native_pose, pose_, 1, window_size_); 
}

void RosettaFragment::collect_carmsd(map<int, list<float> > &rsd_carmsds) {
  list<float> carmsds = rsd_carmsds[rsd_pos_]; 
  carmsds.push_back(carmsd_);
  rsd_carmsds[rsd_pos_] = carmsds;
}

void RosettaFragment::collect_top25_carmsd(map<int, list<float> > &rsd_carmsds) {
  if(!is_in_top25_) return;
  list<float> carmsds = rsd_carmsds[rsd_pos_]; 
  carmsds.push_back(carmsd_);
  rsd_carmsds[rsd_pos_] = carmsds;
} 

Residues::Residues() {
}

Residues::Residues(int residue_pos, RosettaFragment rosetta_fragment) {
  list<RosettaFragment> fragments;
  fragments.clear();
  fragments.push_back(rosetta_fragment);
  residues_[residue_pos] = fragments;
}

void Residues::add(int residue_pos, RosettaFragment rosetta_fragment) {
  list<RosettaFragment> fragments = residues_[residue_pos];
  if(!fragments.empty()) 
    fragments.push_back(rosetta_fragment);
  else 
    fragments.push_back(rosetta_fragment);
  residues_[residue_pos] = fragments;
} 

map<string, Pose> Residues::get_poses(int rsd_pos) {
  map<string, Pose> poses;
  map<int, list<RosettaFragment> > residues; 
  map<int, list<RosettaFragment> >::iterator milt = residues.begin(); 
  while(milt != residues.end()) {
    if(milt->first == rsd_pos) {
      list<RosettaFragment> rosetta_fragments = milt->second;
      list<RosettaFragment>::iterator loit = rosetta_fragments.begin();
      while(loit != rosetta_fragments.end()) {
        Pose pose = loit->get_pose();
        string pdb_id = loit->get_pdb_id_with_chain();
        poses[pdb_id] = pose;
      }
    }
    milt++;
  }  
}

void Residues::do_clustering(float cluster_radius) {
  // map<int, list<RosettaFragment> > residues; 
  // map<int, list<RosettaFragment> >::iterator milt = residues.begin(); 
  // while(milt != residues.end()) {
  //   int rsd_pos = milt->first;
  //   map<string, Pose> fragment_poses = get_poses(rsd_pos); 

  //   int mode_  = 1; 
  //   int smode_ = 2;
  //   int number_of_fragments_  = 0;
  //   int fragment_per_cluster_ = 0; 

  //   //do clustering...
  //   // FragmentCluster clustering_;
  //   // clustering_.set_total_residue(total_residue_);
  //   // clustering_.set_poses(fragment_poses);
  //   // clustering_.do_clustering_using_durandal(cluster_radius); 
  //   // map<string, Pose> selected_templates = clustering_.get_selected_templates(mode_, 
  //   //                                                                           smode_,
  //   //                                                                           number_of_fragment_,
  //   //                                                                           fragment_per_cluster_); 
  //   // write_clustered_templates(f_output_fragment, clustering_, selected_templates);
  //   // write_template_fragments_pdb(f_output_template, clustering_);
  //   // write_cluster_info(f_output, clustering_);
  //   // write_selected_templates_native_RMSD(f_output, clustering_);

  //   //delete after valgrid check 
  //   selected_fragment_poses_.clear();
  // }
  //   mlit++;
  // }  
}

// void Residues::write_fragments_selected_using_clustering(ofstream &f_output_fragment,
//                                                          FragmentCluster clustering_, 
//                                                          map<string, Pose> selected_templates) {
//   try {
//     if(selected_templates.empty()) 
//       throw  " templates are not found for residue. " + rsd_pos_;
//     map<string, Pose> selected_extended_templates;
//     // selected_extended_templates = get_extended_selected_templates(selected_templates); 
//     if(selected_templates.empty()) 
//       throw  " templates are not found for residue. " + rsd_pos_;
//     clustering_.write_fragments(f_output_fragment, selected_extended_templates, rsd_pos_); 
//     //delete after valgrid check 
//     // selected_fragment_poses_.clear();
//   }
//   catch(const char* e) {
//     std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<e<<std::endl;
//     exit(1);
//   }
//   catch(...) {
//     std::cerr<<"fatal: "<<__FILE__<<" "<<__LINE__<<" unknown error."<<std::endl;
//     exit(1);
//   }
// }


