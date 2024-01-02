#ifndef _CLUS_STR_H_
#define _CLUS_STR_H_

#include "calccorr.h"
#include "stringo.h"
#include "lattype.h"
//#include <fstream>
#include "mpiinterf.h"

extern Real complexity_exp;

typedef Array<rVector3d> Cluster;
typedef Array<Cluster> ArrayCluster;
typedef Array<int> Arrayint;

inline Real get_length_quick(const Array<rVector3d> &c) {
  if (c.get_size()==0) return 0;
  return norm(c(c.get_size()-1)-c(0));
}

Real get_length_slow(const Array<rVector3d> &c);

class ClusterBank {
  LinkedList<Cluster> cluster_list;
  LinkedListIterator<Cluster> current_cluster;
  int current_index;
  Real previous_length;
  AtomMultipletIterator multiplet;
  const SpaceGroup &equivalent_by_symmetry;
 public:
  ClusterBank(const rMatrix3d &cell,
              const Array<rVector3d> &atom_pos,
              int ntuple,
              const SpaceGroup &_equivalent_by_symmetry);
  ClusterBank(const ClusterBank &clusterbank,int ntuple);
  void reset(void);
  void operator++(int);
  int get_current_index(void) {
    return current_index;
  }
  Real get_previous_length(void) {
    return previous_length;
  }
  Real get_current_length(void) {
    return get_length_quick(*current_cluster);
  }
  operator Cluster& (void) {
    return *current_cluster;
  }
  LinkedList<Cluster>& get_cluster_list(void) {
    return cluster_list;
  }
};

class MultiClusterBank {
  Structure lat;
  ClusterBank cluster_bank;
  LinkedList<MultiCluster> cluster_list;
  LinkedListIterator<MultiCluster> current_cluster;
  int current_index;
  Real previous_length;
  const SpaceGroup &equivalent_by_symmetry;

  void make_new_clusters(void);
 public:
  MultiClusterBank(const Structure &_lat,
              int ntuple,
              const SpaceGroup &_equivalent_by_symmetry);
  MultiClusterBank(const MultiClusterBank &clusterbank,int ntuple);
  void reset(void);
  void operator++(int);
  int get_current_index(void) {
    return current_index;
  }
  Real get_previous_length(void) {
    return previous_length;
  }
  Real get_current_length(void) {
    return get_length_quick(current_cluster->clus);
  }
  operator MultiCluster& (void) {
    return *current_cluster;
  }
  MultiCluster* operator -> () {
    return current_cluster;
  }
  LinkedList<MultiCluster>& get_cluster_list(void) {
    return cluster_list;
  }
};


// The best algorithm is invoked by not specifying any -D... in the makefile;
// If you want to use MPI, use -DSLOWENUMALGO ;
// For debugging, you can use the slowest algorithm with -DSLOWENUMALGO -DOLD_STR_ALGO ;

#ifndef SLOWENUMALGO
// This is the best algorithm;
template<class T>
class StructureBank {
  int current_index;
  const SpaceGroup &equivalent_by_symmetry;
  Structure basic_structure;
  int current_volume;
  LinkedList<T> structure_list;
  LinkedListIterator<T> current_structure;
  int do2D;

  MultiDimIterator<Arrayint> curnew_config;
  Array<rMatrix3d> supercells;
  int curnew_supercell;
  LinkedListIterator<T> beginning_of_cell;
  Structure blank_superstructure;
  SpaceGroup subgroup;
  int useradded;

  void init(void) {
    current_volume=0;
    current_index=0;
    curnew_supercell=0;
    useradded=0;

    // special case to put pure structures first;
    LinkedListIterator<T> insert_at(structure_list);
    blank_superstructure.cell=find_symmetric_cell(basic_structure.cell);
    find_all_atom_in_supercell(&blank_superstructure.atom_pos,
                               &blank_superstructure.atom_type,basic_structure.atom_pos,
                               basic_structure.atom_type,
                               basic_structure.cell, blank_superstructure.cell);
    Array<int> savetype=blank_superstructure.atom_type;
    for (int t=0; t<max(savetype); t++) {
      for (int s=0; s<blank_superstructure.atom_type.get_size(); s++) {
	blank_superstructure.atom_type(s)=MIN(t,savetype(s)-1);
      }
      LinkedListIterator<T> i(structure_list);
      for ( ; i ; i++) {
	if (equivalent_by_symmetry(*i,blank_superstructure)) break;
      }
      if (!i) {
	structure_list.add(new T(blank_superstructure),insert_at);
	insert_at++;
      }
    }
    reset();
  }
  int find_new_structure(void);
 public:
  StructureBank(const Structure &_basic_structure, const SpaceGroup &_equivalent_by_symmetry, int _do2D=0):
    structure_list(), current_structure(),
    equivalent_by_symmetry(_equivalent_by_symmetry),
    basic_structure(_basic_structure), do2D(_do2D), supercells(), curnew_config(), beginning_of_cell(), subgroup() {init();}
  StructureBank(const Structure &_basic_structure,
		const Array<Arrayint> &site_type_list,
		const SpaceGroup &_equivalent_by_symmetry, int _do2D=0):
    structure_list(), current_structure(),
    equivalent_by_symmetry(_equivalent_by_symmetry),
    basic_structure(_basic_structure), do2D(_do2D), supercells(), curnew_config(), beginning_of_cell(), subgroup() {
      for (int i=0; i<basic_structure.atom_type.get_size(); i++) {
        basic_structure.atom_type(i)=site_type_list(basic_structure.atom_type(i)).get_size();
      }
      init();
  }
  void set_2D_mode(int _do2D) {
    do2D=_do2D;
  }
  void reset(void) {
    current_index=0;
    current_structure.init(structure_list);
  }
  T& get_current_structure(void) {
    find_new_structure();
    return *current_structure;
  }
  operator T& (void) {
    find_new_structure();
    return *current_structure;
  }
  void operator++(int) {
    //    if (find_new_structure()) curnew_config++;
    find_new_structure();
    current_structure++;
    current_index++;
  }
  LinkedList<T>& get_structure_list(void) {
    return structure_list;
  }
  int get_current_index(void) {
    return current_index;
  }
  int get_max_volume(void) {
    return current_volume;
  }
  /*
  int end_of_this_volume(void) { // caution: this routine only works if add_structure is not called.
    if (supercells.get_size()==0) {return 0;}
    if (curnew_config) {return 0;} else {return 1;}
  }
  */
  int find_new_structures(int new_volume) {
    while (1) {
      if (current_volume>new_volume) break;
      (*this)++;
    }
    return 0;
  } 
  int add_structure(const Structure &str, T **pp_str=NULL) {
    LinkedListIterator<T> i(structure_list);
    for ( ; i ; i++) {
      if (equivalent_by_symmetry(*i,str)) break;
    }
    if (!i) {
      LinkedListIterator<T> insert_at(structure_list);
      while (insert_at && insert_at->atom_pos.get_size() < str.atom_pos.get_size()+1) insert_at++;
      T *p_str=new T(str);
      structure_list.add(p_str,insert_at);
      useradded=1;
      if (pp_str) *pp_str=p_str;
      return 1;
    }
    else {
      if (pp_str) *pp_str=i;
      return 0;
    }
  }
};

#define ENUMFLIPTRICK

template<class T>
int StructureBank<T>::find_new_structure(void) {
  while (1) {
    if (!curnew_config) {
      curnew_supercell++;
      if (curnew_supercell>=supercells.get_size()) {
	current_volume=current_volume+1;
	//cerr << "volume= " << current_volume << endl;
	Array<rMatrix3d> pointgroup;
	pointgroup_from_spacegroup(&pointgroup,equivalent_by_symmetry.point_op);
	if (do2D) {
	  find_supercells_2D(&supercells, current_volume, current_volume, basic_structure.cell, pointgroup);
	}
	else {
	  find_supercells(&supercells, current_volume, current_volume, basic_structure.cell, pointgroup);
	}
	curnew_supercell=0;
	int target_nbatom=basic_structure.atom_pos.get_size()*current_volume;
	beginning_of_cell.init(structure_list);
	//cerr << "begincell [\n";
	while (1) {
	  if (!beginning_of_cell) break;
	  if (beginning_of_cell->atom_pos.get_size()>=target_nbatom) break;
	  beginning_of_cell++;
	}
	//cerr << "begincell ]\n";
      }
      if (!useradded && current_volume>1) {beginning_of_cell=current_structure;}
      //cerr << "findsymcell [\n";
      supercells(curnew_supercell)=find_symmetric_cell(supercells(curnew_supercell));
      //cerr << "findsymcell ]\n";
      blank_superstructure.cell=supercells(curnew_supercell);
      find_all_atom_in_supercell(&blank_superstructure.atom_pos,
				 &blank_superstructure.atom_type,basic_structure.atom_pos,
				 basic_structure.atom_type,
				 basic_structure.cell, blank_superstructure.cell);
      curnew_config.init(blank_superstructure.atom_type);
      subfactorgroup_from_supercell_and_spacegroup(&subgroup.point_op,&subgroup.trans, blank_superstructure.cell,equivalent_by_symmetry.point_op,equivalent_by_symmetry.trans);
      //      cerr << "point = " << subgroup.point_op << endl;
      //      cerr << "trans = " << subgroup.trans << endl;

    }

    int need_to_gen=0;
    if (!current_structure) {
      need_to_gen=1;
    }
    else if (current_structure->atom_pos.get_size()/basic_structure.atom_pos.get_size() > MAX(1,current_volume)) {
      need_to_gen=1;
    }
    if (!need_to_gen) return 0;

    for (int s=0; s<((Arrayint)curnew_config).get_size(); s++) {
      blank_superstructure.atom_type(s)=((Arrayint)curnew_config)(s);
    }
    curnew_config++;
#ifdef ENUMFLIPTRICK
    if (!contains_pure_translations(blank_superstructure,basic_structure.cell)) {
      if (!equiv_to_lexico_successor(blank_superstructure,basic_structure.cell,subgroup.point_op,subgroup.trans)) {
	if (!useradded && current_volume>1) {
	  structure_list.add(new T(blank_superstructure),current_structure);
	  break; // exit routine;
	}
	else {
	  LinkedListIterator<T> i=beginning_of_cell;
	  for ( ; i ; i++) {
	    if (equivalent_by_symmetry(*i,blank_superstructure)) break;
	  }
	  if (!i) {
	    structure_list.add(new T(blank_superstructure),current_structure);
	    break; // exit routine;
	  }
	}
      }
    }
#else
    if (!contains_pure_translations_or_lexico_successor(blank_superstructure,basic_structure.cell)) {
      LinkedListIterator<T> i=beginning_of_cell;
      for ( ; i ; i++) {
	if (equivalent_by_symmetry(*i,blank_superstructure)) break;
      }
      if (!i) {
	structure_list.add(new T(blank_superstructure),current_structure);
	break; // exit routine;
      }
    }
#endif
    //    curnew_config++;
  }
  return 0;
}
#else
template<class T>
class StructureBank {
  int current_index;
  const SpaceGroup &equivalent_by_symmetry;
  Structure basic_structure;
  int current_volume;
  LinkedList<T> structure_list;
  LinkedListIterator<T> current_structure;
  int do2D;
  void init(void) {
    current_volume=0;
    current_index=0;

    // special case to put pure structures first;
    LinkedListIterator<T> insert_at(structure_list);
    Structure blank_superstructure;
    blank_superstructure.cell=find_symmetric_cell(basic_structure.cell);
    find_all_atom_in_supercell(&blank_superstructure.atom_pos,
                               &blank_superstructure.atom_type,basic_structure.atom_pos,
                               basic_structure.atom_type,
                               basic_structure.cell, blank_superstructure.cell);
    Array<int> savetype=blank_superstructure.atom_type;
    for (int t=0; t<max(savetype); t++) {
      for (int s=0; s<blank_superstructure.atom_type.get_size(); s++) {
	blank_superstructure.atom_type(s)=MIN(t,savetype(s)-1);
      }
      LinkedListIterator<T> i(structure_list);
      for ( ; i ; i++) {
	if (equivalent_by_symmetry(*i,blank_superstructure)) break;
      }
      if (!i) {
	structure_list.add(new T(blank_superstructure),insert_at);
	insert_at++;
      }
    }
    // find_new_structures(current_volume+1);
    reset();
  }
 public:
  StructureBank(const Structure &_basic_structure, const SpaceGroup &_equivalent_by_symmetry, int _do2D=0):
    structure_list(), current_structure(),
    equivalent_by_symmetry(_equivalent_by_symmetry),
    basic_structure(_basic_structure), do2D(_do2D) {init();}
  StructureBank(const Structure &_basic_structure,
		const Array<Arrayint> &site_type_list,
		const SpaceGroup &_equivalent_by_symmetry, int _do2D=0):
    structure_list(), current_structure(),
    equivalent_by_symmetry(_equivalent_by_symmetry),
    basic_structure(_basic_structure), do2D(_do2D) {
      for (int i=0; i<basic_structure.atom_type.get_size(); i++) {
        basic_structure.atom_type(i)=site_type_list(basic_structure.atom_type(i)).get_size();
      }
      init();
  }
  void set_2D_mode(int _do2D) {
    do2D=_do2D;
  }
  void reset(void) {
    current_index=0;
    current_structure.init(structure_list);
  }
  void operator++(int) {
    current_structure++;
    current_index++;
    if (!current_structure) {
      while (!current_structure) find_new_structures(current_volume+1);
    }
    else if (current_structure->atom_pos.get_size()/basic_structure.atom_pos.get_size() > MAX(1,current_volume)) {
      find_new_structures(current_volume+1);
    }
  }
  operator T& (void) {
    return *current_structure;
  }
  T& get_current_structure(void) {
    return *current_structure;
  }
  LinkedList<T>& get_structure_list(void) {
    return structure_list;
  }
  int get_current_index(void) {
    return current_index;
  }
  int get_max_volume(void) {
    return current_volume;
  }
  int end_of_this_volume(void) {
    if (current_volume==0) {
      find_new_structures(current_volume+1);
    }
    LinkedListIterator<T> tmp_structure=current_structure;
    tmp_structure++;
    if (!tmp_structure) return 1;
    if (tmp_structure->atom_pos.get_size()/basic_structure.atom_pos.get_size() > current_volume) return 1;
    return 0;
  }
  void find_new_structures(int new_max_volume);
  int add_structure(const Structure &str, T **pp_str=NULL) {
    LinkedListIterator<T> i(structure_list);
    for ( ; i ; i++) {
      if (equivalent_by_symmetry(*i,str)) break;
    }
    if (!i) {
      LinkedListIterator<T> insert_at(structure_list);
      while (insert_at && insert_at->atom_pos.get_size() < str.atom_pos.get_size()+1) insert_at++;
      T *p_str=new T(str);
      structure_list.add(p_str,insert_at);
      if (pp_str) *pp_str=p_str;
      return 1;
    }
    else {
      if (pp_str) *pp_str=i;
      return 0;
    }
  }
};

//#define OLD_STR_ALGO
#ifdef OLD_STR_ALGO
template<class T>
void StructureBank<T>::find_new_structures(int new_max_volume) {
  LinkedListIterator<T> insert_at(structure_list);
  int nb_atom=(current_volume+1)*basic_structure.atom_pos.get_size();
  while (insert_at && insert_at->atom_pos.get_size() <= nb_atom) insert_at++;

  Array<rMatrix3d> pointgroup;
  pointgroup_from_spacegroup(&pointgroup,equivalent_by_symmetry.point_op);
  Array<rMatrix3d> supercell;
  if (do2D) {
    find_supercells_2D(&supercell, current_volume+1, new_max_volume, basic_structure.cell, pointgroup);
  }
  else {
    find_supercells(&supercell, current_volume+1, new_max_volume, basic_structure.cell, pointgroup);
  }
  current_volume=MAX(current_volume,new_max_volume);
  for (int c=0; c<supercell.get_size(); c++) {
    Structure blank_superstructure;
    blank_superstructure.cell=find_symmetric_cell(supercell(c));
    find_all_atom_in_supercell(&blank_superstructure.atom_pos,
			       &blank_superstructure.atom_type,basic_structure.atom_pos,
			       basic_structure.atom_type,
			       basic_structure.cell, blank_superstructure.cell);
    MultiDimIterator<Arrayint> config(blank_superstructure.atom_type);
    for (; config; config++) {
      for (int s=0; s<((Arrayint)config).get_size(); s++) {
        blank_superstructure.atom_type(s)=((Arrayint)config)(s);
      }
      LinkedListIterator<T> i(structure_list);
      for ( ; i ; i++) {
	if (equivalent_by_symmetry(*i,blank_superstructure)) break;
      }
      if (!i) {
	structure_list.add(new T(blank_superstructure),insert_at);
	insert_at++;
      }
    }
  }
/*
  {
    ofstream file("strdebug.out");
    LinkedListIterator<T> ii(structure_list);
    for (; ii; ii++) {
      file << ii->cell << endl;
      file << ii->atom_pos << endl;
      file << ii->atom_type << endl;
      file << "end" << endl << endl;
    }
  }
*/
}
#else
template<class T>
void StructureBank<T>::find_new_structures(int new_max_volume) {
  Array<rMatrix3d> pointgroup;
  pointgroup_from_spacegroup(&pointgroup,equivalent_by_symmetry.point_op);
  Array<rMatrix3d> supercell;
  if (do2D) {
    find_supercells_2D(&supercell, current_volume+1, new_max_volume, basic_structure.cell, pointgroup);
  }
  else {
    find_supercells(&supercell, current_volume+1, new_max_volume, basic_structure.cell, pointgroup);
  }
  current_volume=MAX(current_volume,new_max_volume);
  LinkedListIterator<T> insert_at(structure_list);
  int nb_atom=current_volume*basic_structure.atom_pos.get_size();
  while (insert_at && insert_at->atom_pos.get_size() < nb_atom) insert_at++;

  Array<LinkedList<T> > partial_str_list(supercell.get_size());
  {
      MPISynchronizer<LinkedList<T> > sync;
      for (int c=0; c<supercell.get_size(); c++) {
	  if (sync.is_my_job()) {
//cerr << "Process " << MyMPIobj.id << " is doing cell " << c << endl;
	      supercell(c)=find_symmetric_cell(supercell(c));
	      Structure blank_superstructure;
	      blank_superstructure.cell=supercell(c);
	      find_all_atom_in_supercell(&blank_superstructure.atom_pos,
					 &blank_superstructure.atom_type,basic_structure.atom_pos,
					 basic_structure.atom_type,
					 basic_structure.cell, blank_superstructure.cell);
	      MultiDimIterator<Arrayint> config(blank_superstructure.atom_type);
	      for (; config; config++) {
		  for (int s=0; s<((Arrayint)config).get_size(); s++) {
		      blank_superstructure.atom_type(s)=((Arrayint)config)(s);
		  }
		  if (!contains_pure_translations_or_lexico_successor(blank_superstructure,basic_structure.cell)) {
		      LinkedListIterator<T> i(partial_str_list(c));
		      for ( ; i ; i++) {
			  if (equivalent_by_symmetry(*i,blank_superstructure)) break;
		      }
		      if (!i) {
			  partial_str_list(c) << new T(blank_superstructure);
		      }
		  }
	      }
	  
	      LinkedListIterator<T> i(partial_str_list(c));
	      while (i) {
		LinkedListIterator<T> j(insert_at);
		int found=0;
		while (j) {
		  if (j->atom_pos.get_size()!=i->atom_pos.get_size()) break;
		  if (equivalent_by_symmetry(*i,*j)) {
		    found=1;
		    break;
		  }
		  j++;
		}
		if (found) {
		  delete partial_str_list(c).detach(i);
		}
		else {
		  i++;
		}
	      }

	  }
	  sync.sync(&(partial_str_list(c)));
      }
  }
  
  while (insert_at && insert_at->atom_pos.get_size() == nb_atom) insert_at++;
  for (int c=0; c<supercell.get_size(); c++) {
    LinkedListIterator<T> i(partial_str_list(c));
    for ( ; i; insert_at++) {
      structure_list.add(partial_str_list(c).detach(i),insert_at);
    }
  }
/*
  {
    ofstream file("strdebug.out");
    LinkedListIterator<T> ii(structure_list);
    for (; ii; ii++) {
      file << ii->cell << endl;
      file << ii->atom_pos << endl;
      file << ii->atom_type << endl;
      file << "end" << endl << endl;
    }
  }
*/
}
#endif

#endif


#endif
