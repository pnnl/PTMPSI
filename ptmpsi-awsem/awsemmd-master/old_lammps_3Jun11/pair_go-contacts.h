/* ----------------------------------------------------------------------
Copyright (2010) Aram Davtyan and Garegin Papoian

Papoian's Group, University of Maryland at Collage Park
http://papoian.chem.umd.edu/

Last Update: 12/01/2010
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(gocontacts,PairGoContacts)

#else

#ifndef LMP_PAIR_GO_CONTACTS_H
#define LMP_PAIR_GO_CONTACTS_H

#include "pair.h"
#include "random_park.h"

namespace LAMMPS_NS {

class PairGoContacts : public Pair {
  friend class Pair;

 public:
  PairGoContacts(class LAMMPS *);
  ~PairGoContacts();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  
 private:
  char *parfile;
  double **cut, cut_global;
  double epsilon, epsilon2;
  int nres, ncont;
  int **isNative;
  int dev_type;
  int seed;
  double sigma0; // Non-native sigma distance
  double **sigma, **sigma_sq; // Equal either native distance seperation or sigma0
  double rnd;
  int lj_contacts_flag, contacts_flag, contacts_dev_flag;
  
  double dev, devA, devB, devC;
  double sdivf, tcorr, dev0;
  double w, xi;
  
  RanPark *rand;
  
  enum ContactsDevType{DT_NONE=0, DT_CORR, DT_SIN, DT_CONST};
  
  inline void print_log(char *line);
  void compute_contact_deviation();

  void allocate();
  
  FILE *fout;
};

}

#endif
#endif
