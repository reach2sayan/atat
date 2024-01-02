#ifndef __GETVALUE_H__
#define __GETVALUE_H__

#include <iostream>
#include "stringo.h"

enum VarType {TITLEVAL,INTVAL,REALVAL,BOOLVAL,STRINGVAL,CHOICEVAL,ARRAYRVAL,LISTVAL};

struct AskStruct {     // contains the caracteristics of a command-line option;
	const char *shortname; // command line option name;
	const char *longname;  // description string;
	VarType vartype; // type of the variable to initialize (see VarType);
	void *outvar;    // pointer to the variable to initialize;
};

// use the stream s to initialize the variables pointed to in the array questions of lenght nb;
int get_values(istream &s, int nb, AskStruct *label);
// use argc and argv to initialize the variables pointed to in the array questions of lenght nb;
int get_values(int argc, char* argv[], int nb, AskStruct *label);

// Displays the description strings;
void display_help(int nb, AskStruct *label);

#define CMDLINEPARAM(VAR,TYP) {"-"#VAR,#VAR,TYP,&VAR}

void chdir_robust(const char *dir);

int get_string(AutoString *ps, istream &file, const char *delim=" \t\n");
int skip_delim(istream &file, const char *delim=" \t\n");
int get_row_numbers(Array<Real> *pa, istream &file);
void read_table(Array<Array<Real> > *pa, istream &file, int keepempty);
void get_atat_root(AutoString *patatroot);

#endif
