#include <strstream>
#include <fstream>
#include "getvalue.h"
#include "linklist.h"

const char *equal_delim=":= \t\n";
const char *blanks_delim=" \t\n";

LinkedList<AutoString> getvalue_string_buffer;

int get_values(istream &s, int nb, AskStruct *label) {
  char c;
  while (1) {
    while (1) {
      c=s.get();
      if (!strchr(blanks_delim,c) || s.eof()) break;
    }
    if (s.eof()) return 1;
    s.putback(c);
    ostrstream cur_label;
    while (1) {
      c=s.get();
      if (strchr(equal_delim,c) || s.eof()) break;
      cur_label << c;
    }
    cur_label << '\0';
    while (1) {
      c=s.get();
      if (!strchr(equal_delim,c) || s.eof()) break;
    }
    if (!s.eof()) {s.putback(c);}
    int i;
    for (i=0; i<nb; i++) {
      if (strcmp(label[i].shortname,cur_label.str())==0) {
	ostrstream value_str;
	if (label[i].vartype!=BOOLVAL) {
	  while (1) {
	    c=s.get();
	    if (strchr(blanks_delim,c) || s.eof()) break;
	    value_str << c;
	  }
	}
	value_str << '\0';
	istrstream value(value_str.str());
	switch (label[i].vartype) {
        case INTVAL:
          value >> (*(int *)(label[i].outvar));
	  break;
        case REALVAL:
          value >> (*(Real *)(label[i].outvar));
	  break;
        case BOOLVAL:
          (*(int *)(label[i].outvar))=1;
	  break;
        case STRINGVAL:
	  {
	    AutoString *pstr=new AutoString(value_str.str());
	    getvalue_string_buffer << pstr;
	    *((const char **)(label[i].outvar))=(const char *)(*pstr);
	  }
	  break;
        case CHOICEVAL:
          value >> (*(char *)(label[i].outvar));
	  break;
	case ARRAYRVAL:
	  {
	    LinkedList<Real> list;
	    while (1) {
	      skip_delim(value,",;:");
	      Real tmp=MAXFLOAT;
	      value >> tmp;
	      if (tmp==MAXFLOAT) break;
	      list << new Real(tmp);
	    }
	    LinkedList_to_Array((Array<Real> *)(label[i].outvar),list);
	  }
	  break;
	}
	break;
      }
    }
    if (i==nb) return 0;
  }
}

int get_values(int argc, char* argv[], int nb, AskStruct *label) {
  if (argc==1) return 0;
  ostrstream os;
  for (int i=1; i<argc; i++) {
    os << argv[i] << ' ';
  }
  os << '\0';
  istrstream is(os.str());
  return get_values(is, nb, label);
}

int trivial_grep(const char *line, const char *pat) {
  if (strlen(line)<strlen(pat)) {return 0;}
  for (int i=0; i<strlen(line)-strlen(pat); i++) {
    int j;
    for (j=0; j<strlen(pat); j++) {
      if (line[i+j]!=pat[j]) break;
    }
    if (j==strlen(pat)) return 1;
  }
  return 0;
}

void display_help(int nb, AskStruct *label) {
  int maxw=0;
  for (int i=0; i<nb; i++) {
    maxw=MAX(maxw,strlen(label[i].shortname));
  }
  for (int i=0; i<nb; i++) {
    if (label[i].vartype==TITLEVAL) {
      cerr << label[i].longname << endl;
    }
    else {
      Real r;
      cerr.width(maxw);
      cerr << label[i].shortname;
      int dodef=!trivial_grep(label[i].longname,"efault");
      switch (label[i].vartype){
      case INTVAL:
	cerr << "=[int]      " << label[i].longname;
	if (dodef) {
	  cerr <<" (Default: " << *(int *)label[i].outvar << ")";
	}
	cerr << endl;
	break;
      case REALVAL:
	r=*(Real *)label[i].outvar;
	cerr << "=[real]     " << label[i].longname;
	if (dodef) {
	  cerr << " (Default: ";
	  if (r==MAXFLOAT) {
	    cerr << "Off";
	  } else {
	    cerr << r;
	  }
	  cerr << ")";
	}
	cerr << endl;
	break;
      case BOOLVAL:
	cerr << "            " << label[i].longname << " (Default: " << (*(int *)label[i].outvar!=0 ? "On" : "Off") << ")" << endl;
	break;
      case STRINGVAL:
	cerr << "=[string]   " << label[i].longname;
	if (dodef) {
	  cerr << " (Default: " << *(char **)label[i].outvar << ")";
	}
	cerr<< endl;
	break;
      case CHOICEVAL: cerr << "=[choice]   " << label[i].longname;
	if (dodef) {
	  cerr << " (Default: " << *(char *)label[i].outvar << ")";
	}
	cerr << endl;
	break;
      case ARRAYRVAL:
	cerr << "=[real],... " << label[i].longname;
	if (dodef) {
	  cerr << " (Default: ";
	  for (int j=0; j<((Array<Real> *)(label[i].outvar))->get_size(); j++) {
	    if (j>0) {cerr << ",";}
	    cerr << (*(Array<Real> *)(label[i].outvar))(j);
	  }
	  cerr << ")";
	}
	cerr << endl;
      }
    }
  }
}

void chdir_robust(const char *dir) {
  if (chdir(dir)!=0) {
    cerr << "Cannot cd into " << dir << endl;
    ERRORQUIT("Aborting");
  }
}

int get_string(AutoString *ps, istream &file, const char *delim) {
  AutoString accum;
  char c;
  while (1) {
     file.get(c);
     if (strchr(delim,c)!=NULL || file.eof()) break;
     accum+=c;
  }
  *ps=accum;
  if (file.eof()) return 0;
  file.putback(c);
  return 1;
}

int skip_delim(istream &file, const char *delim) {
  char c;
  do {
     file.get(c);
  } while (strchr(delim,c)!=NULL && !file.eof());
  if (file.eof()) return 0;
  file.putback(c);
  return 1;
}

int get_row_numbers(Array<Real> *pa, istream &file) {
  if (file.eof()) {return 0;}
  char buf[MAX_LINE_LEN];
  buf[0]=0;
  file.get(buf,MAX_LINE_LEN-1);
  file.get();
  istrstream line(buf);
  LinkedList<Real> list;
  while (!line.eof()) {
    Real r=MAXFLOAT;
    line >> r;
    if (r==MAXFLOAT) break;
    list << new Real(r);
  }
  LinkedList_to_Array(pa,list);
  return 1;
}

void read_table(Array<Array<Real> > *pa, istream &file, int keepempty) {
  LinkedList<Array<Real> > alist;
  while (1) {
    char buf[MAX_LINE_LEN];
    buf[0]=0;
    file.clear();
    file.get(buf,MAX_LINE_LEN-1);
    if (file.eof()) break;
    istrstream line(buf);
    {
      LinkedList<Real> rlist;
      while (1) {
	Real r=MAXFLOAT;
	line >> r;
	if (r==MAXFLOAT) break;
	rlist << new Real(r);
      }
      if (rlist.get_size()>0 || keepempty) {
	Array<Real> *prs=new Array<Real>();
	LinkedList_to_Array(prs,rlist);
	alist << prs;
      }
    }
    file.clear();
    file.get();
    if (file.eof()) break;
  }
  LinkedList_to_Array(pa,alist);
}

void get_atat_root(AutoString *patatroot) {
  AutoString configfilename(getenv("HOME"));
  configfilename+="/.atat.rc";
  ifstream configfile(configfilename);
  if (!configfile) {
    ERRORQUIT("$HOME/.atat.rc was not found.");
  }
  while (configfile.get()!='=') {};
  skip_delim(configfile," \t");
  get_string(patatroot,configfile);
}
