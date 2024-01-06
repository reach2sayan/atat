char *helpstring =
    ""
    "->What does this program do?\n"
    "\n"
    "1) It reads the lattice file (specified by the -l option).\n"
    "\n"
    "2) It determines the space group of this lattice and\n"
    "   writes it to the sym.out file.\n"
    "\n"
    "3) It finds all symmetrically distinct clusters that satisfy the\n"
    "   conditions specified by the options -2 through -6.\n"
    "   For instance, if -2=2.1 -3=1.1 is specified,\n"
    "   only pairs shorter than 2.1 units and triplets containing\n"
    "   no pairs longer than 1.1 will be selected.\n"
    "\n"
    "4) It writes all clusters found to clusters.out.\n"
    "   If the -c option is specified, clusters are read from clusters.out "
    "instead.\n"
    "\n"
    "5) It reads the structure file (specified by the -s option).\n"
    "\n"
    "6) It determines, for that structure, the correlations associated with "
    "all\n"
    "   the clusters chosen earlier.\n"
    "   This information is then output on one line, in the same order as in "
    "the\n"
    "   clusters.out file. See below for conventions used to calculate "
    "correlations.\n"
    "\n"
    "7) It writes the files corrdump.log containting the list of all "
    "adjustements\n"
    "   needed to map the (possibly relaxed) structure onto the ideal "
    "lattice.\n"
    "\n"
    "->File formats\n"
    "\n"
    "Lattice and structure files\n"
    "\n"
    "Both the lattice and the structure files have a similar structure.\n"
    "First, the coordinate system a,b,c is specified, either as\n"
    "[a] [b] [c] [alpha] [beta] [gamma]\n"
    "or as:\n"
    "[ax] [ay] [az]\n"
    "[bx] [by] [bz]\n"
    "[cx] [cy] [cz]\n"
    "Then the lattice vectors u,v,w are listed, expressed in the coordinate "
    "system\n"
    "just defined:\n"
    "[ua] [ub] [uc]\n"
    "[va] [vb] [vc]\n"
    "[wa] [wb] [wc]\n"
    "Finally, atom positions and types are given, expressed in the same "
    "coordinate system\n"
    "as the lattice vectors:\n"
    "[atom1a] [atom1b] [atom1c] [atom1type]\n"
    "[atom2a] [atom2b] [atom2c] [atom1type]\n"
    "etc.\n"
    "\n"
    "In the lattice file:\n"
    "-The atom type is a comma-separated list of the atomic\n"
    " symbols of the atoms that can sit the lattice site.\n"
    "-The first symbol listed is assigned a spin of -1.\n"
    "-When only one symbol is listed, this site is ignored for the purpose\n"
    " of calculating correlations, but not for determining symmetry.\n"
    "-The atomic symbol 'Vac' is used to indicate a vacancy.\n"
    "\n"
    "In the structure file:\n"
    "-The atom type is just a single atomic symbol\n"
    " (which, of course, has to be among the atomic symbols given in the\n"
    " lattice file).\n"
    "-Vacancies do not need to be specified.\n"
    "\n"
    "Examples\n"
    "\n"
    "The fcc lattice of the Cu-Au system:\n"
    "1 1 1 90 90 90\n"
    "0   0.5 0.5\n"
    "0.5 0   0.5\n"
    "0.5 0.5 0\n"
    "0 0 0 Cu,Au\n"
    "\n"
    "The Cu3Au L1_2 structure:\n"
    "1 1 1 90 90 90\n"
    "1 0 0\n"
    "0 1 0\n"
    "0 0 1\n"
    "0   0   0   Au\n"
    "0.5 0.5 0   Cu\n"
    "0.5 0   0.5 Cu\n"
    "0   0.5 0.5 Cu\n"
    "\n"
    "A lattice for the Li_x Co_y Al_(1-y) O_2 system:\n"
    "0.707 0.707 6.928 90 90 120\n"
    " 0.3333  0.6667 0.3333\n"
    "-0.6667 -0.3333 0.3333\n"
    " 0.3333 -0.3333 0.3333\n"
    " 0       0      0       Li,Vac\n"
    " 0.3333  0.6667 0.0833  O\n"
    " 0.6667  0.3333 0.1667  Co,Al\n"
    " 0       0      0.25    O\n"
    "\n"
    "Symmetry file format (sym.out)\n"
    "\n"
    "[number of symmetry operations]\n"
    "\n"
    "3x3 matrix: point operation\n"
    "\n"
    "1x3 matrix: translation\n"
    "\n"
    "repeat, etc.\n"
    "\n"
    "Note that if you enter more than one unit cell of the lattice,\n"
    "sym.out will contain some pure translations as symmetry elements.\n"
    "\n"
    "Cluster file format (clusters.out)\n"
    "\n"
    "for each cluster:\n"
    "[multiplicity]\n"
    "[length of the longest pair within the cluster]\n"
    "[number of points in cluster]\n"
    "[coordinates of 1st point] [number of possible species-2] [cluster "
    "function]\n"
    "[coordinates of 2nd point] [number of possible species-2] [cluster "
    "function]\n"
    "etc.\n"
    "\n"
    "repeat, etc.\n"
    "\n"
    "(Multiplicity and length are ignored when reading in the clusters.out "
    "file.)\n"
    "For each 'point' the following convention apply\n"
    "-The coordinates are expressed in the coordinate system given in\n"
    " the first line (or the first 3 lines) of the lat.in file.\n"
    "-The 'number of possible species' distinguishes between binaries, "
    "ternaries, etc...\n"
    " Since each site can accomodate any number of atom types,\n"
    " this is specified for each point of the cluster.\n"
    "-In multicomponent system, the cluster function are numbered from 0 to "
    "number of possible species-2.\n"
    "In the simple of a binary system [number of possible species-2] [cluster "
    "function] are just 0 0.\n"
    "For a ternary, the possible values are 1 0 and 1 1.\n"
    "All the utilities that are not yet multicomponent-ready just ignore the "
    "entries [number of possible species-2] [cluster function].\n"
    "\n"
    "Convention used to calculate the correlations:\n"
    "  The cluster functions in a m-component system are defined as\n"
    "   function '0' : -cos(2*PI*  1  *s/m)\n"
    "   function '1' : -sin(2*PI*  1  *s/m)\n"
    "                             .\n"
    "                             .\n"
    "                             .\n"
    "                  -cos(2*PI*[m/2]*s/m)\n"
    "                  -sin(2*PI*[m/2]*s/m)   <--- the last sin( ) is omitted "
    "if m is even\n"
    "  where the occupation variable s can take any values in {0,...,m-1}\n"
    "  and [...] denotes the 'round down' operation.\n"
    "  Note that, these functions reduce to the single function (-1)^s in the "
    "binary case.\n"
    "\n"
    "Special options:\n"
    "\n"
    "-sym:  Just find the space group and then abort.\n"
    "-clus: Just find space group and clusters and then abort.\n"
    "-z:    To find symmetry operations, atoms are considered to lie on\n"
    "       top of one another when they are less than this much apart.\n"
    "-sig:  Number of significant digits printed.\n"
    "\n"
    "\n"
    "->Cautions\n"
    "\n"
    "When vacancies are specified, the program may not be able to warn\n"
    "you that the structure and the lattice just don't fit.\n"
    "Carefully inspect the corrdump.log file!\n"
    "\n"
    "If the structure has significant cell shape relaxations, the program\n"
    "will be unable to find how it relates to the ideal lattice.\n"
    "The problem gets worse as the supercell size of the structure gets\n"
    "bigger.\n"
    "\n"
    "There is no limit on how far an atom in a structure can be from\n"
    "the ideal lattice site. The program first finds the atom that can\n"
    "be the most unambiguously assigned to a lattice site. It then\n"
    "finds the next best assignement and so on. This is actually a\n"
    "pretty robust way to do this. But keep in mind that the -z option\n"
    "does NOT control this process.\n";
