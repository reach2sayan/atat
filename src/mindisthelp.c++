const char *helpstring=""
"Returns 1 if some atoms are within r of each other and returns 0 otherwise.\n"
"This convention allows one to write:\n"
"  mindist -r=1.5 -s=str.out && runstruct_xxxx\n"
"to mean that the command runstruct_xxxx will run only if the minimum distance\n"
"between atoms is at least r.\n"
;
