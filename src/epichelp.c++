const char *helpstring = "";
/*
    "
    "Summary
    "
    "This code find a saddle using a modification of the dimer method called the
   epicycle method. "An epicycle differs from a dimer in that one of the system
   images lies in the center and the "other one is rotating around it while in a
   dimer the two images lie symmetrically about the center. "The algorithm is
   described in: "A. van de Walle, S. Kadkhodaei, R. Sun, and Q.-J. Hong.,
   'Epicycle method for elasticity limit calculations.', "Phys. Rev. B,
   95:144113, 2017. http://dx.doi.org/10.1103/PhysRevB.95.144113
    "
    "It is implemented as two nested optimization problems:
    "an inner optimization problem to find the softest phonon mode using a
   modification of the dimer method
    "(called the epicycle method) and an outer optimization problem to find the
   actual atomic geometry.
    "
    "All options that pertain to the inner 'epicycle' optimization start with
   '-e' while "those that control ionic (and cell parameter) movement start with
   '-i'. "Both optimization are implemented via a conjugate gradient algorithm.
    "
    "Parameters controlling the optimization
    "
    "-el  (Epicycle Length, in Angstrom) specifies the distance between the two
   images of the epicycle.
    "-egt (Epicycle Gradient Tolerance, in eV/Ang^2) is the stopping criterion
   for the soft mode search algorithm.
    "-ets (Epicycle Trial Step, in radiant) trial step for the line
   minimizations. "     If this trial step fails to bracket the minimum, the
   step is size multiplied by the...
    "-eml (Epicycle Multiplier in Line minimization) until the minimum is
   bracketed. "     But the step size will never exceed the...
    "-ems (Epicycle Maximum Step multiplier).
    "     The actual maximum step is equal to trial step times the multiplier.
    "     Bracketing stops when a bracket is found or the
    "-ebl (Epicycle maximum number of Bracketing steps in Line minimizations) is
   reached. "     Once the bracketing step is done, the algorithm does a sequant
   search, finding the point most likely to have zero derivative "     based on
   the derivatives at the end points of the bracket. If the predicted point
   falls too close to the endpoints "     (not in the middle half of the
   bracket) a bisection step is taken instead, except the first 't' times this
   occurs, "     where 't' is the...
    "-ebf (Epicycle motion Bad step Forgiveness).
    "     This is repeated until either the gradiant tolerance is met (-egt
   divided by sqrt(number of degrees of freedom ) or the...
    "-eil (Epicycle maximum number of Iterations in Line minimizations) is
   reached. "     The line minimizations are embedded in a conjugate gradient
   proceedure which is repeated until the gradient "     tolerance (-egt) is met
   or the ...
    "-eic (Epicycle maximum number of Iterations in Conjugate gradient) is
   reached. "     As in any conjugate gradient method, the conjugate direction
   reset to steepest descent every...
    "-eir (Epicycle iterations between conjugate gradient reset) steps.
    "
    "For the outer optimization problem (ion motion), the same parameters are
   available: with the initial '-e' replaced by '-i':
    "-igt,-its,-iml,-ims,-ibl,-ibf,-iil,-iic,iir.
    "     One difference is that -igt is in eV/Ang and -its and -ims are in
   Angstrom. "     In addition, the outer optimization start with a few steps of
   steepest descent, indicated by the
    "-isi (Ion motion Steepest maximum number of Iterations) parameter.
    "     This helps because the outer optimization problem can be far from
   quadratic in behavior initially. "     The step size is equal to the force
   times the
    "-ism (Ion motion Steepest descent Multiplier, in Angstrom^2/eV).
    "     The outer optimization problem also optimizes cell shape in addition
   to ionic coordinates. "     To make cell degress of freedom comparable to
   atomic coordinates, the following scaling is used: "     The 'force' on the
   cell is Omega^(2/3) * (stress tensor) * f where Omega is the average atomic
   volume for the "     initial configuration and f is a constant set with: "
   -ff (Force scale Factor, dimensionless). "     The cell shape is parametrized
   by n*Omega^(1/3)*strain/f where n is the number of atoms. "     This
   convention ensures that stress and cell parameters scale approriately with
   cell size and have units and magnitudes "     comparable to the atomic
   coordinates. Note that the atomic coordinates are stored internally as the
   coordinates for the "     initial cell shape (and are distorded according the
   strain acting on the cell).
    "-ds  indicates whether and how to relax the cell shape. 0 (the default)
   means no cell relaxations, "     3 means fully general 3d strains are
   allowed, 2 means means fully general 2d strains, etc. "     The unique axis
   is x be default but can be changed by adding 4 (for y axis) or 8 (for z
   axis).
    "-st (Sleep Time between read access) is the time between each attempt to
   check if the calculation lauched has completed.
    "
    "Input files
    "  str.in  The initial structure geometry in standard ATAT format (see mmaps
   -h) " Optionally the following files can be read to continue a previous run
   (see below for a description): "  epidir.out "  epipos.out
    "
    "Final output files
    "  cstr_relax.out  The optimized inflection point geometry at the end of the
   calculations. "  cenergy.out     The corresponding final energy. "  stdout
   Log file (written to infdet.log if called by robustrelax_xxxx)
    "
    "A few notes about the log file (assuming it is redirected to infdet.log):
    "  -The inner optimization problem output is bracketted by the phrases
   'begin on_sphere' and 'end on_sphere'. "  -For a quick overview of the
   progress, use:      grep curv infdet.log "  -For a more detailed view (each
   c.g. step), use: grep norm infdet.log "   The inner optimization steps are
   indicated by s_gnorm while the outer by l_gnorm. "  -To see each function
   evaluation, use:           grep dfx infdet.log "  -To see the whole history
   of epicycle positions and directions: " grep -e epipos -e epidir infdet.log
    "
    "Intermediate input/output files
    " The infdet code relies on an external program to obtain its energy and
   force data. At each optimization step, it writes: "  busy            A
   semaphore file indicating that the external program should run now. " str.out
   The current geometry to be run by the external code " When done, the external
   program deletes the busy file and infdet knows it has to read: "  force.out
   the forces acting on each atom (one 3d vector per line) "  stress.out the
   stress tensor acting on the cell "  energy          the total energy "
   str_relax.out   the atomic coordinate (same as str.out but perhaps re-ordered
   - forces will be re-ordered accordingly)
    "
    " As the optimization progresses, the current status is written to various
   files: "  str_current.out   the current best approximation to the inflection
   point (in standard ATAT format) "  epipos.out        the current best
   approximation to the inflection point (as a vector, in internal format) "
   epidir.out        the current direction of the epicycle (softest phonon mode)
    "                    Note that the whole history of these values are
   available in the log file, prefixed by "                    epipos= and
   epidir= , which can be used to restart an aborted or misbehaving run.
    "
    ;
		*/
