# EryNeu
EryNeu package includes code for modeling cellular differentiation using gene circuit methodology. The code was used in the manuscript by Handzlik and Manu (Department of Biology, University of North Dakota) for modeling the differentiation of erythrocytes and neutrophils.

### DEPENDENCIES

1. GNU compiler selection (gcc)
2. GNU scientific library (gsl)

### INSTALLATION/COMPILATION (Unix/Linux)

Commands for building the executables files:

make -e veryclean (optional)

make -e deps

make -e fly

### COMMAND REFERENCE

    Usage: fly_sa [-a <accuracy>] [-b <bkup_freq>] [-B] [-e <freeze_crit>] [-E]
              [-f <param_prec>] [-g <g(u)>] [-h] [-i <stepsize>]    
              [-l] [-L] [-n] [-N] [-p] [-Q] [-s <solver>] [-t] [-v]
              [-w <out_file>] [-y <log_freq>]
              <datafile>
                  
    Usage: fly_sa [options] <datafile>

    Argument:
      <datafile>          input data file

    Options:
      -a <accuracy>       solver accuracy for adaptive stepsize ODE solvers
      -b <bkup_freq>      write state file every <bkup_freq> * tau moves
      -B                  run in benchmark mode (only do fixed initial steps)
      -D                  debugging mode, prints all kinds of debugging info
      -e <freeze_crit>    set annealing freeze criterion to <freeze_crit>
      -E                  run in equilibration mode
      -f <param_prec>     float precision of parameters is <param_prec>
      -g                  chooses g(u) in the right hand side of the equation:                                              
                          -e:  g(u) = 1 / (1 + exp(-2u))     (exponential) 
                          -h:  g(u) = -1 if x<0, 0 if x>0      (heaviside) 
                          -s:  g(u) = 1/2 * ( u / sqrt(1 + u^2) )   (sqrt)
                          -t:  g(u) = 1/2 * ( tanh(u) + 1)          (tanh) 
      -h                  prints this help message
      -i <stepsize>       sets ODE solver stepsize (in minutes)
      -k <cost_function> select error function for calculation the differences between
                          equation and data: sq=squared differences, abs=absolute differences
                          default=sq
      -l                  echo log to the terminal
      -n                  nofile: don't print .log or .state files
      -N                  generates landscape to .landscape file in equilibrate mode 
      -o                  use oldstyle cell division times (3 div only)
      -p                  prints move acceptance stats to .prolix file
      -Q                  quenchit mode, T is lowered immediately to zero
      -s <solver>         choose ODE solver
      -t                  write timing information to .times file
      -v                  print version and compilation date
      -w <out_file>       write output to <out_file> instead of <datafile>
      -y <log_freq>       write log every <log_freq> * tau moves


### INPUT FILE FORMAT

The 'examples' directory of this package includes examples of EryNeu input data files.
