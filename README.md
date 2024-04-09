mpic++ -o convex_hull_mpi main.cpp  
mpirun --oversubscribe -np 8  ./convex_hull_mpi
