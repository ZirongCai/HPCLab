
void GEMM(double* A, double* B, double* C, double** A_pack, double** B_pack, int S, int nTeams = 1, int threadsPerTeam = 1);

void GEBP(double* A, double* B_pack, double* C, double* A_pack, int S, int threadsPerTeam);

void GEPP(double* A, double* B, double* C, double** A_pack, double** B_pack, int S, int nTeams, int threadsPerTeam);

void microkernel(double* A_pack, double* B_pack, double* C, int S);
