#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <assert.h>
#include <time.h>
#include <vector>
#include <math.h>
#include <fstream>

#define NumsubR	2
#define Nnot	20
#define Nodv	100
#define NumMacro 100000
#define NumSys	1
#define NumConstraint	1
#define NumThreshold	2

#define norm 2.328306549295728e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0


double MRG32k3a(void);  //Generate R(0,1) by L'ecuyer (1997)
double  s10 = 12345, s11 = 12345, s12 = 12345, s20 = 12345, s21 = 12345, s22 = 12345;
//double  s10 = 43, s11 =54, s12 =65, s20 =43, s21 =54, s22 =65;
//double  s10 = 4321111, s11 =1115432, s12 =1116543, s20 =4321111, s21 =1115432, s22 =6543111;
//double  s10 = 43221, s11 =54332, s12 =65443, s20 =43321, s21 =54532, s22 =61543;
//double  s10 = 1010, s11 =10, s12 =101, s20 =2001, s21 = 202, s22 = 202;

double minfn(double x1, double x2);
double maxfn(double x, double y);

double normal(double rmean, double rvar);
double configuration(void);
double generate_multiNormal(int numConstraint, double rho, int case_index);
int read_chol_matrix(void);

double mean_value[NumSys][NumConstraint];
double variance_value[NumSys][NumConstraint];
double chol_matrix[3][NumConstraint][NumConstraint];
int system_info[NumSys][2];
double obs_for_each_sys[NumSys];
double pcd_for_each_sys[NumSys];

double observations[NumSys][NumConstraint];
double q[1][NumThreshold][NumConstraint];
double epsilon[NumConstraint];
double num_obs[NumSys][NumConstraint];
double v_min[NumSys][NumConstraint];
double v_max[NumSys][NumConstraint];
double R[NumSys][NumConstraint];
int ON[NumSys][NumConstraint];
int ON_l[NumSys][NumConstraint][NumThreshold];
int Z[NumSys][NumConstraint][NumThreshold];

int true_feasibility[NumSys][NumThreshold];
int feasibility[NumSys][NumThreshold];

//double eta[NumSys];
double eta[NumConstraint];

double rho = -0.25;

double total_obs = 0;
double final_cd = 0;
double correct_decision = 0;
double overall_obs = 0;

int main()
{
  read_chol_matrix();

  //eta[0] = 0.1371;
  eta[0] = 0.1854;
  //eta[1] = 0.2930;
  //eta[2] = 0.3119;

  for (int i=0; i<NumSys; i++) {
    obs_for_each_sys[i] = 0;
    pcd_for_each_sys[i] = 0;
  }

  for (int l=0; l<NumMacro; l++) {

      configuration();
      total_obs = 0;
      final_cd = 1;
      double Sil2[NumSys][NumConstraint];

      for (int i=0; i<NumSys; i++) {
   			for (int j=0; j<NumConstraint; j++) {
                Sil2[i][j] = 0;
                num_obs[i][j] = 0;
   			}
        }

      for (int i=0; i<NumSys; i++) {

            double obs_for_this_sys = 0;
            double correct_decision_for_this_sys = 0;

   			// generate initial samples for obs
   			double sumY[NumConstraint];
   			double sum_squareY[NumConstraint];

            int surviveConstraint = NumConstraint;
            int surviveThreshold[NumConstraint];

   			for (int j=0; j<NumConstraint; j++) {
   				sumY[j] = 0;
   				sum_squareY[j] = 0;
                surviveThreshold[j] = NumThreshold;
   			}

   			for (int n=0; n<Nnot; n++) {
   				generate_multiNormal(NumConstraint, rho, 1);
                obs_for_each_sys[i] += 1;
   				total_obs += 1;

   				for (int j=0; j<NumConstraint; j++) {
   					sumY[j] += observations[i][j];
   					sum_squareY[j] += observations[i][j] * observations[i][j];
                    num_obs[i][j] += 1;
   				}

   			}

   			// find continuation region
   			for (int j=0; j<NumConstraint; j++) {
   				Sil2[i][j] = (sum_squareY[j]/(Nnot-1)) - (sumY[j]/(Nnot-1))*(sumY[j]/Nnot);
   				R[i][j] = maxfn(0, (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j]-epsilon[j]*num_obs[i][j]/2);
                //R[i][j] = (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j];
   			}

   			while (surviveConstraint != 0) {

   				for (int j=0; j<NumConstraint; j++) {

   					if (ON[i][j] == 1) {
   						v_min[i][j] = minfn(v_min[i][j], (sumY[j]+R[i][j])/num_obs[i][j]);
   						v_max[i][j] = maxfn(v_max[i][j], (sumY[j]-R[i][j])/num_obs[i][j]);

   						for (int d=0; d<NumThreshold; d++) {

   							if (ON_l[i][j][d] == 1) {

                                if (v_min[i][j] <= q[system_info[i][0]][d][j]) {
   									Z[i][j][d] = 1;
   									ON_l[i][j][d] = 0;
   									surviveThreshold[j] -= 1;
                                }

                                if (v_max[i][j] >= q[system_info[i][0]][d][j]) {
      							    Z[i][j][d] = 0;
   									ON_l[i][j][d] = 0;
   				                    surviveThreshold[j] -= 1;
   								}
   							}

   					     }

                         if (surviveThreshold[j] == 0) {
                             ON[i][j] = 0;
                             surviveConstraint -= 1;
                         }

   					}

   				}

   				if (surviveConstraint == 0) break;

   				generate_multiNormal(NumConstraint, rho, 1);
                obs_for_each_sys[i] += 1;
   				total_obs += 1;

   				for (int j=0; j<NumConstraint; j++) {
   					sumY[j] += observations[i][j];
                    num_obs[i][j] += 1;
                    R[i][j] = maxfn(0, (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j]-epsilon[j]*num_obs[i][j]/2);
                    //R[i][j] = (Nnot-1)*Sil2[i][j]*(eta[j])/epsilon[j];
   				}
   			}

   			// declare feasibility of a system
   			int feasibility_for_threshold = 1;
   			for (int d=0; d<NumThreshold; d++) {
   				for (int j=0; j<NumConstraint; j++) {
   					feasibility_for_threshold *= Z[i][j][d];
   				}
   				if (feasibility_for_threshold == 1) {
   					feasibility[i][d] = 1;
   				} else {
   					feasibility[i][d] = 0;
   				}
   			}

   			// check whether the decision is correct
   			int cd_for_one_threshold = 1;
   			for (int d=0; d<NumThreshold; d++) {
   				for (int j=0; j<NumConstraint; j++) {
   					if (mean_value[i][j] <= q[system_info[i][0]][d][j]-epsilon[j]) {
   						if (Z[i][j][d] == 1) {
   							cd_for_one_threshold *= 1;
   						} else {
   							cd_for_one_threshold *= 0;
   						}
   					} else if (mean_value[i][j] >= q[system_info[i][0]][d][j]+epsilon[j]) {
   						if (Z[i][j][d] == 0) {
   							cd_for_one_threshold *= 1;
   						} else {
   							cd_for_one_threshold *= 0;
   						}
   					} else {
   						cd_for_one_threshold *= 1;
   					}
   				}
   			}

   			if (cd_for_one_threshold == 1) {
                final_cd *= 1;
                pcd_for_each_sys[i] += 1;
   			} else {

            final_cd *= 0;
            }

   		}

        if (final_cd == 1) {
            correct_decision++;
        }
   		overall_obs += total_obs;

    }

    printf("Overall: %.10f\n", correct_decision/NumMacro);
    printf("Overall: %.4f\n", overall_obs/NumMacro);

	return 0;
}

double MRG32k3a() //L'ecuyer Random number generator(0,1)
{
    long   k;
    double p1, p2;
    // Component 1
    p1 = a12 * s11 - a13n * s10;
    k = p1 / m1;   p1 -= k * m1;   if (p1 < 0.0) p1 += m1;
    s10 = s11;   s11 = s12;   s12 = p1;

    // Component 2
    p2 = a21 * s22 - a23n * s20;
    k  = p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
    s20 = s21;   s21 = s22;   s22 = p2;
    // Combination
    if (p1 <= p2) return ((p1 - p2 + m1) * norm);
    else return ((p1 - p2) * norm)+0.000001;
}

double normal(double rmean, double rvar)
/* return normal random variable with mean rmean and variance rvar. */
// this is modified for Fully Sequential Procedure with CRN
{
	double V1 = 0, V2 = 0, W = 2, Y = 0, X1 = 0;
	do {
		V1 = 2 * MRG32k3a() - 1;
		V2 = 2 * MRG32k3a() - 1;
     	W = pow(V1,2) + pow(V2,2);

	} while (W > 1);

	Y = sqrt( (-2.00 * log(W))/W );
	X1 = rmean + sqrt(rvar) * V1 * Y;
	return X1;
}

int read_chol_matrix() {
  double c1, c2, c3, c4, c5;
  char ch;

  std::ifstream myfile ("cholMatrix_rho0.txt");
  if (myfile.is_open()) {
    int case_counter = 0;
    int pair_counter = 0;
    while ( myfile >> c1 >> ch >> c2 >> ch >> c3 >> ch >> c4 >> ch >> c5 )
    {
      chol_matrix[case_counter][pair_counter][0] = c1;
      //chol_matrix[case_counter][pair_counter][1] = c2;
      //chol_matrix[case_counter][pair_counter][2] = c3;
      //chol_matrix[case_counter][pair_counter][3] = c4;
      //chol_matrix[case_counter][pair_counter][4] = c5;
      pair_counter++;
      if (pair_counter == NumConstraint) {
        case_counter++;
        pair_counter = 0;
      }
    }
    myfile.close();
  }
  return 0;
}

double generate_multiNormal(int numConstraint, double rho, int case_index) {
	double std_normal[numConstraint];
	for (int i=0; i<numConstraint; i++) {
		std_normal[i] = normal(0, 1);
	}

	// input the C matrix first calculated from MATLAB
	// observation = mu + CZ
	double C[numConstraint][numConstraint];

    switch (case_index) {
        case 1:

        for (int i=0; i<numConstraint; i++) {
            for (int j=0; j<numConstraint; j++) {
                C[i][j] = chol_matrix[0][i][j];
            }
        }
        break;

        case 2:

        for (int i=0; i<numConstraint; i++) {
            for (int j=0; j<numConstraint; j++) {
                C[i][j] = chol_matrix[1][i][j];
            }
        }
        break;

        case 3:

        for (int i=0; i<numConstraint; i++) {
            for (int j=0; j<numConstraint; j++) {
                C[i][j] = chol_matrix[2][i][j];
            }
        }
        break;
  }


	for (int i=0; i<NumSys; i++) {
		for (int j=0; j<numConstraint; j++) {
			double temp = 0;
			for (int k=0; k<numConstraint; k++) {
			 	temp += C[j][k]*std_normal[k];
			}
			observations[i][j] = mean_value[i][j]+temp;
		}
	}

	return 0;
}

double configuration(void) {

  //int sys_counter = 0;
  //for (int settingnum = 0; settingnum<4; settingnum++) {
  // for (int casenum = 1; casenum<4; casenum++) {
  //    system_info[sys_counter][0] = settingnum;
  //    system_info[sys_counter][1] = casenum;
  //    sys_counter++;
  //   }
  // }

  system_info[0][0] = 0;
  //system_info[1][0] = 0;
  //system_info[2][0] = 0;
  system_info[0][1] = 1;
  //system_info[1][1] = 2;
  //system_info[2][1] = 3;

	for (int j=0; j<NumConstraint; j++) {
		//epsilon[j] = 1/sqrt(Nnot);
		epsilon[j] = 0.1;
	}

	// single system
	for (int i=0; i<NumSys; i++) {
		for (int j=0; j<NumConstraint; j++) {
			mean_value[i][j] = 0;
			variance_value[i][j] = 1;

			ON[i][j] = 1;
			v_min[i][j] = 100000000000;
			v_max[i][j] = -10000000000;
			for (int d=0; d<NumThreshold; d++) {
				ON_l[i][j][d] = 1;
			}
		}
	}

	// q[0][0][0] = -1.25*epsilon[0];
 //  q[0][1][0] = -epsilon[0];
 //  q[0][2][0] = epsilon[0];
 //  q[0][3][0] = 1.25*epsilon[0];   // redundant

 //  q[0][0][1] = -1.25*epsilon[1];
 //  q[0][1][1] = -epsilon[1];
 //  q[0][2][1] = epsilon[1];
 //  q[0][3][1] = 1.25*epsilon[1];    // redundant

 //  q[0][0][2] = -1.25*epsilon[2];
 //  q[0][1][2] = -epsilon[2];
 //  q[0][2][2] = epsilon[2];
 //  q[0][3][2] = 1.25*epsilon[2];

 //  q[0][0][3] = -1.25*epsilon[3];
 //  q[0][1][3] = -epsilon[3];
 //  q[0][2][3] = epsilon[3];   // redundant
 //  q[0][3][3] = 1.25*epsilon[3];   // redundant

 //  q[0][0][4] = -1.25*epsilon[4];
 //  q[0][1][4] = -epsilon[4];
 //  q[0][2][4] = epsilon[4];   // redundant
 //  q[0][3][4] = 1.25*epsilon[4];   // redundant

  // for (int j=0; j<NumConstraint; j++) {
  //  q[0][0][j] = -epsilon[j];
  //   q[0][1][j] = epsilon[j];
  //   q[0][2][j] = -1.5*epsilon[j];
  //   q[0][3][j] = -1.25*epsilon[j];
  //   q[0][4][j] = -epsilon[j];
  //   q[0][5][j] = epsilon[j];
  //   q[0][6][j] = 1.25*epsilon[j];
  //   q[0][7][j] = 1.5*epsilon[j];
  //   q[0][8][j] = 1.75*epsilon[j];
  //   q[0][9][j] = 2*epsilon[j];
  // }

  q[0][0][0] = -epsilon[0];
  q[0][1][0] = epsilon[0];
  //q[0][2][0] = epsilon[0];
  //q[0][3][0] = 1.5*epsilon[0];

	return 0;
}

double maxfn(double x, double y)
{
	if(x>y) return x;
	else return y;
}

double minfn(double x, double y)
{
	if(x<y) return x;
	else return y;
}
