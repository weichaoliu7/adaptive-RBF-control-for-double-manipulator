#include <stdio.h>
#include <stdlib.h>
#include "sine.h"
#include "cosine.h"
#include "math.h"
#include "inv_matrix.h"

#define PI 3.14159

// reference: [1]Liu JinKun. Robot Control System Design and MATLAB Simulation[M]. Tsinghua University Press, 2008.
// [2]Feng G. A compensating scheme for robot tracking based on neural networks[J]. Robotics and Autonomous Systems, 1995, 15(3): 199-206.

void matrix_Multi(double C[][4], double A[][4], double B[][4], int rows1, int cols1, int cols2){
    for (int j = 0; j < rows1; j++){
        for (int k = 0; k < cols2; k++){
            C[j][k] = 0.0;
            for (int g = 0; g < cols1; g++){
                C[j][k] += A[j][g] * B[g][k];
            }
        }
    }
}

void matrix_Multi1(double C[][2], double A[][4], double B[][2], int rows1, int cols1, int cols2){
    for (int j = 0; j < rows1; j++){
        for (int k = 0; k < cols2; k++){
            C[j][k] = 0.0;
            for (int g = 0; g < cols1; g++){
                C[j][k] += A[j][g] * B[g][k];
            }
        }
    }
}

void matrix_Multi2(double C[2], double A[][2], double B[2], int rows1, int cols1, int cols2){
    for (int j = 0; j < rows1; j++){
        C[j] = 0.0;
        for (int k = 0; k < cols1; k++){
            C[j] += A[j][k] * B[k];
        }
    }
}

typedef struct{
    double qd1_p;
    double qd2_p;
    double dqd1_p;
    double dqd2_p;
    double ddqd1_p;
    double ddqd2_p;
} Variable;


int main(){

    double Ts = 0.01; // sampling period
    double t0 = 0.0;
    double t1 = 20.0;
    Variable a = {0.2, -0.2, 0.2 * 0.5 * PI, 0.2 * 0.5 * PI, -0.2 * pow(0.5 * PI, 2), 0.2 * pow(0.5 * PI, 2)};
    Variable m = {0.5 * PI, 0.5 * PI, 0.5 * PI, 0.5 * PI, 0.5 * PI, 0.5 * PI};
    Variable b1 = {1, 1, 0, 0, 0, 0};

    Data qd1, dqd2, ddqd1;
    Data1 qd2, dqd1, ddqd2;
    sine(&qd1, a.qd1_p, m.qd1_p, b1.qd1_p, Ts, t0, t1);           // desired angular displacement of link 1
    cosine(&qd2, a.qd2_p, m.qd2_p, b1.qd2_p, Ts, t0, t1);         // desired angular displacement of link 2
    cosine(&dqd1, a.dqd1_p, m.dqd1_p, b1.dqd1_p, Ts, t0, t1);     // desired angular velocity of link 1
    sine(&dqd2, a.dqd2_p, m.dqd2_p, b1.dqd2_p, Ts, t0, t1);       // desired angular velocity of link 2
    sine(&ddqd1, a.ddqd1_p, m.ddqd1_p, b1.ddqd1_p, Ts, t0, t1);   // desired angular acceleration of link 1
    cosine(&ddqd2, a.ddqd2_p, m.ddqd2_p, b1.ddqd2_p, Ts, t0, t1); // desired angular acceleration of link 2

    int ARRAY_SIZE = (t1 - t0) / Ts;

    // parameter
    int IN = 4, H = 5, OUT = 2; // input, hidden and output layer cell's number
    int M1 = 2;                 // adaptive law selection
    int M = 3;                  // control scheme selection for precise model
    double gamma = 20;
    double alpha1 = 3, alpha2 = 3;
    double d1 = 2, d2 = 3, d3 = 6; // external perturbation factor
    double kp[2][2] = {{pow(alpha1, 2), 0}, {0, pow(alpha1, 2)}}; // proportional gain
    double kv[2][2] = {{2 * alpha2, 0}, {0, 2 * alpha2}};         // proportional gain
    double A[IN][IN], B[IN][OUT], Q[IN][IN];
    
    double ctrl_u1[ARRAY_SIZE], ctrl_u2[ARRAY_SIZE];
    double ctrl_u3, ctrl_u4, ctrl_u5, ctrl_u6, ctrl_u7, ctrl_u8;
    double q1[ARRAY_SIZE], dq1[ARRAY_SIZE], q2[ARRAY_SIZE], dq2[ARRAY_SIZE], ddq1[ARRAY_SIZE], ddq2[ARRAY_SIZE];
    double f1[ARRAY_SIZE], fn1[ARRAY_SIZE], f2[ARRAY_SIZE], fn2[ARRAY_SIZE], tol_1[ARRAY_SIZE], tol_2[ARRAY_SIZE];

    double time;
    double e1, e2, de1, de2;
    double f[2], fn[2], tol1[2], tol2[2], tol[2];
    double theta[H][OUT][ARRAY_SIZE], d_theta[H][OUT][ARRAY_SIZE];
    double x[IN], h[H], b[H], c[IN][H], h_x[H][IN], h_x_P[H][IN], h_x_P_B[H][OUT];
    double kv_de[2], kp_e[2], diff[2], diff1[2], D0_diff[2], C0_dq[2], dD_ddq[2], dC_dq[2];
    double sum1[2], D0_f[2], D0_fn[2], inv_D0[2][2], inv_D0_diff1[2];

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            if (j - i == 2){
                A[i][j] = 1.0;
            }else if (i >= 2 && j < 2){
                A[i][j] = -kp[i - 2][j];
            }else if (i >= 2 && j >= 2){
                A[i][j] = -kv[i - 2][j - 2];
            }else{
                A[i][j] = 0.0;
            }
        }
    }

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 2; j++){
            if (i - j == 2){
                B[i][j] = 1.0;
            }else{
                B[i][j] = 0.0;              
            }
        }
    }

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 2; j++){
            if (i == j){
                Q[i][j] = 50;
            }else{
                Q[i][j] = 0.0;
            }
        }
    }

    for (int i = 0; i < ARRAY_SIZE; i++){

        time = (i + 1) * Ts + t0;
        printf("time at step %d: %f\n", i, time);

        if (i == 0){
            for (int j = 0; j < H; j++){
                for (int k = 0; k < OUT; k++){
                    theta[j][k][i] = 0.1; // weight of RBF network
                    d_theta[j][k][i] = 0.0; // weight derivative of RBF network
                }
            }
            q1[i] = 0.6;
            dq1[i] = 0.3;
            ddq1[i] = 0.0;
            q2[i] = 0.5;
            dq2[i] = 0.5;
            ddq2[i] = 0.0;
        }

        // controller
        ctrl_u1[i] = qd1.y[i]; // desired angular displacement of link 1
        ctrl_u2[i] = qd2.y[i]; // desired angular displacement of link 2

        ctrl_u3 = q1[i];   // angular displacement of link 1
        ctrl_u4 = dq1[i];  // angular velocity of link 1
        ctrl_u5 = q2[i];   // angular displacement of link 2
        ctrl_u6 = dq2[i];  // angular velocity of link 2
        ctrl_u7 = ddq1[i]; // angular acceleration of link 1
        ctrl_u8 = ddq2[i]; // angular acceleration of link 2

        e1 = ctrl_u3 - ctrl_u1[i]; // angular displacement error of link 1
        e2 = ctrl_u5 - ctrl_u2[i]; // angular displacement error of link 2
        de1 = ctrl_u4 - dqd1.y[i]; // angular velocity error of link 1
        de2 = ctrl_u6 - dqd2.y[i]; // angular velocity error of link 2

        x[0] = e1; // state variable
        x[1] = e2;
        x[2] = de1;
        x[3] = de2;

        double dq[2] = {dq1[i], dq2[i]};           // angular velocity
        double ddq[2] = {ddq1[i], ddq2[i]};        // angular acceleration
        double ddqd[2] = {ddqd1.y[i], ddqd2.y[i]}; // desired angular acceleration

        double column_c[] = {-2, -1, 0, 1, 2};
        for (int j = 0; j < IN; j++){
            for (int k = 0; k < H; k++){
                c[j][k] = column_c[k]; // basis function center
            }
        }

        for (int j = 0; j < H; j++){
            b[j] = 3.0; // basis function width
        }

        for (int j = 0; j < H; j++){
            double sum = 0.0;
            for (int k = 0; k < IN; k++){
                sum += pow(x[k] - c[k][j], 2);
            }
            h[j] = exp(-sum / (2 * b[j] * b[j])); // RBF function's output
        }

        double P[4][4] = {{58.3333, 0, 2.7778, 0},
                          {0, 58.3333, 0, 2.7778},
                          {2.7778, 0, 4.6296, 0},
                          {0, 2.7778, 0, 4.6296}}; // use matlab get the solution of Lyapunov equation: P=lyap(A',Q);

        // calculate weight estimate's derivative
        for (int j = 0; j < H; j++){
            for (int k = 0; k < IN; k++){
                h_x[j][k] = h[j] * x[k];
                //printf("h_x[%d][%d]: %f\n", j, k, h_x[j][k]);
            }
        }

        matrix_Multi(h_x_P, h_x, P, H, IN, IN);
        // for (int j = 0; j < H; j++){
        //     for (int k = 0; k < IN; k++){
        //         printf("h_x_P[%d][%d]: %f\n", j, k, h_x_P[j][k]);
        //     }
        // }

        matrix_Multi1(h_x_P_B, h_x_P, B, H, IN, OUT);

        if (M1 == 1){ // adaptive law
            for (int j = 0; j < H; j++){
                for (int k = 0; k < OUT; k++){
                    d_theta[j][k][i] = gamma * h_x_P_B[j][k]; // derivative of weight estimate
                }
            }
        }
        else if (M1 == 2){ // adaptive law with UUB
            double k1 = 0.001;
            double sum = 0.0;
            for (int j = 0; j < H; j++){
                for (int k = 0; k < OUT; k++){
                    sum += (pow(theta[j][k][i], 2));
                }
            }

            // for (int j = 0; j < H; j++){
            //     for (int k = 0; k < OUT; k++){
            //         printf("theta[%d][%d][%d]: %f\n", j, k, i, theta[j][k][i]);
            //     }
            // }  

            for (int j = 0; j < H; j++){
                for (int k = 0; k < OUT; k++){
                    d_theta[j][k][i] = gamma * h_x_P_B[j][k] + k1 * gamma * sqrt(sum) * theta[j][k][i]; // network weight estimate derivative
                    // printf("d_theta[%d][%d][%d]: %f\n", j, k, i, d_theta[j][k][i]);
                }
            }

            for (int j = 0; j < H; j++){
                for (int k = 0; k < OUT; k++){
                    theta[j][k][i + 1] = theta[j][k][i] + (d_theta[j][k][i] * Ts);
                }
            }

        }

        double g = 9.8; // gravitational acceleration
        double v = 13.33;
        double q01 = 8.98;
        double q02 = 8.75;

        double D0[2][2] = {{v + q01 + 2 * q02 * cos(ctrl_u5), q01 + q02 * cos(ctrl_u5)},
                           {q01 + q02 * cos(ctrl_u5), q01}};  // inertia matrix of nominal model
        double C0[2][2] = {{-q02 * ctrl_u6 * sin(ctrl_u5), -q02 * (ctrl_u4 + ctrl_u6) * sin(ctrl_u5)},
                           {q02 * ctrl_u4 * sin(ctrl_u5), 0}};  // coriolis term of the nominal model
        double G0[2] = {15 * g * cos(ctrl_u3) + 8.75 * g * cos(ctrl_u3 + ctrl_u5), 8.75 * g * cos(ctrl_u3 + ctrl_u5)}; // gravity term of the nominal model
        double d_D[2][2], d_C[2][2], d_G[2];

        for (int j = 0; j < 2; j++){
            for (int k = 0; k < 2; k++){
                d_D[j][k] = 0.2 * D0[j][k];
                d_C[j][k] = 0.2 * C0[j][k];
            }
            d_G[j] = 0.2 * G0[j];
        }

        double D = D0 - d_D;  // inertia matrix of precise model
        double C = C0 - d_C;  // coriolis term of precise model
        double G = G0 - d_G;  // gravity term of precise model
        double d = d1 + d2 * sqrt(pow(e1, 2) + pow(e2, 2)) + d3 * sqrt(pow(de1, 2) + pow(de2, 2)); // external disturbance

        double e[2] = {e1, e2}, de[2] = {de1, de2};

        matrix_Multi2(kv_de, kv, de, 2, 2, 1);
        matrix_Multi2(kp_e, kp, e, 2, 2, 1);

        for (int j = 0; j < OUT; j++){
            diff[j] = ddqd[j] - kv_de[j] - kp_e[j];
        }

        matrix_Multi2(D0_diff, D0, diff, 2, 2, 1);
        matrix_Multi2(C0_dq, C0, dq, 2, 2, 1);

        inv_matrix(inv_D0, D0, 2);

        matrix_Multi2(dD_ddq, d_D, ddq, 2, 2, 1);
        matrix_Multi2(dC_dq, d_C, dq, 2, 2, 1);

        for (int j = 0; j < OUT; j++){
            sum1[j] = dD_ddq[j] + dC_dq[j] + d_G[j] + d;
        }

        matrix_Multi2(f, inv_D0, sum1, 2, 2, 1); // model uncertainty term

        // for (int j = 0; j < OUT; j++){
        //     printf("f[%d]: %f\n", j, f[j]);
        // }

        if (M == 1){ // nominal model based controller
            for (int j = 0; j < OUT; j++){
                tol1[j] = D0_diff[j] + C0_dq[j] + G0[j]; // computed torque controller
                tol2[j] = 0;
                tol[j] = tol1[j] + tol2[j];
            }
        }
        else if (M == 2){ // control with precise nonlinear compensation
            matrix_Multi2(D0_f, D0, f, 2, 2, 1);
            for (int j = 0; j < OUT; j++){
                tol1[j] = D0_diff[j] + C0_dq[j] + G0[j];
                tol2[j] = -D0_f[j];
                tol[j] = tol1[j] + tol2[j]; // improved computed torque controller
            }
        }
        else if (M == 3){ // control with neural compensation
            for (int k = 0; k < OUT; k++){
                fn[k] = 0.0;
                for (int j = 0; j < H; j++){
                    fn[k] += theta[j][k][i] * h[j]; // NN's output = transposition of weight estimate * RBF function's output
                }
                // printf("fn[%d]: %f\n", k, fn[k]);
            }

            matrix_Multi2(D0_fn, D0, fn, 2, 2, 1);

            for (int j = 0; j < OUT; j++){
                tol1[j] = D0_diff[j] + C0_dq[j] + G0[j];                
                tol2[j] = -D0_fn[j];
                tol[j] = tol1[j] + tol2[j];
                // printf("tol1[%d]: %f\n", j, tol1[j]);
            }
        }

        tol_1[i] = tol[0];
        tol_2[i] = tol[1];
        f1[i] = f[0];
        f2[i] = f[1];
        fn1[i] = fn[0];
        fn2[i] = fn[1];

        // plant
        for (int j = 0; j < OUT; j++){
            diff1[j] = tol[j] - C0_dq[j] - G0[j];
        }

        matrix_Multi2(inv_D0_diff1, inv_D0, diff1, 2, 2, 1);

        for (int j = 0; j < OUT; j++){
            ddq[j] = inv_D0_diff1[j] + f[j];
        }

        dq1[i + 1] = dq1[i] + ddq[0] * Ts;
        dq2[i + 1] = dq2[i] + ddq[1] * Ts;
        q1[i + 1] = q1[i] + dq1[i] * Ts;
        q2[i + 1] = q2[i] + dq2[i] * Ts;
        ddq1[i + 1] = ddq[0];
        ddq2[i + 1] = ddq[1];

        q1[0] = 0.6;
        dq1[0] = 0.3;
        ddq1[0] = 0.0;
        q2[0] = 0.5;
        dq2[0] = 0.5;
        ddq2[0] = 0.0;
    }

    const char* filename[] = {
        "qd1.txt", "qd2.txt", "q1.txt", "q2.txt",
        "tol1.txt", "tol2.txt", "f1.txt", "f2.txt",
        "fn1.txt", "fn2.txt", "dqd1.txt", "dqd2.txt",
        "dq1.txt", "dq2.txt", "ddqd1.txt", "ddqd2.txt",
        "ddq1.txt", "ddq2.txt"
    };

    double* dataArray[] = {
        qd1.y, qd2.y, q1, q2,
        tol_1, tol_2, f1, f2,
        fn1, fn2, dqd1.y, dqd2.y,
        dq1, dq2, ddqd1.y, ddqd2.y,
        ddq1, ddq2
    };

    for (int i = 0; i < 18; i++) {
        FILE* file = fopen(filename[i], "w");
        if (file != NULL) {
            for (int j = 0; j < ARRAY_SIZE; j++) {
                fprintf(file, "%lf\n", dataArray[i][j]);
            }
            fclose(file);
        } else {
            printf("Unable to open %s\n", filename[i]);
        }
    }
  
    return 0;

}
