#include "stdio.h"
#include "stdlib.h"
#include <fstream>
#include <iostream>
#include "math.h"
#include "string.h"

#define X0 0
#define Y0 1
#define Z0 2
#define X1 3
#define Y1 4
#define Z1 5

#define ROW 45
#define COL 91
#define X_FACTOR 0.555e-10
#define Y_FACTOR 1.111e-10

#define MIN 1.e-15
#define MAX 2.5e-9

#define AVO_NUM 6.0221367e23
#define DENSITY 2.888
#define MOLAR_M 131.293

#define DIM 3
#define FORCES 2

double rand_uniform();

void ClearTheScreen();

double VonNeumann(double xMin, double xMax, double yMin, double yMax);

using namespace std;

int main(int argc, char **argv) {

    const long BODIES =
            long(floor((DENSITY / MOLAR_M) * AVO_NUM * pow(1e2, 3.) * pow(2. * MAX, 3.) + 0.5));
    double r1[BODIES][DIM], a1[BODIES][DIM], v2[BODIES][DIM];
    double del, dt, max, mass[BODIES];
    int grid[ROW][COL];
    int coord[BODIES][DIM];
    double Norms[BODIES];

    cout << "time_step [s] = ";
    cin >> dt;
    cout << "stop-time [s] = ";
    cin >> max;

    double equil[FORCES] = {0., 0.}; // Not used here, spring equilib for springy forces
    double G_Newton[FORCES] = {-3.18e-132, 3.35e-76}; // Repulsive and attractive terms in LJ for Xenon
    double powerLaw[FORCES] = {-13., -7.};
    int i, j, k, q;
    double cost, sint, phi, sinp, cosp, vx, vy, vz, vNorm;
    double rMin[BODIES][6]; // Find distance of 1 Xe atom to its nearest neighbors
    for (i = 0; i < BODIES; i++) {
        rMin[i][0] = MAX;
        rMin[i][1] = MAX;
        rMin[i][2] = MAX;
        rMin[i][3] = MAX;
        rMin[i][4] = MAX;
        rMin[i][5] = MAX;
        mass[i] = 2.18e-25; //in kg (Xe)
        // Produce a random position in the box, i over bodies, j over dimensions (pretty much everywhere)
        for (j = 0; j < DIM; j++) { r1[i][j] = (2. * rand_uniform() - 1.) * MAX; }
        cost = 1. - 2. * rand_uniform();            // Random trajectory
        sint = sqrt((1. - cost) * (1. + cost));
        phi = 2. * M_PI * rand_uniform();
        sinp = sin(phi);
        cosp = cos(phi);
        vx = sint * cosp;                           // Random starting speed
        vy = sint * sinp;
        vz = cost;
        vNorm = rand_uniform() * 500.; // meters per sec.
        //vNorm = VonNeumann ( 4e1, 5e2, 0., 15e3 );
        v2[i][X0] = vx * vNorm;
        v2[i][Y0] = vy * vNorm;
        v2[i][Z0] = vz * vNorm;
    }

    double radius, cosT, sinT, cosP, sinP, delta[2], flip[DIM];
    double rAvg = 0.;
    // Loop over all the atoms, O(n**2)
    for (i = 0; i < BODIES; i++) {
        a1[i][X0] = 0.;
        a1[i][Y0] = 0.;
        a1[i][Z0] = 0.;
        for (j = 0; j < BODIES; j++) {
            radius = 0.;
            for (k = 0; k < DIM; k++) {
                delta[0] = fabs(r1[i][k] - r1[j][k]);
                flip[k] = 1.;
                delta[1] = 2. * MAX - delta[0];
                if (delta[1] > 2. * MAX || delta[1] <= 0) delta[1] = 2. * MAX;
                if (delta[1] < delta[0]) {
                    radius += pow(delta[1], 2.);
                    flip[k] = -1.; //wrap-around
                } else
                    radius += pow(delta[0], 2.);
            }
            radius = sqrt(radius);
            if (DIM == 3) {
                cosP = (r1[j][2] - r1[i][2]) /
                       sqrt(pow(r1[i][0] - r1[j][0], 2.) + pow(r1[i][1] - r1[j][1], 2.) + pow(r1[i][2] - r1[j][2], 2.));
                sinP = sqrt(1. - pow(cosP, 2.));
            } else { sinP = 1.; }
            cosT = (r1[j][0] - r1[i][0]) / sqrt(pow(r1[i][0] - r1[j][0], 2.) + pow(r1[i][1] - r1[j][1], 2.));
            sinT = (r1[j][1] - r1[i][1]) / sqrt(pow(r1[i][0] - r1[j][0], 2.) + pow(r1[i][1] - r1[j][1], 2.));
            if (r1[i][0] == r1[j][0] && r1[i][1] == r1[j][1]) {
                cosT = 1.0;
                sinT = 0.0;
            }
            if (i != j) {
                if (radius < rMin[i][0] && (r1[i][0] - r1[j][0]) > 0.) rMin[i][0] = radius;
                if (radius < rMin[i][1] && (r1[i][0] - r1[j][0]) < 0.) rMin[i][1] = radius;
                if (radius < rMin[i][2] && (r1[i][1] - r1[j][1]) > 0.) rMin[i][2] = radius;
                if (radius < rMin[i][3] && (r1[i][1] - r1[j][1]) < 0.) rMin[i][3] = radius;
                if (radius < rMin[i][4] && (r1[i][2] - r1[j][2]) > 0.) rMin[i][4] = radius;
                if (radius < rMin[i][5] && (r1[i][2] - r1[j][2]) < 0.) rMin[i][5] = radius;
                for (q = 0; q < FORCES; q++) {
                    a1[i][0] += G_Newton[q] * pow(radius - equil[q], powerLaw[q]) * cosT * sinP / mass[i] * flip[0];
                    a1[i][1] += G_Newton[q] * pow(radius - equil[q], powerLaw[q]) * sinT * sinP / mass[i] * flip[1];
                    if (DIM > 2)
                        a1[i][2] += G_Newton[q] * pow(radius - equil[q], powerLaw[q]) * cosP / mass[i] * flip[2];
                    if (G_Newton[q] == 6.673889e-11) {
                        a1[i][0] *= mass[j] * mass[i];
                        a1[i][1] *= mass[j] * mass[i];
                        a1[i][2] *= mass[j] * mass[i];
                    }
                } //F=m*a across all forces
            } //complete acceleration
        } //loop over masses
    } //loop over masses higher-level
    double time = 0.;
    long counter = 0;
    for (i = 0; i < BODIES; i++) {
        if (rMin[i][0] != MAX) {
            rAvg += rMin[i][0];
            counter++;
        }
        if (rMin[i][1] != MAX) {
            rAvg += rMin[i][1];
            counter++;
        }
        if (rMin[i][2] != MAX) {
            rAvg += rMin[i][2];
            counter++;
        }
        if (rMin[i][3] != MAX) {
            rAvg += rMin[i][3];
            counter++;
        }
        if (rMin[i][4] != MAX) {
            rAvg += rMin[i][4];
            counter++;
        }
        if (rMin[i][5] != MAX) {
            rAvg += rMin[i][5];
            counter++;
        }
        //cout << rMin[i][0] << "\t" << rMin[i][1] << "\t" << rMin[i][2] << "\t" << rMin[i][3] << "\t" << rMin[i][4] << "\t" << rMin[i][5] << endl;
    } //cout << counter << " " << rAvg/double(counter) << endl;

    while (time < max) { //dt = del * pow(radius[0][1],7.);

        for (i = 0; i < ROW; i++) {
            for (j = 0; j < COL; j++) {
                grid[i][j] = 0;
            }
        } //end of the initial grid setup

        for (i = 0; i < BODIES; i++) {
            rMin[i][0] = MAX;
            rMin[i][1] = MAX;
            rMin[i][2] = MAX;
            rMin[i][3] = MAX;
            rMin[i][4] = MAX;
            rMin[i][5] = MAX;
            for (j = 0; j < DIM; j++) {
                if (!time) r1[i][j] += 0.5 * v2[i][j] * dt + 0.25 * a1[i][j] * pow(dt, 2.);
                else r1[i][j] += 0.5 * v2[i][j] * dt;
            } //loop over dimensions
        } //loop over all of the bodies

        rAvg = 0.;
        for (i = 0; i < BODIES; i++) {
            a1[i][0] = 0.;
            a1[i][1] = 0.;
            a1[i][2] = 0.;
            for (j = 0; j < BODIES; j++) {
                radius = 0.;
                for (k = 0; k < DIM; k++) {
                    delta[0] = fabs(r1[i][k] - r1[j][k]);
                    flip[k] = 1.;
                    delta[1] = 2. * MAX - delta[0];
                    if (delta[1] > 2. * MAX || delta[1] <= 0) delta[1] = 2. * MAX;
                    if (delta[1] < delta[0]) {
                        radius += pow(delta[1], 2.);
                        flip[k] = -1.; //wrap-around
                    } else
                        radius += pow(delta[0], 2.);
                }
                radius = sqrt(radius);
                if (radius > 1e-9) continue;
                if (DIM == 3) {
                    cosP = (r1[j][2] - r1[i][2]) / sqrt(pow(r1[i][0] - r1[j][0], 2.) + pow(r1[i][1] - r1[j][1], 2.) +
                                                        pow(r1[i][2] - r1[j][2], 2.));
                    sinP = sqrt(1. - pow(cosP, 2.));
                } else { sinP = 1.; }
                cosT = (r1[j][0] - r1[i][0]) / sqrt(pow(r1[i][0] - r1[j][0], 2.) + pow(r1[i][1] - r1[j][1], 2.));
                sinT = (r1[j][1] - r1[i][1]) / sqrt(pow(r1[i][0] - r1[j][0], 2.) + pow(r1[i][1] - r1[j][1], 2.));
                if (r1[i][0] == r1[j][0] && r1[i][1] == r1[j][1]) {
                    cosT = 1.0;
                    sinT = 0.0;
                }
                if (i != j) {
                    if (radius < rMin[i][0] && (r1[i][0] - r1[j][0]) > 0. && fabs(r1[i][0]) < MAX &&
                        fabs(r1[j][0]) < MAX)
                        rMin[i][0] = radius;
                    if (radius < rMin[i][1] && (r1[i][0] - r1[j][0]) < 0. && fabs(r1[i][0]) < MAX &&
                        fabs(r1[j][0]) < MAX)
                        rMin[i][1] = radius;
                    if (radius < rMin[i][2] && (r1[i][1] - r1[j][1]) > 0. && fabs(r1[i][1]) < MAX &&
                        fabs(r1[j][1]) < MAX)
                        rMin[i][2] = radius;
                    if (radius < rMin[i][3] && (r1[i][1] - r1[j][1]) < 0. && fabs(r1[i][1]) < MAX &&
                        fabs(r1[j][1]) < MAX)
                        rMin[i][3] = radius;
                    if (radius < rMin[i][4] && (r1[i][2] - r1[j][2]) > 0. && fabs(r1[i][2]) < MAX &&
                        fabs(r1[j][2]) < MAX)
                        rMin[i][4] = radius;
                    if (radius < rMin[i][5] && (r1[i][2] - r1[j][2]) < 0. && fabs(r1[i][2]) < MAX &&
                        fabs(r1[j][2]) < MAX)
                        rMin[i][5] = radius;
                    if (!radius || radius < MIN || radius > 1e100 * MAX) {
                        cout << radius << endl;
                        return -1;
                    } // CRASH!
                    for (q = 0; q < FORCES; q++) {
                        if (DIM > 2)
                            a1[i][2] += G_Newton[q] * pow(radius - equil[q], powerLaw[q]) * cosP / mass[i] * flip[2];
                        a1[i][0] += G_Newton[q] * pow(radius - equil[q], powerLaw[q]) * cosT * sinP / mass[i] * flip[0];
                        a1[i][1] += G_Newton[q] * pow(radius - equil[q], powerLaw[q]) * sinT * sinP / mass[i] * flip[1];
                        if (G_Newton[q] == 6.673889e-11) {
                            a1[i][0] *= mass[j] * mass[i];
                            a1[i][1] *= mass[j] * mass[i];
                            a1[i][2] *= mass[j] * mass[i];
                        }
                    } //F=m*a across all forces
                } //complete acceleration
            } //loop over masses
        } //loop over masses higher-level

        double aAvg = 0.;
        double aNorm = 0.;
        bool out;
        long ext = 0;
        double vAvg = 0.;
        vNorm = 0.;
        aAvg /= double(BODIES);
        counter = 0;
        for (i = 0; i < BODIES; i++) {
            if (rMin[i][0] != MAX) {
                rAvg += rMin[i][0];
                counter++;
            }
            if (rMin[i][1] != MAX) {
                rAvg += rMin[i][1];
                counter++;
            }
            if (rMin[i][2] != MAX) {
                rAvg += rMin[i][2];
                counter++;
            }
            if (rMin[i][3] != MAX) {
                rAvg += rMin[i][3];
                counter++;
            }
            if (rMin[i][4] != MAX) {
                rAvg += rMin[i][4];
                counter++;
            }
            if (rMin[i][5] != MAX) {
                rAvg += rMin[i][5];
                counter++;
            }
            out = false;
            for (j = 0; j < DIM; j++) {
                v2[i][j] += a1[i][j] * dt;
                if (v2[i][j] > 295.) v2[i][j] = 295.;
                if (v2[i][j] < -295.) v2[i][j] = -295.;
                r1[i][j] += 0.5 * v2[i][j] * dt;
                if (r1[i][j] > MAX || r1[i][j] < -MAX) {
                    if (!out) ext++;
                    out = true;
                }
            }
            aNorm = sqrt(pow(a1[i][0], 2.) + pow(a1[i][1], 2.) + pow(a1[i][2], 2.));
            aAvg += aNorm;
            vNorm = sqrt(pow(v2[i][0], 2.) + pow(v2[i][1], 2.) + pow(v2[i][2], 2.));
            vAvg += vNorm;
            Norms[i] = vNorm;
            if (vNorm > 5e2 && isnan(vNorm) && out && vNorm < 4e1) { //change and's to or's to activate
                for (j = 0; j < DIM; j++) { r1[i][j] = (2. * rand_uniform() - 1.) * MAX; }
                cost = 1. - 2. * rand_uniform();
                sint = sqrt((1. - cost) * (1. + cost));
                phi = 2. * M_PI * rand_uniform();
                sinp = sin(phi);
                cosp = cos(phi);
                vx = sint * cosp;
                vy = sint * sinp;
                vz = cost;
                vNorm = VonNeumann(4e1, 5e2, 0., 15e3); // meters per sec.
                v2[i][0] = vx * vNorm;
                v2[i][1] = vy * vNorm;
                v2[i][2] = vz * vNorm;
            } //vAvg += vNorm; Norms[i]=vNorm;
            coord[i][0] = int(floor((r1[i][0] / X_FACTOR) + 0.5)) + int(floor(double(COL) / 2.));
            coord[i][1] = -int(floor((r1[i][1] / Y_FACTOR) + 0.5)) + int(floor(double(ROW) / 2.));
            double signZ;
            if (r1[i][2] != 0.) signZ = fabs(r1[i][2]) / r1[i][2]; else signZ = 1.;
            if (coord[i][0] < COL && coord[i][1] < ROW && coord[i][0] >= 0 && coord[i][1] >= 0)
                grid[coord[i][1]][coord[i][0]] = (i + 1) * signZ;
        }
        vAvg /= double(BODIES);
        rAvg /= double(counter); //cout << double(ext)/double(BODIES) << endl;

        //system("sleep 0.1");
        /*ClearTheScreen();
        for ( i = 0; i < ROW; i++ ) {
          for ( j = 0; j < COL; j++ ) {
        if ( grid[i][j] > 0 ) cout << "O";
        else if ( grid[i][j] < 0 ) cout << "o";
        else cout << " ";
          }
          cout << endl;
        }*/ //draw the new grid to the screen by row

        if (!time) printf("time [s]\tseparation [m]\tmean speed [m/s]\n");
        printf("%e\t%e\t%e\n", time + dt, rAvg, vAvg);
        //printf("%e,%f %f,%f %f,%f\t%e,%f %f,%f %f,%f\t%e\n",r1[0][0],r1[0][1],v2[0][0],v2[0][1],a1[0][0],a1[0][1],r1[1][0],r1[1][1],v2[1][0],v2[1][1],a1[1][0],a1[1][1],time);

        time += dt;
        for (i = 0; i < BODIES; i++) {
            for (j = 0; j < DIM; j++) {
                if (r1[i][j] > MAX) r1[i][j] = -MAX + (r1[i][j] - MAX);
                if (r1[i][j] < -MAX) r1[i][j] = MAX + (MAX - r1[i][j]);
                //if ( fabs(r1[i][j]) > MAX ) v2[i][j] *= -1.;
            } //over dims
        } //over N-bodies

    }

    cout << "Speeds" << endl;
    for (i = 0; i < BODIES; i++) cout << Norms[i] << endl;
    return 1;

}

double rand_uniform() {

    return (double) rand() / (double) RAND_MAX;

}

void ClearTheScreen() {

    printf("\033[%d;%dH", 1, 1);

    return;

} //hopefully this func works for all OSes!!

double VonNeumann(double xMin, double xMax, double yMin, double yMax) {
    /* Method for integration and random sampling. Using this for randomly sampling a Maxwellian distribution
     * Samples some random normal velocity
     */

    double xTry = xMin + (xMax - xMin) * rand_uniform();
    double yTry = yMin + (yMax - yMin) * rand_uniform();

    double FuncValue = pow(xTry, 2.) * exp(-2.70e-5 * pow(xTry, 2.));

    while (yTry > FuncValue) {
        xTry = xMin + (xMax - xMin) * rand_uniform();
        yTry = yMin + (yMax - yMin) * rand_uniform();
        FuncValue =
                pow(xTry, 2.) * exp(-2.70e-5 * pow(xTry, 2.));
    }

    return xTry; //selection under curve made

}