#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define CHAINS 36
#define BEADS 20
#define STEPS 21

// Box dimensions
#define XMAX 7.0
#define YMAX 7.0
#define ZMAX 20.0

int main() {
    FILE *f_in, *f_out;
    char infname[50], outfname[50];
    double x[CHAINS * BEADS], y[CHAINS * BEADS], z[CHAINS * BEADS];

    for (int step = 0; step < STEPS; step++) {
        sprintf(infname, "mchains%d.csv", step);
        f_in = fopen(infname, "r");
        if (f_in == NULL) {
            printf("Error opening %s\n", infname);
            continue;
        }

        sprintf(outfname, "rend%d.csv", step);
        f_out = fopen(outfname, "w");
        if (f_out == NULL) {
            printf("Error opening %s\n", outfname);
            fclose(f_in);
            continue;
        }

        // Skip header
        char header[200];
        fgets(header, sizeof(header), f_in);

        // Read bead positions (add 0.5 shift)
        for (int i = 0; i < CHAINS * BEADS; i++) {
            double tx, ty, tz;
            fscanf(f_in, "%*d,%lf,%lf,%lf\n", &tx, &ty, &tz);
            x[i] = tx + 0.5;
            y[i] = ty + 0.5;
            z[i] = tz + 0.5;
        }
        fclose(f_in);

        // Compute end-to-end distances with PBC
        for (int c = 0; c < CHAINS; c++) {
            int start = c * BEADS;
            int end = start + BEADS - 1;

            double dx = x[end] - x[start];
            double dy = y[end] - y[start];
            double dz = z[end] - z[start];

            // Apply periodic boundary correction( wraps the coordinates to get the shortest distance)
            if (dx >  XMAX / 2.0) dx -= XMAX;
            if (dx <= -XMAX / 2.0) dx += XMAX;

            if (dy >  YMAX / 2.0) dy -= YMAX;
            if (dy <= -YMAX / 2.0) dy += YMAX;

            if (dz >  ZMAX / 2.0) dz -= ZMAX;
            if (dz <= -ZMAX / 2.0) dz += ZMAX;

            double Re = sqrt(dx * dx + dy * dy + dz * dz);
            fprintf(f_out, "%f\n", Re);
        }

        fclose(f_out);
        printf("Finished step %d\n", step);
    }

    printf("All Re values written.\n");
    return 0;
}
