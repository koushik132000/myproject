#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define CHAINS 36
#define BEADS_PER_CHAIN 20
#define FILE_COUNT 21

int main() {
    char ipname[] = "mchains";
    char opname[] = "rend";

    for (int i = 0; i < FILE_COUNT; i++) {
        char posfile[100], outfile[100];
        sprintf(posfile, "%s%d.csv", ipname, i);
        sprintf(outfile, "%s%d.csv", opname, i);

        FILE *fp = fopen(posfile, "r");
        if (!fp) {
            perror("Error reading input file");
            return 1;
        }

        // Skip header
        char line[200];
        fgets(line, sizeof(line), fp);

        // Store all x, y, z (translated +0.5)
        double xs[CHAINS * BEADS_PER_CHAIN];
        double ys[CHAINS * BEADS_PER_CHAIN];
        double zs[CHAINS * BEADS_PER_CHAIN];

        for (int j = 0; j < CHAINS * BEADS_PER_CHAIN; j++) {
            double x, y, z;
            fscanf(fp, "%lf,%lf,%lf\n", &x, &y, &z);
            xs[j] = x + 0.5;
            ys[j] = y + 0.5;
            zs[j] = z + 0.5;
        }
        fclose(fp);

        // Open output file
        FILE *fo = fopen(outfile, "w");
        if (!fo) {
            perror("Error writing output file");
            return 1;
        }

        // Compute end-to-end distances
        for (int c = 0; c < CHAINS; c++) {
            int start = c * BEADS_PER_CHAIN;
            int end = start + BEADS_PER_CHAIN - 1;

            double x0 = xs[start];
            double y0 = ys[start];
            double z0 = zs[start];

            double xl = xs[end];
            double yl = ys[end];
            double zl = zs[end];

            double dx = xl - x0;
            double dy = yl - y0;
            double dz = zl - z0;

            double dist = sqrt(dx*dx + dy*dy + dz*dz);

            fprintf(fo, "%.6lf\n", dist);
            printf("Chain %d: End-to-End Distance = %.6lf\n", c, dist);
        }

        fclose(fo);
        printf("Finished snapshot %d\n", i);
    }

    printf("All done.\n");
    return 0;
}
