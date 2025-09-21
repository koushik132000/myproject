#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define CHAINS 36
#define BEADS_PER_CHAIN 20
#define FILE_COUNT 21

int main() {
    char ipname[] = "mchains";
    char icom[] = "com";
    char opname[] = "rg";

    for (int i = 0; i < FILE_COUNT; i++) {
        double sxs[CHAINS] = {0}, sys[CHAINS] = {0}, szs[CHAINS] = {0}, crg[CHAINS] = {0};
        double comx[CHAINS], comy[CHAINS], comz[CHAINS];
        double xs[CHAINS * BEADS_PER_CHAIN], ys[CHAINS * BEADS_PER_CHAIN], zs[CHAINS * BEADS_PER_CHAIN];

        char comfile[100], posfile[100], rgfile[100];
        sprintf(comfile, "%s%d.csv", icom, i);
        sprintf(posfile, "%s%d.csv", ipname, i);
        sprintf(rgfile, "%s%d.csv", opname, i);

        // Read center of mass file
        FILE *fc = fopen(comfile, "r");
        if (!fc) {
            perror("COM file read error");
            return 1;
        }
        char line[200];
        fgets(line, sizeof(line), fc); // Skip header

        for (int c = 0; c < CHAINS; c++) {
            int id;
            fscanf(fc, "%d,%lf,%lf,%lf\n", &id, &comx[c], &comy[c], &comz[c]);
        }
        fclose(fc);

        // Read position file and shift by 0.5
        FILE *fp = fopen(posfile, "r");
        if (!fp) {
            perror("Position file read error");
            return 1;
        }
        fgets(line, sizeof(line), fp); // Skip header

        for (int j = 0; j < CHAINS * BEADS_PER_CHAIN; j++) {
            double x, y, z;
            fscanf(fp, "%lf,%lf,%lf\n", &x, &y, &z);
            xs[j] = x + 0.5;
            ys[j] = y + 0.5;
            zs[j] = z + 0.5;
        }
        fclose(fp);

        // Compute radius of gyration
        for (int c = 0; c < CHAINS; c++) {
            for (int b = 0; b < BEADS_PER_CHAIN; b++) {
                int idx = c * BEADS_PER_CHAIN + b;
                sxs[c] += (comx[c] - xs[idx]) * (comx[c] - xs[idx]);
                sys[c] += (comy[c] - ys[idx]) * (comy[c] - ys[idx]);
                szs[c] += (comz[c] - zs[idx]) * (comz[c] - zs[idx]);
            }
        }

        for (int c = 0; c < CHAINS; c++) {
            crg[c] = sqrt((sxs[c] + sys[c] + szs[c]) / BEADS_PER_CHAIN);
        }

        // Write results to output file
        FILE *fo = fopen(rgfile, "w");
        if (!fo) {
            perror("Output file write error");
            return 1;
        }
        fprintf(fo, "c,Rg\n");
        for (int c = 0; c < CHAINS; c++) {
            fprintf(fo, "%d,%.6lf\n", c, crg[c]);
        }
        fclose(fo);

        printf("Processed file set %d\n", i);
    }

    printf("All done.\n");
    return 0;
}
