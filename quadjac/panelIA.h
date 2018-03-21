
typedef struct {
	double *v0;
	double *v1;
        double len;
        double vel;
        int idxv0;
        int idxv1;
} panel;

void panelIA(int d, int p, double h, panel *pX0, panel *pX1, panel *pY0, panel *pY1, double xLeg[], double wLeg[], double solution[]);
